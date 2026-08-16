//! High-throughput single-end orchestration over the shared combined-index search core.

use bsbit_core::alphabet::Base;
use bsbit_index::reference::ReferenceIndex;
use bsbit_index::storage::fm::ProjectedBase;

use super::mapq::{SingleMapqEvidence, single_mapping_quality_from_evidence};
use crate::AlignmentError;
use crate::placement::{ReadPlacement, placement_origin_key};
use crate::read_mapping_limits::{
    INITIAL_EDIT_DISTANCE, MAX_EDIT_DISTANCE, MAX_READ_BASES, MIN_SUFFIX_BASES,
};
use crate::search::combined_query::{CombinedSearchReferenceExt, CombinedSeedMatches};

use crate::read_mapping::{ReadAlignmentMetrics, ReadCandidate, ReadWorkspace, ungapped_distance};
use crate::search::combined_adaptive::{
    CombinedTwoLaneSearchState, DEFAULT_MAXIMUM_SEED_ROUNDS, DEFAULT_SEARCH_LIMITS,
    DIRECT_SINGLETON_PROOF, DeferredCombinedSeed, EMPTY_SEED_STEP, FLEXIBLE_NOMINAL_PROOF,
    INITIAL_SEARCH_LIMITS, combined_seed_round_is_locatable, continue_combined_two_lane_search,
    prepare_combined_projection,
};

/// Final classification for one directional single read.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SingleMappingStatus {
    /// No verified placement survived the bounded search.
    Unmapped,
    /// Exactly one best biological origin survived.
    Unique,
    /// Multiple equally good biological origins survived.
    Ambiguous,
}

/// Final mapping facts for one directional single read.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SingleAlignmentResult {
    status: SingleMappingStatus,
    placement: Option<ReadPlacement>,
    mapping_quality: u8,
    located_rows: u64,
    verified_placements: u64,
}

impl SingleAlignmentResult {
    /// Returns the final mapping class.
    #[must_use]
    pub const fn status(self) -> SingleMappingStatus {
        self.status
    }

    /// Returns the deterministic representative placement, when mapped.
    #[must_use]
    pub const fn placement(self) -> Option<ReadPlacement> {
        self.placement
    }

    /// Returns the calibrated SAM mapping quality, or zero when not unique.
    #[must_use]
    pub const fn mapping_quality(self) -> u8 {
        self.mapping_quality
    }

    /// Returns the suffix rows located by the candidate search.
    #[must_use]
    pub const fn located_rows(self) -> u64 {
        self.located_rows
    }

    /// Returns the complete verified placement count before best-tier selection.
    #[must_use]
    pub const fn verified_placements(self) -> u64 {
        self.verified_placements
    }

    const fn unmapped(located_rows: u64, verified_placements: u64) -> Self {
        Self {
            status: SingleMappingStatus::Unmapped,
            placement: None,
            mapping_quality: 0,
            located_rows,
            verified_placements,
        }
    }
}

#[derive(Clone, Copy)]
struct PreparedSingleRead<'a> {
    read: &'a [Base],
    projection: &'a [ProjectedBase],
    first_seed: Option<CombinedSeedMatches>,
    search: &'a CombinedTwoLaneSearchState,
    initial_candidates: &'a [ReadCandidate],
}

#[derive(Clone, Copy)]
struct SingleResultEvidence<'a> {
    read_length: usize,
    metrics: ReadAlignmentMetrics,
    verified_distance_limit: u8,
    first_seed: Option<CombinedSeedMatches>,
    search: &'a CombinedTwoLaneSearchState,
    lane: usize,
}

/// Worker-owned batch storage for directional single-read alignment.
pub struct SingleBatchAligner {
    reads: [ReadWorkspace; 2],
    projections: Vec<[ProjectedBase; MAX_READ_BASES]>,
    searchable_reads: Vec<bool>,
    first_seeds: Vec<Option<CombinedSeedMatches>>,
    round_matches: Vec<Option<CombinedSeedMatches>>,
    search_states: Vec<CombinedTwoLaneSearchState>,
    initial_candidates: Vec<Vec<ReadCandidate>>,
    origins: Vec<(u64, bsbit_core::bisulfite::BisulfiteStrand, i128)>,
    results: Vec<SingleAlignmentResult>,
}

impl SingleBatchAligner {
    /// Allocates reusable storage for at least `read_capacity` reads.
    #[must_use]
    pub fn with_capacity(read_capacity: usize) -> Self {
        Self {
            reads: core::array::from_fn(|_| ReadWorkspace::with_capacity(4096, 1024)),
            projections: Vec::with_capacity(read_capacity),
            searchable_reads: Vec::with_capacity(read_capacity),
            first_seeds: Vec::with_capacity(read_capacity),
            round_matches: Vec::with_capacity(read_capacity),
            search_states: Vec::with_capacity(read_capacity),
            initial_candidates: Vec::with_capacity(read_capacity),
            origins: Vec::with_capacity(64),
            results: Vec::with_capacity(read_capacity),
        }
    }

    /// Maps a batch through the same persisted combined index, adaptive seed
    /// schedule, and d3/d5 verifier used by paired-end alignment.
    ///
    /// # Errors
    ///
    /// Returns an unsupported read/edit domain or combined-index failure.
    #[allow(clippy::too_many_lines)]
    pub fn map_reads<'a>(
        &'a mut self,
        reference: &ReferenceIndex,
        reads: &[&[Base]],
        maximum_edit_distance: u8,
    ) -> Result<&'a [SingleAlignmentResult], AlignmentError> {
        if maximum_edit_distance > MAX_EDIT_DISTANCE {
            return Err(AlignmentError::UnsupportedEditDistance {
                requested: maximum_edit_distance,
                maximum: MAX_EDIT_DISTANCE,
            });
        }
        if reads.len() > 64 {
            return Err(AlignmentError::LocatedCountOverflow);
        }
        self.projections.clear();
        self.projections
            .resize(reads.len(), [ProjectedBase::A; MAX_READ_BASES]);
        self.searchable_reads.clear();
        for (projection, read) in self.projections.iter_mut().zip(reads) {
            let searchable = read.len() >= MIN_SUFFIX_BASES
                && read.iter().filter(|base| base.is_unknown()).count()
                    <= usize::from(maximum_edit_distance);
            self.searchable_reads.push(searchable);
            if searchable {
                prepare_combined_projection(read, false, projection)?;
            }
        }
        let searchable = self
            .projections
            .iter()
            .zip(reads)
            .zip(&self.searchable_reads)
            .filter(|(_, searchable)| **searchable)
            .map(|((projection, read), _)| &projection[..read.len()])
            .collect::<Vec<_>>();
        let searchable_seeds = if searchable.is_empty() {
            Vec::new()
        } else {
            reference
                .combined_maximal_suffix_projected_wavefront(&searchable, MIN_SUFFIX_BASES)
                .map_err(|_| AlignmentError::CombinedIndex)?
        };
        self.first_seeds.clear();
        let mut searchable_ordinal = 0_usize;
        for &searchable in &self.searchable_reads {
            if searchable {
                self.first_seeds.push(searchable_seeds[searchable_ordinal]);
                searchable_ordinal += 1;
            } else {
                self.first_seeds.push(None);
            }
        }
        self.prepare_search_wavefront(reference, reads)?;
        self.results.clear();
        let mut ordinal = 0_usize;
        while ordinal < reads.len() {
            if !self.searchable_reads[ordinal] {
                self.results.push(SingleAlignmentResult::unmapped(0, 0));
                ordinal += 1;
                continue;
            }
            if ordinal + 1 < reads.len() && self.searchable_reads[ordinal + 1] {
                let pair = [reads[ordinal], reads[ordinal + 1]];
                let projections = [self.projections[ordinal], self.projections[ordinal + 1]];
                let first_seeds = [self.first_seeds[ordinal], self.first_seeds[ordinal + 1]];
                let searches = [self.search_states[ordinal], self.search_states[ordinal + 1]];
                let initial_candidates = [
                    core::mem::take(&mut self.initial_candidates[ordinal]),
                    core::mem::take(&mut self.initial_candidates[ordinal + 1]),
                ];
                let prepared = core::array::from_fn(|lane| PreparedSingleRead {
                    read: pair[lane],
                    projection: &projections[lane][..pair[lane].len()],
                    first_seed: first_seeds[lane],
                    search: &searches[lane],
                    initial_candidates: &initial_candidates[lane],
                });
                let mapped = self.map_two(reference, prepared, maximum_edit_distance);
                let [first_candidates, second_candidates] = initial_candidates;
                self.initial_candidates[ordinal] = first_candidates;
                self.initial_candidates[ordinal + 1] = second_candidates;
                let mapped = mapped?;
                self.results.extend(mapped);
                ordinal += 2;
            } else {
                let read = reads[ordinal];
                let projection = self.projections[ordinal];
                let first_seed = self.first_seeds[ordinal];
                let search = self.search_states[ordinal];
                let initial_candidates = core::mem::take(&mut self.initial_candidates[ordinal]);
                let result = self.map_one(
                    reference,
                    PreparedSingleRead {
                        read,
                        projection: &projection[..read.len()],
                        first_seed,
                        search: &search,
                        initial_candidates: &initial_candidates,
                    },
                    maximum_edit_distance,
                );
                self.initial_candidates[ordinal] = initial_candidates;
                let result = result?;
                self.results.push(result);
                ordinal += 1;
            }
        }
        Ok(&self.results)
    }

    #[allow(clippy::too_many_lines)]
    fn prepare_search_wavefront(
        &mut self,
        reference: &ReferenceIndex,
        reads: &[&[Base]],
    ) -> Result<(), AlignmentError> {
        const MAX_LANES: usize = 64;
        self.round_matches.clear();
        self.round_matches.resize(reads.len(), None);
        self.search_states.clear();
        self.search_states
            .resize(reads.len(), CombinedTwoLaneSearchState::new());
        self.initial_candidates.resize_with(reads.len(), Vec::new);
        self.initial_candidates.truncate(reads.len());
        for (ordinal, candidates) in self.initial_candidates.iter_mut().enumerate() {
            candidates.clear();
            let state = &mut self.search_states[ordinal];
            state.initialized = true;
            state.active = [self.searchable_reads[ordinal], false];
        }

        let mut active_ordinals = [0_usize; MAX_LANES];
        let mut available_bases = [0_usize; MAX_LANES];
        let mut compact_matches = [None; MAX_LANES];
        let mut locatable_ordinals = [0_usize; MAX_LANES];

        for round in 0..INITIAL_SEARCH_LIMITS.maximum_seed_rounds {
            let mut active_count = 0_usize;
            for (ordinal, state) in self.search_states.iter_mut().enumerate() {
                let available = reads[ordinal].len().saturating_sub(state.offsets[0]);
                available_bases[ordinal] = available;
                state.active[0] &= available >= MIN_SUFFIX_BASES;
                state.completed_rounds = round + 1;
                if state.active[0] {
                    active_ordinals[active_count] = ordinal;
                    active_count += 1;
                }
            }
            if active_count == 0 {
                break;
            }

            self.round_matches.fill(None);
            if round == 0 {
                for &ordinal in &active_ordinals[..active_count] {
                    self.round_matches[ordinal] = self.first_seeds[ordinal];
                }
            } else {
                let first = active_ordinals[0];
                let mut patterns = [&self.projections[first][..available_bases[first]]; MAX_LANES];
                for (slot, &ordinal) in active_ordinals[..active_count].iter().enumerate() {
                    patterns[slot] = &self.projections[ordinal][..available_bases[ordinal]];
                }
                compact_matches[..active_count].fill(None);
                reference
                    .combined_maximal_suffix_projected_wavefront_into(
                        &patterns[..active_count],
                        MIN_SUFFIX_BASES,
                        &mut compact_matches[..active_count],
                    )
                    .map_err(|_| AlignmentError::CombinedIndex)?;
                for (slot, &ordinal) in active_ordinals[..active_count].iter().enumerate() {
                    self.round_matches[ordinal] = compact_matches[slot];
                }
            }

            if INITIAL_SEARCH_LIMITS.maximum_seed_rounds < DEFAULT_MAXIMUM_SEED_ROUNDS {
                let default_limits = DEFAULT_SEARCH_LIMITS;
                for &ordinal in &active_ordinals[..active_count] {
                    if let Some(seed) = self.round_matches[ordinal]
                        && !combined_seed_round_is_locatable(seed, INITIAL_SEARCH_LIMITS)
                        && combined_seed_round_is_locatable(seed, default_limits)
                    {
                        let offset = self.search_states[ordinal].offsets[0];
                        self.search_states[ordinal].defer(
                            0,
                            DeferredCombinedSeed {
                                matches: seed,
                                offset,
                                round,
                            },
                        );
                    }
                }
            }

            let mut locatable_count = 0_usize;
            for &ordinal in &active_ordinals[..active_count] {
                if self.round_matches[ordinal].is_some_and(|seed| {
                    combined_seed_round_is_locatable(seed, INITIAL_SEARCH_LIMITS)
                }) {
                    locatable_ordinals[locatable_count] = ordinal;
                    locatable_count += 1;
                }
            }

            for group in locatable_ordinals[..locatable_count].chunks(2) {
                if group.len() == 2 {
                    let ordinals = [group[0], group[1]];
                    let matches = ordinals
                        .map(|ordinal| self.round_matches[ordinal].expect("locatable seed exists"));
                    let mut direct = [false; 2];
                    let metrics = reference
                        .visit_combined_seed_two_lanes_complete(
                            matches,
                            ordinals.map(|ordinal| {
                                u64::try_from(self.search_states[ordinal].offsets[0])
                                    .unwrap_or(u64::MAX)
                            }),
                            ordinals.map(|ordinal| {
                                u64::try_from(reads[ordinal].len()).unwrap_or(u64::MAX)
                            }),
                            &mut |lane, hit| {
                                let ordinal = ordinals[lane];
                                let mut candidate = ReadCandidate {
                                    contig_ordinal: hit.contig_ordinal(),
                                    start: hit.start(),
                                    strand: hit.strand(),
                                    proof_mask: FLEXIBLE_NOMINAL_PROOF | (1_u8 << round),
                                };
                                if round == 0
                                    && matches[lane].exact_hit_count() == 1
                                    && let Some(distance) =
                                        ungapped_distance(reference, reads[ordinal], candidate)
                                {
                                    candidate.proof_mask = DIRECT_SINGLETON_PROOF | distance;
                                    direct[lane] = true;
                                    self.initial_candidates[ordinal].clear();
                                }
                                self.initial_candidates[ordinal].push(candidate);
                            },
                        )
                        .map_err(|_| AlignmentError::CombinedIndex)?;
                    for lane in 0..2 {
                        let ordinal = ordinals[lane];
                        self.search_states[ordinal].located[0] = self.search_states[ordinal]
                            .located[0]
                            .checked_add(metrics[lane].located_coordinates())
                            .ok_or(AlignmentError::LocatedCountOverflow)?;
                        self.search_states[ordinal].active[0] &= !direct[lane];
                        self.search_states[ordinal].direct[0] |= direct[lane];
                    }
                } else {
                    let ordinal = group[0];
                    let matches =
                        self.round_matches[ordinal].expect("locatable singleton seed exists");
                    let mut direct = false;
                    let metrics = reference
                        .visit_combined_seed(
                            matches,
                            u64::try_from(self.search_states[ordinal].offsets[0])
                                .unwrap_or(u64::MAX),
                            u64::try_from(reads[ordinal].len()).unwrap_or(u64::MAX),
                            &mut |hit| {
                                let mut candidate = ReadCandidate {
                                    contig_ordinal: hit.contig_ordinal(),
                                    start: hit.start(),
                                    strand: hit.strand(),
                                    proof_mask: FLEXIBLE_NOMINAL_PROOF | (1_u8 << round),
                                };
                                if round == 0
                                    && matches.exact_hit_count() == 1
                                    && let Some(distance) =
                                        ungapped_distance(reference, reads[ordinal], candidate)
                                {
                                    candidate.proof_mask = DIRECT_SINGLETON_PROOF | distance;
                                    direct = true;
                                    self.initial_candidates[ordinal].clear();
                                }
                                self.initial_candidates[ordinal].push(candidate);
                                true
                            },
                        )
                        .map_err(|_| AlignmentError::CombinedIndex)?;
                    self.search_states[ordinal].located[0] = self.search_states[ordinal].located[0]
                        .checked_add(metrics.located_coordinates())
                        .ok_or(AlignmentError::LocatedCountOverflow)?;
                    self.search_states[ordinal].active[0] &= !direct;
                    self.search_states[ordinal].direct[0] |= direct;
                }
            }

            for &ordinal in &active_ordinals[..active_count] {
                if let Some(seed) = self.round_matches[ordinal] {
                    let matched = usize::try_from(seed.matched_bases())
                        .map_err(|_| AlignmentError::LocatedCountOverflow)?;
                    self.search_states[ordinal].offsets[0] = self.search_states[ordinal].offsets[0]
                        .saturating_add((matched.saturating_mul(3) / 4).max(1));
                } else {
                    self.search_states[ordinal].offsets[0] =
                        self.search_states[ordinal].offsets[0].saturating_add(EMPTY_SEED_STEP);
                }
            }
        }
        Ok(())
    }

    fn map_one(
        &mut self,
        reference: &ReferenceIndex,
        prepared: PreparedSingleRead<'_>,
        maximum_edit_distance: u8,
    ) -> Result<SingleAlignmentResult, AlignmentError> {
        let PreparedSingleRead {
            read,
            projection,
            first_seed,
            search,
            initial_candidates,
        } = prepared;
        let mut search = *search;
        let (workspaces, origins) = (&mut self.reads, &mut self.origins);
        let [read_workspace, unused_workspace] = workspaces;
        read_workspace.begin_verification_cache_read();
        read_workspace.candidates.clear();
        read_workspace.candidate_nominals.clear();
        read_workspace
            .candidate_nominals
            .extend_from_slice(initial_candidates);
        read_workspace.placements.clear();
        unused_workspace.candidate_nominals.clear();
        let mut metrics = ReadAlignmentMetrics {
            located_rows: search.located[0],
            ..ReadAlignmentMetrics::default()
        };
        let initial_budget = maximum_edit_distance.min(INITIAL_EDIT_DISTANCE);
        let (_, observed) = read_workspace.verify_candidates_with_budget(
            reference,
            read,
            metrics,
            initial_budget,
        )?;
        metrics = observed;
        if !read_workspace.placements.is_empty() {
            return Ok(Self::finish_result(
                read_workspace,
                origins,
                SingleResultEvidence {
                    read_length: read.len(),
                    metrics,
                    verified_distance_limit: initial_budget,
                    first_seed,
                    search: &search,
                    lane: 0,
                },
            ));
        }

        if maximum_edit_distance > initial_budget {
            read_workspace.candidates.clear();
            read_workspace.placements.clear();
            let (_, observed) = read_workspace.verify_candidates_with_budget(
                reference,
                read,
                metrics,
                maximum_edit_distance,
            )?;
            metrics = observed;
            if !read_workspace.placements.is_empty() {
                return Ok(Self::finish_result(
                    read_workspace,
                    origins,
                    SingleResultEvidence {
                        read_length: read.len(),
                        metrics,
                        verified_distance_limit: maximum_edit_distance,
                        first_seed,
                        search: &search,
                        lane: 0,
                    },
                ));
            }
        }

        read_workspace.candidates.clear();
        read_workspace.placements.clear();
        unused_workspace.candidate_nominals.clear();
        let additional = continue_combined_two_lane_search(
            reference,
            [read, &[]],
            [projection, &[]],
            false,
            &mut search,
            &mut read_workspace.candidate_nominals,
            &mut unused_workspace.candidate_nominals,
        )?;
        metrics.located_rows = metrics.located_rows.saturating_add(additional[0]);
        let (_, observed) = read_workspace.verify_candidates_with_budget(
            reference,
            read,
            metrics,
            maximum_edit_distance,
        )?;
        Ok(Self::finish_result(
            read_workspace,
            origins,
            SingleResultEvidence {
                read_length: read.len(),
                metrics: observed,
                verified_distance_limit: maximum_edit_distance,
                first_seed,
                search: &search,
                lane: 0,
            },
        ))
    }

    #[allow(clippy::too_many_lines)]
    fn map_two(
        &mut self,
        reference: &ReferenceIndex,
        prepared: [PreparedSingleRead<'_>; 2],
        maximum_edit_distance: u8,
    ) -> Result<[SingleAlignmentResult; 2], AlignmentError> {
        let reads = prepared.map(|input| input.read);
        let projections = prepared.map(|input| input.projection);
        let first_seeds = prepared.map(|input| input.first_seed);
        let prepared_searches = prepared.map(|input| input.search);
        let initial_candidates = prepared.map(|input| input.initial_candidates);
        let (workspaces, origins) = (&mut self.reads, &mut self.origins);
        for (lane, workspace) in workspaces.iter_mut().enumerate() {
            workspace.begin_verification_cache_read();
            workspace.candidates.clear();
            workspace.candidate_nominals.clear();
            workspace
                .candidate_nominals
                .extend_from_slice(initial_candidates[lane]);
            workspace.placements.clear();
        }
        let mut search = CombinedTwoLaneSearchState::new();
        search.initialized = true;
        search.completed_rounds = prepared_searches
            .iter()
            .map(|state| state.completed_rounds)
            .max()
            .unwrap_or(0);
        for (lane, prepared_search) in prepared_searches.iter().enumerate() {
            search.located[lane] = prepared_search.located[0];
            search.offsets[lane] = prepared_search.offsets[0];
            search.active[lane] = prepared_search.active[0];
            search.direct[lane] = prepared_search.direct[0];
            search.deferred[lane] = prepared_search.deferred[0];
            search.deferred_len[lane] = prepared_search.deferred_len[0];
        }
        let mut metrics = search.located.map(|located_rows| ReadAlignmentMetrics {
            located_rows,
            ..ReadAlignmentMetrics::default()
        });
        let initial_budget = maximum_edit_distance.min(INITIAL_EDIT_DISTANCE);
        let mut results = [None; 2];

        for lane in 0..2 {
            let (_, observed) = workspaces[lane].verify_candidates_with_budget(
                reference,
                reads[lane],
                metrics[lane],
                initial_budget,
            )?;
            metrics[lane] = observed;
            if !workspaces[lane].placements.is_empty() {
                results[lane] = Some(Self::finish_result(
                    &workspaces[lane],
                    origins,
                    SingleResultEvidence {
                        read_length: reads[lane].len(),
                        metrics: metrics[lane],
                        verified_distance_limit: initial_budget,
                        first_seed: first_seeds[lane],
                        search: &search,
                        lane,
                    },
                ));
            }
        }

        if maximum_edit_distance > initial_budget {
            for lane in 0..2 {
                if results[lane].is_some() {
                    continue;
                }
                workspaces[lane].candidates.clear();
                workspaces[lane].placements.clear();
                let (_, observed) = workspaces[lane].verify_candidates_with_budget(
                    reference,
                    reads[lane],
                    metrics[lane],
                    maximum_edit_distance,
                )?;
                metrics[lane] = observed;
                if !workspaces[lane].placements.is_empty() {
                    results[lane] = Some(Self::finish_result(
                        &workspaces[lane],
                        origins,
                        SingleResultEvidence {
                            read_length: reads[lane].len(),
                            metrics: metrics[lane],
                            verified_distance_limit: maximum_edit_distance,
                            first_seed: first_seeds[lane],
                            search: &search,
                            lane,
                        },
                    ));
                }
            }
        }

        if results.iter().all(Option::is_some) {
            return Ok(results.map(|result| result.expect("both single reads resolved")));
        }
        for lane in 0..2 {
            if results[lane].is_some() {
                search.active[lane] = false;
                search.deferred_len[lane] = 0;
            } else {
                workspaces[lane].candidates.clear();
                workspaces[lane].placements.clear();
            }
        }
        let [first_workspace, second_workspace] = workspaces;
        let additional = continue_combined_two_lane_search(
            reference,
            reads,
            projections,
            false,
            &mut search,
            &mut first_workspace.candidate_nominals,
            &mut second_workspace.candidate_nominals,
        )?;
        for lane in 0..2 {
            if results[lane].is_some() {
                continue;
            }
            metrics[lane].located_rows =
                metrics[lane].located_rows.saturating_add(additional[lane]);
            let (_, observed) = workspaces[lane].verify_candidates_with_budget(
                reference,
                reads[lane],
                metrics[lane],
                maximum_edit_distance,
            )?;
            results[lane] = Some(Self::finish_result(
                &workspaces[lane],
                origins,
                SingleResultEvidence {
                    read_length: reads[lane].len(),
                    metrics: observed,
                    verified_distance_limit: maximum_edit_distance,
                    first_seed: first_seeds[lane],
                    search: &search,
                    lane,
                },
            ));
        }
        Ok(results.map(|result| result.expect("single-read continuation resolves every lane")))
    }

    fn finish_result(
        workspace: &ReadWorkspace,
        origins: &mut Vec<(u64, bsbit_core::bisulfite::BisulfiteStrand, i128)>,
        evidence: SingleResultEvidence<'_>,
    ) -> SingleAlignmentResult {
        let SingleResultEvidence {
            read_length,
            metrics,
            verified_distance_limit,
            first_seed,
            search,
            lane,
        } = evidence;
        let Some(best_distance) = workspace
            .placements
            .iter()
            .map(|value| value.distance())
            .min()
        else {
            return SingleAlignmentResult::unmapped(
                metrics.located_rows,
                metrics.verified_placements,
            );
        };
        origins.clear();
        let mut representative = None;
        for placement in workspace
            .placements
            .iter()
            .copied()
            .filter(|placement| placement.distance() == best_distance)
        {
            representative = Some(
                representative.map_or(placement, |current: ReadPlacement| current.min(placement)),
            );
            origins.push(placement_origin_key(placement, read_length));
        }
        origins.sort_unstable();
        origins.dedup();
        let status = if origins.len() == 1 {
            SingleMappingStatus::Unique
        } else {
            SingleMappingStatus::Ambiguous
        };
        let mapping_quality = if matches!(status, SingleMappingStatus::Unique) {
            let best_origin = origins[0];
            let second_best_distance = workspace
                .placements
                .iter()
                .copied()
                .filter(|placement| placement_origin_key(*placement, read_length) != best_origin)
                .map(ReadPlacement::distance)
                .min();
            let (first_seed_hits, first_seed_bases) = first_seed.map_or((0, 0), |seed| {
                (seed.exact_hit_count(), seed.matched_bases())
            });
            single_mapping_quality_from_evidence(SingleMapqEvidence {
                best_distance,
                second_best_distance,
                verified_distance_limit,
                located_rows: metrics.located_rows,
                distinct_candidate_starts: metrics.distinct_candidate_starts,
                verified_placements: metrics.verified_placements,
                first_seed_hits,
                first_seed_bases,
                direct_singleton: search.direct[lane],
            })
        } else {
            0
        };
        SingleAlignmentResult {
            status,
            placement: representative,
            mapping_quality,
            located_rows: metrics.located_rows,
            verified_placements: metrics.verified_placements,
        }
    }
}
