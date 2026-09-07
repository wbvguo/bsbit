//! High-throughput single-end orchestration over the shared combined-index search core.

use bsbit_core::alphabet::Base;
use bsbit_core::bisulfite::BisulfiteStrand;
use bsbit_index::reference::ReferenceIndex;
use bsbit_index::storage::fm::ProjectedBase;
use core::cmp::Reverse;

use super::SINGLE_ALIGNMENT_BATCH_SIZE;
use super::mapq::{
    SingleMapqEvidence, adapter_consensus_mapping_quality, affine_rerank_origin_count_supported,
    affine_unique_mapping_quality, default_confidence_cross_check_required,
    local_locus_mapping_quality, sensitive_replacement_certified,
};
use crate::adapter::supported_three_prime_adapter_start;
use crate::alignment_policy::{
    ADAPTER_STABILITY_DELTA, CombinedSearchLimits, DEFAULT_HIGH_CONFIDENCE_AUDIT_SEARCH_LIMITS,
    EMPTY_SEED_STEP, INITIAL_SEARCH_LIMITS, MIN_ADAPTER_RETAINED_BASES,
};
use crate::library::{ConversionPass, LibraryProfile};
#[cfg(test)]
use crate::placement::placement_conversion_counts;
use crate::placement::{ReadPlacement, placement_net_gap_bases, placement_origin_key};
use crate::read_mapping_limits::{
    INITIAL_EDIT_DISTANCE, MAX_EDIT_DISTANCE, MAX_READ_BASES, MIN_SUFFIX_BASES,
};
use crate::search::combined_query::{CombinedSearchReferenceExt, CombinedSeedMatches};
use crate::{AlignmentError, AlignmentOutputPolicy};

use crate::read_mapping::{
    ReadAlignmentMetrics, ReadCandidate, ReadWorkspace, placement_proof_mask, ungapped_distance,
};
use crate::reporting_tie_break::ReportingTieBreak;
use crate::search::combined_adaptive::{
    CombinedTwoLaneSearchState, DIRECT_SINGLETON_PROOF, DeferredCombinedSeed,
    FLEXIBLE_NOMINAL_PROOF, combined_seed_round_is_locatable, continue_combined_two_lane_search,
    continue_combined_two_lane_search_with_limits, direct_singleton_proof,
    prepare_combined_projection, seed_round_support,
};
use crate::verification::affine::{AffineScoreWorkspace, affine_placement_score};

use super::merge::{
    local_origin_locus_count, merge_non_directional_results_with_tie_break,
    origins_share_local_locus, prefer_fair_ambiguous_representative,
    prefer_minimum_net_gap_representative,
};
#[cfg(test)]
use super::merge::{merge_non_directional_completed_frontiers, merge_non_directional_results};
use super::options::SingleSearchMode;
use super::result::{SingleAlignmentResult, SingleMappingStatus};

#[derive(Clone, Copy)]
struct PreparedSingleRead<'a> {
    read: &'a [Base],
    projection: &'a [ProjectedBase],
    first_seed: Option<CombinedSeedMatches>,
    search: &'a CombinedTwoLaneSearchState,
    initial_candidates: &'a [ReadCandidate],
}

#[derive(Clone, Copy)]
struct SingleResultEvidence {
    read_length: usize,
    metrics: ReadAlignmentMetrics,
    verified_distance_limit: u8,
    first_seed: Option<CombinedSeedMatches>,
    frontier_complete: bool,
    confidence_audit_complete: bool,
}

#[derive(Clone, Copy)]
struct SingleAdapterFallback {
    result: SingleAlignmentResult,
    stability_result: Option<SingleAlignmentResult>,
    final_status: SingleMappingStatus,
    retained_end: usize,
}

#[derive(Clone, Copy)]
struct ReportingTieBreakBatch<'a> {
    seed: u64,
    read_keys: &'a [u64],
}

impl ReportingTieBreakBatch<'_> {
    fn for_read(self, ordinal: usize) -> ReportingTieBreak {
        ReportingTieBreak {
            seed: self.seed,
            read_key: self.read_keys[ordinal],
        }
    }
}

/// Worker-owned batch storage for single-read alignment.
pub struct SingleBatchAligner {
    reads: [ReadWorkspace; 2],
    projections: Vec<[ProjectedBase; MAX_READ_BASES]>,
    searchable_reads: Vec<bool>,
    first_seeds: Vec<Option<CombinedSeedMatches>>,
    round_matches: Vec<Option<CombinedSeedMatches>>,
    search_states: Vec<CombinedTwoLaneSearchState>,
    initial_candidates: Vec<Vec<ReadCandidate>>,
    origins: Vec<(u64, BisulfiteStrand, i128)>,
    results: Vec<SingleAlignmentResult>,
    primary_pass_results: Vec<SingleAlignmentResult>,
    output_results: Vec<SingleAlignmentResult>,
    output_policy: AlignmentOutputPolicy,
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
            primary_pass_results: Vec::with_capacity(read_capacity),
            output_results: Vec::with_capacity(read_capacity),
            output_policy: AlignmentOutputPolicy::default(),
        }
    }

    /// Configures adapter recognition and endpoint soft clipping for this worker.
    #[must_use]
    pub fn with_output_policy(mut self, policy: AlignmentOutputPolicy) -> Self {
        self.output_policy = policy;
        self
    }

    /// Maps a batch through the conversion passes selected by one shared
    /// single-end/paired-end library profile.
    ///
    /// Directional mode executes the original OT/OB pass. Non-directional mode
    /// additionally executes the complementary CTOT/CTOB pass and reduces both
    /// result sets under the single-end global-placement policy.
    ///
    /// # Errors
    ///
    /// Returns an unsupported read/edit domain or combined-index failure.
    pub fn map_reads_with_mode<'a>(
        &'a mut self,
        reference: &ReferenceIndex,
        reads: &[&[Base]],
        maximum_edit_distance: u8,
        library_profile: LibraryProfile,
        search_mode: SingleSearchMode,
    ) -> Result<&'a [SingleAlignmentResult], AlignmentError> {
        self.map_reads_with_mode_impl(
            reference,
            reads,
            maximum_edit_distance,
            library_profile,
            search_mode,
            None,
        )
    }

    fn map_reads_with_mode_impl<'a>(
        &'a mut self,
        reference: &ReferenceIndex,
        reads: &[&[Base]],
        maximum_edit_distance: u8,
        library_profile: LibraryProfile,
        search_mode: SingleSearchMode,
        tie_break: Option<ReportingTieBreakBatch<'_>>,
    ) -> Result<&'a [SingleAlignmentResult], AlignmentError> {
        let (first, second) = match library_profile.conversion_passes() {
            [first] => (*first, None),
            [first, second] => (*first, Some(*second)),
            _ => unreachable!("a library profile has one or two conversion passes"),
        };
        self.map_reads_pass(
            reference,
            reads,
            maximum_edit_distance,
            search_mode,
            first,
            tie_break,
        )?;
        let Some(second) = second else {
            return Ok(&self.results);
        };

        std::mem::swap(&mut self.primary_pass_results, &mut self.results);
        self.map_reads_pass(
            reference,
            reads,
            maximum_edit_distance,
            search_mode,
            second,
            tie_break,
        )?;
        for (ordinal, (complementary, original)) in self
            .results
            .iter_mut()
            .zip(&self.primary_pass_results)
            .enumerate()
        {
            *complementary = merge_non_directional_results_with_tie_break(
                tie_break.map(|batch| (reference, reads[ordinal], batch.for_read(ordinal))),
                *original,
                *complementary,
                search_mode.completes_candidate_frontier(),
            )?;
        }
        Ok(&self.results)
    }

    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn map_reads_pass<'a>(
        &'a mut self,
        reference: &ReferenceIndex,
        reads: &[&[Base]],
        maximum_edit_distance: u8,
        search_mode: SingleSearchMode,
        conversion_pass: ConversionPass,
        tie_break: Option<ReportingTieBreakBatch<'_>>,
    ) -> Result<&'a [SingleAlignmentResult], AlignmentError> {
        if maximum_edit_distance > MAX_EDIT_DISTANCE {
            return Err(AlignmentError::UnsupportedEditDistance {
                requested: maximum_edit_distance,
                maximum: MAX_EDIT_DISTANCE,
            });
        }
        if reads.len() > SINGLE_ALIGNMENT_BATCH_SIZE {
            return Err(AlignmentError::SearchBatchSize {
                observed: reads.len(),
                maximum: SINGLE_ALIGNMENT_BATCH_SIZE,
            });
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
                prepare_combined_projection(read, conversion_pass, projection)?;
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
        self.prepare_search_wavefront(
            reference,
            reads,
            search_mode.limits(),
            conversion_pass,
            false,
        )?;
        self.results.clear();
        let mut ordinal = 0_usize;
        while ordinal < reads.len() {
            if !self.searchable_reads[ordinal] {
                self.results
                    .push(SingleAlignmentResult::unmapped(reads[ordinal].len(), 0, 0));
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
                let mapped = self.map_two(
                    reference,
                    prepared,
                    maximum_edit_distance,
                    search_mode,
                    conversion_pass,
                );
                let [first_candidates, second_candidates] = initial_candidates;
                self.initial_candidates[ordinal] = first_candidates;
                self.initial_candidates[ordinal + 1] = second_candidates;
                let mut mapped = mapped?;
                if let Some(batch) = tie_break {
                    for (lane, result) in mapped.iter_mut().enumerate() {
                        *result = prefer_fair_ambiguous_representative(
                            reference,
                            pair[lane],
                            &self.reads[lane],
                            *result,
                            batch.for_read(ordinal + lane),
                        )?;
                    }
                }
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
                    search_mode,
                    conversion_pass,
                );
                self.initial_candidates[ordinal] = initial_candidates;
                let mut result = result?;
                if let Some(batch) = tie_break {
                    result = prefer_fair_ambiguous_representative(
                        reference,
                        read,
                        &self.reads[0],
                        result,
                        batch.for_read(ordinal),
                    )?;
                }
                self.results.push(result);
                ordinal += 1;
            }
        }
        Ok(&self.results)
    }

    /// Maps complete reads under one shared library profile and applies the
    /// qualified single-end output policy.
    ///
    /// Reads with exact supported Illumina adapter evidence enter the compact
    /// trimmed remap. A tentative unique recovery must remain unique at the
    /// same strand-aware biological origin after an additional eight-base
    /// shortening. An otherwise-unmapped read may recover with MAPQ capped at
    /// 20. An already mapped read may only change its reported endpoint at the
    /// same biological origin, with classification and MAPQ frozen. Adapter
    /// Recovery uses the selected library profile for the complete, trimmed,
    /// and stability passes. Non-directional reads therefore make the same
    /// global four-strand decision at every stage before the endpoint policy
    /// compares biological origins.
    ///
    /// # Errors
    ///
    /// Returns an unsupported read/edit domain or combined-index failure from
    /// any primary, trimmed, or stability mapping phase.
    ///
    /// # Panics
    ///
    /// Panics only if internally generated stability metadata loses its
    /// matching adapter result, which violates this method's construction
    /// invariant.
    #[allow(clippy::too_many_lines)]
    pub fn map_reads_for_output<'a>(
        &'a mut self,
        reference: &ReferenceIndex,
        reads: &[&[Base]],
        maximum_edit_distance: u8,
        library_profile: LibraryProfile,
        search_mode: SingleSearchMode,
    ) -> Result<&'a [SingleAlignmentResult], AlignmentError> {
        self.map_reads_for_output_impl(
            reference,
            reads,
            maximum_edit_distance,
            library_profile,
            search_mode,
            None,
        )
    }

    /// Maps complete reads and resolves only otherwise-equal MAPQ-zero
    /// reporting coordinates with a stable, caller-seeded hash lottery.
    ///
    /// `read_keys` must contain one stable input identity per read. The hash
    /// uses the seed, read key, contig name, strand, and biological coordinate;
    /// it never uses a reference contig ordinal or candidate enumeration order.
    ///
    /// # Errors
    ///
    /// Returns [`AlignmentError::ReportingTieBreakKeyCount`] if the key count
    /// differs from the read count, plus the mapping errors documented by
    /// [`Self::map_reads_for_output`].
    #[allow(clippy::too_many_arguments)]
    pub fn map_reads_for_output_with_tie_break_keys<'a>(
        &'a mut self,
        reference: &ReferenceIndex,
        reads: &[&[Base]],
        maximum_edit_distance: u8,
        library_profile: LibraryProfile,
        search_mode: SingleSearchMode,
        tie_break_seed: u64,
        read_keys: &[u64],
    ) -> Result<&'a [SingleAlignmentResult], AlignmentError> {
        if reads.len() != read_keys.len() {
            return Err(AlignmentError::ReportingTieBreakKeyCount {
                reads: reads.len(),
                keys: read_keys.len(),
            });
        }
        self.map_reads_for_output_impl(
            reference,
            reads,
            maximum_edit_distance,
            library_profile,
            search_mode,
            Some(ReportingTieBreakBatch {
                seed: tie_break_seed,
                read_keys,
            }),
        )
    }

    #[allow(clippy::too_many_lines)]
    fn map_reads_for_output_impl<'a>(
        &'a mut self,
        reference: &ReferenceIndex,
        reads: &[&[Base]],
        maximum_edit_distance: u8,
        library_profile: LibraryProfile,
        search_mode: SingleSearchMode,
        tie_break: Option<ReportingTieBreakBatch<'_>>,
    ) -> Result<&'a [SingleAlignmentResult], AlignmentError> {
        self.map_reads_with_mode_impl(
            reference,
            reads,
            maximum_edit_distance,
            library_profile,
            search_mode,
            tie_break,
        )?;
        if !self.output_policy.adapter_clipping_enabled() {
            return Ok(&self.results);
        }
        let mut clipped_reads = Vec::new();
        let mut clipped_metadata = Vec::new();

        for (offset, read) in reads.iter().enumerate() {
            let Some(retained_end) = supported_three_prime_adapter_start(read, &self.output_policy)
                .filter(|&start| start >= MIN_ADAPTER_RETAINED_BASES)
            else {
                continue;
            };
            clipped_reads.push(&read[..retained_end]);
            clipped_metadata.push((offset, retained_end));
        }

        if clipped_reads.is_empty() {
            return Ok(&self.results);
        }

        let primary = self.results.clone();
        let mut adapter_results = vec![None; reads.len()];
        {
            let clipped_read_keys = tie_break.map(|batch| {
                clipped_metadata
                    .iter()
                    .map(|(offset, _)| batch.read_keys[*offset])
                    .collect::<Vec<_>>()
            });
            let clipped_tie_break =
                tie_break
                    .zip(clipped_read_keys.as_deref())
                    .map(|(batch, read_keys)| ReportingTieBreakBatch {
                        seed: batch.seed,
                        read_keys,
                    });
            let remapped = self
                .map_reads_with_mode_impl(
                    reference,
                    &clipped_reads,
                    maximum_edit_distance,
                    library_profile,
                    search_mode,
                    clipped_tie_break,
                )?
                .to_vec();
            for ((offset, retained_end), result) in clipped_metadata.iter().copied().zip(remapped) {
                adapter_results[offset] = Some(SingleAdapterFallback {
                    result,
                    stability_result: None,
                    final_status: result.status(),
                    retained_end,
                });
            }

            let mut stability_reads = Vec::with_capacity(clipped_reads.len());
            let mut stability_metadata = Vec::with_capacity(clipped_reads.len());
            for (offset, fallback) in adapter_results.iter_mut().enumerate() {
                let Some(fallback) = fallback else {
                    continue;
                };
                if !matches!(fallback.final_status, SingleMappingStatus::Unique) {
                    continue;
                }
                if fallback.retained_end
                    < MIN_ADAPTER_RETAINED_BASES.saturating_add(ADAPTER_STABILITY_DELTA)
                {
                    fallback.final_status = SingleMappingStatus::Ambiguous;
                    continue;
                }
                let stability_end = fallback.retained_end - ADAPTER_STABILITY_DELTA;
                stability_reads.push(&reads[offset][..stability_end]);
                stability_metadata.push(offset);
            }

            if !stability_reads.is_empty() {
                let stability_read_keys = tie_break.map(|batch| {
                    stability_metadata
                        .iter()
                        .map(|offset| batch.read_keys[*offset])
                        .collect::<Vec<_>>()
                });
                let stability_tie_break =
                    tie_break
                        .zip(stability_read_keys.as_deref())
                        .map(|(batch, read_keys)| ReportingTieBreakBatch {
                            seed: batch.seed,
                            read_keys,
                        });
                let stability = self
                    .map_reads_with_mode_impl(
                        reference,
                        &stability_reads,
                        maximum_edit_distance,
                        library_profile,
                        search_mode,
                        stability_tie_break,
                    )?
                    .to_vec();
                for (offset, stability_result) in stability_metadata.into_iter().zip(stability) {
                    let fallback = adapter_results[offset]
                        .as_mut()
                        .expect("stability metadata refers to an adapter result");
                    fallback.stability_result = Some(stability_result);
                    let same_origin = fallback
                        .result
                        .placement()
                        .zip(stability_result.placement())
                        .is_some_and(|(candidate, stability)| {
                            placement_origin_key(candidate, fallback.retained_end)
                                == placement_origin_key(stability, fallback.retained_end)
                        });
                    if !matches!(stability_result.status(), SingleMappingStatus::Unique)
                        || !same_origin
                    {
                        fallback.final_status = SingleMappingStatus::Ambiguous;
                    }
                }
            }
        }

        self.output_results.clear();
        for (offset, strict) in primary.into_iter().enumerate() {
            let Some(fallback) = adapter_results[offset] else {
                self.output_results.push(strict);
                continue;
            };
            if !matches!(strict.status(), SingleMappingStatus::Unmapped) {
                let mut result = strict;
                result.adapter_attempted = true;
                result.adapter_status = Some(fallback.final_status);
                result.adapter_clipped_bases =
                    reads[offset].len().saturating_sub(fallback.retained_end);
                result.located_rows = result
                    .located_rows
                    .saturating_add(fallback.result.located_rows)
                    .saturating_add(
                        fallback
                            .stability_result
                            .map_or(0, SingleAlignmentResult::located_rows),
                    );
                result.verified_placements = result
                    .verified_placements
                    .saturating_add(fallback.result.verified_placements)
                    .saturating_add(
                        fallback
                            .stability_result
                            .map_or(0, SingleAlignmentResult::verified_placements),
                    );
                let same_origin = strict
                    .placement()
                    .zip(fallback.result.placement())
                    .is_some_and(|(selected, endpoint)| {
                        placement_origin_key(selected, reads[offset].len())
                            == placement_origin_key(endpoint, fallback.retained_end)
                    });
                if matches!(fallback.final_status, SingleMappingStatus::Unique) && same_origin {
                    result.placement = fallback.result.placement;
                    result.retained_query_end = fallback.retained_end;
                }
                self.output_results.push(result);
                continue;
            }
            let mut result = fallback.result;
            result.status = fallback.final_status;
            result.adapter_attempted = true;
            result.adapter_status = Some(fallback.final_status);
            result.adapter_clipped_bases =
                reads[offset].len().saturating_sub(fallback.retained_end);
            result.located_rows = strict
                .located_rows
                .saturating_add(result.located_rows)
                .saturating_add(
                    fallback
                        .stability_result
                        .map_or(0, SingleAlignmentResult::located_rows),
                );
            result.verified_placements = strict
                .verified_placements
                .saturating_add(result.verified_placements)
                .saturating_add(
                    fallback
                        .stability_result
                        .map_or(0, SingleAlignmentResult::verified_placements),
                );
            if result.placement.is_some() {
                result.retained_query_end = fallback.retained_end;
            } else {
                result.retained_query_end = reads[offset].len();
            }
            result.mapping_quality = adapter_consensus_mapping_quality(
                matches!(fallback.final_status, SingleMappingStatus::Unique),
                result.mapping_quality,
                fallback
                    .stability_result
                    .map(SingleAlignmentResult::mapping_quality),
            );
            self.output_results.push(result);
        }
        Ok(&self.output_results)
    }

    #[allow(clippy::too_many_lines)]
    fn prepare_search_wavefront(
        &mut self,
        reference: &ReferenceIndex,
        reads: &[&[Base]],
        completion_limits: CombinedSearchLimits,
        conversion_pass: ConversionPass,
        complete_candidate_frontier: bool,
    ) -> Result<(), AlignmentError> {
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

        let mut active_ordinals = [0_usize; SINGLE_ALIGNMENT_BATCH_SIZE];
        let mut available_bases = [0_usize; SINGLE_ALIGNMENT_BATCH_SIZE];
        let mut compact_matches = [None; SINGLE_ALIGNMENT_BATCH_SIZE];
        let mut locatable_ordinals = [0_usize; SINGLE_ALIGNMENT_BATCH_SIZE];

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
                let mut patterns = [&self.projections[first][..available_bases[first]];
                    SINGLE_ALIGNMENT_BATCH_SIZE];
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

            if INITIAL_SEARCH_LIMITS.maximum_seed_rounds < completion_limits.maximum_seed_rounds {
                for &ordinal in &active_ordinals[..active_count] {
                    if let Some(seed) = self.round_matches[ordinal]
                        && !combined_seed_round_is_locatable(seed, INITIAL_SEARCH_LIMITS)
                        && combined_seed_round_is_locatable(seed, completion_limits)
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
                                let Some(strand) =
                                    conversion_pass.relabel_combined_hit(hit.strand())
                                else {
                                    return;
                                };
                                let mut candidate = ReadCandidate {
                                    contig_ordinal: hit.contig_ordinal(),
                                    start: hit.start(),
                                    strand,
                                    proof_mask: FLEXIBLE_NOMINAL_PROOF | (1_u16 << round),
                                };
                                if round == 0
                                    && matches[lane].exact_hit_count() == 1
                                    && let Some(distance) =
                                        ungapped_distance(reference, reads[ordinal], candidate)
                                {
                                    candidate.proof_mask = direct_singleton_proof(distance);
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
                        self.search_states[ordinal].active[0] &=
                            !direct[lane] || complete_candidate_frontier;
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
                                let Some(strand) =
                                    conversion_pass.relabel_combined_hit(hit.strand())
                                else {
                                    return true;
                                };
                                let mut candidate = ReadCandidate {
                                    contig_ordinal: hit.contig_ordinal(),
                                    start: hit.start(),
                                    strand,
                                    proof_mask: FLEXIBLE_NOMINAL_PROOF | (1_u16 << round),
                                };
                                if round == 0
                                    && matches.exact_hit_count() == 1
                                    && let Some(distance) =
                                        ungapped_distance(reference, reads[ordinal], candidate)
                                {
                                    candidate.proof_mask = direct_singleton_proof(distance);
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
                    self.search_states[ordinal].active[0] &= !direct || complete_candidate_frontier;
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

    #[allow(clippy::too_many_lines)]
    fn map_one(
        &mut self,
        reference: &ReferenceIndex,
        prepared: PreparedSingleRead<'_>,
        maximum_edit_distance: u8,
        search_mode: SingleSearchMode,
        conversion_pass: ConversionPass,
    ) -> Result<SingleAlignmentResult, AlignmentError> {
        let PreparedSingleRead {
            read,
            projection,
            first_seed,
            search,
            initial_candidates,
        } = prepared;
        let mut search = *search;
        let completion_search_start = search;
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
        let mut verified_distance_limit = initial_budget;
        if read_workspace.placements.is_empty() && maximum_edit_distance > initial_budget {
            read_workspace.candidates.clear();
            read_workspace.placements.clear();
            let (_, observed) = read_workspace.verify_candidates_with_budget(
                reference,
                read,
                metrics,
                maximum_edit_distance,
            )?;
            metrics = observed;
            verified_distance_limit = maximum_edit_distance;
        }

        if read_workspace.placements.is_empty() {
            read_workspace.candidates.clear();
            read_workspace.placements.clear();
            unused_workspace.candidate_nominals.clear();
            let additional = continue_combined_two_lane_search(
                reference,
                [read, &[]],
                [projection, &[]],
                [conversion_pass, ConversionPass::Original],
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
            metrics = observed;
            verified_distance_limit = maximum_edit_distance;
        }

        let incumbent = Self::finish_result_with_affine_rerank(
            reference,
            read,
            read_workspace,
            origins,
            SingleResultEvidence {
                read_length: read.len(),
                metrics,
                verified_distance_limit,
                first_seed,
                frontier_complete: false,
                confidence_audit_complete: false,
            },
        )?;
        if matches!(search_mode, SingleSearchMode::Default) {
            if !Self::default_high_confidence_cross_check_required(incumbent) {
                return Ok(incumbent);
            }
            let mut cross_check_search = completion_search_start;
            read_workspace.candidates.clear();
            unused_workspace.candidate_nominals.clear();
            continue_combined_two_lane_search_with_limits(
                reference,
                [read, &[]],
                [projection, &[]],
                [conversion_pass, ConversionPass::Original],
                DEFAULT_HIGH_CONFIDENCE_AUDIT_SEARCH_LIMITS,
                true,
                &mut cross_check_search,
                &mut read_workspace.candidate_nominals,
                &mut unused_workspace.candidate_nominals,
            )?;
            metrics.located_rows = cross_check_search.located[0];
            let (_, observed) = read_workspace.verify_candidates_with_budget(
                reference,
                read,
                metrics,
                maximum_edit_distance,
            )?;
            return Self::finish_result_with_affine_rerank(
                reference,
                read,
                read_workspace,
                origins,
                SingleResultEvidence {
                    read_length: read.len(),
                    metrics: observed,
                    verified_distance_limit: maximum_edit_distance,
                    first_seed,
                    frontier_complete: false,
                    confidence_audit_complete: true,
                },
            );
        }
        if !Self::sensitive_audit_required(incumbent) {
            return Ok(incumbent);
        }

        let mut completion_search = completion_search_start;
        read_workspace.candidates.clear();
        unused_workspace.candidate_nominals.clear();
        continue_combined_two_lane_search_with_limits(
            reference,
            [read, &[]],
            [projection, &[]],
            [conversion_pass, ConversionPass::Original],
            search_mode.limits(),
            true,
            &mut completion_search,
            &mut read_workspace.candidate_nominals,
            &mut unused_workspace.candidate_nominals,
        )?;
        let audit_distance_limit =
            Self::sensitive_audit_distance_limit(incumbent, maximum_edit_distance);
        metrics.located_rows = completion_search.located[0];
        let (_, observed) = read_workspace.verify_candidates_with_budget(
            reference,
            read,
            metrics,
            audit_distance_limit,
        )?;
        let completed = Self::finish_result_with_affine_rerank(
            reference,
            read,
            read_workspace,
            origins,
            SingleResultEvidence {
                read_length: read.len(),
                metrics: observed,
                verified_distance_limit: audit_distance_limit,
                first_seed,
                frontier_complete: true,
                confidence_audit_complete: true,
            },
        )?;
        Ok(Self::reconcile_sensitive_result(
            incumbent,
            completed,
            read.len(),
        ))
    }

    #[allow(clippy::too_many_lines)]
    fn map_two(
        &mut self,
        reference: &ReferenceIndex,
        prepared: [PreparedSingleRead<'_>; 2],
        maximum_edit_distance: u8,
        search_mode: SingleSearchMode,
        conversion_pass: ConversionPass,
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
        let completion_search_start = search;
        let mut metrics = search.located.map(|located_rows| ReadAlignmentMetrics {
            located_rows,
            ..ReadAlignmentMetrics::default()
        });
        let initial_budget = maximum_edit_distance.min(INITIAL_EDIT_DISTANCE);
        let mut results = [None; 2];
        let mut verified_distance_limits = [initial_budget; 2];

        for lane in 0..2 {
            let (_, observed) = workspaces[lane].verify_candidates_with_budget(
                reference,
                reads[lane],
                metrics[lane],
                initial_budget,
            )?;
            metrics[lane] = observed;
            if !workspaces[lane].placements.is_empty() {
                results[lane] = Some(Self::finish_result_with_affine_rerank(
                    reference,
                    reads[lane],
                    &mut workspaces[lane],
                    origins,
                    SingleResultEvidence {
                        read_length: reads[lane].len(),
                        metrics: metrics[lane],
                        verified_distance_limit: initial_budget,
                        first_seed: first_seeds[lane],
                        frontier_complete: false,
                        confidence_audit_complete: false,
                    },
                )?);
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
                verified_distance_limits[lane] = maximum_edit_distance;
                if !workspaces[lane].placements.is_empty() {
                    results[lane] = Some(Self::finish_result_with_affine_rerank(
                        reference,
                        reads[lane],
                        &mut workspaces[lane],
                        origins,
                        SingleResultEvidence {
                            read_length: reads[lane].len(),
                            metrics: metrics[lane],
                            verified_distance_limit: maximum_edit_distance,
                            first_seed: first_seeds[lane],
                            frontier_complete: false,
                            confidence_audit_complete: false,
                        },
                    )?);
                }
            }
        }

        if results.iter().any(Option::is_none) {
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
                [conversion_pass; 2],
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
                results[lane] = Some(Self::finish_result_with_affine_rerank(
                    reference,
                    reads[lane],
                    &mut workspaces[lane],
                    origins,
                    SingleResultEvidence {
                        read_length: reads[lane].len(),
                        metrics: observed,
                        verified_distance_limit: maximum_edit_distance,
                        first_seed: first_seeds[lane],
                        frontier_complete: false,
                        confidence_audit_complete: false,
                    },
                )?);
            }
        }
        let incumbents = results
            .map(|result| result.expect("default single-read continuation resolves every lane"));
        if matches!(search_mode, SingleSearchMode::Default) {
            let cross_check_required =
                incumbents.map(Self::default_high_confidence_cross_check_required);
            if !cross_check_required[0] && !cross_check_required[1] {
                return Ok(incumbents);
            }
            let mut cross_check_search = completion_search_start;
            for (lane, &required) in cross_check_required.iter().enumerate() {
                if !required {
                    cross_check_search.active[lane] = false;
                    cross_check_search.direct[lane] = false;
                    cross_check_search.deferred_len[lane] = 0;
                }
            }
            let [first_workspace, second_workspace] = workspaces;
            continue_combined_two_lane_search_with_limits(
                reference,
                reads,
                projections,
                [conversion_pass; 2],
                DEFAULT_HIGH_CONFIDENCE_AUDIT_SEARCH_LIMITS,
                true,
                &mut cross_check_search,
                &mut first_workspace.candidate_nominals,
                &mut second_workspace.candidate_nominals,
            )?;
            let mut completed = incumbents;
            for (lane, &required) in cross_check_required.iter().enumerate() {
                if !required {
                    continue;
                }
                metrics[lane].located_rows = cross_check_search.located[lane];
                let (_, observed) = workspaces[lane].verify_candidates_with_budget(
                    reference,
                    reads[lane],
                    metrics[lane],
                    maximum_edit_distance,
                )?;
                completed[lane] = Self::finish_result_with_affine_rerank(
                    reference,
                    reads[lane],
                    &mut workspaces[lane],
                    origins,
                    SingleResultEvidence {
                        read_length: reads[lane].len(),
                        metrics: observed,
                        verified_distance_limit: maximum_edit_distance,
                        first_seed: first_seeds[lane],
                        frontier_complete: false,
                        confidence_audit_complete: true,
                    },
                )?;
            }
            return Ok(completed);
        }

        let audit_required = incumbents.map(Self::sensitive_audit_required);
        if !audit_required[0] && !audit_required[1] {
            return Ok(incumbents);
        }
        let mut completion_search = completion_search_start;
        for (lane, &required) in audit_required.iter().enumerate() {
            if !required {
                completion_search.active[lane] = false;
                completion_search.direct[lane] = false;
                completion_search.deferred_len[lane] = 0;
            }
        }
        let [first_workspace, second_workspace] = workspaces;
        continue_combined_two_lane_search_with_limits(
            reference,
            reads,
            projections,
            [conversion_pass; 2],
            search_mode.limits(),
            true,
            &mut completion_search,
            &mut first_workspace.candidate_nominals,
            &mut second_workspace.candidate_nominals,
        )?;
        let mut completed = [SingleAlignmentResult::unmapped(0, 0, 0); 2];
        for (lane, &required) in audit_required.iter().enumerate() {
            if !required {
                completed[lane] = incumbents[lane];
                continue;
            }
            let audit_distance_limit =
                Self::sensitive_audit_distance_limit(incumbents[lane], maximum_edit_distance);
            metrics[lane].located_rows = completion_search.located[lane];
            let (_, observed) = workspaces[lane].verify_candidates_with_budget(
                reference,
                reads[lane],
                metrics[lane],
                audit_distance_limit,
            )?;
            completed[lane] = Self::finish_result_with_affine_rerank(
                reference,
                reads[lane],
                &mut workspaces[lane],
                origins,
                SingleResultEvidence {
                    read_length: reads[lane].len(),
                    metrics: observed,
                    verified_distance_limit: audit_distance_limit
                        .max(verified_distance_limits[lane]),
                    first_seed: first_seeds[lane],
                    frontier_complete: true,
                    confidence_audit_complete: true,
                },
            )?;
        }
        Ok(core::array::from_fn(|lane| {
            if audit_required[lane] {
                Self::reconcile_sensitive_result(
                    incumbents[lane],
                    completed[lane],
                    reads[lane].len(),
                )
            } else {
                incumbents[lane]
            }
        }))
    }

    const fn default_high_confidence_cross_check_required(
        incumbent: SingleAlignmentResult,
    ) -> bool {
        if !matches!(incumbent.status, SingleMappingStatus::Unique) {
            return false;
        }
        match incumbent.placement {
            Some(placement) => default_confidence_cross_check_required(
                incumbent.mapping_quality,
                placement.distance(),
            ),
            None => false,
        }
    }

    const fn sensitive_audit_required(_incumbent: SingleAlignmentResult) -> bool {
        // Sensitive mode is a complete bounded-search contract, not merely a
        // higher-confidence rescore of provisional unique mappings. Complete
        // ambiguous and initially unmapped frontiers as well so every
        // conversion pass is compared against the same search boundary.
        true
    }

    fn reconcile_sensitive_result(
        incumbent: SingleAlignmentResult,
        completed: SingleAlignmentResult,
        read_length: usize,
    ) -> SingleAlignmentResult {
        let Some(completed_placement) = completed.placement else {
            return incumbent
                .placement
                .map_or(completed, |placement| SingleAlignmentResult {
                    status: incumbent.status,
                    placement: Some(placement),
                    retained_query_end: read_length,
                    mapping_quality: incumbent.mapping_quality,
                    located_rows: completed.located_rows,
                    distinct_candidate_starts: completed.distinct_candidate_starts,
                    verified_placements: completed.verified_placements,
                    best_origin_count: completed.best_origin_count,
                    adapter_attempted: false,
                    adapter_status: None,
                    adapter_clipped_bases: 0,
                });
        };
        let Some(incumbent_placement) = incumbent.placement else {
            // Preserve a completed ambiguous representative at MAPQ 0. SAM
            // can report the best known coordinate without claiming that the
            // read is uniquely placed.
            return if matches!(completed.status, SingleMappingStatus::Ambiguous)
                || sensitive_replacement_certified(completed.mapping_quality)
            {
                completed
            } else {
                SingleAlignmentResult::unmapped_with_evidence(
                    read_length,
                    completed.located_rows,
                    completed.distinct_candidate_starts,
                    completed.verified_placements,
                )
            };
        };
        if placement_origin_key(incumbent_placement, read_length)
            == placement_origin_key(completed_placement, read_length)
        {
            return completed;
        }
        if matches!(completed.status, SingleMappingStatus::Unique)
            && sensitive_replacement_certified(completed.mapping_quality)
        {
            return completed;
        }
        SingleAlignmentResult {
            status: SingleMappingStatus::Ambiguous,
            placement: Some(incumbent_placement),
            retained_query_end: read_length,
            mapping_quality: 0,
            located_rows: completed.located_rows,
            distinct_candidate_starts: completed.distinct_candidate_starts,
            verified_placements: completed.verified_placements,
            best_origin_count: completed.best_origin_count,
            adapter_attempted: false,
            adapter_status: None,
            adapter_clipped_bases: 0,
        }
    }

    fn sensitive_audit_distance_limit(
        _incumbent: SingleAlignmentResult,
        maximum_edit_distance: u8,
    ) -> u8 {
        // Audit the complete configured verification radius so a farther
        // runner-up contributes to the same score-gap model as every other
        // verified competitor.
        maximum_edit_distance
    }

    #[allow(clippy::too_many_lines)]
    fn finish_result(
        workspace: &ReadWorkspace,
        origins: &mut Vec<(u64, BisulfiteStrand, i128)>,
        evidence: SingleResultEvidence,
    ) -> SingleAlignmentResult {
        let SingleResultEvidence {
            read_length,
            metrics,
            verified_distance_limit,
            first_seed,
            frontier_complete,
            confidence_audit_complete,
        } = evidence;
        let Some(best_distance) = workspace
            .placements
            .iter()
            .map(|value| value.distance())
            .min()
        else {
            return SingleAlignmentResult::unmapped_with_evidence(
                read_length,
                metrics.located_rows,
                metrics.distinct_candidate_starts,
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
        // The canonical representative is already optimal when it spans the
        // complete retained query. Pay the endpoint-parsimony scan only for
        // the sparse residual whose canonical endpoint contains a net gap.
        if best_distance != 0
            && origins.len() > 1
            && representative
                .is_some_and(|placement| placement_net_gap_bases(placement, read_length) != 0)
        {
            representative = prefer_minimum_net_gap_representative(
                &workspace.placements,
                best_distance,
                read_length,
            );
        }
        origins.sort_unstable();
        origins.dedup();
        let exact_origin_count = origins.len();
        let best_locus_count = if frontier_complete {
            local_origin_locus_count(origins, best_distance)
        } else {
            exact_origin_count
        };
        let local_origin_collapsed = exact_origin_count > 1 && best_locus_count == 1;
        let status = if best_locus_count == 1 {
            SingleMappingStatus::Unique
        } else {
            SingleMappingStatus::Ambiguous
        };
        let mapping_quality = if matches!(status, SingleMappingStatus::Unique) {
            let best_origin = representative.map_or(origins[0], |placement| {
                placement_origin_key(placement, read_length)
            });
            let second_best_distance = workspace
                .placements
                .iter()
                .copied()
                .filter(|placement| {
                    let origin = placement_origin_key(*placement, read_length);
                    if local_origin_collapsed {
                        !origins_share_local_locus(origin, best_origin, best_distance)
                    } else {
                        origin != best_origin
                    }
                })
                .map(ReadPlacement::distance)
                .min();
            let first_seed_hits = first_seed.map_or(0, CombinedSeedMatches::exact_hit_count);
            let first_seed_bases = first_seed.map_or(0, CombinedSeedMatches::matched_bases);
            // A one-row first seed certifies only the coordinate it emitted.
            // Sensitive completion may select a different coordinate, so the
            // direct flag must come from the retained placement's own proof.
            let retained_proof =
                representative.map_or(0, |placement| placement_proof_mask(workspace, placement));
            let direct_singleton = retained_proof & DIRECT_SINGLETON_PROOF != 0;
            let seed_round_support = if retained_proof & FLEXIBLE_NOMINAL_PROOF == 0 {
                0
            } else {
                seed_round_support(retained_proof)
            };
            let mapq_evidence = SingleMapqEvidence {
                read_length,
                best_distance,
                second_best_distance,
                verified_distance_limit,
                located_rows: metrics.located_rows,
                distinct_candidate_starts: metrics.distinct_candidate_starts,
                verified_placements: metrics.verified_placements,
                first_seed_hits,
                first_seed_bases,
                direct_singleton,
                frontier_complete,
                confidence_audit_complete,
                seed_round_support,
            };
            local_locus_mapping_quality(mapq_evidence, local_origin_collapsed)
        } else {
            0
        };
        SingleAlignmentResult {
            status,
            placement: representative,
            retained_query_end: read_length,
            mapping_quality,
            located_rows: metrics.located_rows,
            distinct_candidate_starts: metrics.distinct_candidate_starts,
            verified_placements: metrics.verified_placements,
            best_origin_count: u64::try_from(best_locus_count).unwrap_or(u64::MAX),
            adapter_attempted: false,
            adapter_status: None,
            adapter_clipped_bases: 0,
        }
    }

    // Whole-read edit distance remains the primary objective.  Within a
    // bounded equal-distance origin set, affine score is the secondary
    // objective. A unique affine winner may contribute only the policy's
    // conservative affine certificate, and only after the sensitive frontier
    // completed.
    fn finish_result_with_affine_rerank(
        reference: &ReferenceIndex,
        read: &[Base],
        workspace: &mut ReadWorkspace,
        origins: &mut Vec<(u64, BisulfiteStrand, i128)>,
        evidence: SingleResultEvidence,
    ) -> Result<SingleAlignmentResult, AlignmentError> {
        workspace.affine_scores.clear();
        let frontier_complete = evidence.frontier_complete;
        let mut result = Self::finish_result(workspace, origins, evidence);
        if !matches!(result.status, SingleMappingStatus::Ambiguous)
            || !affine_rerank_origin_count_supported(
                u64::try_from(origins.len()).unwrap_or(u64::MAX),
            )
        {
            return Ok(result);
        }
        let Some(current) = result.placement else {
            return Ok(result);
        };
        let best_distance = current.distance();
        let mut affine_workspace = AffineScoreWorkspace::default();
        let mut ranked_origins = Vec::with_capacity(origins.len());
        for placement in workspace
            .placements
            .iter()
            .copied()
            .filter(|placement| placement.distance() == best_distance)
        {
            let origin = placement_origin_key(placement, read.len());
            let score = if let Some(score) = workspace
                .affine_score_cache
                .iter()
                .find_map(|(cached, score)| (*cached == placement).then_some(*score))
            {
                score
            } else {
                let score =
                    affine_placement_score(reference, read, placement, 0, &mut affine_workspace)?;
                workspace.affine_score_cache.push((placement, score));
                score
            };
            workspace.affine_scores.push((placement, score));
            if let Some((_, retained, retained_score)) = ranked_origins
                .iter_mut()
                .find(|(retained_origin, _, _)| *retained_origin == origin)
            {
                if score > *retained_score || score == *retained_score && placement < *retained {
                    *retained = placement;
                    *retained_score = score;
                }
            } else {
                ranked_origins.push((origin, placement, score));
            }
        }
        ranked_origins.sort_unstable_by_key(|(_, placement, score)| (Reverse(*score), *placement));
        let Some((_, representative, best_score)) = ranked_origins.first().copied() else {
            return Ok(result);
        };
        result.placement = Some(representative);
        if let Some(mapping_quality) = affine_unique_mapping_quality(
            frontier_complete,
            best_score,
            ranked_origins.get(1).map(|(_, _, score)| *score),
        ) {
            result.status = SingleMappingStatus::Unique;
            result.mapping_quality = mapping_quality;
            result.best_origin_count = 1;
        }
        Ok(result)
    }
}

#[cfg(test)]
#[path = "../../tests/whitebox/single_end.rs"]
mod whitebox;
