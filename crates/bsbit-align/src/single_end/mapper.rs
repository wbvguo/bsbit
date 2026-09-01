//! High-throughput single-end orchestration over the shared combined-index search core.

use bsbit_core::alphabet::Base;
use bsbit_core::bisulfite::BisulfiteStrand;
use bsbit_index::reference::ReferenceIndex;
use bsbit_index::storage::fm::ProjectedBase;

use super::mapq::{SingleMapqEvidence, single_mapping_quality_from_evidence};
use crate::AlignmentError;
use crate::adapter::{
    ADAPTER_STABILITY_DELTA, MIN_ADAPTER_RETAINED_BASES, supported_three_prime_adapter_start,
};
use crate::library::{ConversionPass, LibraryProfile};
use crate::placement::{ReadPlacement, placement_origin_key};
use crate::read_mapping_limits::{
    INITIAL_EDIT_DISTANCE, MAX_EDIT_DISTANCE, MAX_READ_BASES, MIN_SUFFIX_BASES,
};
use crate::search::combined_query::{CombinedSearchReferenceExt, CombinedSeedMatches};

use crate::read_mapping::{ReadAlignmentMetrics, ReadCandidate, ReadWorkspace, ungapped_distance};
use crate::search::combined_adaptive::{
    CombinedSearchLimits, CombinedTwoLaneSearchState, DEFAULT_SEARCH_LIMITS,
    DIRECT_SINGLETON_PROOF, DeferredCombinedSeed, EMPTY_SEED_STEP, FLEXIBLE_NOMINAL_PROOF,
    INITIAL_SEARCH_LIMITS, SENSITIVE_SEARCH_LIMITS, combined_seed_round_is_locatable,
    continue_combined_two_lane_search, continue_combined_two_lane_search_with_limits,
    prepare_combined_projection,
};

const SINGLE_ADAPTER_MAX_MAPQ: u8 = 20;

/// Final classification for one single read.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SingleMappingStatus {
    /// No verified placement survived the bounded search.
    Unmapped,
    /// Exactly one best biological origin survived.
    Unique,
    /// Multiple plausible biological origins survived the confidence policy.
    Ambiguous,
}

const SENSITIVE_REPLACEMENT_MIN_MAPQ: u8 = 20;

/// Candidate-search effort for single-end alignment.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum SingleSearchMode {
    /// Qualified low-latency alignment with an incremental fallback.
    #[default]
    Default,
    /// Default mapping followed by a bounded confidence audit.
    Sensitive,
}

impl SingleSearchMode {
    const fn limits(self) -> CombinedSearchLimits {
        match self {
            Self::Default => DEFAULT_SEARCH_LIMITS,
            Self::Sensitive => SENSITIVE_SEARCH_LIMITS,
        }
    }

    const fn completes_candidate_frontier(self) -> bool {
        matches!(self, Self::Sensitive)
    }
}

/// Final mapping facts for one single read.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SingleAlignmentResult {
    status: SingleMappingStatus,
    placement: Option<ReadPlacement>,
    retained_query_end: usize,
    mapping_quality: u8,
    located_rows: u64,
    distinct_candidate_starts: u64,
    verified_placements: u64,
    adapter_attempted: bool,
    adapter_status: Option<SingleMappingStatus>,
    adapter_clipped_bases: usize,
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

    /// Returns the retained sequencing-orientation query interval.
    #[must_use]
    pub const fn retained_query_interval(self) -> core::ops::Range<usize> {
        0..self.retained_query_end
    }

    /// Returns the calibrated SAM mapping quality, or zero when not unique.
    #[must_use]
    pub const fn mapping_quality(self) -> u8 {
        self.mapping_quality
    }

    /// Returns suffix rows located across every executed mapping phase.
    #[must_use]
    pub const fn located_rows(self) -> u64 {
        self.located_rows
    }

    /// Returns verified placements across every executed mapping phase before
    /// per-phase best-tier selection.
    #[must_use]
    pub const fn verified_placements(self) -> u64 {
        self.verified_placements
    }

    /// Reports whether exact adapter support triggered a trimmed remap.
    #[must_use]
    pub const fn adapter_attempted(self) -> bool {
        self.adapter_attempted
    }

    /// Returns the adapter-remap class after stability verification.
    #[must_use]
    pub const fn adapter_status(self) -> Option<SingleMappingStatus> {
        self.adapter_status
    }

    /// Returns the number of bases omitted at the supported 3' adapter boundary.
    #[must_use]
    pub const fn adapter_clipped_bases(self) -> usize {
        self.adapter_clipped_bases
    }

    const fn unmapped(read_length: usize, located_rows: u64, verified_placements: u64) -> Self {
        Self {
            status: SingleMappingStatus::Unmapped,
            placement: None,
            retained_query_end: read_length,
            mapping_quality: 0,
            located_rows,
            distinct_candidate_starts: 0,
            verified_placements,
            adapter_attempted: false,
            adapter_status: None,
            adapter_clipped_bases: 0,
        }
    }

    const fn unmapped_with_evidence(
        read_length: usize,
        located_rows: u64,
        distinct_candidate_starts: u64,
        verified_placements: u64,
    ) -> Self {
        Self {
            status: SingleMappingStatus::Unmapped,
            placement: None,
            retained_query_end: read_length,
            mapping_quality: 0,
            located_rows,
            distinct_candidate_starts,
            verified_placements,
            adapter_attempted: false,
            adapter_status: None,
            adapter_clipped_bases: 0,
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

#[derive(Clone, Copy)]
struct SingleAdapterFallback {
    result: SingleAlignmentResult,
    stability_result: Option<SingleAlignmentResult>,
    final_status: SingleMappingStatus,
    retained_end: usize,
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
        self.map_reads_with_mode(
            reference,
            reads,
            maximum_edit_distance,
            SingleSearchMode::Default,
        )
    }

    /// Maps a batch with the selected single-end candidate-search policy.
    ///
    /// Sensitive mode first obtains the default result and then completes the
    /// bounded seed frontier as a confidence audit. Q20 incumbents at edit
    /// distance two or better need verification only through distance three;
    /// weaker results retain the full distance-five verification boundary. A
    /// different-origin result or a new rescue must be unique at Q20 or above;
    /// lower-confidence conflicts retain the default representative at Q0. It
    /// does not apply paired-only mate rescue, template geometry, or
    /// adapter-pair recovery. This search-only entry point does not apply the
    /// single-end output adapter policy; callers producing complete directional
    /// records should use [`Self::map_reads_for_output`].
    ///
    /// # Errors
    ///
    /// Returns an unsupported read/edit domain or combined-index failure.
    #[allow(clippy::too_many_lines)]
    pub fn map_reads_with_mode<'a>(
        &'a mut self,
        reference: &ReferenceIndex,
        reads: &[&[Base]],
        maximum_edit_distance: u8,
        search_mode: SingleSearchMode,
    ) -> Result<&'a [SingleAlignmentResult], AlignmentError> {
        self.map_reads_with_profile_and_mode(
            reference,
            reads,
            maximum_edit_distance,
            LibraryProfile::Directional,
            search_mode,
        )
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
    pub fn map_reads_with_profile_and_mode<'a>(
        &'a mut self,
        reference: &ReferenceIndex,
        reads: &[&[Base]],
        maximum_edit_distance: u8,
        library_profile: LibraryProfile,
        search_mode: SingleSearchMode,
    ) -> Result<&'a [SingleAlignmentResult], AlignmentError> {
        let passes = library_profile.conversion_passes();
        let first = passes[0];
        self.map_reads_pass(reference, reads, maximum_edit_distance, search_mode, first)?;
        if passes.len() == 1 {
            return Ok(&self.results);
        }

        self.primary_pass_results.clear();
        self.primary_pass_results.extend_from_slice(&self.results);
        self.map_reads_pass(
            reference,
            reads,
            maximum_edit_distance,
            search_mode,
            passes[1],
        )?;
        for (complementary, original) in self.results.iter_mut().zip(&self.primary_pass_results) {
            *complementary = merge_non_directional_results(*original, *complementary);
        }
        Ok(&self.results)
    }

    /// Maps each read against OT, OB, CTOT, and CTOB, then makes one global
    /// single-end placement decision across both directional passes.
    ///
    /// Equal-best evidence in the original and complementary passes is
    /// ambiguous. A weaker cross-pass placement contributes to MAPQ separation
    /// and repeat pressure for the selected result.
    ///
    /// # Errors
    ///
    /// Returns an unsupported read/edit domain or combined-index failure.
    pub fn map_reads_non_directional_with_mode<'a>(
        &'a mut self,
        reference: &ReferenceIndex,
        reads: &[&[Base]],
        maximum_edit_distance: u8,
        search_mode: SingleSearchMode,
    ) -> Result<&'a [SingleAlignmentResult], AlignmentError> {
        self.map_reads_with_profile_and_mode(
            reference,
            reads,
            maximum_edit_distance,
            LibraryProfile::NonDirectional,
            search_mode,
        )
    }

    #[allow(clippy::too_many_lines)]
    fn map_reads_pass<'a>(
        &'a mut self,
        reference: &ReferenceIndex,
        reads: &[&[Base]],
        maximum_edit_distance: u8,
        search_mode: SingleSearchMode,
        conversion_pass: ConversionPass,
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
                    search_mode,
                    conversion_pass,
                );
                self.initial_candidates[ordinal] = initial_candidates;
                let result = result?;
                self.results.push(result);
                ordinal += 1;
            }
        }
        Ok(&self.results)
    }

    /// Maps directional complete reads and applies qualified single-end 3'
    /// adapter recovery.
    ///
    /// This compatibility entry point preserves the original directional API.
    /// New orchestration that already owns a library profile should use
    /// [`Self::map_reads_for_output_with_profile`].
    ///
    /// # Errors
    ///
    /// Returns an unsupported read/edit domain or combined-index failure from
    /// any primary, trimmed, or stability mapping phase.
    pub fn map_reads_for_output<'a>(
        &'a mut self,
        reference: &ReferenceIndex,
        reads: &[&[Base]],
        maximum_edit_distance: u8,
        search_mode: SingleSearchMode,
    ) -> Result<&'a [SingleAlignmentResult], AlignmentError> {
        self.map_reads_for_output_with_profile(
            reference,
            reads,
            maximum_edit_distance,
            LibraryProfile::Directional,
            search_mode,
        )
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
    /// recovery remains explicitly directional until a separate
    /// non-directional endpoint policy is qualified; non-directional mapping
    /// still enters through this output-policy boundary rather than bypassing
    /// it in the CLI.
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
    pub fn map_reads_for_output_with_profile<'a>(
        &'a mut self,
        reference: &ReferenceIndex,
        reads: &[&[Base]],
        maximum_edit_distance: u8,
        library_profile: LibraryProfile,
        search_mode: SingleSearchMode,
    ) -> Result<&'a [SingleAlignmentResult], AlignmentError> {
        self.map_reads_with_profile_and_mode(
            reference,
            reads,
            maximum_edit_distance,
            library_profile,
            search_mode,
        )?;
        if matches!(library_profile, LibraryProfile::NonDirectional) {
            return Ok(&self.results);
        }
        let mut clipped_reads = Vec::new();
        let mut clipped_metadata = Vec::new();

        for (offset, read) in reads.iter().enumerate() {
            let Some(retained_end) = supported_three_prime_adapter_start(read)
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
            let remapped = self
                .map_reads_with_profile_and_mode(
                    reference,
                    &clipped_reads,
                    maximum_edit_distance,
                    library_profile,
                    search_mode,
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
                let stability = self
                    .map_reads_with_profile_and_mode(
                        reference,
                        &stability_reads,
                        maximum_edit_distance,
                        library_profile,
                        search_mode,
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
            result.mapping_quality = if matches!(fallback.final_status, SingleMappingStatus::Unique)
            {
                result
                    .mapping_quality
                    .min(
                        fallback
                            .stability_result
                            .map_or(0, SingleAlignmentResult::mapping_quality),
                    )
                    .min(SINGLE_ADAPTER_MAX_MAPQ)
            } else {
                0
            };
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

        let incumbent = Self::finish_result(
            read_workspace,
            origins,
            SingleResultEvidence {
                read_length: read.len(),
                metrics,
                verified_distance_limit,
                first_seed,
                search: &search,
                lane: 0,
            },
        );
        if !search_mode.completes_candidate_frontier() {
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
        let completed = Self::finish_result(
            read_workspace,
            origins,
            SingleResultEvidence {
                read_length: read.len(),
                metrics: observed,
                verified_distance_limit: audit_distance_limit,
                first_seed,
                search: &completion_search,
                lane: 0,
            },
        );
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
                verified_distance_limits[lane] = maximum_edit_distance;
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
        }
        if !search_mode.completes_candidate_frontier() {
            return Ok(
                results.map(|result| result.expect("single-read continuation resolves every lane"))
            );
        }

        let incumbents = results
            .map(|result| result.expect("default single-read continuation resolves every lane"));
        let mut completion_search = completion_search_start;
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
        for lane in 0..2 {
            let audit_distance_limit =
                Self::sensitive_audit_distance_limit(incumbents[lane], maximum_edit_distance);
            metrics[lane].located_rows = completion_search.located[lane];
            let (_, observed) = workspaces[lane].verify_candidates_with_budget(
                reference,
                reads[lane],
                metrics[lane],
                audit_distance_limit,
            )?;
            completed[lane] = Self::finish_result(
                &workspaces[lane],
                origins,
                SingleResultEvidence {
                    read_length: reads[lane].len(),
                    metrics: observed,
                    verified_distance_limit: audit_distance_limit
                        .max(verified_distance_limits[lane]),
                    first_seed: first_seeds[lane],
                    search: &completion_search,
                    lane,
                },
            );
        }
        Ok(core::array::from_fn(|lane| {
            Self::reconcile_sensitive_result(incumbents[lane], completed[lane], reads[lane].len())
        }))
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
                    adapter_attempted: false,
                    adapter_status: None,
                    adapter_clipped_bases: 0,
                });
        };
        let Some(incumbent_placement) = incumbent.placement else {
            return if matches!(completed.status, SingleMappingStatus::Unique)
                && completed.mapping_quality >= SENSITIVE_REPLACEMENT_MIN_MAPQ
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
            && completed.mapping_quality >= SENSITIVE_REPLACEMENT_MIN_MAPQ
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
            adapter_attempted: false,
            adapter_status: None,
            adapter_clipped_bases: 0,
        }
    }

    fn sensitive_audit_distance_limit(
        incumbent: SingleAlignmentResult,
        maximum_edit_distance: u8,
    ) -> u8 {
        // A Q20 incumbent at distance two or better can only be tied, beaten,
        // or lose its Q20 separation to another origin through distance three.
        // Keep the full distance-five audit for weaker and distance-three
        // incumbents, where a distance-four runner-up can still change MAPQ.
        if incumbent
            .placement
            .is_some_and(|placement| placement.distance() <= 2)
            && incumbent.mapping_quality >= SENSITIVE_REPLACEMENT_MIN_MAPQ
        {
            maximum_edit_distance.min(INITIAL_EDIT_DISTANCE)
        } else {
            maximum_edit_distance
        }
    }

    fn finish_result(
        workspace: &ReadWorkspace,
        origins: &mut Vec<(u64, BisulfiteStrand, i128)>,
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
            retained_query_end: read_length,
            mapping_quality,
            located_rows: metrics.located_rows,
            distinct_candidate_starts: metrics.distinct_candidate_starts,
            verified_placements: metrics.verified_placements,
            adapter_attempted: false,
            adapter_status: None,
            adapter_clipped_bases: 0,
        }
    }
}

fn merge_non_directional_results(
    original: SingleAlignmentResult,
    complementary: SingleAlignmentResult,
) -> SingleAlignmentResult {
    let original_distance = original.placement.map(ReadPlacement::distance);
    let complementary_distance = complementary.placement.map(ReadPlacement::distance);
    let (mut selected, other, tied) = match (original_distance, complementary_distance) {
        (Some(left), Some(right)) if left < right => (original, complementary, false),
        (Some(left), Some(right)) if right < left => (complementary, original, false),
        (Some(_), Some(_)) => (original, complementary, true),
        (None, Some(_)) => (complementary, original, false),
        (_, None) => (original, complementary, false),
    };
    selected.located_rows = original
        .located_rows
        .saturating_add(complementary.located_rows);
    selected.distinct_candidate_starts = original
        .distinct_candidate_starts
        .saturating_add(complementary.distinct_candidate_starts);
    selected.verified_placements = original
        .verified_placements
        .saturating_add(complementary.verified_placements);

    if tied {
        selected.status = SingleMappingStatus::Ambiguous;
        selected.mapping_quality = 0;
    } else if matches!(selected.status, SingleMappingStatus::Unique) {
        if let (Some(best), Some(runner_up)) = (
            selected.placement.map(ReadPlacement::distance),
            other.placement.map(ReadPlacement::distance),
        ) {
            selected.mapping_quality = selected
                .mapping_quality
                .min(cross_pass_mapping_quality_cap(best, runner_up));
        }
        if selected.located_rows > 256
            || selected.distinct_candidate_starts > 64
            || selected.verified_placements > 64
        {
            selected.mapping_quality = selected.mapping_quality.min(10);
        }
    } else {
        selected.mapping_quality = 0;
    }
    selected
}

const fn cross_pass_mapping_quality_cap(best_distance: u8, runner_up_distance: u8) -> u8 {
    let separation = runner_up_distance.saturating_sub(best_distance);
    match (best_distance, separation) {
        (_, 0) => 0,
        (5.., _) => 10,
        (4, _) | (_, 1) => 15,
        _ => u8::MAX,
    }
}

#[cfg(test)]
#[path = "../../tests/whitebox/single_end.rs"]
mod whitebox;
