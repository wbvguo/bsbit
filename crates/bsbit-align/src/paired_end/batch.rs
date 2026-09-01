//! Paired-batch orchestration and reusable worker state.

use super::adapter::supported_three_prime_adapter_start;
use super::frontier::{
    append_ranked_block_candidates, collect_ranked_block_seeds,
    selective_unmapped_frontier_deepening_required,
};
use super::rescue::{
    append_ungapped_semi_global_placements, exact_compatible_pair, exact_retained_placement,
    relabel_exact_retained_hit, rescue_from_combined_exact_blocks,
    rescue_from_ranked_anchor_windows,
};
use super::selection::{
    affine_placement_score, collapse_equivalent_pair_origins, pair_net_gap_profile,
    pair_origin_key, prefer_minimum_net_gap_representative,
    select_best_pair_origins_with_affine_score, select_best_pair_origins_with_endpoint_policy,
    select_best_pairs, select_best_pairs_with_fallback_score, select_reported_origin_endpoint,
    spatial_key,
};
use super::{
    ADAPTER_STABILITY_DELTA, AdapterFallbackResult, AffineScoreWorkspace, AlignmentError,
    BWA_MATCH_SCORE, BWA_MISMATCH_PENALTY, BWA_NEAR_SUBOPTIMAL_DELTA, Base, BisulfiteStrand,
    CombinedSearchLimits, CombinedSearchReferenceExt, CombinedSeedMatches,
    CombinedTwoLaneSearchState, ConversionPass, DIRECT_SINGLETON_PROOF, INITIAL_EDIT_DISTANCE,
    INITIAL_SEARCH_LIMITS, MAX_EDIT_DISTANCE, MAX_READ_BASES, MIN_SUFFIX_BASES,
    PAIRED_MAX_EDIT_DISTANCE, PairAlignmentMetrics, PairMappingStatus, PairWorkspace,
    PairedAlignmentOptions, PairedAlignmentResult, PairedAlignmentWorkMetrics, PairedBatchAligner,
    PairedBatchResult, PairedPlacement, PairedSearchMode, ProjectedBase, RankedBlockSeed,
    RankedBlockSelection, ReadAlignmentMetrics, ReadCandidate, ReadWorkspace, ReferenceIndex,
    SEMI_GLOBAL_CLIP_PENALTY, SEMI_GLOBAL_MAX_EXACT_ANCHOR_HITS, SEMI_GLOBAL_MIN_ALIGNED_BASES,
    SENSITIVE_MIN_EVENT_PENALTY, SENSITIVE_POSITIVE_MAPQ_REPORTING_MAX_RETAINED_HITS,
    SENSITIVE_POSITIVE_MAPQ_REPORTING_MIN_RETAINED_HITS, SENSITIVE_PROOF_BLOCKS,
    SENSITIVE_RANKED_BLOCK_HITS, SENSITIVE_REPEAT_RECHECK_ROWS,
    SENSITIVE_SELECTIVE_UNMAPPED_RANKED_BLOCK_HITS, SENSITIVE_UNMAPPED_RANKED_BLOCK_HITS,
    SearchBase, continue_combined_two_lane_search, paired_mapping_quality, placement_net_gap_bases,
    prepare_combined_projection, prepare_combined_search_projection,
    sensitive_ambiguity_q10_certified, sensitive_effective_mapping_quality,
    sensitive_incomplete_sparse_completion_required, sensitive_stable_rescue_q20_certified,
    sensitive_two_way_parsimony_q20_certified, sort_nominal_candidates,
    start_combined_two_lane_search,
};

impl PairedBatchAligner {
    /// Allocates reusable mapping storage for at least `pair_capacity` pairs.
    #[must_use]
    pub fn with_capacity(pair_capacity: usize) -> Self {
        Self {
            pair: PairWorkspace::with_capacity(4096, 1024, 32),
            projections: Vec::with_capacity(pair_capacity),
            first_seeds: Vec::with_capacity(pair_capacity.saturating_mul(2)),
            results: Vec::with_capacity(pair_capacity),
            collect_work_metrics: false,
            last_work_metrics: PairedAlignmentWorkMetrics::default(),
        }
    }

    /// Allocates reusable mapping storage and enables optional work counters.
    #[doc(hidden)]
    #[must_use]
    pub fn with_capacity_and_work_metrics(pair_capacity: usize) -> Self {
        Self {
            collect_work_metrics: true,
            ..Self::with_capacity(pair_capacity)
        }
    }

    /// Returns work performed by the most recent complete output mapping call.
    #[doc(hidden)]
    #[must_use]
    pub const fn last_work_metrics(&self) -> PairedAlignmentWorkMetrics {
        self.last_work_metrics
    }

    pub(super) fn observe_work_metrics(
        &mut self,
        results: &[PairedBatchResult],
        directional_passes_per_pair: u64,
    ) {
        if !self.collect_work_metrics {
            return;
        }
        for result in results {
            let metrics = result.metrics();
            let emitted_candidate_starts = metrics
                .mate1
                .emitted_candidate_starts
                .saturating_add(metrics.mate2.emitted_candidate_starts);
            let distinct_candidate_starts = metrics
                .mate1
                .distinct_candidate_starts
                .saturating_add(metrics.mate2.distinct_candidate_starts);
            let verified_placements = metrics
                .mate1
                .verified_placements
                .saturating_add(metrics.mate2.verified_placements);
            self.last_work_metrics.merge(PairedAlignmentWorkMetrics {
                pair_mapping_passes: directional_passes_per_pair,
                emitted_candidate_starts,
                distinct_candidate_starts,
                verified_placements,
                compatible_pairs: metrics.compatible_pairs,
                best_pair_placements: metrics.best_pair_placements,
            });
        }
    }

    /// Maps a directional or non-directional paired-end batch with one of the
    /// qualified paired-end strategies.
    ///
    /// # Errors
    ///
    /// Returns [`AlignmentError`] for an unsupported library profile,
    /// invalid template spans, or any mapping failure.
    fn map_pairs_combined<'a>(
        &'a mut self,
        reference: &ReferenceIndex,
        reads: &[[&[Base]; 2]],
        options: PairedAlignmentOptions,
    ) -> Result<&'a [PairedBatchResult], AlignmentError> {
        let (maximum_edit_distance, window_rescue, semi_global) = options.derived_policy();
        let passes = options.library_profile.conversion_passes();
        self.map_pair_conversion_pass(
            reference,
            reads,
            maximum_edit_distance,
            options.minimum_template_span,
            options.maximum_template_span,
            window_rescue,
            semi_global,
            options.search_mode,
            passes[0],
        )?;
        if passes.len() == 1 {
            return Ok(&self.results);
        }

        let original = self.results.clone();
        self.map_pair_conversion_pass(
            reference,
            reads,
            maximum_edit_distance,
            options.minimum_template_span,
            options.maximum_template_span,
            window_rescue,
            semi_global,
            options.search_mode,
            passes[1],
        )?;
        for (complementary, original) in self.results.iter_mut().zip(&original) {
            *complementary = merge_non_directional_batch_results(original, complementary);
        }
        Ok(&self.results)
    }

    /// Maps a paired-read batch through the complete qualified output policy.
    ///
    /// Adapter-supported trimming, stability remapping, MAPQ certificates,
    /// positive-MAPQ admission, and endpoint representation are resolved here
    /// so serialization callers receive facts rather than policy controls.
    ///
    /// # Errors
    ///
    /// Returns [`AlignmentError`] when the template span is invalid
    /// or any qualified mapping phase fails.
    ///
    /// # Panics
    ///
    /// Panics only if internally generated adapter-stability metadata loses
    /// its matching adapter result, which would violate this method's local
    /// construction invariant.
    // Adapter repair, stability proof, MAPQ certification, and reporting
    // admission form one ordered output-policy transaction.
    #[allow(clippy::too_many_lines)]
    pub fn map_pairs_for_output(
        &mut self,
        reference: &ReferenceIndex,
        reads: &[[&[Base]; 2]],
        options: PairedAlignmentOptions,
    ) -> Result<Vec<PairedAlignmentResult>, AlignmentError> {
        self.last_work_metrics = PairedAlignmentWorkMetrics::default();
        let directional_passes_per_pair =
            u64::from(options.library_profile.conversion_pass_count());
        let primary = self.map_pairs_combined(reference, reads, options)?.to_vec();
        self.observe_work_metrics(&primary, directional_passes_per_pair);
        let mut adapter_results = vec![None; reads.len()];
        let mut adapter_classes = vec![None; reads.len()];
        let mut adapter_attempted = vec![false; reads.len()];
        let mut adapter_clipped_mates = vec![0_u8; reads.len()];
        let mut adapter_clipped_bases = vec![0_usize; reads.len()];
        let mut clipped_reads = Vec::with_capacity(reads.len());
        let mut clipped_metadata = Vec::with_capacity(reads.len());

        for (offset, (pair, result)) in reads.iter().zip(&primary).enumerate() {
            let should_attempt = matches!(result.class(), PairMappingStatus::Unmapped)
                || (options.search_mode.is_sensitive()
                    && result.metrics().window_rescue_attempted
                    && matches!(result.class(), PairMappingStatus::Ambiguous));
            if !should_attempt {
                continue;
            }
            let retained = [
                supported_three_prime_adapter_start(pair[0])
                    .filter(|&start| start >= SEMI_GLOBAL_MIN_ALIGNED_BASES)
                    .unwrap_or(pair[0].len()),
                supported_three_prime_adapter_start(pair[1])
                    .filter(|&start| start >= SEMI_GLOBAL_MIN_ALIGNED_BASES)
                    .unwrap_or(pair[1].len()),
            ];
            if retained == [pair[0].len(), pair[1].len()] {
                continue;
            }
            adapter_attempted[offset] = true;
            adapter_clipped_mates[offset] = u8::from(retained[0] != pair[0].len())
                .saturating_add(u8::from(retained[1] != pair[1].len()));
            adapter_clipped_bases[offset] = pair[0]
                .len()
                .saturating_sub(retained[0])
                .saturating_add(pair[1].len().saturating_sub(retained[1]));
            clipped_reads.push([&pair[0][..retained[0]], &pair[1][..retained[1]]]);
            clipped_metadata.push((offset, retained));
        }

        if !clipped_reads.is_empty() {
            let adapter_options = PairedAlignmentOptions::adapter_trimmed(
                options.library_profile,
                options.search_mode,
                options.minimum_template_span,
                options.maximum_template_span,
            );
            let remapped = self
                .map_pairs_combined(reference, &clipped_reads, adapter_options)?
                .to_vec();
            self.observe_work_metrics(&remapped, directional_passes_per_pair);
            for ((offset, retained_bases), result) in clipped_metadata.iter().copied().zip(remapped)
            {
                adapter_results[offset] = Some(AdapterFallbackResult {
                    result,
                    stability_result: None,
                    final_class: result.class(),
                    retained_bases,
                });
            }

            let mut stability_reads = Vec::with_capacity(clipped_reads.len());
            let mut stability_metadata = Vec::with_capacity(clipped_reads.len());
            for (offset, fallback) in adapter_results.iter_mut().enumerate() {
                let Some(fallback) = fallback else {
                    continue;
                };
                if !matches!(fallback.final_class, PairMappingStatus::Unique) {
                    continue;
                }
                let full_lengths = [reads[offset][0].len(), reads[offset][1].len()];
                let mut retained = fallback.retained_bases;
                let mut stable_domain = true;
                for mate in 0..2 {
                    if retained[mate] == full_lengths[mate] {
                        continue;
                    }
                    if retained[mate]
                        < SEMI_GLOBAL_MIN_ALIGNED_BASES.saturating_add(ADAPTER_STABILITY_DELTA)
                    {
                        stable_domain = false;
                        break;
                    }
                    retained[mate] -= ADAPTER_STABILITY_DELTA;
                }
                if !stable_domain {
                    fallback.final_class = PairMappingStatus::Ambiguous;
                    continue;
                }
                stability_reads.push([
                    &reads[offset][0][..retained[0]],
                    &reads[offset][1][..retained[1]],
                ]);
                stability_metadata.push(offset);
            }

            if !stability_reads.is_empty() {
                let stability = self
                    .map_pairs_combined(reference, &stability_reads, adapter_options)?
                    .to_vec();
                self.observe_work_metrics(&stability, directional_passes_per_pair);
                for (offset, stability_result) in stability_metadata.into_iter().zip(stability) {
                    let fallback = adapter_results[offset]
                        .as_mut()
                        .expect("stability metadata refers to an adapter result");
                    fallback.stability_result = Some(stability_result);
                    let same_origin = fallback
                        .result
                        .best_pair()
                        .zip(stability_result.best_pair())
                        .is_some_and(|(primary, stability)| {
                            pair_origin_key(
                                primary,
                                fallback.retained_bases[0],
                                fallback.retained_bases[1],
                            ) == pair_origin_key(
                                stability,
                                fallback.retained_bases[0],
                                fallback.retained_bases[1],
                            )
                        });
                    if !matches!(stability_result.class(), PairMappingStatus::Unique)
                        || !same_origin
                    {
                        fallback.final_class = PairMappingStatus::Ambiguous;
                    }
                }
            }

            for (offset, fallback) in adapter_results.iter_mut().enumerate() {
                let Some(candidate) = *fallback else {
                    continue;
                };
                if matches!(primary[offset].class(), PairMappingStatus::Ambiguous)
                    && primary[offset].metrics().window_rescue_attempted
                    && !matches!(candidate.final_class, PairMappingStatus::Unique)
                {
                    adapter_classes[offset] = Some(PairMappingStatus::Ambiguous);
                    *fallback = None;
                } else {
                    adapter_classes[offset] = Some(candidate.final_class);
                }
            }
        }

        let mut outputs = Vec::with_capacity(reads.len());
        for (offset, (pair, strict_result)) in reads.iter().zip(primary).enumerate() {
            let adapter = adapter_results[offset];
            let result = adapter.map_or(strict_result, |fallback| fallback.result);
            let mut class = adapter.map_or(result.class(), |fallback| fallback.final_class);
            let mate_rescue_attempted = result.metrics().window_rescue_attempted
                || adapter.is_some_and(|fallback| {
                    fallback
                        .stability_result
                        .is_some_and(|stability| stability.metrics().window_rescue_attempted)
                });
            let semi_global_attempted = result.metrics().semi_global_attempted;
            let mut semi_global_clipped_mates = 0_u8;
            let mut semi_global_clipped_bases = 0_usize;
            if semi_global_attempted && let Some(selected) = result.best_pair() {
                for (placement, read) in [(selected.mate1(), pair[0]), (selected.mate2(), pair[1])]
                {
                    let retained = placement.retained_query_interval(read.len());
                    let clipped = read
                        .len()
                        .saturating_sub(retained.end.saturating_sub(retained.start));
                    semi_global_clipped_mates =
                        semi_global_clipped_mates.saturating_add(u8::from(clipped != 0));
                    semi_global_clipped_bases = semi_global_clipped_bases.saturating_add(clipped);
                }
            }

            let report_ambiguous =
                matches!(class, PairMappingStatus::Ambiguous) && result.best_pair().is_some();
            let Some(selected) = result
                .best_pair()
                .filter(|_| matches!(class, PairMappingStatus::Unique) || report_ambiguous)
            else {
                outputs.push(PairedAlignmentResult {
                    class,
                    placement: None,
                    retained_query_intervals: [0..pair[0].len(), 0..pair[1].len()],
                    mapping_quality: 0,
                    adapter_attempted: adapter_attempted[offset],
                    adapter_class: adapter_classes[offset],
                    adapter_clipped_mates: adapter_clipped_mates[offset],
                    adapter_clipped_bases: adapter_clipped_bases[offset],
                    semi_global_attempted,
                    semi_global_clipped_mates,
                    semi_global_clipped_bases,
                    mate_rescue_attempted,
                });
                continue;
            };
            let mut retained_query_intervals = adapter.map_or_else(
                || {
                    [
                        selected.mate1().retained_query_interval(pair[0].len()),
                        selected.mate2().retained_query_interval(pair[1].len()),
                    ]
                },
                |fallback| [0..fallback.retained_bases[0], 0..fallback.retained_bases[1]],
            );
            let mapping_quality = paired_mapping_quality(
                result,
                adapter.and_then(|fallback| fallback.stability_result),
                class,
                options.search_mode,
                [pair[0].len(), pair[1].len()],
                [&retained_query_intervals[0], &retained_query_intervals[1]],
            );
            let requires_positive_mapq = strict_result.requires_positive_mapq_for_reporting
                || result.requires_positive_mapq_for_reporting;
            if requires_positive_mapq && mapping_quality == 0 {
                class = PairMappingStatus::Unmapped;
                outputs.push(PairedAlignmentResult {
                    class,
                    placement: None,
                    retained_query_intervals,
                    mapping_quality,
                    adapter_attempted: adapter_attempted[offset],
                    adapter_class: adapter_classes[offset],
                    adapter_clipped_mates: adapter_clipped_mates[offset],
                    adapter_clipped_bases: adapter_clipped_bases[offset],
                    semi_global_attempted,
                    semi_global_clipped_mates,
                    semi_global_clipped_bases,
                    mate_rescue_attempted,
                });
                continue;
            }
            let selected = if adapter.is_none() && options.search_mode.is_sensitive() {
                select_reported_origin_endpoint(
                    reference,
                    *pair,
                    selected,
                    PAIRED_MAX_EDIT_DISTANCE,
                    options.minimum_template_span,
                    options.maximum_template_span,
                )
            } else {
                selected
            };
            retained_query_intervals = adapter.map_or_else(
                || {
                    [
                        selected.mate1().retained_query_interval(pair[0].len()),
                        selected.mate2().retained_query_interval(pair[1].len()),
                    ]
                },
                |fallback| [0..fallback.retained_bases[0], 0..fallback.retained_bases[1]],
            );
            outputs.push(PairedAlignmentResult {
                class,
                placement: Some(selected),
                retained_query_intervals,
                mapping_quality,
                adapter_attempted: adapter_attempted[offset],
                adapter_class: adapter_classes[offset],
                adapter_clipped_mates: adapter_clipped_mates[offset],
                adapter_clipped_bases: adapter_clipped_bases[offset],
                semi_global_attempted,
                semi_global_clipped_mates,
                semi_global_clipped_bases,
                mate_rescue_attempted,
            });
        }
        Ok(outputs)
    }

    #[allow(clippy::too_many_arguments)]
    fn map_pair_conversion_pass<'a>(
        &'a mut self,
        reference: &ReferenceIndex,
        reads: &[[&[Base]; 2]],
        maximum_edit_distance: u8,
        minimum_template_span: u64,
        maximum_template_span: u64,
        window_rescue: bool,
        semi_global: bool,
        search_mode: PairedSearchMode,
        conversion_pass: ConversionPass,
    ) -> Result<&'a [PairedBatchResult], AlignmentError> {
        if !conversion_pass.swaps_mates() {
            return self.map_directional_pairs_combined_inner(
                reference,
                reads,
                maximum_edit_distance,
                minimum_template_span,
                maximum_template_span,
                window_rescue,
                semi_global,
                search_mode,
            );
        }

        let swapped_reads = reads
            .iter()
            .map(|pair| [pair[1], pair[0]])
            .collect::<Vec<_>>();
        self.map_directional_pairs_combined_inner(
            reference,
            &swapped_reads,
            maximum_edit_distance,
            minimum_template_span,
            maximum_template_span,
            window_rescue,
            semi_global,
            search_mode,
        )?;
        for result in &mut self.results {
            *result = swap_batch_result_mates(*result);
        }
        Ok(&self.results)
    }

    #[allow(clippy::too_many_arguments)]
    // Cross-pair wavefront seeding and per-pair completion intentionally share
    // one batch workspace so projected reads and seed states are not copied.
    #[allow(clippy::too_many_lines)]
    fn map_directional_pairs_combined_inner<'a>(
        &'a mut self,
        reference: &ReferenceIndex,
        reads: &[[&[Base]; 2]],
        maximum_edit_distance: u8,
        minimum_template_span: u64,
        maximum_template_span: u64,
        window_rescue: bool,
        semi_global: bool,
        search_mode: PairedSearchMode,
    ) -> Result<&'a [PairedBatchResult], AlignmentError> {
        if maximum_edit_distance > MAX_EDIT_DISTANCE {
            return Err(AlignmentError::UnsupportedEditDistance {
                requested: maximum_edit_distance,
                maximum: MAX_EDIT_DISTANCE,
            });
        }
        if minimum_template_span > maximum_template_span {
            return Err(AlignmentError::InvertedTemplateSpan {
                minimum: minimum_template_span,
                maximum: maximum_template_span,
            });
        }
        if reads.len().saturating_mul(2) > 64 {
            return Err(AlignmentError::LocatedCountOverflow);
        }
        self.projections.clear();
        self.projections
            .resize(reads.len(), [[ProjectedBase::A; MAX_READ_BASES]; 2]);
        for (projection, pair) in self.projections.iter_mut().zip(reads) {
            prepare_combined_projection(pair[0], ConversionPass::Original, &mut projection[0])?;
            prepare_combined_projection(
                pair[1],
                ConversionPass::Complementary,
                &mut projection[1],
            )?;
        }
        let patterns = self
            .projections
            .iter()
            .zip(reads)
            .flat_map(|(projection, pair)| {
                [
                    &projection[0][..pair[0].len()],
                    &projection[1][..pair[1].len()],
                ]
            })
            .collect::<Vec<_>>();
        self.first_seeds = reference
            .combined_maximal_suffix_projected_wavefront(&patterns, MIN_SUFFIX_BASES)
            .map_err(|_| AlignmentError::CombinedIndex)?;
        self.pair.semi_global_clip_penalty = search_mode.semi_global_clip_penalty();
        self.pair.prefer_minimum_net_gap = search_mode.is_sensitive();
        // Sensitive mode uses semi-global alignment as a confidence repair,
        // not as an eager replacement objective. Run the proof-oriented
        // strict search first, preserve every already-high-confidence result,
        // and revisit only the small residual low-confidence frontier below.
        let eager_semi_global = semi_global && !search_mode.is_sensitive();
        self.results.clear();
        for (ordinal, pair) in reads.iter().enumerate() {
            let mut repeat_risk_q20_certified = false;
            let mut parsimony_q20_certified = false;
            let mut ambiguity_q10_certified = false;
            // The final paired-end frontier may append bounded completion work.
            let mut requires_positive_mapq_for_reporting = false;
            let first_seeds = [
                self.first_seeds[ordinal * 2],
                self.first_seeds[ordinal * 2 + 1],
            ];
            let projection = &self.projections[ordinal];
            let initial_edit_distance = maximum_edit_distance.min(INITIAL_EDIT_DISTANCE);
            let (mut class, mut metrics, mut second_best_distance) =
                self.pair.map_directional_pair_combined_prepared(
                    reference,
                    pair[0],
                    pair[1],
                    [
                        &projection[0][..pair[0].len()],
                        &projection[1][..pair[1].len()],
                    ],
                    first_seeds,
                    initial_edit_distance,
                    minimum_template_span,
                    maximum_template_span,
                    window_rescue,
                    eager_semi_global,
                    INITIAL_SEARCH_LIMITS,
                    true,
                )?;
            if matches!(class, PairMappingStatus::Unmapped) {
                if maximum_edit_distance > initial_edit_distance {
                    let reverified = self.pair.reverify_directional_pair_combined_candidates(
                        reference,
                        pair[0],
                        pair[1],
                        maximum_edit_distance,
                        minimum_template_span,
                        maximum_template_span,
                        eager_semi_global,
                        metrics,
                    )?;
                    if !matches!(reverified.0, PairMappingStatus::Unmapped) {
                        (class, metrics, second_best_distance) = reverified;
                    }
                }
                if matches!(class, PairMappingStatus::Unmapped) {
                    (class, metrics, second_best_distance) =
                        self.pair.continue_directional_pair_combined_incremental(
                            reference,
                            pair[0],
                            pair[1],
                            [
                                &projection[0][..pair[0].len()],
                                &projection[1][..pair[1].len()],
                            ],
                            first_seeds,
                            maximum_edit_distance,
                            minimum_template_span,
                            maximum_template_span,
                            window_rescue,
                            eager_semi_global,
                            metrics,
                        )?;
                }
            }
            let complete_suspicious_unique = sensitive_repeat_recheck_required(class, metrics);
            if search_mode.is_sensitive()
                && window_rescue
                && (matches!(class, PairMappingStatus::Unmapped) || complete_suspicious_unique)
            {
                let original = (class, metrics, second_best_distance);
                let original_best = complete_suspicious_unique
                    .then(|| self.pair.best_pairs().first().copied())
                    .flatten();
                let mut completed = self.pair.extend_directional_pair_from_ranked_blocks(
                    reference,
                    pair[0],
                    pair[1],
                    [
                        &projection[0][..pair[0].len()],
                        &projection[1][..pair[1].len()],
                    ],
                    maximum_edit_distance,
                    minimum_template_span,
                    maximum_template_span,
                    eager_semi_global,
                    SENSITIVE_RANKED_BLOCK_HITS,
                )?;
                if matches!(class, PairMappingStatus::Unmapped)
                    && matches!(
                        completed.as_ref().map(|candidate| candidate.0),
                        None | Some(PairMappingStatus::Unmapped)
                    )
                {
                    completed = self.pair.extend_directional_pair_from_ranked_blocks(
                        reference,
                        pair[0],
                        pair[1],
                        [
                            &projection[0][..pair[0].len()],
                            &projection[1][..pair[1].len()],
                        ],
                        maximum_edit_distance,
                        minimum_template_span,
                        maximum_template_span,
                        eager_semi_global,
                        SENSITIVE_UNMAPPED_RANKED_BLOCK_HITS,
                    )?;
                    if matches!(
                        completed.as_ref().map(|candidate| candidate.0),
                        None | Some(PairMappingStatus::Unmapped)
                    ) && self.pair.selective_unmapped_frontier_deepening_required(
                        completed.as_ref().map(|candidate| candidate.1),
                    ) {
                        let extended_frontier_requires_positive_mapq = self
                            .pair
                            .selective_unmapped_frontier_requires_positive_mapq();
                        completed = self.pair.extend_directional_pair_from_ranked_blocks(
                            reference,
                            pair[0],
                            pair[1],
                            [
                                &projection[0][..pair[0].len()],
                                &projection[1][..pair[1].len()],
                            ],
                            maximum_edit_distance,
                            minimum_template_span,
                            maximum_template_span,
                            eager_semi_global,
                            SENSITIVE_SELECTIVE_UNMAPPED_RANKED_BLOCK_HITS,
                        )?;
                        {
                            requires_positive_mapq_for_reporting =
                                completed.as_ref().is_some_and(|candidate| {
                                    !matches!(candidate.0, PairMappingStatus::Unmapped)
                                }) && extended_frontier_requires_positive_mapq;
                        }
                    }
                }
                match (original_best, completed) {
                    (_, Some(completed)) if !matches!(completed.0, PairMappingStatus::Unmapped) => {
                        let completed_best = self.pair.best_pairs().first().copied();
                        if original_best
                            .zip(completed_best)
                            .is_some_and(|(original, completed)| {
                                completed.score() > original.score()
                            })
                        {
                            self.pair.best_pairs.clear();
                            self.pair
                                .best_pairs
                                .push(original_best.expect("rescued unique has a best pair"));
                            (class, metrics, second_best_distance) = original;
                        } else {
                            (class, metrics, second_best_distance) = completed;
                            if let Some((original, completed)) = original_best.zip(completed_best)
                                && original.score() == completed.score()
                                && pair_origin_key(original, pair[0].len(), pair[1].len())
                                    != pair_origin_key(completed, pair[0].len(), pair[1].len())
                            {
                                self.pair.best_pairs.push(original);
                                collapse_equivalent_pair_origins(
                                    &mut self.pair.best_pairs,
                                    pair[0].len(),
                                    pair[1].len(),
                                );
                                prefer_minimum_net_gap_representative(
                                    &mut self.pair.best_pairs,
                                    pair[0].len(),
                                    pair[1].len(),
                                );
                                class = PairMappingStatus::Ambiguous;
                                metrics.best_pair_placements = metrics
                                    .best_pair_placements
                                    .max(
                                        u64::try_from(self.pair.best_pairs.len())
                                            .unwrap_or(u64::MAX),
                                    )
                                    .max(2);
                                metrics.near_best_pairings = metrics.near_best_pairings.max(1);
                            }
                        }
                    }
                    (Some(original_best), _) => {
                        self.pair.best_pairs.clear();
                        self.pair.best_pairs.push(original_best);
                        (class, metrics, second_best_distance) = original;
                    }
                    (None, _) => {}
                }
            }
            if search_mode.is_sensitive()
                && self
                    .pair
                    .should_affine_rescore(class, metrics, pair[0].len(), pair[1].len())
            {
                (class, metrics, second_best_distance) =
                    self.pair.affine_rescore_directional_pair(
                        reference,
                        pair[0],
                        pair[1],
                        class,
                        maximum_edit_distance,
                        minimum_template_span,
                        maximum_template_span,
                        metrics,
                    )?;
            }
            if search_mode.is_sensitive()
                && semi_global
                && sensitive_targeted_semi_global_required(
                    class,
                    metrics,
                    self.pair
                        .best_pairs()
                        .first()
                        .copied()
                        .map(PairedPlacement::distance),
                )
            {
                let original = (class, metrics, second_best_distance);
                let original_best = self.pair.best_pairs().first().copied();
                let original_confidence = sensitive_effective_mapping_quality(class, metrics);
                let mut candidate = self.pair.finish_directional_pair_combined(
                    reference,
                    pair[0],
                    pair[1],
                    maximum_edit_distance,
                    minimum_template_span,
                    maximum_template_span,
                    true,
                    metrics.mate1,
                    metrics.mate2,
                    metrics.window_rescue_attempted,
                )?;
                let candidate_best = self.pair.best_pairs().first().copied();
                let candidate_net_gap_profile =
                    pair_net_gap_profile(self.pair.best_pairs(), pair[0].len(), pair[1].len());
                let same_origin =
                    original_best
                        .zip(candidate_best)
                        .is_some_and(|(original, candidate)| {
                            pair_origin_key(original, pair[0].len(), pair[1].len())
                                == pair_origin_key(candidate, pair[0].len(), pair[1].len())
                        });
                ambiguity_q10_certified = sensitive_ambiguity_q10_certified(
                    original.0,
                    original.1,
                    original_best.map(PairedPlacement::score),
                    original_best.map(PairedPlacement::distance),
                    candidate.0,
                    candidate.1,
                    candidate_best.map(PairedPlacement::distance),
                    candidate_net_gap_profile,
                    same_origin,
                );
                let bounded_parsimony = sensitive_two_way_parsimony_q20_certified(
                    original.0,
                    original.1,
                    original_best,
                    candidate.0,
                    candidate.1,
                    candidate_best,
                    self.pair.best_pairs().len(),
                    candidate_net_gap_profile,
                    pair[0].len(),
                    pair[1].len(),
                    same_origin,
                );
                if bounded_parsimony {
                    // `prefer_minimum_net_gap_representative` has already put
                    // the uniquely parsimonious member first.  Keep the
                    // primary-score tie in the score metrics, but make the
                    // secondary decision explicit in class/cardinality.
                    self.pair.best_pairs.truncate(1);
                    candidate.0 = PairMappingStatus::Unique;
                    candidate.1.best_pair_placements = 1;
                }
                let candidate_confidence = if bounded_parsimony {
                    20
                } else {
                    sensitive_effective_mapping_quality(candidate.0, candidate.1)
                };
                repeat_risk_q20_certified = sensitive_stable_rescue_q20_certified(
                    original.0,
                    original.1,
                    original_confidence,
                    candidate.0,
                    candidate.1,
                    candidate_confidence,
                    same_origin,
                );
                // A first coordinate discovered only after endpoint clipping
                // has no independent full-read origin to stabilize it. Keep
                // it as unresolved; targeted semi-global is allowed to repair
                // confidence only at an origin already supported by the
                // strict/affine search.
                let completed_incomplete_q10 = !original.1.frontier_complete
                    && ambiguity_q10_certified
                    && candidate.1.frontier_complete
                    && same_origin;
                if completed_incomplete_q10 {
                    // Completion proved a stable representative, but this
                    // The qualified cell is only a Q10 ambiguity
                    // certificate.  Do not let the candidate's temporary
                    // Unique class enter the unique-result Q30/Q40 LUT.
                    candidate.0 = PairMappingStatus::Ambiguous;
                    candidate.1.best_pair_placements = candidate.1.best_pair_placements.max(2);
                }
                let accepted = same_origin
                    && ((original.1.frontier_complete
                        && (bounded_parsimony || candidate_confidence > original_confidence))
                        || completed_incomplete_q10);
                parsimony_q20_certified = bounded_parsimony && accepted;
                if accepted {
                    if bounded_parsimony {
                        // The targeted pass is secondary confidence evidence
                        // at the already supported biological origin.  Keep
                        // the strict representative for traceback/CIGAR: its
                        // full-read placement remains the primary alignment,
                        // while the uniquely parsimonious endpoint tie-break
                        // only certifies class and the Q20 boundary.
                        self.pair.best_pairs.clear();
                        self.pair
                            .best_pairs
                            .push(original_best.expect("certified parsimony has a strict origin"));
                    }
                    (class, metrics, second_best_distance) = candidate;
                } else {
                    self.pair.best_pairs.clear();
                    if let Some(original_best) = original_best {
                        self.pair.best_pairs.push(original_best);
                    }
                    (class, metrics, second_best_distance) = original;
                }
            }
            self.results.push(PairedBatchResult {
                class,
                metrics,
                best_pair: self.pair.best_pairs().first().copied(),
                second_best_distance,
                repeat_risk_q20_certified,
                parsimony_q20_certified,
                ambiguity_q10_certified,
                requires_positive_mapq_for_reporting,
            });
        }
        Ok(&self.results)
    }
}

pub(super) fn swap_batch_result_mates(mut result: PairedBatchResult) -> PairedBatchResult {
    result.metrics = PairAlignmentMetrics {
        mate1: result.metrics.mate2,
        mate2: result.metrics.mate1,
        ..result.metrics
    };
    result.best_pair = result.best_pair.map(|pair| PairedPlacement {
        mate1: pair.mate2,
        mate2: pair.mate1,
        ..pair
    });
    result
}

// Classification, score confidence, and MAPQ evidence must be merged under
// the same directional tie decision.
#[allow(clippy::too_many_lines)]
pub(super) fn merge_non_directional_batch_results(
    original: &PairedBatchResult,
    complementary: &PairedBatchResult,
) -> PairedBatchResult {
    let original_score = batch_best_score(original);
    let complementary_score = batch_best_score(complementary);
    let (mut selected, other, tied) = match (original_score, complementary_score) {
        (Some(left), Some(right)) if left > right => (*original, complementary, false),
        (Some(left), Some(right)) if right > left => (*complementary, original, false),
        (Some(_), Some(_)) => (*original, complementary, true),
        (None, Some(_)) => (*complementary, original, false),
        (_, None) => (*original, complementary, false),
    };
    let selected_score = batch_best_score(&selected);
    let other_score = batch_best_score(other);
    let selected_metrics = selected.metrics;
    let other_metrics = other.metrics;

    selected.metrics = PairAlignmentMetrics {
        mate1: merge_read_metrics(selected_metrics.mate1, other_metrics.mate1),
        mate2: merge_read_metrics(selected_metrics.mate2, other_metrics.mate2),
        compatible_pairs: selected_metrics
            .compatible_pairs
            .saturating_add(other_metrics.compatible_pairs),
        best_pair_placements: if tied {
            selected_metrics
                .best_pair_placements
                .saturating_add(other_metrics.best_pair_placements)
                .max(2)
        } else {
            selected_metrics.best_pair_placements
        },
        window_rescue_attempted: selected_metrics.window_rescue_attempted
            || other_metrics.window_rescue_attempted,
        semi_global_attempted: selected_metrics.semi_global_attempted
            || other_metrics.semi_global_attempted,
        best_pair_score: selected_score,
        second_best_pair_score: if tied {
            selected_score
        } else {
            maximum_optional_score(selected_metrics.second_best_pair_score, other_score)
        },
        near_best_pairings: selected_metrics
            .near_best_pairings
            .saturating_add(other_score.map_or(0, |other_score| {
                selected_score.map_or(0, |selected_score| {
                    if selected_score.saturating_sub(other_score) <= BWA_NEAR_SUBOPTIMAL_DELTA {
                        other_metrics.near_best_pairings.saturating_add(1)
                    } else {
                        0
                    }
                })
            })),
        mapq_compatible_pairs: selected_metrics
            .mapq_compatible_pairs
            .saturating_add(other_metrics.mapq_compatible_pairs),
        mapq_best_pair_score: maximum_optional_score(
            selected_metrics.mapq_best_pair_score,
            other_metrics.mapq_best_pair_score,
        ),
        mapq_second_best_pair_score: if tied {
            maximum_optional_score(
                selected_metrics.mapq_best_pair_score,
                other_metrics.mapq_best_pair_score,
            )
        } else {
            maximum_optional_score(
                selected_metrics.mapq_second_best_pair_score,
                other_metrics.mapq_best_pair_score,
            )
        },
        mapq_near_best_pairings: selected_metrics.mapq_near_best_pairings.saturating_add(
            other_metrics.mapq_best_pair_score.map_or(0, |other_score| {
                selected_metrics
                    .mapq_best_pair_score
                    .map_or(0, |selected_score| {
                        if selected_score.saturating_sub(other_score) <= BWA_NEAR_SUBOPTIMAL_DELTA {
                            other_metrics.mapq_near_best_pairings.saturating_add(1)
                        } else {
                            0
                        }
                    })
            }),
        ),
        frontier_complete: selected_metrics.frontier_complete && other_metrics.frontier_complete,
    };
    selected.second_best_distance = if tied {
        selected.best_pair.map(PairedPlacement::score)
    } else {
        minimum_optional_distance(
            selected.second_best_distance,
            other.best_pair.map(PairedPlacement::score),
        )
    };
    if tied {
        selected.class = PairMappingStatus::Ambiguous;
        selected.repeat_risk_q20_certified = false;
        selected.parsimony_q20_certified = false;
    } else if !selected.metrics.frontier_complete
        && matches!(selected.class, PairMappingStatus::Unique)
    {
        selected.class = PairMappingStatus::Ambiguous;
        selected.metrics.best_pair_placements = selected.metrics.best_pair_placements.max(2);
        selected.repeat_risk_q20_certified = false;
        selected.parsimony_q20_certified = false;
    } else if other_score.is_some() {
        // Cross-conversion evidence did not participate in the directional
        // pass's specialized confidence repair, so do not carry that repair
        // certificate across the global four-strand decision.
        selected.repeat_risk_q20_certified = false;
        selected.parsimony_q20_certified = false;
    }
    selected
}

fn batch_best_score(result: &PairedBatchResult) -> Option<i16> {
    result.metrics.best_pair_score.or_else(|| {
        result
            .best_pair
            .map(|pair| -i16::from(pair.score()).saturating_mul(BWA_MISMATCH_PENALTY))
    })
}

const fn maximum_optional_score(left: Option<i16>, right: Option<i16>) -> Option<i16> {
    match (left, right) {
        (Some(left), Some(right)) => Some(if left > right { left } else { right }),
        (Some(value), None) | (None, Some(value)) => Some(value),
        (None, None) => None,
    }
}

const fn minimum_optional_distance(left: Option<u8>, right: Option<u8>) -> Option<u8> {
    match (left, right) {
        (Some(left), Some(right)) => Some(if left < right { left } else { right }),
        (Some(value), None) | (None, Some(value)) => Some(value),
        (None, None) => None,
    }
}

const fn merge_read_metrics(
    left: ReadAlignmentMetrics,
    right: ReadAlignmentMetrics,
) -> ReadAlignmentMetrics {
    ReadAlignmentMetrics {
        located_rows: left.located_rows.saturating_add(right.located_rows),
        emitted_candidate_starts: left
            .emitted_candidate_starts
            .saturating_add(right.emitted_candidate_starts),
        distinct_candidate_starts: left
            .distinct_candidate_starts
            .saturating_add(right.distinct_candidate_starts),
        verified_placements: left
            .verified_placements
            .saturating_add(right.verified_placements),
    }
}

pub(super) fn sensitive_repeat_recheck_required(
    class: PairMappingStatus,
    metrics: PairAlignmentMetrics,
) -> bool {
    matches!(class, PairMappingStatus::Unique)
        && (metrics.window_rescue_attempted
            || metrics.mate1.located_rows.max(metrics.mate2.located_rows)
                >= SENSITIVE_REPEAT_RECHECK_ROWS)
}

pub(super) fn sensitive_targeted_semi_global_required(
    class: PairMappingStatus,
    metrics: PairAlignmentMetrics,
    pair_distance: Option<u8>,
) -> bool {
    !matches!(class, PairMappingStatus::Unmapped)
        && ((metrics.frontier_complete && sensitive_effective_mapping_quality(class, metrics) < 20)
            || sensitive_incomplete_sparse_completion_required(class, metrics, pair_distance))
}

pub(super) fn conservatively_mark_incomplete_frontier(
    result: &mut (PairMappingStatus, PairAlignmentMetrics, Option<u8>),
    complete: bool,
) {
    result.1.frontier_complete = complete;
    if !complete && matches!(result.0, PairMappingStatus::Unique) {
        result.0 = PairMappingStatus::Ambiguous;
        result.1.best_pair_placements = result.1.best_pair_placements.max(2);
        result.2 = None;
    }
}

impl PairWorkspace {
    #[must_use]
    pub(super) fn with_capacity(
        mate_candidate_capacity: usize,
        mate_placement_capacity: usize,
        pair_capacity: usize,
    ) -> Self {
        Self {
            mate1: ReadWorkspace::with_capacity(mate_candidate_capacity, mate_placement_capacity),
            mate2: ReadWorkspace::with_capacity(mate_candidate_capacity, mate_placement_capacity),
            rescue_windows: Vec::with_capacity(mate_placement_capacity),
            best_pairs: Vec::with_capacity(pair_capacity),
            exact_anchor_candidates: Vec::with_capacity(
                usize::try_from(SEMI_GLOBAL_MAX_EXACT_ANCHOR_HITS)
                    .expect("exact anchor limit fits usize"),
            ),
            ranked_anchor_placements: Vec::with_capacity(mate_placement_capacity),
            mate1_affine_scores: Vec::with_capacity(mate_placement_capacity),
            mate2_affine_scores: Vec::with_capacity(mate_placement_capacity),
            affine: AffineScoreWorkspace::default(),
            semi_global_clip_penalty: SEMI_GLOBAL_CLIP_PENALTY,
            prefer_minimum_net_gap: false,
            origin_pair_evidence: std::collections::HashMap::with_capacity(pair_capacity),
            combined_search_state: CombinedTwoLaneSearchState::new(),
            fallback_mate1_nominals: Vec::with_capacity(mate_candidate_capacity / 5 + 1),
            fallback_mate2_nominals: Vec::with_capacity(mate_candidate_capacity / 5 + 1),
            ranked_extension_selections: [None; 2],
        }
    }

    fn selective_unmapped_frontier_deepening_required(
        &self,
        ordinary_frontier_metrics: Option<PairAlignmentMetrics>,
    ) -> bool {
        selective_unmapped_frontier_deepening_required(
            self.ranked_extension_selections,
            ordinary_frontier_metrics,
        )
    }

    fn selective_unmapped_frontier_requires_positive_mapq(&self) -> bool {
        let [Some(first), Some(second)] = self.ranked_extension_selections else {
            return false;
        };
        !(SENSITIVE_POSITIVE_MAPQ_REPORTING_MIN_RETAINED_HITS
            ..SENSITIVE_POSITIVE_MAPQ_REPORTING_MAX_RETAINED_HITS)
            .contains(&first.retained_hits.saturating_add(second.retained_hits))
    }

    #[allow(clippy::too_many_arguments)]
    fn collect_ranked_block_seeds_for_lane(
        reference: &ReferenceIndex,
        read: &[Base],
        reversed_projected: &[ProjectedBase],
        lane: usize,
        maximum_edit_distance: u8,
        maximum_ranked_block_hits: u64,
        output: &mut [Option<RankedBlockSeed>; SENSITIVE_PROOF_BLOCKS],
    ) -> Result<Option<RankedBlockSelection>, AlignmentError> {
        let budget = usize::from(maximum_edit_distance);
        debug_assert!(lane < 2);
        debug_assert!(budget < SENSITIVE_PROOF_BLOCKS);
        let selection = collect_ranked_block_seeds(
            reference,
            read,
            reversed_projected,
            lane == 0,
            maximum_edit_distance,
            maximum_ranked_block_hits,
            output,
        )?;
        Ok(selection)
    }

    fn append_ranked_block_candidates_for_lane(
        &mut self,
        reference: &ReferenceIndex,
        read_len: usize,
        lane: usize,
        maximum_edit_distance: u8,
        seeds: &[Option<RankedBlockSeed>; SENSITIVE_PROOF_BLOCKS],
    ) -> Result<u64, AlignmentError> {
        let budget = usize::from(maximum_edit_distance);
        debug_assert!(lane < 2);
        debug_assert!(budget < SENSITIVE_PROOF_BLOCKS);
        let candidates = if lane == 0 {
            &mut self.mate1.candidate_nominals
        } else {
            &mut self.mate2.candidate_nominals
        };
        let located_rows =
            append_ranked_block_candidates(reference, read_len, lane == 0, seeds, candidates)?;
        Ok(located_rows)
    }

    // This is the internal handoff of one fully prepared pair; a parameter
    // object would obscure which search proof each value belongs to.
    #[allow(clippy::too_many_arguments)]
    fn map_directional_pair_combined_prepared(
        &mut self,
        reference: &ReferenceIndex,
        read1: &[Base],
        read2: &[Base],
        projected: [&[ProjectedBase]; 2],
        first_seeds: [Option<CombinedSeedMatches>; 2],
        maximum_edit_distance: u8,
        minimum_template_span: u64,
        maximum_template_span: u64,
        window_rescue: bool,
        semi_global: bool,
        search_limits: CombinedSearchLimits,
        preserve_fallback_frontier: bool,
    ) -> Result<(PairMappingStatus, PairAlignmentMetrics, Option<u8>), AlignmentError> {
        {
            self.mate1.begin_verification_cache_read();
            self.mate2.begin_verification_cache_read();
        }
        self.best_pairs.clear();
        self.mate1.candidates.clear();
        self.mate1.candidate_nominals.clear();
        self.mate1.placements.clear();
        self.mate2.candidates.clear();
        self.mate2.candidate_nominals.clear();
        self.mate2.placements.clear();
        self.combined_search_state = CombinedTwoLaneSearchState::new();
        self.fallback_mate1_nominals.clear();
        self.fallback_mate2_nominals.clear();
        if !semi_global
            && (read1.iter().filter(|base| base.is_unknown()).count()
                > usize::from(maximum_edit_distance)
                || read2.iter().filter(|base| base.is_unknown()).count()
                    > usize::from(maximum_edit_distance))
        {
            return Ok((PairMappingStatus::Unmapped, empty_pair_metrics(), None));
        }
        self.combined_search_state = start_combined_two_lane_search(
            reference,
            [read1, read2],
            projected,
            first_seeds,
            [ConversionPass::Original, ConversionPass::Complementary],
            search_limits,
            &mut self.mate1.candidate_nominals,
            &mut self.mate2.candidate_nominals,
        )?;
        let located_rows = self.combined_search_state.located;
        let mate1_metrics = ReadAlignmentMetrics {
            located_rows: located_rows[0],
            ..ReadAlignmentMetrics::default()
        };
        let mate2_metrics = ReadAlignmentMetrics {
            located_rows: located_rows[1],
            ..ReadAlignmentMetrics::default()
        };
        self.verify_directional_pair_combined_frontier(
            reference,
            read1,
            read2,
            projected,
            first_seeds,
            maximum_edit_distance,
            minimum_template_span,
            maximum_template_span,
            window_rescue,
            semi_global,
            search_limits,
            preserve_fallback_frontier,
            mate1_metrics,
            mate2_metrics,
        )
    }

    // Geometry filtering, rescue selection, and verification share one
    // candidate frontier and must update its metrics atomically.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn verify_directional_pair_combined_frontier(
        &mut self,
        reference: &ReferenceIndex,
        read1: &[Base],
        read2: &[Base],
        projected: [&[ProjectedBase]; 2],
        first_seeds: [Option<CombinedSeedMatches>; 2],
        maximum_edit_distance: u8,
        minimum_template_span: u64,
        maximum_template_span: u64,
        window_rescue: bool,
        semi_global: bool,
        search_limits: CombinedSearchLimits,
        preserve_fallback_frontier: bool,
        mate1_metrics: ReadAlignmentMetrics,
        mate2_metrics: ReadAlignmentMetrics,
    ) -> Result<(PairMappingStatus, PairAlignmentMetrics, Option<u8>), AlignmentError> {
        sort_nominal_candidates(&mut self.mate1.candidate_nominals);
        sort_nominal_candidates(&mut self.mate2.candidate_nominals);
        let nominal_geometry = nominal_pair_geometry_exists(
            &self.mate1.candidate_nominals,
            &self.mate2.candidate_nominals,
            read1.len(),
            read2.len(),
            maximum_template_span,
            maximum_edit_distance,
        );
        if preserve_fallback_frontier && !nominal_geometry {
            self.fallback_mate1_nominals
                .clone_from(&self.mate1.candidate_nominals);
            self.fallback_mate2_nominals
                .clone_from(&self.mate2.candidate_nominals);
        }
        let rescue_anchor = if window_rescue && !nominal_geometry {
            select_combined_window_rescue_anchor(
                first_seeds,
                &self.mate1.candidate_nominals,
                &self.mate2.candidate_nominals,
            )
        } else {
            None
        };
        let (mate1_metrics, mate2_metrics, window_rescue_attempted) = match rescue_anchor {
            Some(0) => {
                let (anchors, mate1_metrics) = self.mate1.verify_sorted_candidates_with_budget(
                    reference,
                    read1,
                    mate1_metrics,
                    maximum_edit_distance,
                )?;
                if anchors.is_empty() {
                    self.mate2.candidate_nominals.clear();
                    self.mate2.candidates.clear();
                    self.mate2.placements.clear();
                    (mate1_metrics, mate2_metrics, false)
                } else {
                    let mut rescued_metrics = rescue_from_combined_exact_blocks(
                        &mut self.mate2,
                        &mut self.rescue_windows,
                        reference,
                        read2,
                        projected[1],
                        anchors,
                        false,
                        maximum_template_span,
                        maximum_edit_distance,
                        preserve_fallback_frontier,
                        search_limits.maximum_combined_rescue_hits,
                    )?;
                    rescued_metrics.located_rows = mate2_metrics.located_rows;
                    (mate1_metrics, rescued_metrics, true)
                }
            }
            Some(1) => {
                let (anchors, mate2_metrics) = self.mate2.verify_sorted_candidates_with_budget(
                    reference,
                    read2,
                    mate2_metrics,
                    maximum_edit_distance,
                )?;
                if anchors.is_empty() {
                    self.mate1.candidate_nominals.clear();
                    self.mate1.candidates.clear();
                    self.mate1.placements.clear();
                    (mate1_metrics, mate2_metrics, false)
                } else {
                    let mut rescued_metrics = rescue_from_combined_exact_blocks(
                        &mut self.mate1,
                        &mut self.rescue_windows,
                        reference,
                        read1,
                        projected[0],
                        anchors,
                        true,
                        maximum_template_span,
                        maximum_edit_distance,
                        preserve_fallback_frontier,
                        search_limits.maximum_combined_rescue_hits,
                    )?;
                    rescued_metrics.located_rows = mate1_metrics.located_rows;
                    (rescued_metrics, mate2_metrics, true)
                }
            }
            Some(_) => unreachable!("a pair has exactly two mates"),
            None => {
                retain_nominal_pair_geometry(
                    &mut self.mate1.candidate_nominals,
                    &mut self.mate2.candidate_nominals,
                    read1.len(),
                    read2.len(),
                    maximum_template_span,
                    maximum_edit_distance,
                );
                let (_, mate1_metrics) = self.mate1.verify_sorted_candidates_with_budget(
                    reference,
                    read1,
                    mate1_metrics,
                    maximum_edit_distance,
                )?;
                let (_, mate2_metrics) = self.mate2.verify_sorted_candidates_with_budget(
                    reference,
                    read2,
                    mate2_metrics,
                    maximum_edit_distance,
                )?;
                (mate1_metrics, mate2_metrics, false)
            }
        };
        self.finish_directional_pair_combined(
            reference,
            read1,
            read2,
            maximum_edit_distance,
            minimum_template_span,
            maximum_template_span,
            semi_global,
            mate1_metrics,
            mate2_metrics,
            window_rescue_attempted,
        )
    }

    /// Reuses the candidate frontier left by the initial d3 pass and only reruns
    /// bounded verification at the requested edit budget. A pair that remains
    /// unmapped can still fall through to the deeper incremental seed search.
    #[allow(clippy::too_many_arguments)]
    fn reverify_directional_pair_combined_candidates(
        &mut self,
        reference: &ReferenceIndex,
        read1: &[Base],
        read2: &[Base],
        maximum_edit_distance: u8,
        minimum_template_span: u64,
        maximum_template_span: u64,
        semi_global: bool,
        previous_metrics: PairAlignmentMetrics,
    ) -> Result<(PairMappingStatus, PairAlignmentMetrics, Option<u8>), AlignmentError> {
        self.best_pairs.clear();
        self.mate1.candidates.clear();
        self.mate1.placements.clear();
        self.mate2.candidates.clear();
        self.mate2.placements.clear();
        self.mate1
            .candidate_nominals
            .extend_from_slice(&self.fallback_mate1_nominals);
        self.mate2
            .candidate_nominals
            .extend_from_slice(&self.fallback_mate2_nominals);
        sort_nominal_candidates(&mut self.mate1.candidate_nominals);
        sort_nominal_candidates(&mut self.mate2.candidate_nominals);
        retain_nominal_pair_geometry(
            &mut self.mate1.candidate_nominals,
            &mut self.mate2.candidate_nominals,
            read1.len(),
            read2.len(),
            maximum_template_span,
            maximum_edit_distance,
        );
        let mate1_metrics = ReadAlignmentMetrics {
            located_rows: previous_metrics.mate1.located_rows,
            ..ReadAlignmentMetrics::default()
        };
        let mate2_metrics = ReadAlignmentMetrics {
            located_rows: previous_metrics.mate2.located_rows,
            ..ReadAlignmentMetrics::default()
        };
        let (_, mate1_metrics) = self.mate1.verify_sorted_candidates_with_budget(
            reference,
            read1,
            mate1_metrics,
            maximum_edit_distance,
        )?;
        let (_, mate2_metrics) = self.mate2.verify_sorted_candidates_with_budget(
            reference,
            read2,
            mate2_metrics,
            maximum_edit_distance,
        )?;
        self.finish_directional_pair_combined(
            reference,
            read1,
            read2,
            maximum_edit_distance,
            minimum_template_span,
            maximum_template_span,
            semi_global,
            mate1_metrics,
            mate2_metrics,
            previous_metrics.window_rescue_attempted,
        )
    }

    /// Replays only initial-pass seeds that become admissible in the incremental
    /// fallback,
    /// then continues from the saved seed offset into the additional round.
    /// This preserves the deeper frontier without repeating rounds 0 through 4.
    #[allow(clippy::too_many_arguments)]
    fn continue_directional_pair_combined_incremental(
        &mut self,
        reference: &ReferenceIndex,
        read1: &[Base],
        read2: &[Base],
        projected: [&[ProjectedBase]; 2],
        first_seeds: [Option<CombinedSeedMatches>; 2],
        maximum_edit_distance: u8,
        minimum_template_span: u64,
        maximum_template_span: u64,
        window_rescue: bool,
        semi_global: bool,
        previous_metrics: PairAlignmentMetrics,
    ) -> Result<(PairMappingStatus, PairAlignmentMetrics, Option<u8>), AlignmentError> {
        if !self.combined_search_state.initialized {
            return self.map_directional_pair_combined_prepared(
                reference,
                read1,
                read2,
                projected,
                first_seeds,
                maximum_edit_distance,
                minimum_template_span,
                maximum_template_span,
                window_rescue,
                semi_global,
                PairedSearchMode::Default.limits(),
                false,
            );
        }
        self.best_pairs.clear();
        self.mate1.candidates.clear();
        self.mate1.placements.clear();
        self.mate2.candidates.clear();
        self.mate2.placements.clear();
        self.mate1
            .candidate_nominals
            .append(&mut self.fallback_mate1_nominals);
        self.mate2
            .candidate_nominals
            .append(&mut self.fallback_mate2_nominals);
        let additional_located = continue_combined_two_lane_search(
            reference,
            [read1, read2],
            projected,
            [ConversionPass::Original, ConversionPass::Complementary],
            &mut self.combined_search_state,
            &mut self.mate1.candidate_nominals,
            &mut self.mate2.candidate_nominals,
        )?;
        let mate1_metrics = ReadAlignmentMetrics {
            located_rows: previous_metrics
                .mate1
                .located_rows
                .saturating_add(additional_located[0]),
            ..ReadAlignmentMetrics::default()
        };
        let mate2_metrics = ReadAlignmentMetrics {
            located_rows: previous_metrics
                .mate2
                .located_rows
                .saturating_add(additional_located[1]),
            ..ReadAlignmentMetrics::default()
        };
        self.verify_directional_pair_combined_frontier(
            reference,
            read1,
            read2,
            projected,
            first_seeds,
            maximum_edit_distance,
            minimum_template_span,
            maximum_template_span,
            window_rescue,
            semi_global,
            PairedSearchMode::Default.limits(),
            false,
            mate1_metrics,
            mate2_metrics,
        )
    }

    /// Sensitive failed-pair completion.
    ///
    /// The `d + 1` disjoint exact blocks are ranked by occurrence count and the
    /// rarest intervals fitting the global enumeration budget form the anchor
    /// frontier. Every verified anchor then induces a bounded window in which
    /// the partner is completed with the full disjoint-block proof. A result
    /// recovered from an incomplete global block frontier is reported only as
    /// ambiguous; no truth coordinate participates in paired-end
    /// classification.
    // The ranked proof pipeline intentionally retains its seed, verification,
    // and fallback frontiers in one worker-owned transaction.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn extend_directional_pair_from_ranked_blocks(
        &mut self,
        reference: &ReferenceIndex,
        read1: &[Base],
        read2: &[Base],
        projected: [&[ProjectedBase]; 2],
        maximum_edit_distance: u8,
        minimum_template_span: u64,
        maximum_template_span: u64,
        semi_global: bool,
        maximum_ranked_block_hits: u64,
    ) -> Result<Option<(PairMappingStatus, PairAlignmentMetrics, Option<u8>)>, AlignmentError> {
        let mut seed_sets = [[None; SENSITIVE_PROOF_BLOCKS]; 2];
        let first_selection = Self::collect_ranked_block_seeds_for_lane(
            reference,
            read1,
            projected[0],
            0,
            maximum_edit_distance,
            maximum_ranked_block_hits,
            &mut seed_sets[0],
        )?;
        let second_selection = Self::collect_ranked_block_seeds_for_lane(
            reference,
            read2,
            projected[1],
            1,
            maximum_edit_distance,
            maximum_ranked_block_hits,
            &mut seed_sets[1],
        )?;
        {
            self.ranked_extension_selections = [first_selection, second_selection];
        }
        // The paired geometry is much more selective than either short
        // bisulfite block by itself.  Intersect both bounded frontiers before
        // invoking the d5 verifier; the anchor/window path below remains the
        // fallback when only one mate has an informative retained block.
        if let (Some(first_ranked), Some(second_ranked)) = (first_selection, second_selection) {
            self.best_pairs.clear();
            self.mate1.candidates.clear();
            self.mate1.candidate_nominals.clear();
            self.mate1.placements.clear();
            self.mate2.candidates.clear();
            self.mate2.candidate_nominals.clear();
            self.mate2.placements.clear();
            let mate1_rows = self.append_ranked_block_candidates_for_lane(
                reference,
                read1.len(),
                0,
                maximum_edit_distance,
                &seed_sets[0],
            )?;
            let mate2_rows = self.append_ranked_block_candidates_for_lane(
                reference,
                read2.len(),
                1,
                maximum_edit_distance,
                &seed_sets[1],
            )?;
            sort_nominal_candidates(&mut self.mate1.candidate_nominals);
            sort_nominal_candidates(&mut self.mate2.candidate_nominals);
            retain_nominal_pair_geometry(
                &mut self.mate1.candidate_nominals,
                &mut self.mate2.candidate_nominals,
                read1.len(),
                read2.len(),
                maximum_template_span,
                maximum_edit_distance,
            );
            if !self.mate1.candidate_nominals.is_empty()
                && !self.mate2.candidate_nominals.is_empty()
            {
                let mate1_seed_metrics = ReadAlignmentMetrics {
                    located_rows: mate1_rows,
                    ..ReadAlignmentMetrics::default()
                };
                let mate2_seed_metrics = ReadAlignmentMetrics {
                    located_rows: mate2_rows,
                    ..ReadAlignmentMetrics::default()
                };
                let (_, mate1_metrics) = self.mate1.verify_sorted_candidates_with_budget(
                    reference,
                    read1,
                    mate1_seed_metrics,
                    maximum_edit_distance,
                )?;
                let (_, mate2_metrics) = self.mate2.verify_sorted_candidates_with_budget(
                    reference,
                    read2,
                    mate2_seed_metrics,
                    maximum_edit_distance,
                )?;
                let mut completed = self.finish_directional_pair_combined(
                    reference,
                    read1,
                    read2,
                    maximum_edit_distance,
                    minimum_template_span,
                    maximum_template_span,
                    semi_global,
                    mate1_metrics,
                    mate2_metrics,
                    false,
                )?;
                if !matches!(completed.0, PairMappingStatus::Unmapped) {
                    let complete_frontier = first_ranked.complete && second_ranked.complete;
                    if !complete_frontier
                        && let Some(certified) = self.certify_ranked_pair_frontier(
                            reference,
                            read1,
                            read2,
                            projected,
                            maximum_edit_distance,
                            minimum_template_span,
                            maximum_template_span,
                            semi_global,
                            SENSITIVE_RANKED_BLOCK_HITS.saturating_mul(2),
                            completed,
                        )?
                    {
                        return Ok(Some(certified));
                    }
                    conservatively_mark_incomplete_frontier(&mut completed, complete_frontier);
                    return Ok(Some(completed));
                }
            }
        }

        let anchor_lane = match (first_selection, second_selection) {
            (None, None) => return Ok(None),
            (Some(_), None) => 0,
            (None, Some(_)) => 1,
            (Some(first), Some(second)) => match (first.complete, second.complete) {
                (true, false) => 0,
                (false, true) => 1,
                _ => usize::from(second.retained_hits < first.retained_hits),
            },
        };
        let anchor_frontier_complete = [first_selection, second_selection][anchor_lane]
            .expect("selected anchor has retained hits")
            .complete;

        self.best_pairs.clear();
        self.mate1.candidates.clear();
        self.mate1.candidate_nominals.clear();
        self.mate1.placements.clear();
        self.mate2.candidates.clear();
        self.mate2.candidate_nominals.clear();
        self.mate2.placements.clear();
        self.ranked_anchor_placements.clear();

        let located_rows = if anchor_lane == 0 {
            self.append_ranked_block_candidates_for_lane(
                reference,
                read1.len(),
                0,
                maximum_edit_distance,
                &seed_sets[0],
            )?
        } else {
            self.append_ranked_block_candidates_for_lane(
                reference,
                read2.len(),
                1,
                maximum_edit_distance,
                &seed_sets[1],
            )?
        };
        let anchor_metrics = ReadAlignmentMetrics {
            located_rows,
            ..ReadAlignmentMetrics::default()
        };
        let anchor_metrics = if anchor_lane == 0 {
            let (_, metrics) = self.mate1.verify_candidates_with_budget(
                reference,
                read1,
                anchor_metrics,
                maximum_edit_distance,
            )?;
            self.ranked_anchor_placements
                .extend_from_slice(&self.mate1.placements);
            metrics
        } else {
            let (_, metrics) = self.mate2.verify_candidates_with_budget(
                reference,
                read2,
                anchor_metrics,
                maximum_edit_distance,
            )?;
            self.ranked_anchor_placements
                .extend_from_slice(&self.mate2.placements);
            metrics
        };

        let (mate1_metrics, mate2_metrics) = if self.ranked_anchor_placements.is_empty() {
            if anchor_lane == 0 {
                (anchor_metrics, ReadAlignmentMetrics::default())
            } else {
                (ReadAlignmentMetrics::default(), anchor_metrics)
            }
        } else if anchor_lane == 0 {
            let partner_metrics = rescue_from_ranked_anchor_windows(
                &mut self.mate2,
                &mut self.rescue_windows,
                reference,
                read2,
                &self.ranked_anchor_placements,
                false,
                maximum_template_span,
                maximum_edit_distance,
            )?;
            (anchor_metrics, partner_metrics)
        } else {
            let partner_metrics = rescue_from_ranked_anchor_windows(
                &mut self.mate1,
                &mut self.rescue_windows,
                reference,
                read1,
                &self.ranked_anchor_placements,
                true,
                maximum_template_span,
                maximum_edit_distance,
            )?;
            (partner_metrics, anchor_metrics)
        };
        let mut completed = self.finish_directional_pair_combined(
            reference,
            read1,
            read2,
            maximum_edit_distance,
            minimum_template_span,
            maximum_template_span,
            semi_global,
            mate1_metrics,
            mate2_metrics,
            !self.ranked_anchor_placements.is_empty(),
        )?;
        if !anchor_frontier_complete
            && let Some(certified) = self.certify_ranked_pair_frontier(
                reference,
                read1,
                read2,
                projected,
                maximum_edit_distance,
                minimum_template_span,
                maximum_template_span,
                semi_global,
                SENSITIVE_RANKED_BLOCK_HITS.saturating_mul(2),
                completed,
            )?
        {
            return Ok(Some(certified));
        }
        conservatively_mark_incomplete_frontier(&mut completed, anchor_frontier_complete);
        Ok(Some(completed))
    }

    /// Attempts a complete uniqueness proof at the score already discovered by
    /// an incomplete maximum-distance frontier.
    ///
    /// Under strict scoring, a pair with score `k` can disrupt at most `k`
    /// query blocks. Under semi-global scoring every retained mismatch costs
    /// seven and every clipped base has the configured sensitive
    /// penalty, so dividing by the smaller event penalty is a conservative
    /// bound. A `k + 1` partition therefore covers every pair
    /// that could tie or beat the selected score while using longer, much less
    /// repetitive exact blocks than the original distance-five search.
    // Certification replays the complete ranked frontier while preserving the
    // originally selected pair for conservative fallback.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn certify_ranked_pair_frontier(
        &mut self,
        reference: &ReferenceIndex,
        read1: &[Base],
        read2: &[Base],
        projected: [&[ProjectedBase]; 2],
        maximum_edit_distance: u8,
        minimum_template_span: u64,
        maximum_template_span: u64,
        semi_global: bool,
        maximum_combined_proof_hits: u64,
        original: (PairMappingStatus, PairAlignmentMetrics, Option<u8>),
    ) -> Result<Option<(PairMappingStatus, PairAlignmentMetrics, Option<u8>)>, AlignmentError> {
        if !matches!(original.0, PairMappingStatus::Unique) {
            return Ok(None);
        }
        let Some(original_best) = self.best_pairs.first().copied() else {
            return Ok(None);
        };
        let proof_budget = if original.1.semi_global_attempted {
            original_best.score() / SENSITIVE_MIN_EVENT_PENALTY
        } else {
            original_best.score()
        };
        if proof_budget >= maximum_edit_distance {
            return Ok(None);
        }

        let mut seed_sets = [[None; SENSITIVE_PROOF_BLOCKS]; 2];
        let first_selection = Self::collect_ranked_block_seeds_for_lane(
            reference,
            read1,
            projected[0],
            0,
            proof_budget,
            SENSITIVE_RANKED_BLOCK_HITS,
            &mut seed_sets[0],
        )?;
        let second_selection = Self::collect_ranked_block_seeds_for_lane(
            reference,
            read2,
            projected[1],
            1,
            proof_budget,
            SENSITIVE_RANKED_BLOCK_HITS,
            &mut seed_sets[1],
        )?;
        let Some((first_selection, second_selection)) = first_selection.zip(second_selection)
        else {
            return Ok(None);
        };
        if !first_selection.complete || !second_selection.complete {
            return Ok(None);
        }
        if first_selection
            .retained_hits
            .saturating_add(second_selection.retained_hits)
            > maximum_combined_proof_hits
        {
            return Ok(None);
        }

        // The proof partition is complete for every per-mate placement whose
        // distance can contribute to a pair tying or beating `original_best`.
        // Verification still uses the paired-end maximum-distance budget so
        // candidates outside the proof partition retain ordinary semantics.
        let certification_budget = maximum_edit_distance;

        self.best_pairs.clear();
        self.mate1.candidates.clear();
        self.mate1.candidate_nominals.clear();
        self.mate1.placements.clear();
        self.mate2.candidates.clear();
        self.mate2.candidate_nominals.clear();
        self.mate2.placements.clear();
        let mate1_rows = self.append_ranked_block_candidates_for_lane(
            reference,
            read1.len(),
            0,
            proof_budget,
            &seed_sets[0],
        )?;
        let mate2_rows = self.append_ranked_block_candidates_for_lane(
            reference,
            read2.len(),
            1,
            proof_budget,
            &seed_sets[1],
        )?;
        sort_nominal_candidates(&mut self.mate1.candidate_nominals);
        sort_nominal_candidates(&mut self.mate2.candidate_nominals);
        retain_nominal_pair_geometry(
            &mut self.mate1.candidate_nominals,
            &mut self.mate2.candidate_nominals,
            read1.len(),
            read2.len(),
            maximum_template_span,
            certification_budget,
        );
        if self.mate1.candidate_nominals.is_empty() || self.mate2.candidate_nominals.is_empty() {
            self.best_pairs.push(original_best);
            return Ok(None);
        }

        let mate1_seed_metrics = ReadAlignmentMetrics {
            located_rows: mate1_rows,
            ..ReadAlignmentMetrics::default()
        };
        let mate2_seed_metrics = ReadAlignmentMetrics {
            located_rows: mate2_rows,
            ..ReadAlignmentMetrics::default()
        };
        let (_, mate1_metrics) = self.mate1.verify_sorted_candidates_with_budget(
            reference,
            read1,
            mate1_seed_metrics,
            certification_budget,
        )?;
        let (_, mate2_metrics) = self.mate2.verify_sorted_candidates_with_budget(
            reference,
            read2,
            mate2_seed_metrics,
            certification_budget,
        )?;
        let certified = self.finish_directional_pair_combined(
            reference,
            read1,
            read2,
            certification_budget,
            minimum_template_span,
            maximum_template_span,
            semi_global,
            mate1_metrics,
            mate2_metrics,
            false,
        )?;
        if matches!(certified.0, PairMappingStatus::Unmapped)
            || self
                .best_pairs
                .first()
                .is_none_or(|best| best.score() > original_best.score())
        {
            self.best_pairs.clear();
            self.best_pairs.push(original_best);
            return Ok(None);
        }
        Ok(Some(certified))
    }

    // Pair selection, optional rescoring, and confidence aggregation consume
    // the same workspace state and therefore remain one finishing pass.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn finish_directional_pair_combined(
        &mut self,
        reference: &ReferenceIndex,
        read1: &[Base],
        read2: &[Base],
        maximum_edit_distance: u8,
        minimum_template_span: u64,
        maximum_template_span: u64,
        semi_global: bool,
        mut mate1_metrics: ReadAlignmentMetrics,
        mut mate2_metrics: ReadAlignmentMetrics,
        window_rescue_attempted: bool,
    ) -> Result<(PairMappingStatus, PairAlignmentMetrics, Option<u8>), AlignmentError> {
        let mut selection = if self.prefer_minimum_net_gap {
            select_best_pair_origins_with_endpoint_policy(
                &self.mate1.placements,
                &self.mate2.placements,
                [read1, read2],
                maximum_edit_distance,
                minimum_template_span,
                maximum_template_span,
                false,
                &mut self.origin_pair_evidence,
                &mut self.best_pairs,
            )
        } else {
            select_best_pairs(
                &self.mate1.placements,
                &self.mate2.placements,
                maximum_edit_distance,
                minimum_template_span,
                maximum_template_span,
                &mut self.best_pairs,
            )
        };
        collapse_equivalent_pair_origins(&mut self.best_pairs, read1.len(), read2.len());
        if self.prefer_minimum_net_gap && self.best_pairs.len() > 1 {
            prefer_minimum_net_gap_representative(&mut self.best_pairs, read1.len(), read2.len());
        }
        let semi_global_attempted = semi_global
            && (self.best_pairs.len() != 1 || self.best_pairs[0].distance() != 0)
            && !self.mate1.candidate_nominals.is_empty()
            && !self.mate2.candidate_nominals.is_empty();
        if semi_global_attempted {
            append_ungapped_semi_global_placements(
                &mut self.mate1,
                reference,
                read1,
                maximum_edit_distance,
                self.semi_global_clip_penalty,
            );
            mate1_metrics.verified_placements =
                u64::try_from(self.mate1.placements.len()).unwrap_or(u64::MAX);
            append_ungapped_semi_global_placements(
                &mut self.mate2,
                reference,
                read2,
                maximum_edit_distance,
                self.semi_global_clip_penalty,
            );
            mate2_metrics.verified_placements =
                u64::try_from(self.mate2.placements.len()).unwrap_or(u64::MAX);
            // The pair join uses partition points over spatially sorted mate-2
            // placements. Appending a second, individually ordered frontier
            // does not preserve that global order.
            self.mate1
                .placements
                .sort_unstable_by_key(|placement| spatial_key(*placement));
            self.mate2
                .placements
                .sort_unstable_by_key(|placement| spatial_key(*placement));
            selection = if self.prefer_minimum_net_gap {
                select_best_pair_origins_with_endpoint_policy(
                    &self.mate1.placements,
                    &self.mate2.placements,
                    [read1, read2],
                    maximum_edit_distance,
                    minimum_template_span,
                    maximum_template_span,
                    true,
                    &mut self.origin_pair_evidence,
                    &mut self.best_pairs,
                )
            } else {
                select_best_pairs_with_fallback_score(
                    &self.mate1.placements,
                    &self.mate2.placements,
                    maximum_edit_distance,
                    minimum_template_span,
                    maximum_template_span,
                    &mut self.best_pairs,
                )
            };
            collapse_equivalent_pair_origins(&mut self.best_pairs, read1.len(), read2.len());
            if self.prefer_minimum_net_gap && self.best_pairs.len() > 1 {
                prefer_minimum_net_gap_representative(
                    &mut self.best_pairs,
                    read1.len(),
                    read2.len(),
                );
            }
        }
        let exact_retained_alternative = if semi_global_attempted && self.best_pairs.len() == 1 {
            self.exact_retained_pair_has_alternative(
                reference,
                read1,
                read2,
                self.best_pairs[0],
                minimum_template_span,
                maximum_template_span,
            )?
        } else {
            false
        };
        let class = match (self.best_pairs.len(), exact_retained_alternative) {
            (0, _) => PairMappingStatus::Unmapped,
            (1, false) => PairMappingStatus::Unique,
            _ => PairMappingStatus::Ambiguous,
        };
        Ok((
            class,
            PairAlignmentMetrics {
                mate1: mate1_metrics,
                mate2: mate2_metrics,
                compatible_pairs: selection.compatible_pairs,
                best_pair_placements: if exact_retained_alternative {
                    2
                } else {
                    u64::try_from(self.best_pairs.len()).unwrap_or(u64::MAX)
                },
                window_rescue_attempted,
                semi_global_attempted,
                best_pair_score: selection.best_pair_score,
                second_best_pair_score: selection.second_best_pair_score,
                near_best_pairings: selection
                    .near_best_pairings
                    .max(u64::from(exact_retained_alternative)),
                mapq_compatible_pairs: selection.mapq_compatible_pairs,
                mapq_best_pair_score: selection.mapq_best_pair_score,
                mapq_second_best_pair_score: selection.mapq_second_best_pair_score,
                mapq_near_best_pairings: selection
                    .mapq_near_best_pairings
                    .max(u64::from(exact_retained_alternative)),
                frontier_complete: true,
            },
            selection.second_best_distance,
        ))
    }

    // Both retained mates must be re-enumerated under one exact-origin proof;
    // splitting the scan would duplicate its shared anchor state.
    #[allow(clippy::too_many_lines)]
    fn exact_retained_pair_has_alternative(
        &mut self,
        reference: &ReferenceIndex,
        read1: &[Base],
        read2: &[Base],
        selected: PairedPlacement,
        minimum_template_span: u64,
        maximum_template_span: u64,
    ) -> Result<bool, AlignmentError> {
        let placements = [selected.mate1(), selected.mate2()];
        if placements.iter().any(|placement| placement.distance() != 0)
            || !placements[0].is_soft_clipped(read1.len())
                && !placements[1].is_soft_clipped(read2.len())
        {
            return Ok(false);
        }
        let ranges = [
            placements[0].retained_query_interval(read1.len()),
            placements[1].retained_query_interval(read2.len()),
        ];
        let retained = [&read1[ranges[0].clone()], &read2[ranges[1].clone()]];
        let mut first_projected = [SearchBase::A; MAX_READ_BASES];
        let mut second_projected = [SearchBase::A; MAX_READ_BASES];
        prepare_combined_search_projection(
            retained[0],
            ConversionPass::Original,
            &mut first_projected,
        )?;
        prepare_combined_search_projection(
            retained[1],
            ConversionPass::Complementary,
            &mut second_projected,
        )?;
        let projected = [
            &first_projected[..retained[0].len()],
            &second_projected[..retained[1].len()],
        ];
        let Some(first_seed) = reference
            .combined_exact_seed(projected[0])
            .map_err(|_| AlignmentError::CombinedIndex)?
        else {
            return Ok(true);
        };
        let Some(second_seed) = reference
            .combined_exact_seed(projected[1])
            .map_err(|_| AlignmentError::CombinedIndex)?
        else {
            return Ok(true);
        };
        let seeds = [first_seed, second_seed];
        for lane in 0..2 {
            if seeds[lane].matched_bases()
                != u64::try_from(retained[lane].len()).expect("bounded retained length fits u64")
            {
                return Err(AlignmentError::CombinedIndex);
            }
        }
        let hits = seeds.map(CombinedSeedMatches::exact_hit_count);
        if hits == [1, 1] {
            return Ok(false);
        }
        let anchor_lane = usize::from(hits[1] < hits[0]);
        if hits[anchor_lane] > SEMI_GLOBAL_MAX_EXACT_ANCHOR_HITS {
            return Ok(true);
        }
        let other_lane = 1 - anchor_lane;
        self.exact_anchor_candidates.clear();
        reference
            .visit_combined_seed(
                seeds[anchor_lane],
                0,
                u64::try_from(retained[anchor_lane].len())
                    .expect("bounded retained length fits u64"),
                &mut |hit| {
                    if let Some(candidate) = relabel_exact_retained_hit(hit, anchor_lane) {
                        self.exact_anchor_candidates.push(candidate);
                    }
                    true
                },
            )
            .map_err(|_| AlignmentError::CombinedIndex)?;
        if self.exact_anchor_candidates.is_empty() {
            return Ok(true);
        }

        let selected_origin = pair_origin_key(selected, read1.len(), read2.len());
        let anchors = &self.exact_anchor_candidates;
        let mut alternative = false;
        reference
            .visit_combined_seed(
                seeds[other_lane],
                0,
                u64::try_from(retained[other_lane].len())
                    .expect("bounded retained length fits u64"),
                &mut |hit| {
                    let Some(other_candidate) = relabel_exact_retained_hit(hit, other_lane) else {
                        return true;
                    };
                    let Some(other) = exact_retained_placement(
                        other_candidate,
                        placements[other_lane],
                        retained[other_lane].len(),
                    ) else {
                        return true;
                    };
                    for &anchor_candidate in anchors {
                        let Some(anchor) = exact_retained_placement(
                            anchor_candidate,
                            placements[anchor_lane],
                            retained[anchor_lane].len(),
                        ) else {
                            continue;
                        };
                        let pair = if anchor_lane == 0 {
                            exact_compatible_pair(
                                anchor,
                                other,
                                minimum_template_span,
                                maximum_template_span,
                            )
                        } else {
                            exact_compatible_pair(
                                other,
                                anchor,
                                minimum_template_span,
                                maximum_template_span,
                            )
                        };
                        if pair.is_some_and(|pair| {
                            pair_origin_key(pair, read1.len(), read2.len()) != selected_origin
                        }) {
                            alternative = true;
                            return false;
                        }
                    }
                    true
                },
            )
            .map_err(|_| AlignmentError::CombinedIndex)?;
        Ok(alternative)
    }

    pub(super) fn should_affine_rescore(
        &self,
        class: PairMappingStatus,
        metrics: PairAlignmentMetrics,
        read1_len: usize,
        read2_len: usize,
    ) -> bool {
        if !matches!(class, PairMappingStatus::Ambiguous) || self.best_pairs.is_empty() {
            return false;
        }
        metrics.compatible_pairs <= 64
            && self.best_pairs.iter().any(|pair| {
                pair.mate1().is_soft_clipped(read1_len)
                    || pair.mate2().is_soft_clipped(read2_len)
                    || placement_net_gap_bases(pair.mate1(), read1_len) != 0
                    || placement_net_gap_bases(pair.mate2(), read2_len) != 0
            })
    }

    #[allow(clippy::too_many_arguments)]
    fn affine_rescore_directional_pair(
        &mut self,
        reference: &ReferenceIndex,
        read1: &[Base],
        read2: &[Base],
        original_class: PairMappingStatus,
        maximum_edit_distance: u8,
        minimum_template_span: u64,
        maximum_template_span: u64,
        mut metrics: PairAlignmentMetrics,
    ) -> Result<(PairMappingStatus, PairAlignmentMetrics, Option<u8>), AlignmentError> {
        let prior_best_count = metrics.best_pair_placements;
        let prior_retained_count = u64::try_from(self.best_pairs.len()).unwrap_or(u64::MAX);
        let concealed_or_incomplete =
            !metrics.frontier_complete || prior_best_count > prior_retained_count;

        self.mate1_affine_scores.clear();
        for &placement in &self.mate1.placements {
            self.mate1_affine_scores.push(affine_placement_score(
                reference,
                read1,
                placement,
                self.semi_global_clip_penalty,
                &mut self.affine,
            )?);
        }
        self.mate2_affine_scores.clear();
        for &placement in &self.mate2.placements {
            self.mate2_affine_scores.push(affine_placement_score(
                reference,
                read2,
                placement,
                self.semi_global_clip_penalty,
                &mut self.affine,
            )?);
        }

        let selection = select_best_pair_origins_with_affine_score(
            &self.mate1.placements,
            &self.mate1_affine_scores,
            &self.mate2.placements,
            &self.mate2_affine_scores,
            [read1, read2],
            maximum_edit_distance,
            minimum_template_span,
            maximum_template_span,
            &mut self.origin_pair_evidence,
            &mut self.best_pairs,
        );
        collapse_equivalent_pair_origins(&mut self.best_pairs, read1.len(), read2.len());
        if self.best_pairs.len() > 1 {
            prefer_minimum_net_gap_representative(&mut self.best_pairs, read1.len(), read2.len());
        }
        let insufficient_affine_separation = matches!(original_class, PairMappingStatus::Ambiguous)
            && selection
                .best_pair_score
                .zip(selection.second_best_pair_score)
                .is_some_and(|(best, second)| {
                    best.saturating_sub(second)
                        < BWA_MATCH_SCORE.saturating_add(BWA_MISMATCH_PENALTY)
                });
        let must_remain_ambiguous = concealed_or_incomplete || insufficient_affine_separation;
        let class = match (self.best_pairs.len(), must_remain_ambiguous) {
            (0, _) => PairMappingStatus::Unmapped,
            (1, false) => PairMappingStatus::Unique,
            _ => PairMappingStatus::Ambiguous,
        };
        metrics.compatible_pairs = selection.compatible_pairs;
        metrics.best_pair_placements = if must_remain_ambiguous {
            u64::try_from(self.best_pairs.len())
                .unwrap_or(u64::MAX)
                .max(2)
        } else {
            u64::try_from(self.best_pairs.len()).unwrap_or(u64::MAX)
        };
        metrics.best_pair_score = selection.best_pair_score;
        metrics.second_best_pair_score = selection.second_best_pair_score;
        metrics.near_best_pairings = selection
            .near_best_pairings
            .max(u64::from(must_remain_ambiguous));
        metrics.mapq_compatible_pairs = selection.mapq_compatible_pairs;
        metrics.mapq_best_pair_score = selection.mapq_best_pair_score;
        metrics.mapq_second_best_pair_score = selection.mapq_second_best_pair_score;
        metrics.mapq_near_best_pairings = selection
            .mapq_near_best_pairings
            .max(u64::from(must_remain_ambiguous));
        Ok((class, metrics, None))
    }

    #[must_use]
    pub fn best_pairs(&self) -> &[PairedPlacement] {
        &self.best_pairs
    }
}

pub(super) const fn empty_pair_metrics() -> PairAlignmentMetrics {
    PairAlignmentMetrics {
        mate1: ReadAlignmentMetrics {
            located_rows: 0,
            emitted_candidate_starts: 0,
            distinct_candidate_starts: 0,
            verified_placements: 0,
        },
        mate2: ReadAlignmentMetrics {
            located_rows: 0,
            emitted_candidate_starts: 0,
            distinct_candidate_starts: 0,
            verified_placements: 0,
        },
        compatible_pairs: 0,
        best_pair_placements: 0,
        window_rescue_attempted: false,
        semi_global_attempted: false,
        best_pair_score: None,
        second_best_pair_score: None,
        near_best_pairings: 0,
        mapq_compatible_pairs: 0,
        mapq_best_pair_score: None,
        mapq_second_best_pair_score: None,
        mapq_near_best_pairings: 0,
        frontier_complete: false,
    }
}

fn select_combined_window_rescue_anchor(
    first_seeds: [Option<CombinedSeedMatches>; 2],
    mate1: &[ReadCandidate],
    mate2: &[ReadCandidate],
) -> Option<usize> {
    let pools = [mate1, mate2];
    let evidence = core::array::from_fn::<_, 2, _>(|mate| {
        let seed = first_seeds[mate]?;
        if seed.exact_hit_count() != 1 || pools[mate].is_empty() {
            return None;
        }
        let direct = pools[mate]
            .iter()
            .any(|candidate| candidate.proof_mask & DIRECT_SINGLETON_PROOF != 0);
        Some((direct, seed.matched_bases(), pools[mate].len()))
    });
    match (evidence[0], evidence[1]) {
        (None, None) => None,
        (Some(_), None) => Some(0),
        (None, Some(_)) => Some(1),
        (Some(first), Some(second)) => {
            if first.0 != second.0 {
                Some(usize::from(!first.0))
            } else if first.1 != second.1 {
                Some(usize::from(first.1 < second.1))
            } else {
                Some(usize::from(first.2 > second.2))
            }
        }
    }
}

fn nominal_pair_geometry_exists(
    mate1: &[ReadCandidate],
    mate2: &[ReadCandidate],
    read1_len: usize,
    read2_len: usize,
    maximum_span: u64,
    maximum_edit_distance: u8,
) -> bool {
    mate1.iter().any(|&left| {
        nominal_partner_exists(
            left,
            mate2,
            true,
            read1_len,
            read2_len,
            maximum_span,
            maximum_edit_distance,
        )
    })
}

fn retain_nominal_pair_geometry(
    mate1: &mut Vec<ReadCandidate>,
    mate2: &mut Vec<ReadCandidate>,
    read1_len: usize,
    read2_len: usize,
    maximum_span: u64,
    maximum_edit_distance: u8,
) {
    mate1.retain(|left| {
        nominal_partner_exists(
            *left,
            mate2,
            true,
            read1_len,
            read2_len,
            maximum_span,
            maximum_edit_distance,
        )
    });
    mate2.retain(|right| {
        nominal_partner_exists(
            *right,
            mate1,
            false,
            read1_len,
            read2_len,
            maximum_span,
            maximum_edit_distance,
        )
    });
}

fn nominal_partner_exists(
    candidate: ReadCandidate,
    pool: &[ReadCandidate],
    candidate_is_mate1: bool,
    read1_len: usize,
    read2_len: usize,
    maximum_span: u64,
    maximum_edit_distance: u8,
) -> bool {
    let edit_budget = u64::from(maximum_edit_distance);
    let target = match (candidate_is_mate1, candidate.strand()) {
        (true, BisulfiteStrand::OT) => BisulfiteStrand::CTOT,
        (true, BisulfiteStrand::OB) => BisulfiteStrand::CTOB,
        (false, BisulfiteStrand::CTOT) => BisulfiteStrand::OT,
        (false, BisulfiteStrand::CTOB) => BisulfiteStrand::OB,
        _ => return false,
    };
    let lower_start = candidate.start().saturating_sub(maximum_span);
    let upper_start = candidate.start().saturating_add(maximum_span);
    let lower = pool.partition_point(|partner| {
        (partner.strand(), partner.contig_ordinal(), partner.start())
            < (target, candidate.contig_ordinal(), lower_start)
    });
    let upper = pool.partition_point(|partner| {
        (partner.strand(), partner.contig_ordinal(), partner.start())
            <= (target, candidate.contig_ordinal(), upper_start)
    });
    pool[lower..upper].iter().any(|partner| {
        let (left, right) = if candidate_is_mate1 {
            (candidate, *partner)
        } else {
            (*partner, candidate)
        };
        match (left.strand(), right.strand()) {
            (BisulfiteStrand::OT, BisulfiteStrand::CTOT) => {
                left.start()
                    < right
                        .start()
                        .saturating_add(u64::try_from(read2_len).unwrap_or(u64::MAX) + edit_budget)
            }
            (BisulfiteStrand::OB, BisulfiteStrand::CTOB) => {
                right.start()
                    < left
                        .start()
                        .saturating_add(u64::try_from(read1_len).unwrap_or(u64::MAX) + edit_budget)
            }
            _ => false,
        }
    })
}
