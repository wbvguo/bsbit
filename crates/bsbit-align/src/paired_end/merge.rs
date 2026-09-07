//! Cross-pass evidence reduction and conservative frontier reconciliation.

use bsbit_core::alphabet::Base;
use bsbit_index::reference::ReferenceIndex;

use super::evidence::{ExactRetainedPairCheck, PairAlignmentMetrics, PairedBatchResult};
use super::mapq::{
    confidence_score_from_edit_distance,
    effective_mapping_quality as sensitive_effective_mapping_quality, scores_are_near,
    targeted_completion_required,
};
use super::reporting::pair_less_without_reference_order;
use super::result::{PairMappingStatus, PairedPlacement};
use crate::AlignmentError;
use crate::read_mapping::ReadAlignmentMetrics;
use crate::reporting_tie_break::{ReportingTieBreak, pair_origin_hash};

#[allow(clippy::too_many_lines)]
#[cfg(test)]
pub(super) fn merge_non_directional_batch_results(
    original: &PairedBatchResult,
    complementary: &PairedBatchResult,
) -> PairedBatchResult {
    merge_non_directional_batch_results_with_tie_break(None, original, complementary)
        .expect("a disabled tie-break does not access a reference")
}

#[allow(clippy::too_many_lines)]
pub(super) fn merge_non_directional_batch_results_with_tie_break(
    tie_break: Option<(&ReferenceIndex, [&[Base]; 2], ReportingTieBreak)>,
    original: &PairedBatchResult,
    complementary: &PairedBatchResult,
) -> Result<PairedBatchResult, AlignmentError> {
    let original_score = batch_best_score(original);
    let complementary_score = batch_best_score(complementary);
    let (mut selected, other, tied) = match (original_score, complementary_score) {
        (Some(left), Some(right)) if left > right => (*original, complementary, false),
        (Some(left), Some(right)) if right > left => (*complementary, original, false),
        (Some(_), Some(_)) => {
            let select_complementary = if let (
                Some((reference, reads, tie_break)),
                Some(original_pair),
                Some(complementary_pair),
            ) =
                (tie_break, original.best_pair, complementary.best_pair)
            {
                let original_hash = pair_origin_hash(
                    reference,
                    tie_break,
                    original_pair.mate1(),
                    reads[0].len(),
                    original_pair.mate2(),
                    reads[1].len(),
                )?;
                let complementary_hash = pair_origin_hash(
                    reference,
                    tie_break,
                    complementary_pair.mate1(),
                    reads[0].len(),
                    complementary_pair.mate2(),
                    reads[1].len(),
                )?;
                complementary_hash < original_hash
                    || (complementary_hash == original_hash
                        && pair_less_without_reference_order(
                            reference,
                            reads,
                            complementary_pair,
                            original_pair,
                        )?)
            } else {
                false
            };
            if select_complementary {
                (*complementary, original, true)
            } else {
                (*original, complementary, true)
            }
        }
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
        exact_retained_pair_check: merge_exact_retained_pair_checks(
            selected_metrics.exact_retained_pair_check,
            other_metrics.exact_retained_pair_check,
        ),
        resolved_prior_ambiguity: selected_metrics.resolved_prior_ambiguity
            || other_metrics.resolved_prior_ambiguity,
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
                    if scores_are_near(selected_score, other_score) {
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
                        if scores_are_near(selected_score, other_score) {
                            other_metrics.mapq_near_best_pairings.saturating_add(1)
                        } else {
                            0
                        }
                    })
            }),
        ),
        frontier_complete: selected_metrics.frontier_complete && other_metrics.frontier_complete,
        alternative_margin_frontier_complete: selected_metrics.alternative_margin_frontier_complete
            && other_metrics.alternative_margin_frontier_complete,
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
    } else if !selected.metrics.frontier_complete
        && matches!(selected.class, PairMappingStatus::Unique)
    {
        selected.class = PairMappingStatus::Ambiguous;
        selected.metrics.best_pair_placements = selected.metrics.best_pair_placements.max(2);
    }
    Ok(selected)
}

pub(super) fn batch_best_score(result: &PairedBatchResult) -> Option<i16> {
    result.metrics.best_pair_score.or_else(|| {
        result
            .best_pair
            .map(|pair| confidence_score_from_edit_distance(pair.score()))
    })
}

pub(super) const fn merge_exact_retained_pair_checks(
    left: ExactRetainedPairCheck,
    right: ExactRetainedPairCheck,
) -> ExactRetainedPairCheck {
    use ExactRetainedPairCheck::{
        AlternativeFound, InconclusiveAnchorLimit, InconclusiveEmptyAnchorSet,
        InconclusiveMissingSeed, NoAlternative, NotRequired,
    };
    match (left, right) {
        (AlternativeFound, _) | (_, AlternativeFound) => AlternativeFound,
        (InconclusiveAnchorLimit, _) | (_, InconclusiveAnchorLimit) => InconclusiveAnchorLimit,
        (InconclusiveMissingSeed, _) | (_, InconclusiveMissingSeed) => InconclusiveMissingSeed,
        (InconclusiveEmptyAnchorSet, _) | (_, InconclusiveEmptyAnchorSet) => {
            InconclusiveEmptyAnchorSet
        }
        (NoAlternative, _) | (_, NoAlternative) => NoAlternative,
        (NotRequired, NotRequired) => NotRequired,
    }
}

pub(super) const fn maximum_optional_score(left: Option<i16>, right: Option<i16>) -> Option<i16> {
    match (left, right) {
        (Some(left), Some(right)) => Some(if left > right { left } else { right }),
        (Some(value), None) | (None, Some(value)) => Some(value),
        (None, None) => None,
    }
}

pub(super) const fn minimum_optional_distance(left: Option<u8>, right: Option<u8>) -> Option<u8> {
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

pub(super) fn retain_completed_runner_up_evidence(
    mut incumbent: PairAlignmentMetrics,
    completed: PairAlignmentMetrics,
) -> PairAlignmentMetrics {
    // A completed disjoint-block search can return a worse best placement
    // when the fast adaptive incumbent was not rediscovered. Keep only
    // adverse evidence while preserving the incumbent coordinate.
    incumbent.mate1 = merge_read_metrics(incumbent.mate1, completed.mate1);
    incumbent.mate2 = merge_read_metrics(incumbent.mate2, completed.mate2);
    incumbent.compatible_pairs = incumbent
        .compatible_pairs
        .saturating_add(completed.compatible_pairs);
    incumbent.mapq_compatible_pairs = incumbent
        .mapq_compatible_pairs
        .saturating_add(completed.mapq_compatible_pairs);
    incumbent.second_best_pair_score =
        maximum_optional_score(incumbent.second_best_pair_score, completed.best_pair_score);
    incumbent.mapq_second_best_pair_score = maximum_optional_score(
        incumbent.mapq_second_best_pair_score,
        completed.mapq_best_pair_score,
    );
    if incumbent
        .best_pair_score
        .zip(completed.best_pair_score)
        .is_some_and(|(best, runner_up)| scores_are_near(best, runner_up))
    {
        incumbent.near_best_pairings = incumbent
            .near_best_pairings
            .saturating_add(completed.near_best_pairings)
            .saturating_add(1);
    }
    if incumbent
        .mapq_best_pair_score
        .zip(completed.mapq_best_pair_score)
        .is_some_and(|(best, runner_up)| scores_are_near(best, runner_up))
    {
        incumbent.mapq_near_best_pairings = incumbent
            .mapq_near_best_pairings
            .saturating_add(completed.mapq_near_best_pairings)
            .saturating_add(1);
    }
    incumbent.window_rescue_attempted |= completed.window_rescue_attempted;
    incumbent.semi_global_attempted |= completed.semi_global_attempted;
    incumbent.resolved_prior_ambiguity |= completed.resolved_prior_ambiguity;
    incumbent.frontier_complete = completed.frontier_complete;
    incumbent.alternative_margin_frontier_complete &=
        completed.alternative_margin_frontier_complete;
    incumbent
}

pub(super) const fn sensitive_unique_frontier_completion_required(
    class: PairMappingStatus,
) -> bool {
    // A singleton found by one adaptive seed is not a uniqueness certificate.
    matches!(class, PairMappingStatus::Unique)
}

pub(super) fn sensitive_targeted_semi_global_required(
    class: PairMappingStatus,
    metrics: PairAlignmentMetrics,
) -> bool {
    !matches!(class, PairMappingStatus::Unmapped)
        && (!metrics.frontier_complete
            || targeted_completion_required(sensitive_effective_mapping_quality(class, metrics)))
}

pub(super) fn restore_rejected_targeted_frontier(
    best_pairs: &mut Vec<PairedPlacement>,
    mut incumbent: Vec<PairedPlacement>,
    preserve_ambiguous_tie: bool,
) {
    if !preserve_ambiguous_tie {
        incumbent.truncate(1);
    }
    *best_pairs = incumbent;
}

pub(super) fn conservatively_mark_incomplete_frontier(
    result: &mut (PairMappingStatus, PairAlignmentMetrics, Option<u8>),
    complete: bool,
) {
    result.1.frontier_complete = complete;
    if !complete {
        result.1.alternative_margin_frontier_complete = false;
    }
    if !complete && matches!(result.0, PairMappingStatus::Unique) {
        result.0 = PairMappingStatus::Ambiguous;
        result.1.best_pair_placements = result.1.best_pair_placements.max(2);
        result.2 = None;
    }
}
