//! Paired-end mapping-quality evidence and qualified confidence policy.
//!
//! Pair search supplies complete-frontier and score-gap evidence. This module
//! converts that evidence to MAPQ without owning discovery or serialization.

use std::ops::Range;

use crate::placement::placement_net_gap_bases;

use super::{
    PairAlignmentMetrics, PairMappingStatus, PairedBatchResult, PairedPlacement, PairedSearchMode,
};

const AMBIGUITY_Q10_COMPLETE_SPARSE_ORIGINAL_BEST_COUNT_MAX: u64 = 3;
const AMBIGUITY_Q10_COMPLETE_SPARSE_CANDIDATE_ROWS_MIN_MAX: u64 = 3;
const AMBIGUITY_Q10_COMPLETE_SPARSE_CANDIDATE_NEAR_BEST_MAX: u64 = 1;
const AMBIGUITY_Q10_COMPLETE_SPARSE_ORIGINAL_PAIR_SCORE_MIN: u64 = 7;
const AMBIGUITY_Q10_COMPLETE_SPARSE_CANDIDATE_PAIR_DISTANCE_MIN: u64 = 3;
const AMBIGUITY_Q10_INCOMPLETE_BEST_COUNT: u64 = 2;
const AMBIGUITY_Q10_INCOMPLETE_NEAR_BEST_MAX: u64 = 0;
const AMBIGUITY_Q10_INCOMPLETE_ROWS_MAX_MIN: u64 = 4;
const AMBIGUITY_Q10_INCOMPLETE_PAIR_DISTANCE_MIN: u64 = 1;
const AMBIGUITY_Q10_INCOMPLETE_PAIR_DISTANCE_MAX: u64 = 3;

/// Candidate-row pressure above which a completed sensitive result
/// lacks enough uniqueness evidence for MAPQ 20 or greater.
pub const SENSITIVE_MAPQ_REPEAT_RISK_ROWS: u64 = 384;

pub(crate) const PARSIMONY_MAX_LOCATED_ROWS: u64 = 126;
const PARSIMONY_MAX_VERIFIED_PLACEMENTS: u64 = 7;
pub(crate) const PARSIMONY_REQUIRED_SCORE_GAP: i16 = 3;
const PARSIMONY_MAX_PAIR_SCORE: u8 = 38;

/// Computes the pair-level BWA-style score-gap mapping quality used by the
/// paired-end aligner.
///
/// Reporting-only clipping and repeat-risk caps are deliberately separate so
/// search and serialization can consume the same evidence calculation.
#[must_use]
// This implements the established BWA logarithmic repeat penalty. Evidence
// counts are converted to f64 only for ln(), then clamped to the u8 MAPQ range.
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss
)]
pub fn bwa_pair_mapping_quality_from_evidence(
    class: PairMappingStatus,
    frontier_complete: bool,
    best: Option<i16>,
    second_best: Option<i16>,
    near_best_pairings: u64,
) -> u8 {
    if !matches!(class, PairMappingStatus::Unique) || !frontier_complete {
        return 0;
    }
    let Some(best) = best else {
        return 0;
    };
    let raw = second_best.map_or(60_i32, |second| {
        let score_gap = i32::from(best.saturating_sub(second));
        (602_i32.saturating_mul(score_gap).saturating_add(50)) / 100
    });
    let repeat_penalty = if near_best_pairings == 0 {
        0
    } else {
        (4.343_f64 * (near_best_pairings as f64 + 1.0).ln() + 0.499) as i32
    };
    raw.saturating_sub(repeat_penalty).clamp(0, 60) as u8
}

pub(crate) fn incomplete_sparse_completion_required(
    class: PairMappingStatus,
    metrics: PairAlignmentMetrics,
    pair_distance: Option<u8>,
) -> bool {
    matches!(class, PairMappingStatus::Ambiguous)
        && !metrics.frontier_complete
        && metrics.best_pair_placements == AMBIGUITY_Q10_INCOMPLETE_BEST_COUNT
        && metrics.near_best_pairings == AMBIGUITY_Q10_INCOMPLETE_NEAR_BEST_MAX
        && metrics.second_best_pair_score.is_none()
        && metrics.mate1.located_rows.max(metrics.mate2.located_rows)
            >= AMBIGUITY_Q10_INCOMPLETE_ROWS_MAX_MIN
        && pair_distance.is_some_and(|distance| {
            let distance = u64::from(distance);
            (AMBIGUITY_Q10_INCOMPLETE_PAIR_DISTANCE_MIN
                ..=AMBIGUITY_Q10_INCOMPLETE_PAIR_DISTANCE_MAX)
                .contains(&distance)
        })
}

pub(crate) fn effective_mapping_quality(
    class: PairMappingStatus,
    metrics: PairAlignmentMetrics,
) -> u8 {
    let adjusted_mapq = bwa_pair_mapping_quality_from_evidence(
        class,
        metrics.frontier_complete,
        metrics.best_pair_score,
        metrics.second_best_pair_score,
        metrics.near_best_pairings,
    );
    if metrics.window_rescue_attempted
        || metrics.mate1.located_rows.max(metrics.mate2.located_rows)
            >= SENSITIVE_MAPQ_REPEAT_RISK_ROWS
    {
        adjusted_mapq.min(19)
    } else {
        adjusted_mapq
    }
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn stable_rescue_q20_certified(
    original_class: PairMappingStatus,
    original: PairAlignmentMetrics,
    original_confidence: u8,
    candidate_class: PairMappingStatus,
    candidate: PairAlignmentMetrics,
    candidate_confidence: u8,
    same_origin: bool,
) -> bool {
    if !same_origin
        || !matches!(original_class, PairMappingStatus::Unique)
        || !matches!(candidate_class, PairMappingStatus::Unique)
        || !original.frontier_complete
        || !candidate.frontier_complete
        || original_confidence != 19
        || candidate_confidence != 19
    {
        return false;
    }
    let rows = [candidate.mate1.located_rows, candidate.mate2.located_rows];
    if rows[0].min(rows[1]) != 0 || rows[0].max(rows[1]) < 5 {
        return false;
    }
    let verified = candidate
        .mate1
        .verified_placements
        .saturating_add(candidate.mate2.verified_placements);
    if (4..=54).contains(&verified) {
        return true;
    }
    verified <= 3
        && bwa_pair_mapping_quality_from_evidence(
            candidate_class,
            candidate.frontier_complete,
            candidate.best_pair_score,
            candidate.second_best_pair_score,
            candidate.near_best_pairings,
        ) == 60
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn ambiguity_q10_certified(
    original_class: PairMappingStatus,
    original: PairAlignmentMetrics,
    original_pair_score: Option<u8>,
    original_pair_distance: Option<u8>,
    candidate_class: PairMappingStatus,
    candidate: PairAlignmentMetrics,
    candidate_pair_distance: Option<u8>,
    candidate_net_gap_profile: (Option<u64>, Option<u64>, u64),
    same_origin: bool,
) -> bool {
    let unique_minimum_net_gap = original.frontier_complete
        && matches!(candidate_class, PairMappingStatus::Ambiguous)
        && candidate.frontier_complete
        && candidate_net_gap_profile.2 == 1
        && candidate_net_gap_profile.1.is_some();
    let candidate_rows_min = candidate
        .mate1
        .located_rows
        .min(candidate.mate2.located_rows);
    let sparse_distance = same_origin
        && matches!(original_class, PairMappingStatus::Ambiguous)
        && matches!(candidate_class, PairMappingStatus::Ambiguous)
        && original.frontier_complete
        && candidate.frontier_complete
        && original.best_pair_placements <= AMBIGUITY_Q10_COMPLETE_SPARSE_ORIGINAL_BEST_COUNT_MAX
        && candidate_rows_min <= AMBIGUITY_Q10_COMPLETE_SPARSE_CANDIDATE_ROWS_MIN_MAX
        && candidate.near_best_pairings <= AMBIGUITY_Q10_COMPLETE_SPARSE_CANDIDATE_NEAR_BEST_MAX
        && original_pair_score.is_some_and(|score| {
            u64::from(score) >= AMBIGUITY_Q10_COMPLETE_SPARSE_ORIGINAL_PAIR_SCORE_MIN
        })
        && candidate_pair_distance.is_some_and(|distance| {
            u64::from(distance) >= AMBIGUITY_Q10_COMPLETE_SPARSE_CANDIDATE_PAIR_DISTANCE_MIN
        });
    let completed_incomplete_sparse_frontier = same_origin
        && incomplete_sparse_completion_required(original_class, original, original_pair_distance)
        && matches!(candidate_class, PairMappingStatus::Unique)
        && candidate.frontier_complete;
    unique_minimum_net_gap || sparse_distance || completed_incomplete_sparse_frontier
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn two_way_parsimony_q20_certified(
    original_class: PairMappingStatus,
    original: PairAlignmentMetrics,
    original_best: Option<PairedPlacement>,
    candidate_class: PairMappingStatus,
    candidate: PairAlignmentMetrics,
    candidate_best: Option<PairedPlacement>,
    candidate_pair_count: usize,
    candidate_net_gap_profile: (Option<u64>, Option<u64>, u64),
    read1_len: usize,
    read2_len: usize,
    same_origin: bool,
) -> bool {
    if !same_origin
        || !matches!(original_class, PairMappingStatus::Ambiguous)
        || !matches!(candidate_class, PairMappingStatus::Ambiguous)
        || !original.frontier_complete
        || !candidate.frontier_complete
        || original.best_pair_placements != 2
        || candidate.best_pair_placements != 2
        || candidate_pair_count != 2
        || original.window_rescue_attempted
        || original.mate1.located_rows.max(original.mate2.located_rows) > PARSIMONY_MAX_LOCATED_ROWS
    {
        return false;
    }
    if original
        .best_pair_score
        .zip(original.second_best_pair_score)
        .is_none_or(|(best, second)| best.saturating_sub(second) != PARSIMONY_REQUIRED_SCORE_GAP)
    {
        return false;
    }
    let verified = candidate
        .mate1
        .verified_placements
        .saturating_add(candidate.mate2.verified_placements);
    if verified > PARSIMONY_MAX_VERIFIED_PLACEMENTS {
        return false;
    }
    let Some(original_best) = original_best else {
        return false;
    };
    if original_best.score() > PARSIMONY_MAX_PAIR_SCORE {
        return false;
    }
    let Some(candidate_best) = candidate_best else {
        return false;
    };
    if candidate_best.mate1().is_soft_clipped(read1_len)
        || candidate_best.mate2().is_soft_clipped(read2_len)
    {
        return false;
    }
    let (minimum_gap, second_gap, minimum_count) = candidate_net_gap_profile;
    minimum_count == 1
        && minimum_gap
            .zip(second_gap)
            .is_some_and(|(minimum, second)| minimum < second)
        && minimum_gap.is_some_and(|minimum| {
            placement_net_gap_bases(candidate_best.mate1(), read1_len)
                .saturating_add(placement_net_gap_bases(candidate_best.mate2(), read2_len))
                == minimum
        })
}

#[derive(Clone, Copy, Debug)]
struct MapqPlacementContext<'a> {
    selected: PairedPlacement,
    read_lengths: [usize; 2],
    retained_ranges: [&'a Range<usize>; 2],
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct SensitiveMapqEvidence {
    baseline_mapq: u8,
    raw_mapq: u8,
    reported_ambiguous: bool,
    frontier_complete: bool,
    best_pair_placements: u64,
    compatible_pairs: u64,
    best_score: i16,
    second_best_present: bool,
    score_gap: i32,
    near_best_pairings: u64,
    located_rows_min: u64,
    located_rows_max: u64,
    located_rows_sum: u64,
    emitted_candidate_starts_sum: u64,
    distinct_candidate_starts_sum: u64,
    verified_placements_sum: u64,
    pair_distance: u8,
    pair_score: u8,
    mate_distance_max: u8,
    net_gap_sum: u64,
    clipped_bases_sum: u64,
}

impl SensitiveMapqEvidence {
    fn from_result(
        result: PairedBatchResult,
        context: MapqPlacementContext<'_>,
        baseline_mapq: u8,
        raw_mapq: u8,
        reported_ambiguous: bool,
    ) -> Self {
        let MapqPlacementContext {
            selected,
            read_lengths,
            retained_ranges,
        } = context;
        let metrics = result.metrics();
        let located_rows_min = metrics.mate1.located_rows.min(metrics.mate2.located_rows);
        let located_rows_max = metrics.mate1.located_rows.max(metrics.mate2.located_rows);
        let located_rows_sum = metrics
            .mate1
            .located_rows
            .saturating_add(metrics.mate2.located_rows);
        let emitted_candidate_starts_sum = metrics
            .mate1
            .emitted_candidate_starts
            .saturating_add(metrics.mate2.emitted_candidate_starts);
        let distinct_candidate_starts_sum = metrics
            .mate1
            .distinct_candidate_starts
            .saturating_add(metrics.mate2.distinct_candidate_starts);
        let verified_placements_sum = metrics
            .mate1
            .verified_placements
            .saturating_add(metrics.mate2.verified_placements);
        let score_gap = metrics
            .mapq_best_pair_score
            .zip(metrics.mapq_second_best_pair_score)
            .map_or(99, |(best, second)| i32::from(best) - i32::from(second));
        let placements = [selected.mate1(), selected.mate2()];
        let retained_bases = [
            retained_ranges[0]
                .end
                .saturating_sub(retained_ranges[0].start),
            retained_ranges[1]
                .end
                .saturating_sub(retained_ranges[1].start),
        ];
        let net_gap_sum = placements
            .iter()
            .zip(retained_bases)
            .map(|(placement, retained)| {
                let reference_bases = placement.end().saturating_sub(placement.start());
                reference_bases.abs_diff(
                    u64::try_from(retained).expect("bounded retained read length fits u64"),
                )
            })
            .fold(0_u64, u64::saturating_add);
        let clipped_bases_sum = read_lengths
            .iter()
            .zip(retained_bases)
            .map(|(read_length, retained)| read_length.saturating_sub(retained))
            .map(|clipped| u64::try_from(clipped).expect("bounded read length fits u64"))
            .fold(0_u64, u64::saturating_add);

        Self {
            baseline_mapq,
            raw_mapq,
            reported_ambiguous,
            frontier_complete: metrics.frontier_complete,
            best_pair_placements: metrics.best_pair_placements,
            compatible_pairs: metrics.mapq_compatible_pairs,
            best_score: metrics.mapq_best_pair_score.unwrap_or(i16::MIN),
            second_best_present: metrics.mapq_second_best_pair_score.is_some(),
            score_gap,
            near_best_pairings: metrics.mapq_near_best_pairings,
            located_rows_min,
            located_rows_max,
            located_rows_sum,
            emitted_candidate_starts_sum,
            distinct_candidate_starts_sum,
            verified_placements_sum,
            pair_distance: selected.distance(),
            pair_score: selected.score(),
            mate_distance_max: placements[0].distance().max(placements[1].distance()),
            net_gap_sum,
            clipped_bases_sum,
        }
    }

    fn with_endpoint_evidence(mut self, metrics: PairAlignmentMetrics) -> Self {
        self.compatible_pairs = metrics.compatible_pairs;
        self.best_score = metrics.best_pair_score.unwrap_or(i16::MIN);
        self.second_best_present = metrics.second_best_pair_score.is_some();
        self.score_gap = metrics
            .best_pair_score
            .zip(metrics.second_best_pair_score)
            .map_or(99, |(best, second)| i32::from(best) - i32::from(second));
        self.near_best_pairings = metrics.near_best_pairings;
        self
    }
}

fn sensitive_repeat_risk(
    search_mode: PairedSearchMode,
    window_rescue_attempted: bool,
    mate1_located_rows: u64,
    mate2_located_rows: u64,
) -> bool {
    matches!(search_mode, PairedSearchMode::Sensitive)
        && (window_rescue_attempted
            || mate1_located_rows.max(mate2_located_rows) >= SENSITIVE_MAPQ_REPEAT_RISK_ROWS)
}

fn origin_mapq_evidence_matches_endpoints(metrics: PairAlignmentMetrics) -> bool {
    metrics.mapq_compatible_pairs == metrics.compatible_pairs
        && metrics.mapq_best_pair_score == metrics.best_pair_score
        && metrics.mapq_second_best_pair_score == metrics.second_best_pair_score
        && metrics.mapq_near_best_pairings == metrics.near_best_pairings
}

const fn sensitive_q10_certified(evidence: SensitiveMapqEvidence) -> bool {
    evidence.baseline_mapq < 10
        && evidence.reported_ambiguous
        && evidence.best_pair_placements <= 2
        && ((evidence.verified_placements_sum <= 3
            && evidence.located_rows_sum >= 3
            && evidence.pair_distance <= 4
            && (evidence.located_rows_min >= 3 || evidence.emitted_candidate_starts_sum <= 3))
            || (evidence.verified_placements_sum >= 4
                && evidence.near_best_pairings == 0
                && ((evidence.pair_distance <= 3 && evidence.located_rows_sum >= 144)
                    || (evidence.pair_distance >= 4
                        && evidence.located_rows_sum >= 101
                        && evidence.pair_score <= 4)))
            || (evidence.verified_placements_sum >= 4
                && evidence.near_best_pairings >= 1
                && evidence.emitted_candidate_starts_sum <= 2
                && evidence.pair_score <= 38
                && evidence.net_gap_sum == 0))
}

const fn sensitive_q30_low_certified(evidence: SensitiveMapqEvidence) -> bool {
    evidence.baseline_mapq < 20
        && !evidence.reported_ambiguous
        && evidence.frontier_complete
        && evidence.best_pair_placements <= 2
        && ((evidence.raw_mapq <= 8
            && evidence.verified_placements_sum >= 4
            && evidence.near_best_pairings == 0
            && evidence.best_score >= -13
            && evidence.located_rows_sum >= 144
            && evidence.score_gap >= 11)
            || (evidence.raw_mapq >= 9
                && evidence.pair_distance <= 3
                && evidence.located_rows_max <= 534
                && evidence.emitted_candidate_starts_sum >= 6
                && evidence.compatible_pairs <= 355
                && evidence.located_rows_sum <= 622))
}

const fn sensitive_q40_certified(evidence: SensitiveMapqEvidence) -> bool {
    evidence.baseline_mapq >= 20
        && evidence.baseline_mapq < 30
        && !evidence.reported_ambiguous
        && evidence.frontier_complete
        && evidence.best_pair_placements <= 2
        && evidence.raw_mapq >= 9
        && ((evidence.compatible_pairs <= 2
            && (evidence.raw_mapq >= 25
                || (evidence.raw_mapq <= 24
                    && evidence.located_rows_sum <= 165
                    && evidence.located_rows_min >= 2)))
            || (evidence.compatible_pairs >= 3
                && evidence.compatible_pairs <= 7
                && evidence.clipped_bases_sum <= 1
                && evidence.mate_distance_max <= 2
                && evidence.net_gap_sum == 0)
            || (evidence.compatible_pairs >= 11 && evidence.clipped_bases_sum == 0))
}

const fn apply_sensitive_mapq_policy(evidence: SensitiveMapqEvidence) -> u8 {
    let mut adjusted_mapq = evidence.baseline_mapq;
    if sensitive_q10_certified(evidence) {
        adjusted_mapq = 10;
    }
    if evidence.baseline_mapq >= 20 && evidence.baseline_mapq < 30 {
        adjusted_mapq = 30;
    }
    if sensitive_q30_low_certified(evidence) {
        adjusted_mapq = 30;
    }
    if sensitive_q40_certified(evidence) {
        adjusted_mapq = 40;
    }
    adjusted_mapq
}

const fn apply_origin_grouped_mapq_policy(
    endpoint_mapq: u8,
    grouped_mapq: u8,
    grouped: SensitiveMapqEvidence,
) -> u8 {
    if grouped_mapq <= endpoint_mapq {
        return endpoint_mapq;
    }
    let mut adjusted_mapq = endpoint_mapq;
    if adjusted_mapq < 10
        && grouped_mapq >= 10
        && grouped.pair_distance <= 2
        && grouped.frontier_complete
    {
        adjusted_mapq = 10;
    }
    if adjusted_mapq < 30 && grouped_mapq >= 30 && grouped.score_gap >= 16 {
        adjusted_mapq = 30;
    }
    let q40_certified = grouped.pair_distance == 0
        || (grouped.clipped_bases_sum == 0 && grouped.pair_distance <= 3)
        || (grouped.near_best_pairings == 0
            && grouped.second_best_present
            && grouped.score_gap >= 16);
    if adjusted_mapq < 40 && grouped_mapq >= 40 && q40_certified {
        adjusted_mapq = 40;
    }
    adjusted_mapq
}

const fn sensitive_q40_common_evidence(evidence: SensitiveMapqEvidence) -> bool {
    !evidence.reported_ambiguous && evidence.frontier_complete && evidence.best_pair_placements == 1
}

const fn sensitive_q40_v5_evidence(evidence: SensitiveMapqEvidence) -> bool {
    evidence.distinct_candidate_starts_sum == 0
        || evidence.located_rows_max < 16
        || (evidence.distinct_candidate_starts_sum <= 1
            && evidence.located_rows_max < 128
            && evidence.pair_distance <= 3)
}

const fn apply_sensitive_q40_certificate(mapq: u8, evidence: SensitiveMapqEvidence) -> u8 {
    if mapq < 40 {
        return mapq;
    }
    if sensitive_q40_common_evidence(evidence) && sensitive_q40_v5_evidence(evidence) {
        mapq
    } else {
        30
    }
}

pub(crate) fn paired_mapping_quality(
    result: PairedBatchResult,
    stability_result: Option<PairedBatchResult>,
    class: PairMappingStatus,
    search_mode: PairedSearchMode,
    read_lengths: [usize; 2],
    retained_ranges: [&Range<usize>; 2],
) -> u8 {
    let selected = result
        .best_pair()
        .expect("a reportable paired-end result has a selected pair");
    let reported_ambiguous = matches!(class, PairMappingStatus::Ambiguous);
    let soft_clipped = read_lengths
        .iter()
        .zip(retained_ranges)
        .any(|(length, retained)| retained.start != 0 || retained.end != *length);
    let raw_mapq = result.evidence_mapping_quality();
    let mut adjusted_mapq = if reported_ambiguous {
        0
    } else if result.parsimony_q20_certified() {
        20
    } else {
        raw_mapq
    };
    if let Some(stability) = stability_result {
        adjusted_mapq = adjusted_mapq.min(stability.evidence_mapping_quality());
    }
    let repeat_risk = |candidate: PairedBatchResult| {
        let metrics = candidate.metrics();
        sensitive_repeat_risk(
            search_mode,
            metrics.window_rescue_attempted,
            metrics.mate1.located_rows,
            metrics.mate2.located_rows,
        )
    };
    let result_repeat_risk = repeat_risk(result);
    let uncertified_repeat_risk = result_repeat_risk && !result.repeat_risk_q20_certified()
        || stability_result.is_some_and(|stability| {
            repeat_risk(stability) && !stability.repeat_risk_q20_certified()
        });
    if uncertified_repeat_risk {
        adjusted_mapq = adjusted_mapq.min(19);
    } else if result_repeat_risk && adjusted_mapq >= 19 && result.repeat_risk_q20_certified() {
        adjusted_mapq = 20;
    }
    if soft_clipped {
        adjusted_mapq = adjusted_mapq.min(20);
    }
    if !matches!(search_mode, PairedSearchMode::Sensitive) {
        return adjusted_mapq;
    }
    let context = MapqPlacementContext {
        selected,
        read_lengths,
        retained_ranges,
    };
    let evidence = SensitiveMapqEvidence::from_result(
        result,
        context,
        adjusted_mapq,
        raw_mapq,
        reported_ambiguous,
    );
    let grouped_mapq = apply_sensitive_mapq_policy(evidence);
    let grouped_mapq = if origin_mapq_evidence_matches_endpoints(result.metrics()) {
        grouped_mapq
    } else {
        apply_origin_grouped_mapq_policy(
            apply_sensitive_mapq_policy(evidence.with_endpoint_evidence(result.metrics())),
            grouped_mapq,
            evidence,
        )
    };
    apply_sensitive_q40_certificate(grouped_mapq, evidence)
}

#[cfg(test)]
#[path = "../../tests/whitebox/mapq_policy.rs"]
mod whitebox;
