//! Paired-end mapping-quality evidence and qualified confidence policy.
//!
//! Pair search supplies complete-frontier and score-gap evidence. This module
//! converts that evidence to MAPQ without owning discovery or serialization.

use std::ops::Range;

use crate::mapq_policy::{MAPQ_POLICY, MapqCertificate};

use super::evidence::{PairAlignmentMetrics, PairedBatchResult};
use super::result::PairMappingStatus;

/// Candidate-row pressure above which a completed result requires the
/// sensitive completion pass before high-confidence reporting.
const PAIR_MAPQ_REPEAT_RISK_ROWS: u64 = MAPQ_POLICY.paired.repeat_risk_rows;
/// Largest score deficit counted as a near-best alternative origin.
pub(crate) const PAIR_NEAR_SUBOPTIMAL_SCORE_DELTA: i16 =
    MAPQ_POLICY.paired.near_suboptimal_score_delta;

/// Converts a pair edit distance to the affine-free confidence score used
/// when no richer endpoint score is available.
#[must_use]
pub(crate) fn confidence_score_from_edit_distance(distance: u8) -> i16 {
    -i16::from(distance).saturating_mul(MAPQ_POLICY.paired.confidence_mismatch_penalty)
}

/// Reports whether one lower-scoring origin belongs to the near-best MAPQ
/// evidence band.
#[must_use]
pub(crate) const fn scores_are_near(best: i16, runner_up: i16) -> bool {
    best.saturating_sub(runner_up) <= PAIR_NEAR_SUBOPTIMAL_SCORE_DELTA
}

/// Additional edit events enumerated to certify the high-confidence pair
/// margin.
#[must_use]
pub(crate) const fn confidence_proof_extra_edits() -> u8 {
    MAPQ_POLICY.paired.confidence_proof_extra_edits
}

/// Reports whether current pair confidence warrants targeted endpoint
/// completion.
#[must_use]
pub(crate) const fn targeted_completion_required(mapping_quality: u8) -> bool {
    mapping_quality < MAPQ_POLICY.paired.targeted_completion_trigger
}

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
pub(super) fn bwa_pair_mapping_quality_from_evidence(
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
    let raw = second_best.map_or(i32::from(MAPQ_POLICY.paired.maximum_mapq), |second| {
        let score_gap = i32::from(best.saturating_sub(second));
        (MAPQ_POLICY
            .paired
            .score_gap_scale_hundredths
            .saturating_mul(score_gap)
            .saturating_add(MAPQ_POLICY.paired.score_gap_scale_hundredths / 20))
            / 100
    });
    // The best runner-up is already represented by `score_gap`.  Count only
    // additional near-best origins in the logarithmic multiplicity penalty;
    // including the first runner-up here would double-count the same evidence.
    let additional_near_best = near_best_pairings.saturating_sub(1);
    let repeat_penalty = if additional_near_best == 0 {
        0
    } else {
        (4.343_f64 * (additional_near_best as f64 + 1.0).ln() + 0.499) as i32
    };
    let adjusted = raw
        .saturating_sub(repeat_penalty)
        .clamp(0, i32::from(MAPQ_POLICY.paired.maximum_mapq)) as u8;
    if near_best_pairings == 0 {
        adjusted
    } else {
        adjusted.min(MAPQ_POLICY.paired.observed_near_best_cap)
    }
}

pub(super) fn evidence_mapping_quality(result: PairedBatchResult) -> u8 {
    let metrics = result.metrics();
    bwa_pair_mapping_quality_from_evidence(
        result.class(),
        metrics.frontier_complete,
        metrics.mapq_best_pair_score,
        metrics.mapq_second_best_pair_score,
        metrics.mapq_near_best_pairings,
    )
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
    if search_repeat_risk(
        metrics.window_rescue_attempted,
        metrics.resolved_prior_ambiguity,
        metrics.mate1.located_rows,
        metrics.mate2.located_rows,
    ) {
        adjusted_mapq.min(MAPQ_POLICY.paired.search_repeat_risk_cap)
    } else {
        adjusted_mapq
    }
}

fn search_repeat_risk(
    window_rescue_attempted: bool,
    resolved_prior_ambiguity: bool,
    mate1_located_rows: u64,
    mate2_located_rows: u64,
) -> bool {
    window_rescue_attempted
        || resolved_prior_ambiguity
        || mate1_located_rows.max(mate2_located_rows) >= PAIR_MAPQ_REPEAT_RISK_ROWS
}

const fn rescue_risk(window_rescue_attempted: bool) -> bool {
    window_rescue_attempted
}

const fn row_pressure(mate1_located_rows: u64, mate2_located_rows: u64) -> bool {
    if mate1_located_rows > mate2_located_rows {
        mate1_located_rows >= PAIR_MAPQ_REPEAT_RISK_ROWS
    } else {
        mate2_located_rows >= PAIR_MAPQ_REPEAT_RISK_ROWS
    }
}

const fn alternative_margin_certificate(candidate: PairedBatchResult) -> Option<MapqCertificate> {
    if candidate.metrics().alternative_margin_frontier_complete {
        Some(MapqCertificate::CompletedPairAlternativeMargin)
    } else {
        None
    }
}

/// Applies evidence caps without ever increasing confidence.  Keeping this
/// operation separate makes the monotonicity invariant directly testable.
const fn apply_final_mapq_caps(
    mapq: u8,
    stability_mapq: Option<u8>,
    repeat_cap: Option<u8>,
    soft_clipped: bool,
) -> u8 {
    // Preserve the score-gap MAPQ itself.  Confidence reductions below are
    // evidence caps; there is no scientific reason to truncate every complete
    // pair frontier at one fixed tier before considering those caps.
    let mut capped = mapq;
    if let Some(stability) = stability_mapq {
        capped = if capped < stability {
            capped
        } else {
            stability
        };
    }
    if let Some(repeat) = repeat_cap {
        capped = if capped < repeat { capped } else { repeat };
    }
    if soft_clipped && capped > MAPQ_POLICY.paired.soft_clip_cap {
        capped = MAPQ_POLICY.paired.soft_clip_cap;
    }
    capped
}

pub(crate) fn paired_mapping_quality(
    result: PairedBatchResult,
    stability_result: Option<PairedBatchResult>,
    class: PairMappingStatus,
    read_lengths: [usize; 2],
    retained_ranges: [&Range<usize>; 2],
) -> u8 {
    if !matches!(class, PairMappingStatus::Unique) {
        return 0;
    }
    let soft_clipped = read_lengths
        .iter()
        .zip(retained_ranges)
        .any(|(length, retained)| retained.start != 0 || retained.end != *length);
    let baseline_mapq = evidence_mapping_quality(result);
    let rescue_risk = |candidate: PairedBatchResult| {
        let metrics = candidate.metrics();
        rescue_risk(metrics.window_rescue_attempted)
    };
    let resolved_frontier =
        |candidate: PairedBatchResult| candidate.metrics().resolved_prior_ambiguity;
    let row_pressure = |candidate: PairedBatchResult| {
        let metrics = candidate.metrics();
        row_pressure(metrics.mate1.located_rows, metrics.mate2.located_rows)
            && alternative_margin_certificate(candidate).is_none()
    };
    let mut evidence_cap = if rescue_risk(result) || stability_result.is_some_and(rescue_risk) {
        Some(MAPQ_POLICY.paired.rescue_risk_cap)
    } else if resolved_frontier(result) || stability_result.is_some_and(resolved_frontier) {
        Some(MAPQ_POLICY.paired.resolved_frontier_cap)
    } else if row_pressure(result) || stability_result.is_some_and(row_pressure) {
        Some(MAPQ_POLICY.paired.row_pressure_cap)
    } else {
        None
    };
    let alternative_margin_incomplete =
        |candidate: PairedBatchResult| alternative_margin_certificate(candidate).is_none();
    if alternative_margin_incomplete(result)
        || stability_result.is_some_and(alternative_margin_incomplete)
    {
        evidence_cap = Some(
            evidence_cap.map_or(MAPQ_POLICY.paired.incomplete_high_confidence_cap, |cap| {
                cap.min(MAPQ_POLICY.paired.incomplete_high_confidence_cap)
            }),
        );
    }
    apply_final_mapq_caps(
        baseline_mapq,
        stability_result.map(evidence_mapping_quality),
        evidence_cap,
        soft_clipped,
    )
}

#[cfg(test)]
#[path = "../../tests/whitebox/mapq_policy.rs"]
mod whitebox;
