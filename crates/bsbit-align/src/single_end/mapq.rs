//! Single-end mapping-quality evidence and confidence policy.
//!
//! The mapper supplies evidence already observed while selecting one read;
//! MAPQ calculation performs no additional index search or verification.

use crate::mapq_policy::{MAPQ_POLICY, MapqCertificate};

/// Already-observed evidence used to score one unique single-read origin.
///
/// The policy intentionally consumes only the candidate and verification
/// frontier that selected the alignment.  Producing MAPQ therefore adds no
/// second FM-index search and no second verification pass.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct SingleMapqEvidence {
    pub(crate) read_length: usize,
    pub(crate) best_distance: u8,
    pub(crate) second_best_distance: Option<u8>,
    pub(crate) verified_distance_limit: u8,
    pub(crate) located_rows: u64,
    pub(crate) distinct_candidate_starts: u64,
    pub(crate) verified_placements: u64,
    pub(crate) first_seed_hits: u64,
    pub(crate) first_seed_bases: u64,
    pub(crate) direct_singleton: bool,
    /// True only after the configured sensitive candidate frontier completed.
    pub(crate) frontier_complete: bool,
    /// True after a bounded high-confidence audit reached its configured
    /// search boundary.  Default can establish this selectively without
    /// claiming completion of the wider sensitive candidate frontier.
    pub(crate) confidence_audit_complete: bool,
    /// Distinct offset seed rounds that rediscovered the winner.
    pub(crate) seed_round_support: u8,
}

/// Returns the truth-blind evidence supporting the retained coordinate after
/// the bounded sensitive frontier completes.
#[must_use]
pub(crate) const fn completed_coordinate_certificate(
    evidence: SingleMapqEvidence,
) -> Option<MapqCertificate> {
    let policy = MAPQ_POLICY.single;
    if !evidence.frontier_complete {
        None
    } else if evidence.seed_round_support >= policy.strong_seed_rounds {
        Some(MapqCertificate::CompletedMultiSeedCoordinate)
    } else if evidence.read_length >= policy.long_read_two_seed_min_bases
        && evidence.direct_singleton
        && evidence.seed_round_support >= policy.corroborating_seed_rounds
        && evidence.best_distance <= policy.singleton_max_edit_distance
    {
        Some(MapqCertificate::LongReadSingletonCorroboration)
    } else if evidence.read_length >= policy.long_read_two_seed_min_bases
        && evidence.seed_round_support >= policy.corroborating_seed_rounds
        && evidence.best_distance == evidence.verified_distance_limit
        && evidence.second_best_distance.is_none()
    {
        Some(MapqCertificate::LongReadBoundaryCorroboration)
    } else {
        None
    }
}

/// Returns supporting evidence for an already-collapsed local indel-shift
/// locus.  This does not apply to distinct genomic origins.
#[must_use]
pub(crate) const fn completed_local_locus_certificate(
    evidence: SingleMapqEvidence,
) -> Option<MapqCertificate> {
    if let Some(certificate) = completed_coordinate_certificate(evidence) {
        Some(certificate)
    } else if evidence.frontier_complete
        && evidence.best_distance <= MAPQ_POLICY.single.very_low_edit_distance
        && evidence.seed_round_support >= MAPQ_POLICY.single.corroborating_seed_rounds
    {
        Some(MapqCertificate::LowEditLocalLocusCorroboration)
    } else {
        None
    }
}

const fn completed_multiseed_frontier(evidence: SingleMapqEvidence) -> bool {
    evidence.frontier_complete
        && evidence.seed_round_support >= MAPQ_POLICY.single.corroborating_seed_rounds
        && evidence.best_distance <= MAPQ_POLICY.single.maximum_confident_edit_distance
}

const fn completed_high_confidence_certificate(
    evidence: SingleMapqEvidence,
    raw_mapq: u8,
) -> Option<MapqCertificate> {
    let policy = MAPQ_POLICY.single;
    let short_singleton = evidence.direct_singleton
        && evidence.first_seed_hits == 1
        && evidence.first_seed_bases >= policy.short_singleton_min_seed_bases
        && evidence.first_seed_bases <= policy.short_singleton_max_seed_bases;
    if evidence.confidence_audit_complete
        && evidence.best_distance <= policy.maximum_confident_edit_distance
        && short_singleton
    {
        return Some(MapqCertificate::AuditedShortSingletonCoordinate);
    }
    let audited_multiseed = evidence.confidence_audit_complete
        && evidence.seed_round_support >= policy.corroborating_seed_rounds
        && evidence.best_distance <= policy.maximum_confident_edit_distance;
    let strong_multiseed = audited_multiseed
        && (evidence.seed_round_support >= policy.strong_seed_rounds
            || evidence.seed_round_support >= policy.corroborating_seed_rounds
                && (raw_mapq
                    >= policy
                        .high_confidence_boundary
                        .saturating_add(policy.edit_separation_scale)
                    || raw_mapq >= policy.high_confidence_boundary
                        && evidence.best_distance <= policy.very_low_edit_distance));
    if strong_multiseed {
        Some(MapqCertificate::AuditedStrongMultiSeedCoordinate)
    } else {
        None
    }
}

fn raw_separation_mapq(evidence: SingleMapqEvidence) -> u8 {
    let separation = evidence.second_best_distance.map_or_else(
        || {
            evidence
                .verified_distance_limit
                .saturating_add(1)
                .saturating_sub(evidence.best_distance)
        },
        |second| second.saturating_sub(evidence.best_distance),
    );
    separation
        .saturating_mul(MAPQ_POLICY.single.edit_separation_scale)
        .min(MAPQ_POLICY.single.maximum_mapq)
}

fn apply_adverse_caps(raw_mapq: u8, evidence: SingleMapqEvidence) -> u8 {
    let policy = MAPQ_POLICY.single;
    let mut mapq = raw_mapq;
    if mapq >= policy.high_confidence_boundary
        && completed_high_confidence_certificate(evidence, raw_mapq).is_none()
    {
        mapq = mapq.min(policy.uncertified_high_confidence_cap);
    }

    let moderate_edit_multiseed = evidence.frontier_complete
        && evidence.best_distance == policy.moderate_edit_distance
        && evidence.seed_round_support >= policy.corroborating_seed_rounds;
    if evidence.best_distance >= policy.high_edit_distance {
        mapq = mapq.min(policy.high_edit_cap);
    } else if evidence.best_distance == policy.moderate_edit_distance {
        mapq = mapq.min(if moderate_edit_multiseed {
            policy.moderate_edit_multiseed_cap
        } else {
            policy.moderate_edit_cap
        });
    }

    let repeat_risk = evidence.first_seed_hits > policy.repeat_first_seed_hits
        || evidence.located_rows > policy.repeat_located_rows
        || evidence.distinct_candidate_starts > policy.repeat_candidate_starts
        || evidence.verified_placements > policy.repeat_verified_placements;
    if repeat_risk && !evidence.frontier_complete {
        mapq = mapq.min(policy.incomplete_repeat_cap);
    }
    if !evidence.frontier_complete && evidence.seed_round_support < 1 {
        mapq = mapq.min(policy.incomplete_unconfirmed_cap);
    }
    mapq
}

fn apply_positive_certificates(capped_mapq: u8, raw_mapq: u8, evidence: SingleMapqEvidence) -> u8 {
    let policy = MAPQ_POLICY.single;
    let mut mapq = capped_mapq;
    if mapq == raw_mapq
        && policy.multiseed_increment_boundaries.contains(&raw_mapq)
        && completed_multiseed_frontier(evidence)
        && evidence.seed_round_support >= policy.strong_seed_rounds
    {
        mapq = mapq.saturating_add(policy.edit_separation_scale);
    }
    if completed_coordinate_certificate(evidence).is_some() {
        mapq = mapq.max(policy.coordinate_certificate_floor);
    }
    mapq
}

/// Returns an integer-only structural single-read MAPQ from one retained
/// search frontier.  Whole-read edit separation supplies raw confidence;
/// incomplete frontier evidence, unresolved repeat pressure, and high edit
/// burden can only lower it.  Formal probability calibration remains an
/// external step.
#[must_use]
pub(crate) fn single_mapping_quality_from_evidence(evidence: SingleMapqEvidence) -> u8 {
    // The order is part of the scientific contract: score observed origin
    // separation, apply only confidence-reducing evidence, then admit the
    // small set of explicitly named positive certificates.
    let raw_mapq = raw_separation_mapq(evidence);
    let capped_mapq = apply_adverse_caps(raw_mapq, evidence);
    apply_positive_certificates(capped_mapq, raw_mapq, evidence)
}

/// Applies the confidence cap for one adapter-supported endpoint that remained
/// unique across the primary, clipped, and stability searches.
#[must_use]
pub(crate) fn adapter_consensus_mapping_quality(
    endpoint_is_unique: bool,
    clipped_mapq: u8,
    stability_mapq: Option<u8>,
) -> u8 {
    if endpoint_is_unique {
        clipped_mapq
            .min(stability_mapq.unwrap_or(0))
            .min(MAPQ_POLICY.single.adapter_cap)
    } else {
        0
    }
}

/// Reports whether a provisional default result merits the bounded confidence
/// cross-check. Selection remains in the mapper; this function owns the MAPQ
/// threshold and eligible edit-distance policy.
#[must_use]
pub(crate) const fn default_confidence_cross_check_required(
    mapping_quality: u8,
    best_distance: u8,
) -> bool {
    mapping_quality >= MAPQ_POLICY.single.default_cross_check_min_mapq
        && best_distance <= MAPQ_POLICY.single.very_low_edit_distance
}

/// Reports whether a completed different-origin result carries enough
/// confidence to replace a provisional single-end result.
#[must_use]
pub(crate) const fn sensitive_replacement_certified(mapping_quality: u8) -> bool {
    mapping_quality >= MAPQ_POLICY.single.sensitive_replacement_min_mapq
}

/// Scores a unique origin and applies the policy cap for an indel-shifted
/// local-locus collapse, when one occurred.
#[must_use]
pub(crate) fn local_locus_mapping_quality(
    evidence: SingleMapqEvidence,
    local_origin_collapsed: bool,
) -> u8 {
    let mapping_quality = single_mapping_quality_from_evidence(evidence);
    if !local_origin_collapsed {
        return mapping_quality;
    }
    let cap = if completed_local_locus_certificate(evidence).is_some() {
        MAPQ_POLICY.single.certified_local_locus_cap
    } else {
        MAPQ_POLICY.single.local_locus_cap
    };
    mapping_quality.min(cap)
}

/// Reports whether an ambiguous origin count is inside the bounded affine
/// secondary-ranking domain.
#[must_use]
pub(crate) fn affine_rerank_origin_count_supported(origin_count: u64) -> bool {
    (2..=u64::try_from(MAPQ_POLICY.single.affine_rerank_origin_limit).unwrap_or(u64::MAX))
        .contains(&origin_count)
}

/// Returns the conservative confidence assigned to a completed, uniquely best
/// affine origin, or `None` when the evidence does not certify uniqueness.
#[must_use]
pub(crate) fn affine_unique_mapping_quality(
    frontier_complete: bool,
    best_score: i16,
    runner_up_score: Option<i16>,
) -> Option<u8> {
    let separated = runner_up_score.is_some_and(|runner_up| {
        best_score.saturating_sub(runner_up) >= MAPQ_POLICY.single.affine_unique_min_gap
    });
    (frontier_complete && separated).then_some(MAPQ_POLICY.single.affine_unique_mapq)
}

/// Applies cross-pass repeat pressure to a selected unique result when the
/// merged frontiers remain incomplete.
#[must_use]
pub(crate) fn merged_incomplete_repeat_mapping_quality(
    mapping_quality: u8,
    frontiers_complete: bool,
    located_rows: u64,
    distinct_candidate_starts: u64,
    verified_placements: u64,
) -> u8 {
    let policy = MAPQ_POLICY.single;
    let repeat_pressure = located_rows > policy.repeat_located_rows
        || distinct_candidate_starts > policy.repeat_candidate_starts
        || verified_placements > policy.repeat_verified_placements;
    if !frontiers_complete && mapping_quality < policy.high_confidence_boundary && repeat_pressure {
        mapping_quality.min(policy.incomplete_repeat_cap)
    } else {
        mapping_quality
    }
}

/// Returns the confidence ceiling contributed by an alternative conversion
/// pass.  This reuses the same whole-read separation scale and edit-burden
/// caps as within-pass single-end evidence.
#[must_use]
pub(crate) fn cross_pass_mapping_quality_cap(best_distance: u8, runner_up_distance: u8) -> u8 {
    let policy = MAPQ_POLICY.single;
    let separation = runner_up_distance.saturating_sub(best_distance);
    let mut mapq = separation
        .saturating_mul(policy.edit_separation_scale)
        .min(policy.maximum_mapq);
    if best_distance >= policy.high_edit_distance {
        mapq = mapq.min(policy.high_edit_cap);
    } else if best_distance == policy.moderate_edit_distance {
        mapq = mapq.min(policy.moderate_edit_cap);
    }
    mapq
}

#[cfg(test)]
#[path = "../../tests/whitebox/single_mapq.rs"]
mod whitebox;
