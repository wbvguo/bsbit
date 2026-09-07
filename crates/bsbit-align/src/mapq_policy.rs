//! Versioned structural MAPQ policy shared by single- and paired-end mapping.
//!
//! This module is the only source of numeric confidence thresholds.  Search
//! and selection code may collect evidence, but it must not invent local MAPQ
//! cutoffs.  Keeping the policy in one table makes the scientific contract
//! reviewable and lets validation reports pin the exact decision rule.

/// Stable identifier recorded by alignment metrics and validation reports.
///
/// Increment this identifier whenever a numeric threshold or decision rule in
/// `MAPQ_POLICY` changes. Refactors that preserve every decision retain the
/// same identifier.
pub const MAPQ_POLICY_ID: &str = "structural-origin-evidence-v1";

/// Truth-blind evidence that can support a positive confidence declaration.
///
/// Variant names describe observable evidence rather than the MAPQ tier that
/// currently consumes it.  The association between a certificate and a score
/// floor or ceiling belongs to [`MapqPolicy`].
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum MapqCertificate {
    /// Three or more independent offset seeds rediscovered one coordinate.
    CompletedMultiSeedCoordinate,
    /// A long-read direct singleton was rediscovered by another offset seed.
    LongReadSingletonCorroboration,
    /// Independent long-read seeds reached the completed verification bound.
    LongReadBoundaryCorroboration,
    /// Two offsets support one already-collapsed low-edit local locus.
    LowEditLocalLocusCorroboration,
    /// A bounded short suffix became a singleton in a completed confidence audit.
    AuditedShortSingletonCoordinate,
    /// Multiple offsets support the winner in a completed low-edit confidence audit.
    AuditedStrongMultiSeedCoordinate,
    /// The pair search enumerated the configured alternative-score margin.
    CompletedPairAlternativeMargin,
}

/// Single-end structural confidence thresholds.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct SingleMapqPolicy {
    pub(crate) maximum_mapq: u8,
    pub(crate) edit_separation_scale: u8,
    pub(crate) coordinate_certificate_floor: u8,
    pub(crate) high_confidence_boundary: u8,
    pub(crate) uncertified_high_confidence_cap: u8,
    pub(crate) incomplete_unconfirmed_cap: u8,
    pub(crate) incomplete_repeat_cap: u8,
    pub(crate) moderate_edit_distance: u8,
    pub(crate) high_edit_distance: u8,
    pub(crate) moderate_edit_cap: u8,
    pub(crate) moderate_edit_multiseed_cap: u8,
    pub(crate) high_edit_cap: u8,
    pub(crate) maximum_confident_edit_distance: u8,
    pub(crate) corroborating_seed_rounds: u8,
    pub(crate) strong_seed_rounds: u8,
    pub(crate) singleton_max_edit_distance: u8,
    pub(crate) very_low_edit_distance: u8,
    pub(crate) multiseed_increment_boundaries: [u8; 2],
    pub(crate) long_read_two_seed_min_bases: usize,
    pub(crate) short_singleton_min_seed_bases: u64,
    pub(crate) short_singleton_max_seed_bases: u64,
    pub(crate) repeat_first_seed_hits: u64,
    pub(crate) repeat_located_rows: u64,
    pub(crate) repeat_candidate_starts: u64,
    pub(crate) repeat_verified_placements: u64,
    pub(crate) adapter_cap: u8,
    pub(crate) affine_unique_min_gap: i16,
    pub(crate) affine_unique_mapq: u8,
    pub(crate) local_locus_cap: u8,
    pub(crate) certified_local_locus_cap: u8,
    pub(crate) sensitive_replacement_min_mapq: u8,
    pub(crate) default_cross_check_min_mapq: u8,
    pub(crate) affine_rerank_origin_limit: usize,
}

/// Paired-end structural confidence thresholds.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct PairedMapqPolicy {
    pub(crate) maximum_mapq: u8,
    pub(crate) score_gap_scale_hundredths: i32,
    pub(crate) repeat_risk_rows: u64,
    pub(crate) search_repeat_risk_cap: u8,
    pub(crate) rescue_risk_cap: u8,
    pub(crate) resolved_frontier_cap: u8,
    pub(crate) row_pressure_cap: u8,
    pub(crate) incomplete_high_confidence_cap: u8,
    pub(crate) observed_near_best_cap: u8,
    pub(crate) soft_clip_cap: u8,
    pub(crate) confidence_mismatch_penalty: i16,
    pub(crate) near_suboptimal_score_delta: i16,
    pub(crate) confidence_proof_extra_edits: u8,
    pub(crate) targeted_completion_trigger: u8,
}

/// Complete frozen MAPQ decision table.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct MapqPolicy {
    pub(crate) single: SingleMapqPolicy,
    pub(crate) paired: PairedMapqPolicy,
}

pub(crate) const MAPQ_POLICY: MapqPolicy = MapqPolicy {
    single: SingleMapqPolicy {
        maximum_mapq: 60,
        edit_separation_scale: 10,
        coordinate_certificate_floor: 20,
        high_confidence_boundary: 40,
        uncertified_high_confidence_cap: 30,
        incomplete_unconfirmed_cap: 20,
        incomplete_repeat_cap: 10,
        moderate_edit_distance: 4,
        high_edit_distance: 5,
        moderate_edit_cap: 15,
        moderate_edit_multiseed_cap: 20,
        high_edit_cap: 10,
        maximum_confident_edit_distance: 3,
        corroborating_seed_rounds: 2,
        strong_seed_rounds: 3,
        singleton_max_edit_distance: 2,
        very_low_edit_distance: 1,
        multiseed_increment_boundaries: [10, 30],
        long_read_two_seed_min_bases: 128,
        short_singleton_min_seed_bases: 16,
        short_singleton_max_seed_bases: 46,
        repeat_first_seed_hits: 64,
        repeat_located_rows: 256,
        repeat_candidate_starts: 64,
        repeat_verified_placements: 64,
        adapter_cap: 20,
        affine_unique_min_gap: 2,
        affine_unique_mapq: 10,
        local_locus_cap: 10,
        certified_local_locus_cap: 20,
        sensitive_replacement_min_mapq: 10,
        default_cross_check_min_mapq: 20,
        affine_rerank_origin_limit: 15,
    },
    paired: PairedMapqPolicy {
        maximum_mapq: 60,
        score_gap_scale_hundredths: 1_000,
        repeat_risk_rows: 384,
        search_repeat_risk_cap: 19,
        rescue_risk_cap: 20,
        resolved_frontier_cap: 30,
        row_pressure_cap: 30,
        incomplete_high_confidence_cap: 39,
        observed_near_best_cap: 39,
        soft_clip_cap: 20,
        confidence_mismatch_penalty: 4,
        near_suboptimal_score_delta: 7,
        confidence_proof_extra_edits: 2,
        targeted_completion_trigger: 20,
    },
};
