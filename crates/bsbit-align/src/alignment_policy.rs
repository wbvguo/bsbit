//! Versioned, truth-blind alignment search and endpoint policy.
//!
//! This module is the only source of numeric thresholds that can change
//! candidate-frontier coverage, endpoint admission, or representative
//! selection. Low-level implementation dimensions remain with their kernels;
//! MAPQ thresholds remain in `mapq_policy`.

use crate::read_mapping_limits::{INITIAL_EDIT_DISTANCE, MAX_EDIT_DISTANCE};

/// Stable identifier recorded by alignment metrics and validation reports.
///
/// Increment this identifier whenever a numeric value or decision rule in
/// this module changes. A behavior-preserving refactor retains the identifier.
pub const ALIGNMENT_POLICY_ID: &str = "bounded-structural-alignment-v1";

/// Bounded combined-index search effort for one read lane.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct CombinedSearchLimits {
    pub(crate) minimum_multi_hit_seed_bases: usize,
    pub(crate) maximum_seed_hits: u64,
    pub(crate) maximum_combined_rescue_hits: u64,
    pub(crate) maximum_seed_rounds: usize,
}

const INITIAL_MINIMUM_MULTI_HIT_SEED_BASES: usize = 17;
const INITIAL_MAXIMUM_SEED_HITS: u64 = 1_000;
const INITIAL_MAXIMUM_COMBINED_RESCUE_HITS: u64 = 4_096;
pub(crate) const INITIAL_MAXIMUM_SEED_ROUNDS: usize = 5;
pub(crate) const DEFAULT_MINIMUM_MULTI_HIT_SEED_BASES: usize = 16;
pub(crate) const DEFAULT_MAXIMUM_SEED_HITS: u64 = 1_000;
pub(crate) const DEFAULT_MAXIMUM_COMBINED_RESCUE_HITS: u64 = 4_096;
pub(crate) const DEFAULT_MAXIMUM_SEED_ROUNDS: usize = 6;
const SENSITIVE_MAXIMUM_SEED_HITS: u64 = 4_096;
const SENSITIVE_MAXIMUM_COMBINED_RESCUE_HITS: u64 = 4_096;
const SENSITIVE_MAXIMUM_SEED_ROUNDS: usize = 6;
pub(crate) const SEED_PROOF_ROUNDS: usize = 10;
const SENSITIVE_SINGLE_MAXIMUM_SEED_ROUNDS: usize = SEED_PROOF_ROUNDS;

pub(crate) const INITIAL_SEARCH_LIMITS: CombinedSearchLimits = CombinedSearchLimits {
    minimum_multi_hit_seed_bases: INITIAL_MINIMUM_MULTI_HIT_SEED_BASES,
    maximum_seed_hits: INITIAL_MAXIMUM_SEED_HITS,
    maximum_combined_rescue_hits: INITIAL_MAXIMUM_COMBINED_RESCUE_HITS,
    maximum_seed_rounds: INITIAL_MAXIMUM_SEED_ROUNDS,
};

pub(crate) const DEFAULT_SEARCH_LIMITS: CombinedSearchLimits = CombinedSearchLimits {
    minimum_multi_hit_seed_bases: DEFAULT_MINIMUM_MULTI_HIT_SEED_BASES,
    maximum_seed_hits: DEFAULT_MAXIMUM_SEED_HITS,
    maximum_combined_rescue_hits: DEFAULT_MAXIMUM_COMBINED_RESCUE_HITS,
    maximum_seed_rounds: DEFAULT_MAXIMUM_SEED_ROUNDS,
};

/// Narrow default-mode audit used to corroborate provisional confidence.
pub(crate) const DEFAULT_HIGH_CONFIDENCE_AUDIT_SEARCH_LIMITS: CombinedSearchLimits =
    CombinedSearchLimits {
        minimum_multi_hit_seed_bases: DEFAULT_MINIMUM_MULTI_HIT_SEED_BASES,
        maximum_seed_hits: DEFAULT_MAXIMUM_SEED_HITS,
        maximum_combined_rescue_hits: DEFAULT_MAXIMUM_COMBINED_RESCUE_HITS,
        maximum_seed_rounds: DEFAULT_MAXIMUM_SEED_ROUNDS,
    };

pub(crate) const SENSITIVE_SEARCH_LIMITS: CombinedSearchLimits = CombinedSearchLimits {
    minimum_multi_hit_seed_bases: DEFAULT_MINIMUM_MULTI_HIT_SEED_BASES,
    maximum_seed_hits: SENSITIVE_MAXIMUM_SEED_HITS,
    maximum_combined_rescue_hits: SENSITIVE_MAXIMUM_COMBINED_RESCUE_HITS,
    maximum_seed_rounds: SENSITIVE_MAXIMUM_SEED_ROUNDS,
};

pub(crate) const SENSITIVE_SINGLE_SEARCH_LIMITS: CombinedSearchLimits = CombinedSearchLimits {
    minimum_multi_hit_seed_bases: DEFAULT_MINIMUM_MULTI_HIT_SEED_BASES,
    maximum_seed_hits: SENSITIVE_MAXIMUM_SEED_HITS,
    maximum_combined_rescue_hits: SENSITIVE_MAXIMUM_COMBINED_RESCUE_HITS,
    maximum_seed_rounds: SENSITIVE_SINGLE_MAXIMUM_SEED_ROUNDS,
};

/// Search scheduling and bounded local-filter geometry.
pub(crate) const EMPTY_SEED_STEP: usize = 8;
pub(crate) const LOCAL_FILTER_BLOCKS: usize = 8;

/// Adapter-boundary stability and minimum retained alignment domain.
pub(crate) const ADAPTER_STABILITY_DELTA: usize = 8;
pub(crate) const MIN_ADAPTER_RETAINED_BASES: usize = 50;
pub(crate) const ILLUMINA_UNIVERSAL_ADAPTER: &[u8] = b"AGATCGGAAGAGC";
pub(crate) const MIN_ADAPTER_SUPPORT_BASES: usize = 8;
pub(crate) const DEFAULT_ADAPTER_MAX_CLIP_BASES: usize = 30;
pub(crate) const DEFAULT_MAX_SOFT_CLIP_BASES: usize = 30;

/// Edit cost used by the endpoint fallback objective.
pub(crate) const SEMI_GLOBAL_EDIT_PENALTY: u8 = 7;

/// Paired search proof and candidate-enumeration bounds.
pub(crate) const RESCUE_BLOCKS: usize = INITIAL_EDIT_DISTANCE as usize + 1;
pub(crate) const SENSITIVE_RANKED_BLOCK_HITS: u64 = 512;
pub(crate) const SENSITIVE_ALTERNATIVE_MARGIN_BLOCK_HITS: u64 =
    SENSITIVE_RANKED_BLOCK_HITS.saturating_mul(32);
pub(crate) const SENSITIVE_UNMAPPED_RANKED_BLOCK_HITS: u64 =
    SENSITIVE_RANKED_BLOCK_HITS.saturating_mul(2);
pub(crate) const SENSITIVE_DEEP_UNMAPPED_RANKED_BLOCK_HITS: u64 =
    SENSITIVE_RANKED_BLOCK_HITS.saturating_mul(8);
pub(crate) const SENSITIVE_PROOF_BLOCKS: usize = MAX_EDIT_DISTANCE as usize + 1;
pub(crate) const SENSITIVE_ADAPTIVE_MIN_BLOCK_BASES: usize = 19;
pub(crate) const SENSITIVE_BALANCED_BOUNDARY_SHIFTS: [i8; SENSITIVE_PROOF_BLOCKS - 1] =
    [0; SENSITIVE_PROOF_BLOCKS - 1];
pub(crate) const SENSITIVE_ADAPTIVE_BOUNDARY_SHIFTS: [i8; 3] = [-3, 0, 3];
pub(crate) const SEMI_GLOBAL_MAX_EXACT_ANCHOR_HITS: u64 = 256;
pub(crate) const PAIR_ORIGIN_EXACT_SCAN_LIMIT: usize = 64;

/// Paired endpoint admission and representation objective.
pub(crate) const SEMI_GLOBAL_MIN_ALIGNED_BASES: usize = 50;
pub(crate) const SEMI_GLOBAL_ADMISSION_EDIT_PENALTY: u8 = 2;
pub(crate) const SEMI_GLOBAL_CLIP_PENALTY: u8 = 1;
pub(crate) const SENSITIVE_CLIP_PENALTY: u8 = 4;
pub(crate) const ORIGIN_ENDPOINT_CLIP_OPEN_PENALTY: u16 = 8;
pub(crate) const ORIGIN_ENDPOINT_CLIP_EXTENSION_PENALTY: u16 = 7;
pub(crate) const ORIGIN_ENDPOINT_ADAPTER_CLIP_OPEN_PENALTY: u16 = 2;
pub(crate) const ORIGIN_ENDPOINT_ADAPTER_CLIP_EXTENSION_PENALTY: u16 = 0;

pub(crate) const SENSITIVE_MIN_EVENT_PENALTY: u8 =
    if SENSITIVE_CLIP_PENALTY < SEMI_GLOBAL_EDIT_PENALTY {
        SENSITIVE_CLIP_PENALTY
    } else {
        SEMI_GLOBAL_EDIT_PENALTY
    };
