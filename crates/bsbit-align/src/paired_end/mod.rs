//! Canonical paired-end read-to-reference alignment.
//!
//! Information-first maximal-suffix seeds feed an integer-only candidate path
//! and a worker-owned edit-distance-three verifier.

use crate::AlignmentError;
use crate::placement::{
    ReadPlacement, SEMI_GLOBAL_EDIT_PENALTY, placement_net_gap_bases, placement_origin_key,
};
#[cfg(test)]
use crate::read_mapping::{
    LOCAL_FILTER_BLOCKS, LocalCandidateFilter, PlacementVerifier, VerificationCacheEntry,
};
use crate::read_mapping::{
    ReadAlignmentMetrics, ReadCandidate, ReadWorkspace, sort_nominal_candidates, strand_index,
};
use crate::read_mapping_limits::{
    INITIAL_EDIT_DISTANCE, MAX_EDIT_DISTANCE, MAX_READ_BASES, MIN_SUFFIX_BASES,
};
use crate::verification::affine::{AffineScoreWorkspace, banded_affine_score};
use crate::verification::ungapped::UngappedEndpoint;
use crate::verification::ungapped::{BoundedSemiglobalConfig, UngappedProfile};

#[cfg(test)]
use self::mapq::{
    PARSIMONY_MAX_LOCATED_ROWS as SENSITIVE_PARSIMONY_MAX_LOCATED_ROWS,
    PARSIMONY_REQUIRED_SCORE_GAP as SENSITIVE_PARSIMONY_REQUIRED_SCORE_GAP,
};
pub use self::mapq::{SENSITIVE_MAPQ_REPEAT_RISK_ROWS, bwa_pair_mapping_quality_from_evidence};
use self::mapq::{
    ambiguity_q10_certified as sensitive_ambiguity_q10_certified,
    effective_mapping_quality as sensitive_effective_mapping_quality,
    incomplete_sparse_completion_required as sensitive_incomplete_sparse_completion_required,
    paired_mapping_quality, stable_rescue_q20_certified as sensitive_stable_rescue_q20_certified,
    two_way_parsimony_q20_certified as sensitive_two_way_parsimony_q20_certified,
};
use crate::library::PairedLibraryProfile;
use crate::search::combined_adaptive::{
    CombinedSearchLimits, CombinedTwoLaneSearchState, DIRECT_SINGLETON_PROOF,
    FLEXIBLE_NOMINAL_PROOF, INITIAL_SEARCH_LIMITS, continue_combined_two_lane_search,
    prepare_combined_projection, prepare_combined_search_projection,
    start_combined_two_lane_search,
};
use crate::search::combined_query::{
    CombinedSearchReferenceExt, CombinedSeedHit, CombinedSeedMatches,
};
use bsbit_core::alphabet::Base;
use bsbit_core::bisulfite::{AlignmentOrientation, BisulfiteStrand, strand_semantics};
use bsbit_index::reference::ReferenceIndex;
use bsbit_index::storage::fm::{ProjectedBase, SearchBase};

mod adapter;
mod batch;
mod frontier;
mod mapq;
mod options;
mod rescue;
mod result;
mod selection;

#[cfg(test)]
use self::adapter::{
    best_ungapped_semi_global_placement, placement_endpoint_cost,
    sequencing_three_prime_adapter_supported, supported_three_prime_adapter_start,
};
#[cfg(test)]
use self::batch::{
    conservatively_mark_incomplete_frontier, empty_pair_metrics,
    merge_non_directional_batch_results, sensitive_repeat_recheck_required,
    sensitive_targeted_semi_global_required, swap_batch_result_mates,
};
#[cfg(test)]
use self::frontier::{
    append_local_flexible_proof_candidates, ranked_block_boundaries,
    selective_unmapped_frontier_deepening_required,
};
#[cfg(test)]
use self::selection::{
    affine_placement_score, collapse_equivalent_pair_origins, pair_net_gap_profile,
    pair_origin_key, placements_may_share_origin, prefer_minimum_net_gap_representative,
    select_best_pair_origins_with_endpoint_policy, select_best_pairs,
    select_best_pairs_with_fallback_score, select_reported_origin_endpoint,
};

pub use options::{PairedAlignmentOptions, PairedSearchMode};
pub use result::{PairMappingStatus, PairedAlignmentResult, PairedPlacement};

/// Maximum number of paired reads mapped in one paired-end worker batch.
pub const PAIRED_ALIGNMENT_BATCH_SIZE: usize = 32;
/// Fixed edit-distance budget used by every supported paired-end search mode.
pub const PAIRED_MAX_EDIT_DISTANCE: u8 = MAX_EDIT_DISTANCE;
/// Largest per-mate edit-distance budget supported by the paired-end mapper.
/// Distances four and five use the generic narrow-band AVX2 kernel only in an
/// incremental fallback; the distance-three specialization is retained for
/// the common first pass inside both supported modes.
// Four disjoint blocks make the local mate-rescue candidate frontier complete
// for the paired-end edit-distance-three budget: at most three edits can
// disturb at most three blocks, leaving at least one exact proof block.
const RESCUE_BLOCKS: usize = INITIAL_EDIT_DISTANCE as usize + 1;
const SENSITIVE_RANKED_BLOCK_HITS: u64 = 512;
const SENSITIVE_UNMAPPED_RANKED_BLOCK_HITS: u64 = SENSITIVE_RANKED_BLOCK_HITS.saturating_mul(2);
const SENSITIVE_SELECTIVE_UNMAPPED_RANKED_BLOCK_HITS: u64 =
    SENSITIVE_RANKED_BLOCK_HITS.saturating_mul(8);
const SENSITIVE_SELECTIVE_UNMAPPED_MIN_RETAINED_HITS: u64 = 32;
const SENSITIVE_POSITIVE_MAPQ_REPORTING_MIN_RETAINED_HITS: u64 = 128;
const SENSITIVE_SELECTIVE_UNMAPPED_MAX_RETAINED_HITS: u64 = 2_049;
const SENSITIVE_POSITIVE_MAPQ_REPORTING_MAX_RETAINED_HITS: u64 = 512;
/// Candidate-row pressure that triggers a bounded second-best completion.
/// This starts below the MAPQ risk threshold so search can prove away false
/// uniqueness before the reporting layer needs to lower confidence.
const SENSITIVE_REPEAT_RECHECK_ROWS: u64 = 256;
const SENSITIVE_PROOF_BLOCKS: usize = MAX_EDIT_DISTANCE as usize + 1;
const SENSITIVE_ADAPTIVE_MIN_BLOCK_BASES: usize = 19;
const SENSITIVE_BALANCED_BOUNDARY_SHIFTS: [i8; SENSITIVE_PROOF_BLOCKS - 1] =
    [0; SENSITIVE_PROOF_BLOCKS - 1];
const SENSITIVE_ADAPTIVE_BOUNDARY_SHIFTS: [i8; 3] = [-3, 0, 3];
// The qualified endpoint search is bounded to a 30-base terminal domain so it
// remains a candidate-local operation rather than an unrestricted local
// aligner.
const SEMI_GLOBAL_MAX_CLIP_BASES: usize = 30;
const SEMI_GLOBAL_MIN_ALIGNED_BASES: usize = 50;
const ADAPTER_STABILITY_DELTA: usize = 8;
const SEMI_GLOBAL_ADMISSION_EDIT_PENALTY: u8 = 2;
const SEMI_GLOBAL_CLIP_PENALTY: u8 = 1;
// Endpoint representation is selected independently from genomic-locus
// ranking in the origin-grouped policy. Unsupported clipping must not win
// merely because several sequencing errors happen to be terminal: its affine
// extension equals the mismatch penalty. Explicit adapter evidence receives a
// separate, favorable clipping prior below.
const ORIGIN_ENDPOINT_CLIP_OPEN_PENALTY: u16 = 8;
const ORIGIN_ENDPOINT_CLIP_EXTENSION_PENALTY: u16 = 7;
const ORIGIN_ENDPOINT_ADAPTER_CLIP_OPEN_PENALTY: u16 = 2;
const ORIGIN_ENDPOINT_ADAPTER_CLIP_EXTENSION_PENALTY: u16 = 0;
const ORIGIN_ENDPOINT_MIN_ADAPTER_SUPPORT: usize = 8;
const ILLUMINA_UNIVERSAL_ADAPTER: &[u8] = b"AGATCGGAAGAGC";
const SENSITIVE_CLIP_PENALTY: u8 = 4;
// A complete terminal mismatch costs seven. Eight makes endpoint selection
// prefer the full-read placement over removing that mismatch;
// seven would still prefer the clipped placement through the retained-edit
// tie-break in pair selection.
const SENSITIVE_MIN_EVENT_PENALTY: u8 = if SENSITIVE_CLIP_PENALTY < SEMI_GLOBAL_EDIT_PENALTY {
    SENSITIVE_CLIP_PENALTY
} else {
    SEMI_GLOBAL_EDIT_PENALTY
};
const SEMI_GLOBAL_MAX_EXACT_ANCHOR_HITS: u64 = 256;

// BWA-MEM-compatible score units used only by the residual sensitive
// selector.  Candidate discovery and the qualified d3/d5 verifier continue to
// use conversion-aware edit distance, so the common path does not pay dynamic
// programming cost.
const BWA_MATCH_SCORE: i16 = 1;
const BWA_MISMATCH_PENALTY: i16 = 4;
const BWA_NEAR_SUBOPTIMAL_DELTA: i16 = 7;

#[derive(Clone, Copy)]
struct RankedBlockSeed {
    matches: CombinedSeedMatches,
    query_offset: u64,
    proof_mask: u8,
}

#[derive(Clone, Copy)]
struct RankedBlockSelection {
    retained_hits: u64,
    complete: bool,
}

type RankedBlockPartition = Option<(u64, [Option<RankedBlockSeed>; SENSITIVE_PROOF_BLOCKS])>;
type EndpointKey = (u16, usize, u8, usize, usize, u8, usize, usize);

/// One disjoint exact block and its already-counted FM interval.
#[derive(Clone, Copy, Debug)]
struct ProofBlock {
    query_start: u16,
    query_end: u16,
}

impl ProofBlock {
    #[must_use]
    const fn query_start(self) -> u16 {
        self.query_start
    }

    #[must_use]
    const fn query_end(self) -> u16 {
        self.query_end
    }
}

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
struct MateRescueWindow {
    contig_ordinal: u64,
    strand: BisulfiteStrand,
    start: u64,
    end: u64,
}

type OriginPairStorageKey = ((u64, u8, i128), (u64, u8, i128));

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct OriginPairEvidence {
    // Larger is better. Each biological origin contributes only its best
    // endpoint score to MAPQ, while raw endpoint statistics remain available
    // to control the search pipeline.
    mapq_score: i16,
}

/// Result counters for one paired-end alignment.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct PairAlignmentMetrics {
    pub(crate) mate1: ReadAlignmentMetrics,
    pub(crate) mate2: ReadAlignmentMetrics,
    pub(crate) compatible_pairs: u64,
    pub(crate) best_pair_placements: u64,
    pub(crate) window_rescue_attempted: bool,
    pub(crate) semi_global_attempted: bool,
    /// Best compatible pair score in BWA score units (larger is better).
    pub(crate) best_pair_score: Option<i16>,
    /// Best strictly lower compatible pair score, when one was observed.
    pub(crate) second_best_pair_score: Option<i16>,
    /// Number of alternative pairings within the BWA near-suboptimal window.
    pub(crate) near_best_pairings: u64,
    /// Confidence evidence collapsed to distinct biological pair origins.
    /// These fields are consumed only by MAPQ; raw fields above continue to
    /// control affine rescoring, rescue, and candidate-search decisions.
    pub(crate) mapq_compatible_pairs: u64,
    pub(crate) mapq_best_pair_score: Option<i16>,
    pub(crate) mapq_second_best_pair_score: Option<i16>,
    pub(crate) mapq_near_best_pairings: u64,
    /// Whether all candidate work required by the active bounded search ended.
    pub(crate) frontier_complete: bool,
}

/// Reusable state for one paired-end mapping worker.
struct PairWorkspace {
    mate1: ReadWorkspace,
    mate2: ReadWorkspace,
    rescue_windows: Vec<MateRescueWindow>,
    best_pairs: Vec<PairedPlacement>,
    exact_anchor_candidates: Vec<ReadCandidate>,
    ranked_anchor_placements: Vec<ReadPlacement>,
    mate1_affine_scores: Vec<i16>,
    mate2_affine_scores: Vec<i16>,
    affine: AffineScoreWorkspace,
    semi_global_clip_penalty: u8,
    prefer_minimum_net_gap: bool,
    origin_pair_evidence: std::collections::HashMap<OriginPairStorageKey, OriginPairEvidence>,
    combined_search_state: CombinedTwoLaneSearchState,
    fallback_mate1_nominals: Vec<ReadCandidate>,
    fallback_mate2_nominals: Vec<ReadCandidate>,
    ranked_extension_selections: [Option<RankedBlockSelection>; 2],
}

/// One copied result from a cross-pair combined first-seed wavefront.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
// These booleans are independent output certificates, not interchangeable
// state flags; a bitfield would make their meanings less explicit.
#[allow(clippy::struct_excessive_bools)]
pub(crate) struct PairedBatchResult {
    class: PairMappingStatus,
    metrics: PairAlignmentMetrics,
    best_pair: Option<PairedPlacement>,
    second_best_distance: Option<u8>,
    repeat_risk_q20_certified: bool,
    parsimony_q20_certified: bool,
    ambiguity_q10_certified: bool,
    requires_positive_mapq_for_reporting: bool,
}

impl PairedBatchResult {
    #[must_use]
    pub(crate) const fn class(self) -> PairMappingStatus {
        self.class
    }

    #[must_use]
    pub(crate) const fn metrics(self) -> PairAlignmentMetrics {
        self.metrics
    }

    #[must_use]
    pub(crate) const fn best_pair(self) -> Option<PairedPlacement> {
        self.best_pair
    }

    /// Returns the selected pair score in BWA score units (larger is better).
    #[must_use]
    #[cfg(test)]
    pub(crate) const fn best_pair_score(self) -> Option<i16> {
        self.metrics.best_pair_score
    }

    /// Returns the best strictly lower pair score, when observed.
    #[must_use]
    #[cfg(test)]
    pub(crate) const fn second_best_pair_score(self) -> Option<i16> {
        self.metrics.second_best_pair_score
    }

    /// Returns the number of alternative pairings close enough to penalize
    /// BWA-style mapping quality.
    #[must_use]
    #[cfg(test)]
    pub(crate) const fn near_best_pairings(self) -> u64 {
        self.metrics.near_best_pairings
    }

    /// Returns the pair-level BWA-style score-gap mapping quality before any
    /// reporting-layer cap for clipping or repeat-risk provenance.
    #[must_use]
    pub(crate) fn evidence_mapping_quality(self) -> u8 {
        bwa_pair_mapping_quality_from_evidence(
            self.class,
            self.metrics.frontier_complete,
            self.metrics.best_pair_score,
            self.metrics.second_best_pair_score,
            self.metrics.near_best_pairings,
        )
    }

    /// Reports whether an independent, bounded endpoint pass supplied enough
    /// stable one-mate-rescue evidence to clear only the MAPQ-20 boundary.
    #[must_use]
    pub(crate) const fn repeat_risk_q20_certified(self) -> bool {
        self.repeat_risk_q20_certified
    }

    /// Reports whether a complete two-way semi-global tie had one uniquely
    /// parsimonious representative inside the qualified Q20 envelope.
    #[must_use]
    pub(crate) const fn parsimony_q20_certified(self) -> bool {
        self.parsimony_q20_certified
    }
}

#[derive(Clone, Copy)]
struct AdapterFallbackResult {
    result: PairedBatchResult,
    stability_result: Option<PairedBatchResult>,
    final_class: PairMappingStatus,
    retained_bases: [usize; 2],
}

/// Worker-owned storage for one combined cross-read-pair seed wavefront.
pub struct PairedBatchAligner {
    pair: PairWorkspace,
    projections: Vec<[[ProjectedBase; MAX_READ_BASES]; 2]>,
    first_seeds: Vec<Option<CombinedSeedMatches>>,
    results: Vec<PairedBatchResult>,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
struct PairSelection {
    compatible_pairs: u64,
    second_best_distance: Option<u8>,
    best_pair_score: Option<i16>,
    second_best_pair_score: Option<i16>,
    near_best_pairings: u64,
    mapq_compatible_pairs: u64,
    mapq_best_pair_score: Option<i16>,
    mapq_second_best_pair_score: Option<i16>,
    mapq_near_best_pairings: u64,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
struct PairScoreConfidence {
    best: Option<i16>,
    second: Option<i16>,
    counts_by_delta: [u64; BWA_NEAR_SUBOPTIMAL_DELTA as usize + 1],
}

#[cfg(test)]
#[path = "../../tests/whitebox/paired_end.rs"]
mod tests;
