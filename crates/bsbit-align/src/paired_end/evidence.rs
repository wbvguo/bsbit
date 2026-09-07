//! Internal paired-end search evidence shared by mapping, merging, and MAPQ.

use crate::read_mapping::ReadAlignmentMetrics;

use super::result::{PairMappingStatus, PairedPlacement};

/// Outcome of the exact retained-sequence uniqueness check used after
/// sensitive semi-global endpoint completion.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub(super) enum ExactRetainedPairCheck {
    #[default]
    NotRequired,
    NoAlternative,
    AlternativeFound,
    InconclusiveMissingSeed,
    InconclusiveAnchorLimit,
    InconclusiveEmptyAnchorSet,
}

impl ExactRetainedPairCheck {
    pub(super) const fn is_unresolved(self) -> bool {
        !matches!(self, Self::NotRequired | Self::NoAlternative)
    }

    pub(super) const fn found_alternative(self) -> bool {
        matches!(self, Self::AlternativeFound)
    }
}

/// Search and confidence evidence for one paired-end alignment.
// These booleans are independent evidence flags, not mutually exclusive
// states; replacing them with one enum would erase valid combinations.
#[allow(clippy::struct_excessive_bools)]
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) struct PairAlignmentMetrics {
    pub(super) mate1: ReadAlignmentMetrics,
    pub(super) mate2: ReadAlignmentMetrics,
    pub(super) compatible_pairs: u64,
    pub(super) best_pair_placements: u64,
    pub(super) window_rescue_attempted: bool,
    pub(super) semi_global_attempted: bool,
    pub(super) exact_retained_pair_check: ExactRetainedPairCheck,
    /// A complete follow-up search reduced an earlier ambiguous or incomplete
    /// endpoint set to one biological origin.
    pub(super) resolved_prior_ambiguity: bool,
    /// Best compatible pair score in BWA score units (larger is better).
    pub(super) best_pair_score: Option<i16>,
    /// Best strictly lower compatible pair score, when one was observed.
    pub(super) second_best_pair_score: Option<i16>,
    /// Number of alternative pairings within the BWA near-suboptimal window.
    pub(super) near_best_pairings: u64,
    /// Confidence evidence collapsed to distinct biological pair origins.
    pub(super) mapq_compatible_pairs: u64,
    pub(super) mapq_best_pair_score: Option<i16>,
    pub(super) mapq_second_best_pair_score: Option<i16>,
    pub(super) mapq_near_best_pairings: u64,
    /// Whether all candidate work required by the active bounded search ended.
    pub(super) frontier_complete: bool,
    /// Whether the additional alternative-score confidence margin was
    /// completely enumerated.
    pub(super) alternative_margin_frontier_complete: bool,
}

/// One copied result from a cross-pair combined first-seed wavefront.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) struct PairedBatchResult {
    pub(super) class: PairMappingStatus,
    pub(super) metrics: PairAlignmentMetrics,
    pub(super) best_pair: Option<PairedPlacement>,
    pub(super) second_best_distance: Option<u8>,
}

impl PairedBatchResult {
    #[must_use]
    pub(super) const fn class(self) -> PairMappingStatus {
        self.class
    }

    #[must_use]
    pub(super) const fn metrics(self) -> PairAlignmentMetrics {
        self.metrics
    }

    #[must_use]
    pub(super) const fn best_pair(self) -> Option<PairedPlacement> {
        self.best_pair
    }

    #[must_use]
    #[cfg(test)]
    pub(super) const fn best_pair_score(self) -> Option<i16> {
        self.metrics.best_pair_score
    }

    #[must_use]
    #[cfg(test)]
    pub(super) const fn second_best_pair_score(self) -> Option<i16> {
        self.metrics.second_best_pair_score
    }

    #[must_use]
    #[cfg(test)]
    pub(super) const fn near_best_pairings(self) -> u64 {
        self.metrics.near_best_pairings
    }
}
