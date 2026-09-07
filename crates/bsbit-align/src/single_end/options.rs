//! Single-end candidate-search effort.

use crate::alignment_policy::{
    CombinedSearchLimits, DEFAULT_SEARCH_LIMITS, SENSITIVE_SINGLE_SEARCH_LIMITS,
};

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
    pub(super) const fn limits(self) -> CombinedSearchLimits {
        match self {
            Self::Default => DEFAULT_SEARCH_LIMITS,
            Self::Sensitive => SENSITIVE_SINGLE_SEARCH_LIMITS,
        }
    }

    pub(super) const fn completes_candidate_frontier(self) -> bool {
        matches!(self, Self::Sensitive)
    }
}
