use crate::library::LibraryProfile;
use crate::search::combined_adaptive::{
    CombinedSearchLimits, DEFAULT_SEARCH_LIMITS, SENSITIVE_SEARCH_LIMITS,
};

use super::{PAIRED_MAX_EDIT_DISTANCE, SEMI_GLOBAL_CLIP_PENALTY, SENSITIVE_CLIP_PENALTY};

/// Candidate-search effort for paired-end alignment.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum PairedSearchMode {
    /// Qualified low-latency alignment.
    #[default]
    Default,
    /// Qualified alignment with bounded completion and confidence repair.
    Sensitive,
}

/// Internal remapping phase derived by the paired-end aligner.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub(super) enum AlignmentPhase {
    /// Initial mapping of complete input reads.
    #[default]
    Primary,
    /// Recheck of reads trimmed at a supported adapter boundary.
    AdapterTrimmed,
}

/// Stable options for one paired-end alignment batch.
///
/// Edit distance, rescue, and clipping behavior are derived from the selected
/// mode. Callers cannot compose internal alignment stages.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct PairedAlignmentOptions {
    pub(super) library_profile: LibraryProfile,
    pub(super) search_mode: PairedSearchMode,
    pub(super) minimum_template_span: u64,
    pub(super) maximum_template_span: u64,
    pub(super) phase: AlignmentPhase,
}

impl PairedAlignmentOptions {
    /// Creates options for the initial complete-read alignment.
    #[must_use]
    pub const fn primary(
        library_profile: LibraryProfile,
        search_mode: PairedSearchMode,
        minimum_template_span: u64,
        maximum_template_span: u64,
    ) -> Self {
        Self {
            library_profile,
            search_mode,
            minimum_template_span,
            maximum_template_span,
            phase: AlignmentPhase::Primary,
        }
    }

    pub(super) const fn adapter_trimmed(
        library_profile: LibraryProfile,
        search_mode: PairedSearchMode,
        minimum_template_span: u64,
        maximum_template_span: u64,
    ) -> Self {
        Self {
            library_profile,
            search_mode,
            minimum_template_span,
            maximum_template_span,
            phase: AlignmentPhase::AdapterTrimmed,
        }
    }

    pub(super) const fn derived_policy(self) -> (u8, bool, bool) {
        let sensitive = self.search_mode.is_sensitive();
        (
            PAIRED_MAX_EDIT_DISTANCE,
            sensitive,
            sensitive && matches!(self.phase, AlignmentPhase::Primary),
        )
    }
}

impl PairedSearchMode {
    pub(super) const fn limits(self) -> CombinedSearchLimits {
        match self {
            Self::Default => DEFAULT_SEARCH_LIMITS,
            Self::Sensitive => SENSITIVE_SEARCH_LIMITS,
        }
    }

    pub(super) const fn semi_global_clip_penalty(self) -> u8 {
        match self {
            Self::Default => SEMI_GLOBAL_CLIP_PENALTY,
            Self::Sensitive => SENSITIVE_CLIP_PENALTY,
        }
    }

    pub(super) const fn is_sensitive(self) -> bool {
        matches!(self, Self::Sensitive)
    }
}
