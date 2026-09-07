//! Single-end mapping classification and final result values.

use crate::placement::ReadPlacement;

/// Final classification for one single read.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SingleMappingStatus {
    /// No verified placement survived the bounded search.
    Unmapped,
    /// Exactly one best biological origin survived the configured alignment
    /// objective and confidence policy.
    Unique,
    /// Multiple plausible biological origins survived the confidence policy.
    Ambiguous,
}

/// Final mapping facts for one single read.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SingleAlignmentResult {
    pub(super) status: SingleMappingStatus,
    pub(super) placement: Option<ReadPlacement>,
    pub(super) retained_query_end: usize,
    pub(super) mapping_quality: u8,
    pub(super) located_rows: u64,
    pub(super) distinct_candidate_starts: u64,
    pub(super) verified_placements: u64,
    pub(super) best_origin_count: u64,
    pub(super) adapter_attempted: bool,
    pub(super) adapter_status: Option<SingleMappingStatus>,
    pub(super) adapter_clipped_bases: usize,
}

impl SingleAlignmentResult {
    /// Returns the final mapping class.
    #[must_use]
    pub const fn status(self) -> SingleMappingStatus {
        self.status
    }

    /// Returns the deterministic representative placement, when mapped.
    #[must_use]
    pub const fn placement(self) -> Option<ReadPlacement> {
        self.placement
    }

    /// Returns the retained sequencing-orientation query interval.
    #[must_use]
    pub const fn retained_query_interval(self) -> core::ops::Range<usize> {
        0..self.retained_query_end
    }

    /// Returns the evidence-derived SAM mapping quality, or zero when not unique.
    #[must_use]
    pub const fn mapping_quality(self) -> u8 {
        self.mapping_quality
    }

    /// Returns suffix rows located across every executed mapping phase.
    #[must_use]
    pub const fn located_rows(self) -> u64 {
        self.located_rows
    }

    /// Returns verified placements across every executed mapping phase before
    /// per-phase best-tier selection.
    #[must_use]
    pub const fn verified_placements(self) -> u64 {
        self.verified_placements
    }

    /// Reports whether exact adapter support triggered a trimmed remap.
    #[must_use]
    pub const fn adapter_attempted(self) -> bool {
        self.adapter_attempted
    }

    /// Returns the adapter-remap class after stability verification.
    #[must_use]
    pub const fn adapter_status(self) -> Option<SingleMappingStatus> {
        self.adapter_status
    }

    /// Returns the number of bases omitted at the supported 3' adapter boundary.
    #[must_use]
    pub const fn adapter_clipped_bases(self) -> usize {
        self.adapter_clipped_bases
    }

    pub(super) const fn unmapped(
        read_length: usize,
        located_rows: u64,
        verified_placements: u64,
    ) -> Self {
        Self {
            status: SingleMappingStatus::Unmapped,
            placement: None,
            retained_query_end: read_length,
            mapping_quality: 0,
            located_rows,
            distinct_candidate_starts: 0,
            verified_placements,
            best_origin_count: 0,
            adapter_attempted: false,
            adapter_status: None,
            adapter_clipped_bases: 0,
        }
    }

    pub(super) const fn unmapped_with_evidence(
        read_length: usize,
        located_rows: u64,
        distinct_candidate_starts: u64,
        verified_placements: u64,
    ) -> Self {
        Self {
            status: SingleMappingStatus::Unmapped,
            placement: None,
            retained_query_end: read_length,
            mapping_quality: 0,
            located_rows,
            distinct_candidate_starts,
            verified_placements,
            best_origin_count: 0,
            adapter_attempted: false,
            adapter_status: None,
            adapter_clipped_bases: 0,
        }
    }
}
