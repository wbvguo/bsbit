use crate::placement::ReadPlacement;

/// One complete paired-end placement.
#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub struct PairedPlacement {
    pub(super) mate1: ReadPlacement,
    pub(super) mate2: ReadPlacement,
    pub(super) template_start: u64,
    pub(super) template_end: u64,
    pub(super) distance: u8,
    pub(super) score: u8,
}

impl PairedPlacement {
    /// Returns the first mate placement.
    #[must_use]
    pub const fn mate1(self) -> ReadPlacement {
        self.mate1
    }

    /// Returns the second mate placement.
    #[must_use]
    pub const fn mate2(self) -> ReadPlacement {
        self.mate2
    }

    /// Returns the inclusive template start coordinate.
    #[must_use]
    pub const fn template_start(self) -> u64 {
        self.template_start
    }

    /// Returns the exclusive template end coordinate.
    #[must_use]
    pub const fn template_end(self) -> u64 {
        self.template_end
    }

    /// Returns the combined pair edit distance.
    #[must_use]
    pub const fn distance(self) -> u8 {
        self.distance
    }

    /// Returns the active pair-selection score.
    #[must_use]
    pub const fn score(self) -> u8 {
        self.score
    }
}

/// Final classification for one paired read.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum PairMappingStatus {
    /// No reportable compatible placement survived the alignment policy.
    Unmapped,
    /// Exactly one reportable biological pair origin survived.
    Unique,
    /// Multiple equally supported origins remain reportable at MAPQ zero.
    Ambiguous,
}

/// Final alignment facts for one paired read.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PairedAlignmentResult {
    pub(super) class: PairMappingStatus,
    pub(super) placement: Option<PairedPlacement>,
    pub(super) retained_query_intervals: [std::ops::Range<usize>; 2],
    pub(super) mapping_quality: u8,
    pub(super) adapter_attempted: bool,
    pub(super) adapter_class: Option<PairMappingStatus>,
    pub(super) adapter_clipped_mates: u8,
    pub(super) adapter_clipped_bases: usize,
    pub(super) semi_global_attempted: bool,
    pub(super) semi_global_clipped_mates: u8,
    pub(super) semi_global_clipped_bases: usize,
    pub(super) mate_rescue_attempted: bool,
}

impl PairedAlignmentResult {
    /// Returns the final reporting classification.
    #[must_use]
    pub const fn class(&self) -> PairMappingStatus {
        self.class
    }

    /// Returns the final pair placement, if the pair is reportable.
    #[must_use]
    pub const fn placement(&self) -> Option<PairedPlacement> {
        self.placement
    }

    /// Returns the retained sequencing-orientation interval for each mate.
    #[must_use]
    pub fn retained_query_intervals(&self) -> [std::ops::Range<usize>; 2] {
        self.retained_query_intervals.clone()
    }

    /// Returns the final calibrated pair mapping quality.
    #[must_use]
    pub const fn mapping_quality(&self) -> u8 {
        self.mapping_quality
    }

    /// Reports whether exact adapter support triggered a trimmed remap.
    #[must_use]
    pub const fn adapter_attempted(&self) -> bool {
        self.adapter_attempted
    }

    /// Returns the adapter-remap class before final reporting admission.
    #[must_use]
    pub const fn adapter_class(&self) -> Option<PairMappingStatus> {
        self.adapter_class
    }

    /// Returns the number of mates trimmed at an adapter boundary.
    #[must_use]
    pub const fn adapter_clipped_mates(&self) -> u8 {
        self.adapter_clipped_mates
    }

    /// Returns the total number of adapter-trimmed bases.
    #[must_use]
    pub const fn adapter_clipped_bases(&self) -> usize {
        self.adapter_clipped_bases
    }

    /// Reports whether endpoint completion ran.
    #[must_use]
    pub const fn semi_global_attempted(&self) -> bool {
        self.semi_global_attempted
    }

    /// Returns the number of mates clipped by endpoint completion.
    #[must_use]
    pub const fn semi_global_clipped_mates(&self) -> u8 {
        self.semi_global_clipped_mates
    }

    /// Returns the total number of bases clipped by endpoint completion.
    #[must_use]
    pub const fn semi_global_clipped_bases(&self) -> usize {
        self.semi_global_clipped_bases
    }

    /// Reports whether bounded mate rescue was attempted.
    #[must_use]
    pub const fn mate_rescue_attempted(&self) -> bool {
        self.mate_rescue_attempted
    }
}
