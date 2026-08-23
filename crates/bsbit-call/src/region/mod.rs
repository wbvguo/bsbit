//! Public region-selection contract.

mod planner;

pub(crate) use planner::{CallRegion, plan_call_regions};

use std::path::PathBuf;

/// One zero-based, half-open genomic interval.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct GenomicInterval {
    /// BAM reference-sequence name.
    pub contig: String,
    /// Zero-based inclusive start coordinate.
    pub start: u64,
    /// Zero-based exclusive end coordinate.
    pub end: u64,
}

/// Optional restriction of calling to selected genomic intervals.
///
/// An empty selection means the whole BAM dictionary. Direct intervals and a
/// regions file are combined as a union; overlaps are merged before work is
/// scheduled, so one locus is never counted twice. The regions file is BED3+
/// with zero-based, half-open coordinates and may be plain, gzip, or BGZF.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct RegionSelection {
    /// Direct zero-based, half-open intervals.
    pub intervals: Vec<GenomicInterval>,
    /// Optional local BED3+ path.
    pub regions_file: Option<PathBuf>,
}
