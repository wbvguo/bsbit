//! Public region-selection contract and crate-internal execution facade.

mod planner;
mod workers;

pub(crate) use planner::{CallRegion, plan_call_regions};
pub(crate) use workers::{IndexedCallMode, region_bases_for, stream_indexed_region_workers_mode};

use std::path::PathBuf;

use crate::meth::aggregation::DenseMethRegion;
use crate::snp::result::VariantCall;

#[derive(Debug, Default)]
pub(crate) struct RegionAggregation {
    pub(crate) meth: Option<DenseMethRegion>,
    pub(crate) variants: Vec<(u32, VariantCall)>,
}

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
