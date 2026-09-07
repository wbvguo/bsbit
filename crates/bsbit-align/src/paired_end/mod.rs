//! Canonical paired-end alignment API.
//!
//! Public options and results are separated from internal search evidence.
//! Mapping, merging, reporting, rescue, selection, and MAPQ remain focused
//! sibling modules under this layout boundary.

mod evidence;
mod mapper;
mod mapq;
mod merge;
mod options;
mod ranked_blocks;
mod reporting;
mod rescue;
mod result;
mod selection;

pub use mapper::PairedBatchAligner;
pub use options::{PairedAlignmentOptions, PairedSearchMode};
pub use result::{PairMappingStatus, PairedAlignmentResult, PairedPlacement};

/// Maximum number of paired reads accepted by one paired-end search wavefront.
pub const PAIRED_ALIGNMENT_BATCH_SIZE: usize = 32;
