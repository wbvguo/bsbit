//! Canonical single-end read-to-reference alignment.
//!
//! Public options and results are separated from combined-index batching,
//! result merging, and evidence-derived MAPQ. Layout-neutral search and
//! verification mechanisms remain shared with paired-end mapping outside this
//! module.

mod mapper;
mod mapq;
mod merge;
mod options;
mod result;

pub use mapper::SingleBatchAligner;
pub use options::SingleSearchMode;
pub use result::{SingleAlignmentResult, SingleMappingStatus};

/// Maximum number of reads accepted by one single-end search wavefront.
pub const SINGLE_ALIGNMENT_BATCH_SIZE: usize = 64;
