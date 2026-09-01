//! Canonical single-end read-to-reference alignment.
//!
//! The high-throughput mapper owns combined-index batching, single-read
//! classification, and calibrated MAPQ. Layout-neutral search and verification
//! mechanisms remain shared with paired-end mapping outside this module.

mod mapper;
mod mapq;

pub use mapper::{
    SingleAlignmentResult, SingleBatchAligner, SingleMappingStatus, SingleSearchMode,
};
