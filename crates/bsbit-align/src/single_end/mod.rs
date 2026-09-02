//! Canonical single-end read-to-reference alignment.
//!
//! The high-throughput mapper owns combined-index batching, single-read
//! classification, and calibrated MAPQ. Layout-neutral search and verification
//! mechanisms remain shared with paired-end mapping outside this module.

mod mapper;
mod mapq;

/// Largest edit-distance budget supported by the single-end mapper.
pub const SINGLE_MAX_EDIT_DISTANCE: u8 = crate::read_mapping_limits::MAX_EDIT_DISTANCE;

pub use mapper::{
    SingleAlignmentResult, SingleBatchAligner, SingleMappingStatus, SingleSearchMode,
};
