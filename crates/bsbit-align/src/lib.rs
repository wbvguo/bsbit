//! Read-to-reference alignment orchestration.
//!
//! This crate owns sequence-to-sequence algorithms, candidate discovery,
//! seeding, verification scheduling, paired-end selection, search policy, and
//! MAPQ. Architecture-specific kernels are isolated behind private dispatch
//! adapters; reference storage and FM rank/locate primitives live in
//! `bsbit-index`, while stable DNA and chemistry values live in `bsbit-core`.

#![deny(unsafe_code)]

pub mod verification;

mod adapter;
mod alignment_policy;
mod error;
pub mod extension;
pub mod library;
mod mapq_policy;
pub mod materialize;
pub mod paired_end;
pub mod placement;
mod read_mapping;
mod read_mapping_limits;
mod reporting_tie_break;
pub mod score;
pub mod search;
pub mod single_end;

pub use adapter::{AlignmentOutputPolicy, AlignmentOutputPolicyError, SoftClipMode};
pub use alignment_policy::ALIGNMENT_POLICY_ID;
pub use error::AlignmentError;
pub use mapq_policy::MAPQ_POLICY_ID;
pub use read_mapping_limits::{MAX_EDIT_DISTANCE, MAX_READ_BASES, MIN_READ_BASES};
