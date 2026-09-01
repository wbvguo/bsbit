//! Read-to-reference alignment orchestration.
//!
//! This crate owns sequence-to-sequence algorithms, candidate discovery,
//! seeding, verification scheduling, paired-end selection, search policy, and
//! MAPQ. Architecture-specific kernels are isolated in their own public
//! namespace; reference storage and FM rank/locate primitives live in
//! `bsbit-index`, while stable DNA and chemistry values live in `bsbit-core`.

#![deny(unsafe_code)]

// The complete aligner is safe by default. Architecture-specific intrinsics
// and their checked dispatch wrappers are confined to this implementation
// module and audited under `unsafe_op_in_unsafe_fn = "deny"`.
#[allow(unsafe_code)]
pub mod verification;

mod adapter;
mod error;
pub mod extension;
pub mod library;
pub mod materialize;
pub mod paired_end;
pub mod placement;
mod read_mapping;
mod read_mapping_limits;
pub mod score;
pub mod search;
pub mod single_end;

pub use error::AlignmentError;
