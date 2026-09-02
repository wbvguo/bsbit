//! Shared alignment verification and traceback algorithms.
//!
//! Low-level kernels accept caller-provided equality masks, while semantic
//! modules depend only on stable domain values from `bsbit-core`. Global,
//! affine, CIGAR-replay, and bounded ungapped endpoint algorithms live here.
//! This implementation layer deliberately has no reference-index ownership,
//! candidate search, pairing, or MAPQ policy.

#![deny(unsafe_op_in_unsafe_fn)]

pub mod affine;
pub mod cigar;
pub mod distance;
pub(crate) mod prefix_filter;
pub mod ungapped;

/// Maximum query length handled by one machine-word Myers state.
pub(crate) const MAX_QUERY_BASES: usize = 64;
/// Largest edit band represented by one 32-bit narrow-band state.
pub const MAX_NARROW_BAND_DISTANCE: usize = 15;

/// Safe symbol conversion used by gather-style narrow-band verification.
pub trait NarrowReferenceCode: Copy {
    /// Returns a caller-defined code in `0..=4`; larger values mismatch.
    fn narrow_reference_code(self) -> u8;
}

impl NarrowReferenceCode for u8 {
    fn narrow_reference_code(self) -> u8 {
        self
    }
}

mod narrow;
pub use narrow::*;

#[cfg(test)]
mod qualification;
#[cfg(test)]
pub(crate) use qualification::{KernelFlavor, myers_distance_batch};
