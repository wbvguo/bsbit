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

// Production intrinsics and their checked dispatch boundary live in this
// private module; the rest of alignment remains under crate-wide deny(unsafe).
#[allow(unsafe_code)]
mod narrow;
pub(crate) use narrow::narrow_banded_placement_distances_interleaved_batch_d3;
pub use narrow::{
    NarrowBandedError, NarrowBandedResult, NarrowEndpointDistances, NarrowPlacementDistances,
    myers_prefix_distances_u128_batch, narrow_banded_fixed_start_batch,
    narrow_banded_fixed_start_gather_batch, narrow_banded_placement_distances,
    narrow_banded_placement_distances_batch, narrow_banded_placement_distances_batch_d3,
    narrow_banded_placement_distances_batch_d5, narrow_banded_placement_distances_d3,
    narrow_banded_placement_distances_d5, narrow_banded_prefix_batch,
};
#[cfg(test)]
use narrow::{
    NarrowBandedFlavor, myers_prefix_distances_u128_scalar_one, narrow_avx2_available,
    narrow_avx512_available, narrow_banded_fixed_start_batch_with_flavor,
    narrow_banded_prefix_batch_with_flavor, narrow_banded_scalar_one,
    narrow_fixed_start_scalar_one, narrow_placement_distances_scalar,
    narrow_placement_distances_scalar_interleaved, narrow_sse2_available, narrow_sse42_available,
};
#[cfg(all(test, target_arch = "x86_64"))]
use narrow::{
    narrow_placement_distances_avx2, narrow_placement_distances_sse2,
    narrow_placement_distances_sse42,
};

#[cfg(test)]
#[allow(unsafe_code)]
#[path = "../../tests/whitebox/alignment_kernels.rs"]
mod tests;
