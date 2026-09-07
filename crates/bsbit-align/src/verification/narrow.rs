//! Narrow-band scalar and architecture-dispatched alignment kernels.

// Keep the five-byte relation table borrowed across public, function-pointer,
// and target-feature boundaries. Copying it into every dispatch trampoline
// adds work to the batch entry path without simplifying either ABI.
#![allow(clippy::trivially_copy_pass_by_ref)]

use super::{MAX_NARROW_BAND_DISTANCE, NarrowReferenceCode};
#[cfg(all(test, target_arch = "x86_64"))]
use bsbit_cpu::configuration;
use bsbit_cpu::{Backend, selected_backend};
use core::fmt;
use std::sync::OnceLock;

#[cfg(target_arch = "aarch64")]
#[path = "narrow_neon.rs"]
mod neon;

#[cfg(target_arch = "x86_64")]
#[path = "narrow_x86_64.rs"]
mod x86_64;

#[cfg(target_arch = "x86_64")]
#[path = "narrow_x86_64_avx512.rs"]
mod x86_64_avx512;

#[cfg(all(test, target_arch = "x86_64"))]
pub(super) use x86_64::{
    narrow_placement_distances_avx2, narrow_placement_distances_sse2,
    narrow_placement_distances_sse42,
};

/// Compact exact distances for every in-budget narrow-band endpoint.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct NarrowEndpointDistances {
    packed: [u64; 2],
    in_budget_mask: u32,
}

impl NarrowEndpointDistances {
    /// Empty endpoint frontier.
    pub const EMPTY: Self = Self {
        packed: [0; 2],
        in_budget_mask: 0,
    };

    /// Returns the mask of endpoints whose exact distance is within budget.
    #[must_use]
    pub const fn in_budget_mask(self) -> u32 {
        self.in_budget_mask
    }

    /// Returns the distance at `delta`, or `None` when that endpoint is outside
    /// the budget or the 31-endpoint narrow domain.
    #[must_use]
    pub fn distance(self, delta: usize) -> Option<u32> {
        if delta > MAX_NARROW_BAND_DISTANCE * 2 || self.in_budget_mask & (1_u32 << delta) == 0 {
            return None;
        }
        let lane = delta / 16;
        let shift = (delta % 16) * 4;
        Some(((self.packed[lane] >> shift) & 0xf) as u32)
    }

    /// Returns every in-budget endpoint whose exact distance equals `distance`.
    #[must_use]
    pub fn mask_at_distance(self, distance: u32) -> u32 {
        let mut matching = 0_u32;
        let mut remaining = self.in_budget_mask;
        while remaining != 0 {
            let delta = remaining.trailing_zeros() as usize;
            let bit = 1_u32 << delta;
            remaining &= remaining - 1;
            if self.distance(delta) == Some(distance) {
                matching |= bit;
            }
        }
        matching
    }

    pub(super) fn insert(&mut self, delta: usize, distance: u32) {
        debug_assert!(delta < MAX_NARROW_BAND_DISTANCE * 2 + 1);
        debug_assert!(
            usize::try_from(distance).is_ok_and(|value| value <= MAX_NARROW_BAND_DISTANCE)
        );
        let lane = delta / 16;
        let shift = (delta % 16) * 4;
        self.packed[lane] |= u64::from(distance) << shift;
        self.in_budget_mask |= 1_u32 << delta;
    }
}

const MAX_NARROW_ENDPOINTS: usize = MAX_NARROW_BAND_DISTANCE * 2 + 1;

/// Minimum restricted-band distance for every fixed start/end interval.
#[derive(Clone, Copy, Debug)]
pub struct NarrowPlacementDistances {
    distances: [[core::mem::MaybeUninit<u8>; MAX_NARROW_ENDPOINTS]; MAX_NARROW_ENDPOINTS],
    band_length: u8,
    max_distance: u8,
}

impl PartialEq for NarrowPlacementDistances {
    fn eq(&self, other: &Self) -> bool {
        if self.band_length != other.band_length || self.max_distance != other.max_distance {
            return false;
        }
        let band_length = usize::from(self.band_length);
        (0..band_length).all(|start| {
            (0..band_length).all(|endpoint| {
                // SAFETY: every constructor initializes the complete active
                // `band_length * band_length` square before publishing it.
                unsafe {
                    self.distances[start][endpoint].assume_init()
                        == other.distances[start][endpoint].assume_init()
                }
            })
        })
    }
}

impl Eq for NarrowPlacementDistances {}

impl NarrowPlacementDistances {
    /// Empty placeholder overwritten by a placement kernel call.
    pub const EMPTY: Self = Self {
        distances: [[core::mem::MaybeUninit::uninit(); MAX_NARROW_ENDPOINTS]; MAX_NARROW_ENDPOINTS],
        band_length: 0,
        max_distance: 0,
    };

    fn for_band(band_length: usize, max_distance: usize) -> Self {
        Self {
            band_length: u8::try_from(band_length).expect("validated narrow band fits u8"),
            max_distance: u8::try_from(max_distance).expect("validated narrow distance fits u8"),
            ..Self::EMPTY
        }
    }

    fn insert_distance(&mut self, start_delta: usize, endpoint_delta: usize, distance: u8) {
        debug_assert!(start_delta < usize::from(self.band_length));
        debug_assert!(endpoint_delta < usize::from(self.band_length));
        self.distances[start_delta][endpoint_delta].write(distance);
    }

    /// Returns the minimum in-budget distance for one fixed interval.
    ///
    /// `start_delta` is the alignment start relative to the verification
    /// window. `endpoint_delta` selects the exclusive endpoint at
    /// `query_length + endpoint_delta`. `None` means that either coordinate is
    /// outside the band or the interval's minimum distance exceeds the budget.
    #[must_use]
    pub fn distance(&self, start_delta: usize, endpoint_delta: usize) -> Option<u32> {
        if start_delta >= usize::from(self.band_length)
            || endpoint_delta >= usize::from(self.band_length)
        {
            return None;
        }
        // SAFETY: coordinates inside the published band are initialized by
        // every scalar and SIMD constructor before this method can be called.
        let distance = unsafe { self.distances[start_delta][endpoint_delta].assume_init() };
        (distance <= self.max_distance).then_some(u32::from(distance))
    }
}

/// One selected narrow-band prefix and its complete in-budget endpoint
/// frontier for a candidate pattern.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct NarrowBandedResult {
    /// Minimum edit distance, or `u32::MAX` when it exceeds the band.
    pub distance: u32,
    /// Selected exclusive pattern-prefix length, or `usize::MAX` when absent.
    pub prefix_length: usize,
    /// Bit `d` is set when prefix `query_length + d` attains `distance`.
    pub tied_prefix_mask: u32,
    /// Exact distance for every endpoint that is within the edit budget.
    pub endpoint_distances: NarrowEndpointDistances,
}

impl NarrowBandedResult {
    pub(super) const ABSENT: Self = Self {
        distance: u32::MAX,
        prefix_length: usize::MAX,
        tied_prefix_mask: 0,
        endpoint_distances: NarrowEndpointDistances::EMPTY,
    };
}

/// Private implementation selector shared by dispatch trampolines and
/// differential kernel tests.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[cfg_attr(not(target_arch = "x86_64"), allow(dead_code))]
pub(super) enum NarrowBandedFlavor {
    /// Portable one-candidate 32-bit implementation.
    Scalar,
    /// Four independent candidate patterns in SSE2 32-bit lanes.
    Sse2,
    /// Four independent candidate patterns in SSE4.2-qualified 32-bit lanes.
    Sse42,
    /// Eight independent candidate patterns in AVX2 32-bit lanes.
    Avx2,
    /// Sixteen independent candidate patterns in AVX-512 32-bit lanes.
    Avx512,
}

/// Invalid dimensions for a narrow-band candidate batch.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum NarrowBandedError {
    /// The query was empty.
    EmptyQuery,
    /// The supplied query length exceeds the two-word prefix domain.
    QueryLength {
        /// Supplied query bases.
        observed: usize,
    },
    /// The edit-distance band cannot fit one 32-bit state.
    Band {
        /// Supplied maximum edit distance.
        observed: usize,
    },
    /// Flat candidate-pattern storage has the wrong length.
    PatternDimension {
        /// Required flat pattern bytes.
        expected: usize,
        /// Supplied flat pattern bytes.
        observed: usize,
    },
    /// Output count does not equal the candidate count.
    OutputDimension {
        /// Candidate pattern count.
        candidates: usize,
        /// Supplied output slots.
        outputs: usize,
    },
    /// A forced test implementation is unavailable.
    #[cfg(test)]
    UnsupportedFlavor,
    /// A placement batch cannot fit in the 32 byte lanes.
    PlacementBatch {
        /// Supplied candidate count.
        observed: usize,
        /// Maximum candidate count for this band.
        maximum: usize,
    },
}

impl fmt::Display for NarrowBandedError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyQuery => formatter.write_str("narrow-band query is empty"),
            Self::QueryLength { observed } => {
                write!(formatter, "prefix query length {observed} exceeds 128")
            }
            Self::Band { observed } => write!(
                formatter,
                "narrow-band distance {observed} exceeds {MAX_NARROW_BAND_DISTANCE}"
            ),
            Self::PatternDimension { expected, observed } => write!(
                formatter,
                "narrow-band pattern bytes differ: expected {expected}, observed {observed}"
            ),
            Self::OutputDimension {
                candidates,
                outputs,
            } => write!(
                formatter,
                "narrow-band candidate/output counts differ: {candidates}/{outputs}"
            ),
            #[cfg(test)]
            Self::UnsupportedFlavor => {
                formatter.write_str("forced narrow-band test kernel is unavailable")
            }
            Self::PlacementBatch { observed, maximum } => write!(
                formatter,
                "narrow placement batch {observed} exceeds SIMD capacity {maximum}"
            ),
        }
    }
}

impl std::error::Error for NarrowBandedError {}

type PlacementKernel = fn(&[u8; 5], &[u8], &[u8], usize) -> NarrowPlacementDistances;
type SpecializedPlacementKernel = fn(&[u8; 5], &[u8], &[u8]) -> NarrowPlacementDistances;
type PlacementBatchKernel =
    fn(&[u8; 5], &[u8], &[u8], usize, usize, &mut [NarrowPlacementDistances]);
type InterleavedPlacementBatchD3Kernel =
    fn(&[u8; 5], &[u8], &[u8], &mut [NarrowPlacementDistances]);
type PrefixBatchKernel = fn([u8; 5], &[u8], &[u8], usize, usize, &mut [NarrowBandedResult]);
type FixedStartBatchKernel =
    fn([u8; 5], &[u8], &[u8], usize, usize, &mut [NarrowEndpointDistances]);

/// Alignment-owned function table constructed from the process CPU policy.
///
/// Only this adapter maps a `bsbit-cpu` backend to alignment kernels. Callers
/// elsewhere in the aligner invoke safe algorithm functions and never inspect
/// an instruction-set choice.
#[derive(Clone, Copy)]
struct NarrowDispatch {
    placement: PlacementKernel,
    placement_d3: SpecializedPlacementKernel,
    placement_d5: SpecializedPlacementKernel,
    placement_batch: PlacementBatchKernel,
    placement_batch_d3: PlacementBatchKernel,
    placement_interleaved_batch_d3: InterleavedPlacementBatchD3Kernel,
    placement_batch_d5: PlacementBatchKernel,
    prefix_batch: PrefixBatchKernel,
    fixed_start_batch: FixedStartBatchKernel,
    gather_avx2: bool,
    myers_u128_avx2: bool,
}

impl NarrowDispatch {
    const fn scalar() -> Self {
        Self {
            placement: placement_scalar_dispatch,
            placement_d3: placement_d3_scalar_dispatch,
            placement_d5: placement_d5_scalar_dispatch,
            placement_batch: placement_batch_scalar_dispatch,
            placement_batch_d3: placement_batch_scalar_dispatch,
            placement_interleaved_batch_d3: placement_interleaved_batch_d3_scalar_dispatch,
            placement_batch_d5: placement_batch_scalar_dispatch,
            prefix_batch: prefix_batch_scalar_dispatch,
            fixed_start_batch: fixed_start_batch_scalar_dispatch,
            gather_avx2: false,
            myers_u128_avx2: false,
        }
    }

    #[cfg(target_arch = "x86_64")]
    const fn sse2() -> Self {
        Self {
            placement: placement_sse2_dispatch,
            placement_d3: placement_d3_sse2_dispatch,
            placement_d5: placement_d5_sse2_dispatch,
            placement_batch: placement_batch_sse2_dispatch,
            placement_batch_d3: placement_batch_sse2_dispatch,
            placement_interleaved_batch_d3: placement_interleaved_batch_d3_scalar_dispatch,
            placement_batch_d5: placement_batch_sse2_dispatch,
            prefix_batch: prefix_batch_sse2_dispatch,
            fixed_start_batch: fixed_start_batch_sse2_dispatch,
            gather_avx2: false,
            myers_u128_avx2: false,
        }
    }

    #[cfg(target_arch = "x86_64")]
    const fn sse42_popcnt() -> Self {
        Self {
            placement: placement_sse42_dispatch,
            placement_d3: placement_d3_sse42_dispatch,
            placement_d5: placement_d5_sse42_dispatch,
            placement_batch: placement_batch_sse42_dispatch,
            placement_batch_d3: placement_batch_sse42_dispatch,
            placement_interleaved_batch_d3: placement_interleaved_batch_d3_scalar_dispatch,
            placement_batch_d5: placement_batch_sse42_dispatch,
            prefix_batch: prefix_batch_sse42_dispatch,
            fixed_start_batch: fixed_start_batch_sse42_dispatch,
            gather_avx2: false,
            myers_u128_avx2: false,
        }
    }

    #[cfg(target_arch = "x86_64")]
    const fn avx2_popcnt() -> Self {
        Self {
            placement: placement_avx2_dispatch,
            placement_d3: placement_d3_avx2_dispatch,
            placement_d5: placement_d5_avx2_dispatch,
            placement_batch: placement_batch_avx2_dispatch,
            placement_batch_d3: placement_batch_d3_avx2_dispatch,
            placement_interleaved_batch_d3: placement_interleaved_batch_d3_avx2_dispatch,
            placement_batch_d5: placement_batch_d5_avx2_dispatch,
            prefix_batch: prefix_batch_avx2_dispatch,
            fixed_start_batch: fixed_start_batch_avx2_dispatch,
            gather_avx2: true,
            myers_u128_avx2: true,
        }
    }

    #[cfg(target_arch = "x86_64")]
    const fn avx512_bw_popcnt() -> Self {
        Self {
            placement: placement_avx2_dispatch,
            placement_d3: placement_d3_avx2_dispatch,
            placement_d5: placement_d5_avx2_dispatch,
            placement_batch: placement_batch_avx2_dispatch,
            placement_batch_d3: placement_batch_d3_avx2_dispatch,
            placement_interleaved_batch_d3: placement_interleaved_batch_d3_avx2_dispatch,
            placement_batch_d5: placement_batch_d5_avx2_dispatch,
            prefix_batch: prefix_batch_avx512_dispatch,
            fixed_start_batch: fixed_start_batch_avx512_dispatch,
            gather_avx2: true,
            myers_u128_avx2: true,
        }
    }

    #[cfg(target_arch = "aarch64")]
    const fn neon() -> Self {
        Self {
            placement: placement_neon_dispatch,
            placement_d3: placement_d3_neon_dispatch,
            placement_d5: placement_d5_neon_dispatch,
            placement_batch: placement_batch_neon_dispatch,
            placement_batch_d3: placement_batch_neon_dispatch,
            placement_interleaved_batch_d3: placement_interleaved_batch_d3_scalar_dispatch,
            placement_batch_d5: placement_batch_neon_dispatch,
            prefix_batch: prefix_batch_neon_dispatch,
            fixed_start_batch: fixed_start_batch_neon_dispatch,
            gather_avx2: false,
            myers_u128_avx2: false,
        }
    }

    fn for_backend(backend: Backend) -> Self {
        match backend {
            Backend::Scalar => Self::scalar(),
            #[cfg(target_arch = "x86_64")]
            Backend::Sse2 => Self::sse2(),
            #[cfg(target_arch = "x86_64")]
            Backend::Sse42Popcnt => Self::sse42_popcnt(),
            #[cfg(target_arch = "x86_64")]
            Backend::Avx2Popcnt => Self::avx2_popcnt(),
            #[cfg(target_arch = "x86_64")]
            Backend::Avx512BwPopcnt => Self::avx512_bw_popcnt(),
            #[cfg(not(target_arch = "x86_64"))]
            Backend::Sse2
            | Backend::Sse42Popcnt
            | Backend::Avx2Popcnt
            | Backend::Avx512BwPopcnt => Self::scalar(),
            #[cfg(target_arch = "aarch64")]
            Backend::Neon => Self::neon(),
            #[cfg(not(target_arch = "aarch64"))]
            Backend::Neon => Self::scalar(),
        }
    }
}

static NARROW_DISPATCH: OnceLock<NarrowDispatch> = OnceLock::new();

fn narrow_dispatch() -> &'static NarrowDispatch {
    NARROW_DISPATCH.get_or_init(|| NarrowDispatch::for_backend(selected_backend()))
}

/// Computes the complete fixed-start/fixed-end frontier inside the narrow
/// verification band.
///
/// The pattern contains exactly `query.len() + 2 * max_distance` codes. Both
/// start and endpoint deltas range over `0..=2 * max_distance`; paths remain
/// inside that same diagonal band. Each returned cell is the minimum distance
/// for its complete `(start, end)` interval, rather than only the minimum for
/// an endpoint after free-start minimization.
///
/// # Errors
///
/// Rejects an empty query, an unsupported band, or a pattern with the wrong
/// fixed width before evaluating any cell.
pub fn narrow_banded_placement_distances(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    pattern: &[u8],
    max_distance: usize,
) -> Result<NarrowPlacementDistances, NarrowBandedError> {
    validate_narrow_pattern(query, pattern, max_distance)?;
    Ok((narrow_dispatch().placement)(
        reference_masks_by_query,
        query,
        pattern,
        max_distance,
    ))
}

/// Computes a distance-three frontier with a seven-vector state.
///
/// The general kernel reserves the maximum supported number of diagonals so
/// it can serve every band width. This fixed-distance entry point avoids
/// carrying unused vector slots through its hot recurrence.
///
/// # Errors
///
/// Rejects an empty query or a pattern whose fixed distance-three width is
/// inconsistent with the query.
pub fn narrow_banded_placement_distances_d3(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    pattern: &[u8],
) -> Result<NarrowPlacementDistances, NarrowBandedError> {
    const MAX_DISTANCE: usize = 3;

    validate_narrow_pattern(query, pattern, MAX_DISTANCE)?;
    Ok((narrow_dispatch().placement_d3)(
        reference_masks_by_query,
        query,
        pattern,
    ))
}

/// Computes a distance-five frontier with an eleven-vector state.
///
/// This is the exact distance-five specialization of
/// [`narrow_banded_placement_distances`]. It avoids reserving and rotating the
/// unused distance-six-through-fifteen diagonal state in the repeat-audit hot
/// path.
///
/// # Errors
///
/// Rejects an empty query or a pattern whose fixed distance-five width is
/// inconsistent with the query.
pub fn narrow_banded_placement_distances_d5(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    pattern: &[u8],
) -> Result<NarrowPlacementDistances, NarrowBandedError> {
    const MAX_DISTANCE: usize = 5;

    validate_narrow_pattern(query, pattern, MAX_DISTANCE)?;
    Ok((narrow_dispatch().placement_d5)(
        reference_masks_by_query,
        query,
        pattern,
    ))
}

/// Computes complete start/end frontiers for several independent candidate
/// patterns in one SIMD operation. At edit distance three, four candidates
/// occupy 28 of the 32 AVX2 byte lanes.
///
/// # Errors
///
/// Rejects an empty query, an unsupported band, or inconsistent pattern and
/// output dimensions.
pub fn narrow_banded_placement_distances_batch(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    patterns: &[u8],
    max_distance: usize,
    output: &mut [NarrowPlacementDistances],
) -> Result<(), NarrowBandedError> {
    let pattern_len = query
        .len()
        .checked_add(max_distance.saturating_mul(2))
        .ok_or(NarrowBandedError::QueryLength {
            observed: query.len(),
        })?;
    if query.is_empty() {
        return Err(NarrowBandedError::EmptyQuery);
    }
    if max_distance > MAX_NARROW_BAND_DISTANCE {
        return Err(NarrowBandedError::Band {
            observed: max_distance,
        });
    }
    let band_length = max_distance * 2 + 1;
    let maximum = 32 / band_length;
    if output.is_empty() || output.len() > maximum {
        return Err(NarrowBandedError::PlacementBatch {
            observed: output.len(),
            maximum,
        });
    }
    let expected = pattern_len.saturating_mul(output.len());
    if patterns.len() != expected {
        return Err(NarrowBandedError::PatternDimension {
            expected,
            observed: patterns.len(),
        });
    }
    let dispatch = narrow_dispatch();
    (dispatch.placement_batch)(
        reference_masks_by_query,
        query,
        patterns,
        pattern_len,
        max_distance,
        output,
    );
    Ok(())
}

/// Computes distance-three frontiers for up to four candidates
/// while retaining only the seven active diagonal vectors.
///
/// # Errors
///
/// Rejects an empty query or inconsistent pattern and output dimensions.
pub fn narrow_banded_placement_distances_batch_d3(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    patterns: &[u8],
    output: &mut [NarrowPlacementDistances],
) -> Result<(), NarrowBandedError> {
    const BAND_LENGTH: usize = 7;
    const MAX_DISTANCE: usize = 3;

    let pattern_len =
        query
            .len()
            .checked_add(2 * MAX_DISTANCE)
            .ok_or(NarrowBandedError::QueryLength {
                observed: query.len(),
            })?;
    if query.is_empty() {
        return Err(NarrowBandedError::EmptyQuery);
    }
    let maximum = 32 / BAND_LENGTH;
    if output.is_empty() || output.len() > maximum {
        return Err(NarrowBandedError::PlacementBatch {
            observed: output.len(),
            maximum,
        });
    }
    let expected = pattern_len.saturating_mul(output.len());
    if patterns.len() != expected {
        return Err(NarrowBandedError::PatternDimension {
            expected,
            observed: patterns.len(),
        });
    }
    let dispatch = narrow_dispatch();
    (dispatch.placement_batch_d3)(
        reference_masks_by_query,
        query,
        patterns,
        pattern_len,
        MAX_DISTANCE,
        output,
    );
    Ok(())
}

/// Computes distance-three frontiers from a position-major, four-lane pattern
/// slab. Each reference position occupies four adjacent bytes; unused lanes
/// contain the sentinel reference code and are ignored.
///
/// This layout is private to the flexible verifier. Keeping its four reference
/// codes adjacent lets the AVX2 implementation classify every candidate with
/// one unaligned 32-bit load and one byte shuffle per diagonal.
///
/// # Errors
///
/// Rejects an empty query, a batch outside two through four candidates, or a
/// slab whose length is not exactly four times the candidate pattern length.
pub(crate) fn narrow_banded_placement_distances_interleaved_batch_d3(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    interleaved_patterns: &[u8],
    output: &mut [NarrowPlacementDistances],
) -> Result<(), NarrowBandedError> {
    const MAX_DISTANCE: usize = 3;
    const PATTERN_LANES: usize = 4;

    let pattern_len =
        query
            .len()
            .checked_add(2 * MAX_DISTANCE)
            .ok_or(NarrowBandedError::QueryLength {
                observed: query.len(),
            })?;
    if query.is_empty() {
        return Err(NarrowBandedError::EmptyQuery);
    }
    if !(2..=PATTERN_LANES).contains(&output.len()) {
        return Err(NarrowBandedError::PlacementBatch {
            observed: output.len(),
            maximum: PATTERN_LANES,
        });
    }
    let expected = pattern_len.saturating_mul(PATTERN_LANES);
    if interleaved_patterns.len() != expected {
        return Err(NarrowBandedError::PatternDimension {
            expected,
            observed: interleaved_patterns.len(),
        });
    }
    (narrow_dispatch().placement_interleaved_batch_d3)(
        reference_masks_by_query,
        query,
        interleaved_patterns,
        output,
    );
    Ok(())
}

/// Computes distance-five frontiers for up to two candidates while
/// retaining only the eleven active diagonal vectors per candidate.
///
/// # Errors
///
/// Rejects an empty query or inconsistent pattern and output dimensions.
///
/// # Panics
///
/// Panics only if the internally validated one-or-two-candidate batch cannot
/// be represented by its corresponding fixed-size SIMD view.
pub fn narrow_banded_placement_distances_batch_d5(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    patterns: &[u8],
    output: &mut [NarrowPlacementDistances],
) -> Result<(), NarrowBandedError> {
    const BAND_LENGTH: usize = 11;
    const MAX_DISTANCE: usize = 5;

    let pattern_len =
        query
            .len()
            .checked_add(2 * MAX_DISTANCE)
            .ok_or(NarrowBandedError::QueryLength {
                observed: query.len(),
            })?;
    if query.is_empty() {
        return Err(NarrowBandedError::EmptyQuery);
    }
    let maximum = 32 / BAND_LENGTH;
    if output.is_empty() || output.len() > maximum {
        return Err(NarrowBandedError::PlacementBatch {
            observed: output.len(),
            maximum,
        });
    }
    let expected = pattern_len.saturating_mul(output.len());
    if patterns.len() != expected {
        return Err(NarrowBandedError::PatternDimension {
            expected,
            observed: patterns.len(),
        });
    }
    let dispatch = narrow_dispatch();
    (dispatch.placement_batch_d5)(
        reference_masks_by_query,
        query,
        patterns,
        pattern_len,
        MAX_DISTANCE,
        output,
    );
    Ok(())
}

/// Finds one best in-band prefix for each fixed-width candidate pattern.
///
/// Every candidate occupies `query.len() + 2 * max_distance` consecutive
/// bytes in `patterns`. Codes `0..=4` are caller-defined symbols. For each
/// query code, `reference_masks_by_query[code]` names the reference codes with
/// zero substitution cost. Tied minima follow the frozen endpoint policy: prefer the later
/// prefix, except that an equal ungapped-center endpoint wins last.
///
/// # Errors
///
/// Validates the band and every flat dimension before writing output.
pub fn narrow_banded_prefix_batch(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    patterns: &[u8],
    max_distance: usize,
    output: &mut [NarrowBandedResult],
) -> Result<(), NarrowBandedError> {
    let pattern_length = validate_narrow_batch(query, patterns, max_distance, output)?;
    let dispatch = narrow_dispatch();
    (dispatch.prefix_batch)(
        *reference_masks_by_query,
        query,
        patterns,
        pattern_length,
        max_distance,
        output,
    );
    Ok(())
}

/// Computes every in-budget endpoint for fixed-start candidate patterns.
///
/// Each candidate pattern contains `query.len() + 2 * max_distance` symbols.
/// Its fixed alignment start is `max_distance` symbols from the beginning, so
/// endpoint delta `0..=2 * max_distance` represents interval lengths
/// `query.len() - max_distance ..= query.len() + max_distance`. Unlike
/// [`narrow_banded_prefix_batch`], this routine never minimizes over alternate
/// starts. Each SSE backend and `AArch64` NEON evaluate up to 16 independent
/// candidates in byte lanes; AVX2 evaluates up to 32 and AVX-512BW up to 64.
///
/// # Errors
///
/// Rejects an empty query, an unsupported band, or inconsistent flat pattern
/// and output dimensions before writing output.
pub fn narrow_banded_fixed_start_batch(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    patterns: &[u8],
    max_distance: usize,
    output: &mut [NarrowEndpointDistances],
) -> Result<(), NarrowBandedError> {
    let pattern_length =
        validate_narrow_fixed_start_batch(query, patterns, max_distance, output.len())?;
    let dispatch = narrow_dispatch();
    (dispatch.fixed_start_batch)(
        *reference_masks_by_query,
        query,
        patterns,
        pattern_length,
        max_distance,
        output,
    );
    Ok(())
}

/// Runs an explicitly selected fixed-start implementation for differential
/// qualification. `Scalar` is the portable reference implementation.
#[cfg(test)]
pub(super) fn narrow_banded_fixed_start_batch_with_flavor(
    flavor: NarrowBandedFlavor,
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    patterns: &[u8],
    max_distance: usize,
    output: &mut [NarrowEndpointDistances],
) -> Result<(), NarrowBandedError> {
    let pattern_length =
        validate_narrow_fixed_start_batch(query, patterns, max_distance, output.len())?;
    if matches!(flavor, NarrowBandedFlavor::Sse2) && !narrow_sse2_available()
        || matches!(flavor, NarrowBandedFlavor::Sse42) && !narrow_sse42_available()
        || matches!(flavor, NarrowBandedFlavor::Avx2) && !narrow_avx2_available()
        || matches!(flavor, NarrowBandedFlavor::Avx512) && !narrow_avx512_available()
    {
        return Err(NarrowBandedError::UnsupportedFlavor);
    }
    run_narrow_banded_fixed_start_batch(
        flavor,
        *reference_masks_by_query,
        query,
        patterns,
        pattern_length,
        max_distance,
        output,
    );
    Ok(())
}

fn run_narrow_banded_fixed_start_batch(
    flavor: NarrowBandedFlavor,
    reference_masks_by_query: [u8; 5],
    query: &[u8],
    patterns: &[u8],
    pattern_length: usize,
    max_distance: usize,
    output: &mut [NarrowEndpointDistances],
) {
    match flavor {
        NarrowBandedFlavor::Scalar => {
            for (pattern, result) in patterns.chunks_exact(pattern_length).zip(output) {
                *result = narrow_fixed_start_scalar_one(
                    &reference_masks_by_query,
                    query,
                    pattern,
                    max_distance,
                );
            }
        }
        NarrowBandedFlavor::Sse2 => {
            #[cfg(target_arch = "x86_64")]
            // SAFETY: runtime detection above proves SSE2; validation proves
            // every fixed-width pattern/output dimension.
            unsafe {
                x86_64::narrow_fixed_start_sse2(
                    &reference_masks_by_query,
                    query,
                    patterns,
                    pattern_length,
                    max_distance,
                    output,
                );
            }
            #[cfg(not(target_arch = "x86_64"))]
            unreachable!("non-x86_64 cannot report SSE2 available");
        }
        NarrowBandedFlavor::Sse42 => {
            #[cfg(target_arch = "x86_64")]
            // SAFETY: runtime detection above proves SSE4.1, SSE4.2, and
            // POPCNT; validation proves every fixed-width dimension.
            unsafe {
                x86_64::narrow_fixed_start_sse42(
                    &reference_masks_by_query,
                    query,
                    patterns,
                    pattern_length,
                    max_distance,
                    output,
                );
            }
            #[cfg(not(target_arch = "x86_64"))]
            unreachable!("non-x86_64 cannot report SSE4.2 available");
        }
        NarrowBandedFlavor::Avx2 => {
            #[cfg(target_arch = "x86_64")]
            // SAFETY: runtime detection above proves both AVX2 and POPCNT;
            // validation proves every fixed-width pattern/output dimension.
            unsafe {
                x86_64::narrow_fixed_start_avx2(
                    &reference_masks_by_query,
                    query,
                    patterns,
                    pattern_length,
                    max_distance,
                    output,
                );
            }
            #[cfg(not(target_arch = "x86_64"))]
            unreachable!("non-x86_64 cannot report AVX2 available");
        }
        NarrowBandedFlavor::Avx512 => {
            #[cfg(target_arch = "x86_64")]
            // SAFETY: runtime detection above proves AVX-512F, AVX-512BW,
            // AVX2, and POPCNT; validation proves all fixed-width dimensions.
            unsafe {
                x86_64_avx512::narrow_fixed_start_avx512(
                    &reference_masks_by_query,
                    query,
                    patterns,
                    pattern_length,
                    max_distance,
                    output,
                );
            }
            #[cfg(not(target_arch = "x86_64"))]
            unreachable!("non-x86_64 cannot report AVX-512 available");
        }
    }
}

/// Computes fixed-start endpoint frontiers without materializing candidate
/// patterns. All starts address one shared reference slice.
///
/// # Errors
///
/// Rejects an empty query, an unsupported band, mismatched start/output
/// dimensions, or a reference interval that cannot contain a candidate.
pub fn narrow_banded_fixed_start_gather_batch<T: NarrowReferenceCode>(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    reference: &[T],
    starts: &[usize],
    max_distance: usize,
    output: &mut [NarrowEndpointDistances],
) -> Result<(), NarrowBandedError> {
    if query.is_empty() {
        return Err(NarrowBandedError::EmptyQuery);
    }
    if max_distance > MAX_NARROW_BAND_DISTANCE {
        return Err(NarrowBandedError::Band {
            observed: max_distance,
        });
    }
    if starts.len() != output.len() {
        return Err(NarrowBandedError::OutputDimension {
            candidates: starts.len(),
            outputs: output.len(),
        });
    }
    let required_length = query.len().saturating_add(max_distance);
    for &start in starts {
        if start
            .checked_add(required_length)
            .is_none_or(|end| end > reference.len())
        {
            return Err(NarrowBandedError::PatternDimension {
                expected: start.saturating_add(required_length),
                observed: reference.len(),
            });
        }
    }
    if narrow_dispatch().gather_avx2 {
        #[cfg(target_arch = "x86_64")]
        // SAFETY: process initialization proved AVX2 and POPCNT support;
        // validation bounds all reference gathers.
        unsafe {
            x86_64::narrow_fixed_start_gather_avx2(
                reference_masks_by_query,
                query,
                reference,
                starts,
                max_distance,
                output,
            );
        }
        #[cfg(not(target_arch = "x86_64"))]
        unreachable!("non-x86_64 cannot select the AVX2 gather kernel");
    } else {
        for (&start, result) in starts.iter().zip(output) {
            *result = narrow_fixed_start_gather_scalar_one(
                reference_masks_by_query,
                query,
                reference,
                start,
                max_distance,
            );
        }
    }
    Ok(())
}

/// Computes complete fixed-start prefix-distance frontiers for two-word
/// queries, evaluating four independent candidates per AVX2 vector.
///
/// Each candidate pattern starts at its authoritative alignment coordinate and
/// contains `query_length + max_distance` symbols. Output delta
/// `0..=2 * max_distance` corresponds to prefix lengths
/// `query_length - max_distance ..= query_length + max_distance`. The supplied
/// equality masks use one query-position bit per reference code and may span
/// at most 128 query positions.
///
/// # Errors
///
/// Rejects an empty or overlong query, an unsupported distance, or inconsistent
/// flat pattern/output dimensions before writing output.
pub fn myers_prefix_distances_u128_batch(
    equality_masks: &[u128; 5],
    query_length: usize,
    patterns: &[u8],
    max_distance: usize,
    output: &mut [NarrowEndpointDistances],
) -> Result<(), NarrowBandedError> {
    if query_length == 0 {
        return Err(NarrowBandedError::EmptyQuery);
    }
    if query_length > u128::BITS as usize {
        return Err(NarrowBandedError::QueryLength {
            observed: query_length,
        });
    }
    if max_distance > MAX_NARROW_BAND_DISTANCE || max_distance >= query_length {
        return Err(NarrowBandedError::Band {
            observed: max_distance,
        });
    }
    let pattern_length = query_length.saturating_add(max_distance);
    let expected = pattern_length.saturating_mul(output.len());
    if patterns.len() != expected {
        return Err(NarrowBandedError::PatternDimension {
            expected,
            observed: patterns.len(),
        });
    }
    if narrow_dispatch().myers_u128_avx2 && query_length > u64::BITS as usize {
        #[cfg(target_arch = "x86_64")]
        // SAFETY: process initialization proved AVX2 and POPCNT support;
        // all dimensions and the two-word query bound were validated.
        unsafe {
            x86_64::myers_prefix_distances_u128_avx2(
                equality_masks,
                query_length,
                patterns,
                pattern_length,
                max_distance,
                output,
            );
        }
        #[cfg(not(target_arch = "x86_64"))]
        unreachable!("non-x86_64 cannot select the AVX2 u128 kernel");
    } else {
        for (pattern, result) in patterns.chunks_exact(pattern_length).zip(output) {
            *result = myers_prefix_distances_u128_scalar_one(
                equality_masks,
                query_length,
                pattern,
                max_distance,
            );
        }
    }
    Ok(())
}

/// Runs an explicitly selected narrow-band implementation for differential
/// qualification.
// This table is shared by scalar and SIMD dispatch; retaining one borrowed
// signature avoids copying it at every qualification seam.
#[cfg(test)]
pub(super) fn narrow_banded_prefix_batch_with_flavor(
    flavor: NarrowBandedFlavor,
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    patterns: &[u8],
    max_distance: usize,
    output: &mut [NarrowBandedResult],
) -> Result<(), NarrowBandedError> {
    let pattern_length = validate_narrow_batch(query, patterns, max_distance, output)?;
    if matches!(flavor, NarrowBandedFlavor::Sse2) && !narrow_sse2_available()
        || matches!(flavor, NarrowBandedFlavor::Sse42) && !narrow_sse42_available()
        || matches!(flavor, NarrowBandedFlavor::Avx2) && !narrow_avx2_available()
        || matches!(flavor, NarrowBandedFlavor::Avx512) && !narrow_avx512_available()
    {
        return Err(NarrowBandedError::UnsupportedFlavor);
    }
    run_narrow_banded_prefix_batch(
        flavor,
        *reference_masks_by_query,
        query,
        patterns,
        pattern_length,
        max_distance,
        output,
    );
    Ok(())
}

fn run_narrow_banded_prefix_batch(
    flavor: NarrowBandedFlavor,
    reference_masks_by_query: [u8; 5],
    query: &[u8],
    patterns: &[u8],
    pattern_length: usize,
    max_distance: usize,
    output: &mut [NarrowBandedResult],
) {
    match flavor {
        NarrowBandedFlavor::Scalar => {
            for (pattern, result) in patterns.chunks_exact(pattern_length).zip(output) {
                *result = narrow_banded_scalar_one(
                    &reference_masks_by_query,
                    query,
                    pattern,
                    max_distance,
                );
            }
        }
        NarrowBandedFlavor::Sse2 => {
            #[cfg(target_arch = "x86_64")]
            // SAFETY: runtime detection above proves SSE2 and validation
            // proves all fixed-width candidate/output dimensions.
            unsafe {
                x86_64::narrow_banded_sse2(
                    &reference_masks_by_query,
                    query,
                    patterns,
                    pattern_length,
                    max_distance,
                    output,
                );
            }
            #[cfg(not(target_arch = "x86_64"))]
            unreachable!("non-x86_64 cannot report SSE2 available");
        }
        NarrowBandedFlavor::Sse42 => {
            #[cfg(target_arch = "x86_64")]
            // SAFETY: runtime detection above proves SSE4.1, SSE4.2, and
            // POPCNT; validation proves all fixed-width dimensions.
            unsafe {
                x86_64::narrow_banded_sse42(
                    &reference_masks_by_query,
                    query,
                    patterns,
                    pattern_length,
                    max_distance,
                    output,
                );
            }
            #[cfg(not(target_arch = "x86_64"))]
            unreachable!("non-x86_64 cannot report SSE4.2 available");
        }
        NarrowBandedFlavor::Avx2 => {
            #[cfg(target_arch = "x86_64")]
            // SAFETY: runtime detection above proves both AVX2 and POPCNT;
            // validation proves all fixed-width candidate/output dimensions.
            unsafe {
                x86_64::narrow_banded_avx2(
                    &reference_masks_by_query,
                    query,
                    patterns,
                    pattern_length,
                    max_distance,
                    output,
                );
            }
            #[cfg(not(target_arch = "x86_64"))]
            unreachable!("non-x86_64 cannot report AVX2 available");
        }
        NarrowBandedFlavor::Avx512 => {
            #[cfg(target_arch = "x86_64")]
            // SAFETY: runtime detection above proves AVX-512F, AVX-512BW,
            // AVX2, and POPCNT; validation proves all fixed-width dimensions.
            unsafe {
                x86_64_avx512::narrow_banded_avx512(
                    &reference_masks_by_query,
                    query,
                    patterns,
                    pattern_length,
                    max_distance,
                    output,
                );
            }
            #[cfg(not(target_arch = "x86_64"))]
            unreachable!("non-x86_64 cannot report AVX-512 available");
        }
    }
}

fn validate_narrow_pattern(
    query: &[u8],
    pattern: &[u8],
    max_distance: usize,
) -> Result<usize, NarrowBandedError> {
    if query.is_empty() {
        return Err(NarrowBandedError::EmptyQuery);
    }
    if max_distance > MAX_NARROW_BAND_DISTANCE {
        return Err(NarrowBandedError::Band {
            observed: max_distance,
        });
    }
    let expected =
        query
            .len()
            .checked_add(max_distance.checked_mul(2).ok_or(
                NarrowBandedError::PatternDimension {
                    expected: usize::MAX,
                    observed: pattern.len(),
                },
            )?)
            .ok_or(NarrowBandedError::PatternDimension {
                expected: usize::MAX,
                observed: pattern.len(),
            })?;
    if pattern.len() != expected {
        return Err(NarrowBandedError::PatternDimension {
            expected,
            observed: pattern.len(),
        });
    }
    Ok(expected)
}

fn validate_narrow_batch(
    query: &[u8],
    patterns: &[u8],
    max_distance: usize,
    output: &[NarrowBandedResult],
) -> Result<usize, NarrowBandedError> {
    if query.is_empty() {
        return Err(NarrowBandedError::EmptyQuery);
    }
    if max_distance > MAX_NARROW_BAND_DISTANCE {
        return Err(NarrowBandedError::Band {
            observed: max_distance,
        });
    }
    let pattern_length =
        query
            .len()
            .checked_add(max_distance.checked_mul(2).ok_or(
                NarrowBandedError::PatternDimension {
                    expected: usize::MAX,
                    observed: patterns.len(),
                },
            )?)
            .ok_or(NarrowBandedError::PatternDimension {
                expected: usize::MAX,
                observed: patterns.len(),
            })?;
    let expected =
        pattern_length
            .checked_mul(output.len())
            .ok_or(NarrowBandedError::PatternDimension {
                expected: usize::MAX,
                observed: patterns.len(),
            })?;
    if patterns.len() != expected {
        return Err(NarrowBandedError::PatternDimension {
            expected,
            observed: patterns.len(),
        });
    }
    if pattern_length == 0 {
        return Err(NarrowBandedError::OutputDimension {
            candidates: 0,
            outputs: output.len(),
        });
    }
    Ok(pattern_length)
}

fn validate_narrow_fixed_start_batch(
    query: &[u8],
    patterns: &[u8],
    max_distance: usize,
    outputs: usize,
) -> Result<usize, NarrowBandedError> {
    let pattern_length = query.len().saturating_add(max_distance.saturating_mul(2));
    if query.is_empty() {
        return Err(NarrowBandedError::EmptyQuery);
    }
    if max_distance > MAX_NARROW_BAND_DISTANCE {
        return Err(NarrowBandedError::Band {
            observed: max_distance,
        });
    }
    let expected = pattern_length.saturating_mul(outputs);
    if patterns.len() != expected {
        return Err(NarrowBandedError::PatternDimension {
            expected,
            observed: patterns.len(),
        });
    }
    Ok(pattern_length)
}

// This table is shared by scalar and SIMD dispatch; retain one borrowed
// signature for the paired differential-test surface.
pub(super) fn narrow_banded_scalar_one(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    pattern: &[u8],
    max_distance: usize,
) -> NarrowBandedResult {
    let band_length = max_distance * 2 + 1;
    let high_bit = 1_u32 << (band_length - 1);
    let mut peq = [0_u32; 5];
    for (position, &reference_code) in pattern[..band_length].iter().enumerate() {
        if usize::from(reference_code) < peq.len() {
            peq[usize::from(reference_code)] |= 1_u32 << position;
        }
    }
    let mut positive = 0_u32;
    let mut negative = 0_u32;
    let mut error = 0_u32;
    for (query_position, &query_code) in query.iter().enumerate() {
        let mut equal = 0_u32;
        let reference_mask = reference_masks_by_query
            .get(usize::from(query_code))
            .copied()
            .unwrap_or(0);
        for (reference_code, &bits) in peq.iter().enumerate() {
            if reference_mask & (1_u8 << reference_code) != 0 {
                equal |= bits;
            }
        }
        let horizontal_input = equal | negative;
        let horizontal =
            (positive.wrapping_add(horizontal_input & positive) ^ positive) | horizontal_input;
        let negative_horizontal = positive & horizontal;
        let positive_horizontal = negative | !(positive | horizontal);
        let shifted = horizontal >> 1;
        negative = shifted & positive_horizontal;
        positive = negative_horizontal | !(shifted | positive_horizontal);
        error = error.wrapping_add(1_u32.wrapping_sub(horizontal & 1));
        if query_position + 1 != query.len() {
            for bits in &mut peq {
                *bits >>= 1;
            }
            let entering = pattern[band_length + query_position];
            if let Some(bits) = peq.get_mut(usize::from(entering)) {
                *bits |= high_bit;
            }
        }
    }
    finish_narrow_result(error, positive, negative, query.len(), max_distance)
}

fn finish_narrow_result(
    mut error: u32,
    positive: u32,
    negative: u32,
    query_length: usize,
    max_distance: usize,
) -> NarrowBandedResult {
    let mut best = NarrowBandedResult::ABSENT;
    let mut endpoint_distances = NarrowEndpointDistances::EMPTY;
    let mut center_error = error;
    for delta in 0..=max_distance * 2 {
        if delta != 0 {
            let bit = delta - 1;
            error = error
                .wrapping_add((positive >> bit) & 1)
                .wrapping_sub((negative >> bit) & 1);
        }
        if delta == max_distance {
            center_error = error;
        }
        if usize::try_from(error).unwrap_or(usize::MAX) <= max_distance {
            endpoint_distances.insert(delta, error);
            if error < best.distance {
                best = NarrowBandedResult {
                    distance: error,
                    prefix_length: query_length + delta,
                    tied_prefix_mask: 1_u32 << delta,
                    endpoint_distances: NarrowEndpointDistances::EMPTY,
                };
            } else if error == best.distance {
                best.prefix_length = query_length + delta;
                best.tied_prefix_mask |= 1_u32 << delta;
            }
        }
    }
    if center_error <= best.distance
        && usize::try_from(center_error).unwrap_or(usize::MAX) <= max_distance
    {
        best = NarrowBandedResult {
            distance: center_error,
            prefix_length: query_length + max_distance,
            tied_prefix_mask: best.tied_prefix_mask,
            endpoint_distances: NarrowEndpointDistances::EMPTY,
        };
    }
    best.endpoint_distances = endpoint_distances;
    best
}

// This table is shared by scalar and SIMD dispatch; retain one borrowed
// signature for the paired differential-test surface.
pub(super) fn narrow_placement_distances_scalar(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    pattern: &[u8],
    max_distance: usize,
) -> NarrowPlacementDistances {
    let band_length = max_distance * 2 + 1;
    let capped = u8::try_from(max_distance)
        .unwrap_or(u8::MAX - 1)
        .saturating_add(1);
    let mut distances = NarrowPlacementDistances::for_band(band_length, max_distance);
    let mut previous = [capped; MAX_NARROW_ENDPOINTS];
    let mut current = [capped; MAX_NARROW_ENDPOINTS];
    for start_delta in 0..band_length {
        previous[..band_length].fill(capped);
        previous[start_delta] = 0;
        for diagonal in start_delta + 1..band_length {
            previous[diagonal] = previous[diagonal - 1].saturating_add(1).min(capped);
        }
        for (query_position, &query_code) in query.iter().enumerate() {
            current[..band_length].fill(capped);
            let reference_mask = reference_masks_by_query
                .get(usize::from(query_code))
                .copied()
                .unwrap_or(0);
            for diagonal in 0..band_length {
                let reference_position = query_position + diagonal;
                let substitution = u8::from(
                    reference_mask & (1_u8 << usize::from(pattern[reference_position])) == 0,
                );
                let mut best = previous[diagonal].saturating_add(substitution).min(capped);
                if diagonal + 1 < band_length {
                    best = best.min(previous[diagonal + 1].saturating_add(1).min(capped));
                }
                if diagonal != 0 {
                    best = best.min(current[diagonal - 1].saturating_add(1).min(capped));
                }
                current[diagonal] = best;
            }
            core::mem::swap(&mut previous, &mut current);
        }
        for (endpoint_delta, &distance) in previous[..band_length].iter().enumerate() {
            distances.insert_distance(start_delta, endpoint_delta, distance);
        }
    }
    distances
}

pub(super) fn narrow_placement_distances_scalar_interleaved(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    patterns: &[u8],
    candidate: usize,
    pattern_lanes: usize,
    max_distance: usize,
) -> NarrowPlacementDistances {
    let band_length = max_distance * 2 + 1;
    let capped = u8::try_from(max_distance)
        .unwrap_or(u8::MAX - 1)
        .saturating_add(1);
    let mut distances = NarrowPlacementDistances::for_band(band_length, max_distance);
    let mut previous = [capped; MAX_NARROW_ENDPOINTS];
    let mut current = [capped; MAX_NARROW_ENDPOINTS];
    for start_delta in 0..band_length {
        previous[..band_length].fill(capped);
        previous[start_delta] = 0;
        for diagonal in start_delta + 1..band_length {
            previous[diagonal] = previous[diagonal - 1].saturating_add(1).min(capped);
        }
        for (query_position, &query_code) in query.iter().enumerate() {
            current[..band_length].fill(capped);
            let reference_mask = reference_masks_by_query
                .get(usize::from(query_code))
                .copied()
                .unwrap_or(0);
            for diagonal in 0..band_length {
                let reference_position = query_position + diagonal;
                let reference_code = patterns[reference_position * pattern_lanes + candidate];
                let substitution =
                    u8::from(reference_mask & (1_u8 << usize::from(reference_code)) == 0);
                let mut best = previous[diagonal].saturating_add(substitution).min(capped);
                if diagonal + 1 < band_length {
                    best = best.min(previous[diagonal + 1].saturating_add(1).min(capped));
                }
                if diagonal != 0 {
                    best = best.min(current[diagonal - 1].saturating_add(1).min(capped));
                }
                current[diagonal] = best;
            }
            core::mem::swap(&mut previous, &mut current);
        }
        for (endpoint_delta, &distance) in previous[..band_length].iter().enumerate() {
            distances.insert_distance(start_delta, endpoint_delta, distance);
        }
    }
    distances
}

// Keep the relation table borrowed to match its SIMD counterpart exactly.
pub(super) fn narrow_fixed_start_scalar_one(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    pattern: &[u8],
    max_distance: usize,
) -> NarrowEndpointDistances {
    let placements =
        narrow_placement_distances_scalar(reference_masks_by_query, query, pattern, max_distance);
    let mut result = NarrowEndpointDistances::EMPTY;
    for endpoint in 0..=max_distance * 2 {
        if let Some(distance) = placements.distance(max_distance, endpoint) {
            result.insert(endpoint, distance);
        }
    }
    result
}

// Keep the relation table borrowed to match its SIMD counterpart exactly.
fn narrow_fixed_start_gather_scalar_one<T: NarrowReferenceCode>(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    reference: &[T],
    start: usize,
    max_distance: usize,
) -> NarrowEndpointDistances {
    let band_length = max_distance * 2 + 1;
    let capped = u8::try_from(max_distance)
        .unwrap_or(u8::MAX - 1)
        .saturating_add(1);
    let mut previous = [capped; MAX_NARROW_ENDPOINTS];
    let mut current = [capped; MAX_NARROW_ENDPOINTS];
    for (diagonal, slot) in previous[..band_length].iter_mut().enumerate() {
        *slot = if diagonal < max_distance {
            capped
        } else {
            u8::try_from(diagonal - max_distance)
                .unwrap_or(capped)
                .min(capped)
        };
    }
    for (query_position, &query_code) in query.iter().enumerate() {
        current[..band_length].fill(capped);
        let reference_mask = reference_masks_by_query
            .get(usize::from(query_code))
            .copied()
            .unwrap_or(0);
        for diagonal in 0..band_length {
            let shifted = query_position + diagonal;
            let reference_code = if shifted < max_distance {
                u8::MAX
            } else {
                reference[start + shifted - max_distance].narrow_reference_code()
            };
            let substitution =
                u8::from(reference_code >= 8 || reference_mask & (1_u8 << reference_code) == 0);
            let mut best = previous[diagonal].saturating_add(substitution).min(capped);
            if diagonal + 1 < band_length {
                best = best.min(previous[diagonal + 1].saturating_add(1).min(capped));
            }
            if diagonal != 0 {
                best = best.min(current[diagonal - 1].saturating_add(1).min(capped));
            }
            current[diagonal] = best;
        }
        core::mem::swap(&mut previous, &mut current);
    }
    let mut result = NarrowEndpointDistances::EMPTY;
    for (endpoint, &distance) in previous[..band_length].iter().enumerate() {
        if usize::from(distance) <= max_distance {
            result.insert(endpoint, u32::from(distance));
        }
    }
    result
}

// Both conversions are guarded by the validated 128-base query and
// distance-at-most-fifteen kernel domain.
#[allow(clippy::cast_possible_truncation)]
pub(super) fn myers_prefix_distances_u128_scalar_one(
    equality_masks: &[u128; 5],
    query_length: usize,
    pattern: &[u8],
    max_distance: usize,
) -> NarrowEndpointDistances {
    let mut positive = !0_u128;
    let mut negative = 0_u128;
    let mut score = u64::try_from(query_length).expect("query length fits u64");
    let high_bit = 1_u128 << (query_length - 1);
    let minimum_end = query_length - max_distance;
    let mut result = NarrowEndpointDistances::EMPTY;
    for (position, &code) in pattern.iter().enumerate() {
        let equal = match code {
            0 => equality_masks[0],
            1 => equality_masks[1],
            2 => equality_masks[2],
            3 => equality_masks[3],
            _ => 0,
        };
        let horizontal_input = equal | negative;
        let horizontal = (((equal & positive).wrapping_add(positive)) ^ positive) | equal;
        let positive_horizontal = negative | !(horizontal | positive);
        let negative_horizontal = positive & horizontal;
        if positive_horizontal & high_bit != 0 {
            score = score.saturating_add(1);
        } else if negative_horizontal & high_bit != 0 {
            score = score.saturating_sub(1);
        }
        let shifted_positive = (positive_horizontal << 1) | 1;
        let shifted_negative = negative_horizontal << 1;
        positive = shifted_negative | !(horizontal_input | shifted_positive);
        negative = shifted_positive & horizontal_input;
        let prefix_length = position + 1;
        if prefix_length >= minimum_end && score <= max_distance as u64 {
            result.insert(prefix_length - minimum_end, score as u32);
        }
    }
    result
}

#[cfg(target_arch = "x86_64")]
#[cfg(test)]
pub(super) fn narrow_sse2_available() -> bool {
    configuration().features().sse2()
}

#[cfg(not(target_arch = "x86_64"))]
#[cfg(test)]
pub(super) const fn narrow_sse2_available() -> bool {
    false
}

#[cfg(target_arch = "x86_64")]
#[cfg(test)]
pub(super) fn narrow_sse42_available() -> bool {
    Backend::Sse42Popcnt.is_supported_by(configuration().features())
}

#[cfg(not(target_arch = "x86_64"))]
#[cfg(test)]
pub(super) const fn narrow_sse42_available() -> bool {
    false
}

#[cfg(target_arch = "x86_64")]
#[cfg(test)]
pub(super) fn narrow_avx2_available() -> bool {
    let features = configuration().features();
    features.avx2() && features.popcnt()
}

#[cfg(not(target_arch = "x86_64"))]
#[cfg(test)]
pub(super) const fn narrow_avx2_available() -> bool {
    false
}

#[cfg(target_arch = "x86_64")]
#[cfg(test)]
pub(super) fn narrow_avx512_available() -> bool {
    Backend::Avx512BwPopcnt.is_supported_by(configuration().features())
}

#[cfg(not(target_arch = "x86_64"))]
#[cfg(test)]
pub(super) const fn narrow_avx512_available() -> bool {
    false
}

fn placement_scalar_dispatch(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    pattern: &[u8],
    max_distance: usize,
) -> NarrowPlacementDistances {
    narrow_placement_distances_scalar(reference_masks_by_query, query, pattern, max_distance)
}

fn placement_d3_scalar_dispatch(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    pattern: &[u8],
) -> NarrowPlacementDistances {
    narrow_placement_distances_scalar(reference_masks_by_query, query, pattern, 3)
}

fn placement_d5_scalar_dispatch(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    pattern: &[u8],
) -> NarrowPlacementDistances {
    narrow_placement_distances_scalar(reference_masks_by_query, query, pattern, 5)
}

fn placement_batch_scalar_dispatch(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    patterns: &[u8],
    pattern_length: usize,
    max_distance: usize,
    output: &mut [NarrowPlacementDistances],
) {
    for (pattern, destination) in patterns.chunks_exact(pattern_length).zip(output) {
        *destination = narrow_placement_distances_scalar(
            reference_masks_by_query,
            query,
            pattern,
            max_distance,
        );
    }
}

fn placement_interleaved_batch_d3_scalar_dispatch(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    interleaved_patterns: &[u8],
    output: &mut [NarrowPlacementDistances],
) {
    const PATTERN_LANES: usize = 4;
    const MAX_DISTANCE: usize = 3;

    for (candidate, destination) in output.iter_mut().enumerate() {
        *destination = narrow_placement_distances_scalar_interleaved(
            reference_masks_by_query,
            query,
            interleaved_patterns,
            candidate,
            PATTERN_LANES,
            MAX_DISTANCE,
        );
    }
}

fn prefix_batch_scalar_dispatch(
    reference_masks_by_query: [u8; 5],
    query: &[u8],
    patterns: &[u8],
    pattern_length: usize,
    max_distance: usize,
    output: &mut [NarrowBandedResult],
) {
    run_narrow_banded_prefix_batch(
        NarrowBandedFlavor::Scalar,
        reference_masks_by_query,
        query,
        patterns,
        pattern_length,
        max_distance,
        output,
    );
}

fn fixed_start_batch_scalar_dispatch(
    reference_masks_by_query: [u8; 5],
    query: &[u8],
    patterns: &[u8],
    pattern_length: usize,
    max_distance: usize,
    output: &mut [NarrowEndpointDistances],
) {
    run_narrow_banded_fixed_start_batch(
        NarrowBandedFlavor::Scalar,
        reference_masks_by_query,
        query,
        patterns,
        pattern_length,
        max_distance,
        output,
    );
}

#[cfg(target_arch = "aarch64")]
fn placement_neon_dispatch(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    pattern: &[u8],
    max_distance: usize,
) -> NarrowPlacementDistances {
    // SAFETY: AArch64 process initialization validates NEON support before
    // this private trampoline enters the dispatch table.
    unsafe {
        neon::narrow_placement_distances(reference_masks_by_query, query, pattern, max_distance)
    }
}

#[cfg(target_arch = "aarch64")]
fn placement_d3_neon_dispatch(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    pattern: &[u8],
) -> NarrowPlacementDistances {
    placement_neon_dispatch(reference_masks_by_query, query, pattern, 3)
}

#[cfg(target_arch = "aarch64")]
fn placement_d5_neon_dispatch(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    pattern: &[u8],
) -> NarrowPlacementDistances {
    placement_neon_dispatch(reference_masks_by_query, query, pattern, 5)
}

#[cfg(target_arch = "aarch64")]
fn placement_batch_neon_dispatch(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    patterns: &[u8],
    pattern_length: usize,
    max_distance: usize,
    output: &mut [NarrowPlacementDistances],
) {
    for (pattern, destination) in patterns.chunks_exact(pattern_length).zip(output) {
        *destination =
            placement_neon_dispatch(reference_masks_by_query, query, pattern, max_distance);
    }
}

#[cfg(target_arch = "aarch64")]
fn prefix_batch_neon_dispatch(
    reference_masks_by_query: [u8; 5],
    query: &[u8],
    patterns: &[u8],
    pattern_length: usize,
    max_distance: usize,
    output: &mut [NarrowBandedResult],
) {
    // SAFETY: installed only for the AArch64 NEON process backend; callers
    // validated every flattened pattern and output dimension.
    unsafe {
        neon::narrow_banded_prefix_batch(
            &reference_masks_by_query,
            query,
            patterns,
            pattern_length,
            max_distance,
            output,
        );
    }
}

#[cfg(target_arch = "aarch64")]
fn fixed_start_batch_neon_dispatch(
    reference_masks_by_query: [u8; 5],
    query: &[u8],
    patterns: &[u8],
    pattern_length: usize,
    max_distance: usize,
    output: &mut [NarrowEndpointDistances],
) {
    // SAFETY: installed only for the AArch64 NEON process backend; callers
    // validated every flattened pattern and output dimension.
    unsafe {
        neon::narrow_fixed_start_batch(
            &reference_masks_by_query,
            query,
            patterns,
            pattern_length,
            max_distance,
            output,
        );
    }
}

#[cfg(target_arch = "x86_64")]
fn placement_sse2_dispatch(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    pattern: &[u8],
    max_distance: usize,
) -> NarrowPlacementDistances {
    // SAFETY: the validated SSE2 backend installs this private trampoline.
    unsafe {
        x86_64::narrow_placement_distances_sse2(
            reference_masks_by_query,
            query,
            pattern,
            max_distance,
        )
    }
}

#[cfg(target_arch = "x86_64")]
fn placement_d3_sse2_dispatch(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    pattern: &[u8],
) -> NarrowPlacementDistances {
    placement_sse2_dispatch(reference_masks_by_query, query, pattern, 3)
}

#[cfg(target_arch = "x86_64")]
fn placement_d5_sse2_dispatch(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    pattern: &[u8],
) -> NarrowPlacementDistances {
    placement_sse2_dispatch(reference_masks_by_query, query, pattern, 5)
}

#[cfg(target_arch = "x86_64")]
fn placement_batch_sse2_dispatch(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    patterns: &[u8],
    pattern_length: usize,
    max_distance: usize,
    output: &mut [NarrowPlacementDistances],
) {
    for (pattern, destination) in patterns.chunks_exact(pattern_length).zip(output) {
        *destination =
            placement_sse2_dispatch(reference_masks_by_query, query, pattern, max_distance);
    }
}

#[cfg(target_arch = "x86_64")]
fn prefix_batch_sse2_dispatch(
    reference_masks_by_query: [u8; 5],
    query: &[u8],
    patterns: &[u8],
    pattern_length: usize,
    max_distance: usize,
    output: &mut [NarrowBandedResult],
) {
    run_narrow_banded_prefix_batch(
        NarrowBandedFlavor::Sse2,
        reference_masks_by_query,
        query,
        patterns,
        pattern_length,
        max_distance,
        output,
    );
}

#[cfg(target_arch = "x86_64")]
fn fixed_start_batch_sse2_dispatch(
    reference_masks_by_query: [u8; 5],
    query: &[u8],
    patterns: &[u8],
    pattern_length: usize,
    max_distance: usize,
    output: &mut [NarrowEndpointDistances],
) {
    run_narrow_banded_fixed_start_batch(
        NarrowBandedFlavor::Sse2,
        reference_masks_by_query,
        query,
        patterns,
        pattern_length,
        max_distance,
        output,
    );
}

#[cfg(target_arch = "x86_64")]
fn placement_sse42_dispatch(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    pattern: &[u8],
    max_distance: usize,
) -> NarrowPlacementDistances {
    // SAFETY: the validated SSE4.2 backend installs this private trampoline.
    unsafe {
        x86_64::narrow_placement_distances_sse42(
            reference_masks_by_query,
            query,
            pattern,
            max_distance,
        )
    }
}

#[cfg(target_arch = "x86_64")]
fn placement_d3_sse42_dispatch(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    pattern: &[u8],
) -> NarrowPlacementDistances {
    placement_sse42_dispatch(reference_masks_by_query, query, pattern, 3)
}

#[cfg(target_arch = "x86_64")]
fn placement_d5_sse42_dispatch(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    pattern: &[u8],
) -> NarrowPlacementDistances {
    placement_sse42_dispatch(reference_masks_by_query, query, pattern, 5)
}

#[cfg(target_arch = "x86_64")]
fn placement_batch_sse42_dispatch(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    patterns: &[u8],
    pattern_length: usize,
    max_distance: usize,
    output: &mut [NarrowPlacementDistances],
) {
    for (pattern, destination) in patterns.chunks_exact(pattern_length).zip(output) {
        *destination =
            placement_sse42_dispatch(reference_masks_by_query, query, pattern, max_distance);
    }
}

#[cfg(target_arch = "x86_64")]
fn prefix_batch_sse42_dispatch(
    reference_masks_by_query: [u8; 5],
    query: &[u8],
    patterns: &[u8],
    pattern_length: usize,
    max_distance: usize,
    output: &mut [NarrowBandedResult],
) {
    run_narrow_banded_prefix_batch(
        NarrowBandedFlavor::Sse42,
        reference_masks_by_query,
        query,
        patterns,
        pattern_length,
        max_distance,
        output,
    );
}

#[cfg(target_arch = "x86_64")]
fn fixed_start_batch_sse42_dispatch(
    reference_masks_by_query: [u8; 5],
    query: &[u8],
    patterns: &[u8],
    pattern_length: usize,
    max_distance: usize,
    output: &mut [NarrowEndpointDistances],
) {
    run_narrow_banded_fixed_start_batch(
        NarrowBandedFlavor::Sse42,
        reference_masks_by_query,
        query,
        patterns,
        pattern_length,
        max_distance,
        output,
    );
}

#[cfg(target_arch = "x86_64")]
fn placement_avx2_dispatch(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    pattern: &[u8],
    max_distance: usize,
) -> NarrowPlacementDistances {
    // SAFETY: this private function can enter the dispatch table only for the
    // validated AVX2+POPCNT process backend.
    unsafe {
        x86_64::narrow_placement_distances_avx2(
            reference_masks_by_query,
            query,
            pattern,
            max_distance,
        )
    }
}

#[cfg(target_arch = "x86_64")]
fn placement_d3_avx2_dispatch(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    pattern: &[u8],
) -> NarrowPlacementDistances {
    // SAFETY: installed only by the validated AVX2+POPCNT dispatch table.
    unsafe { x86_64::narrow_placement_distances_d3_avx2(reference_masks_by_query, query, pattern) }
}

#[cfg(target_arch = "x86_64")]
fn placement_d5_avx2_dispatch(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    pattern: &[u8],
) -> NarrowPlacementDistances {
    // SAFETY: installed only by the validated AVX2+POPCNT dispatch table.
    unsafe { x86_64::narrow_placement_distances_d5_avx2(reference_masks_by_query, query, pattern) }
}

#[cfg(target_arch = "x86_64")]
fn placement_batch_avx2_dispatch(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    patterns: &[u8],
    pattern_length: usize,
    max_distance: usize,
    output: &mut [NarrowPlacementDistances],
) {
    // SAFETY: installed only by the validated AVX2+POPCNT dispatch table;
    // safe callers validate all flat dimensions before invoking this pointer.
    unsafe {
        x86_64::narrow_placement_distances_batch_avx2(
            reference_masks_by_query,
            query,
            patterns,
            pattern_length,
            max_distance,
            output,
        );
    }
}

#[cfg(target_arch = "x86_64")]
fn placement_batch_d3_avx2_dispatch(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    patterns: &[u8],
    pattern_length: usize,
    max_distance: usize,
    output: &mut [NarrowPlacementDistances],
) {
    debug_assert_eq!(max_distance, 3);
    // SAFETY: installed only by the validated AVX2+POPCNT dispatch table;
    // the public distance-three wrapper validates all dimensions.
    unsafe {
        x86_64::narrow_placement_distances_batch_d3_avx2(
            reference_masks_by_query,
            query,
            patterns,
            pattern_length,
            output,
        );
    }
}

#[cfg(target_arch = "x86_64")]
fn placement_interleaved_batch_d3_avx2_dispatch(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    interleaved_patterns: &[u8],
    output: &mut [NarrowPlacementDistances],
) {
    // SAFETY: installed only by the validated AVX2+POPCNT dispatch table;
    // the private wrapper proves a complete four-byte group exists for every
    // diagonal access and limits the active candidate count to four.
    unsafe {
        x86_64::narrow_placement_distances_interleaved_batch_d3_avx2(
            reference_masks_by_query,
            query,
            interleaved_patterns,
            output,
        );
    }
}

#[cfg(target_arch = "x86_64")]
fn placement_batch_d5_avx2_dispatch(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    patterns: &[u8],
    pattern_length: usize,
    max_distance: usize,
    output: &mut [NarrowPlacementDistances],
) {
    debug_assert_eq!(max_distance, 5);
    // SAFETY: installed only by the validated AVX2+POPCNT dispatch table;
    // validation limits this batch to one or two complete candidates.
    unsafe {
        if let [destination] = output {
            *destination = x86_64::narrow_placement_distances_d5_avx2(
                reference_masks_by_query,
                query,
                patterns,
            );
        } else if output.len() == 2 {
            let pair: &mut [NarrowPlacementDistances; 2] = output
                .try_into()
                .expect("validated distance-five SIMD batch has two candidates");
            x86_64::narrow_placement_distances_batch_d5_avx2(
                reference_masks_by_query,
                query,
                patterns,
                pattern_length,
                pair,
            );
        } else {
            unreachable!("validated distance-five batch has one or two candidates");
        }
    }
}

#[cfg(target_arch = "x86_64")]
fn prefix_batch_avx2_dispatch(
    reference_masks_by_query: [u8; 5],
    query: &[u8],
    patterns: &[u8],
    pattern_length: usize,
    max_distance: usize,
    output: &mut [NarrowBandedResult],
) {
    run_narrow_banded_prefix_batch(
        NarrowBandedFlavor::Avx2,
        reference_masks_by_query,
        query,
        patterns,
        pattern_length,
        max_distance,
        output,
    );
}

#[cfg(target_arch = "x86_64")]
fn fixed_start_batch_avx2_dispatch(
    reference_masks_by_query: [u8; 5],
    query: &[u8],
    patterns: &[u8],
    pattern_length: usize,
    max_distance: usize,
    output: &mut [NarrowEndpointDistances],
) {
    run_narrow_banded_fixed_start_batch(
        NarrowBandedFlavor::Avx2,
        reference_masks_by_query,
        query,
        patterns,
        pattern_length,
        max_distance,
        output,
    );
}

#[cfg(target_arch = "x86_64")]
fn prefix_batch_avx512_dispatch(
    reference_masks_by_query: [u8; 5],
    query: &[u8],
    patterns: &[u8],
    pattern_length: usize,
    max_distance: usize,
    output: &mut [NarrowBandedResult],
) {
    run_narrow_banded_prefix_batch(
        NarrowBandedFlavor::Avx512,
        reference_masks_by_query,
        query,
        patterns,
        pattern_length,
        max_distance,
        output,
    );
}

#[cfg(target_arch = "x86_64")]
fn fixed_start_batch_avx512_dispatch(
    reference_masks_by_query: [u8; 5],
    query: &[u8],
    patterns: &[u8],
    pattern_length: usize,
    max_distance: usize,
    output: &mut [NarrowEndpointDistances],
) {
    run_narrow_banded_fixed_start_batch(
        NarrowBandedFlavor::Avx512,
        reference_masks_by_query,
        query,
        patterns,
        pattern_length,
        max_distance,
        output,
    );
}
