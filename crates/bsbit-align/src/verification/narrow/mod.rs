//! Narrow-band scalar and architecture-dispatched alignment kernels.

use super::{MAX_NARROW_BAND_DISTANCE, NarrowReferenceCode};
use core::fmt;

mod scalar;
mod x86;

pub(super) use self::scalar::{
    myers_prefix_distances_u128_scalar_one, narrow_banded_scalar_one,
    narrow_fixed_start_scalar_one, narrow_placement_distances_scalar,
};
use self::scalar::{
    narrow_fixed_start_gather_scalar_one, narrow_placement_distances_scalar_interleaved,
};
#[cfg(target_arch = "x86_64")]
use self::x86::{
    myers_prefix_distances_u128_avx2, narrow_banded_avx2, narrow_banded_sse42,
    narrow_fixed_start_avx2, narrow_fixed_start_gather_avx2, narrow_placement_distances_batch_avx2,
    narrow_placement_distances_batch_d3_avx2, narrow_placement_distances_batch_d5_avx2,
    narrow_placement_distances_d3_avx2, narrow_placement_distances_d5_avx2,
    narrow_placement_distances_interleaved_batch_d3_avx2,
};
pub(super) use self::x86::{narrow_avx2_available, narrow_sse42_available};
#[cfg(target_arch = "x86_64")]
pub(super) use self::x86::{
    narrow_fixed_start_sse42, narrow_placement_distances_avx2, narrow_placement_distances_sse42,
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

/// Runtime implementation used by the 32-bit narrow-band candidate kernel.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum NarrowBandedFlavor {
    /// Portable one-candidate 32-bit implementation.
    Scalar,
    /// Four independent candidate patterns in SSE4.2 32-bit lanes.
    Sse42,
    /// Eight independent candidate patterns in AVX2 32-bit lanes.
    Avx2,
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
    /// A forced SIMD implementation is unavailable.
    UnsupportedFlavor {
        /// Requested runtime implementation.
        flavor: NarrowBandedFlavor,
    },
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
            Self::UnsupportedFlavor { flavor } => {
                write!(formatter, "narrow-band kernel {flavor:?} is unavailable")
            }
            Self::PlacementBatch { observed, maximum } => write!(
                formatter,
                "narrow placement batch {observed} exceeds SIMD capacity {maximum}"
            ),
        }
    }
}

impl std::error::Error for NarrowBandedError {}

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
    #[cfg(target_arch = "x86_64")]
    if narrow_avx2_available() {
        // SAFETY: runtime detection proves AVX2 support, and validation above
        // bounds every pattern and matrix access.
        return Ok(unsafe {
            narrow_placement_distances_avx2(reference_masks_by_query, query, pattern, max_distance)
        });
    }
    #[cfg(target_arch = "x86_64")]
    if narrow_sse42_available() {
        // SAFETY: runtime detection proves SSE4.2 support, and validation above
        // bounds every pattern and matrix access.
        return Ok(unsafe {
            narrow_placement_distances_sse42(reference_masks_by_query, query, pattern, max_distance)
        });
    }
    Ok(narrow_placement_distances_scalar(
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
    #[cfg(target_arch = "x86_64")]
    if narrow_avx2_available() {
        // SAFETY: runtime detection proves AVX2 support, and validation above
        // bounds every pattern and matrix access.
        return Ok(unsafe {
            narrow_placement_distances_d3_avx2(reference_masks_by_query, query, pattern)
        });
    }
    #[cfg(target_arch = "x86_64")]
    if narrow_sse42_available() {
        // SAFETY: runtime detection proves SSE4.2 support, and validation above
        // bounds every pattern and matrix access.
        return Ok(unsafe {
            narrow_placement_distances_sse42(reference_masks_by_query, query, pattern, MAX_DISTANCE)
        });
    }
    Ok(narrow_placement_distances_scalar(
        reference_masks_by_query,
        query,
        pattern,
        MAX_DISTANCE,
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
    #[cfg(target_arch = "x86_64")]
    if narrow_avx2_available() {
        // SAFETY: runtime detection proves AVX2 support, and validation above
        // bounds every pattern and matrix access.
        return Ok(unsafe {
            narrow_placement_distances_d5_avx2(reference_masks_by_query, query, pattern)
        });
    }
    #[cfg(target_arch = "x86_64")]
    if narrow_sse42_available() {
        // SAFETY: runtime detection proves SSE4.2 support, and validation above
        // bounds every pattern and matrix access.
        return Ok(unsafe {
            narrow_placement_distances_sse42(reference_masks_by_query, query, pattern, MAX_DISTANCE)
        });
    }
    Ok(narrow_placement_distances_scalar(
        reference_masks_by_query,
        query,
        pattern,
        MAX_DISTANCE,
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
) -> Result<NarrowBandedFlavor, NarrowBandedError> {
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
    if narrow_avx2_available() {
        #[cfg(target_arch = "x86_64")]
        // SAFETY: runtime detection proves AVX2 and dimensions above keep all
        // active candidate/start lanes inside one vector.
        unsafe {
            narrow_placement_distances_batch_avx2(
                reference_masks_by_query,
                query,
                patterns,
                pattern_len,
                max_distance,
                output,
            );
        }
        return Ok(NarrowBandedFlavor::Avx2);
    }
    if narrow_sse42_available() {
        #[cfg(target_arch = "x86_64")]
        for (pattern, destination) in patterns.chunks_exact(pattern_len).zip(output) {
            // SAFETY: runtime detection proves SSE4.2 and the batch validation
            // proves every fixed-width candidate dimension.
            *destination = unsafe {
                narrow_placement_distances_sse42(
                    reference_masks_by_query,
                    query,
                    pattern,
                    max_distance,
                )
            };
        }
        return Ok(NarrowBandedFlavor::Sse42);
    }
    for (pattern, destination) in patterns.chunks_exact(pattern_len).zip(output) {
        *destination = narrow_placement_distances_scalar(
            reference_masks_by_query,
            query,
            pattern,
            max_distance,
        );
    }
    Ok(NarrowBandedFlavor::Scalar)
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
) -> Result<NarrowBandedFlavor, NarrowBandedError> {
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
    if narrow_avx2_available() {
        #[cfg(target_arch = "x86_64")]
        // SAFETY: runtime detection proves AVX2 and dimensions above keep all
        // active candidate/start lanes inside one vector.
        unsafe {
            narrow_placement_distances_batch_d3_avx2(
                reference_masks_by_query,
                query,
                patterns,
                pattern_len,
                output,
            );
        }
        return Ok(NarrowBandedFlavor::Avx2);
    }
    if narrow_sse42_available() {
        #[cfg(target_arch = "x86_64")]
        for (pattern, destination) in patterns.chunks_exact(pattern_len).zip(output) {
            // SAFETY: runtime detection proves SSE4.2 and the batch validation
            // proves every fixed distance-three candidate dimension.
            *destination = unsafe {
                narrow_placement_distances_sse42(
                    reference_masks_by_query,
                    query,
                    pattern,
                    MAX_DISTANCE,
                )
            };
        }
        return Ok(NarrowBandedFlavor::Sse42);
    }
    for (pattern, destination) in patterns.chunks_exact(pattern_len).zip(output) {
        *destination = narrow_placement_distances_scalar(
            reference_masks_by_query,
            query,
            pattern,
            MAX_DISTANCE,
        );
    }
    Ok(NarrowBandedFlavor::Scalar)
}

/// Computes distance-three frontiers from a position-major, four-lane pattern
/// slab. Each reference position occupies four adjacent bytes; unused lanes
/// must contain any valid reference code and are ignored.
///
/// # Errors
///
/// Rejects an empty query, a batch outside two through four candidates, or a
/// slab whose length is not exactly four times the candidate pattern length.
#[allow(clippy::trivially_copy_pass_by_ref)]
pub(crate) fn narrow_banded_placement_distances_interleaved_batch_d3(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    interleaved_patterns: &[u8],
    output: &mut [NarrowPlacementDistances],
) -> Result<NarrowBandedFlavor, NarrowBandedError> {
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
    if narrow_avx2_available() {
        #[cfg(target_arch = "x86_64")]
        // SAFETY: runtime detection proves AVX2/SSSE3, and the validated slab
        // provides one complete four-byte code group for every DP access.
        unsafe {
            narrow_placement_distances_interleaved_batch_d3_avx2(
                reference_masks_by_query,
                query,
                interleaved_patterns,
                output,
            );
        }
        return Ok(NarrowBandedFlavor::Avx2);
    }
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
    Ok(NarrowBandedFlavor::Scalar)
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
) -> Result<NarrowBandedFlavor, NarrowBandedError> {
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
    if narrow_avx2_available() {
        #[cfg(target_arch = "x86_64")]
        // SAFETY: runtime detection proves AVX2 and dimensions above keep all
        // active candidate/start lanes inside one vector.
        unsafe {
            if let [destination] = output {
                *destination =
                    narrow_placement_distances_d5_avx2(reference_masks_by_query, query, patterns);
            } else if output.len() == 2 {
                let pair: &mut [NarrowPlacementDistances; 2] = output
                    .try_into()
                    .expect("validated distance-five SIMD batch has two candidates");
                narrow_placement_distances_batch_d5_avx2(
                    reference_masks_by_query,
                    query,
                    patterns,
                    pattern_len,
                    pair,
                );
            } else {
                unreachable!("validated non-quad distance-five batch has at most two candidates");
            }
        }
        return Ok(NarrowBandedFlavor::Avx2);
    }
    if narrow_sse42_available() {
        #[cfg(target_arch = "x86_64")]
        for (pattern, destination) in patterns.chunks_exact(pattern_len).zip(output) {
            // SAFETY: runtime detection proves SSE4.2 and the batch validation
            // proves every fixed distance-five candidate dimension.
            *destination = unsafe {
                narrow_placement_distances_sse42(
                    reference_masks_by_query,
                    query,
                    pattern,
                    MAX_DISTANCE,
                )
            };
        }
        return Ok(NarrowBandedFlavor::Sse42);
    }
    for (pattern, destination) in patterns.chunks_exact(pattern_len).zip(output) {
        *destination = narrow_placement_distances_scalar(
            reference_masks_by_query,
            query,
            pattern,
            MAX_DISTANCE,
        );
    }
    Ok(NarrowBandedFlavor::Scalar)
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
) -> Result<NarrowBandedFlavor, NarrowBandedError> {
    validate_narrow_batch(query, patterns, max_distance, output)?;
    let flavor = if narrow_avx2_available() {
        NarrowBandedFlavor::Avx2
    } else if narrow_sse42_available() {
        NarrowBandedFlavor::Sse42
    } else {
        NarrowBandedFlavor::Scalar
    };
    narrow_banded_prefix_batch_with_flavor(
        flavor,
        reference_masks_by_query,
        query,
        patterns,
        max_distance,
        output,
    )?;
    Ok(flavor)
}

/// Computes every in-budget endpoint for fixed-start candidate patterns.
///
/// Each candidate pattern contains `query.len() + 2 * max_distance` symbols.
/// Its fixed alignment start is `max_distance` symbols from the beginning, so
/// endpoint delta `0..=2 * max_distance` represents interval lengths
/// `query.len() - max_distance ..= query.len() + max_distance`. Unlike
/// [`narrow_banded_prefix_batch`], this routine never minimizes over alternate
/// starts. SSE4.2 evaluates up to 16 independent candidates in byte lanes;
/// AVX2 evaluates up to 32.
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
) -> Result<NarrowBandedFlavor, NarrowBandedError> {
    let pattern_length =
        validate_narrow_fixed_start_batch(query, patterns, max_distance, output.len())?;
    let flavor = if narrow_avx2_available() {
        NarrowBandedFlavor::Avx2
    } else if narrow_sse42_available() {
        NarrowBandedFlavor::Sse42
    } else {
        NarrowBandedFlavor::Scalar
    };
    match flavor {
        NarrowBandedFlavor::Scalar => {
            for (pattern, result) in patterns.chunks_exact(pattern_length).zip(output) {
                *result = narrow_fixed_start_scalar_one(
                    reference_masks_by_query,
                    query,
                    pattern,
                    max_distance,
                );
            }
        }
        NarrowBandedFlavor::Sse42 => {
            #[cfg(target_arch = "x86_64")]
            // SAFETY: runtime detection proves SSE4.2, and validation above
            // proves every fixed-width pattern and output dimension.
            unsafe {
                narrow_fixed_start_sse42(
                    reference_masks_by_query,
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
            // SAFETY: runtime detection proves AVX2, and validation above
            // proves every fixed-width pattern and output dimension.
            unsafe {
                narrow_fixed_start_avx2(
                    reference_masks_by_query,
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
    }
    Ok(flavor)
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
) -> Result<NarrowBandedFlavor, NarrowBandedError> {
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
    let flavor = if narrow_avx2_available() {
        NarrowBandedFlavor::Avx2
    } else {
        NarrowBandedFlavor::Scalar
    };
    match flavor {
        NarrowBandedFlavor::Scalar => {
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
        NarrowBandedFlavor::Sse42 => {
            unreachable!("gather dispatch does not select SSE4.2")
        }
        NarrowBandedFlavor::Avx2 => {
            #[cfg(target_arch = "x86_64")]
            // SAFETY: runtime detection proves AVX2 and validation bounds all
            // reference gathers.
            unsafe {
                narrow_fixed_start_gather_avx2(
                    reference_masks_by_query,
                    query,
                    reference,
                    starts,
                    max_distance,
                    output,
                );
            }
            #[cfg(not(target_arch = "x86_64"))]
            unreachable!("non-x86_64 cannot report AVX2 available");
        }
    }
    Ok(flavor)
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
) -> Result<NarrowBandedFlavor, NarrowBandedError> {
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
    let flavor = if narrow_avx2_available() && query_length > u64::BITS as usize {
        NarrowBandedFlavor::Avx2
    } else {
        NarrowBandedFlavor::Scalar
    };
    match flavor {
        NarrowBandedFlavor::Scalar => {
            for (pattern, result) in patterns.chunks_exact(pattern_length).zip(output) {
                *result = myers_prefix_distances_u128_scalar_one(
                    equality_masks,
                    query_length,
                    pattern,
                    max_distance,
                );
            }
        }
        NarrowBandedFlavor::Sse42 => {
            unreachable!("u128 dispatch does not select SSE4.2")
        }
        NarrowBandedFlavor::Avx2 => {
            #[cfg(target_arch = "x86_64")]
            // SAFETY: runtime detection proves AVX2, and all flat dimensions
            // and the two-word query bound were validated above.
            unsafe {
                myers_prefix_distances_u128_avx2(
                    equality_masks,
                    query_length,
                    patterns,
                    pattern_length,
                    max_distance,
                    output,
                );
            }
            #[cfg(not(target_arch = "x86_64"))]
            unreachable!("non-x86_64 cannot report AVX2 available");
        }
    }
    Ok(flavor)
}

/// Runs an explicitly selected narrow-band implementation for differential
/// qualification.
// This table is shared by scalar and SIMD dispatch; retaining one borrowed
// signature avoids copying it at every qualification seam.
#[allow(clippy::trivially_copy_pass_by_ref)]
pub(crate) fn narrow_banded_prefix_batch_with_flavor(
    flavor: NarrowBandedFlavor,
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    patterns: &[u8],
    max_distance: usize,
    output: &mut [NarrowBandedResult],
) -> Result<(), NarrowBandedError> {
    let pattern_length = validate_narrow_batch(query, patterns, max_distance, output)?;
    if matches!(flavor, NarrowBandedFlavor::Sse42) && !narrow_sse42_available()
        || matches!(flavor, NarrowBandedFlavor::Avx2) && !narrow_avx2_available()
    {
        return Err(NarrowBandedError::UnsupportedFlavor { flavor });
    }
    match flavor {
        NarrowBandedFlavor::Scalar => {
            for (pattern, result) in patterns.chunks_exact(pattern_length).zip(output) {
                *result = narrow_banded_scalar_one(
                    reference_masks_by_query,
                    query,
                    pattern,
                    max_distance,
                );
            }
        }
        NarrowBandedFlavor::Sse42 => {
            #[cfg(target_arch = "x86_64")]
            // SAFETY: runtime detection above proves SSE4.2, and validation
            // proves all fixed-width candidate and output dimensions.
            unsafe {
                narrow_banded_sse42(
                    reference_masks_by_query,
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
            // SAFETY: runtime detection above proves AVX2, and validation
            // proves all fixed-width candidate and output dimensions.
            unsafe {
                narrow_banded_avx2(
                    reference_masks_by_query,
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
    }
    Ok(())
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
