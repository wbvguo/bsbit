//! Narrow-band scalar and architecture-dispatched alignment kernels.

use super::{MAX_NARROW_BAND_DISTANCE, NarrowReferenceCode};
use core::fmt;

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

// This table is shared by scalar and SIMD dispatch; retain one borrowed
// signature for the paired differential-test surface.
#[allow(clippy::trivially_copy_pass_by_ref)]
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
#[allow(clippy::trivially_copy_pass_by_ref)]
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

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "sse4.2")]
// Unaligned SSE loads/stores intentionally accept byte-array pointers; the
// intrinsic provides the alignment guarantee that a typed dereference would.
#[allow(clippy::cast_ptr_alignment, clippy::needless_range_loop)]
pub(super) unsafe fn narrow_placement_distances_sse42(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    pattern: &[u8],
    max_distance: usize,
) -> NarrowPlacementDistances {
    use std::arch::x86_64::{
        __m128i, _mm_adds_epu8, _mm_loadu_si128, _mm_min_epu8, _mm_set1_epi8, _mm_storeu_si128,
    };

    let band_length = max_distance * 2 + 1;
    let capped = u8::try_from(max_distance)
        .unwrap_or(u8::MAX - 1)
        .saturating_add(1);
    let cap_vector = _mm_set1_epi8(capped.cast_signed());
    let one = _mm_set1_epi8(1);
    let mut distances = NarrowPlacementDistances::for_band(band_length, max_distance);

    for start_base in (0..band_length).step_by(16) {
        let active_starts = (band_length - start_base).min(16);
        let mut previous = [cap_vector; MAX_NARROW_ENDPOINTS];
        let mut current = [cap_vector; MAX_NARROW_ENDPOINTS];
        for (diagonal, slot) in previous[..band_length].iter_mut().enumerate() {
            let mut lanes = [capped; 16];
            for (lane, value) in lanes[..active_starts].iter_mut().enumerate() {
                let start = start_base + lane;
                if start <= diagonal {
                    *value = u8::try_from(diagonal - start).unwrap_or(capped).min(capped);
                }
            }
            // SAFETY: one vector reads exactly the 16-byte local array.
            *slot = unsafe { _mm_loadu_si128(lanes.as_ptr().cast::<__m128i>()) };
        }

        for (query_position, &query_code) in query.iter().enumerate() {
            current[..band_length].fill(cap_vector);
            let reference_mask = reference_masks_by_query
                .get(usize::from(query_code))
                .copied()
                .unwrap_or(0);
            for diagonal in 0..band_length {
                let reference_position = query_position + diagonal;
                let substitution = u8::from(
                    reference_mask & (1_u8 << usize::from(pattern[reference_position])) == 0,
                );
                let substitution = _mm_set1_epi8(substitution.cast_signed());
                let diagonal_score = _mm_adds_epu8(previous[diagonal], substitution);
                let query_gap = if diagonal + 1 < band_length {
                    _mm_adds_epu8(previous[diagonal + 1], one)
                } else {
                    cap_vector
                };
                let reference_gap = if diagonal != 0 {
                    _mm_adds_epu8(current[diagonal - 1], one)
                } else {
                    cap_vector
                };
                current[diagonal] = _mm_min_epu8(
                    _mm_min_epu8(_mm_min_epu8(diagonal_score, query_gap), reference_gap),
                    cap_vector,
                );
            }
            core::mem::swap(&mut previous, &mut current);
        }

        for (endpoint, &endpoint_distances) in previous[..band_length].iter().enumerate() {
            let mut lanes = [capped; 16];
            // SAFETY: one vector writes exactly the 16-byte local array.
            unsafe {
                _mm_storeu_si128(lanes.as_mut_ptr().cast::<__m128i>(), endpoint_distances);
            }
            for (lane, &distance) in lanes[..active_starts].iter().enumerate() {
                distances.insert_distance(start_base + lane, endpoint, distance);
            }
        }
    }

    distances
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
// Unaligned AVX2 loads/stores intentionally accept byte-array pointers; the
// intrinsic provides the alignment guarantee that a typed dereference would.
#[allow(clippy::cast_ptr_alignment, clippy::needless_range_loop)]
pub(super) unsafe fn narrow_placement_distances_avx2(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    pattern: &[u8],
    max_distance: usize,
) -> NarrowPlacementDistances {
    use std::arch::x86_64::{
        __m256i, _mm256_adds_epu8, _mm256_loadu_si256, _mm256_min_epu8, _mm256_set1_epi8,
        _mm256_storeu_si256,
    };

    let band_length = max_distance * 2 + 1;
    let capped = u8::try_from(max_distance)
        .unwrap_or(u8::MAX - 1)
        .saturating_add(1);
    let cap_vector = _mm256_set1_epi8(capped.cast_signed());
    let mut previous = [cap_vector; MAX_NARROW_ENDPOINTS];
    let mut current = [cap_vector; MAX_NARROW_ENDPOINTS];
    for (diagonal, slot) in previous[..band_length].iter_mut().enumerate() {
        let mut lanes = [capped; 32];
        for (start, lane) in lanes[..band_length].iter_mut().enumerate() {
            if start <= diagonal {
                *lane = u8::try_from(diagonal - start).unwrap_or(capped).min(capped);
            }
        }
        // SAFETY: one vector reads exactly the 32-byte local array.
        *slot = unsafe { _mm256_loadu_si256(lanes.as_ptr().cast::<__m256i>()) };
    }
    let one = _mm256_set1_epi8(1);
    for (query_position, &query_code) in query.iter().enumerate() {
        current[..band_length].fill(cap_vector);
        let reference_mask = reference_masks_by_query
            .get(usize::from(query_code))
            .copied()
            .unwrap_or(0);
        for diagonal in 0..band_length {
            let reference_position = query_position + diagonal;
            let substitution =
                u8::from(reference_mask & (1_u8 << usize::from(pattern[reference_position])) == 0);
            let substitution = _mm256_set1_epi8(substitution.cast_signed());
            let diagonal_score = _mm256_adds_epu8(previous[diagonal], substitution);
            let query_gap = if diagonal + 1 < band_length {
                _mm256_adds_epu8(previous[diagonal + 1], one)
            } else {
                cap_vector
            };
            let reference_gap = if diagonal != 0 {
                _mm256_adds_epu8(current[diagonal - 1], one)
            } else {
                cap_vector
            };
            current[diagonal] = _mm256_min_epu8(
                _mm256_min_epu8(_mm256_min_epu8(diagonal_score, query_gap), reference_gap),
                cap_vector,
            );
        }
        core::mem::swap(&mut previous, &mut current);
    }

    let mut distances = NarrowPlacementDistances::for_band(band_length, max_distance);
    for (endpoint, &endpoint_distances) in previous[..band_length].iter().enumerate() {
        let mut lanes = [capped; 32];
        // SAFETY: one vector writes exactly the 32-byte local array.
        unsafe {
            _mm256_storeu_si256(lanes.as_mut_ptr().cast::<__m256i>(), endpoint_distances);
        }
        for start in 0..band_length {
            distances.insert_distance(start, endpoint, lanes[start]);
        }
    }
    distances
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
// Fixed constants bound every narrowing cast, and the unaligned intrinsics
// intentionally accept byte-array pointers.
#[allow(clippy::cast_ptr_alignment, clippy::cast_possible_truncation)]
unsafe fn narrow_placement_distances_d3_avx2(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    pattern: &[u8],
) -> NarrowPlacementDistances {
    use std::arch::x86_64::{
        __m256i, _mm256_adds_epu8, _mm256_loadu_si256, _mm256_min_epu8, _mm256_set1_epi8,
        _mm256_storeu_si256,
    };

    const BAND_LENGTH: usize = 7;
    const MAX_DISTANCE: usize = 3;
    const CAPPED: u8 = 4;

    let cap_vector = _mm256_set1_epi8(CAPPED.cast_signed());
    let mut previous_storage = [cap_vector; BAND_LENGTH];
    let mut current_storage = [cap_vector; BAND_LENGTH];
    for (diagonal, slot) in previous_storage.iter_mut().enumerate() {
        let mut lanes = [CAPPED; 32];
        for (start, lane) in lanes[..BAND_LENGTH].iter_mut().enumerate() {
            if start <= diagonal {
                *lane = (diagonal - start) as u8;
            }
        }
        // SAFETY: one vector reads exactly the 32-byte local array.
        *slot = unsafe { _mm256_loadu_si256(lanes.as_ptr().cast::<__m256i>()) };
    }
    let mut previous = &mut previous_storage;
    let mut current = &mut current_storage;
    let one = _mm256_set1_epi8(1);
    for (query_position, &query_code) in query.iter().enumerate() {
        current.fill(cap_vector);
        let reference_mask = reference_masks_by_query
            .get(usize::from(query_code))
            .copied()
            .unwrap_or(0);
        for diagonal in 0..BAND_LENGTH {
            let reference_position = query_position + diagonal;
            let substitution =
                u8::from(reference_mask & (1_u8 << usize::from(pattern[reference_position])) == 0);
            let substitution = _mm256_set1_epi8(substitution.cast_signed());
            let diagonal_score = _mm256_adds_epu8(previous[diagonal], substitution);
            let query_gap = if diagonal + 1 < BAND_LENGTH {
                _mm256_adds_epu8(previous[diagonal + 1], one)
            } else {
                cap_vector
            };
            let reference_gap = if diagonal != 0 {
                _mm256_adds_epu8(current[diagonal - 1], one)
            } else {
                cap_vector
            };
            current[diagonal] = _mm256_min_epu8(
                _mm256_min_epu8(_mm256_min_epu8(diagonal_score, query_gap), reference_gap),
                cap_vector,
            );
        }
        core::mem::swap(&mut previous, &mut current);
    }

    let mut distances = NarrowPlacementDistances::for_band(BAND_LENGTH, MAX_DISTANCE);
    for (endpoint, &endpoint_distances) in previous.iter().enumerate() {
        let mut lanes = [CAPPED; 32];
        // SAFETY: one vector writes exactly the 32-byte local array.
        unsafe {
            _mm256_storeu_si256(lanes.as_mut_ptr().cast::<__m256i>(), endpoint_distances);
        }
        for (start, &distance) in lanes[..BAND_LENGTH].iter().enumerate() {
            distances.insert_distance(start, endpoint, distance);
        }
    }
    distances
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
// Fixed constants bound every narrowing cast, and the unaligned intrinsics
// intentionally accept byte-array pointers.
#[allow(clippy::cast_ptr_alignment, clippy::cast_possible_truncation)]
unsafe fn narrow_placement_distances_d5_avx2(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    pattern: &[u8],
) -> NarrowPlacementDistances {
    use std::arch::x86_64::{
        __m256i, _mm256_adds_epu8, _mm256_loadu_si256, _mm256_min_epu8, _mm256_set1_epi8,
        _mm256_storeu_si256,
    };

    const BAND_LENGTH: usize = 11;
    const MAX_DISTANCE: usize = 5;
    const CAPPED: u8 = 6;

    let cap_vector = _mm256_set1_epi8(CAPPED.cast_signed());
    let mut previous_storage = [cap_vector; BAND_LENGTH];
    let mut current_storage = [cap_vector; BAND_LENGTH];
    for (diagonal, slot) in previous_storage.iter_mut().enumerate() {
        let mut lanes = [CAPPED; 32];
        for (start, lane) in lanes[..BAND_LENGTH].iter_mut().enumerate() {
            if start <= diagonal {
                *lane = (diagonal - start) as u8;
            }
        }
        // SAFETY: one vector reads exactly the 32-byte local array.
        *slot = unsafe { _mm256_loadu_si256(lanes.as_ptr().cast::<__m256i>()) };
    }
    let mut previous = &mut previous_storage;
    let mut current = &mut current_storage;
    let one = _mm256_set1_epi8(1);
    for (query_position, &query_code) in query.iter().enumerate() {
        current.fill(cap_vector);
        let reference_mask = reference_masks_by_query
            .get(usize::from(query_code))
            .copied()
            .unwrap_or(0);
        for diagonal in 0..BAND_LENGTH {
            let reference_position = query_position + diagonal;
            let substitution =
                u8::from(reference_mask & (1_u8 << usize::from(pattern[reference_position])) == 0);
            let substitution = _mm256_set1_epi8(substitution.cast_signed());
            let diagonal_score = _mm256_adds_epu8(previous[diagonal], substitution);
            let query_gap = if diagonal + 1 < BAND_LENGTH {
                _mm256_adds_epu8(previous[diagonal + 1], one)
            } else {
                cap_vector
            };
            let reference_gap = if diagonal != 0 {
                _mm256_adds_epu8(current[diagonal - 1], one)
            } else {
                cap_vector
            };
            current[diagonal] = _mm256_min_epu8(
                _mm256_min_epu8(_mm256_min_epu8(diagonal_score, query_gap), reference_gap),
                cap_vector,
            );
        }
        core::mem::swap(&mut previous, &mut current);
    }

    let mut distances = NarrowPlacementDistances::for_band(BAND_LENGTH, MAX_DISTANCE);
    for (endpoint, &endpoint_distances) in previous.iter().enumerate() {
        let mut lanes = [CAPPED; 32];
        // SAFETY: one vector writes exactly the 32-byte local array.
        unsafe {
            _mm256_storeu_si256(lanes.as_mut_ptr().cast::<__m256i>(), endpoint_distances);
        }
        for (start, &distance) in lanes[..BAND_LENGTH].iter().enumerate() {
            distances.insert_distance(start, endpoint, distance);
        }
    }
    distances
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
// Unaligned AVX2 loads/stores intentionally accept byte-array pointers.
#[allow(clippy::cast_ptr_alignment)]
unsafe fn narrow_placement_distances_batch_avx2(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    patterns: &[u8],
    pattern_len: usize,
    max_distance: usize,
    output: &mut [NarrowPlacementDistances],
) {
    use std::arch::x86_64::{
        __m256i, _mm256_adds_epu8, _mm256_loadu_si256, _mm256_min_epu8, _mm256_or_si256,
        _mm256_set1_epi8, _mm256_setzero_si256, _mm256_storeu_si256,
    };

    let band_length = max_distance * 2 + 1;
    let capped = u8::try_from(max_distance)
        .unwrap_or(u8::MAX - 1)
        .saturating_add(1);
    let cap_vector = _mm256_set1_epi8(capped.cast_signed());
    let mut previous = [cap_vector; MAX_NARROW_ENDPOINTS];
    let mut current = [cap_vector; MAX_NARROW_ENDPOINTS];
    for (diagonal, slot) in previous[..band_length].iter_mut().enumerate() {
        let mut lanes = [capped; 32];
        for candidate in 0..output.len() {
            for start in 0..band_length {
                if start <= diagonal {
                    lanes[candidate * band_length + start] =
                        u8::try_from(diagonal - start).unwrap_or(capped).min(capped);
                }
            }
        }
        // SAFETY: one vector reads the complete 32-byte local array.
        *slot = unsafe { _mm256_loadu_si256(lanes.as_ptr().cast::<__m256i>()) };
    }
    let one = _mm256_set1_epi8(1);
    let lane_masks: [__m256i; 4] = core::array::from_fn(|candidate| {
        let mut lanes = [0_u8; 32];
        if candidate < output.len() {
            lanes[candidate * band_length..(candidate + 1) * band_length].fill(1);
        }
        // SAFETY: one vector reads the complete 32-byte local array.
        unsafe { _mm256_loadu_si256(lanes.as_ptr().cast::<__m256i>()) }
    });
    for (query_position, &query_code) in query.iter().enumerate() {
        current[..band_length].fill(cap_vector);
        let reference_mask = reference_masks_by_query
            .get(usize::from(query_code))
            .copied()
            .unwrap_or(0);
        for diagonal in 0..band_length {
            let mut substitution = _mm256_setzero_si256();
            for candidate in 0..output.len() {
                let reference_code = patterns[candidate * pattern_len + query_position + diagonal];
                if reference_mask & (1_u8 << usize::from(reference_code)) == 0 {
                    substitution = _mm256_or_si256(substitution, lane_masks[candidate]);
                }
            }
            let diagonal_score = _mm256_adds_epu8(previous[diagonal], substitution);
            let query_gap = if diagonal + 1 < band_length {
                _mm256_adds_epu8(previous[diagonal + 1], one)
            } else {
                cap_vector
            };
            let reference_gap = if diagonal != 0 {
                _mm256_adds_epu8(current[diagonal - 1], one)
            } else {
                cap_vector
            };
            current[diagonal] = _mm256_min_epu8(
                _mm256_min_epu8(_mm256_min_epu8(diagonal_score, query_gap), reference_gap),
                cap_vector,
            );
        }
        core::mem::swap(&mut previous, &mut current);
    }
    for destination in output.iter_mut() {
        destination.band_length = u8::try_from(band_length).expect("validated narrow band fits u8");
        destination.max_distance = u8::try_from(max_distance).expect("validated distance fits u8");
    }
    for (endpoint, &endpoint_distances) in previous[..band_length].iter().enumerate() {
        let mut lanes = [capped; 32];
        // SAFETY: one vector writes the complete 32-byte local array.
        unsafe {
            _mm256_storeu_si256(lanes.as_mut_ptr().cast::<__m256i>(), endpoint_distances);
        }
        for (candidate, destination) in output.iter_mut().enumerate() {
            for start in 0..band_length {
                destination.insert_distance(
                    start,
                    endpoint,
                    lanes[candidate * band_length + start],
                );
            }
        }
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
// Fixed distance-three constants bound the byte casts; unaligned AVX2
// intrinsics intentionally accept byte-array pointers.
#[allow(clippy::cast_ptr_alignment, clippy::cast_possible_truncation)]
unsafe fn narrow_placement_distances_batch_d3_avx2(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    patterns: &[u8],
    pattern_len: usize,
    output: &mut [NarrowPlacementDistances],
) {
    use std::arch::x86_64::{
        __m256i, _mm256_adds_epu8, _mm256_loadu_si256, _mm256_min_epu8, _mm256_or_si256,
        _mm256_set1_epi8, _mm256_setzero_si256, _mm256_storeu_si256,
    };

    const BAND_LENGTH: usize = 7;
    const MAX_DISTANCE: usize = 3;
    const CAPPED: u8 = 4;

    let cap_vector = _mm256_set1_epi8(CAPPED.cast_signed());
    let mut previous_storage = [cap_vector; BAND_LENGTH];
    let mut current_storage = [cap_vector; BAND_LENGTH];
    for (diagonal, slot) in previous_storage.iter_mut().enumerate() {
        let mut lanes = [CAPPED; 32];
        for candidate in 0..output.len() {
            for start in 0..BAND_LENGTH {
                if start <= diagonal {
                    lanes[candidate * BAND_LENGTH + start] = (diagonal - start) as u8;
                }
            }
        }
        // SAFETY: one vector reads the complete 32-byte local array.
        *slot = unsafe { _mm256_loadu_si256(lanes.as_ptr().cast::<__m256i>()) };
    }
    let mut previous = &mut previous_storage;
    let mut current = &mut current_storage;
    let one = _mm256_set1_epi8(1);
    let lane_masks: [__m256i; 4] = core::array::from_fn(|candidate| {
        let mut lanes = [0_u8; 32];
        lanes[candidate * BAND_LENGTH..(candidate + 1) * BAND_LENGTH].fill(1);
        // SAFETY: one vector reads the complete 32-byte local array.
        unsafe { _mm256_loadu_si256(lanes.as_ptr().cast::<__m256i>()) }
    });
    let substitution_masks: [__m256i; 16] = core::array::from_fn(|mismatch_mask| {
        let mut substitution = _mm256_setzero_si256();
        for (candidate, &lane_mask) in lane_masks.iter().enumerate() {
            if mismatch_mask & (1 << candidate) != 0 {
                substitution = _mm256_or_si256(substitution, lane_mask);
            }
        }
        substitution
    });
    for (query_position, &query_code) in query.iter().enumerate() {
        current.fill(cap_vector);
        let reference_mask = reference_masks_by_query
            .get(usize::from(query_code))
            .copied()
            .unwrap_or(0);
        for diagonal in 0..BAND_LENGTH {
            let mut mismatch_mask = 0_usize;
            for candidate in 0..output.len() {
                let reference_code = patterns[candidate * pattern_len + query_position + diagonal];
                if reference_mask & (1_u8 << usize::from(reference_code)) == 0 {
                    mismatch_mask |= 1 << candidate;
                }
            }
            let substitution = substitution_masks[mismatch_mask];
            let diagonal_score = _mm256_adds_epu8(previous[diagonal], substitution);
            let query_gap = if diagonal + 1 < BAND_LENGTH {
                _mm256_adds_epu8(previous[diagonal + 1], one)
            } else {
                cap_vector
            };
            let reference_gap = if diagonal != 0 {
                _mm256_adds_epu8(current[diagonal - 1], one)
            } else {
                cap_vector
            };
            current[diagonal] = _mm256_min_epu8(
                _mm256_min_epu8(_mm256_min_epu8(diagonal_score, query_gap), reference_gap),
                cap_vector,
            );
        }
        core::mem::swap(&mut previous, &mut current);
    }
    for destination in output.iter_mut() {
        destination.band_length = BAND_LENGTH as u8;
        destination.max_distance = MAX_DISTANCE as u8;
    }
    for (endpoint, &endpoint_distances) in previous.iter().enumerate() {
        let mut lanes = [CAPPED; 32];
        // SAFETY: one vector writes the complete 32-byte local array.
        unsafe {
            _mm256_storeu_si256(lanes.as_mut_ptr().cast::<__m256i>(), endpoint_distances);
        }
        for (candidate, destination) in output.iter_mut().enumerate() {
            for start in 0..BAND_LENGTH {
                destination.insert_distance(
                    start,
                    endpoint,
                    lanes[candidate * BAND_LENGTH + start],
                );
            }
        }
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
// Fixed distance-five constants bound the byte casts; unaligned AVX2
// intrinsics intentionally accept byte-array pointers.
#[allow(clippy::cast_ptr_alignment, clippy::cast_possible_truncation)]
unsafe fn narrow_placement_distances_batch_d5_avx2(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    patterns: &[u8],
    pattern_len: usize,
    output: &mut [NarrowPlacementDistances; 2],
) {
    use std::arch::x86_64::{
        __m256i, _mm256_adds_epu8, _mm256_loadu_si256, _mm256_min_epu8, _mm256_or_si256,
        _mm256_set1_epi8, _mm256_setzero_si256, _mm256_storeu_si256,
    };

    const BAND_LENGTH: usize = 11;
    const MAX_DISTANCE: usize = 5;
    const CAPPED: u8 = 6;

    let cap_vector = _mm256_set1_epi8(CAPPED.cast_signed());
    let mut previous_storage = [cap_vector; BAND_LENGTH];
    let mut current_storage = [cap_vector; BAND_LENGTH];
    for (diagonal, slot) in previous_storage.iter_mut().enumerate() {
        let mut lanes = [CAPPED; 32];
        for candidate in 0..2 {
            for start in 0..BAND_LENGTH {
                if start <= diagonal {
                    lanes[candidate * BAND_LENGTH + start] = (diagonal - start) as u8;
                }
            }
        }
        // SAFETY: one vector reads the complete 32-byte local array.
        *slot = unsafe { _mm256_loadu_si256(lanes.as_ptr().cast::<__m256i>()) };
    }
    let mut previous = &mut previous_storage;
    let mut current = &mut current_storage;
    let one = _mm256_set1_epi8(1);
    let lane_masks: [__m256i; 2] = core::array::from_fn(|candidate| {
        let mut lanes = [0_u8; 32];
        lanes[candidate * BAND_LENGTH..(candidate + 1) * BAND_LENGTH].fill(1);
        // SAFETY: one vector reads the complete 32-byte local array.
        unsafe { _mm256_loadu_si256(lanes.as_ptr().cast::<__m256i>()) }
    });
    let substitution_masks = [
        _mm256_setzero_si256(),
        lane_masks[0],
        lane_masks[1],
        _mm256_or_si256(lane_masks[0], lane_masks[1]),
    ];
    let (first_pattern, second_pattern) = patterns.split_at(pattern_len);
    for (query_position, &query_code) in query.iter().enumerate() {
        current.fill(cap_vector);
        let reference_mask = reference_masks_by_query
            .get(usize::from(query_code))
            .copied()
            .unwrap_or(0);
        for diagonal in 0..BAND_LENGTH {
            let first_code = first_pattern[query_position + diagonal];
            let second_code = second_pattern[query_position + diagonal];
            let mismatch_mask =
                usize::from(reference_mask & (1_u8 << usize::from(first_code)) == 0)
                    | (usize::from(reference_mask & (1_u8 << usize::from(second_code)) == 0) << 1);
            let substitution = substitution_masks[mismatch_mask];
            let diagonal_score = _mm256_adds_epu8(previous[diagonal], substitution);
            let query_gap = if diagonal + 1 < BAND_LENGTH {
                _mm256_adds_epu8(previous[diagonal + 1], one)
            } else {
                cap_vector
            };
            let reference_gap = if diagonal != 0 {
                _mm256_adds_epu8(current[diagonal - 1], one)
            } else {
                cap_vector
            };
            current[diagonal] = _mm256_min_epu8(
                _mm256_min_epu8(_mm256_min_epu8(diagonal_score, query_gap), reference_gap),
                cap_vector,
            );
        }
        core::mem::swap(&mut previous, &mut current);
    }
    for destination in output.iter_mut() {
        destination.band_length = BAND_LENGTH as u8;
        destination.max_distance = MAX_DISTANCE as u8;
    }
    for (endpoint, &endpoint_distances) in previous.iter().enumerate() {
        let mut lanes = [CAPPED; 32];
        // SAFETY: one vector writes the complete 32-byte local array.
        unsafe {
            _mm256_storeu_si256(lanes.as_mut_ptr().cast::<__m256i>(), endpoint_distances);
        }
        for (candidate, destination) in output.iter_mut().enumerate() {
            for start in 0..BAND_LENGTH {
                destination.insert_distance(
                    start,
                    endpoint,
                    lanes[candidate * BAND_LENGTH + start],
                );
            }
        }
    }
}

// Keep the relation table borrowed to match its SIMD counterpart exactly.
#[allow(clippy::trivially_copy_pass_by_ref)]
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
#[allow(clippy::trivially_copy_pass_by_ref)]
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
#[target_feature(enable = "avx2")]
// Keep the equality table borrowed across scalar/SIMD dispatch.
#[allow(clippy::trivially_copy_pass_by_ref)]
unsafe fn myers_prefix_distances_u128_avx2(
    equality_masks: &[u128; 5],
    query_length: usize,
    patterns: &[u8],
    pattern_length: usize,
    max_distance: usize,
    output: &mut [NarrowEndpointDistances],
) {
    for base in (0..output.len()).step_by(4) {
        let lanes = (output.len() - base).min(4);
        unsafe {
            myers_prefix_distances_u128_avx2_chunk(
                equality_masks,
                query_length,
                &patterns[base * pattern_length..(base + lanes) * pattern_length],
                pattern_length,
                lanes,
                max_distance,
                &mut output[base..base + lanes],
            );
        }
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
// This is one coupled two-word Myers recurrence. Splitting it would duplicate
// carry state; bounds validation makes its narrowing casts exact, and all
// pointer casts are consumed solely by unaligned AVX2 intrinsics.
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_ptr_alignment,
    clippy::too_many_lines
)]
unsafe fn myers_prefix_distances_u128_avx2_chunk(
    equality_masks: &[u128; 5],
    query_length: usize,
    patterns: &[u8],
    pattern_length: usize,
    lanes: usize,
    max_distance: usize,
    output: &mut [NarrowEndpointDistances],
) {
    use std::arch::x86_64::{
        __m256i, _mm256_add_epi64, _mm256_and_si256, _mm256_andnot_si256, _mm256_cmpeq_epi64,
        _mm256_cmpgt_epi64, _mm256_loadu_si256, _mm256_or_si256, _mm256_set1_epi64x,
        _mm256_setzero_si256, _mm256_slli_epi64, _mm256_srli_epi64, _mm256_storeu_si256,
        _mm256_sub_epi64, _mm256_xor_si256,
    };

    debug_assert!(query_length > u64::BITS as usize && query_length <= u128::BITS as usize);
    let all = _mm256_set1_epi64x(-1);
    let one = _mm256_set1_epi64x(1);
    let sign = _mm256_set1_epi64x(i64::MIN);
    let zero = _mm256_setzero_si256();
    let mut positive_low = all;
    let mut positive_high = all;
    let mut negative_low = zero;
    let mut negative_high = zero;
    let mut score = _mm256_set1_epi64x(i64::try_from(query_length).unwrap_or(i64::MAX));
    let high = _mm256_set1_epi64x((1_u64 << (query_length - u64::BITS as usize - 1)).cast_signed());
    let minimum_end = query_length - max_distance;
    output.fill(NarrowEndpointDistances::EMPTY);

    for position in 0..pattern_length {
        let mut equal_low = [0_u64; 4];
        let mut equal_high = [0_u64; 4];
        for lane in 0..lanes {
            let mask = match patterns[lane * pattern_length + position] {
                0 => equality_masks[0],
                1 => equality_masks[1],
                2 => equality_masks[2],
                3 => equality_masks[3],
                _ => 0,
            };
            equal_low[lane] = mask as u64;
            equal_high[lane] = (mask >> u64::BITS) as u64;
        }
        let equal_low = unsafe { _mm256_loadu_si256(equal_low.as_ptr().cast::<__m256i>()) };
        let equal_high = unsafe { _mm256_loadu_si256(equal_high.as_ptr().cast::<__m256i>()) };
        let horizontal_input_low = _mm256_or_si256(equal_low, negative_low);
        let horizontal_input_high = _mm256_or_si256(equal_high, negative_high);
        let addend_low = _mm256_and_si256(equal_low, positive_low);
        let addend_high = _mm256_and_si256(equal_high, positive_high);
        let sum_low = _mm256_add_epi64(addend_low, positive_low);
        let carry_mask = _mm256_cmpgt_epi64(
            _mm256_xor_si256(addend_low, sign),
            _mm256_xor_si256(sum_low, sign),
        );
        let sum_high = _mm256_add_epi64(
            _mm256_add_epi64(addend_high, positive_high),
            _mm256_and_si256(carry_mask, one),
        );
        let horizontal_low = _mm256_or_si256(_mm256_xor_si256(sum_low, positive_low), equal_low);
        let horizontal_high =
            _mm256_or_si256(_mm256_xor_si256(sum_high, positive_high), equal_high);
        let positive_horizontal_low = _mm256_or_si256(
            negative_low,
            _mm256_andnot_si256(_mm256_or_si256(horizontal_low, positive_low), all),
        );
        let positive_horizontal_high = _mm256_or_si256(
            negative_high,
            _mm256_andnot_si256(_mm256_or_si256(horizontal_high, positive_high), all),
        );
        let negative_horizontal_low = _mm256_and_si256(positive_low, horizontal_low);
        let negative_horizontal_high = _mm256_and_si256(positive_high, horizontal_high);
        let positive_hit =
            _mm256_cmpeq_epi64(_mm256_and_si256(positive_horizontal_high, high), high);
        let negative_hit =
            _mm256_cmpeq_epi64(_mm256_and_si256(negative_horizontal_high, high), high);
        score = _mm256_add_epi64(score, _mm256_and_si256(positive_hit, one));
        score = _mm256_sub_epi64(score, _mm256_and_si256(negative_hit, one));

        let shifted_positive_low =
            _mm256_or_si256(_mm256_slli_epi64(positive_horizontal_low, 1), one);
        let shifted_positive_high = _mm256_or_si256(
            _mm256_slli_epi64(positive_horizontal_high, 1),
            _mm256_srli_epi64(positive_horizontal_low, 63),
        );
        let shifted_negative_low = _mm256_slli_epi64(negative_horizontal_low, 1);
        let shifted_negative_high = _mm256_or_si256(
            _mm256_slli_epi64(negative_horizontal_high, 1),
            _mm256_srli_epi64(negative_horizontal_low, 63),
        );
        positive_low = _mm256_or_si256(
            shifted_negative_low,
            _mm256_andnot_si256(
                _mm256_or_si256(horizontal_input_low, shifted_positive_low),
                all,
            ),
        );
        positive_high = _mm256_or_si256(
            shifted_negative_high,
            _mm256_andnot_si256(
                _mm256_or_si256(horizontal_input_high, shifted_positive_high),
                all,
            ),
        );
        negative_low = _mm256_and_si256(shifted_positive_low, horizontal_input_low);
        negative_high = _mm256_and_si256(shifted_positive_high, horizontal_input_high);

        let prefix_length = position + 1;
        if prefix_length >= minimum_end {
            let mut scores = [u64::MAX; 4];
            unsafe {
                _mm256_storeu_si256(scores.as_mut_ptr().cast::<__m256i>(), score);
            }
            let endpoint = prefix_length - minimum_end;
            for lane in 0..lanes {
                if scores[lane] <= max_distance as u64 {
                    output[lane].insert(endpoint, scores[lane] as u32);
                }
            }
        }
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "sse4.2")]
// Keep the relation table borrowed across scalar/SIMD dispatch.
#[allow(clippy::trivially_copy_pass_by_ref)]
pub(super) unsafe fn narrow_fixed_start_sse42(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    patterns: &[u8],
    pattern_length: usize,
    max_distance: usize,
    output: &mut [NarrowEndpointDistances],
) {
    for base in (0..output.len()).step_by(16) {
        let lanes = (output.len() - base).min(16);
        unsafe {
            narrow_fixed_start_sse42_chunk(
                reference_masks_by_query,
                query,
                &patterns[base * pattern_length..(base + lanes) * pattern_length],
                pattern_length,
                lanes,
                max_distance,
                &mut output[base..base + lanes],
            );
        }
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "sse4.2")]
// Bounds validation makes the byte casts exact; the pointer casts are used
// only by unaligned SSE intrinsics.
#[allow(clippy::cast_possible_truncation, clippy::cast_ptr_alignment)]
unsafe fn narrow_fixed_start_sse42_chunk(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    patterns: &[u8],
    pattern_length: usize,
    lanes: usize,
    max_distance: usize,
    output: &mut [NarrowEndpointDistances],
) {
    use std::arch::x86_64::{
        __m128i, _mm_adds_epu8, _mm_loadu_si128, _mm_min_epu8, _mm_set1_epi8, _mm_storeu_si128,
    };

    let band_length = max_distance * 2 + 1;
    let capped = u8::try_from(max_distance)
        .unwrap_or(u8::MAX - 1)
        .saturating_add(1);
    let cap_vector = _mm_set1_epi8(capped.cast_signed());
    let one = _mm_set1_epi8(1);
    let mut previous = [cap_vector; MAX_NARROW_ENDPOINTS];
    let mut current = [cap_vector; MAX_NARROW_ENDPOINTS];
    for (diagonal, slot) in previous[..band_length].iter_mut().enumerate() {
        let initial = if diagonal < max_distance {
            capped
        } else {
            u8::try_from(diagonal - max_distance)
                .unwrap_or(capped)
                .min(capped)
        };
        *slot = _mm_set1_epi8(initial.cast_signed());
    }
    for (query_position, &query_code) in query.iter().enumerate() {
        current[..band_length].fill(cap_vector);
        let reference_mask = reference_masks_by_query
            .get(usize::from(query_code))
            .copied()
            .unwrap_or(0);
        for diagonal in 0..band_length {
            let reference_position = query_position + diagonal;
            let mut substitutions = [capped; 16];
            for lane in 0..lanes {
                let reference_code = patterns[lane * pattern_length + reference_position];
                substitutions[lane] = u8::from(
                    reference_code >= u8::BITS as u8
                        || reference_mask & (1_u8 << reference_code) == 0,
                );
            }
            let substitution = unsafe { _mm_loadu_si128(substitutions.as_ptr().cast::<__m128i>()) };
            let diagonal_score = _mm_adds_epu8(previous[diagonal], substitution);
            let query_gap = if diagonal + 1 < band_length {
                _mm_adds_epu8(previous[diagonal + 1], one)
            } else {
                cap_vector
            };
            let reference_gap = if diagonal != 0 {
                _mm_adds_epu8(current[diagonal - 1], one)
            } else {
                cap_vector
            };
            current[diagonal] = _mm_min_epu8(
                _mm_min_epu8(_mm_min_epu8(diagonal_score, query_gap), reference_gap),
                cap_vector,
            );
        }
        core::mem::swap(&mut previous, &mut current);
    }

    output.fill(NarrowEndpointDistances::EMPTY);
    for (endpoint, distances) in previous[..band_length].iter().copied().enumerate() {
        let mut lanes_out = [capped; 16];
        unsafe {
            _mm_storeu_si128(lanes_out.as_mut_ptr().cast::<__m128i>(), distances);
        }
        for lane in 0..lanes {
            let distance = lanes_out[lane];
            if usize::from(distance) <= max_distance {
                output[lane].insert(endpoint, u32::from(distance));
            }
        }
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
// Keep the relation table borrowed across scalar/SIMD dispatch.
#[allow(clippy::trivially_copy_pass_by_ref)]
unsafe fn narrow_fixed_start_avx2(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    patterns: &[u8],
    pattern_length: usize,
    max_distance: usize,
    output: &mut [NarrowEndpointDistances],
) {
    for base in (0..output.len()).step_by(32) {
        let lanes = (output.len() - base).min(32);
        unsafe {
            narrow_fixed_start_avx2_chunk(
                reference_masks_by_query,
                query,
                &patterns[base * pattern_length..(base + lanes) * pattern_length],
                pattern_length,
                lanes,
                max_distance,
                &mut output[base..base + lanes],
            );
        }
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
// Keep the relation table borrowed across scalar/SIMD dispatch.
#[allow(clippy::trivially_copy_pass_by_ref)]
unsafe fn narrow_fixed_start_gather_avx2<T: NarrowReferenceCode>(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    reference: &[T],
    starts: &[usize],
    max_distance: usize,
    output: &mut [NarrowEndpointDistances],
) {
    for base in (0..output.len()).step_by(32) {
        let lanes = (output.len() - base).min(32);
        unsafe {
            narrow_fixed_start_gather_avx2_chunk(
                reference_masks_by_query,
                query,
                reference,
                &starts[base..base + lanes],
                max_distance,
                &mut output[base..base + lanes],
            );
        }
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
// Bounds validation makes the byte casts exact; the pointer casts are used
// only by unaligned AVX2 intrinsics.
#[allow(clippy::cast_possible_truncation, clippy::cast_ptr_alignment)]
unsafe fn narrow_fixed_start_gather_avx2_chunk<T: NarrowReferenceCode>(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    reference: &[T],
    starts: &[usize],
    max_distance: usize,
    output: &mut [NarrowEndpointDistances],
) {
    use std::arch::x86_64::{
        __m256i, _mm256_adds_epu8, _mm256_loadu_si256, _mm256_min_epu8, _mm256_set1_epi8,
        _mm256_storeu_si256,
    };

    let lanes = starts.len();
    let band_length = max_distance * 2 + 1;
    let capped = u8::try_from(max_distance)
        .unwrap_or(u8::MAX - 1)
        .saturating_add(1);
    let cap_vector = _mm256_set1_epi8(capped.cast_signed());
    let one = _mm256_set1_epi8(1);
    let mut previous = [cap_vector; MAX_NARROW_ENDPOINTS];
    let mut current = [cap_vector; MAX_NARROW_ENDPOINTS];
    for (diagonal, slot) in previous[..band_length].iter_mut().enumerate() {
        let initial = if diagonal < max_distance {
            capped
        } else {
            u8::try_from(diagonal - max_distance)
                .unwrap_or(capped)
                .min(capped)
        };
        *slot = _mm256_set1_epi8(initial.cast_signed());
    }
    for (query_position, &query_code) in query.iter().enumerate() {
        current[..band_length].fill(cap_vector);
        let reference_mask = reference_masks_by_query
            .get(usize::from(query_code))
            .copied()
            .unwrap_or(0);
        for diagonal in 0..band_length {
            let shifted = query_position + diagonal;
            let mut substitutions = [capped; 32];
            if shifted >= max_distance {
                let offset = shifted - max_distance;
                for (lane, &start) in starts.iter().enumerate() {
                    let reference_code = reference[start + offset].narrow_reference_code();
                    substitutions[lane] = u8::from(
                        reference_code >= u8::BITS as u8
                            || reference_mask & (1_u8 << reference_code) == 0,
                    );
                }
            }
            let substitution =
                unsafe { _mm256_loadu_si256(substitutions.as_ptr().cast::<__m256i>()) };
            let diagonal_score = _mm256_adds_epu8(previous[diagonal], substitution);
            let query_gap = if diagonal + 1 < band_length {
                _mm256_adds_epu8(previous[diagonal + 1], one)
            } else {
                cap_vector
            };
            let reference_gap = if diagonal != 0 {
                _mm256_adds_epu8(current[diagonal - 1], one)
            } else {
                cap_vector
            };
            current[diagonal] = _mm256_min_epu8(
                _mm256_min_epu8(_mm256_min_epu8(diagonal_score, query_gap), reference_gap),
                cap_vector,
            );
        }
        core::mem::swap(&mut previous, &mut current);
    }
    output.fill(NarrowEndpointDistances::EMPTY);
    for (endpoint, distances) in previous[..band_length].iter().copied().enumerate() {
        let mut lanes_out = [capped; 32];
        unsafe {
            _mm256_storeu_si256(lanes_out.as_mut_ptr().cast::<__m256i>(), distances);
        }
        for lane in 0..lanes {
            let distance = lanes_out[lane];
            if usize::from(distance) <= max_distance {
                output[lane].insert(endpoint, u32::from(distance));
            }
        }
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
// Bounds validation makes the byte casts exact; the pointer casts are used
// only by unaligned AVX2 intrinsics.
#[allow(clippy::cast_possible_truncation, clippy::cast_ptr_alignment)]
unsafe fn narrow_fixed_start_avx2_chunk(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    patterns: &[u8],
    pattern_length: usize,
    lanes: usize,
    max_distance: usize,
    output: &mut [NarrowEndpointDistances],
) {
    use std::arch::x86_64::{
        __m256i, _mm256_adds_epu8, _mm256_loadu_si256, _mm256_min_epu8, _mm256_set1_epi8,
        _mm256_storeu_si256,
    };

    let band_length = max_distance * 2 + 1;
    let capped = u8::try_from(max_distance)
        .unwrap_or(u8::MAX - 1)
        .saturating_add(1);
    let cap_vector = _mm256_set1_epi8(capped.cast_signed());
    let one = _mm256_set1_epi8(1);
    let mut previous = [cap_vector; MAX_NARROW_ENDPOINTS];
    let mut current = [cap_vector; MAX_NARROW_ENDPOINTS];
    for (diagonal, slot) in previous[..band_length].iter_mut().enumerate() {
        let initial = if diagonal < max_distance {
            capped
        } else {
            u8::try_from(diagonal - max_distance)
                .unwrap_or(capped)
                .min(capped)
        };
        *slot = _mm256_set1_epi8(initial.cast_signed());
    }
    for (query_position, &query_code) in query.iter().enumerate() {
        current[..band_length].fill(cap_vector);
        let reference_mask = reference_masks_by_query
            .get(usize::from(query_code))
            .copied()
            .unwrap_or(0);
        for diagonal in 0..band_length {
            let reference_position = query_position + diagonal;
            let mut substitutions = [capped; 32];
            for lane in 0..lanes {
                let reference_code = patterns[lane * pattern_length + reference_position];
                substitutions[lane] = u8::from(
                    reference_code >= u8::BITS as u8
                        || reference_mask & (1_u8 << reference_code) == 0,
                );
            }
            let substitution =
                unsafe { _mm256_loadu_si256(substitutions.as_ptr().cast::<__m256i>()) };
            let diagonal_score = _mm256_adds_epu8(previous[diagonal], substitution);
            let query_gap = if diagonal + 1 < band_length {
                _mm256_adds_epu8(previous[diagonal + 1], one)
            } else {
                cap_vector
            };
            let reference_gap = if diagonal != 0 {
                _mm256_adds_epu8(current[diagonal - 1], one)
            } else {
                cap_vector
            };
            current[diagonal] = _mm256_min_epu8(
                _mm256_min_epu8(_mm256_min_epu8(diagonal_score, query_gap), reference_gap),
                cap_vector,
            );
        }
        core::mem::swap(&mut previous, &mut current);
    }

    output.fill(NarrowEndpointDistances::EMPTY);
    for (endpoint, distances) in previous[..band_length].iter().copied().enumerate() {
        let mut lanes_out = [capped; 32];
        unsafe {
            _mm256_storeu_si256(lanes_out.as_mut_ptr().cast::<__m256i>(), distances);
        }
        for lane in 0..lanes {
            let distance = lanes_out[lane];
            if usize::from(distance) <= max_distance {
                output[lane].insert(endpoint, u32::from(distance));
            }
        }
    }
}

#[cfg(target_arch = "x86_64")]
pub(super) fn narrow_sse42_available() -> bool {
    std::arch::is_x86_feature_detected!("sse4.2")
}

#[cfg(not(target_arch = "x86_64"))]
pub(super) const fn narrow_sse42_available() -> bool {
    false
}

#[cfg(target_arch = "x86_64")]
pub(super) fn narrow_avx2_available() -> bool {
    std::arch::is_x86_feature_detected!("avx2")
}

#[cfg(not(target_arch = "x86_64"))]
pub(super) const fn narrow_avx2_available() -> bool {
    false
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "sse4.2")]
// Keep the relation table borrowed across scalar/SIMD dispatch.
#[allow(clippy::trivially_copy_pass_by_ref)]
unsafe fn narrow_banded_sse42(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    patterns: &[u8],
    pattern_length: usize,
    max_distance: usize,
    output: &mut [NarrowBandedResult],
) {
    for base in (0..output.len()).step_by(4) {
        let lanes = (output.len() - base).min(4);
        unsafe {
            narrow_banded_sse42_chunk(
                reference_masks_by_query,
                query,
                &patterns[base * pattern_length..(base + lanes) * pattern_length],
                pattern_length,
                lanes,
                max_distance,
                &mut output[base..base + lanes],
            );
        }
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "sse4.2")]
// Unaligned SSE loads/stores intentionally accept the lane-array pointers.
#[allow(clippy::cast_ptr_alignment)]
unsafe fn narrow_banded_sse42_chunk(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    patterns: &[u8],
    pattern_length: usize,
    lanes: usize,
    max_distance: usize,
    output: &mut [NarrowBandedResult],
) {
    use std::arch::x86_64::{
        __m128i, _mm_add_epi32, _mm_and_si128, _mm_andnot_si128, _mm_loadu_si128, _mm_or_si128,
        _mm_set1_epi32, _mm_setzero_si128, _mm_srli_epi32, _mm_storeu_si128, _mm_sub_epi32,
        _mm_xor_si128,
    };

    let band_length = max_distance * 2 + 1;
    let high_bit = 1_u32 << (band_length - 1);
    let mut initial = [[0_u32; 4]; 5];
    for lane in 0..lanes {
        let pattern = &patterns[lane * pattern_length..(lane + 1) * pattern_length];
        for (position, &reference_code) in pattern[..band_length].iter().enumerate() {
            if let Some(bits) = initial.get_mut(usize::from(reference_code)) {
                bits[lane] |= 1_u32 << position;
            }
        }
    }
    let mut peq: [__m128i; 5] =
        core::array::from_fn(|code| unsafe { _mm_loadu_si128(initial[code].as_ptr().cast()) });
    let all = _mm_set1_epi32(-1);
    let one = _mm_set1_epi32(1);
    let mut positive = _mm_setzero_si128();
    let mut negative = _mm_setzero_si128();
    let mut error = _mm_setzero_si128();
    for (query_position, &query_code) in query.iter().enumerate() {
        let reference_mask = reference_masks_by_query
            .get(usize::from(query_code))
            .copied()
            .unwrap_or(0);
        let mut equal = _mm_setzero_si128();
        for (reference_code, &bits) in peq.iter().enumerate() {
            if reference_mask & (1_u8 << reference_code) != 0 {
                equal = _mm_or_si128(equal, bits);
            }
        }
        let horizontal_input = _mm_or_si128(equal, negative);
        let horizontal = _mm_or_si128(
            _mm_xor_si128(
                _mm_add_epi32(_mm_and_si128(horizontal_input, positive), positive),
                positive,
            ),
            horizontal_input,
        );
        let negative_horizontal = _mm_and_si128(positive, horizontal);
        let positive_horizontal = _mm_or_si128(
            negative,
            _mm_andnot_si128(_mm_or_si128(positive, horizontal), all),
        );
        let shifted = _mm_srli_epi32(horizontal, 1);
        negative = _mm_and_si128(shifted, positive_horizontal);
        positive = _mm_or_si128(
            negative_horizontal,
            _mm_andnot_si128(_mm_or_si128(shifted, positive_horizontal), all),
        );
        error = _mm_sub_epi32(_mm_add_epi32(error, one), _mm_and_si128(horizontal, one));
        if query_position + 1 != query.len() {
            for bits in &mut peq {
                *bits = _mm_srli_epi32(*bits, 1);
            }
            let entering_position = band_length + query_position;
            let mut entering = [[0_u32; 4]; 5];
            for lane in 0..lanes {
                let reference_code = patterns[lane * pattern_length + entering_position];
                if let Some(bits) = entering.get_mut(usize::from(reference_code)) {
                    bits[lane] = high_bit;
                }
            }
            for (bits, additions) in peq.iter_mut().zip(&entering) {
                *bits = _mm_or_si128(*bits, unsafe { _mm_loadu_si128(additions.as_ptr().cast()) });
            }
        }
    }
    let mut errors = [0_u32; 4];
    let mut positives = [0_u32; 4];
    let mut negatives = [0_u32; 4];
    unsafe {
        _mm_storeu_si128(errors.as_mut_ptr().cast(), error);
        _mm_storeu_si128(positives.as_mut_ptr().cast(), positive);
        _mm_storeu_si128(negatives.as_mut_ptr().cast(), negative);
    }
    for lane in 0..lanes {
        output[lane] = finish_narrow_result(
            errors[lane],
            positives[lane],
            negatives[lane],
            query.len(),
            max_distance,
        );
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
// Keep the relation table borrowed across scalar/SIMD dispatch.
#[allow(clippy::trivially_copy_pass_by_ref)]
unsafe fn narrow_banded_avx2(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    patterns: &[u8],
    pattern_length: usize,
    max_distance: usize,
    output: &mut [NarrowBandedResult],
) {
    for base in (0..output.len()).step_by(8) {
        let lanes = (output.len() - base).min(8);
        unsafe {
            narrow_banded_avx2_chunk(
                reference_masks_by_query,
                query,
                &patterns[base * pattern_length..(base + lanes) * pattern_length],
                pattern_length,
                lanes,
                max_distance,
                &mut output[base..base + lanes],
            );
        }
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
// Unaligned AVX2 loads/stores intentionally accept the lane-array pointers.
#[allow(clippy::cast_ptr_alignment)]
unsafe fn narrow_banded_avx2_chunk(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    patterns: &[u8],
    pattern_length: usize,
    lanes: usize,
    max_distance: usize,
    output: &mut [NarrowBandedResult],
) {
    use std::arch::x86_64::{
        __m256i, _mm256_add_epi32, _mm256_and_si256, _mm256_andnot_si256, _mm256_loadu_si256,
        _mm256_or_si256, _mm256_set1_epi32, _mm256_setzero_si256, _mm256_srli_epi32,
        _mm256_storeu_si256, _mm256_sub_epi32, _mm256_xor_si256,
    };

    let band_length = max_distance * 2 + 1;
    let high_bit = 1_u32 << (band_length - 1);
    let mut initial = [[0_u32; 8]; 5];
    for lane in 0..lanes {
        let pattern = &patterns[lane * pattern_length..(lane + 1) * pattern_length];
        for (position, &reference_code) in pattern[..band_length].iter().enumerate() {
            if let Some(bits) = initial.get_mut(usize::from(reference_code)) {
                bits[lane] |= 1_u32 << position;
            }
        }
    }
    let mut peq: [__m256i; 5] =
        core::array::from_fn(|code| unsafe { _mm256_loadu_si256(initial[code].as_ptr().cast()) });
    let all = _mm256_set1_epi32(-1);
    let one = _mm256_set1_epi32(1);
    let mut positive = _mm256_setzero_si256();
    let mut negative = _mm256_setzero_si256();
    let mut error = _mm256_setzero_si256();
    for (query_position, &query_code) in query.iter().enumerate() {
        let reference_mask = reference_masks_by_query
            .get(usize::from(query_code))
            .copied()
            .unwrap_or(0);
        let mut equal = _mm256_setzero_si256();
        for (reference_code, &bits) in peq.iter().enumerate() {
            if reference_mask & (1_u8 << reference_code) != 0 {
                equal = _mm256_or_si256(equal, bits);
            }
        }
        let horizontal_input = _mm256_or_si256(equal, negative);
        let horizontal = _mm256_or_si256(
            _mm256_xor_si256(
                _mm256_add_epi32(_mm256_and_si256(horizontal_input, positive), positive),
                positive,
            ),
            horizontal_input,
        );
        let negative_horizontal = _mm256_and_si256(positive, horizontal);
        let positive_horizontal = _mm256_or_si256(
            negative,
            _mm256_andnot_si256(_mm256_or_si256(positive, horizontal), all),
        );
        let shifted = _mm256_srli_epi32(horizontal, 1);
        negative = _mm256_and_si256(shifted, positive_horizontal);
        positive = _mm256_or_si256(
            negative_horizontal,
            _mm256_andnot_si256(_mm256_or_si256(shifted, positive_horizontal), all),
        );
        error = _mm256_sub_epi32(
            _mm256_add_epi32(error, one),
            _mm256_and_si256(horizontal, one),
        );
        if query_position + 1 != query.len() {
            for bits in &mut peq {
                *bits = _mm256_srli_epi32(*bits, 1);
            }
            let entering_position = band_length + query_position;
            let mut entering = [[0_u32; 8]; 5];
            for lane in 0..lanes {
                let reference_code = patterns[lane * pattern_length + entering_position];
                if let Some(bits) = entering.get_mut(usize::from(reference_code)) {
                    bits[lane] = high_bit;
                }
            }
            for (bits, additions) in peq.iter_mut().zip(&entering) {
                *bits = _mm256_or_si256(*bits, unsafe {
                    _mm256_loadu_si256(additions.as_ptr().cast())
                });
            }
        }
    }
    let mut errors = [0_u32; 8];
    let mut positives = [0_u32; 8];
    let mut negatives = [0_u32; 8];
    unsafe {
        _mm256_storeu_si256(errors.as_mut_ptr().cast(), error);
        _mm256_storeu_si256(positives.as_mut_ptr().cast(), positive);
        _mm256_storeu_si256(negatives.as_mut_ptr().cast(), negative);
    }
    for lane in 0..lanes {
        output[lane] = finish_narrow_result(
            errors[lane],
            positives[lane],
            negatives[lane],
            query.len(),
            max_distance,
        );
    }
}
