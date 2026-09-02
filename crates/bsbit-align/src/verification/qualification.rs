//! Test-only scalar/SIMD qualification surface.

use core::fmt;

use super::narrow::{
    myers_prefix_distances_u128_scalar_one, narrow_avx2_available, narrow_banded_scalar_one,
    narrow_fixed_start_scalar_one, narrow_placement_distances_scalar, narrow_sse42_available,
};
#[cfg(target_arch = "x86_64")]
use super::narrow::{
    narrow_fixed_start_sse42, narrow_placement_distances_avx2, narrow_placement_distances_sse42,
};
use super::*;

/// Runtime implementation used for one completed batch.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[cfg(test)]
pub(crate) enum KernelFlavor {
    /// Portable single-candidate machine-word implementation.
    Scalar,
    /// Two candidates evaluated in parallel with SSE4.2 64-bit lanes.
    Sse42,
    /// Four candidates evaluated in parallel with AVX2 64-bit lanes.
    Avx2,
    /// Eight candidates evaluated in parallel with AVX-512 64-bit lanes.
    Avx512,
}

#[cfg(test)]
impl fmt::Display for KernelFlavor {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::Scalar => "scalar-u64",
            Self::Sse42 => "sse4.2-u64x2",
            Self::Avx2 => "avx2-u64x4",
            Self::Avx512 => "avx512-u64x8",
        })
    }
}

/// Invalid safe-wrapper dimensions for a bit-vector batch.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[cfg(test)]
pub(crate) enum BatchError {
    /// Query length was outside the supported one-word domain.
    QueryLength {
        /// Supplied query bases.
        observed: usize,
    },
    /// Candidate start and length arrays differ in size.
    CandidateDimension {
        /// Start count.
        starts: usize,
        /// Length count.
        lengths: usize,
    },
    /// Output storage does not have one slot per candidate.
    OutputDimension {
        /// Candidate count.
        candidates: usize,
        /// Output slot count.
        outputs: usize,
    },
    /// A candidate range overflowed or exceeded the reference window.
    CandidateOutOfBounds {
        /// Candidate ordinal.
        candidate: usize,
        /// Inclusive candidate start.
        start: usize,
        /// Candidate length.
        length: usize,
        /// Available reference bases.
        reference_length: usize,
    },
    /// A forced test/benchmark implementation is unavailable on this CPU.
    UnsupportedFlavor {
        /// Requested runtime implementation.
        flavor: KernelFlavor,
    },
}

#[cfg(test)]
impl fmt::Display for BatchError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::QueryLength { observed } => write!(
                formatter,
                "Myers query length {observed} is outside 1..={MAX_QUERY_BASES}"
            ),
            Self::CandidateDimension { starts, lengths } => write!(
                formatter,
                "candidate start/length counts differ: {starts}/{lengths}"
            ),
            Self::OutputDimension {
                candidates,
                outputs,
            } => write!(
                formatter,
                "candidate/output counts differ: {candidates}/{outputs}"
            ),
            Self::CandidateOutOfBounds {
                candidate,
                start,
                length,
                reference_length,
            } => write!(
                formatter,
                "candidate {candidate} range {start}+{length} exceeds reference length {reference_length}"
            ),
            Self::UnsupportedFlavor { flavor } => {
                write!(
                    formatter,
                    "alignment kernel {flavor} is unavailable on this CPU"
                )
            }
        }
    }
}

#[cfg(test)]
impl std::error::Error for BatchError {}

/// Computes exact global edit distances for several reference intervals.
///
/// `reference_codes` uses `0=A`, `1=C`, `2=G`, `3=T`; every other byte is
/// treated as unknown/mismatch. Candidate intervals may have different
/// lengths. The qualified automatic default is scalar because the available
/// AVX2 implementation is slower on the qualification host; explicit SIMD
/// selection is retained for equivalence and benchmark work. Empty batches are
/// valid and report scalar.
///
/// # Errors
///
/// Returns a checked dimension/range error before writing any output.
#[cfg(test)]
pub(crate) fn myers_distance_batch(
    equality_masks: &[u64; 5],
    query_length: usize,
    reference_codes: &[u8],
    starts: &[usize],
    lengths: &[usize],
    output: &mut [u64],
) -> Result<KernelFlavor, BatchError> {
    validate_batch(
        query_length,
        reference_codes.len(),
        starts,
        lengths,
        output.len(),
    )?;
    let preferred = preferred_flavor();
    let flavor = if flavor_available(preferred) {
        preferred
    } else {
        KernelFlavor::Scalar
    };
    // SAFETY: validation proved every candidate range and output slot. Runtime
    // detection proves the selected target features before the specialized
    // function is entered.
    unsafe {
        dispatch_validated(
            flavor,
            equality_masks,
            query_length,
            reference_codes,
            starts,
            lengths,
            output,
        );
    }
    Ok(flavor)
}

/// Runs an explicitly selected implementation for equivalence/benchmark work.
///
/// Normal product code should use [`myers_distance_batch`].
///
/// # Errors
///
/// Returns a dimension/range failure or rejects a SIMD flavor unavailable on
/// the current CPU before writing output.
#[cfg(test)]
pub(crate) fn myers_distance_batch_with_flavor(
    flavor: KernelFlavor,
    equality_masks: &[u64; 5],
    query_length: usize,
    reference_codes: &[u8],
    starts: &[usize],
    lengths: &[usize],
    output: &mut [u64],
) -> Result<(), BatchError> {
    validate_batch(
        query_length,
        reference_codes.len(),
        starts,
        lengths,
        output.len(),
    )?;
    if !flavor_available(flavor) {
        return Err(BatchError::UnsupportedFlavor { flavor });
    }
    // SAFETY: the same validated range and runtime-feature invariants as the
    // normal dispatcher hold for this explicitly selected flavor.
    unsafe {
        dispatch_validated(
            flavor,
            equality_masks,
            query_length,
            reference_codes,
            starts,
            lengths,
            output,
        );
    }
    Ok(())
}

/// Returns the implementation accepted for automatic product dispatch.
///
/// SSE4.2/AVX2/AVX-512 remain explicitly selectable alternate implementations
/// until each has representative evidence on qualified hardware. CPU feature
/// presence alone is not a performance qualification.
#[must_use]
#[cfg(test)]
pub(crate) const fn preferred_flavor() -> KernelFlavor {
    KernelFlavor::Scalar
}

#[cfg(test)]
fn flavor_available(flavor: KernelFlavor) -> bool {
    match flavor {
        KernelFlavor::Scalar => true,
        KernelFlavor::Sse42 => sse42_available(),
        KernelFlavor::Avx2 => avx2_available(),
        KernelFlavor::Avx512 => avx512_available(),
    }
}

#[cfg(target_arch = "x86_64")]
#[cfg(test)]
fn sse42_available() -> bool {
    std::arch::is_x86_feature_detected!("sse4.2")
}

#[cfg(not(target_arch = "x86_64"))]
#[cfg(test)]
const fn sse42_available() -> bool {
    false
}

#[cfg(target_arch = "x86_64")]
#[cfg(test)]
fn avx2_available() -> bool {
    std::arch::is_x86_feature_detected!("avx2")
}

#[cfg(not(target_arch = "x86_64"))]
#[cfg(test)]
const fn avx2_available() -> bool {
    false
}

#[cfg(target_arch = "x86_64")]
#[cfg(test)]
fn avx512_available() -> bool {
    std::arch::is_x86_feature_detected!("avx512f")
        && std::arch::is_x86_feature_detected!("avx512bw")
}

#[cfg(not(target_arch = "x86_64"))]
#[cfg(test)]
const fn avx512_available() -> bool {
    false
}

#[cfg(test)]
fn validate_batch(
    query_length: usize,
    reference_length: usize,
    starts: &[usize],
    lengths: &[usize],
    output_length: usize,
) -> Result<(), BatchError> {
    if !(1..=MAX_QUERY_BASES).contains(&query_length) {
        return Err(BatchError::QueryLength {
            observed: query_length,
        });
    }
    if starts.len() != lengths.len() {
        return Err(BatchError::CandidateDimension {
            starts: starts.len(),
            lengths: lengths.len(),
        });
    }
    if output_length != starts.len() {
        return Err(BatchError::OutputDimension {
            candidates: starts.len(),
            outputs: output_length,
        });
    }
    for (candidate, (&start, &length)) in starts.iter().zip(lengths).enumerate() {
        if start
            .checked_add(length)
            .is_none_or(|end| end > reference_length)
        {
            return Err(BatchError::CandidateOutOfBounds {
                candidate,
                start,
                length,
                reference_length,
            });
        }
    }
    Ok(())
}

#[cfg(test)]
unsafe fn dispatch_validated(
    flavor: KernelFlavor,
    equality_masks: &[u64; 5],
    query_length: usize,
    reference_codes: &[u8],
    starts: &[usize],
    lengths: &[usize],
    output: &mut [u64],
) {
    match flavor {
        KernelFlavor::Scalar => scalar_batch(
            equality_masks,
            query_length,
            reference_codes,
            starts,
            lengths,
            output,
        ),
        #[cfg(target_arch = "x86_64")]
        KernelFlavor::Sse42 => unsafe {
            sse42_batch(
                equality_masks,
                query_length,
                reference_codes,
                starts,
                lengths,
                output,
            );
        },
        #[cfg(target_arch = "x86_64")]
        KernelFlavor::Avx2 => unsafe {
            avx2_batch(
                equality_masks,
                query_length,
                reference_codes,
                starts,
                lengths,
                output,
            );
        },
        #[cfg(target_arch = "x86_64")]
        KernelFlavor::Avx512 => unsafe {
            avx512_batch(
                equality_masks,
                query_length,
                reference_codes,
                starts,
                lengths,
                output,
            );
        },
        #[cfg(not(target_arch = "x86_64"))]
        KernelFlavor::Sse42 | KernelFlavor::Avx2 | KernelFlavor::Avx512 => {
            unreachable!("feature detection is scalar")
        }
    }
}

#[cfg(test)]
fn scalar_batch(
    equality_masks: &[u64; 5],
    query_length: usize,
    reference_codes: &[u8],
    starts: &[usize],
    lengths: &[usize],
    output: &mut [u64],
) {
    for ((&start, &length), result) in starts.iter().zip(lengths).zip(output) {
        *result = scalar_candidate(
            equality_masks,
            query_length,
            &reference_codes[start..start + length],
        );
    }
}

#[cfg(test)]
fn scalar_candidate(equality_masks: &[u64; 5], query_length: usize, reference: &[u8]) -> u64 {
    let mut positive = !0_u64;
    let mut negative = 0_u64;
    let mut score = u64::try_from(query_length).expect("query length fits u64");
    let high_bit = 1_u64 << (query_length - 1);
    for &code in reference {
        let equal = equality_mask(equality_masks, code);
        let horizontal_input = equal | negative;
        let horizontal = (((equal & positive).wrapping_add(positive)) ^ positive) | equal;
        let positive_horizontal = negative | !(horizontal | positive);
        let negative_horizontal = positive & horizontal;
        if positive_horizontal & high_bit != 0 {
            score += 1;
        } else if negative_horizontal & high_bit != 0 {
            score -= 1;
        }
        let shifted_positive = (positive_horizontal << 1) | 1;
        let shifted_negative = negative_horizontal << 1;
        positive = shifted_negative | !(horizontal_input | shifted_positive);
        negative = shifted_positive & horizontal_input;
    }
    score
}

#[cfg(test)]
const fn equality_mask(masks: &[u64; 5], code: u8) -> u64 {
    match code {
        0 => masks[0],
        1 => masks[1],
        2 => masks[2],
        3 => masks[3],
        _ => 0,
    }
}

#[cfg(target_arch = "x86_64")]
#[cfg(test)]
#[target_feature(enable = "sse4.2")]
unsafe fn sse42_batch(
    equality_masks: &[u64; 5],
    query_length: usize,
    reference_codes: &[u8],
    starts: &[usize],
    lengths: &[usize],
    output: &mut [u64],
) {
    for base in (0..starts.len()).step_by(2) {
        let lanes = (starts.len() - base).min(2);
        unsafe {
            sse42_chunk(
                equality_masks,
                query_length,
                reference_codes,
                &starts[base..base + lanes],
                &lengths[base..base + lanes],
                &mut output[base..base + lanes],
            );
        }
    }
}

#[cfg(target_arch = "x86_64")]
#[cfg(test)]
#[target_feature(enable = "sse4.2")]
unsafe fn sse42_chunk(
    equality_masks: &[u64; 5],
    query_length: usize,
    reference_codes: &[u8],
    starts: &[usize],
    lengths: &[usize],
    output: &mut [u64],
) {
    use std::arch::x86_64::{
        _mm_add_epi64, _mm_and_si128, _mm_andnot_si128, _mm_cmpeq_epi64, _mm_or_si128,
        _mm_set_epi64x, _mm_set1_epi64x, _mm_setzero_si128, _mm_slli_epi64, _mm_storeu_si128,
        _mm_sub_epi64, _mm_xor_si128,
    };

    if lengths.windows(2).all(|pair| pair[0] == pair[1]) {
        unsafe {
            sse42_equal_length_chunk(
                equality_masks,
                query_length,
                reference_codes,
                starts,
                lengths.first().copied().unwrap_or(0),
                output,
            );
        }
        return;
    }

    let all = _mm_set1_epi64x(-1);
    let ones = _mm_set1_epi64x(1);
    let high = _mm_set1_epi64x((1_u64 << (query_length - 1)).cast_signed());
    let mut positive = all;
    let mut negative = _mm_setzero_si128();
    let mut score = _mm_set1_epi64x(i64::try_from(query_length).expect("query length fits i64"));
    let maximum_length = lengths.iter().copied().max().unwrap_or(0);
    for position in 0..maximum_length {
        let mut equal = [0_u64; 2];
        let mut active = [0_i64; 2];
        for lane in 0..starts.len() {
            if position < lengths[lane] {
                equal[lane] =
                    equality_mask(equality_masks, reference_codes[starts[lane] + position]);
                active[lane] = -1;
            }
        }
        let equal = _mm_set_epi64x(equal[1].cast_signed(), equal[0].cast_signed());
        let active = _mm_set_epi64x(active[1], active[0]);
        let horizontal_input = _mm_or_si128(equal, negative);
        let horizontal = _mm_or_si128(
            _mm_xor_si128(
                _mm_add_epi64(_mm_and_si128(equal, positive), positive),
                positive,
            ),
            equal,
        );
        let positive_horizontal = _mm_or_si128(
            negative,
            _mm_andnot_si128(_mm_or_si128(horizontal, positive), all),
        );
        let negative_horizontal = _mm_and_si128(positive, horizontal);
        let positive_high = _mm_and_si128(
            active,
            _mm_cmpeq_epi64(_mm_and_si128(positive_horizontal, high), high),
        );
        let negative_high = _mm_and_si128(
            active,
            _mm_cmpeq_epi64(_mm_and_si128(negative_horizontal, high), high),
        );
        score = _mm_add_epi64(score, _mm_and_si128(positive_high, ones));
        score = _mm_sub_epi64(score, _mm_and_si128(negative_high, ones));
        let shifted_positive = _mm_or_si128(_mm_slli_epi64(positive_horizontal, 1), ones);
        let shifted_negative = _mm_slli_epi64(negative_horizontal, 1);
        let next_positive = _mm_or_si128(
            shifted_negative,
            _mm_andnot_si128(_mm_or_si128(horizontal_input, shifted_positive), all),
        );
        let next_negative = _mm_and_si128(shifted_positive, horizontal_input);
        positive = _mm_or_si128(
            _mm_and_si128(active, next_positive),
            _mm_andnot_si128(active, positive),
        );
        negative = _mm_or_si128(
            _mm_and_si128(active, next_negative),
            _mm_andnot_si128(active, negative),
        );
    }
    let mut scores = [0_u64; 2];
    unsafe {
        _mm_storeu_si128(scores.as_mut_ptr().cast(), score);
    }
    output.copy_from_slice(&scores[..output.len()]);
}

#[cfg(target_arch = "x86_64")]
#[cfg(test)]
#[target_feature(enable = "sse4.2")]
unsafe fn sse42_equal_length_chunk(
    equality_masks: &[u64; 5],
    query_length: usize,
    reference_codes: &[u8],
    starts: &[usize],
    length: usize,
    output: &mut [u64],
) {
    use std::arch::x86_64::{
        _mm_add_epi64, _mm_and_si128, _mm_andnot_si128, _mm_cmpeq_epi64, _mm_or_si128,
        _mm_set_epi64x, _mm_set1_epi64x, _mm_setzero_si128, _mm_slli_epi64, _mm_storeu_si128,
        _mm_sub_epi64, _mm_xor_si128,
    };

    let all = _mm_set1_epi64x(-1);
    let ones = _mm_set1_epi64x(1);
    let high = _mm_set1_epi64x((1_u64 << (query_length - 1)).cast_signed());
    let mut positive = all;
    let mut negative = _mm_setzero_si128();
    let mut score = _mm_set1_epi64x(i64::try_from(query_length).expect("query length fits i64"));
    for position in 0..length {
        let mut equal = [0_u64; 2];
        for lane in 0..starts.len() {
            equal[lane] = equality_mask(equality_masks, reference_codes[starts[lane] + position]);
        }
        let equal = _mm_set_epi64x(equal[1].cast_signed(), equal[0].cast_signed());
        let horizontal_input = _mm_or_si128(equal, negative);
        let horizontal = _mm_or_si128(
            _mm_xor_si128(
                _mm_add_epi64(_mm_and_si128(equal, positive), positive),
                positive,
            ),
            equal,
        );
        let positive_horizontal = _mm_or_si128(
            negative,
            _mm_andnot_si128(_mm_or_si128(horizontal, positive), all),
        );
        let negative_horizontal = _mm_and_si128(positive, horizontal);
        let positive_high = _mm_cmpeq_epi64(_mm_and_si128(positive_horizontal, high), high);
        let negative_high = _mm_cmpeq_epi64(_mm_and_si128(negative_horizontal, high), high);
        score = _mm_add_epi64(score, _mm_and_si128(positive_high, ones));
        score = _mm_sub_epi64(score, _mm_and_si128(negative_high, ones));
        let shifted_positive = _mm_or_si128(_mm_slli_epi64(positive_horizontal, 1), ones);
        let shifted_negative = _mm_slli_epi64(negative_horizontal, 1);
        positive = _mm_or_si128(
            shifted_negative,
            _mm_andnot_si128(_mm_or_si128(horizontal_input, shifted_positive), all),
        );
        negative = _mm_and_si128(shifted_positive, horizontal_input);
    }
    let mut scores = [0_u64; 2];
    unsafe {
        _mm_storeu_si128(scores.as_mut_ptr().cast(), score);
    }
    output.copy_from_slice(&scores[..output.len()]);
}

#[cfg(target_arch = "x86_64")]
#[cfg(test)]
#[target_feature(enable = "avx2")]
unsafe fn avx2_batch(
    equality_masks: &[u64; 5],
    query_length: usize,
    reference_codes: &[u8],
    starts: &[usize],
    lengths: &[usize],
    output: &mut [u64],
) {
    for base in (0..starts.len()).step_by(4) {
        let lanes = (starts.len() - base).min(4);
        unsafe {
            avx2_chunk(
                equality_masks,
                query_length,
                reference_codes,
                &starts[base..base + lanes],
                &lengths[base..base + lanes],
                &mut output[base..base + lanes],
            );
        }
    }
}

#[cfg(target_arch = "x86_64")]
#[cfg(test)]
#[target_feature(enable = "avx2")]
unsafe fn avx2_chunk(
    equality_masks: &[u64; 5],
    query_length: usize,
    reference_codes: &[u8],
    starts: &[usize],
    lengths: &[usize],
    output: &mut [u64],
) {
    use std::arch::x86_64::{
        _mm256_add_epi64, _mm256_and_si256, _mm256_andnot_si256, _mm256_cmpeq_epi64,
        _mm256_or_si256, _mm256_set_epi64x, _mm256_set1_epi64x, _mm256_setzero_si256,
        _mm256_slli_epi64, _mm256_storeu_si256, _mm256_sub_epi64, _mm256_xor_si256,
    };

    if lengths.windows(2).all(|pair| pair[0] == pair[1]) {
        unsafe {
            avx2_equal_length_chunk(
                equality_masks,
                query_length,
                reference_codes,
                starts,
                lengths.first().copied().unwrap_or(0),
                output,
            );
        }
        return;
    }

    let mut positive = _mm256_set1_epi64x(-1);
    let mut negative = _mm256_setzero_si256();
    let mut score = _mm256_set1_epi64x(i64::try_from(query_length).expect("query length fits i64"));
    let high = _mm256_set1_epi64x((1_u64 << (query_length - 1)).cast_signed());
    let ones = _mm256_set1_epi64x(1);
    let maximum_length = lengths.iter().copied().max().unwrap_or(0);
    for position in 0..maximum_length {
        let mut equal = [0_u64; 4];
        let mut active = [0_i64; 4];
        for lane in 0..starts.len() {
            if position < lengths[lane] {
                equal[lane] =
                    equality_mask(equality_masks, reference_codes[starts[lane] + position]);
                active[lane] = -1;
            }
        }
        let equal = _mm256_set_epi64x(
            equal[3].cast_signed(),
            equal[2].cast_signed(),
            equal[1].cast_signed(),
            equal[0].cast_signed(),
        );
        let active = _mm256_set_epi64x(active[3], active[2], active[1], active[0]);
        let horizontal_input = _mm256_or_si256(equal, negative);
        let horizontal = _mm256_or_si256(
            _mm256_xor_si256(
                _mm256_add_epi64(_mm256_and_si256(equal, positive), positive),
                positive,
            ),
            equal,
        );
        let positive_horizontal = _mm256_or_si256(
            negative,
            _mm256_andnot_si256(
                _mm256_or_si256(horizontal, positive),
                _mm256_set1_epi64x(-1),
            ),
        );
        let negative_horizontal = _mm256_and_si256(positive, horizontal);
        let positive_high = _mm256_and_si256(
            active,
            _mm256_cmpeq_epi64(_mm256_and_si256(positive_horizontal, high), high),
        );
        let negative_high = _mm256_and_si256(
            active,
            _mm256_cmpeq_epi64(_mm256_and_si256(negative_horizontal, high), high),
        );
        score = _mm256_add_epi64(score, _mm256_and_si256(positive_high, ones));
        score = _mm256_sub_epi64(score, _mm256_and_si256(negative_high, ones));
        let shifted_positive = _mm256_or_si256(_mm256_slli_epi64(positive_horizontal, 1), ones);
        let shifted_negative = _mm256_slli_epi64(negative_horizontal, 1);
        let next_positive = _mm256_or_si256(
            shifted_negative,
            _mm256_andnot_si256(
                _mm256_or_si256(horizontal_input, shifted_positive),
                _mm256_set1_epi64x(-1),
            ),
        );
        let next_negative = _mm256_and_si256(shifted_positive, horizontal_input);
        positive = _mm256_or_si256(
            _mm256_and_si256(active, next_positive),
            _mm256_andnot_si256(active, positive),
        );
        negative = _mm256_or_si256(
            _mm256_and_si256(active, next_negative),
            _mm256_andnot_si256(active, negative),
        );
    }
    let mut scores = [0_u64; 4];
    unsafe {
        _mm256_storeu_si256(scores.as_mut_ptr().cast(), score);
    }
    output.copy_from_slice(&scores[..output.len()]);
}

#[cfg(target_arch = "x86_64")]
#[cfg(test)]
#[target_feature(enable = "avx2")]
unsafe fn avx2_equal_length_chunk(
    equality_masks: &[u64; 5],
    query_length: usize,
    reference_codes: &[u8],
    starts: &[usize],
    length: usize,
    output: &mut [u64],
) {
    use std::arch::x86_64::{
        _mm256_add_epi64, _mm256_and_si256, _mm256_andnot_si256, _mm256_cmpeq_epi64,
        _mm256_or_si256, _mm256_set_epi64x, _mm256_set1_epi64x, _mm256_setzero_si256,
        _mm256_slli_epi64, _mm256_storeu_si256, _mm256_sub_epi64, _mm256_xor_si256,
    };

    let all = _mm256_set1_epi64x(-1);
    let ones = _mm256_set1_epi64x(1);
    let high = _mm256_set1_epi64x((1_u64 << (query_length - 1)).cast_signed());
    let mut positive = all;
    let mut negative = _mm256_setzero_si256();
    let mut score = _mm256_set1_epi64x(i64::try_from(query_length).expect("query length fits i64"));
    for position in 0..length {
        let mut equal = [0_u64; 4];
        for lane in 0..starts.len() {
            equal[lane] = equality_mask(equality_masks, reference_codes[starts[lane] + position]);
        }
        let equal = _mm256_set_epi64x(
            equal[3].cast_signed(),
            equal[2].cast_signed(),
            equal[1].cast_signed(),
            equal[0].cast_signed(),
        );
        let horizontal_input = _mm256_or_si256(equal, negative);
        let horizontal = _mm256_or_si256(
            _mm256_xor_si256(
                _mm256_add_epi64(_mm256_and_si256(equal, positive), positive),
                positive,
            ),
            equal,
        );
        let positive_horizontal = _mm256_or_si256(
            negative,
            _mm256_andnot_si256(_mm256_or_si256(horizontal, positive), all),
        );
        let negative_horizontal = _mm256_and_si256(positive, horizontal);
        let positive_high = _mm256_cmpeq_epi64(_mm256_and_si256(positive_horizontal, high), high);
        let negative_high = _mm256_cmpeq_epi64(_mm256_and_si256(negative_horizontal, high), high);
        score = _mm256_add_epi64(score, _mm256_and_si256(positive_high, ones));
        score = _mm256_sub_epi64(score, _mm256_and_si256(negative_high, ones));
        let shifted_positive = _mm256_or_si256(_mm256_slli_epi64(positive_horizontal, 1), ones);
        let shifted_negative = _mm256_slli_epi64(negative_horizontal, 1);
        positive = _mm256_or_si256(
            shifted_negative,
            _mm256_andnot_si256(_mm256_or_si256(horizontal_input, shifted_positive), all),
        );
        negative = _mm256_and_si256(shifted_positive, horizontal_input);
    }
    let mut scores = [0_u64; 4];
    unsafe {
        _mm256_storeu_si256(scores.as_mut_ptr().cast(), score);
    }
    output.copy_from_slice(&scores[..output.len()]);
}

#[cfg(target_arch = "x86_64")]
#[cfg(test)]
#[target_feature(enable = "avx512f,avx512bw")]
unsafe fn avx512_batch(
    equality_masks: &[u64; 5],
    query_length: usize,
    reference_codes: &[u8],
    starts: &[usize],
    lengths: &[usize],
    output: &mut [u64],
) {
    for base in (0..starts.len()).step_by(8) {
        let lanes = (starts.len() - base).min(8);
        unsafe {
            avx512_chunk(
                equality_masks,
                query_length,
                reference_codes,
                &starts[base..base + lanes],
                &lengths[base..base + lanes],
                &mut output[base..base + lanes],
            );
        }
    }
}

#[cfg(target_arch = "x86_64")]
#[cfg(test)]
#[target_feature(enable = "avx512f,avx512bw")]
unsafe fn avx512_chunk(
    equality_masks: &[u64; 5],
    query_length: usize,
    reference_codes: &[u8],
    starts: &[usize],
    lengths: &[usize],
    output: &mut [u64],
) {
    use std::arch::x86_64::{
        _mm512_add_epi64, _mm512_and_si512, _mm512_cmpeq_epi64_mask, _mm512_mask_add_epi64,
        _mm512_mask_mov_epi64, _mm512_mask_sub_epi64, _mm512_or_si512, _mm512_set_epi64,
        _mm512_set1_epi64, _mm512_setzero_si512, _mm512_slli_epi64, _mm512_storeu_si512,
        _mm512_xor_si512,
    };

    let mut positive = _mm512_set1_epi64(-1);
    let mut negative = _mm512_setzero_si512();
    let mut score = _mm512_set1_epi64(i64::try_from(query_length).expect("query length fits i64"));
    let high = _mm512_set1_epi64((1_u64 << (query_length - 1)).cast_signed());
    let ones = _mm512_set1_epi64(1);
    let maximum_length = lengths.iter().copied().max().unwrap_or(0);
    for position in 0..maximum_length {
        let mut equal = [0_u64; 8];
        let mut active = 0_u8;
        for lane in 0..starts.len() {
            if position < lengths[lane] {
                equal[lane] =
                    equality_mask(equality_masks, reference_codes[starts[lane] + position]);
                active |= 1 << lane;
            }
        }
        let equal = _mm512_set_epi64(
            equal[7].cast_signed(),
            equal[6].cast_signed(),
            equal[5].cast_signed(),
            equal[4].cast_signed(),
            equal[3].cast_signed(),
            equal[2].cast_signed(),
            equal[1].cast_signed(),
            equal[0].cast_signed(),
        );
        let horizontal_input = _mm512_or_si512(equal, negative);
        let horizontal = _mm512_or_si512(
            _mm512_xor_si512(
                _mm512_add_epi64(_mm512_and_si512(equal, positive), positive),
                positive,
            ),
            equal,
        );
        let positive_horizontal = _mm512_or_si512(
            negative,
            _mm512_xor_si512(_mm512_or_si512(horizontal, positive), _mm512_set1_epi64(-1)),
        );
        let negative_horizontal = _mm512_and_si512(positive, horizontal);
        let positive_high =
            _mm512_cmpeq_epi64_mask(_mm512_and_si512(positive_horizontal, high), high) & active;
        let negative_high =
            _mm512_cmpeq_epi64_mask(_mm512_and_si512(negative_horizontal, high), high) & active;
        score = _mm512_mask_add_epi64(score, positive_high, score, ones);
        score = _mm512_mask_sub_epi64(score, negative_high, score, ones);
        let shifted_positive = _mm512_or_si512(_mm512_slli_epi64(positive_horizontal, 1), ones);
        let shifted_negative = _mm512_slli_epi64(negative_horizontal, 1);
        let next_positive = _mm512_or_si512(
            shifted_negative,
            _mm512_xor_si512(
                _mm512_or_si512(horizontal_input, shifted_positive),
                _mm512_set1_epi64(-1),
            ),
        );
        let next_negative = _mm512_and_si512(shifted_positive, horizontal_input);
        positive = _mm512_mask_mov_epi64(positive, active, next_positive);
        negative = _mm512_mask_mov_epi64(negative, active, next_negative);
    }
    let mut scores = [0_u64; 8];
    unsafe {
        _mm512_storeu_si512(scores.as_mut_ptr().cast(), score);
    }
    output.copy_from_slice(&scores[..output.len()]);
}

#[cfg(test)]
#[path = "../../tests/whitebox/alignment_kernels.rs"]
mod tests;
