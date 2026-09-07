//! x86-64 AVX-512BW candidate-parallel narrow-band kernels.
//!
//! The AVX-512 backend deliberately reuses the mature AVX2 kernels for
//! operations that are not candidate-parallel. These kernels widen the two
//! batch primitives whose independent candidate lanes benefit directly from
//! 512-bit byte and word vectors.

#![allow(unsafe_code)]
#![deny(unsafe_op_in_unsafe_fn)]

use super::{
    MAX_NARROW_ENDPOINTS, NarrowBandedResult, NarrowEndpointDistances, finish_narrow_result,
};

#[target_feature(enable = "avx512f,avx512bw,popcnt")]
// Keep the relation table borrowed across scalar/SIMD dispatch.
#[allow(clippy::trivially_copy_pass_by_ref)]
pub(super) unsafe fn narrow_fixed_start_avx512(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    patterns: &[u8],
    pattern_length: usize,
    max_distance: usize,
    output: &mut [NarrowEndpointDistances],
) {
    for base in (0..output.len()).step_by(64) {
        let lanes = (output.len() - base).min(64);
        unsafe {
            narrow_fixed_start_avx512_chunk(
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

#[target_feature(enable = "avx512f,avx512bw,popcnt")]
// Bounds validation makes the byte casts exact; the pointer casts are used
// only by unaligned AVX-512 intrinsics.
#[allow(clippy::cast_possible_truncation, clippy::cast_ptr_alignment)]
unsafe fn narrow_fixed_start_avx512_chunk(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    patterns: &[u8],
    pattern_length: usize,
    lanes: usize,
    max_distance: usize,
    output: &mut [NarrowEndpointDistances],
) {
    use std::arch::x86_64::{
        __m512i, _mm512_adds_epu8, _mm512_loadu_si512, _mm512_min_epu8, _mm512_set1_epi8,
        _mm512_storeu_si512,
    };

    let band_length = max_distance * 2 + 1;
    let capped = u8::try_from(max_distance)
        .unwrap_or(u8::MAX - 1)
        .saturating_add(1);
    let cap_vector = _mm512_set1_epi8(capped.cast_signed());
    let one = _mm512_set1_epi8(1);
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
        *slot = _mm512_set1_epi8(initial.cast_signed());
    }
    for (query_position, &query_code) in query.iter().enumerate() {
        current[..band_length].fill(cap_vector);
        let reference_mask = reference_masks_by_query
            .get(usize::from(query_code))
            .copied()
            .unwrap_or(0);
        for diagonal in 0..band_length {
            let reference_position = query_position + diagonal;
            let mut substitutions = [capped; 64];
            let mut lane = 0;
            while lane + 4 <= lanes {
                for offset in 0..4 {
                    let lane = lane + offset;
                    let reference_code = patterns[lane * pattern_length + reference_position];
                    substitutions[lane] = u8::from(
                        reference_code >= u8::BITS as u8
                            || reference_mask & (1_u8 << reference_code) == 0,
                    );
                }
                lane += 4;
            }
            while lane < lanes {
                let reference_code = patterns[lane * pattern_length + reference_position];
                substitutions[lane] = u8::from(
                    reference_code >= u8::BITS as u8
                        || reference_mask & (1_u8 << reference_code) == 0,
                );
                lane += 1;
            }
            let substitution =
                unsafe { _mm512_loadu_si512(substitutions.as_ptr().cast::<__m512i>()) };
            let diagonal_score = _mm512_adds_epu8(previous[diagonal], substitution);
            let query_gap = if diagonal + 1 < band_length {
                _mm512_adds_epu8(previous[diagonal + 1], one)
            } else {
                cap_vector
            };
            let reference_gap = if diagonal != 0 {
                _mm512_adds_epu8(current[diagonal - 1], one)
            } else {
                cap_vector
            };
            current[diagonal] = _mm512_min_epu8(
                _mm512_min_epu8(_mm512_min_epu8(diagonal_score, query_gap), reference_gap),
                cap_vector,
            );
        }
        core::mem::swap(&mut previous, &mut current);
    }

    output.fill(NarrowEndpointDistances::EMPTY);
    for (endpoint, distances) in previous[..band_length].iter().copied().enumerate() {
        let mut lanes_out = [capped; 64];
        unsafe {
            _mm512_storeu_si512(lanes_out.as_mut_ptr().cast::<__m512i>(), distances);
        }
        for lane in 0..lanes {
            let distance = lanes_out[lane];
            if usize::from(distance) <= max_distance {
                output[lane].insert(endpoint, u32::from(distance));
            }
        }
    }
}

#[target_feature(enable = "avx512f,avx512bw,popcnt")]
// Keep the relation table borrowed across scalar/SIMD dispatch.
#[allow(clippy::trivially_copy_pass_by_ref)]
pub(super) unsafe fn narrow_banded_avx512(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    patterns: &[u8],
    pattern_length: usize,
    max_distance: usize,
    output: &mut [NarrowBandedResult],
) {
    for base in (0..output.len()).step_by(16) {
        let lanes = (output.len() - base).min(16);
        unsafe {
            narrow_banded_avx512_chunk(
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

#[target_feature(enable = "avx512f,avx512bw,popcnt")]
// Unaligned AVX-512 loads/stores intentionally accept lane-array pointers.
#[allow(clippy::cast_ptr_alignment)]
unsafe fn narrow_banded_avx512_chunk(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    patterns: &[u8],
    pattern_length: usize,
    lanes: usize,
    max_distance: usize,
    output: &mut [NarrowBandedResult],
) {
    use std::arch::x86_64::{
        __m512i, _mm512_add_epi32, _mm512_and_si512, _mm512_andnot_si512, _mm512_loadu_si512,
        _mm512_or_si512, _mm512_set1_epi32, _mm512_setzero_si512, _mm512_srli_epi32,
        _mm512_storeu_si512, _mm512_sub_epi32, _mm512_xor_si512,
    };

    let band_length = max_distance * 2 + 1;
    let high_bit = 1_u32 << (band_length - 1);
    let mut initial = [[0_u32; 16]; 5];
    for lane in 0..lanes {
        let pattern = &patterns[lane * pattern_length..(lane + 1) * pattern_length];
        for (position, &reference_code) in pattern[..band_length].iter().enumerate() {
            if let Some(bits) = initial.get_mut(usize::from(reference_code)) {
                bits[lane] |= 1_u32 << position;
            }
        }
    }
    let mut peq: [__m512i; 5] = core::array::from_fn(|code| unsafe {
        _mm512_loadu_si512(initial[code].as_ptr().cast::<__m512i>())
    });
    let all = _mm512_set1_epi32(-1);
    let one = _mm512_set1_epi32(1);
    let mut positive = _mm512_setzero_si512();
    let mut negative = _mm512_setzero_si512();
    let mut error = _mm512_setzero_si512();
    for (query_position, &query_code) in query.iter().enumerate() {
        let reference_mask = reference_masks_by_query
            .get(usize::from(query_code))
            .copied()
            .unwrap_or(0);
        let mut equal = _mm512_setzero_si512();
        for (reference_code, &bits) in peq.iter().enumerate() {
            if reference_mask & (1_u8 << reference_code) != 0 {
                equal = _mm512_or_si512(equal, bits);
            }
        }
        let horizontal_input = _mm512_or_si512(equal, negative);
        let horizontal = _mm512_or_si512(
            _mm512_xor_si512(
                _mm512_add_epi32(_mm512_and_si512(horizontal_input, positive), positive),
                positive,
            ),
            horizontal_input,
        );
        let negative_horizontal = _mm512_and_si512(positive, horizontal);
        let positive_horizontal = _mm512_or_si512(
            negative,
            _mm512_andnot_si512(_mm512_or_si512(positive, horizontal), all),
        );
        let shifted = _mm512_srli_epi32(horizontal, 1);
        negative = _mm512_and_si512(shifted, positive_horizontal);
        positive = _mm512_or_si512(
            negative_horizontal,
            _mm512_andnot_si512(_mm512_or_si512(shifted, positive_horizontal), all),
        );
        error = _mm512_sub_epi32(
            _mm512_add_epi32(error, one),
            _mm512_and_si512(horizontal, one),
        );
        if query_position + 1 != query.len() {
            for bits in &mut peq {
                *bits = _mm512_srli_epi32(*bits, 1);
            }
            let entering_position = band_length + query_position;
            let mut entering = [[0_u32; 16]; 5];
            for lane in 0..lanes {
                let reference_code = patterns[lane * pattern_length + entering_position];
                if let Some(bits) = entering.get_mut(usize::from(reference_code)) {
                    bits[lane] = high_bit;
                }
            }
            for (bits, additions) in peq.iter_mut().zip(&entering) {
                *bits = _mm512_or_si512(*bits, unsafe {
                    _mm512_loadu_si512(additions.as_ptr().cast::<__m512i>())
                });
            }
        }
    }
    let mut errors = [0_u32; 16];
    let mut positives = [0_u32; 16];
    let mut negatives = [0_u32; 16];
    unsafe {
        _mm512_storeu_si512(errors.as_mut_ptr().cast::<__m512i>(), error);
        _mm512_storeu_si512(positives.as_mut_ptr().cast::<__m512i>(), positive);
        _mm512_storeu_si512(negatives.as_mut_ptr().cast::<__m512i>(), negative);
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
