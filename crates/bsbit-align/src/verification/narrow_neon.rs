//! `AArch64` NEON implementations for narrow-band verification.
//!
//! This child module owns only algorithm-specific intrinsics. Safe callers in
//! `narrow` validate dimensions and install these functions only for the
//! process-wide `AArch64` NEON backend.

#![allow(unsafe_code)]
#![deny(unsafe_op_in_unsafe_fn)]

use super::{
    MAX_NARROW_ENDPOINTS, NarrowBandedResult, NarrowEndpointDistances, NarrowPlacementDistances,
    finish_narrow_result,
};

/// Computes the complete start/end frontier in sixteen NEON byte lanes.
///
/// # Safety
///
/// The caller must run on `AArch64` with NEON available and must satisfy the
/// validated narrow-band dimensions enforced by the public wrapper.
#[target_feature(enable = "neon")]
#[allow(clippy::needless_range_loop)]
pub(super) unsafe fn narrow_placement_distances(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    pattern: &[u8],
    max_distance: usize,
) -> NarrowPlacementDistances {
    use std::arch::aarch64::{uint8x16_t, vdupq_n_u8, vld1q_u8, vminq_u8, vqaddq_u8, vst1q_u8};

    let band_length = max_distance * 2 + 1;
    let capped = u8::try_from(max_distance)
        .unwrap_or(u8::MAX - 1)
        .saturating_add(1);
    let cap_vector = vdupq_n_u8(capped);
    let one = vdupq_n_u8(1);
    let mut distances = NarrowPlacementDistances::for_band(band_length, max_distance);

    for start_base in (0..band_length).step_by(16) {
        let active_starts = (band_length - start_base).min(16);
        let mut previous: [uint8x16_t; MAX_NARROW_ENDPOINTS] = [cap_vector; MAX_NARROW_ENDPOINTS];
        let mut current: [uint8x16_t; MAX_NARROW_ENDPOINTS] = [cap_vector; MAX_NARROW_ENDPOINTS];
        for (diagonal, slot) in previous[..band_length].iter_mut().enumerate() {
            let mut lanes = [capped; 16];
            for (lane, value) in lanes[..active_starts].iter_mut().enumerate() {
                let start = start_base + lane;
                if start <= diagonal {
                    *value = u8::try_from(diagonal - start).unwrap_or(capped).min(capped);
                }
            }
            // SAFETY: one vector reads exactly the live sixteen-byte array.
            *slot = unsafe { vld1q_u8(lanes.as_ptr()) };
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
                let diagonal_score = vqaddq_u8(previous[diagonal], vdupq_n_u8(substitution));
                let query_gap = if diagonal + 1 < band_length {
                    vqaddq_u8(previous[diagonal + 1], one)
                } else {
                    cap_vector
                };
                let reference_gap = if diagonal != 0 {
                    vqaddq_u8(current[diagonal - 1], one)
                } else {
                    cap_vector
                };
                current[diagonal] = vminq_u8(
                    vminq_u8(vminq_u8(diagonal_score, query_gap), reference_gap),
                    cap_vector,
                );
            }
            core::mem::swap(&mut previous, &mut current);
        }

        for (endpoint, &endpoint_distances) in previous[..band_length].iter().enumerate() {
            let mut lanes = [capped; 16];
            // SAFETY: one vector writes exactly the live sixteen-byte array.
            unsafe { vst1q_u8(lanes.as_mut_ptr(), endpoint_distances) };
            for (lane, &distance) in lanes[..active_starts].iter().enumerate() {
                distances.insert_distance(start_base + lane, endpoint, distance);
            }
        }
    }

    distances
}

/// Computes fixed-start endpoint frontiers for sixteen candidates at a time.
///
/// # Safety
///
/// The caller must run on `AArch64` with NEON available and provide dimensions
/// already validated by `narrow_banded_fixed_start_batch`.
#[target_feature(enable = "neon")]
pub(super) unsafe fn narrow_fixed_start_batch(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    patterns: &[u8],
    pattern_length: usize,
    max_distance: usize,
    output: &mut [NarrowEndpointDistances],
) {
    for base in (0..output.len()).step_by(16) {
        let lanes = (output.len() - base).min(16);
        // SAFETY: the caller validated every flattened candidate range and
        // this chunk contains exactly `lanes` complete patterns and outputs.
        unsafe {
            narrow_fixed_start_chunk(
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

#[target_feature(enable = "neon")]
unsafe fn narrow_fixed_start_chunk(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    patterns: &[u8],
    pattern_length: usize,
    lanes: usize,
    max_distance: usize,
    output: &mut [NarrowEndpointDistances],
) {
    use std::arch::aarch64::{uint8x16_t, vdupq_n_u8, vld1q_u8, vminq_u8, vqaddq_u8, vst1q_u8};

    let band_length = max_distance * 2 + 1;
    let capped = u8::try_from(max_distance)
        .unwrap_or(u8::MAX - 1)
        .saturating_add(1);
    let cap_vector = vdupq_n_u8(capped);
    let one = vdupq_n_u8(1);
    let mut previous: [uint8x16_t; MAX_NARROW_ENDPOINTS] = [cap_vector; MAX_NARROW_ENDPOINTS];
    let mut current: [uint8x16_t; MAX_NARROW_ENDPOINTS] = [cap_vector; MAX_NARROW_ENDPOINTS];
    for (diagonal, slot) in previous[..band_length].iter_mut().enumerate() {
        let initial = if diagonal < max_distance {
            capped
        } else {
            u8::try_from(diagonal - max_distance)
                .unwrap_or(capped)
                .min(capped)
        };
        *slot = vdupq_n_u8(initial);
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
                substitutions[lane] =
                    u8::from(reference_code >= 8 || reference_mask & (1_u8 << reference_code) == 0);
            }
            // SAFETY: the local array supplies one complete vector.
            let substitution = unsafe { vld1q_u8(substitutions.as_ptr()) };
            let diagonal_score = vqaddq_u8(previous[diagonal], substitution);
            let query_gap = if diagonal + 1 < band_length {
                vqaddq_u8(previous[diagonal + 1], one)
            } else {
                cap_vector
            };
            let reference_gap = if diagonal != 0 {
                vqaddq_u8(current[diagonal - 1], one)
            } else {
                cap_vector
            };
            current[diagonal] = vminq_u8(
                vminq_u8(vminq_u8(diagonal_score, query_gap), reference_gap),
                cap_vector,
            );
        }
        core::mem::swap(&mut previous, &mut current);
    }

    output.fill(NarrowEndpointDistances::EMPTY);
    for (endpoint, distances) in previous[..band_length].iter().copied().enumerate() {
        let mut lanes_out = [capped; 16];
        // SAFETY: the local array receives one complete vector.
        unsafe { vst1q_u8(lanes_out.as_mut_ptr(), distances) };
        for lane in 0..lanes {
            let distance = lanes_out[lane];
            if usize::from(distance) <= max_distance {
                output[lane].insert(endpoint, u32::from(distance));
            }
        }
    }
}

/// Computes prefix frontiers for four candidates at a time in 32-bit lanes.
///
/// # Safety
///
/// The caller must run on `AArch64` with NEON available and provide dimensions
/// already validated by `narrow_banded_prefix_batch`.
#[target_feature(enable = "neon")]
pub(super) unsafe fn narrow_banded_prefix_batch(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    patterns: &[u8],
    pattern_length: usize,
    max_distance: usize,
    output: &mut [NarrowBandedResult],
) {
    for base in (0..output.len()).step_by(4) {
        let lanes = (output.len() - base).min(4);
        // SAFETY: the caller validated every flattened candidate range and
        // this chunk contains exactly `lanes` complete patterns and outputs.
        unsafe {
            narrow_banded_prefix_chunk(
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

#[target_feature(enable = "neon")]
unsafe fn narrow_banded_prefix_chunk(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    patterns: &[u8],
    pattern_length: usize,
    lanes: usize,
    max_distance: usize,
    output: &mut [NarrowBandedResult],
) {
    use std::arch::aarch64::{
        uint32x4_t, vaddq_u32, vandq_u32, vdupq_n_u32, veorq_u32, vld1q_u32, vmvnq_u32, vorrq_u32,
        vshrq_n_u32, vst1q_u32, vsubq_u32,
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
    let mut peq: [uint32x4_t; 5] = core::array::from_fn(|code| {
        // SAFETY: each row is one complete four-lane vector.
        unsafe { vld1q_u32(initial[code].as_ptr()) }
    });
    let one = vdupq_n_u32(1);
    let mut positive = vdupq_n_u32(0);
    let mut negative = vdupq_n_u32(0);
    let mut error = vdupq_n_u32(0);
    for (query_position, &query_code) in query.iter().enumerate() {
        let reference_mask = reference_masks_by_query
            .get(usize::from(query_code))
            .copied()
            .unwrap_or(0);
        let mut equal = vdupq_n_u32(0);
        for (reference_code, &bits) in peq.iter().enumerate() {
            if reference_mask & (1_u8 << reference_code) != 0 {
                equal = vorrq_u32(equal, bits);
            }
        }
        let horizontal_input = vorrq_u32(equal, negative);
        let horizontal = vorrq_u32(
            veorq_u32(
                vaddq_u32(vandq_u32(horizontal_input, positive), positive),
                positive,
            ),
            horizontal_input,
        );
        let negative_horizontal = vandq_u32(positive, horizontal);
        let positive_horizontal = vorrq_u32(negative, vmvnq_u32(vorrq_u32(positive, horizontal)));
        let shifted = vshrq_n_u32::<1>(horizontal);
        negative = vandq_u32(shifted, positive_horizontal);
        positive = vorrq_u32(
            negative_horizontal,
            vmvnq_u32(vorrq_u32(shifted, positive_horizontal)),
        );
        error = vsubq_u32(vaddq_u32(error, one), vandq_u32(horizontal, one));
        if query_position + 1 != query.len() {
            for bits in &mut peq {
                *bits = vshrq_n_u32::<1>(*bits);
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
                // SAFETY: each row is one complete four-lane vector.
                *bits = vorrq_u32(*bits, unsafe { vld1q_u32(additions.as_ptr()) });
            }
        }
    }

    let mut errors = [0_u32; 4];
    let mut positives = [0_u32; 4];
    let mut negatives = [0_u32; 4];
    // SAFETY: each output row receives one complete four-lane vector.
    unsafe {
        vst1q_u32(errors.as_mut_ptr(), error);
        vst1q_u32(positives.as_mut_ptr(), positive);
        vst1q_u32(negatives.as_mut_ptr(), negative);
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
