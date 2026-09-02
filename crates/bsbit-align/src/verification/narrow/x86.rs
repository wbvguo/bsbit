//! Runtime-detected x86 narrow-verification kernels.

use super::{
    MAX_NARROW_ENDPOINTS, NarrowBandedResult, NarrowEndpointDistances, NarrowPlacementDistances,
    NarrowReferenceCode, finish_narrow_result,
};

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "sse4.2")]
// Unaligned SSE loads/stores intentionally accept byte-array pointers; the
// intrinsic provides the alignment guarantee that a typed dereference would.
#[allow(clippy::cast_ptr_alignment, clippy::needless_range_loop)]
pub(in crate::verification) unsafe fn narrow_placement_distances_sse42(
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
pub(in crate::verification) unsafe fn narrow_placement_distances_avx2(
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
pub(super) unsafe fn narrow_placement_distances_d3_avx2(
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
        // Every current diagonal is overwritten in dependency order; filling
        // the inactive buffer here only adds seven vector stores per base.
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
pub(super) unsafe fn narrow_placement_distances_d5_avx2(
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
        // Every current diagonal is overwritten in dependency order; filling
        // the inactive buffer here only adds eleven vector stores per base.
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
pub(super) unsafe fn narrow_placement_distances_batch_avx2(
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
pub(super) unsafe fn narrow_placement_distances_batch_d3_avx2(
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
        // Every current diagonal is overwritten in dependency order; filling
        // the inactive buffer here only adds seven vector stores per base.
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
#[allow(
    clippy::cast_ptr_alignment,
    clippy::cast_possible_truncation,
    clippy::too_many_lines
)]
pub(super) unsafe fn narrow_placement_distances_interleaved_batch_d3_avx2(
    reference_masks_by_query: &[u8; 5],
    query: &[u8],
    patterns: &[u8],
    output: &mut [NarrowPlacementDistances],
) {
    use std::arch::x86_64::{
        __m128i, __m256i, _mm_cvtsi32_si128, _mm_loadu_si128, _mm_movemask_epi8, _mm_setzero_si128,
        _mm_shuffle_epi8, _mm256_adds_epu8, _mm256_loadu_si256, _mm256_min_epu8, _mm256_or_si256,
        _mm256_set1_epi8, _mm256_setzero_si256, _mm256_storeu_si256,
    };

    const BAND_LENGTH: usize = 7;
    const MAX_DISTANCE: usize = 3;
    const CAPPED: u8 = 4;
    const PATTERN_LANES: usize = 4;

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
        *slot = unsafe { _mm256_loadu_si256(lanes.as_ptr().cast::<__m256i>()) };
    }
    let mut previous = &mut previous_storage;
    let mut current = &mut current_storage;
    let one = _mm256_set1_epi8(1);
    let lane_masks: [__m256i; 4] = core::array::from_fn(|candidate| {
        let mut lanes = [0_u8; 32];
        lanes[candidate * BAND_LENGTH..(candidate + 1) * BAND_LENGTH].fill(1);
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
    let relation_tables: [__m128i; 5] = core::array::from_fn(|query_code| {
        let reference_mask = reference_masks_by_query[query_code];
        let mut table = [0_i8; 16];
        for (reference_code, relation) in table[..5].iter_mut().enumerate() {
            if reference_mask & (1_u8 << reference_code) != 0 {
                *relation = -1;
            }
        }
        unsafe { _mm_loadu_si128(table.as_ptr().cast::<__m128i>()) }
    });
    let active_candidates = (1_usize << output.len()) - 1;
    for (query_position, &query_code) in query.iter().enumerate() {
        current.fill(cap_vector);
        let relation_table = relation_tables
            .get(usize::from(query_code))
            .copied()
            .unwrap_or_else(|| _mm_setzero_si128());
        for diagonal in 0..BAND_LENGTH {
            let offset = (query_position + diagonal) * PATTERN_LANES;
            let packed_codes =
                unsafe { patterns.as_ptr().add(offset).cast::<i32>().read_unaligned() };
            let codes = _mm_cvtsi32_si128(packed_codes);
            let matches = _mm_shuffle_epi8(relation_table, codes);
            let mismatch_mask = (!(u32::try_from(_mm_movemask_epi8(matches)).unwrap_or(0))
                as usize)
                & active_candidates;
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
pub(super) unsafe fn narrow_placement_distances_batch_d5_avx2(
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
        // Every current diagonal is overwritten in dependency order; filling
        // the inactive buffer here only adds eleven vector stores per base.
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

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
// Keep the equality table borrowed across scalar/SIMD dispatch.
#[allow(clippy::trivially_copy_pass_by_ref)]
pub(super) unsafe fn myers_prefix_distances_u128_avx2(
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
pub(in crate::verification) unsafe fn narrow_fixed_start_sse42(
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
pub(super) unsafe fn narrow_fixed_start_avx2(
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
pub(super) unsafe fn narrow_fixed_start_gather_avx2<T: NarrowReferenceCode>(
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
pub(in crate::verification) fn narrow_sse42_available() -> bool {
    std::arch::is_x86_feature_detected!("sse4.2")
}

#[cfg(not(target_arch = "x86_64"))]
pub(in crate::verification) const fn narrow_sse42_available() -> bool {
    false
}

#[cfg(target_arch = "x86_64")]
pub(in crate::verification) fn narrow_avx2_available() -> bool {
    std::arch::is_x86_feature_detected!("avx2")
}

#[cfg(not(target_arch = "x86_64"))]
pub(in crate::verification) const fn narrow_avx2_available() -> bool {
    false
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "sse4.2")]
// Keep the relation table borrowed across scalar/SIMD dispatch.
#[allow(clippy::trivially_copy_pass_by_ref)]
pub(super) unsafe fn narrow_banded_sse42(
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
pub(super) unsafe fn narrow_banded_avx2(
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
