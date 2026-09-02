//! Portable scalar implementations for narrow verification.

use super::{
    MAX_NARROW_ENDPOINTS, NarrowBandedResult, NarrowEndpointDistances, NarrowPlacementDistances,
    NarrowReferenceCode, finish_narrow_result,
};

// This table is shared by scalar and SIMD dispatch; retain one borrowed
// signature for the paired differential-test surface.
#[allow(clippy::trivially_copy_pass_by_ref)]
pub(in crate::verification) fn narrow_banded_scalar_one(
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

// This table is shared by scalar and SIMD dispatch; retain one borrowed
// signature for the paired differential-test surface.
#[allow(clippy::trivially_copy_pass_by_ref)]
pub(in crate::verification) fn narrow_placement_distances_scalar(
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

#[allow(clippy::trivially_copy_pass_by_ref)]
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
#[allow(clippy::trivially_copy_pass_by_ref)]
pub(in crate::verification) fn narrow_fixed_start_scalar_one(
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
pub(super) fn narrow_fixed_start_gather_scalar_one<T: NarrowReferenceCode>(
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
pub(in crate::verification) fn myers_prefix_distances_u128_scalar_one(
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
