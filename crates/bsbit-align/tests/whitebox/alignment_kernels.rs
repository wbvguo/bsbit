//! White-box tests for the alignment-kernel implementation.
//!
//! Kept outside implementation `src/` while remaining a child module so private
//! invariants can be tested without widening the crate API.

use super::*;

const LITERAL_MASKS: [u64; 5] = [
    0x1111_1111_1111_1111,
    0x2222_2222_2222_2222,
    0x4444_4444_4444_4444,
    0x8888_8888_8888_8888,
    0,
];

fn brute_narrow_result(
    reference_masks_by_query: [u8; 5],
    query: &[u8],
    pattern: &[u8],
    max_distance: usize,
) -> NarrowBandedResult {
    let infinity = u32::MAX / 4;
    let maximum_lead = max_distance * 2;
    let mut matrix = vec![vec![infinity; pattern.len() + 1]; query.len() + 1];
    // The narrow-band contract gives the reference a free 0..=2k lead, then retains
    // only cells whose reference/query prefix displacement stays in that
    // same band. This matrix is intentionally independent of the bit
    // recurrence used by both optimized implementations above.
    for cell in matrix[0].iter_mut().take(maximum_lead + 1) {
        *cell = 0;
    }
    for query_position in 1..=query.len() {
        let first_reference = query_position;
        let last_reference = (query_position + maximum_lead).min(pattern.len());
        let query_code = query[query_position - 1];
        let reference_mask = reference_masks_by_query
            .get(usize::from(query_code))
            .copied()
            .unwrap_or(0);
        for reference_position in first_reference..=last_reference {
            let reference_code = pattern[reference_position - 1];
            let substitution =
                u32::from(reference_mask & (1_u8 << usize::from(reference_code)) == 0);
            matrix[query_position][reference_position] = matrix[query_position - 1]
                [reference_position]
                .saturating_add(1)
                .min(matrix[query_position][reference_position - 1].saturating_add(1))
                .min(
                    matrix[query_position - 1][reference_position - 1].saturating_add(substitution),
                );
        }
    }
    let mut best = NarrowBandedResult::ABSENT;
    let mut endpoint_distances = NarrowEndpointDistances::EMPTY;
    let mut center_distance = u32::MAX;
    for delta in 0..=max_distance * 2 {
        let prefix_length = query.len() + delta;
        let distance = matrix[query.len()][prefix_length];
        if delta == max_distance {
            center_distance = distance;
        }
        if usize::try_from(distance).unwrap_or(usize::MAX) <= max_distance {
            endpoint_distances.insert(delta, distance);
            if distance < best.distance {
                best = NarrowBandedResult {
                    distance,
                    prefix_length,
                    tied_prefix_mask: 1_u32 << delta,
                    endpoint_distances: NarrowEndpointDistances::EMPTY,
                };
            } else if distance == best.distance {
                best.prefix_length = prefix_length;
                best.tied_prefix_mask |= 1_u32 << delta;
            }
        }
    }
    if center_distance <= best.distance
        && usize::try_from(center_distance).unwrap_or(usize::MAX) <= max_distance
    {
        best = NarrowBandedResult {
            distance: center_distance,
            prefix_length: query.len() + max_distance,
            tied_prefix_mask: best.tied_prefix_mask,
            endpoint_distances: NarrowEndpointDistances::EMPTY,
        };
    }
    best.endpoint_distances = endpoint_distances;
    best
}

fn brute_unrestricted_semiglobal_distance(query: &[u8], pattern: &[u8]) -> u32 {
    let mut matrix = vec![vec![0_u32; pattern.len() + 1]; query.len() + 1];
    for (query_position, row) in matrix.iter_mut().enumerate().skip(1) {
        row[0] = u32::try_from(query_position).unwrap_or(u32::MAX);
    }
    for query_position in 1..=query.len() {
        for reference_position in 1..=pattern.len() {
            let substitution =
                u32::from(query[query_position - 1] != pattern[reference_position - 1]);
            matrix[query_position][reference_position] = matrix[query_position - 1]
                [reference_position]
                .saturating_add(1)
                .min(matrix[query_position][reference_position - 1].saturating_add(1))
                .min(
                    matrix[query_position - 1][reference_position - 1].saturating_add(substitution),
                );
        }
    }
    matrix[query.len()][query.len()..]
        .iter()
        .copied()
        .min()
        .unwrap_or(u32::MAX)
}

#[test]
fn validation_is_fail_closed_before_output() {
    let reference = [0_u8, 1, 2, 3];
    let mut output = [99_u64];
    assert!(matches!(
        myers_distance_batch(&LITERAL_MASKS, 0, &reference, &[0], &[1], &mut output),
        Err(BatchError::QueryLength { observed: 0 })
    ));
    assert_eq!(output, [99]);
    assert!(matches!(
        myers_distance_batch(&LITERAL_MASKS, 4, &reference, &[3], &[2], &mut output),
        Err(BatchError::CandidateOutOfBounds { candidate: 0, .. })
    ));
    assert_eq!(output, [99]);
}

#[test]
fn runtime_and_available_simd_equal_scalar_for_variable_candidates() {
    let reference = (0_u16..193)
        .map(|index| u8::try_from(index * 3 % 5).expect("modulo five fits u8"))
        .collect::<Vec<_>>();
    let starts = [0, 3, 9, 17, 35, 67, 101, 140, 192];
    let lengths = [7, 7, 7, 16, 31, 48, 64, 53, 1];
    let mut scalar = [0_u64; 9];
    myers_distance_batch_with_flavor(
        KernelFlavor::Scalar,
        &LITERAL_MASKS,
        64,
        &reference,
        &starts,
        &lengths,
        &mut scalar,
    )
    .unwrap();
    let mut runtime = [0_u64; 9];
    myers_distance_batch(
        &LITERAL_MASKS,
        64,
        &reference,
        &starts,
        &lengths,
        &mut runtime,
    )
    .unwrap();
    assert_eq!(runtime, scalar);
    for flavor in [
        KernelFlavor::Sse42,
        KernelFlavor::Avx2,
        KernelFlavor::Avx512,
    ] {
        if flavor_available(flavor) {
            let mut observed = [0_u64; 9];
            myers_distance_batch_with_flavor(
                flavor,
                &LITERAL_MASKS,
                64,
                &reference,
                &starts,
                &lengths,
                &mut observed,
            )
            .unwrap();
            assert_eq!(observed, scalar, "{flavor}");
        }
    }
}

#[test]
fn empty_reference_and_unknown_codes_have_global_semantics() {
    let reference = [9_u8; 8];
    let mut output = [0_u64; 2];
    let mut masks = LITERAL_MASKS;
    masks[4] = u64::MAX;
    myers_distance_batch(&masks, 4, &reference, &[0, 0], &[0, 8], &mut output).unwrap();
    assert_eq!(output, [4, 8]);
}

#[test]
fn narrow_available_simd_equals_scalar_for_full_and_partial_chunks() {
    let masks = [1_u8, 2, 4, 8, 0];
    let query = (0..100)
        .map(|position| u8::try_from(position % 4).unwrap())
        .collect::<Vec<_>>();
    let distance = 2;
    let pattern_length = query.len() + 2 * distance;
    let mut patterns = Vec::new();
    for lane in 0..13 {
        let mut pattern = vec![u8::try_from(lane % 4).unwrap(); pattern_length];
        let shift = lane % 5;
        let copied = query.len().min(pattern_length - shift);
        pattern[shift..shift + copied].copy_from_slice(&query[..copied]);
        if lane % 3 == 1 {
            pattern[37] = (pattern[37] + 1) % 4;
        }
        patterns.extend(pattern);
    }
    let mut scalar = vec![NarrowBandedResult::ABSENT; 13];
    narrow_banded_prefix_batch_with_flavor(
        NarrowBandedFlavor::Scalar,
        &masks,
        &query,
        &patterns,
        distance,
        &mut scalar,
    )
    .unwrap();
    let mut runtime = vec![NarrowBandedResult::ABSENT; 13];
    narrow_banded_prefix_batch(&masks, &query, &patterns, distance, &mut runtime).unwrap();
    assert_eq!(runtime, scalar);
    if narrow_sse42_available() {
        let mut sse42 = vec![NarrowBandedResult::ABSENT; 13];
        narrow_banded_prefix_batch_with_flavor(
            NarrowBandedFlavor::Sse42,
            &masks,
            &query,
            &patterns,
            distance,
            &mut sse42,
        )
        .unwrap();
        assert_eq!(sse42, scalar);
    }
    if narrow_avx2_available() {
        let mut avx2 = vec![NarrowBandedResult::ABSENT; 13];
        narrow_banded_prefix_batch_with_flavor(
            NarrowBandedFlavor::Avx2,
            &masks,
            &query,
            &patterns,
            distance,
            &mut avx2,
        )
        .unwrap();
        assert_eq!(avx2, scalar);
    }
}

#[test]
fn fixed_start_batch_equals_center_row_for_full_and_partial_chunks() {
    let masks = [1_u8, 2, 4, 2 | 8, 0];
    let mut state = 0x3c6e_f372_fe94_f82b_u64;
    for max_distance in 0..=4 {
        let query_length = 23 + max_distance;
        let pattern_length = query_length + max_distance * 2;
        let mut query = vec![0_u8; query_length];
        for code in &mut query {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1);
            *code = u8::try_from((state >> 61) % 5).unwrap();
        }
        let mut patterns = Vec::new();
        let mut expected = Vec::new();
        for _ in 0..37 {
            let mut pattern = vec![0_u8; pattern_length];
            for code in &mut pattern {
                state = state
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1);
                *code = u8::try_from((state >> 61) % 5).unwrap();
            }
            expected.push(narrow_fixed_start_scalar_one(
                &masks,
                &query,
                &pattern,
                max_distance,
            ));
            patterns.extend(pattern);
        }
        let mut observed = vec![NarrowEndpointDistances::EMPTY; expected.len()];
        narrow_banded_fixed_start_batch(&masks, &query, &patterns, max_distance, &mut observed)
            .unwrap();
        assert_eq!(observed, expected, "distance {max_distance}");
        #[cfg(target_arch = "x86_64")]
        if narrow_sse42_available() {
            let mut sse42 = vec![NarrowEndpointDistances::EMPTY; expected.len()];
            // SAFETY: the runtime check proves SSE4.2 support, and the local
            // construction supplies fixed-width patterns and output slots.
            unsafe {
                narrow_fixed_start_sse42(
                    &masks,
                    &query,
                    &patterns,
                    pattern_length,
                    max_distance,
                    &mut sse42,
                );
            }
            assert_eq!(sse42, expected, "SSE4.2 distance {max_distance}");
        }
    }
}

#[test]
fn fixed_start_gather_equals_materialized_patterns() {
    let masks = [1_u8, 2, 4, 8, 0];
    let query = (0..151)
        .map(|position| u8::try_from(position % 4).unwrap())
        .collect::<Vec<_>>();
    let reference = (0..4096)
        .map(|position| u8::try_from((position * 3 + position / 7) % 5).unwrap())
        .collect::<Vec<_>>();
    let starts = (0..37).map(|lane| lane * 83).collect::<Vec<_>>();
    let distance = 2;
    let pattern_length = query.len() + 2 * distance;
    let mut patterns = Vec::new();
    for &start in &starts {
        let mut pattern = vec![4_u8; pattern_length];
        pattern[distance..].copy_from_slice(&reference[start..start + query.len() + distance]);
        patterns.extend(pattern);
    }
    let mut materialized = vec![NarrowEndpointDistances::EMPTY; starts.len()];
    narrow_banded_fixed_start_batch(&masks, &query, &patterns, distance, &mut materialized)
        .unwrap();
    let mut gathered = vec![NarrowEndpointDistances::EMPTY; starts.len()];
    narrow_banded_fixed_start_gather_batch(
        &masks,
        &query,
        &reference,
        &starts,
        distance,
        &mut gathered,
    )
    .unwrap();
    assert_eq!(gathered, materialized);
}

#[test]
fn u128_prefix_batch_equals_scalar_for_full_and_partial_chunks() {
    let query_length = 100;
    let max_distance = 2;
    let pattern_length = query_length + max_distance;
    let equality_masks = [
        (0..query_length)
            .filter(|position| position % 4 == 0)
            .fold(0_u128, |mask, position| mask | (1_u128 << position)),
        (0..query_length)
            .filter(|position| position % 4 == 1)
            .fold(0_u128, |mask, position| mask | (1_u128 << position)),
        (0..query_length)
            .filter(|position| position % 4 == 2)
            .fold(0_u128, |mask, position| mask | (1_u128 << position)),
        (0..query_length)
            .filter(|position| position % 4 == 3)
            .fold(0_u128, |mask, position| mask | (1_u128 << position)),
        0,
    ];
    let mut state = 0xa54f_f53a_5f1d_36f1_u64;
    let mut patterns = Vec::new();
    let mut expected = Vec::new();
    for _ in 0..13 {
        let mut pattern = vec![0_u8; pattern_length];
        for code in &mut pattern {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1);
            *code = u8::try_from((state >> 61) % 5).unwrap();
        }
        expected.push(myers_prefix_distances_u128_scalar_one(
            &equality_masks,
            query_length,
            &pattern,
            max_distance,
        ));
        patterns.extend(pattern);
    }
    let mut observed = vec![NarrowEndpointDistances::EMPTY; expected.len()];
    myers_prefix_distances_u128_batch(
        &equality_masks,
        query_length,
        &patterns,
        max_distance,
        &mut observed,
    )
    .unwrap();
    assert_eq!(observed, expected);
}

#[test]
fn narrow_placement_available_simd_equals_scalar_for_every_supported_band() {
    let masks = [1_u8, 2, 4, 2 | 8, 0];
    let mut state = 0xbb67_ae85_84ca_a73b_u64;
    for max_distance in 0..=MAX_NARROW_BAND_DISTANCE {
        let query_length = 17 + max_distance;
        let pattern_length = query_length + max_distance * 2;
        for _ in 0..64 {
            let mut query = vec![0_u8; query_length];
            let mut pattern = vec![0_u8; pattern_length];
            for code in query.iter_mut().chain(&mut pattern) {
                state = state
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);
                *code = u8::try_from((state >> 61) & 3).unwrap();
            }
            let scalar = narrow_placement_distances_scalar(&masks, &query, &pattern, max_distance);
            let runtime = narrow_banded_placement_distances(&masks, &query, &pattern, max_distance)
                .expect("validated placement dimensions");
            assert_eq!(runtime, scalar);
            #[cfg(target_arch = "x86_64")]
            if narrow_sse42_available() {
                // SAFETY: the runtime check proves SSE4.2 support and all
                // dimensions satisfy the private kernel preconditions.
                let sse42 = unsafe {
                    narrow_placement_distances_sse42(&masks, &query, &pattern, max_distance)
                };
                assert_eq!(sse42, scalar);
            }
            #[cfg(target_arch = "x86_64")]
            if narrow_avx2_available() {
                // SAFETY: the runtime check proves AVX2 support and all
                // dimensions satisfy the private kernel preconditions.
                let avx2 = unsafe {
                    narrow_placement_distances_avx2(&masks, &query, &pattern, max_distance)
                };
                assert_eq!(avx2, scalar);
            }
        }
    }
}

#[test]
fn narrow_placement_distance_three_compact_state_equals_general_kernel() {
    let masks = [1_u8, 2, 4, 2 | 8, 0];
    let mut state = 0x510e_527f_ade6_82d1_u64;
    for query_length in [1, 17, 100, 192] {
        for _ in 0..64 {
            let mut query = vec![0_u8; query_length];
            let mut pattern = vec![0_u8; query_length + 6];
            for code in query.iter_mut().chain(&mut pattern) {
                state = state
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);
                *code = u8::try_from((state >> 61) % 5).unwrap();
            }
            let expected = narrow_banded_placement_distances(&masks, &query, &pattern, 3).unwrap();
            let observed = narrow_banded_placement_distances_d3(&masks, &query, &pattern).unwrap();
            assert_eq!(observed, expected);
        }
    }
}

#[test]
fn narrow_placement_distance_five_compact_state_equals_general_kernel() {
    let masks = [1_u8, 2, 4, 2 | 8, 0];
    let mut state = 0xa54f_f53a_5f1d_36f1_u64;
    for query_length in [1, 17, 100, 192] {
        for _ in 0..64 {
            let mut query = vec![0_u8; query_length];
            let mut pattern = vec![0_u8; query_length + 10];
            for code in query.iter_mut().chain(&mut pattern) {
                state = state
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);
                *code = u8::try_from((state >> 61) % 5).unwrap();
            }
            let expected = narrow_banded_placement_distances(&masks, &query, &pattern, 5).unwrap();
            let observed = narrow_banded_placement_distances_d5(&masks, &query, &pattern).unwrap();
            assert_eq!(observed, expected);
        }
    }
}

#[test]
fn narrow_placement_batch_equals_independent_candidates() {
    let masks = [1_u8, 2, 4, 2 | 8, 0];
    let mut state = 0x3c6e_f372_fe94_f82b_u64;
    for max_distance in 0..=4 {
        let query_length = 19 + max_distance;
        let pattern_length = query_length + max_distance * 2;
        let mut query = vec![0_u8; query_length];
        for code in &mut query {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            *code = u8::try_from((state >> 61) & 3).unwrap();
        }
        let maximum_candidates = (32 / (2 * max_distance + 1)).min(4);
        for candidates in 1..=maximum_candidates {
            let mut patterns = vec![0_u8; pattern_length * candidates];
            for code in &mut patterns {
                state = state
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);
                *code = u8::try_from((state >> 61) & 3).unwrap();
            }
            let expected = patterns
                .chunks_exact(pattern_length)
                .map(|pattern| {
                    narrow_banded_placement_distances(&masks, &query, pattern, max_distance)
                        .unwrap()
                })
                .collect::<Vec<_>>();
            let mut observed = vec![NarrowPlacementDistances::EMPTY; candidates];
            narrow_banded_placement_distances_batch(
                &masks,
                &query,
                &patterns,
                max_distance,
                &mut observed,
            )
            .unwrap();
            assert_eq!(observed, expected);
        }
    }
}

#[test]
fn narrow_placement_distance_three_compact_batch_equals_general_batch() {
    let masks = [1_u8, 2, 4, 2 | 8, 0];
    let mut state = 0x1f83_d9ab_fb41_bd6b_u64;
    for query_length in [1, 17, 100, 192] {
        let pattern_length = query_length + 6;
        for candidates in 1..=4 {
            let mut query = vec![0_u8; query_length];
            let mut patterns = vec![0_u8; pattern_length * candidates];
            for code in query.iter_mut().chain(&mut patterns) {
                state = state
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);
                *code = u8::try_from((state >> 61) % 5).unwrap();
            }
            let mut expected = vec![NarrowPlacementDistances::EMPTY; candidates];
            narrow_banded_placement_distances_batch(&masks, &query, &patterns, 3, &mut expected)
                .unwrap();
            let mut observed = vec![NarrowPlacementDistances::EMPTY; candidates];
            narrow_banded_placement_distances_batch_d3(&masks, &query, &patterns, &mut observed)
                .unwrap();
            assert_eq!(observed, expected);
        }
    }
}

#[test]
fn narrow_placement_distance_five_compact_batch_equals_independent_candidates() {
    let masks = [1_u8, 2, 4, 2 | 8, 0];
    let mut state = 0x5be0_cd19_137e_2179_u64;
    for query_length in [1, 17, 100, 192] {
        let pattern_length = query_length + 10;
        let maximum = 2;
        for candidates in 1..=maximum {
            let mut query = vec![0_u8; query_length];
            let mut patterns = vec![0_u8; pattern_length * candidates];
            for code in query.iter_mut().chain(&mut patterns) {
                state = state
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);
                *code = u8::try_from((state >> 61) % 5).unwrap();
            }
            let mut expected = vec![NarrowPlacementDistances::EMPTY; candidates];
            for (pattern, destination) in patterns.chunks_exact(pattern_length).zip(&mut expected) {
                *destination =
                    narrow_banded_placement_distances(&masks, &query, pattern, 5).unwrap();
            }
            let mut observed = vec![NarrowPlacementDistances::EMPTY; candidates];
            narrow_banded_placement_distances_batch_d5(&masks, &query, &patterns, &mut observed)
                .unwrap();
            assert_eq!(observed, expected);
        }
    }
}

#[test]
fn narrow_scalar_matches_independent_banded_semiglobal_dp() {
    let masks = [1_u8, 2, 4, 2 | 8, 0];
    let mut state = 0x6a09_e667_f3bc_c909_u64;
    for max_distance in 0..=3 {
        let query_length = 7 + max_distance;
        let pattern_length = query_length + 2 * max_distance;
        for case in 0..256 {
            let mut query = vec![0_u8; query_length];
            let mut pattern = vec![0_u8; pattern_length];
            for code in query.iter_mut().chain(&mut pattern) {
                state = state
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);
                *code = u8::try_from((state >> 61) & 3).unwrap();
            }
            let observed = narrow_banded_scalar_one(&masks, &query, &pattern, max_distance);
            let expected = brute_narrow_result(masks, &query, &pattern, max_distance);
            assert_eq!(
                observed, expected,
                "case {case}, k={max_distance}, query={query:?}, pattern={pattern:?}"
            );
        }
    }
}

#[test]
fn narrow_center_tie_and_bisulfite_relation_are_explicit() {
    let query = vec![3_u8; 80];
    let pattern = vec![1_u8; 84];
    let mut output = [NarrowBandedResult::ABSENT];
    narrow_banded_prefix_batch(&[1, 2, 4, 2 | 8, 0], &query, &pattern, 2, &mut output).unwrap();
    assert_eq!(
        output[0],
        NarrowBandedResult {
            distance: 0,
            prefix_length: 82,
            tied_prefix_mask: 0b1_1111,
            endpoint_distances: {
                let mut distances = NarrowEndpointDistances::EMPTY;
                for delta in 0..5 {
                    distances.insert(delta, 0);
                }
                distances
            },
        }
    );
}

#[test]
fn narrow_frontier_preserves_higher_in_budget_endpoint_distances() {
    let result = narrow_banded_scalar_one(&[1, 2, 4, 8, 0], &[0], &[0, 1, 1], 1);
    assert_eq!(result.distance, 0);
    assert_eq!(result.prefix_length, 1);
    assert_eq!(result.tied_prefix_mask, 0b001);
    assert_eq!(result.endpoint_distances.in_budget_mask(), 0b111);
    assert_eq!(result.endpoint_distances.distance(0), Some(0));
    assert_eq!(result.endpoint_distances.distance(1), Some(1));
    assert_eq!(result.endpoint_distances.distance(2), Some(1));
    assert_eq!(result.endpoint_distances.mask_at_distance(0), 0b001);
    assert_eq!(result.endpoint_distances.mask_at_distance(1), 0b110);
}

#[test]
fn unanchored_paths_are_outside_the_fixed_narrow_contract() {
    let query = [0_u8, 1, 2, 1];
    let pattern = [0_u8, 2, 1, 0, 0, 3, 3, 2];
    assert_eq!(brute_unrestricted_semiglobal_distance(&query, &pattern), 2);
    assert_eq!(
        narrow_banded_scalar_one(&[1, 2, 4, 8, 0], &query, &pattern, 2),
        NarrowBandedResult::ABSENT
    );
}
