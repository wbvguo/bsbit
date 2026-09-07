//! White-box tests for canonical paired-end mapping.
//!
//! Kept outside implementation `src/` while remaining a child module so private
//! invariants can be tested without widening the crate API.

use super::*;
use crate::alignment_policy::{
    DEFAULT_MAX_SOFT_CLIP_BASES, LOCAL_FILTER_BLOCKS, ORIGIN_ENDPOINT_ADAPTER_CLIP_OPEN_PENALTY,
    ORIGIN_ENDPOINT_CLIP_EXTENSION_PENALTY, ORIGIN_ENDPOINT_CLIP_OPEN_PENALTY,
    SEMI_GLOBAL_ADMISSION_EDIT_PENALTY, SEMI_GLOBAL_CLIP_PENALTY, SEMI_GLOBAL_EDIT_PENALTY,
    SEMI_GLOBAL_MIN_ALIGNED_BASES, SENSITIVE_ADAPTIVE_BOUNDARY_SHIFTS,
    SENSITIVE_ADAPTIVE_MIN_BLOCK_BASES, SENSITIVE_BALANCED_BOUNDARY_SHIFTS, SENSITIVE_CLIP_PENALTY,
    SENSITIVE_PROOF_BLOCKS,
};
use crate::paired_end::rescue::append_local_flexible_proof_candidates;
use crate::placement::placement_origin_key;
use crate::search::combined_adaptive::FLEXIBLE_NOMINAL_PROOF;
use bsbit_core::bisulfite::BisulfiteStrand;
use bsbit_core::sequence::NormalizedSequence;
use bsbit_index::reference::{ContigInput, ReferenceBuildLimits};

fn reference(bases: &[Base]) -> ReferenceIndex {
    ReferenceIndex::build(
        vec![ContigInput::new(
            b"test".to_vec(),
            NormalizedSequence::from_bases(bases.iter().copied()),
        )],
        ReferenceBuildLimits::MAX,
    )
    .expect("test reference builds")
}

fn repeated_named_reference(names: [&[u8]; 2]) -> ReferenceIndex {
    ReferenceIndex::build(
        names
            .into_iter()
            .map(|name| {
                ContigInput::new(
                    name.to_vec(),
                    NormalizedSequence::from_bases([Base::A; 128]),
                )
            })
            .collect(),
        ReferenceBuildLimits::MAX,
    )
    .expect("repeated named reference builds")
}

#[test]
fn paired_reporting_lottery_is_seeded_and_reference_order_independent() {
    let forward = repeated_named_reference([b"alpha", b"beta"]);
    let reversed = repeated_named_reference([b"beta", b"alpha"]);
    let reads = [&[Base::A; 20][..], &[Base::A; 20][..]];
    let pair = |contig_ordinal| PairedPlacement {
        mate1: ReadPlacement::strict(contig_ordinal, 10, 30, BisulfiteStrand::OT, 0),
        mate2: ReadPlacement::strict(contig_ordinal, 80, 100, BisulfiteStrand::OB, 0),
        template_start: 10,
        template_end: 100,
        distance: 0,
        score: 0,
    };
    let choose_name = |reference: &ReferenceIndex, mut pairs: Vec<PairedPlacement>, seed| {
        let retained = pairs.clone();
        prefer_fair_pair_representative(
            reference,
            reads,
            &mut pairs,
            ReportingTieBreak { seed, read_key: 91 },
            false,
        )
        .expect("paired fair tie-break succeeds");
        assert_eq!(pairs.len(), retained.len());
        assert!(pairs.iter().all(|candidate| retained.contains(candidate)));
        reference
            .contig_by_ordinal(pairs[0].mate1().contig_ordinal())
            .expect("selected contig exists")
            .name()
            .to_vec()
    };

    for seed in 0..32 {
        assert_eq!(
            choose_name(&forward, vec![pair(0), pair(1)], seed),
            choose_name(&reversed, vec![pair(1), pair(0)], seed),
        );
    }
    let first = choose_name(&forward, vec![pair(0), pair(1)], 0);
    let alternate = (1..256)
        .map(|seed| choose_name(&forward, vec![pair(0), pair(1)], seed))
        .find(|name| *name != first)
        .expect("a seed changes a fair paired two-way lottery");
    assert_ne!(first, alternate);
}

#[test]
fn nondirectional_second_pass_uses_final_mate_order_for_the_lottery() {
    let reference = repeated_named_reference([b"alpha", b"beta"]);
    let read1 = [Base::A; 20];
    let read2 = [Base::A; 30];
    // A non-directional second pass is searched internally as [R2, R1].
    let internal_reads = [&read2[..], &read1[..]];
    let candidate = |mate1_contig, mate2_contig| PairedPlacement {
        mate1: ReadPlacement::strict(mate1_contig, 70, 100, BisulfiteStrand::CTOT, 0),
        mate2: ReadPlacement::strict(mate2_contig, 10, 30, BisulfiteStrand::OT, 0),
        template_start: 10,
        template_end: 100,
        distance: 0,
        score: 0,
    };
    let left = candidate(0, 1);
    let right = candidate(1, 0);

    let seed = (0..10_000)
        .find(|&seed| {
            let tie_break = ReportingTieBreak { seed, read_key: 17 };
            let final_left = pair_origin_hash_in_reporting_order(
                &reference,
                tie_break,
                left,
                internal_reads,
                true,
            )
            .expect("final-order hash succeeds");
            let final_right = pair_origin_hash_in_reporting_order(
                &reference,
                tie_break,
                right,
                internal_reads,
                true,
            )
            .expect("final-order hash succeeds");
            let internal_left = pair_origin_hash_in_reporting_order(
                &reference,
                tie_break,
                left,
                internal_reads,
                false,
            )
            .expect("internal-order hash succeeds");
            let internal_right = pair_origin_hash_in_reporting_order(
                &reference,
                tie_break,
                right,
                internal_reads,
                false,
            )
            .expect("internal-order hash succeeds");
            (final_left < final_right) != (internal_left < internal_right)
        })
        .expect("a seed distinguishes final and internal mate order");
    let tie_break = ReportingTieBreak { seed, read_key: 17 };
    let expected =
        if pair_origin_hash_in_reporting_order(&reference, tie_break, left, internal_reads, true)
            .expect("final-order hash succeeds")
            < pair_origin_hash_in_reporting_order(
                &reference,
                tie_break,
                right,
                internal_reads,
                true,
            )
            .expect("final-order hash succeeds")
        {
            left
        } else {
            right
        };
    let mut pairs = vec![left, right];
    prefer_fair_pair_representative(&reference, internal_reads, &mut pairs, tie_break, true)
        .expect("swapped-order fair tie-break succeeds");
    assert_eq!(pairs[0], expected);
}

#[test]
fn sensitive_profile_is_separate_and_prefers_whole_read_edits() {
    let default = PairedSearchMode::Default.limits();
    let sensitive = PairedSearchMode::Sensitive.limits();
    assert!(sensitive.maximum_seed_hits > default.maximum_seed_hits);
    assert_eq!(
        sensitive.maximum_combined_rescue_hits,
        default.maximum_combined_rescue_hits
    );
    assert_eq!(sensitive.maximum_seed_rounds, default.maximum_seed_rounds);
    assert!(
        PairedSearchMode::Sensitive.semi_global_clip_penalty()
            > PairedSearchMode::Default.semi_global_clip_penalty()
    );
    assert!(PairedSearchMode::Sensitive.is_sensitive());
}

#[test]
fn mapping_options_fix_primary_and_adapter_trimmed_policies() {
    let default = PairedAlignmentOptions::primary(
        LibraryProfile::Directional,
        PairedSearchMode::Default,
        0,
        1_000,
    );
    let sensitive = PairedAlignmentOptions::primary(
        LibraryProfile::Directional,
        PairedSearchMode::Sensitive,
        0,
        1_000,
    );
    let trimmed = PairedAlignmentOptions::adapter_trimmed(
        LibraryProfile::Directional,
        PairedSearchMode::Sensitive,
        0,
        1_000,
    );
    assert_eq!(default.derived_policy(), (MAX_EDIT_DISTANCE, false, false));
    assert_eq!(sensitive.derived_policy(), (MAX_EDIT_DISTANCE, true, true));
    assert_eq!(trimmed.derived_policy(), (MAX_EDIT_DISTANCE, true, false));
    let bounded = sensitive
        .with_maximum_edit_distance(3)
        .expect("PE edit budget inside the fixed domain");
    assert_eq!(bounded.maximum_edit_distance(), 3);
    assert_eq!(bounded.derived_policy(), (3, true, true));
    assert!(matches!(
        sensitive.with_maximum_edit_distance(MAX_EDIT_DISTANCE + 1),
        Err(AlignmentError::UnsupportedEditDistance { .. })
    ));
}

#[test]
fn adaptive_ranked_partitions_are_disjoint_and_cover_the_read() {
    for block_count in 2..=SENSITIVE_PROOF_BLOCKS {
        let balanced = ranked_block_boundaries(
            150,
            block_count,
            SENSITIVE_BALANCED_BOUNDARY_SHIFTS,
            SENSITIVE_ADAPTIVE_MIN_BLOCK_BASES,
        )
        .expect("qualified adaptive partition");
        assert_eq!(balanced[0], 0);
        assert_eq!(balanced[block_count], 150);

        let partition_count = SENSITIVE_ADAPTIVE_BOUNDARY_SHIFTS
            .len()
            .pow(u32::try_from(block_count - 1).expect("bounded exponent"));
        for encoded in 0..partition_count {
            let mut remainder = encoded;
            let mut boundary_shifts = [0_i8; SENSITIVE_PROOF_BLOCKS - 1];
            for shift in &mut boundary_shifts[..block_count - 1] {
                *shift = SENSITIVE_ADAPTIVE_BOUNDARY_SHIFTS
                    [remainder % SENSITIVE_ADAPTIVE_BOUNDARY_SHIFTS.len()];
                remainder /= SENSITIVE_ADAPTIVE_BOUNDARY_SHIFTS.len();
            }
            let boundaries = ranked_block_boundaries(
                150,
                block_count,
                boundary_shifts,
                SENSITIVE_ADAPTIVE_MIN_BLOCK_BASES,
            )
            .expect("qualified adaptive partition");
            assert_eq!(boundaries[0], 0);
            assert_eq!(boundaries[block_count], 150);
            assert!(
                boundaries[..=block_count]
                    .windows(2)
                    .all(|window| { window[1] - window[0] >= SENSITIVE_ADAPTIVE_MIN_BLOCK_BASES })
            );
        }
    }
}

#[test]
fn sensitive_frontier_completion_covers_every_provisional_unique_result() {
    assert!(sensitive_unique_frontier_completion_required(
        PairMappingStatus::Unique,
    ));
    assert!(!sensitive_unique_frontier_completion_required(
        PairMappingStatus::Ambiguous,
    ));
    assert!(!sensitive_unique_frontier_completion_required(
        PairMappingStatus::Unmapped,
    ));
}

#[test]
fn completed_worse_frontier_is_retained_as_adverse_mapq_evidence() {
    let incumbent = PairAlignmentMetrics {
        compatible_pairs: 1,
        best_pair_placements: 1,
        best_pair_score: Some(-12),
        mapq_compatible_pairs: 1,
        mapq_best_pair_score: Some(-12),
        frontier_complete: true,
        alternative_margin_frontier_complete: true,
        ..empty_pair_metrics()
    };
    let completed = PairAlignmentMetrics {
        compatible_pairs: 1,
        best_pair_placements: 1,
        best_pair_score: Some(-16),
        mapq_compatible_pairs: 1,
        mapq_best_pair_score: Some(-16),
        frontier_complete: true,
        alternative_margin_frontier_complete: true,
        ..empty_pair_metrics()
    };

    let merged = retain_completed_runner_up_evidence(incumbent, completed);
    assert_eq!(merged.best_pair_score, Some(-12));
    assert_eq!(merged.second_best_pair_score, Some(-16));
    assert_eq!(merged.mapq_second_best_pair_score, Some(-16));
    assert_eq!(merged.near_best_pairings, 1);
    assert_eq!(merged.mapq_near_best_pairings, 1);
    assert_eq!(
        sensitive_effective_mapping_quality(PairMappingStatus::Unique, merged),
        39,
    );
}

#[test]
fn targeted_semi_global_uses_frontier_state_instead_of_benchmark_cells() {
    let complete_ambiguous = PairAlignmentMetrics {
        best_pair_placements: 2,
        best_pair_score: Some(100),
        frontier_complete: true,
        alternative_margin_frontier_complete: true,
        ..empty_pair_metrics()
    };
    assert!(sensitive_targeted_semi_global_required(
        PairMappingStatus::Ambiguous,
        complete_ambiguous,
    ));

    let incomplete = PairAlignmentMetrics {
        frontier_complete: false,
        alternative_margin_frontier_complete: false,
        ..complete_ambiguous
    };
    assert!(sensitive_targeted_semi_global_required(
        PairMappingStatus::Ambiguous,
        incomplete,
    ));
    assert!(!sensitive_targeted_semi_global_required(
        PairMappingStatus::Unmapped,
        complete_ambiguous,
    ));

    let high_confidence = PairAlignmentMetrics {
        best_pair_placements: 1,
        ..complete_ambiguous
    };
    assert!(!sensitive_targeted_semi_global_required(
        PairMappingStatus::Unique,
        high_confidence,
    ));
}

#[test]
fn rejected_targeted_semi_global_restores_every_incumbent_origin() {
    let first = PairedPlacement {
        mate1: ReadPlacement::strict(0, 100, 250, BisulfiteStrand::OT, 1),
        mate2: ReadPlacement::strict(0, 300, 450, BisulfiteStrand::CTOT, 0),
        template_start: 100,
        template_end: 450,
        distance: 1,
        score: 1,
    };
    let second = PairedPlacement {
        mate1: ReadPlacement::strict(1, 500, 650, BisulfiteStrand::OT, 1),
        mate2: ReadPlacement::strict(1, 700, 850, BisulfiteStrand::CTOT, 0),
        template_start: 500,
        template_end: 850,
        ..first
    };
    let incumbent = vec![first, second];
    let mut speculative = vec![second];

    restore_rejected_targeted_frontier(&mut speculative, incumbent.clone(), true);

    assert_eq!(speculative, incumbent);

    let mut unique = vec![second];
    restore_rejected_targeted_frontier(&mut unique, incumbent, false);
    assert_eq!(unique, [first]);
}

#[test]
fn incomplete_sensitive_frontier_cannot_claim_unique() {
    let mut result = (
        PairMappingStatus::Unique,
        PairAlignmentMetrics {
            best_pair_placements: 1,
            ..empty_pair_metrics()
        },
        Some(4),
    );
    conservatively_mark_incomplete_frontier(&mut result, false);
    assert_eq!(result.0, PairMappingStatus::Ambiguous);
    assert_eq!(result.1.best_pair_placements, 2);
    assert_eq!(result.2, None);

    let mut complete = (
        PairMappingStatus::Unique,
        PairAlignmentMetrics {
            best_pair_placements: 1,
            ..empty_pair_metrics()
        },
        Some(4),
    );
    conservatively_mark_incomplete_frontier(&mut complete, true);
    assert_eq!(complete.0, PairMappingStatus::Unique);
    assert_eq!(complete.1.best_pair_placements, 1);
    assert_eq!(complete.2, Some(4));
}

#[test]
fn non_directional_result_merge_swaps_mates_and_resolves_global_evidence() {
    let directional_pair = PairedPlacement {
        mate1: placement(0, 100, 151, BisulfiteStrand::OT, 0),
        mate2: placement(0, 250, 301, BisulfiteStrand::CTOT, 0),
        template_start: 100,
        template_end: 301,
        distance: 0,
        score: 0,
    };
    let swapped_pass_pair = PairedPlacement {
        mate1: placement(0, 80, 131, BisulfiteStrand::OT, 1),
        mate2: placement(0, 280, 331, BisulfiteStrand::CTOT, 0),
        template_start: 80,
        template_end: 331,
        distance: 1,
        score: 1,
    };
    let result = |pair, score, mate1_rows, mate2_rows| {
        let metrics = PairAlignmentMetrics {
            mate1: ReadAlignmentMetrics {
                located_rows: mate1_rows,
                ..ReadAlignmentMetrics::default()
            },
            mate2: ReadAlignmentMetrics {
                located_rows: mate2_rows,
                ..ReadAlignmentMetrics::default()
            },
            compatible_pairs: 1,
            best_pair_placements: 1,
            best_pair_score: Some(score),
            frontier_complete: true,
            alternative_margin_frontier_complete: true,
            ..empty_pair_metrics()
        };
        let best_pair = Some(pair);
        PairedBatchResult {
            class: PairMappingStatus::Unique,
            metrics,
            best_pair,
            second_best_distance: None,
        }
    };

    let complementary = swap_batch_result_mates(result(swapped_pass_pair, -4, 7, 11));
    let complementary_pair = complementary.best_pair().expect("swapped pair retained");
    assert_eq!(complementary_pair.mate1().strand(), BisulfiteStrand::CTOT);
    assert_eq!(complementary_pair.mate2().strand(), BisulfiteStrand::OT);
    assert_eq!(complementary.metrics().mate1.located_rows, 11);
    assert_eq!(complementary.metrics().mate2.located_rows, 7);

    let merged =
        merge_non_directional_batch_results(&result(directional_pair, -8, 3, 5), &complementary);
    assert_eq!(merged.class(), PairMappingStatus::Unique);
    assert_eq!(
        merged.best_pair().expect("global winner").mate1().strand(),
        BisulfiteStrand::CTOT
    );
    assert_eq!(merged.best_pair_score(), Some(-4));
    assert_eq!(merged.second_best_pair_score(), Some(-8));
    assert_eq!(merged.near_best_pairings(), 1);
    assert_eq!(merged.metrics().mate1.located_rows, 14);
    assert_eq!(merged.metrics().mate2.located_rows, 12);
    let tied =
        merge_non_directional_batch_results(&result(directional_pair, -4, 1, 1), &complementary);
    assert_eq!(tied.class(), PairMappingStatus::Ambiguous);
    assert_eq!(tied.best_pair_score(), Some(-4));
    assert_eq!(tied.second_best_pair_score(), Some(-4));
    assert!(tied.metrics().best_pair_placements >= 2);
    assert_eq!(super::super::mapq::evidence_mapping_quality(tied), 0);
    assert_eq!(
        tied.best_pair()
            .expect("original-strand tie preference")
            .mate1()
            .strand(),
        BisulfiteStrand::OT
    );
}

#[test]
fn pre_sorted_candidate_verification_matches_sorting_entrypoint() {
    let bases = (0..320)
        .map(|position| [Base::A, Base::C, Base::G, Base::T][position % 4])
        .collect::<Vec<_>>();
    let index = reference(&bases);
    let read = &bases[40..191];
    let candidates = vec![
        ReadCandidate {
            contig_ordinal: 0,
            start: 44,
            strand: BisulfiteStrand::OT,
            proof_mask: FLEXIBLE_NOMINAL_PROOF,
        },
        ReadCandidate {
            contig_ordinal: 0,
            start: 40,
            strand: BisulfiteStrand::OT,
            proof_mask: FLEXIBLE_NOMINAL_PROOF,
        },
    ];

    let mut ordinary = ReadWorkspace::with_capacity(8, 8);
    ordinary.candidate_nominals.clone_from(&candidates);
    let (ordinary_placements, ordinary_metrics) = ordinary
        .verify_candidates_with_budget(
            &index,
            read,
            ReadAlignmentMetrics::default(),
            INITIAL_EDIT_DISTANCE,
        )
        .expect("ordinary verification");
    let ordinary_placements = ordinary_placements.to_vec();

    let mut sorted = ReadWorkspace::with_capacity(8, 8);
    sorted.candidate_nominals = candidates;
    sort_nominal_candidates(&mut sorted.candidate_nominals);
    let (sorted_placements, sorted_metrics) = sorted
        .verify_sorted_candidates_with_budget(
            &index,
            read,
            ReadAlignmentMetrics::default(),
            INITIAL_EDIT_DISTANCE,
        )
        .expect("pre-sorted verification");
    assert_eq!(sorted_placements, ordinary_placements);
    assert_eq!(
        (
            sorted_metrics.located_rows,
            sorted_metrics.emitted_candidate_starts,
            sorted_metrics.distinct_candidate_starts,
            sorted_metrics.verified_placements,
        ),
        (
            ordinary_metrics.located_rows,
            ordinary_metrics.emitted_candidate_starts,
            ordinary_metrics.distinct_candidate_starts,
            ordinary_metrics.verified_placements,
        )
    );
}

#[test]
fn local_flexible_rescue_keeps_a_d5_candidate_from_six_disjoint_blocks() {
    let read = (0..151)
        .map(|position| [Base::A, Base::C, Base::G, Base::T][position % 4])
        .collect::<Vec<_>>();
    let mut reference_bases = read.clone();
    // Six d5 proof blocks start at 0, 26, 51, 76, 101, and 126.
    // Disturb the first five and leave the sixth as the exact proof.
    for position in [5_usize, 31, 56, 81, 106] {
        reference_bases[position] = Base::N;
    }
    let mut candidates = Vec::new();
    append_local_flexible_proof_candidates(
        &read,
        &reference_bases,
        MateRescueWindow {
            contig_ordinal: 0,
            strand: BisulfiteStrand::OT,
            start: 0,
            end: 0,
        },
        5,
        &mut candidates,
    );
    assert!(candidates.iter().any(|candidate| {
        candidate.start() == 0 && candidate.proof_mask & FLEXIBLE_NOMINAL_PROOF != 0
    }));
}

#[test]
fn verifier_accepts_a_151_base_exact_alignment_without_allocation() {
    let bases = (0..151)
        .map(|position| [Base::A, Base::C, Base::G, Base::T][position % 4])
        .collect::<Vec<_>>();
    let reference = reference(&bases);
    let mut verifier = PlacementVerifier::new(&bases).expect("verifier");
    let observed = verifier
        .verify(
            &reference,
            ReadCandidate {
                contig_ordinal: 0,
                start: 0,
                strand: BisulfiteStrand::OT,
                proof_mask: 1,
            },
        )
        .expect("verification")
        .expect("exact hit");
    assert_eq!(observed.distance(), 0);
    assert_eq!(observed.tied_lengths(), 1 << 3);
}

#[test]
fn verifier_recovers_a_one_base_reference_insertion_endpoint() {
    let read = (0..151)
        .map(|position| [Base::A, Base::C, Base::G, Base::T][position % 4])
        .collect::<Vec<_>>();
    let mut inserted = read.clone();
    inserted.insert(73, Base::T);
    let reference = reference(&inserted);
    let mut verifier = PlacementVerifier::new(&read).expect("verifier");
    let observed = verifier
        .verify(
            &reference,
            ReadCandidate {
                contig_ordinal: 0,
                start: 0,
                strand: BisulfiteStrand::OT,
                proof_mask: 1,
            },
        )
        .expect("verification")
        .expect("one-edit hit");
    assert_eq!(observed.distance(), 1);
    assert_ne!(observed.tied_lengths() & (1 << 4), 0);
}

#[test]
fn verifier_applies_top_and_reverse_bottom_bisulfite_semantics() {
    let top_reference = reference(&[Base::C, Base::C, Base::C, Base::C]);
    let mut top = PlacementVerifier::new(&[Base::T, Base::T, Base::T, Base::T]).expect("top");
    assert_eq!(
        top.verify(
            &top_reference,
            ReadCandidate {
                contig_ordinal: 0,
                start: 0,
                strand: BisulfiteStrand::OT,
                proof_mask: 1,
            },
        )
        .expect("top verification")
        .expect("top hit")
        .distance(),
        0
    );

    let bottom_reference = reference(&[Base::G, Base::G, Base::G, Base::G]);
    let mut bottom = PlacementVerifier::new(&[Base::T, Base::T, Base::T, Base::T]).expect("bottom");
    assert_eq!(
        bottom
            .verify(
                &bottom_reference,
                ReadCandidate {
                    contig_ordinal: 0,
                    start: 0,
                    strand: BisulfiteStrand::OB,
                    proof_mask: 1,
                },
            )
            .expect("bottom verification")
            .expect("bottom hit")
            .distance(),
        0
    );
}

#[test]
fn ungapped_semi_global_chooses_bounded_five_prime_clip() {
    let reference_bases = vec![Base::A; 151];
    let reference = reference(&reference_bases);
    let mut read = reference_bases;
    read[..3].fill(Base::C);
    let placement = best_ungapped_semi_global_placement(
        &reference,
        &read,
        ReadCandidate {
            contig_ordinal: 0,
            start: 0,
            strand: BisulfiteStrand::OT,
            proof_mask: 1,
        },
        MAX_EDIT_DISTANCE,
        SEMI_GLOBAL_CLIP_PENALTY,
    )
    .expect("terminal mismatch run is clipped");
    assert_eq!(placement.start(), 3);
    assert_eq!(placement.end(), 151);
    assert_eq!(placement.distance(), 0);
    assert_eq!(placement.fallback_score, 3);
    assert_eq!(placement.retained_query_interval(151), 3..151);
}

#[test]
fn ungapped_semi_global_converts_reverse_endpoints_to_sequencing_coordinates() {
    let reference = reference(&[Base::A; 151]);
    let mut read = vec![Base::T; 151];
    read[..4].fill(Base::C);
    let placement = best_ungapped_semi_global_placement(
        &reference,
        &read,
        ReadCandidate {
            contig_ordinal: 0,
            start: 0,
            strand: BisulfiteStrand::OB,
            proof_mask: 1,
        },
        MAX_EDIT_DISTANCE,
        SEMI_GLOBAL_CLIP_PENALTY,
    )
    .expect("reverse sequencing 5-prime mismatch run is clipped");
    assert_eq!(placement.start(), 0);
    assert_eq!(placement.end(), 147);
    assert_eq!(placement.retained_query_interval(151), 4..151);
}

#[test]
fn linear_semi_global_endpoint_choice_matches_exhaustive_grid() {
    let reference = reference(&[Base::A; 151]);
    let candidate = ReadCandidate {
        contig_ordinal: 0,
        start: 0,
        strand: BisulfiteStrand::OT,
        proof_mask: 1,
    };
    let mut state = 0x6a09_e667_f3bc_c909_u64;
    for _ in 0..512 {
        let mut read = vec![Base::A; 151];
        for base in &mut read {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            if state.is_multiple_of(23) {
                *base = Base::C;
            }
        }
        let mut expected = None;
        for left in 0..=DEFAULT_MAX_SOFT_CLIP_BASES {
            for right in 0..=DEFAULT_MAX_SOFT_CLIP_BASES {
                let clipped = left + right;
                if clipped == 0
                    || read.len().saturating_sub(clipped) < SEMI_GLOBAL_MIN_ALIGNED_BASES
                {
                    continue;
                }
                let distance = u8::try_from(
                    read[left..read.len() - right]
                        .iter()
                        .filter(|&&base| base != Base::A)
                        .count(),
                )
                .expect("bounded distance");
                let score = distance
                    .saturating_mul(SEMI_GLOBAL_EDIT_PENALTY)
                    .saturating_add(
                        u8::try_from(clipped)
                            .expect("bounded clips")
                            .saturating_mul(SEMI_GLOBAL_CLIP_PENALTY),
                    );
                let admission_score = distance
                    .saturating_mul(SEMI_GLOBAL_ADMISSION_EDIT_PENALTY)
                    .saturating_add(
                        u8::try_from(clipped)
                            .expect("bounded clips")
                            .saturating_mul(SEMI_GLOBAL_CLIP_PENALTY),
                    );
                if distance <= MAX_EDIT_DISTANCE
                    && admission_score <= u8::try_from(read.len() / 5).expect("bounded score")
                {
                    let key = (score, clipped, distance, left, right);
                    if expected.is_none_or(|current| key < current) {
                        expected = Some(key);
                    }
                }
            }
        }
        let observed = best_ungapped_semi_global_placement(
            &reference,
            &read,
            candidate,
            MAX_EDIT_DISTANCE,
            SEMI_GLOBAL_CLIP_PENALTY,
        )
        .map(|placement| {
            let retained = placement.retained_query_interval(read.len());
            (
                placement.fallback_score,
                read.len() - (retained.end - retained.start),
                placement.distance(),
                retained.start,
                read.len() - retained.end,
            )
        });
        assert_eq!(
            observed,
            expected,
            "mismatch positions: {:?}",
            read.iter()
                .enumerate()
                .filter_map(|(position, &base)| (base != Base::A).then_some(position))
                .collect::<Vec<_>>()
        );
    }
}

#[test]
fn pair_selector_uses_direction_span_and_complete_best_score() {
    let mate1 = [
        placement(0, 100, 251, BisulfiteStrand::OT, 1),
        placement(0, 500, 651, BisulfiteStrand::OT, 0),
    ];
    let mate2 = [
        placement(0, 200, 351, BisulfiteStrand::CTOT, 1),
        placement(0, 600, 751, BisulfiteStrand::CTOT, 0),
    ];
    let mut best = Vec::new();
    let selection = select_best_pairs(&mate1, &mate2, MAX_EDIT_DISTANCE, 100, 300, &mut best);
    assert_eq!(selection.compatible_pairs, 2);
    assert_eq!(selection.second_best_distance, Some(2));
    assert_eq!(best.len(), 1);
    assert_eq!(best[0].template_start(), 500);
    assert_eq!(best[0].template_end(), 751);
    assert_eq!(best[0].distance(), 0);

    let exact_only = select_best_pairs(&mate1, &mate2, 0, 100, 300, &mut best);
    assert_eq!(exact_only.compatible_pairs, 1);
    assert_eq!(exact_only.second_best_distance, None);
    assert_eq!(best.len(), 1);
    assert_eq!(best[0].distance(), 0);
}

#[test]
fn ambiguous_representative_prefers_the_smallest_net_gap() {
    let ungapped = PairedPlacement {
        mate1: ReadPlacement::strict(0, 100, 250, BisulfiteStrand::OT, 1),
        mate2: ReadPlacement::strict(0, 300, 450, BisulfiteStrand::CTOT, 1),
        template_start: 100,
        template_end: 450,
        distance: 2,
        score: 2,
    };
    let net_deletion = PairedPlacement {
        mate1: ReadPlacement::strict(0, 99, 250, BisulfiteStrand::OT, 1),
        ..ungapped
    };
    let mut tied = [net_deletion, ungapped];
    prefer_minimum_net_gap_representative(&mut tied, 150, 150);
    assert_eq!(tied[0], ungapped);
    assert_eq!(tied[1], net_deletion);
}

#[test]
fn fallback_score_resolves_endpoint_ambiguity_with_a_better_clipped_origin() {
    let mate1 = [
        ReadPlacement::strict(0, 100, 251, BisulfiteStrand::OT, 2),
        ReadPlacement::strict(0, 101, 252, BisulfiteStrand::OT, 2),
        ReadPlacement {
            contig_ordinal: 0,
            start: 102,
            end: 251,
            strand: BisulfiteStrand::OT,
            distance: 0,
            query_start: 2,
            query_end: 151,
            fallback_score: 2,
        },
    ];
    let mate2 = [ReadPlacement::strict(0, 300, 451, BisulfiteStrand::CTOT, 0)];
    let mut best = Vec::new();
    select_best_pairs_with_fallback_score(&mate1, &mate2, 3, 0, 500, &mut best);
    assert_eq!(best.len(), 1);
    assert_eq!(best[0].mate1().start(), 102);
    assert_eq!(best[0].score(), 2);
}

#[test]
fn fallback_score_prefers_fewer_retained_edits_on_an_exact_score_tie() {
    let mate1 = [
        ReadPlacement::strict(0, 100, 251, BisulfiteStrand::OT, 1),
        ReadPlacement {
            contig_ordinal: 0,
            start: 107,
            end: 251,
            strand: BisulfiteStrand::OT,
            distance: 0,
            query_start: 7,
            query_end: 151,
            fallback_score: 7,
        },
    ];
    let mate2 = [ReadPlacement::strict(0, 300, 451, BisulfiteStrand::CTOT, 0)];
    let mut best = Vec::new();
    let selection = select_best_pairs_with_fallback_score(&mate1, &mate2, 3, 0, 500, &mut best);
    assert_eq!(best.len(), 1);
    assert_eq!(best[0].mate1().retained_query_interval(151), 7..151);
    assert_eq!(selection.second_best_distance, Some(7));
}

#[test]
fn origin_grouping_fast_path_proves_distinct_small_frontiers() {
    let complete = ReadPlacement::strict(0, 100, 251, BisulfiteStrand::OT, 1);
    let distinct = ReadPlacement::strict(0, 200, 351, BisulfiteStrand::OT, 1);
    let same_origin_endpoint = ReadPlacement {
        contig_ordinal: 0,
        start: 101,
        end: 251,
        strand: BisulfiteStrand::OT,
        distance: 0,
        query_start: 1,
        query_end: 151,
        fallback_score: SENSITIVE_CLIP_PENALTY,
    };

    assert!(!placements_may_share_origin(&[complete], 151));
    assert!(!placements_may_share_origin(&[complete, distinct], 151));
    assert!(placements_may_share_origin(
        &[complete, same_origin_endpoint],
        151
    ));
}

#[test]
fn origin_grouping_preserves_raw_selection_and_groups_only_mapq_evidence() {
    let mut read1 = vec![Base::A; 151];
    read1[0] = Base::C;
    let read2 = vec![Base::T; 151];
    let complete = ReadPlacement::strict(0, 100, 251, BisulfiteStrand::OT, 1);
    let clipped = ReadPlacement {
        contig_ordinal: 0,
        start: 101,
        end: 251,
        strand: BisulfiteStrand::OT,
        distance: 0,
        query_start: 1,
        query_end: 151,
        fallback_score: SENSITIVE_CLIP_PENALTY,
    };
    let mate2 = [ReadPlacement::strict(0, 300, 451, BisulfiteStrand::CTOT, 0)];
    let mut origins = std::collections::HashMap::new();
    let mut best = Vec::new();
    let selection = select_best_pair_origins_with_endpoint_policy(
        &[complete, clipped],
        &mate2,
        [&read1, &read2],
        MAX_EDIT_DISTANCE,
        0,
        500,
        true,
        &mut origins,
        &mut best,
    );

    assert_eq!(selection.compatible_pairs, 2);
    assert_eq!(selection.mapq_compatible_pairs, 1);
    assert_eq!(selection.near_best_pairings, 0);
    assert_eq!(selection.mapq_near_best_pairings, 0);
    assert_eq!(best.len(), 1);
    assert_eq!(best[0].mate1(), clipped);
    assert_eq!(best[0].distance(), 0);
    assert_eq!(best[0].score(), SENSITIVE_CLIP_PENALTY);
}

#[test]
fn origin_grouping_clips_a_terminal_mismatch_run_but_preserves_distinct_loci() {
    let read1 = vec![Base::A; 151];
    let read2 = vec![Base::T; 151];
    let endpoint_pair = |origin| {
        [
            ReadPlacement::strict(origin, 100, 251, BisulfiteStrand::OT, 2),
            ReadPlacement {
                contig_ordinal: origin,
                start: 102,
                end: 251,
                strand: BisulfiteStrand::OT,
                distance: 0,
                query_start: 2,
                query_end: 151,
                fallback_score: 2 * SENSITIVE_CLIP_PENALTY,
            },
        ]
    };
    let first = endpoint_pair(0);
    let second = endpoint_pair(1);
    let mate1 = [first[0], first[1], second[0], second[1]];
    let mate2 = [
        ReadPlacement::strict(0, 300, 451, BisulfiteStrand::CTOT, 0),
        ReadPlacement::strict(1, 300, 451, BisulfiteStrand::CTOT, 0),
    ];
    let mut origins = std::collections::HashMap::new();
    let mut best = Vec::new();
    let selection = select_best_pair_origins_with_endpoint_policy(
        &mate1,
        &mate2,
        [&read1, &read2],
        MAX_EDIT_DISTANCE,
        0,
        500,
        true,
        &mut origins,
        &mut best,
    );

    assert_eq!(selection.compatible_pairs, 4);
    assert_eq!(selection.mapq_compatible_pairs, 2);
    assert_eq!(selection.near_best_pairings, 1);
    assert_eq!(selection.mapq_near_best_pairings, 1);
    assert_eq!(best.len(), 2);
    assert!(best.iter().all(|pair| {
        pair.mate1().retained_query_interval(151) == (2..151) && pair.score() == 8
    }));
    assert_ne!(
        pair_origin_key(best[0], 151, 151),
        pair_origin_key(best[1], 151, 151)
    );
}

#[test]
fn reported_origin_endpoint_does_not_clip_unsupported_terminal_errors() {
    let reference = reference(&vec![Base::A; 600]);
    let mut read1 = vec![Base::A; 151];
    read1[..2].fill(Base::C);
    let read2 = vec![Base::T; 151];
    let selected = PairedPlacement {
        mate1: ReadPlacement::strict(0, 100, 251, BisulfiteStrand::OT, 2),
        mate2: ReadPlacement::strict(0, 300, 451, BisulfiteStrand::CTOT, 0),
        template_start: 100,
        template_end: 451,
        distance: 2,
        score: 14,
    };

    let reported = select_reported_origin_endpoint(
        &reference,
        [&read1, &read2],
        selected,
        MAX_EDIT_DISTANCE,
        0,
        500,
        &AlignmentOutputPolicy::default(),
    );

    assert_eq!(
        pair_origin_key(reported, read1.len(), read2.len()),
        pair_origin_key(selected, read1.len(), read2.len())
    );
    assert_eq!(reported.mate1(), selected.mate1());
    assert_eq!(reported.score(), selected.score());
}

#[test]
fn reported_origin_endpoint_keeps_an_isolated_terminal_mismatch_aligned() {
    let reference = reference(&vec![Base::A; 600]);
    let mut read1 = vec![Base::A; 151];
    read1[0] = Base::C;
    let read2 = vec![Base::T; 151];
    let clipped = ReadPlacement {
        contig_ordinal: 0,
        start: 101,
        end: 251,
        strand: BisulfiteStrand::OT,
        distance: 0,
        query_start: 1,
        query_end: 151,
        fallback_score: SENSITIVE_CLIP_PENALTY,
    };
    let selected = PairedPlacement {
        mate1: clipped,
        mate2: ReadPlacement::strict(0, 300, 451, BisulfiteStrand::CTOT, 0),
        template_start: 101,
        template_end: 451,
        distance: 0,
        score: SENSITIVE_CLIP_PENALTY,
    };

    let reported = select_reported_origin_endpoint(
        &reference,
        [&read1, &read2],
        selected,
        MAX_EDIT_DISTANCE,
        0,
        500,
        &AlignmentOutputPolicy::default(),
    );

    assert_eq!(
        pair_origin_key(reported, read1.len(), read2.len()),
        pair_origin_key(selected, read1.len(), read2.len())
    );
    assert_eq!(reported.mate1().start(), 100);
    assert_eq!(reported.mate1().retained_query_interval(151), 0..151);
    assert_eq!(reported.mate1().distance(), 1);
    assert_eq!(reported.score(), selected.score());
}

#[test]
fn endpoint_policy_recognizes_supported_three_prime_adapter_sequence() {
    let mut adapter_read = vec![Base::A; 100];
    let adapter = [
        Base::A,
        Base::G,
        Base::A,
        Base::T,
        Base::C,
        Base::G,
        Base::G,
        Base::A,
        Base::A,
        Base::G,
        Base::A,
        Base::G,
        Base::C,
    ];
    adapter_read[87..].copy_from_slice(&adapter);
    let supported = ReadPlacement {
        contig_ordinal: 0,
        start: 0,
        end: 87,
        strand: BisulfiteStrand::OT,
        distance: 0,
        query_start: 0,
        query_end: 87,
        fallback_score: u8::MAX,
    };
    let unsupported_read = vec![Base::A; 100];

    let output_policy = AlignmentOutputPolicy::default();
    assert!(sequencing_three_prime_adapter_supported(
        &adapter_read,
        87,
        &output_policy,
    ));
    assert_eq!(
        supported_three_prime_adapter_start(&adapter_read, &output_policy),
        Some(87),
    );
    assert_eq!(
        placement_endpoint_cost(&adapter_read, supported),
        ORIGIN_ENDPOINT_ADAPTER_CLIP_OPEN_PENALTY
    );
    assert_eq!(
        placement_endpoint_cost(&unsupported_read, supported),
        ORIGIN_ENDPOINT_CLIP_OPEN_PENALTY + 12 * ORIGIN_ENDPOINT_CLIP_EXTENSION_PENALTY
    );

    let mut partial_adapter_read = vec![Base::A; 100];
    partial_adapter_read[90..].copy_from_slice(&adapter[..10]);
    assert_eq!(
        supported_three_prime_adapter_start(&partial_adapter_read, &output_policy),
        Some(90)
    );
    partial_adapter_read[90] = Base::C;
    assert_eq!(
        supported_three_prime_adapter_start(&partial_adapter_read, &output_policy),
        None
    );
}

#[test]
fn reported_origin_endpoint_clips_an_explicit_three_prime_adapter() {
    let reference = reference(&vec![Base::A; 600]);
    let mut read1 = vec![Base::A; 151];
    let adapter = [
        Base::A,
        Base::G,
        Base::A,
        Base::T,
        Base::C,
        Base::G,
        Base::G,
        Base::A,
        Base::A,
        Base::G,
        Base::A,
        Base::G,
        Base::C,
    ];
    read1[138..].copy_from_slice(&adapter);
    let read2 = vec![Base::T; 151];
    let selected = PairedPlacement {
        mate1: ReadPlacement::strict(0, 100, 251, BisulfiteStrand::OT, 13),
        mate2: ReadPlacement::strict(0, 300, 451, BisulfiteStrand::CTOT, 0),
        template_start: 100,
        template_end: 451,
        distance: 13,
        score: 91,
    };

    let reported = select_reported_origin_endpoint(
        &reference,
        [&read1, &read2],
        selected,
        MAX_EDIT_DISTANCE,
        0,
        500,
        &AlignmentOutputPolicy::default(),
    );

    assert_eq!(
        pair_origin_key(reported, read1.len(), read2.len()),
        pair_origin_key(selected, read1.len(), read2.len())
    );
    assert_eq!(reported.mate1().retained_query_interval(151), 0..138);
    assert_eq!(reported.mate1().distance(), 0);
    assert_eq!(reported.score(), selected.score());
}

#[test]
fn adapter_only_endpoint_policy_never_adds_generic_five_prime_clipping() {
    let reference = reference(&vec![Base::A; 600]);
    let mut read1 = vec![Base::A; 151];
    read1[..2].fill(Base::C);
    read1[138..].copy_from_slice(&[
        Base::A,
        Base::G,
        Base::A,
        Base::T,
        Base::C,
        Base::G,
        Base::G,
        Base::A,
        Base::A,
        Base::G,
        Base::A,
        Base::G,
        Base::C,
    ]);
    let read2 = vec![Base::T; 151];
    let selected = PairedPlacement {
        mate1: ReadPlacement::strict(0, 100, 251, BisulfiteStrand::OT, 15),
        mate2: ReadPlacement::strict(0, 300, 451, BisulfiteStrand::CTOT, 0),
        template_start: 100,
        template_end: 451,
        distance: 15,
        score: 105,
    };
    let policy = AlignmentOutputPolicy::new(
        Some(b"AGATCGGAAGAGC"),
        8,
        30,
        crate::SoftClipMode::Adapter,
        30,
    )
    .expect("adapter-only policy");

    let reported = select_reported_origin_endpoint(
        &reference,
        [&read1, &read2],
        selected,
        MAX_EDIT_DISTANCE,
        0,
        500,
        &policy,
    );

    assert_eq!(reported.mate1().retained_query_interval(151), 0..138);
    assert_eq!(reported.mate1().distance(), 2);
}

#[test]
fn fallback_score_keeps_equal_biological_loci_ambiguous() {
    let clipped = |start, strand| ReadPlacement {
        contig_ordinal: 0,
        start,
        end: start + 149,
        strand,
        distance: 0,
        query_start: 2,
        query_end: 151,
        fallback_score: 2,
    };
    let mate1 = [
        clipped(100, BisulfiteStrand::OT),
        clipped(500, BisulfiteStrand::OT),
    ];
    let mate2 = [
        clipped(300, BisulfiteStrand::CTOT),
        clipped(700, BisulfiteStrand::CTOT),
    ];
    let mut best = Vec::new();
    select_best_pairs_with_fallback_score(&mate1, &mate2, 3, 0, 500, &mut best);
    assert_eq!(best.len(), 2);
    assert_ne!(best[0].template_start(), best[1].template_start());
}

#[test]
fn equivalent_cigar_endpoints_share_one_five_prime_origin() {
    let clipped = |start, query_start| ReadPlacement {
        contig_ordinal: 0,
        start,
        end: start + 149,
        strand: BisulfiteStrand::OT,
        distance: 0,
        query_start,
        query_end: 151,
        fallback_score: 2,
    };
    let mate1 = [clipped(102, 2), clipped(103, 3)];
    let mate2 = [ReadPlacement::strict(0, 300, 451, BisulfiteStrand::CTOT, 0)];
    let mut best = Vec::new();
    select_best_pairs_with_fallback_score(&mate1, &mate2, 3, 0, 500, &mut best);
    assert_eq!(best.len(), 2);
    collapse_equivalent_pair_origins(&mut best, 151, 151, false);
    assert_eq!(best.len(), 1);
    assert_eq!(placement_origin_key(best[0].mate1(), 151).2, 100);
}

#[test]
fn equivalent_pair_origins_retain_the_minimum_net_gap_endpoint() {
    let mate2 = ReadPlacement::strict(0, 300, 450, BisulfiteStrand::CTOT, 0);
    let gapped = PairedPlacement {
        mate1: ReadPlacement::strict(0, 100, 249, BisulfiteStrand::OT, 1),
        mate2,
        template_start: 100,
        template_end: 450,
        distance: 1,
        score: 1,
    };
    let ungapped = PairedPlacement {
        mate1: ReadPlacement::strict(0, 100, 250, BisulfiteStrand::OT, 1),
        ..gapped
    };
    let mut best = vec![gapped, ungapped];

    collapse_equivalent_pair_origins(&mut best, 150, 150, true);

    assert_eq!(best, [ungapped]);
}

#[test]
fn distinct_pair_origins_can_report_the_minimum_net_gap_representative() {
    let gapped = PairedPlacement {
        mate1: ReadPlacement::strict(0, 100, 249, BisulfiteStrand::OT, 1),
        mate2: ReadPlacement::strict(0, 300, 450, BisulfiteStrand::CTOT, 0),
        template_start: 100,
        template_end: 450,
        distance: 1,
        score: 1,
    };
    let ungapped = PairedPlacement {
        mate1: ReadPlacement::strict(0, 500, 650, BisulfiteStrand::OT, 1),
        mate2: ReadPlacement::strict(0, 700, 850, BisulfiteStrand::CTOT, 0),
        template_start: 500,
        template_end: 850,
        ..gapped
    };
    let mut best = vec![gapped, ungapped];

    collapse_equivalent_pair_origins(&mut best, 150, 150, true);

    assert_eq!(best.len(), 2);
    assert_eq!(best[0], ungapped);
}

#[test]
fn local_eight_block_filter_preserves_three_edits_and_rejects_four_blocks() {
    let read = vec![Base::A; 151];
    let candidate = ReadCandidate {
        contig_ordinal: 0,
        start: 0,
        strand: BisulfiteStrand::OT,
        proof_mask: 1,
    };
    let mut reference_bases = read.clone();
    reference_bases[5] = Base::C;
    reference_bases[35] = Base::C;
    let reference_index = reference(&reference_bases);
    let filter = LocalCandidateFilter::new(&read, candidate.strand());
    assert!(filter.supports(&reference_index, candidate));
    reference_bases[65] = Base::C;
    let reference_index = reference(&reference_bases);
    assert!(filter.supports(&reference_index, candidate));
    let mut four_destroyed_blocks = read.clone();
    for ordinal in [0_usize, 2, 4, 6] {
        let start = ordinal * read.len() / LOCAL_FILTER_BLOCKS;
        let end = (ordinal + 1) * read.len() / LOCAL_FILTER_BLOCKS;
        four_destroyed_blocks[start..end].fill(Base::C);
    }
    let reference_index = reference(&four_destroyed_blocks);
    assert!(!filter.supports(&reference_index, candidate));

    let mut inserted = read.clone();
    inserted.insert(73, Base::C);
    let reference_index = reference(&inserted);
    assert!(filter.supports(&reference_index, candidate));
}

#[test]
fn affine_reranking_is_limited_to_structural_ambiguous_pairs() {
    let mut workspace = PairWorkspace::with_capacity(8, 8, 4);
    let clipped = ReadPlacement {
        contig_ordinal: 0,
        start: 100,
        end: 250,
        strand: BisulfiteStrand::OT,
        distance: 0,
        query_start: 1,
        query_end: 151,
        fallback_score: SENSITIVE_CLIP_PENALTY,
    };
    workspace.best_pairs.push(PairedPlacement {
        mate1: clipped,
        mate2: ReadPlacement::strict(0, 300, 451, BisulfiteStrand::CTOT, 0),
        template_start: 100,
        template_end: 451,
        distance: 0,
        score: SENSITIVE_CLIP_PENALTY,
    });
    let metrics = PairAlignmentMetrics {
        compatible_pairs: 2,
        ..empty_pair_metrics()
    };

    assert!(workspace.should_affine_rescore(PairMappingStatus::Ambiguous, metrics, 151, 151,));
    assert!(!workspace.should_affine_rescore(PairMappingStatus::Unique, metrics, 151, 151,));

    workspace.best_pairs[0].mate1 = ReadPlacement::strict(0, 100, 251, BisulfiteStrand::OT, 0);
    assert!(!workspace.should_affine_rescore(PairMappingStatus::Ambiguous, metrics, 151, 151,));
}

#[test]
fn affine_score_uses_bwa_penalties_and_bisulfite_zero_cost_matches() {
    let reference_index = reference(&[Base::C, Base::C, Base::C, Base::C]);
    let placement = ReadPlacement::strict(0, 0, 4, BisulfiteStrand::OT, 0);
    let mut workspace = AffineScoreWorkspace::default();
    assert_eq!(
        affine_placement_score(
            &reference_index,
            &[Base::T, Base::T, Base::T, Base::T],
            placement,
            SENSITIVE_CLIP_PENALTY,
            &mut workspace,
        )
        .expect("conversion-aware exact affine score"),
        4
    );
    assert_eq!(
        affine_placement_score(
            &reference_index,
            &[Base::T, Base::T, Base::A, Base::T],
            placement,
            SENSITIVE_CLIP_PENALTY,
            &mut workspace,
        )
        .expect("one-mismatch affine score"),
        -1
    );
}

#[test]
fn selective_unmapped_deepening_requires_two_incomplete_frontiers() {
    let selection = |retained_hits, complete| RankedBlockSelection {
        retained_hits,
        complete,
    };
    let required = |first, second| incomplete_unmapped_frontier_deepening_required([first, second]);
    assert!(required(
        Some(selection(1, false)),
        Some(selection(u64::MAX, false))
    ));
    assert!(!required(
        Some(selection(64, true)),
        Some(selection(64, false))
    ));
    assert!(!required(Some(selection(128, false)), None));
}

#[test]
fn verification_cache_is_exact_and_read_scoped() {
    assert_eq!(core::mem::size_of::<VerificationCacheEntry>(), 32);
    let mut workspace = ReadWorkspace::with_capacity(8, 8);
    let candidate = ReadCandidate {
        contig_ordinal: 3,
        start: 100,
        strand: BisulfiteStrand::OT,
        proof_mask: FLEXIBLE_NOMINAL_PROOF,
    };
    let expected = ReadPlacement::strict(3, 99, 250, BisulfiteStrand::OT, 2);

    workspace.begin_verification_cache_read();
    workspace.placements.push(expected);
    ReadWorkspace::cache_candidate_verification(
        &mut workspace.verification_cache,
        &mut workspace.verification_cache_placements,
        workspace.verification_cache_generation,
        &mut workspace.verification_cache_population,
        &workspace.placements,
        candidate,
        INITIAL_EDIT_DISTANCE,
        true,
        0,
    );

    workspace.placements.clear();
    workspace.candidates.push(candidate);
    workspace.retain_uncached_candidates(INITIAL_EDIT_DISTANCE, true);
    assert!(workspace.candidates.is_empty());
    assert_eq!(workspace.placements, [expected]);

    workspace.placements.clear();
    workspace.candidates.push(candidate);
    workspace.retain_uncached_candidates(MAX_EDIT_DISTANCE, true);
    assert_eq!(workspace.candidates, [candidate]);
    assert!(workspace.placements.is_empty());

    workspace.candidates.clear();
    workspace.begin_verification_cache_read();
    workspace.candidates.push(candidate);
    workspace.retain_uncached_candidates(INITIAL_EDIT_DISTANCE, true);
    assert_eq!(workspace.candidates, [candidate]);
    assert!(workspace.placements.is_empty());
}

fn placement(
    contig_ordinal: u64,
    start: u64,
    end: u64,
    strand: BisulfiteStrand,
    distance: u8,
) -> ReadPlacement {
    ReadPlacement::strict(contig_ordinal, start, end, strand, distance)
}
