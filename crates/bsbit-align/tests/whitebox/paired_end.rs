//! White-box tests for canonical paired-end mapping.
//!
//! Kept outside implementation `src/` while remaining a child module so private
//! invariants can be tested without widening the crate API.

use super::*;
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
        PairedLibraryProfile::Directional,
        PairedSearchMode::Default,
        0,
        1_000,
    );
    let sensitive = PairedAlignmentOptions::primary(
        PairedLibraryProfile::Directional,
        PairedSearchMode::Sensitive,
        0,
        1_000,
    );
    let trimmed = PairedAlignmentOptions::adapter_trimmed(
        PairedLibraryProfile::Directional,
        PairedSearchMode::Sensitive,
        0,
        1_000,
    );
    assert_eq!(
        default.derived_policy(),
        (PAIRED_MAX_EDIT_DISTANCE, false, false)
    );
    assert_eq!(
        sensitive.derived_policy(),
        (PAIRED_MAX_EDIT_DISTANCE, true, true)
    );
    assert_eq!(
        trimmed.derived_policy(),
        (PAIRED_MAX_EDIT_DISTANCE, true, false)
    );
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
fn sensitive_repeat_recheck_is_limited_to_suspicious_unique_results() {
    let below_threshold = PairAlignmentMetrics {
        mate1: ReadAlignmentMetrics {
            located_rows: SENSITIVE_REPEAT_RECHECK_ROWS - 1,
            ..ReadAlignmentMetrics::default()
        },
        ..empty_pair_metrics()
    };
    assert!(!sensitive_repeat_recheck_required(
        PairMappingStatus::Unique,
        below_threshold,
    ));

    let at_threshold = PairAlignmentMetrics {
        mate2: ReadAlignmentMetrics {
            located_rows: SENSITIVE_REPEAT_RECHECK_ROWS,
            ..ReadAlignmentMetrics::default()
        },
        ..below_threshold
    };
    assert!(sensitive_repeat_recheck_required(
        PairMappingStatus::Unique,
        at_threshold,
    ));
    assert!(!sensitive_repeat_recheck_required(
        PairMappingStatus::Ambiguous,
        at_threshold,
    ));

    let rescued = PairAlignmentMetrics {
        window_rescue_attempted: true,
        ..below_threshold
    };
    assert!(sensitive_repeat_recheck_required(
        PairMappingStatus::Unique,
        rescued,
    ));
}

#[test]
fn targeted_semi_global_admits_only_the_frozen_incomplete_sparse_cell() {
    let complete_ambiguous = PairAlignmentMetrics {
        best_pair_placements: 2,
        best_pair_score: Some(100),
        frontier_complete: true,
        ..empty_pair_metrics()
    };
    assert!(sensitive_targeted_semi_global_required(
        PairMappingStatus::Ambiguous,
        complete_ambiguous,
        None,
    ));

    let incomplete = PairAlignmentMetrics {
        frontier_complete: false,
        ..complete_ambiguous
    };
    assert!(!sensitive_targeted_semi_global_required(
        PairMappingStatus::Ambiguous,
        incomplete,
        Some(2),
    ));

    let incomplete_sparse = PairAlignmentMetrics {
        mate1: ReadAlignmentMetrics {
            located_rows: 4,
            ..ReadAlignmentMetrics::default()
        },
        near_best_pairings: 0,
        second_best_pair_score: None,
        ..incomplete
    };
    assert!(sensitive_targeted_semi_global_required(
        PairMappingStatus::Ambiguous,
        incomplete_sparse,
        Some(2),
    ));
    assert!(!sensitive_targeted_semi_global_required(
        PairMappingStatus::Ambiguous,
        incomplete_sparse,
        Some(4),
    ));
    assert!(!sensitive_targeted_semi_global_required(
        PairMappingStatus::Unmapped,
        complete_ambiguous,
        None,
    ));

    let high_confidence = PairAlignmentMetrics {
        best_pair_placements: 1,
        ..complete_ambiguous
    };
    assert!(!sensitive_targeted_semi_global_required(
        PairMappingStatus::Unique,
        high_confidence,
        None,
    ));
}

#[test]
fn ambiguity_q10_uses_only_the_three_frozen_integer_cells() {
    let original = PairAlignmentMetrics {
        best_pair_placements: 3,
        frontier_complete: true,
        ..empty_pair_metrics()
    };
    let candidate = PairAlignmentMetrics {
        mate1: ReadAlignmentMetrics {
            located_rows: 3,
            ..ReadAlignmentMetrics::default()
        },
        mate2: ReadAlignmentMetrics {
            located_rows: 9,
            ..ReadAlignmentMetrics::default()
        },
        near_best_pairings: 1,
        frontier_complete: true,
        ..empty_pair_metrics()
    };
    assert!(sensitive_ambiguity_q10_certified(
        PairMappingStatus::Ambiguous,
        original,
        Some(7),
        Some(2),
        PairMappingStatus::Ambiguous,
        candidate,
        Some(3),
        (Some(0), Some(1), 1),
        false,
    ));
    assert!(sensitive_ambiguity_q10_certified(
        PairMappingStatus::Ambiguous,
        original,
        Some(7),
        Some(2),
        PairMappingStatus::Ambiguous,
        candidate,
        Some(3),
        (Some(0), None, 2),
        true,
    ));
    assert!(!sensitive_ambiguity_q10_certified(
        PairMappingStatus::Ambiguous,
        original,
        Some(7),
        Some(2),
        PairMappingStatus::Ambiguous,
        candidate,
        Some(3),
        (Some(0), None, 2),
        false,
    ));

    let incomplete_sparse = PairAlignmentMetrics {
        mate1: ReadAlignmentMetrics {
            located_rows: 4,
            ..ReadAlignmentMetrics::default()
        },
        best_pair_placements: 2,
        near_best_pairings: 0,
        second_best_pair_score: None,
        frontier_complete: false,
        ..empty_pair_metrics()
    };
    let completed_unique = PairAlignmentMetrics {
        best_pair_placements: 1,
        ..candidate
    };
    assert!(sensitive_ambiguity_q10_certified(
        PairMappingStatus::Ambiguous,
        incomplete_sparse,
        Some(4),
        Some(2),
        PairMappingStatus::Unique,
        completed_unique,
        Some(2),
        (Some(0), None, 2),
        true,
    ));
    assert!(!sensitive_ambiguity_q10_certified(
        PairMappingStatus::Ambiguous,
        incomplete_sparse,
        Some(4),
        Some(2),
        PairMappingStatus::Unique,
        completed_unique,
        Some(2),
        (Some(0), None, 2),
        false,
    ));
}

#[test]
// This table-driven policy test keeps all Q20 boundary counterexamples beside
// the positive case so omissions are visible in one assertion group.
#[allow(clippy::too_many_lines)]
fn stable_one_mate_rescue_certifies_only_the_q20_boundary() {
    let original = PairAlignmentMetrics {
        window_rescue_attempted: true,
        best_pair_score: Some(100),
        frontier_complete: true,
        ..empty_pair_metrics()
    };
    let moderate_candidate = PairAlignmentMetrics {
        mate1: ReadAlignmentMetrics {
            located_rows: 0,
            verified_placements: 2,
            ..ReadAlignmentMetrics::default()
        },
        mate2: ReadAlignmentMetrics {
            located_rows: 8,
            verified_placements: 4,
            ..ReadAlignmentMetrics::default()
        },
        window_rescue_attempted: true,
        best_pair_score: Some(100),
        frontier_complete: true,
        ..empty_pair_metrics()
    };
    assert!(sensitive_stable_rescue_q20_certified(
        PairMappingStatus::Unique,
        original,
        19,
        PairMappingStatus::Unique,
        moderate_candidate,
        19,
        true,
    ));
    assert!(!sensitive_stable_rescue_q20_certified(
        PairMappingStatus::Unique,
        original,
        19,
        PairMappingStatus::Unique,
        moderate_candidate,
        19,
        false,
    ));

    let two_informative_mates = PairAlignmentMetrics {
        mate1: ReadAlignmentMetrics {
            located_rows: 1,
            ..moderate_candidate.mate1
        },
        ..moderate_candidate
    };
    assert!(!sensitive_stable_rescue_q20_certified(
        PairMappingStatus::Unique,
        original,
        19,
        PairMappingStatus::Unique,
        two_informative_mates,
        19,
        true,
    ));

    let below_row_floor = PairAlignmentMetrics {
        mate2: ReadAlignmentMetrics {
            located_rows: 4,
            ..moderate_candidate.mate2
        },
        ..moderate_candidate
    };
    assert!(!sensitive_stable_rescue_q20_certified(
        PairMappingStatus::Unique,
        original,
        19,
        PairMappingStatus::Unique,
        below_row_floor,
        19,
        true,
    ));

    let above_verified_ceiling = PairAlignmentMetrics {
        mate1: ReadAlignmentMetrics {
            verified_placements: 27,
            ..moderate_candidate.mate1
        },
        mate2: ReadAlignmentMetrics {
            verified_placements: 28,
            ..moderate_candidate.mate2
        },
        ..moderate_candidate
    };
    assert!(!sensitive_stable_rescue_q20_certified(
        PairMappingStatus::Unique,
        original,
        19,
        PairMappingStatus::Unique,
        above_verified_ceiling,
        19,
        true,
    ));
    assert!(!sensitive_stable_rescue_q20_certified(
        PairMappingStatus::Unique,
        original,
        19,
        PairMappingStatus::Unique,
        moderate_candidate,
        18,
        true,
    ));

    let sparse_candidate = PairAlignmentMetrics {
        mate1: ReadAlignmentMetrics {
            located_rows: 0,
            verified_placements: 1,
            ..ReadAlignmentMetrics::default()
        },
        mate2: ReadAlignmentMetrics {
            located_rows: 8,
            verified_placements: 2,
            ..ReadAlignmentMetrics::default()
        },
        ..moderate_candidate
    };
    assert!(sensitive_stable_rescue_q20_certified(
        PairMappingStatus::Unique,
        original,
        19,
        PairMappingStatus::Unique,
        sparse_candidate,
        19,
        true,
    ));
    let weak_gap = PairAlignmentMetrics {
        second_best_pair_score: Some(95),
        ..sparse_candidate
    };
    assert!(!sensitive_stable_rescue_q20_certified(
        PairMappingStatus::Unique,
        original,
        19,
        PairMappingStatus::Unique,
        weak_gap,
        19,
        true,
    ));
}

#[test]
// The complete and incomplete parsimony counterexamples intentionally share
// one fixture to make the certificate boundary explicit.
#[allow(clippy::too_many_lines)]
fn bounded_two_way_parsimony_certifies_only_qualified_complete_ties() {
    let minimum_gap_pair = PairedPlacement {
        mate1: placement(0, 100, 251, BisulfiteStrand::OT, 1),
        mate2: placement(0, 300, 451, BisulfiteStrand::CTOT, 1),
        template_start: 100,
        template_end: 451,
        distance: 2,
        score: 14,
    };
    let larger_gap_pair = PairedPlacement {
        mate1: placement(0, 100, 252, BisulfiteStrand::OT, 1),
        ..minimum_gap_pair
    };
    let candidate_pairs = [minimum_gap_pair, larger_gap_pair];
    let gap_profile = pair_net_gap_profile(&candidate_pairs, 151, 151);
    assert_eq!(gap_profile, (Some(0), Some(1), 1));

    let original = PairAlignmentMetrics {
        mate1: ReadAlignmentMetrics {
            located_rows: SENSITIVE_PARSIMONY_MAX_LOCATED_ROWS,
            ..ReadAlignmentMetrics::default()
        },
        best_pair_placements: 2,
        best_pair_score: Some(100),
        second_best_pair_score: Some(100 - SENSITIVE_PARSIMONY_REQUIRED_SCORE_GAP),
        frontier_complete: true,
        ..empty_pair_metrics()
    };
    let candidate = PairAlignmentMetrics {
        mate1: ReadAlignmentMetrics {
            verified_placements: 3,
            ..ReadAlignmentMetrics::default()
        },
        mate2: ReadAlignmentMetrics {
            verified_placements: 4,
            ..ReadAlignmentMetrics::default()
        },
        best_pair_placements: 2,
        frontier_complete: true,
        ..empty_pair_metrics()
    };
    let certified = |original, candidate, candidate_best, profile, same_origin| {
        sensitive_two_way_parsimony_q20_certified(
            PairMappingStatus::Ambiguous,
            original,
            Some(minimum_gap_pair),
            PairMappingStatus::Ambiguous,
            candidate,
            candidate_best,
            2,
            profile,
            151,
            151,
            same_origin,
        )
    };
    assert!(certified(
        original,
        candidate,
        Some(minimum_gap_pair),
        gap_profile,
        true,
    ));
    assert!(!certified(
        original,
        candidate,
        Some(minimum_gap_pair),
        gap_profile,
        false,
    ));
    assert!(!certified(
        PairAlignmentMetrics {
            mate1: ReadAlignmentMetrics {
                located_rows: SENSITIVE_PARSIMONY_MAX_LOCATED_ROWS + 1,
                ..original.mate1
            },
            ..original
        },
        candidate,
        Some(minimum_gap_pair),
        gap_profile,
        true,
    ));
    assert!(!certified(
        original,
        PairAlignmentMetrics {
            mate2: ReadAlignmentMetrics {
                verified_placements: 5,
                ..candidate.mate2
            },
            ..candidate
        },
        Some(minimum_gap_pair),
        gap_profile,
        true,
    ));
    assert!(!certified(
        original,
        candidate,
        Some(larger_gap_pair),
        gap_profile,
        true,
    ));
    assert!(!certified(
        original,
        candidate,
        Some(minimum_gap_pair),
        (Some(0), None, 2),
        true,
    ));
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
            ..empty_pair_metrics()
        };
        let best_pair = Some(pair);
        PairedBatchResult {
            class: PairMappingStatus::Unique,
            metrics,
            best_pair,
            second_best_distance: None,
            repeat_risk_q20_certified: true,
            parsimony_q20_certified: true,
            ambiguity_q10_certified: false,
            requires_positive_mapq_for_reporting: false,
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
    assert!(!merged.repeat_risk_q20_certified());
    assert!(!merged.parsimony_q20_certified());

    let tied =
        merge_non_directional_batch_results(&result(directional_pair, -4, 1, 1), &complementary);
    assert_eq!(tied.class(), PairMappingStatus::Ambiguous);
    assert_eq!(tied.best_pair_score(), Some(-4));
    assert_eq!(tied.second_best_pair_score(), Some(-4));
    assert!(tied.metrics().best_pair_placements >= 2);
    assert_eq!(tied.evidence_mapping_quality(), 0);
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
        for left in 0..=SEMI_GLOBAL_MAX_CLIP_BASES {
            for right in 0..=SEMI_GLOBAL_MAX_CLIP_BASES {
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

    assert!(sequencing_three_prime_adapter_supported(&adapter_read, 87));
    assert_eq!(supported_three_prime_adapter_start(&adapter_read), Some(87));
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
        supported_three_prime_adapter_start(&partial_adapter_read),
        Some(90)
    );
    partial_adapter_read[90] = Base::C;
    assert_eq!(
        supported_three_prime_adapter_start(&partial_adapter_read),
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
    collapse_equivalent_pair_origins(&mut best, 151, 151);
    assert_eq!(best.len(), 1);
    assert_eq!(placement_origin_key(best[0].mate1(), 151).2, 100);
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
fn selective_unmapped_deepening_uses_a_bounded_reciprocal_hit_window() {
    let selection = |retained_hits, complete| RankedBlockSelection {
        retained_hits,
        complete,
    };
    let required =
        |first, second| selective_unmapped_frontier_deepening_required([first, second], None);
    let required_sum = |retained_hits| {
        let first = retained_hits / 2;
        required(
            Some(selection(first, false)),
            Some(selection(retained_hits - first, false)),
        )
    };

    if let Some(below_minimum) = SENSITIVE_SELECTIVE_UNMAPPED_MIN_RETAINED_HITS.checked_sub(1) {
        assert!(!required_sum(below_minimum));
    }
    assert!(required_sum(SENSITIVE_SELECTIVE_UNMAPPED_MIN_RETAINED_HITS));
    assert!(required_sum(
        SENSITIVE_SELECTIVE_UNMAPPED_MAX_RETAINED_HITS - 1
    ));
    assert!(!required_sum(
        SENSITIVE_SELECTIVE_UNMAPPED_MAX_RETAINED_HITS
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
    workspace.retain_uncached_candidates(PAIRED_MAX_EDIT_DISTANCE, true);
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
