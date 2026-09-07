//! White-box tests for canonical single-end mapping.

use super::*;
use crate::mapq_policy::MAPQ_POLICY;
use bsbit_core::alphabet::Base;
use bsbit_core::bisulfite::BisulfiteStrand;
use bsbit_core::sequence::NormalizedSequence;
use bsbit_index::reference::{ContigInput, ReferenceBuildLimits, ReferenceIndex};

fn test_reference(bases: &[Base]) -> ReferenceIndex {
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

fn mapped(start: u64, mapping_quality: u8) -> SingleAlignmentResult {
    mapped_at_distance(start, mapping_quality, 1)
}

fn mapped_at_distance(start: u64, mapping_quality: u8, edit_distance: u8) -> SingleAlignmentResult {
    mapped_on_strand(start, mapping_quality, edit_distance, BisulfiteStrand::OT)
}

fn mapped_on_strand(
    start: u64,
    mapping_quality: u8,
    edit_distance: u8,
    strand: BisulfiteStrand,
) -> SingleAlignmentResult {
    SingleAlignmentResult {
        status: SingleMappingStatus::Unique,
        placement: Some(ReadPlacement::strict(
            0,
            start,
            start + 100,
            strand,
            edit_distance,
        )),
        retained_query_end: 100,
        mapping_quality,
        located_rows: 1,
        distinct_candidate_starts: 1,
        verified_placements: 1,
        best_origin_count: 1,
        adapter_attempted: false,
        adapter_status: None,
        adapter_clipped_bases: 0,
    }
}

#[test]
fn non_directional_cross_pass_tie_remains_ambiguous() {
    let original = mapped_at_distance(100, 20, 2);
    let complementary = mapped_at_distance(500, 20, 2);
    let merged = merge_non_directional_results(original, complementary);
    assert_eq!(merged.status(), SingleMappingStatus::Ambiguous);
}

#[test]
fn sensitive_single_profile_completes_a_wider_bounded_frontier() {
    let default = SingleSearchMode::Default.limits();
    let sensitive = SingleSearchMode::Sensitive.limits();
    assert!(sensitive.maximum_seed_hits > default.maximum_seed_hits);
    assert_eq!(
        sensitive.maximum_combined_rescue_hits,
        default.maximum_combined_rescue_hits
    );
    assert!(sensitive.maximum_seed_rounds > default.maximum_seed_rounds);
    assert!(!SingleSearchMode::Default.completes_candidate_frontier());
    assert!(SingleSearchMode::Sensitive.completes_candidate_frontier());
}

#[test]
fn sensitive_low_confidence_conflict_preserves_the_incumbent_as_ambiguous() {
    let incumbent = mapped(100, 30);
    let completed = mapped(500, MAPQ_POLICY.single.sensitive_replacement_min_mapq - 1);
    let reconciled = SingleBatchAligner::reconcile_sensitive_result(incumbent, completed, 100);
    assert_eq!(reconciled.status(), SingleMappingStatus::Ambiguous);
    assert_eq!(reconciled.placement(), incumbent.placement());
    assert_eq!(reconciled.mapping_quality(), 0);
}

#[test]
fn sensitive_declarable_conflict_may_replace_the_incumbent() {
    let incumbent = mapped(100, 30);
    let completed = mapped(500, MAPQ_POLICY.single.sensitive_replacement_min_mapq);
    let reconciled = SingleBatchAligner::reconcile_sensitive_result(incumbent, completed, 100);
    assert_eq!(reconciled, completed);
}

#[test]
fn sensitive_low_confidence_rescue_remains_unmapped() {
    let incumbent = SingleAlignmentResult::unmapped(1, 0, 0);
    let completed = mapped(500, MAPQ_POLICY.single.sensitive_replacement_min_mapq - 1);
    let reconciled = SingleBatchAligner::reconcile_sensitive_result(incumbent, completed, 100);
    assert_eq!(reconciled.status(), SingleMappingStatus::Unmapped);
    assert_eq!(reconciled.placement(), None);
}

#[test]
fn sensitive_ambiguous_rescue_retains_a_zero_mapq_coordinate() {
    let incumbent = SingleAlignmentResult::unmapped(1, 0, 0);
    let mut completed = mapped(500, 0);
    completed.status = SingleMappingStatus::Ambiguous;
    completed.best_origin_count = 2;
    let reconciled = SingleBatchAligner::reconcile_sensitive_result(incumbent, completed, 100);
    assert_eq!(reconciled.status(), SingleMappingStatus::Ambiguous);
    assert_eq!(reconciled.placement(), completed.placement());
    assert_eq!(reconciled.mapping_quality(), 0);
}

#[test]
fn sensitive_frontier_audit_uses_the_complete_configured_distance_boundary() {
    assert_eq!(
        SingleBatchAligner::sensitive_audit_distance_limit(mapped(100, 20), 5),
        5
    );
    assert_eq!(
        SingleBatchAligner::sensitive_audit_distance_limit(mapped(100, 15), 5),
        5
    );
    assert_eq!(
        SingleBatchAligner::sensitive_audit_distance_limit(mapped_at_distance(100, 20, 3), 5),
        5
    );
    assert_eq!(
        SingleBatchAligner::sensitive_audit_distance_limit(
            SingleAlignmentResult::unmapped(0, 0, 0),
            5,
        ),
        5
    );
}

#[test]
fn sensitive_frontier_audits_every_search_outcome() {
    assert!(SingleBatchAligner::sensitive_audit_required(mapped(
        100, 15
    )));
    assert!(SingleBatchAligner::sensitive_audit_required(mapped(
        100, 20
    )));

    let mut ambiguous = mapped(100, 0);
    ambiguous.status = SingleMappingStatus::Ambiguous;
    assert!(SingleBatchAligner::sensitive_audit_required(ambiguous));
    assert!(SingleBatchAligner::sensitive_audit_required(
        SingleAlignmentResult::unmapped(0, 0, 0)
    ));
}

#[test]
fn highly_repetitive_best_distance_set_retains_a_zero_mapq_representative() {
    let mut workspace = ReadWorkspace::with_capacity(32, 32);
    for start in 0..=MAPQ_POLICY.single.affine_rerank_origin_limit {
        let start = u64::try_from(start).expect("small fixture start fits u64") * 1_000;
        workspace.placements.push(ReadPlacement::strict(
            0,
            start,
            start + 100,
            BisulfiteStrand::OT,
            1,
        ));
    }
    let mut origins = Vec::new();
    let completed = SingleBatchAligner::finish_result(
        &workspace,
        &mut origins,
        SingleResultEvidence {
            read_length: 100,
            metrics: ReadAlignmentMetrics {
                verified_placements: u64::try_from(workspace.placements.len())
                    .expect("small fixture length fits u64"),
                ..ReadAlignmentMetrics::default()
            },
            verified_distance_limit: INITIAL_EDIT_DISTANCE,
            first_seed: None,
            frontier_complete: false,
            confidence_audit_complete: false,
        },
    );
    let reconciled =
        SingleBatchAligner::reconcile_sensitive_result(mapped(10_000, 15), completed, 100);
    assert_eq!(reconciled.status(), SingleMappingStatus::Ambiguous);
    assert!(reconciled.placement().is_some());
    assert_eq!(reconciled.mapping_quality(), 0);

    let result = completed;
    assert_eq!(result.status(), SingleMappingStatus::Ambiguous);
    assert!(result.placement().is_some());
    assert_eq!(result.mapping_quality(), 0);
    assert_eq!(
        result.best_origin_count,
        u64::try_from(MAPQ_POLICY.single.affine_rerank_origin_limit + 1).unwrap()
    );

    workspace.placements.pop();
    let boundary = SingleBatchAligner::finish_result(
        &workspace,
        &mut origins,
        SingleResultEvidence {
            read_length: 100,
            metrics: ReadAlignmentMetrics::default(),
            verified_distance_limit: INITIAL_EDIT_DISTANCE,
            first_seed: None,
            frontier_complete: false,
            confidence_audit_complete: false,
        },
    );
    assert_eq!(boundary.status(), SingleMappingStatus::Ambiguous);
    assert!(boundary.placement().is_some());
    assert_eq!(boundary.mapping_quality(), 0);
}

#[test]
fn fair_reporting_tie_is_seeded_and_reference_order_independent() {
    let read = [Base::A; 20];
    let forward = repeated_named_reference([b"alpha", b"beta"]);
    let reversed = repeated_named_reference([b"beta", b"alpha"]);
    let workspace = |first_ordinal, second_ordinal| {
        let mut workspace = ReadWorkspace::with_capacity(2, 2);
        workspace.placements = vec![
            ReadPlacement::strict(first_ordinal, 11, 31, BisulfiteStrand::OT, 0),
            ReadPlacement::strict(second_ordinal, 11, 31, BisulfiteStrand::OT, 0),
        ];
        workspace
    };
    let ambiguous = |placement| SingleAlignmentResult {
        status: SingleMappingStatus::Ambiguous,
        placement: Some(placement),
        retained_query_end: read.len(),
        mapping_quality: 0,
        located_rows: 2,
        distinct_candidate_starts: 2,
        verified_placements: 2,
        best_origin_count: 2,
        adapter_attempted: false,
        adapter_status: None,
        adapter_clipped_bases: 0,
    };
    let forward_workspace = workspace(0, 1);
    let reversed_workspace = workspace(1, 0);
    let mut cached_forward_workspace = workspace(0, 1);
    cached_forward_workspace.affine_scores = cached_forward_workspace
        .placements
        .iter()
        .copied()
        .map(|placement| (placement, 20))
        .collect();
    let choose_name =
        |reference: &ReferenceIndex, workspace: &ReadWorkspace, seed: u64| -> Vec<u8> {
            let selected = prefer_fair_ambiguous_representative(
                reference,
                &read,
                workspace,
                ambiguous(workspace.placements[0]),
                ReportingTieBreak { seed, read_key: 77 },
            )
            .expect("fair tie-break succeeds")
            .placement()
            .expect("ambiguous representative is retained");
            reference
                .contig_by_ordinal(selected.contig_ordinal())
                .expect("selected contig exists")
                .name()
                .to_vec()
        };

    for seed in 0..32 {
        assert_eq!(
            choose_name(&forward, &forward_workspace, seed),
            choose_name(&reversed, &reversed_workspace, seed),
        );
        assert_eq!(
            choose_name(&forward, &forward_workspace, seed),
            choose_name(&forward, &cached_forward_workspace, seed),
        );
    }
    let first = choose_name(&forward, &forward_workspace, 0);
    let alternate = (1..256)
        .map(|seed| choose_name(&forward, &forward_workspace, seed))
        .find(|name| *name != first)
        .expect("a seed changes a fair two-way lottery");
    assert_ne!(first, alternate);

    let selected = prefer_fair_ambiguous_representative(
        &forward,
        &read,
        &forward_workspace,
        ambiguous(forward_workspace.placements[0]),
        ReportingTieBreak {
            seed: 0,
            read_key: 77,
        },
    )
    .expect("fair tie-break succeeds");
    assert_eq!(selected.status(), SingleMappingStatus::Ambiguous);
    assert_eq!(selected.mapping_quality(), 0);
    assert_eq!(selected.best_origin_count, 2);
}

#[test]
fn equally_scoring_endpoints_prefer_the_minimum_net_gap() {
    let mut workspace = ReadWorkspace::with_capacity(2, 2);
    let shorter_reference_span = ReadPlacement::strict(0, 100, 199, BisulfiteStrand::OT, 1);
    let ungapped_span = ReadPlacement::strict(0, 100, 200, BisulfiteStrand::OT, 1);
    workspace.placements.push(shorter_reference_span);
    workspace.placements.push(ungapped_span);

    let mut origins = Vec::new();
    let result = SingleBatchAligner::finish_result(
        &workspace,
        &mut origins,
        SingleResultEvidence {
            read_length: 100,
            metrics: ReadAlignmentMetrics {
                verified_placements: 2,
                ..ReadAlignmentMetrics::default()
            },
            verified_distance_limit: INITIAL_EDIT_DISTANCE,
            first_seed: None,
            frontier_complete: false,
            confidence_audit_complete: false,
        },
    );

    assert_eq!(result.status(), SingleMappingStatus::Unique);
    assert_eq!(result.placement(), Some(ungapped_span));
    assert_eq!(result.best_origin_count, 1);
}

#[test]
fn representative_reranking_preserves_distinct_origin_ambiguity() {
    let mut workspace = ReadWorkspace::with_capacity(2, 2);
    let canonical_but_gapped = ReadPlacement::strict(0, 100, 199, BisulfiteStrand::OT, 1);
    let later_ungapped = ReadPlacement::strict(0, 500, 600, BisulfiteStrand::OT, 1);
    workspace.placements.push(canonical_but_gapped);
    workspace.placements.push(later_ungapped);

    let mut origins = Vec::new();
    let result = SingleBatchAligner::finish_result(
        &workspace,
        &mut origins,
        SingleResultEvidence {
            read_length: 100,
            metrics: ReadAlignmentMetrics::default(),
            verified_distance_limit: INITIAL_EDIT_DISTANCE,
            first_seed: None,
            frontier_complete: false,
            confidence_audit_complete: false,
        },
    );

    assert_eq!(result.status(), SingleMappingStatus::Ambiguous);
    assert_eq!(result.placement(), Some(later_ungapped));
    assert_eq!(result.mapping_quality(), 0);
    assert_eq!(result.best_origin_count, 2);
}

#[test]
fn completed_frontier_collapses_indel_shifted_endpoints_into_one_q10_locus() {
    let mut workspace = ReadWorkspace::with_capacity(2, 2);
    let first = ReadPlacement::strict(0, 100, 200, BisulfiteStrand::OT, 2);
    let shifted = ReadPlacement::strict(0, 102, 202, BisulfiteStrand::OT, 2);
    workspace.placements.extend([first, shifted]);

    let mut origins = Vec::new();
    let incomplete = SingleBatchAligner::finish_result(
        &workspace,
        &mut origins,
        SingleResultEvidence {
            read_length: 100,
            metrics: ReadAlignmentMetrics::default(),
            verified_distance_limit: MAX_EDIT_DISTANCE,
            first_seed: None,
            frontier_complete: false,
            confidence_audit_complete: false,
        },
    );
    assert_eq!(incomplete.status(), SingleMappingStatus::Ambiguous);
    assert_eq!(incomplete.mapping_quality(), 0);
    assert_eq!(incomplete.best_origin_count, 2);

    let complete = SingleBatchAligner::finish_result(
        &workspace,
        &mut origins,
        SingleResultEvidence {
            read_length: 100,
            metrics: ReadAlignmentMetrics::default(),
            verified_distance_limit: MAX_EDIT_DISTANCE,
            first_seed: None,
            frontier_complete: true,
            confidence_audit_complete: true,
        },
    );
    assert_eq!(complete.status(), SingleMappingStatus::Unique);
    assert_eq!(
        complete.mapping_quality(),
        MAPQ_POLICY.single.local_locus_cap
    );
    assert_eq!(complete.best_origin_count, 1);
}

#[test]
fn completed_frontier_keeps_endpoints_beyond_edit_radius_ambiguous() {
    let mut workspace = ReadWorkspace::with_capacity(2, 2);
    workspace.placements.extend([
        ReadPlacement::strict(0, 100, 200, BisulfiteStrand::OT, 2),
        ReadPlacement::strict(0, 103, 203, BisulfiteStrand::OT, 2),
    ]);

    let mut origins = Vec::new();
    let complete = SingleBatchAligner::finish_result(
        &workspace,
        &mut origins,
        SingleResultEvidence {
            read_length: 100,
            metrics: ReadAlignmentMetrics::default(),
            verified_distance_limit: MAX_EDIT_DISTANCE,
            first_seed: None,
            frontier_complete: true,
            confidence_audit_complete: true,
        },
    );
    assert_eq!(complete.status(), SingleMappingStatus::Ambiguous);
    assert_eq!(complete.mapping_quality(), 0);
    assert_eq!(complete.best_origin_count, 2);
}

#[test]
fn bounded_affine_rerank_prefers_the_simpler_equal_distance_path() {
    let read = [
        Base::A,
        Base::C,
        Base::G,
        Base::T,
        Base::A,
        Base::C,
        Base::G,
        Base::T,
        Base::A,
        Base::C,
    ];
    let mut reference_bases = vec![Base::A; 30];
    reference_bases[..10].copy_from_slice(&[
        Base::C,
        Base::G,
        Base::T,
        Base::A,
        Base::C,
        Base::G,
        Base::T,
        Base::A,
        Base::C,
        Base::A,
    ]);
    reference_bases[20..].copy_from_slice(&[
        Base::G,
        Base::C,
        Base::G,
        Base::T,
        Base::G,
        Base::C,
        Base::G,
        Base::T,
        Base::A,
        Base::C,
    ]);
    let reference = test_reference(&reference_bases);
    let shifted = ReadPlacement::strict(0, 0, 10, BisulfiteStrand::OT, 2);
    let two_substitutions = ReadPlacement::strict(0, 20, 30, BisulfiteStrand::OT, 2);
    let mut workspace = ReadWorkspace::with_capacity(2, 2);
    workspace.placements.push(shifted);
    workspace.placements.push(two_substitutions);

    let mut origins = Vec::new();
    let result = SingleBatchAligner::finish_result_with_affine_rerank(
        &reference,
        &read,
        &mut workspace,
        &mut origins,
        SingleResultEvidence {
            read_length: read.len(),
            metrics: ReadAlignmentMetrics::default(),
            verified_distance_limit: INITIAL_EDIT_DISTANCE,
            first_seed: None,
            frontier_complete: false,
            confidence_audit_complete: false,
        },
    )
    .expect("bounded fixture reranking succeeds");

    assert_eq!(result.status(), SingleMappingStatus::Ambiguous);
    assert_eq!(result.placement(), Some(two_substitutions));
    assert_eq!(result.mapping_quality(), 0);
    assert_eq!(result.best_origin_count, 2);
    assert_eq!(workspace.affine_score_cache.len(), 2);

    origins.clear();
    let completed = SingleBatchAligner::finish_result_with_affine_rerank(
        &reference,
        &read,
        &mut workspace,
        &mut origins,
        SingleResultEvidence {
            read_length: read.len(),
            metrics: ReadAlignmentMetrics::default(),
            verified_distance_limit: INITIAL_EDIT_DISTANCE,
            first_seed: None,
            frontier_complete: true,
            confidence_audit_complete: true,
        },
    )
    .expect("completed bounded fixture reranking succeeds");

    assert_eq!(completed.status(), SingleMappingStatus::Unique);
    assert_eq!(completed.placement(), Some(two_substitutions));
    assert_eq!(
        completed.mapping_quality(),
        MAPQ_POLICY.single.affine_unique_mapq
    );
    assert_eq!(completed.best_origin_count, 1);
    assert_eq!(workspace.affine_score_cache.len(), 2);

    workspace.begin_verification_cache_read();
    assert!(workspace.affine_scores.is_empty());
    assert!(workspace.affine_score_cache.is_empty());
}

#[test]
fn conversion_counts_distinguish_cg_chg_and_chh_on_both_strands() {
    let top_reference = test_reference(&[
        Base::C,
        Base::G,
        Base::C,
        Base::A,
        Base::G,
        Base::C,
        Base::T,
        Base::C,
        Base::C,
    ]);
    let top_read = [
        Base::T,
        Base::G,
        Base::C,
        Base::A,
        Base::G,
        Base::T,
        Base::T,
        Base::C,
        Base::C,
    ];
    assert_eq!(
        placement_conversion_counts(
            &top_reference,
            &top_read,
            ReadPlacement::strict(0, 0, 9, BisulfiteStrand::OT, 0),
        ),
        Some(([1, 0, 1], [0, 1, 2]))
    );

    let bottom_bases = [
        Base::G,
        Base::C,
        Base::T,
        Base::G,
        Base::A,
        Base::C,
        Base::G,
        Base::G,
    ];
    let bottom_reference = test_reference(&bottom_bases);
    let aligned_query = [
        Base::A,
        Base::C,
        Base::T,
        Base::G,
        Base::A,
        Base::C,
        Base::A,
        Base::G,
    ];
    let bottom_read = aligned_query
        .iter()
        .rev()
        .map(|base| base.complement())
        .collect::<Vec<_>>();
    assert_eq!(
        placement_conversion_counts(
            &bottom_reference,
            &bottom_read,
            ReadPlacement::strict(0, 0, 8, BisulfiteStrand::OB, 0),
        ),
        Some(([1, 0, 1], [0, 2, 0]))
    );
}

#[test]
fn representative_reranking_never_promotes_an_ambiguous_origin() {
    let read = [Base::A; 10];
    let reference = test_reference(&[Base::A; 30]);
    let gapped = ReadPlacement::strict(0, 0, 9, BisulfiteStrand::OT, 2);
    let ungapped = ReadPlacement::strict(0, 20, 30, BisulfiteStrand::OT, 2);
    let mut workspace = ReadWorkspace::with_capacity(2, 2);
    workspace.placements.push(gapped);
    workspace.placements.push(ungapped);
    let mut origins = Vec::new();

    let result = SingleBatchAligner::finish_result_with_affine_rerank(
        &reference,
        &read,
        &mut workspace,
        &mut origins,
        SingleResultEvidence {
            read_length: read.len(),
            metrics: ReadAlignmentMetrics {
                located_rows: 4,
                verified_placements: 2,
                ..ReadAlignmentMetrics::default()
            },
            verified_distance_limit: INITIAL_EDIT_DISTANCE,
            first_seed: None,
            frontier_complete: false,
            confidence_audit_complete: false,
        },
    )
    .expect("bounded fixture reranking succeeds");
    assert_eq!(result.status(), SingleMappingStatus::Ambiguous);
    assert_eq!(result.placement(), Some(ungapped));
    assert_eq!(result.mapping_quality(), 0);
}

#[test]
fn non_directional_merge_selects_the_global_best_pass() {
    let original = mapped_on_strand(100, 30, 2, BisulfiteStrand::OT);
    let complementary = mapped_on_strand(500, 40, 0, BisulfiteStrand::CTOT);
    let merged = merge_non_directional_results(original, complementary);
    assert_eq!(merged.status(), SingleMappingStatus::Unique);
    assert_eq!(merged.placement(), complementary.placement());
    assert_eq!(merged.mapping_quality(), 20);
    assert_eq!(merged.located_rows(), 2);
    assert_eq!(merged.verified_placements(), 2);
}

#[test]
fn non_directional_equal_best_passes_are_ambiguous() {
    let original = mapped_on_strand(100, 40, 0, BisulfiteStrand::OT);
    let complementary = mapped_on_strand(500, 40, 0, BisulfiteStrand::CTOT);
    let merged = merge_non_directional_results(original, complementary);
    assert_eq!(merged.status(), SingleMappingStatus::Ambiguous);
    assert_eq!(merged.placement(), original.placement());
    assert_eq!(merged.mapping_quality(), 0);
}

#[test]
fn non_directional_equal_best_passes_combine_their_origin_counts() {
    let mut original = mapped_on_strand(100, 0, 0, BisulfiteStrand::OT);
    original.status = SingleMappingStatus::Ambiguous;
    original.best_origin_count = 8;
    let mut complementary = mapped_on_strand(500, 0, 0, BisulfiteStrand::CTOT);
    complementary.status = SingleMappingStatus::Ambiguous;
    complementary.best_origin_count = 8;

    let merged = merge_non_directional_results(original, complementary);
    assert_eq!(merged.best_origin_count, 16);
    assert_eq!(merged.status(), SingleMappingStatus::Ambiguous);
    assert_eq!(merged.placement(), original.placement());
    assert_eq!(merged.mapping_quality(), 0);
}

#[test]
fn non_directional_runner_up_caps_mapq_and_combines_repeat_pressure() {
    let original = mapped_on_strand(100, 40, 0, BisulfiteStrand::OT);
    let mut complementary = mapped_on_strand(500, 30, 1, BisulfiteStrand::CTOT);
    complementary.distinct_candidate_starts = 65;
    let merged = merge_non_directional_results(original, complementary);
    assert_eq!(merged.status(), SingleMappingStatus::Unique);
    assert_eq!(merged.mapping_quality(), 10);
}

#[test]
fn completed_q40_evidence_supersedes_merged_row_pressure() {
    let mut original = mapped_on_strand(100, 40, 0, BisulfiteStrand::OT);
    original.located_rows = 257;
    let complementary = SingleAlignmentResult::unmapped(0, 0, 0);
    let merged = merge_non_directional_results(original, complementary);
    assert_eq!(merged.status(), SingleMappingStatus::Unique);
    assert_eq!(merged.mapping_quality(), 40);
}

#[test]
fn complete_non_directional_frontiers_do_not_reapply_row_pressure_proxy() {
    let mut original = mapped_on_strand(100, 30, 0, BisulfiteStrand::OT);
    original.located_rows = 257;
    let complementary = SingleAlignmentResult::unmapped(0, 0, 0);

    let incomplete = merge_non_directional_results(original, complementary);
    assert_eq!(incomplete.mapping_quality(), 10);

    let complete = merge_non_directional_completed_frontiers(original, complementary);
    assert_eq!(complete.status(), SingleMappingStatus::Unique);
    assert_eq!(complete.mapping_quality(), 30);
}

#[test]
fn complementary_pass_relabels_both_search_strands() {
    assert_eq!(
        ConversionPass::Complementary.relabel_combined_hit(BisulfiteStrand::OT),
        Some(BisulfiteStrand::CTOT)
    );
    assert_eq!(
        ConversionPass::Complementary.relabel_combined_hit(BisulfiteStrand::OB),
        Some(BisulfiteStrand::CTOB)
    );
}
