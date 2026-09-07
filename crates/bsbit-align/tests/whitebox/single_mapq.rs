//! White-box tests for single-end MAPQ evidence gates.

use super::*;

fn baseline_single_mapq_evidence() -> SingleMapqEvidence {
    SingleMapqEvidence {
        read_length: 100,
        best_distance: 0,
        second_best_distance: None,
        verified_distance_limit: 3,
        located_rows: 1,
        distinct_candidate_starts: 1,
        verified_placements: 1,
        first_seed_hits: 1,
        first_seed_bases: 80,
        direct_singleton: false,
        frontier_complete: false,
        confidence_audit_complete: false,
        seed_round_support: 1,
    }
}

#[test]
fn single_mapq_uses_whole_read_edit_separation() {
    let mut evidence = baseline_single_mapq_evidence();
    evidence.frontier_complete = true;
    evidence.confidence_audit_complete = true;
    evidence.seed_round_support = 3;
    assert_eq!(single_mapping_quality_from_evidence(evidence), 40);

    evidence.verified_distance_limit = 5;
    assert_eq!(single_mapping_quality_from_evidence(evidence), 60);
}

#[test]
fn single_mapq_adverse_evidence_can_only_lower_confidence() {
    let mut evidence = baseline_single_mapq_evidence();
    evidence.second_best_distance = Some(1);
    assert_eq!(single_mapping_quality_from_evidence(evidence), 10);

    evidence = baseline_single_mapq_evidence();
    evidence.located_rows = 257;
    assert_eq!(single_mapping_quality_from_evidence(evidence), 10);

    evidence = baseline_single_mapq_evidence();
    evidence.best_distance = 5;
    evidence.verified_distance_limit = 5;
    assert_eq!(single_mapping_quality_from_evidence(evidence), 10);
}

#[test]
fn single_mapq_uses_the_verified_boundary_without_an_observed_runner_up() {
    let mut evidence = baseline_single_mapq_evidence();
    evidence.best_distance = 3;
    evidence.verified_distance_limit = 3;
    assert_eq!(single_mapping_quality_from_evidence(evidence), 10);

    evidence.best_distance = 2;
    assert_eq!(single_mapping_quality_from_evidence(evidence), 20);
}

#[test]
fn incomplete_frontier_needs_an_offset_seed_for_q30() {
    let mut evidence = baseline_single_mapq_evidence();
    evidence.verified_distance_limit = 5;
    evidence.seed_round_support = 0;
    assert_eq!(single_mapping_quality_from_evidence(evidence), 20);

    evidence.seed_round_support = 1;
    assert_eq!(single_mapping_quality_from_evidence(evidence), 30);
}

#[test]
fn completed_default_confidence_audit_can_reach_q40_without_full_frontier() {
    let mut evidence = baseline_single_mapq_evidence();
    evidence.seed_round_support = MAPQ_POLICY.single.strong_seed_rounds;

    assert_eq!(single_mapping_quality_from_evidence(evidence), 30);

    evidence.confidence_audit_complete = true;
    assert!(!evidence.frontier_complete);
    assert_eq!(single_mapping_quality_from_evidence(evidence), 40);

    evidence.located_rows = MAPQ_POLICY.single.repeat_located_rows + 1;
    assert_eq!(single_mapping_quality_from_evidence(evidence), 10);
}

#[test]
fn completed_frontier_uses_three_distinct_seed_rounds_for_q40() {
    let mut evidence = baseline_single_mapq_evidence();
    evidence.best_distance = 2;
    evidence.verified_distance_limit = 3;
    evidence.frontier_complete = true;
    evidence.confidence_audit_complete = true;
    evidence.seed_round_support = 2;
    assert_eq!(single_mapping_quality_from_evidence(evidence), 20);

    evidence.located_rows = 257;
    assert_eq!(single_mapping_quality_from_evidence(evidence), 20);

    evidence.best_distance = 0;
    evidence.seed_round_support = 3;
    assert_eq!(single_mapping_quality_from_evidence(evidence), 40);

    evidence.verified_distance_limit = 4;
    evidence.seed_round_support = 2;
    assert_eq!(single_mapping_quality_from_evidence(evidence), 50);

    evidence.best_distance = 1;
    evidence.verified_distance_limit = 4;
    evidence.seed_round_support = 2;
    assert_eq!(single_mapping_quality_from_evidence(evidence), 40);

    evidence.verified_distance_limit = 5;
    assert_eq!(single_mapping_quality_from_evidence(evidence), 50);

    evidence.seed_round_support = 3;
    assert_eq!(single_mapping_quality_from_evidence(evidence), 50);

    evidence.best_distance = 2;
    evidence.verified_distance_limit = 5;
    evidence.seed_round_support = 3;
    assert_eq!(single_mapping_quality_from_evidence(evidence), 40);

    evidence.best_distance = 3;
    evidence.verified_distance_limit = 5;
    assert_eq!(single_mapping_quality_from_evidence(evidence), 40);

    evidence.best_distance = 1;
    evidence.verified_distance_limit = 4;
    evidence.seed_round_support = 0;
    evidence.direct_singleton = true;
    evidence.first_seed_bases = 32;
    assert_eq!(single_mapping_quality_from_evidence(evidence), 40);
}

#[test]
fn completed_frontier_applies_multiseed_evidence_consistently_at_q20() {
    let mut evidence = baseline_single_mapq_evidence();
    evidence.best_distance = 3;
    evidence.second_best_distance = Some(4);
    evidence.verified_distance_limit = 5;
    evidence.frontier_complete = true;
    evidence.confidence_audit_complete = true;
    evidence.seed_round_support = 3;
    assert_eq!(single_mapping_quality_from_evidence(evidence), 20);

    evidence.seed_round_support = 2;
    assert_eq!(single_mapping_quality_from_evidence(evidence), 10);
}

#[test]
fn complete_runner_up_gap_replaces_incomplete_row_pressure_proxy() {
    let mut evidence = baseline_single_mapq_evidence();
    evidence.second_best_distance = Some(2);
    evidence.located_rows = 257;
    assert_eq!(single_mapping_quality_from_evidence(evidence), 10);

    evidence.frontier_complete = true;
    evidence.confidence_audit_complete = true;
    assert_eq!(single_mapping_quality_from_evidence(evidence), 20);
}

#[test]
fn completed_d4_multiseed_frontier_can_reach_q20() {
    let mut evidence = baseline_single_mapq_evidence();
    evidence.best_distance = 4;
    evidence.second_best_distance = None;
    evidence.verified_distance_limit = 5;
    evidence.seed_round_support = 2;
    assert_eq!(single_mapping_quality_from_evidence(evidence), 15);

    evidence.frontier_complete = true;
    evidence.confidence_audit_complete = true;
    assert_eq!(single_mapping_quality_from_evidence(evidence), 20);
}

#[test]
fn completed_coordinate_certificate_adapts_to_read_length() {
    let mut evidence = baseline_single_mapq_evidence();
    evidence.best_distance = 2;
    evidence.second_best_distance = Some(3);
    evidence.frontier_complete = true;
    evidence.confidence_audit_complete = true;
    evidence.direct_singleton = true;
    evidence.seed_round_support = 2;
    assert_eq!(single_mapping_quality_from_evidence(evidence), 10);

    evidence.read_length = 150;
    assert_eq!(
        completed_coordinate_certificate(evidence),
        Some(MapqCertificate::LongReadSingletonCorroboration)
    );
    assert_eq!(single_mapping_quality_from_evidence(evidence), 20);

    evidence.read_length = 100;
    evidence.direct_singleton = false;
    evidence.seed_round_support = 3;
    assert_eq!(
        completed_coordinate_certificate(evidence),
        Some(MapqCertificate::CompletedMultiSeedCoordinate)
    );
    assert_eq!(single_mapping_quality_from_evidence(evidence), 20);

    evidence.read_length = 150;
    evidence.best_distance = 5;
    evidence.verified_distance_limit = 5;
    evidence.second_best_distance = None;
    evidence.seed_round_support = 2;
    assert_eq!(
        completed_coordinate_certificate(evidence),
        Some(MapqCertificate::LongReadBoundaryCorroboration)
    );
    assert_eq!(single_mapping_quality_from_evidence(evidence), 20);
}

#[test]
fn completed_local_locus_certificate_is_narrower_than_a_global_floor() {
    let mut evidence = baseline_single_mapq_evidence();
    evidence.best_distance = 1;
    evidence.second_best_distance = Some(2);
    evidence.frontier_complete = true;
    evidence.confidence_audit_complete = true;
    evidence.seed_round_support = 2;

    assert_eq!(completed_coordinate_certificate(evidence), None);
    assert_eq!(
        completed_local_locus_certificate(evidence),
        Some(MapqCertificate::LowEditLocalLocusCorroboration)
    );
    assert_eq!(single_mapping_quality_from_evidence(evidence), 10);

    evidence.best_distance = 2;
    assert_eq!(completed_local_locus_certificate(evidence), None);
}

#[test]
fn single_mapq_layers_preserve_the_frozen_decision_order() {
    let mut evidence = baseline_single_mapq_evidence();
    evidence.frontier_complete = true;
    evidence.confidence_audit_complete = true;
    evidence.seed_round_support = MAPQ_POLICY.single.strong_seed_rounds;

    let raw = raw_separation_mapq(evidence);
    let capped = apply_adverse_caps(raw, evidence);
    let certified = apply_positive_certificates(capped, raw, evidence);

    assert_eq!(raw, 40);
    assert_eq!(capped, 40);
    assert_eq!(certified, 40);
}

fn frozen_structural_origin_evidence_v1(evidence: SingleMapqEvidence) -> u8 {
    let separation = evidence.second_best_distance.map_or_else(
        || {
            evidence
                .verified_distance_limit
                .saturating_add(1)
                .saturating_sub(evidence.best_distance)
        },
        |second| second.saturating_sub(evidence.best_distance),
    );
    let mut mapq = separation.saturating_mul(10).min(60);
    let short_singleton = evidence.direct_singleton
        && evidence.first_seed_hits == 1
        && (16..=46).contains(&evidence.first_seed_bases);
    let completed_multiseed_frontier = evidence.frontier_complete
        && evidence.seed_round_support >= 2
        && evidence.best_distance <= 3;
    let strong_multiseed_certificate = completed_multiseed_frontier
        && (evidence.seed_round_support >= 3
            || evidence.seed_round_support >= 2
                && (mapq >= 50 || mapq >= 40 && evidence.best_distance <= 1));
    let high_confidence_certificate = evidence.frontier_complete
        && evidence.best_distance <= 3
        && (short_singleton || strong_multiseed_certificate);
    if matches!(mapq, 10 | 30) && completed_multiseed_frontier && evidence.seed_round_support >= 3 {
        mapq = mapq.saturating_add(10);
    }
    if mapq >= 40 && !high_confidence_certificate {
        mapq = mapq.min(30);
    }
    let moderate_edit_multiseed = evidence.frontier_complete
        && evidence.best_distance == 4
        && evidence.seed_round_support >= 2;
    if evidence.best_distance >= 5 {
        mapq = mapq.min(10);
    } else if evidence.best_distance == 4 {
        mapq = mapq.min(if moderate_edit_multiseed { 20 } else { 15 });
    }
    let repeat_risk = evidence.first_seed_hits > 64
        || evidence.located_rows > 256
        || evidence.distinct_candidate_starts > 64
        || evidence.verified_placements > 64;
    if repeat_risk && !evidence.frontier_complete {
        mapq = mapq.min(10);
    }
    if !evidence.frontier_complete && evidence.seed_round_support < 1 {
        mapq = mapq.min(20);
    }
    let coordinate_certificate = evidence.frontier_complete
        && (evidence.seed_round_support >= 3
            || evidence.read_length >= 128
                && evidence.direct_singleton
                && evidence.seed_round_support >= 2
                && evidence.best_distance <= 2
            || evidence.read_length >= 128
                && evidence.seed_round_support >= 2
                && evidence.best_distance == evidence.verified_distance_limit
                && evidence.second_best_distance.is_none());
    if coordinate_certificate {
        mapq = mapq.max(20);
    }
    mapq
}

#[test]
fn refactored_layers_match_the_frozen_v1_boundary_oracle() {
    const READ_LENGTHS: [usize; 4] = [100, 127, 128, 192];
    const SECOND_BEST: [Option<u8>; 8] = [
        None,
        Some(0),
        Some(1),
        Some(2),
        Some(3),
        Some(4),
        Some(5),
        Some(6),
    ];
    const COUNTS: [u64; 3] = [1, 64, 65];
    const LOCATED_ROWS: [u64; 3] = [1, 256, 257];
    const SEED_BASES: [u64; 4] = [15, 16, 46, 47];

    for read_length in READ_LENGTHS {
        for best_distance in 0..=5 {
            for second_best_distance in SECOND_BEST {
                for verified_distance_limit in [3, 5] {
                    for located_rows in LOCATED_ROWS {
                        for distinct_candidate_starts in COUNTS {
                            for verified_placements in COUNTS {
                                for first_seed_hits in COUNTS {
                                    for first_seed_bases in SEED_BASES {
                                        for direct_singleton in [false, true] {
                                            for frontier_complete in [false, true] {
                                                for seed_round_support in 0..=3 {
                                                    let evidence = SingleMapqEvidence {
                                                        read_length,
                                                        best_distance,
                                                        second_best_distance,
                                                        verified_distance_limit,
                                                        located_rows,
                                                        distinct_candidate_starts,
                                                        verified_placements,
                                                        first_seed_hits,
                                                        first_seed_bases,
                                                        direct_singleton,
                                                        frontier_complete,
                                                        confidence_audit_complete:
                                                            frontier_complete,
                                                        seed_round_support,
                                                    };
                                                    assert_eq!(
                                                        single_mapping_quality_from_evidence(
                                                            evidence
                                                        ),
                                                        frozen_structural_origin_evidence_v1(
                                                            evidence
                                                        ),
                                                        "evidence={evidence:?}"
                                                    );
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
