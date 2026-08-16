//! White-box tests for single-end MAPQ evidence gates.

use super::*;

fn baseline_single_mapq_evidence() -> SingleMapqEvidence {
    SingleMapqEvidence {
        best_distance: 0,
        second_best_distance: None,
        verified_distance_limit: 3,
        located_rows: 1,
        distinct_candidate_starts: 1,
        verified_placements: 1,
        first_seed_hits: 1,
        first_seed_bases: 30,
        direct_singleton: true,
    }
}

#[test]
fn single_mapq_uses_early_interval_uniqueness_in_the_right_direction() {
    let mut evidence = baseline_single_mapq_evidence();
    assert_eq!(single_mapping_quality_from_evidence(evidence), 40);

    evidence.first_seed_bases = 42;
    assert_eq!(single_mapping_quality_from_evidence(evidence), 30);

    evidence.first_seed_bases = 47;
    assert_eq!(single_mapping_quality_from_evidence(evidence), 20);
}

#[test]
fn single_mapq_reserves_q40_for_an_exact_short_singleton() {
    let mut evidence = baseline_single_mapq_evidence();
    evidence.best_distance = 1;
    assert_eq!(single_mapping_quality_from_evidence(evidence), 30);

    evidence.direct_singleton = false;
    assert_eq!(single_mapping_quality_from_evidence(evidence), 20);
}

#[test]
fn single_mapq_adverse_evidence_can_only_lower_confidence() {
    let mut evidence = baseline_single_mapq_evidence();
    evidence.second_best_distance = Some(1);
    assert_eq!(single_mapping_quality_from_evidence(evidence), 15);

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
    assert_eq!(single_mapping_quality_from_evidence(evidence), 15);

    evidence.best_distance = 2;
    assert_eq!(single_mapping_quality_from_evidence(evidence), 30);
}
