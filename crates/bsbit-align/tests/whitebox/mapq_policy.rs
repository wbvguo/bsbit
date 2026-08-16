//! White-box tests for the fixed implementation MAPQ certificates.

use super::*;

fn baseline_sensitive_mapq_evidence() -> SensitiveMapqEvidence {
    SensitiveMapqEvidence {
        baseline_mapq: 0,
        raw_mapq: 0,
        reported_ambiguous: true,
        frontier_complete: false,
        best_pair_placements: 1,
        compatible_pairs: 1,
        best_score: 0,
        second_best_present: false,
        score_gap: 99,
        near_best_pairings: 0,
        located_rows_min: 1,
        located_rows_max: 2,
        located_rows_sum: 3,
        emitted_candidate_starts_sum: 3,
        distinct_candidate_starts_sum: 3,
        verified_placements_sum: 3,
        pair_distance: 4,
        pair_score: 4,
        mate_distance_max: 2,
        net_gap_sum: 0,
        clipped_bases_sum: 0,
    }
}

#[test]
fn q10_uses_bounded_pair_evidence() {
    let accepted = baseline_sensitive_mapq_evidence();
    assert!(sensitive_q10_certified(accepted));
    assert_eq!(apply_sensitive_mapq_policy(accepted), 10);

    let mut too_many_best_placements = accepted;
    too_many_best_placements.best_pair_placements = 3;
    assert!(!sensitive_q10_certified(too_many_best_placements));
    assert_eq!(apply_sensitive_mapq_policy(too_many_best_placements), 0);

    let mut unreported_ambiguity = accepted;
    unreported_ambiguity.reported_ambiguous = false;
    assert!(!sensitive_q10_certified(unreported_ambiguity));
}

#[test]
fn q30_accepts_complete_low_mapq_evidence() {
    let mut evidence = baseline_sensitive_mapq_evidence();
    evidence.baseline_mapq = 19;
    evidence.raw_mapq = 9;
    evidence.reported_ambiguous = false;
    evidence.frontier_complete = true;
    evidence.pair_distance = 3;
    evidence.located_rows_max = 400;
    evidence.located_rows_sum = 500;
    evidence.emitted_candidate_starts_sum = 6;
    evidence.compatible_pairs = 100;

    assert!(sensitive_q30_low_certified(evidence));
    assert_eq!(apply_sensitive_mapq_policy(evidence), 30);
}

#[test]
fn q30_promotes_the_qualified_q20_bin() {
    let mut evidence = baseline_sensitive_mapq_evidence();
    evidence.baseline_mapq = 25;
    evidence.raw_mapq = 25;
    evidence.reported_ambiguous = false;
    assert_eq!(apply_sensitive_mapq_policy(evidence), 30);
}

#[test]
fn q40_requires_complete_strong_evidence() {
    let mut evidence = baseline_sensitive_mapq_evidence();
    evidence.baseline_mapq = 20;
    evidence.raw_mapq = 25;
    evidence.reported_ambiguous = false;
    evidence.frontier_complete = true;
    evidence.compatible_pairs = 2;

    assert!(sensitive_q40_certified(evidence));
    assert_eq!(apply_sensitive_mapq_policy(evidence), 40);

    evidence.frontier_complete = false;
    assert!(!sensitive_q40_certified(evidence));
    assert_eq!(apply_sensitive_mapq_policy(evidence), 30);
}

#[test]
fn origin_grouping_is_monotone_and_requires_boundary_evidence() {
    let mut grouped = baseline_sensitive_mapq_evidence();
    grouped.pair_distance = 2;
    assert_eq!(apply_origin_grouped_mapq_policy(0, 10, grouped), 0);
    grouped.frontier_complete = true;
    assert_eq!(apply_origin_grouped_mapq_policy(0, 10, grouped), 10);
    grouped.pair_distance = 3;
    assert_eq!(apply_origin_grouped_mapq_policy(0, 10, grouped), 0);

    grouped.score_gap = 16;
    assert_eq!(apply_origin_grouped_mapq_policy(19, 30, grouped), 30);
    grouped.score_gap = 15;
    assert_eq!(apply_origin_grouped_mapq_policy(19, 30, grouped), 19);

    grouped.pair_distance = 0;
    assert_eq!(apply_origin_grouped_mapq_policy(30, 40, grouped), 40);
    grouped.pair_distance = 4;
    grouped.clipped_bases_sum = 0;
    grouped.second_best_present = false;
    assert_eq!(apply_origin_grouped_mapq_policy(30, 40, grouped), 30);
    assert_eq!(apply_origin_grouped_mapq_policy(40, 30, grouped), 40);
}

#[test]
fn q40_certificate_only_demotes_diffuse_high_pressure_evidence() {
    let mut evidence = baseline_sensitive_mapq_evidence();
    evidence.reported_ambiguous = false;
    evidence.frontier_complete = true;
    evidence.distinct_candidate_starts_sum = 0;
    evidence.located_rows_max = 512;
    assert_eq!(apply_sensitive_q40_certificate(40, evidence), 40);

    evidence.distinct_candidate_starts_sum = 2;
    evidence.located_rows_max = 15;
    assert_eq!(apply_sensitive_q40_certificate(60, evidence), 60);

    evidence.distinct_candidate_starts_sum = 1;
    evidence.located_rows_max = 127;
    evidence.pair_distance = 3;
    assert_eq!(apply_sensitive_q40_certificate(40, evidence), 40);

    evidence.pair_distance = 4;
    assert_eq!(apply_sensitive_q40_certificate(40, evidence), 30);
    assert_eq!(apply_sensitive_q40_certificate(30, evidence), 30);

    evidence.pair_distance = 3;
    evidence.frontier_complete = false;
    assert_eq!(apply_sensitive_q40_certificate(40, evidence), 30);
}
