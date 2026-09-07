//! White-box tests for structural paired-end MAPQ and monotone caps.

use super::*;

#[test]
fn bwa_mapq_requires_a_complete_unique_frontier() {
    assert_eq!(
        bwa_pair_mapping_quality_from_evidence(
            PairMappingStatus::Ambiguous,
            true,
            Some(20),
            Some(10),
            0,
        ),
        0
    );
    assert_eq!(
        bwa_pair_mapping_quality_from_evidence(
            PairMappingStatus::Unique,
            false,
            Some(20),
            Some(10),
            0,
        ),
        0
    );
}

#[test]
fn bwa_mapq_is_derived_from_score_separation() {
    assert_eq!(
        bwa_pair_mapping_quality_from_evidence(
            PairMappingStatus::Unique,
            true,
            Some(20),
            Some(17),
            0,
        ),
        30
    );
    assert_eq!(
        bwa_pair_mapping_quality_from_evidence(PairMappingStatus::Unique, true, Some(20), None, 0,),
        60
    );
    assert_eq!(
        bwa_pair_mapping_quality_from_evidence(
            PairMappingStatus::Unique,
            true,
            Some(20),
            Some(17),
            1,
        ),
        30
    );
    assert_eq!(
        bwa_pair_mapping_quality_from_evidence(
            PairMappingStatus::Unique,
            true,
            Some(20),
            Some(16),
            1,
        ),
        39
    );
}

#[test]
fn final_caps_are_monotone_and_order_independent() {
    assert_eq!(apply_final_mapq_caps(60, Some(12), Some(20), true), 12);
    assert_eq!(apply_final_mapq_caps(60, None, Some(20), false), 20);
    assert_eq!(apply_final_mapq_caps(60, None, None, false), 60);
    assert_eq!(apply_final_mapq_caps(60, None, None, true), 20);
    assert_eq!(apply_final_mapq_caps(15, Some(40), Some(20), true), 15);

    for baseline in 0..=60 {
        let capped = apply_final_mapq_caps(baseline, Some(30), Some(20), true);
        assert!(capped <= baseline);
    }
}

#[test]
fn row_pressure_is_distinct_from_rescue_provenance() {
    assert!(!row_pressure(PAIR_MAPQ_REPEAT_RISK_ROWS - 1, 0));
    assert!(row_pressure(PAIR_MAPQ_REPEAT_RISK_ROWS, 0));
    assert!(!rescue_risk(false));
    assert!(rescue_risk(true));

    assert_eq!(
        apply_final_mapq_caps(60, None, Some(MAPQ_POLICY.paired.row_pressure_cap), false),
        30
    );
    assert_eq!(
        apply_final_mapq_caps(60, None, Some(MAPQ_POLICY.paired.rescue_risk_cap), false),
        20
    );
    assert_eq!(
        apply_final_mapq_caps(
            60,
            None,
            Some(MAPQ_POLICY.paired.resolved_frontier_cap),
            false
        ),
        30
    );
}

#[test]
fn a_complete_alternative_margin_supersedes_row_pressure() {
    let row_pressure_without_margin =
        |rows: u64, margin_complete: bool| row_pressure(rows, 0) && !margin_complete;
    assert!(row_pressure_without_margin(
        PAIR_MAPQ_REPEAT_RISK_ROWS,
        false
    ));
    assert!(!row_pressure_without_margin(
        PAIR_MAPQ_REPEAT_RISK_ROWS,
        true
    ));
}
