//! Public mapping-quality evidence contract.

use bsbit_align::paired_end::{PairMappingStatus, bwa_pair_mapping_quality_from_evidence};

#[test]
fn unique_complete_frontier_uses_score_gap() {
    assert_eq!(
        bwa_pair_mapping_quality_from_evidence(
            PairMappingStatus::Unique,
            true,
            Some(100),
            Some(90),
            0,
        ),
        60,
    );
}

#[test]
fn incomplete_or_ambiguous_frontier_has_zero_quality() {
    assert_eq!(
        bwa_pair_mapping_quality_from_evidence(
            PairMappingStatus::Unique,
            false,
            Some(100),
            Some(90),
            0,
        ),
        0,
    );
    assert_eq!(
        bwa_pair_mapping_quality_from_evidence(
            PairMappingStatus::Ambiguous,
            true,
            Some(100),
            Some(90),
            0,
        ),
        0,
    );
}

#[test]
fn near_best_pairings_reduce_quality() {
    let without_repeats = bwa_pair_mapping_quality_from_evidence(
        PairMappingStatus::Unique,
        true,
        Some(100),
        Some(95),
        0,
    );
    let with_repeats = bwa_pair_mapping_quality_from_evidence(
        PairMappingStatus::Unique,
        true,
        Some(100),
        Some(95),
        3,
    );
    assert!(with_repeats < without_repeats);
}
