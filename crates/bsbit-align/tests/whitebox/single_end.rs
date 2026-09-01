//! White-box tests for canonical single-end mapping.

use super::*;
use bsbit_core::bisulfite::BisulfiteStrand;

fn mapped(start: u64, mapping_quality: u8) -> SingleAlignmentResult {
    mapped_at_distance(start, mapping_quality, 1)
}

fn mapped_at_distance(start: u64, mapping_quality: u8, edit_distance: u8) -> SingleAlignmentResult {
    SingleAlignmentResult {
        status: SingleMappingStatus::Unique,
        placement: Some(ReadPlacement::strict(
            0,
            start,
            start + 100,
            BisulfiteStrand::OT,
            edit_distance,
        )),
        mapping_quality,
        located_rows: 1,
        verified_placements: 1,
    }
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
    assert_eq!(sensitive.maximum_seed_rounds, default.maximum_seed_rounds);
    assert!(!SingleSearchMode::Default.completes_candidate_frontier());
    assert!(SingleSearchMode::Sensitive.completes_candidate_frontier());
}

#[test]
fn sensitive_low_confidence_conflict_preserves_the_incumbent_as_ambiguous() {
    let incumbent = mapped(100, 30);
    let completed = mapped(500, SENSITIVE_REPLACEMENT_MIN_MAPQ - 1);
    let reconciled = SingleBatchAligner::reconcile_sensitive_result(incumbent, completed, 100);
    assert_eq!(reconciled.status(), SingleMappingStatus::Ambiguous);
    assert_eq!(reconciled.placement(), incumbent.placement());
    assert_eq!(reconciled.mapping_quality(), 0);
}

#[test]
fn sensitive_q20_conflict_may_replace_the_incumbent() {
    let incumbent = mapped(100, 30);
    let completed = mapped(500, SENSITIVE_REPLACEMENT_MIN_MAPQ);
    let reconciled = SingleBatchAligner::reconcile_sensitive_result(incumbent, completed, 100);
    assert_eq!(reconciled, completed);
}

#[test]
fn sensitive_low_confidence_rescue_remains_unmapped() {
    let incumbent = SingleAlignmentResult::unmapped(1, 0);
    let completed = mapped(500, SENSITIVE_REPLACEMENT_MIN_MAPQ - 1);
    let reconciled = SingleBatchAligner::reconcile_sensitive_result(incumbent, completed, 100);
    assert_eq!(reconciled.status(), SingleMappingStatus::Unmapped);
    assert_eq!(reconciled.placement(), None);
}

#[test]
fn sensitive_strong_incumbent_uses_the_d3_confidence_boundary() {
    assert_eq!(
        SingleBatchAligner::sensitive_audit_distance_limit(mapped(100, 20), 5),
        3
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
            SingleAlignmentResult::unmapped(0, 0),
            5,
        ),
        5
    );
}
