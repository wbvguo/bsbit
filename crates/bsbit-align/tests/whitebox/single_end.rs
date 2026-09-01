//! White-box tests for canonical single-end mapping.

use super::*;

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
