//! White-box tests for the reference implementation.
//!
//! Kept outside production `src/` while remaining a child module so private
//! invariants can be tested without widening the crate API.

use super::*;
use crate::reference::validate_reference_catalog;
#[cfg(feature = "combined-index")]
use crate::storage::fm::ProjectedBase;
use bsbit_core::sequence::normalize_dna;

#[cfg(feature = "combined-index")]
struct LengthOnlyCombined(u64);

#[cfg(feature = "combined-index")]
impl PrivateCombinedIndex for LengthOnlyCombined {
    fn reference_len(&self) -> u64 {
        self.0
    }

    fn exact_interval(
        &self,
        _reversed_projected_pattern: &[SearchBase],
    ) -> Result<Option<FmInterval>, CombinedIndexBackendError> {
        Ok(None)
    }

    fn backward_extend_interval(
        &self,
        _interval: FmInterval,
        _symbol: SearchBase,
    ) -> Result<FmInterval, CombinedIndexBackendError> {
        Err(CombinedIndexBackendError::Structure)
    }

    fn backward_extend_projected_interval(
        &self,
        _interval: FmInterval,
        _symbol: ProjectedBase,
    ) -> Result<FmInterval, CombinedIndexBackendError> {
        Err(CombinedIndexBackendError::Structure)
    }

    fn visit_interval(
        &self,
        _interval: FmInterval,
        _visitor: &mut dyn FnMut(u64) -> bool,
    ) -> Result<PrivateCombinedLocateMetrics, CombinedIndexBackendError> {
        Ok(PrivateCombinedLocateMetrics::new(0, 0, 0, 0))
    }
}

#[cfg(feature = "combined-index")]
struct DiagnosticCombined;

#[cfg(feature = "combined-index")]
impl PrivateCombinedIndex for DiagnosticCombined {
    fn reference_len(&self) -> u64 {
        64
    }

    fn exact_interval(
        &self,
        _reversed_projected_pattern: &[SearchBase],
    ) -> Result<Option<FmInterval>, CombinedIndexBackendError> {
        Ok(None)
    }

    fn resolve_projected_suffix_intervals(
        &self,
        patterns: &[&[ProjectedBase]],
        _minimum_suffix_bases: usize,
        _stop_interval_length: u64,
        output: &mut [Option<(FmInterval, u64)>],
    ) -> Result<(), CombinedIndexBackendError> {
        if patterns.len() != 2 || output.len() != 2 {
            return Err(CombinedIndexBackendError::Structure);
        }
        output[0] = Some((
            FmInterval::private_checked(10, 14, 129)
                .map_err(|_| CombinedIndexBackendError::Structure)?,
            18,
        ));
        output[1] = Some((
            FmInterval::private_checked(20, 21, 129)
                .map_err(|_| CombinedIndexBackendError::Structure)?,
            20,
        ));
        Ok(())
    }

    fn backward_extend_interval(
        &self,
        _interval: FmInterval,
        _symbol: SearchBase,
    ) -> Result<FmInterval, CombinedIndexBackendError> {
        Err(CombinedIndexBackendError::Structure)
    }

    fn backward_extend_projected_interval(
        &self,
        _interval: FmInterval,
        _symbol: ProjectedBase,
    ) -> Result<FmInterval, CombinedIndexBackendError> {
        Err(CombinedIndexBackendError::Structure)
    }

    fn visit_interval(
        &self,
        interval: FmInterval,
        _visitor: &mut dyn FnMut(u64) -> bool,
    ) -> Result<PrivateCombinedLocateMetrics, CombinedIndexBackendError> {
        Ok(PrivateCombinedLocateMetrics::new(interval.len(), 0, 0, 0))
    }

    fn visit_intervals_two_lanes_complete(
        &self,
        intervals: [FmInterval; 2],
        _visitor: &mut dyn FnMut(usize, u64),
    ) -> Result<[PrivateCombinedLocateMetrics; 2], CombinedIndexBackendError> {
        Ok([
            PrivateCombinedLocateMetrics::new(intervals[0].len(), 7, 9, 2),
            PrivateCombinedLocateMetrics::new(intervals[1].len(), 11, 13, 3),
        ])
    }
}

fn catalog_fixture() -> Vec<ContigInput> {
    vec![
        ContigInput::new(
            b"fragmented".to_vec(),
            normalize_dna(b"NACGTNN").expect("valid fixture"),
        ),
        ContigInput::new(
            b"all-n".to_vec(),
            normalize_dna(b"NNN").expect("valid fixture"),
        ),
    ]
}

#[cfg(feature = "combined-index")]
#[test]
fn combined_owner_validates_length_and_retains_reference_metrics() {
    let sequence = normalize_dna(b"NCGTNNGTACGTACN").unwrap();
    let reference = ReferenceIndex::from_private_combined(
        vec![ContigInput::new(b"chr".to_vec(), sequence)],
        PrivateCombinedReference::new(Box::new(LengthOnlyCombined(15))),
    )
    .unwrap();
    assert_eq!(reference.metrics().total_reference_bases(), 15);
    assert!(reference.combined_index_query().is_some());

    let error = ReferenceIndex::from_private_combined(
        vec![ContigInput::new(
            b"chr".to_vec(),
            normalize_dna(b"ACGT").unwrap(),
        )],
        PrivateCombinedReference::new(Box::new(LengthOnlyCombined(3))),
    )
    .expect_err("combined index length must match the reference");
    assert_eq!(
        error,
        PrivateCombinedReferenceError::CombinedDimensions {
            expected_reference_len: 4,
            observed_reference_len: 3,
        }
    );
}

#[cfg(feature = "combined-index")]
#[test]
fn optional_query_diagnostics_count_suffix_rank_and_complete_locate_work() {
    let reference = ReferenceIndex::from_private_combined(
        vec![ContigInput::new(
            b"chr".to_vec(),
            normalize_dna(&[b'A'; 64]).unwrap(),
        )],
        PrivateCombinedReference::new(Box::new(DiagnosticCombined)),
    )
    .unwrap();
    reference.enable_query_diagnostics();
    let query = reference.combined_index_query().unwrap();
    let first = [ProjectedBase::A; 20];
    let second = [ProjectedBase::G; 20];
    let mut intervals = [None; 2];
    query
        .resolve_projected_suffix_intervals(&[&first, &second], 16, 1, &mut intervals)
        .unwrap();
    let intervals = intervals.map(|interval| interval.unwrap().0);
    query
        .visit_raw_intervals_two_lanes_complete(intervals, &mut |_, _| {})
        .unwrap();

    let diagnostics = reference.disable_and_take_query_diagnostics();
    assert_eq!(diagnostics.suffix_search_lanes(), 2);
    assert_eq!(diagnostics.suffix_search_rank_operations(), 14);
    assert_eq!(diagnostics.locate_calls(), 2);
    assert_eq!(diagnostics.singleton_locate_calls(), 1);
    assert_eq!(diagnostics.multi_hit_locate_calls(), 1);
    assert_eq!(diagnostics.located_rows(), 5);
    assert_eq!(diagnostics.locate_lf_steps(), 18);
    assert_eq!(diagnostics.locate_rank_operations(), 22);
    assert_eq!(diagnostics.locate_interval_nodes(), 5);

    reference.enable_query_diagnostics();
    assert_eq!(
        reference.disable_and_take_query_diagnostics(),
        ReferenceQueryDiagnostics::default()
    );
}

#[test]
fn public_catalog_validation_reports_only_catalog_dimensions() {
    let contigs = catalog_fixture();
    let metrics = validate_reference_catalog(&contigs, ReferenceCatalogLimits::MAX)
        .expect("catalog validates without FM construction");
    assert_eq!(metrics.contig_count(), 2);
    assert_eq!(metrics.total_name_bytes(), 15);
    assert_eq!(metrics.total_reference_bases(), 10);

    let build_error = ReferenceIndex::build(
        contigs.clone(),
        ReferenceBuildLimits::MAX.with_max_canonical_runs(0),
    )
    .expect_err("FM-specific build limit rejects the same catalog");
    assert_eq!(
        build_error,
        ReferenceBuildError::LimitExceeded {
            resource: ReferenceResource::CanonicalRuns,
            requested: 1,
            maximum: 0,
        }
    );
    assert!(validate_reference_catalog(&contigs, ReferenceCatalogLimits::MAX).is_ok());
}

#[test]
fn public_catalog_limits_accept_exact_values_and_reject_the_next_catalog() {
    let contigs = catalog_fixture();
    let exact = ReferenceCatalogLimits::MAX
        .with_max_contigs(2)
        .with_max_total_name_bytes(15)
        .with_max_total_reference_bases(10);
    assert!(validate_reference_catalog(&contigs, exact).is_ok());

    let cases = [
        (exact.with_max_contigs(1), ReferenceResource::Contigs, 2, 1),
        (
            exact.with_max_total_name_bytes(14),
            ReferenceResource::TotalNameBytes,
            15,
            14,
        ),
        (
            exact.with_max_total_reference_bases(9),
            ReferenceResource::TotalReferenceBases,
            10,
            9,
        ),
    ];
    for (limits, resource, requested, maximum) in cases {
        assert_eq!(
            validate_reference_catalog(&contigs, limits)
                .expect_err("next catalog exceeds one exact limit"),
            ReferenceBuildError::LimitExceeded {
                resource,
                requested,
                maximum,
            }
        );
    }
}

#[test]
fn public_catalog_validation_matches_build_catalog_error_priority() {
    let invalid_catalogs = [
        Vec::new(),
        vec![ContigInput::new(
            Vec::new(),
            normalize_dna(b"A").expect("valid fixture"),
        )],
        vec![
            ContigInput::new(
                b"duplicate".to_vec(),
                normalize_dna(b"A").expect("valid fixture"),
            ),
            ContigInput::new(
                b"duplicate".to_vec(),
                normalize_dna(b"").expect("empty normalized fixture"),
            ),
        ],
        vec![ContigInput::new(
            b"empty-sequence".to_vec(),
            normalize_dna(b"").expect("empty normalized fixture"),
        )],
    ];
    for contigs in invalid_catalogs {
        let catalog_error = validate_reference_catalog(&contigs, ReferenceCatalogLimits::MAX)
            .expect_err("catalog is invalid");
        let build_error = ReferenceIndex::build(contigs, ReferenceBuildLimits::MAX)
            .expect_err("build rejects the same catalog prefix");
        assert_eq!(catalog_error, build_error);
    }
}

#[test]
fn raw_and_reference_lane_conversions_are_the_accepted_four_lane_table() {
    let expected = [
        (BisulfiteStrand::OT, ThreeLetterConversion::CToT),
        (BisulfiteStrand::OB, ThreeLetterConversion::CToT),
        (BisulfiteStrand::CTOT, ThreeLetterConversion::GToA),
        (BisulfiteStrand::CTOB, ThreeLetterConversion::GToA),
    ];
    for (strand, conversion) in expected {
        assert_eq!(raw_view_conversion(strand), conversion);
        assert_eq!(lane_reference_conversion(strand), conversion);
    }
    assert_eq!(
        strand_semantics(BisulfiteStrand::OB).search_conversion(),
        ThreeLetterConversion::GToA
    );
    assert_eq!(
        strand_semantics(BisulfiteStrand::CTOT).search_conversion(),
        ThreeLetterConversion::CToT
    );
}

#[test]
fn canonical_run_iterator_is_maximal_ordered_and_all_n_safe() {
    let sequence = normalize_dna(b"NACNNTGN").expect("fixture is normalized");
    assert_eq!(
        CanonicalRuns::new(sequence.bases()).collect::<Vec<_>>(),
        vec![(1, 3), (5, 7)]
    );
    let all_n = normalize_dna(b"NNN").expect("fixture is normalized");
    assert_eq!(CanonicalRuns::new(all_n.bases()).next(), None);
}

#[test]
fn retained_formula_and_synthetic_overflows_are_exact() {
    let observed = retained_fm_bytes(10, 3, 12).expect("tiny dimensions fit");
    let fm_width = u64::try_from(size_of::<FmIndex>()).expect("width fits u64");
    let usize_width = u64::try_from(size_of::<usize>()).expect("width fits u64");
    let rank_width = u64::try_from(size_of::<[u64; 4]>()).expect("width fits u64");
    let expected = 12 * fm_width + 4 * (10 + 3) * (usize_width + 1) + 4 * (10 + 2 * 3) * rank_width;
    assert_eq!(observed, expected);

    assert_eq!(
        checked_build_add(ReferenceResource::CanonicalBases, u64::MAX, 1),
        Err(ReferenceBuildError::ArithmeticOverflow {
            resource: ReferenceResource::CanonicalBases,
            operation: ReferenceArithmetic::Add,
            lhs: u64::MAX,
            rhs: 1,
        })
    );
    assert_eq!(
        checked_build_mul(ReferenceResource::Lanes, u64::MAX, 4),
        Err(ReferenceBuildError::ArithmeticOverflow {
            resource: ReferenceResource::Lanes,
            operation: ReferenceArithmetic::Multiply,
            lhs: u64::MAX,
            rhs: 4,
        })
    );
    assert!(matches!(
        retained_fm_bytes(0, u64::MAX, u64::MAX),
        Err(ReferenceBuildError::ArithmeticOverflow {
            resource: ReferenceResource::EstimatedRetainedFmBytes,
            operation: ReferenceArithmetic::Multiply,
            ..
        })
    ));
}

#[test]
fn capacity_and_count_guards_fail_with_exact_public_context() {
    assert_eq!(
        ensure_build_capacity(3, 2, 2),
        Err(ReferenceBuildError::InternalInvariant {
            expected: 3,
            observed: 2,
        })
    );
    assert_eq!(
        ensure_query_capacity(5, 5),
        Err(ReferenceQueryError::CapacityInvariant {
            reserved: 5,
            materialized: 5,
        })
    );
    assert_eq!(
        ensure_query_count(ReferenceQueryCounter::ExactHits, 7, 6),
        Err(ReferenceQueryError::InvariantMismatch {
            counter: ReferenceQueryCounter::ExactHits,
            expected: 7,
            observed: 6,
        })
    );
    assert_eq!(
        ensure_locate_capacity(11, 11),
        Err(ReferenceLocateError::Invariant {
            invariant: ReferenceLocateInvariant::FinalHitCapacity,
            expected: 11,
            observed: 11,
        })
    );
    assert_eq!(
        ensure_locate_equal(ReferenceLocateInvariant::OffsetCount, 13, 12),
        Err(ReferenceLocateError::Invariant {
            invariant: ReferenceLocateInvariant::OffsetCount,
            expected: 13,
            observed: 12,
        })
    );
    assert_eq!(
        ensure_locate_equal(ReferenceLocateInvariant::FinalHitCount, 17, 16),
        Err(ReferenceLocateError::Invariant {
            invariant: ReferenceLocateInvariant::FinalHitCount,
            expected: 17,
            observed: 16,
        })
    );

    let mut scratch = Vec::new();
    assert_eq!(
        push_lane_base(&mut scratch, 0, Base::A, ThreeLetterConversion::CToT,),
        Err(ReferenceBuildError::InternalInvariant {
            expected: 0,
            observed: 0,
        })
    );
}

#[test]
fn query_storage_length_conversion_is_checked_on_supported_targets() {
    assert_eq!(query_storage_to_logical(0), Ok(0));
    let largest_storage = usize::MAX;
    let largest_logical = u64::try_from(largest_storage)
        .expect("crate compile-time target guard limits usize to at most u64");
    assert_eq!(
        query_storage_to_logical(largest_storage),
        Ok(largest_logical)
    );
}

#[test]
fn allocation_preflights_preserve_their_public_contexts() {
    let element_size = u64::try_from(size_of::<u64>()).expect("width fits u64");
    assert_eq!(
        preflight_query_allocation::<u64>(ReferenceAllocation::ProjectedPattern, u64::MAX,),
        Err(ReferenceQueryError::AllocationSizeOverflow {
            allocation: ReferenceAllocation::ProjectedPattern,
            elements: u64::MAX,
            element_size,
        })
    );
    assert_eq!(
        preflight_query_allocation::<u64>(ReferenceAllocation::OpaqueMatches, u64::MAX,),
        Err(ReferenceQueryError::AllocationSizeOverflow {
            allocation: ReferenceAllocation::OpaqueMatches,
            elements: u64::MAX,
            element_size,
        })
    );
    assert_eq!(
        preflight_locate_allocation::<u64>(ReferenceAllocation::FinalHits, u64::MAX,),
        Err(ReferenceLocateError::AllocationSizeOverflow {
            allocation: ReferenceAllocation::FinalHits,
            elements: u64::MAX,
            element_size,
        })
    );

    let run_width = u64::try_from(size_of::<RunMetadata>()).expect("run width fits u64");
    assert_eq!(
        preflight_build_allocation::<RunMetadata>(ReferenceAllocation::RunMetadata, u64::MAX,),
        Err(ReferenceBuildError::AllocationSizeOverflow {
            allocation: ReferenceAllocation::RunMetadata,
            elements: u64::MAX,
            element_size: run_width,
        })
    );
    assert_eq!(
        preflight_build_allocation::<u64>(ReferenceAllocation::ProjectionScratch, u64::MAX,),
        Err(ReferenceBuildError::AllocationSizeOverflow {
            allocation: ReferenceAllocation::ProjectionScratch,
            elements: u64::MAX,
            element_size,
        })
    );
}
