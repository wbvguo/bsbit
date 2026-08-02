//! Public-contract tests for logical CIGAR values and validation.

use bsbit_core::cigar::{
    CigarDomain, CigarError, CoreCigar, CoreCigarOp, RawCigarRun, RawCoreCigar, canonicalize_cigar,
    canonicalize_operations, parse_core_cigar, try_core_cigar, validate_cigar,
};

#[test]
fn cigar_operations_and_canonicalization_have_exact_semantics() {
    assert_eq!(
        CoreCigarOp::LEXICOGRAPHIC,
        [CoreCigarOp::D, CoreCigarOp::I, CoreCigarOp::M]
    );
    assert!(CoreCigarOp::D.consumes_reference());
    assert!(!CoreCigarOp::D.consumes_query());
    assert!(!CoreCigarOp::I.consumes_reference());
    assert!(CoreCigarOp::I.consumes_query());
    assert!(CoreCigarOp::M.consumes_reference());
    assert!(CoreCigarOp::M.consumes_query());
    assert!(CoreCigarOp::D.is_gap());
    assert!(CoreCigarOp::I.is_gap());
    assert!(!CoreCigarOp::M.is_gap());

    let raw_run = RawCigarRun::new(CoreCigarOp::D, 2);
    assert_eq!(raw_run.operation(), CoreCigarOp::D);
    assert_eq!(raw_run.length(), 2);

    let raw = RawCoreCigar::new([
        RawCigarRun::new(CoreCigarOp::D, 1),
        RawCigarRun::new(CoreCigarOp::D, 2),
        RawCigarRun::new(CoreCigarOp::I, 1),
        RawCigarRun::new(CoreCigarOp::M, 2),
        RawCigarRun::new(CoreCigarOp::M, 3),
    ]);
    assert_eq!(raw.runs().len(), 5);
    let canonical = canonicalize_cigar(&raw).expect("positive runs can be coalesced");
    assert_eq!(canonical.to_string(), "3D1I5M");
    assert_eq!(canonical.run_count(), 3);
    assert_canonical_runs(&canonical);

    let from_operations = canonicalize_operations([
        CoreCigarOp::D,
        CoreCigarOp::D,
        CoreCigarOp::I,
        CoreCigarOp::M,
        CoreCigarOp::M,
    ])
    .expect("finite expanded operation stream");
    assert_eq!(from_operations.to_string(), "2D1I2M");

    let empty = canonicalize_operations([]).expect("empty operation stream is valid");
    assert!(empty.is_empty());
    assert_eq!(empty.to_string(), "");
}

#[test]
fn strict_cigar_construction_rejects_zero_adjacent_and_length_mismatch() {
    let zero = RawCoreCigar::new([RawCigarRun::new(CoreCigarOp::I, 0)]);
    assert_eq!(
        try_core_cigar(&zero, 0, 0),
        Err(CigarError::ZeroLengthCigarRun {
            run_index: 0,
            operation: CoreCigarOp::I,
        })
    );
    assert_eq!(
        canonicalize_cigar(&zero),
        Err(CigarError::ZeroLengthCigarRun {
            run_index: 0,
            operation: CoreCigarOp::I,
        })
    );

    let adjacent = RawCoreCigar::new([
        RawCigarRun::new(CoreCigarOp::M, 1),
        RawCigarRun::new(CoreCigarOp::M, 1),
    ]);
    assert_eq!(
        try_core_cigar(&adjacent, 2, 2),
        Err(CigarError::NonCanonicalCigarRuns {
            previous_run_index: 0,
            run_index: 1,
            operation: CoreCigarOp::M,
        })
    );

    let two_matches = parse_core_cigar("2M").expect("canonical CIGAR");
    assert_eq!(
        validate_cigar(&two_matches, 3, 4),
        Err(CigarError::CigarLengthMismatch {
            expected_reference: 3,
            observed_reference: 2,
            expected_query: 4,
            observed_query: 2,
        })
    );
}

#[test]
fn cigar_consumption_and_coalescing_overflow_report_exact_fields() {
    let coalescing = RawCoreCigar::new([
        RawCigarRun::new(CoreCigarOp::M, u64::MAX),
        RawCigarRun::new(CoreCigarOp::M, 1),
    ]);
    assert_eq!(
        canonicalize_cigar(&coalescing),
        Err(CigarError::CigarConsumptionOverflow {
            run_index: 1,
            operation: CoreCigarOp::M,
            domain: CigarDomain::Reference,
            accumulated: u64::MAX,
            run_length: 1,
        })
    );

    let reference_overflow = RawCoreCigar::new([
        RawCigarRun::new(CoreCigarOp::D, u64::MAX),
        RawCigarRun::new(CoreCigarOp::I, 1),
        RawCigarRun::new(CoreCigarOp::D, 1),
    ]);
    assert_eq!(
        try_core_cigar(&reference_overflow, 0, 0),
        Err(CigarError::CigarConsumptionOverflow {
            run_index: 2,
            operation: CoreCigarOp::D,
            domain: CigarDomain::Reference,
            accumulated: u64::MAX,
            run_length: 1,
        })
    );

    let query_overflow = RawCoreCigar::new([
        RawCigarRun::new(CoreCigarOp::I, u64::MAX),
        RawCigarRun::new(CoreCigarOp::D, 1),
        RawCigarRun::new(CoreCigarOp::I, 1),
    ]);
    assert_eq!(
        try_core_cigar(&query_overflow, 0, 0),
        Err(CigarError::CigarConsumptionOverflow {
            run_index: 2,
            operation: CoreCigarOp::I,
            domain: CigarDomain::Query,
            accumulated: u64::MAX,
            run_length: 1,
        })
    );
}

#[test]
fn logical_cigar_parser_is_strict_and_round_trips() {
    for logical in ["", "1M", "2M1D3I4M", "18446744073709551615D"] {
        let parsed = parse_core_cigar(logical).expect("valid canonical logical CIGAR");
        assert_eq!(parsed.to_string(), logical);
        let from_trait: CoreCigar = logical.parse().expect("FromStr uses the strict parser");
        assert_eq!(from_trait, parsed);
    }

    assert_eq!(
        parse_core_cigar("0M"),
        Err(CigarError::ZeroLengthCigarRun {
            run_index: 0,
            operation: CoreCigarOp::M,
        })
    );
    assert_eq!(
        parse_core_cigar("1M2M"),
        Err(CigarError::NonCanonicalCigarRuns {
            previous_run_index: 0,
            run_index: 1,
            operation: CoreCigarOp::M,
        })
    );
    assert_eq!(
        parse_core_cigar("M"),
        Err(CigarError::ExpectedCigarRunLength {
            byte_offset: 0,
            found: b'M',
        })
    );
    assert_eq!(
        parse_core_cigar("1"),
        Err(CigarError::MissingCigarOperation {
            run_index: 0,
            byte_offset: 1,
        })
    );
    assert_eq!(
        parse_core_cigar("1X"),
        Err(CigarError::UnknownCigarOperation {
            run_index: 0,
            byte_offset: 1,
            found: b'X',
        })
    );
    assert_eq!(
        parse_core_cigar("1Mgarbage"),
        Err(CigarError::ExpectedCigarRunLength {
            byte_offset: 2,
            found: b'g',
        })
    );
    assert_eq!(
        parse_core_cigar("18446744073709551616M"),
        Err(CigarError::CigarRunLengthOverflow {
            run_index: 0,
            byte_offset: 19,
        })
    );
    for (logical, run_index, byte_offset) in [("00M", 0, 0), ("01M", 0, 0), ("1M01D", 1, 2)] {
        assert_eq!(
            parse_core_cigar(logical),
            Err(CigarError::NonCanonicalCigarRunLength {
                run_index,
                byte_offset,
            }),
            "leading-zero run length must be rejected: {logical}"
        );
    }
}

#[test]
fn canonicalize_operations_does_not_trust_an_unbounded_size_hint() {
    struct EmptyWithAdversarialHint;

    impl Iterator for EmptyWithAdversarialHint {
        type Item = CoreCigarOp;

        fn next(&mut self) -> Option<Self::Item> {
            None
        }

        fn size_hint(&self) -> (usize, Option<usize>) {
            (usize::MAX, None)
        }
    }

    let result = std::panic::catch_unwind(|| canonicalize_operations(EmptyWithAdversarialHint));
    let value = result.expect("an iterator hint must not trigger an allocation panic");
    assert!(
        value
            .expect("an empty operation stream is valid")
            .is_empty()
    );
}

fn assert_canonical_runs(cigar: &CoreCigar) {
    let mut previous = None;
    for run in cigar.runs() {
        assert!(run.length() > 0);
        assert_ne!(previous, Some(run.operation()));
        previous = Some(run.operation());
    }
}
