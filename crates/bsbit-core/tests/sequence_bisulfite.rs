//! Exhaustive and ground-truth tests for sequence and bisulfite primitives.

use std::collections::BTreeSet;
use std::sync::Arc;
use std::thread;

use bsbit_core::alphabet::Base;
use bsbit_core::bisulfite::{
    AlignmentOrientation, BaseRelation, BisulfiteStrand, CytosineStrand,
    InconsistentStrandAssignment, StrandAssignmentAxis, ThreeLetterConversion, classify_bases,
    strand_semantics, validate_strand_assignment,
};
use bsbit_core::sequence::{
    NormalizationError, NormalizedSequence, normalize_dna, reverse_complement, three_letter_convert,
};

const BASES: [Base; 5] = [Base::A, Base::C, Base::G, Base::T, Base::N];

fn seq(text: &str) -> NormalizedSequence {
    normalize_dna(text.as_bytes()).expect("test sequence is valid")
}

fn opposite_orientation(orientation: AlignmentOrientation) -> AlignmentOrientation {
    match orientation {
        AlignmentOrientation::Forward => AlignmentOrientation::Reverse,
        AlignmentOrientation::Reverse => AlignmentOrientation::Forward,
    }
}

fn opposite_cytosine_strand(strand: CytosineStrand) -> CytosineStrand {
    match strand {
        CytosineStrand::Top => CytosineStrand::Bottom,
        CytosineStrand::Bottom => CytosineStrand::Top,
    }
}

fn expected_relation(reference: Base, query: Base, strand: CytosineStrand) -> BaseRelation {
    if reference == Base::N || query == Base::N {
        BaseRelation::Unknown
    } else if reference == query {
        BaseRelation::LiteralMatch
    } else if (strand == CytosineStrand::Top && reference == Base::C && query == Base::T)
        || (strand == CytosineStrand::Bottom && reference == Base::G && query == Base::A)
    {
        BaseRelation::BisulfiteCompatible
    } else {
        BaseRelation::Mismatch
    }
}

fn oracle_complement(base: Base) -> Base {
    match base {
        Base::A => Base::T,
        Base::C => Base::G,
        Base::G => Base::C,
        Base::T => Base::A,
        _ => Base::N,
    }
}

fn oracle_reverse_complement(input: &NormalizedSequence) -> NormalizedSequence {
    NormalizedSequence::from_bases(input.bases().iter().rev().copied().map(oracle_complement))
}

fn oracle_conversion(
    input: &NormalizedSequence,
    conversion: ThreeLetterConversion,
) -> NormalizedSequence {
    NormalizedSequence::from_bases(input.bases().iter().copied().map(|base| {
        match (conversion, base) {
            (ThreeLetterConversion::CToT, Base::C) => Base::T,
            (ThreeLetterConversion::GToA, Base::G) => Base::A,
            _ => base,
        }
    }))
}

fn for_each_sequence_through_len(max_len: usize, mut visit: impl FnMut(&NormalizedSequence)) {
    fn recurse(
        remaining: usize,
        prefix: &mut Vec<Base>,
        visit: &mut impl FnMut(&NormalizedSequence),
    ) {
        if remaining == 0 {
            visit(&NormalizedSequence::from_bases(prefix.iter().copied()));
            return;
        }
        for base in BASES {
            prefix.push(base);
            recurse(remaining - 1, prefix, visit);
            prefix.pop();
        }
    }

    for length in 0..=max_len {
        recurse(length, &mut Vec::with_capacity(length), &mut visit);
    }
}

#[test]
fn base_ascii_and_complement_tables_are_complete() {
    let expected = [
        (Base::A, b'A', Base::T),
        (Base::C, b'C', Base::G),
        (Base::G, b'G', Base::C),
        (Base::T, b'T', Base::A),
        (Base::N, b'N', Base::N),
    ];
    assert_eq!(Base::ALL, BASES);
    assert_eq!(Base::CANONICAL, [Base::A, Base::C, Base::G, Base::T]);
    for (base, ascii, complement) in expected {
        assert_eq!(base.as_ascii(), ascii);
        assert_eq!(base.complement(), complement);
        assert_eq!(base.complement().complement(), base);
        assert_eq!(base.is_unknown(), base == Base::N);
        assert_eq!(base.to_string().as_bytes(), &[ascii]);
    }
}

#[test]
fn normalization_accepts_every_supported_case_and_empty_input() {
    assert_eq!(normalize_dna(b"").unwrap(), NormalizedSequence::default());
    assert_eq!(normalize_dna(b"acgtn").unwrap(), seq("ACGTN"));
    assert_eq!(normalize_dna(b"AaCcGgTtNn").unwrap(), seq("AACCGGTTNN"));
    assert_eq!(seq("ACGTN").to_ascii(), b"ACGTN");
    assert_eq!(seq("ACGTN").to_string(), "ACGTN");
    assert_eq!(seq("ACGTN").get(0), Some(Base::A));
    assert_eq!(seq("ACGTN").get(4), Some(Base::N));
    assert_eq!(seq("ACGTN").get(5), None);
}

#[test]
fn canonical_ascii_snapshot_round_trips_strictly() {
    for_each_sequence_through_len(6, |original| {
        let encoded = original.to_ascii();
        assert_eq!(normalize_dna(&encoded).unwrap(), *original);
        assert!(encoded.iter().all(u8::is_ascii_uppercase));

        let mut with_trailing_garbage = encoded;
        with_trailing_garbage.push(b'?');
        assert_eq!(
            normalize_dna(&with_trailing_garbage).unwrap_err(),
            NormalizationError::InvalidBaseByte {
                byte: b'?',
                offset: original.len(),
            }
        );
    });
}

#[test]
fn every_recognized_iupac_code_has_a_position_aware_error() {
    for &byte in b"RYSWKMBDHVryswkmbdhv" {
        let input = [b'A', b'C', byte, b'G'];
        let error = normalize_dna(&input).unwrap_err();
        assert_eq!(
            error,
            NormalizationError::UnsupportedIupac { byte, offset: 2 }
        );
        assert_eq!(error.byte(), byte);
        assert_eq!(error.offset(), 2);
        assert!(error.to_string().contains("offset 2"));
    }
}

#[test]
fn invalid_non_iupac_bytes_are_not_silently_trimmed_or_normalized() {
    for &byte in &[
        b'U', b'u', b' ', b'\t', b'\n', b'-', b'.', b'=', b'0', b'?', 0, 0xFF,
    ] {
        let input = [b'A', byte, b'C'];
        let error = normalize_dna(&input).unwrap_err();
        assert_eq!(
            error,
            NormalizationError::InvalidBaseByte { byte, offset: 1 }
        );
        assert_eq!(error.byte(), byte);
        assert_eq!(error.offset(), 1);
    }
}

#[test]
fn normalization_reports_the_first_error_and_returns_no_partial_sequence() {
    assert_eq!(
        normalize_dna(b"AR?").unwrap_err(),
        NormalizationError::UnsupportedIupac {
            byte: b'R',
            offset: 1,
        }
    );
    assert_eq!(
        normalize_dna(b"A?R").unwrap_err(),
        NormalizationError::InvalidBaseByte {
            byte: b'?',
            offset: 1,
        }
    );
    assert_eq!(
        normalize_dna(b"?AC").unwrap_err(),
        NormalizationError::InvalidBaseByte {
            byte: b'?',
            offset: 0,
        }
    );
    assert_eq!(
        normalize_dna(b"AC?").unwrap_err(),
        NormalizationError::InvalidBaseByte {
            byte: b'?',
            offset: 2,
        }
    );
}

#[test]
fn reverse_complement_and_named_conversion_fixtures_are_exact() {
    let input = seq("ACGTN");
    let retained = input.clone();
    assert_eq!(reverse_complement(&input), seq("NACGT"));
    assert_eq!(input.reverse_complement(), seq("NACGT"));
    assert_eq!(input, retained);

    assert_eq!(
        three_letter_convert(&seq("ACGTCN"), ThreeLetterConversion::CToT),
        seq("ATGTTN")
    );
    assert_eq!(
        three_letter_convert(&seq("AGGTAN"), ThreeLetterConversion::GToA),
        seq("AAATAN")
    );
}

#[test]
fn conversion_maps_only_its_target_base_and_preserves_length() {
    let input = NormalizedSequence::from_bases(BASES);
    assert_eq!(
        input.three_letter_convert(ThreeLetterConversion::CToT),
        NormalizedSequence::from_bases([Base::A, Base::T, Base::G, Base::T, Base::N])
    );
    assert_eq!(
        input.three_letter_convert(ThreeLetterConversion::GToA),
        NormalizedSequence::from_bases([Base::A, Base::C, Base::A, Base::T, Base::N])
    );
    assert_eq!(input.len(), 5);
    assert!(!input.is_empty());
}

#[test]
fn transforms_obey_exhaustive_properties_through_length_six() {
    for_each_sequence_through_len(6, |input| {
        let original = input.clone();
        let rc = reverse_complement(input);
        assert_eq!(rc, oracle_reverse_complement(input));
        assert_eq!(reverse_complement(&rc), *input);

        for conversion in [ThreeLetterConversion::CToT, ThreeLetterConversion::GToA] {
            let converted = three_letter_convert(input, conversion);
            assert_eq!(converted, oracle_conversion(input, conversion));
            assert_eq!(converted.len(), input.len());
            assert_eq!(three_letter_convert(&converted, conversion), converted);
        }

        let ct_then_rc =
            reverse_complement(&three_letter_convert(input, ThreeLetterConversion::CToT));
        let rc_then_ga = three_letter_convert(&rc, ThreeLetterConversion::GToA);
        assert_eq!(ct_then_rc, rc_then_ga);

        let ga_then_rc =
            reverse_complement(&three_letter_convert(input, ThreeLetterConversion::GToA));
        let rc_then_ct = three_letter_convert(&rc, ThreeLetterConversion::CToT);
        assert_eq!(ga_then_rc, rc_then_ct);
        assert_eq!(*input, original);
    });
}

#[test]
fn ten_thousand_base_boundary_matches_independent_oracles() {
    let input =
        NormalizedSequence::from_bases((0..10_000_usize).map(|index| BASES[index % BASES.len()]));
    let original = input.clone();
    assert_eq!(input.len(), 10_000);
    assert_eq!(normalize_dna(&input.to_ascii()).unwrap(), input);
    assert_eq!(
        reverse_complement(&input),
        oracle_reverse_complement(&input)
    );
    for conversion in [ThreeLetterConversion::CToT, ThreeLetterConversion::GToA] {
        assert_eq!(
            three_letter_convert(&input, conversion),
            oracle_conversion(&input, conversion)
        );
    }
    assert_eq!(input, original);
}

#[test]
fn strand_semantics_match_the_exact_four_strand_table() {
    let table = [
        (
            BisulfiteStrand::OT,
            AlignmentOrientation::Forward,
            CytosineStrand::Top,
            ThreeLetterConversion::CToT,
        ),
        (
            BisulfiteStrand::OB,
            AlignmentOrientation::Reverse,
            CytosineStrand::Bottom,
            ThreeLetterConversion::GToA,
        ),
        (
            BisulfiteStrand::CTOT,
            AlignmentOrientation::Reverse,
            CytosineStrand::Top,
            ThreeLetterConversion::CToT,
        ),
        (
            BisulfiteStrand::CTOB,
            AlignmentOrientation::Forward,
            CytosineStrand::Bottom,
            ThreeLetterConversion::GToA,
        ),
    ];

    assert_eq!(
        BisulfiteStrand::ALL.map(|strand| strand.to_string()),
        ["OT", "OB", "CTOT", "CTOB"]
    );
    assert_eq!(
        [
            AlignmentOrientation::Forward.to_string(),
            AlignmentOrientation::Reverse.to_string(),
        ],
        ["Forward", "Reverse"]
    );
    assert_eq!(
        [
            CytosineStrand::Top.to_string(),
            CytosineStrand::Bottom.to_string(),
        ],
        ["Top", "Bottom"]
    );
    assert_eq!(
        [
            ThreeLetterConversion::CToT.to_string(),
            ThreeLetterConversion::GToA.to_string(),
        ],
        ["CToT", "GToA"]
    );
    for (strand, orientation, cytosine_strand, conversion) in table {
        let semantics = strand_semantics(strand);
        assert_eq!(semantics.strand(), strand);
        assert_eq!(semantics.orientation(), orientation);
        assert_eq!(semantics.cytosine_strand(), cytosine_strand);
        assert_eq!(semantics.search_conversion(), conversion);
        assert_eq!(
            validate_strand_assignment(strand, orientation, cytosine_strand),
            Ok(semantics)
        );
    }
    assert_eq!(
        BisulfiteStrand::ALL
            .into_iter()
            .collect::<BTreeSet<_>>()
            .len(),
        4
    );
}

#[test]
fn every_inconsistent_strand_axis_is_rejected_in_defined_order() {
    for strand in BisulfiteStrand::ALL {
        let expected = strand_semantics(strand);
        let wrong_orientation = opposite_orientation(expected.orientation());
        let wrong_cytosine = opposite_cytosine_strand(expected.cytosine_strand());

        let orientation_error =
            validate_strand_assignment(strand, wrong_orientation, expected.cytosine_strand())
                .unwrap_err();
        assert_eq!(
            orientation_error,
            InconsistentStrandAssignment {
                strand,
                axis: StrandAssignmentAxis::Orientation,
                supplied_orientation: wrong_orientation,
                expected_orientation: expected.orientation(),
                supplied_cytosine_strand: expected.cytosine_strand(),
                expected_cytosine_strand: expected.cytosine_strand(),
            }
        );
        let diagnostic = orientation_error.to_string();
        assert!(diagnostic.contains(&strand.to_string()));
        assert!(diagnostic.contains("Orientation"));

        let cytosine_error =
            validate_strand_assignment(strand, expected.orientation(), wrong_cytosine).unwrap_err();
        assert_eq!(
            cytosine_error,
            InconsistentStrandAssignment {
                strand,
                axis: StrandAssignmentAxis::CytosineStrand,
                supplied_orientation: expected.orientation(),
                expected_orientation: expected.orientation(),
                supplied_cytosine_strand: wrong_cytosine,
                expected_cytosine_strand: expected.cytosine_strand(),
            }
        );

        let both_error =
            validate_strand_assignment(strand, wrong_orientation, wrong_cytosine).unwrap_err();
        assert_eq!(
            both_error,
            InconsistentStrandAssignment {
                strand,
                axis: StrandAssignmentAxis::Orientation,
                supplied_orientation: wrong_orientation,
                expected_orientation: expected.orientation(),
                supplied_cytosine_strand: wrong_cytosine,
                expected_cytosine_strand: expected.cytosine_strand(),
            }
        );
    }
}

#[test]
fn all_twenty_five_pairs_on_each_strand_have_the_exact_relation_and_cost() {
    for strand in [CytosineStrand::Top, CytosineStrand::Bottom] {
        for reference in BASES {
            for query in BASES {
                let expected = expected_relation(reference, query, strand);
                let actual = classify_bases(reference, query, strand);
                assert_eq!(
                    actual, expected,
                    "strand={strand:?}, ref={reference}, query={query}"
                );
                assert_eq!(actual.cost(), u64::from(!actual.is_zero_cost()));
            }
        }
    }
}

#[test]
fn conversion_relation_is_asymmetric_and_unknown_always_costs_one() {
    assert_eq!(
        classify_bases(Base::C, Base::T, CytosineStrand::Top),
        BaseRelation::BisulfiteCompatible
    );
    assert_eq!(
        classify_bases(Base::T, Base::C, CytosineStrand::Top),
        BaseRelation::Mismatch
    );
    assert_eq!(
        classify_bases(Base::G, Base::A, CytosineStrand::Bottom),
        BaseRelation::BisulfiteCompatible
    );
    assert_eq!(
        classify_bases(Base::A, Base::G, CytosineStrand::Bottom),
        BaseRelation::Mismatch
    );
    for strand in [CytosineStrand::Top, CytosineStrand::Bottom] {
        for base in BASES {
            assert_eq!(classify_bases(Base::N, base, strand), BaseRelation::Unknown);
            assert_eq!(classify_bases(base, Base::N, strand), BaseRelation::Unknown);
        }
        assert_eq!(classify_bases(Base::N, Base::N, strand).cost(), 1);
    }
}

#[test]
fn immutable_sequences_and_stateless_operations_are_thread_safe() {
    let input = Arc::new(seq("ACCGTNNGCACT"));
    let expected_rc = input.reverse_complement();
    let expected_ct = input.three_letter_convert(ThreeLetterConversion::CToT);
    let mut handles = Vec::new();
    for _ in 0..8 {
        let input = Arc::clone(&input);
        let expected_rc = expected_rc.clone();
        let expected_ct = expected_ct.clone();
        handles.push(thread::spawn(move || {
            for _ in 0..100 {
                assert_eq!(input.reverse_complement(), expected_rc);
                assert_eq!(
                    input.three_letter_convert(ThreeLetterConversion::CToT),
                    expected_ct
                );
                assert_eq!(
                    classify_bases(Base::C, Base::T, CytosineStrand::Top),
                    BaseRelation::BisulfiteCompatible
                );
            }
        }));
    }
    for handle in handles {
        handle.join().expect("worker did not panic");
    }
}

fn orient_fixture_query(
    raw: &NormalizedSequence,
    orientation: AlignmentOrientation,
) -> NormalizedSequence {
    match orientation {
        AlignmentOrientation::Forward => raw.clone(),
        AlignmentOrientation::Reverse => reverse_complement(raw),
    }
}

fn fixture_relations(
    reference: &NormalizedSequence,
    query: &NormalizedSequence,
    strand: CytosineStrand,
) -> Vec<BaseRelation> {
    reference
        .bases()
        .iter()
        .copied()
        .zip(query.bases().iter().copied())
        .map(|(reference_base, query_base)| classify_bases(reference_base, query_base, strand))
        .collect()
}

fn assert_four_strand_fixture(
    strand: BisulfiteStrand,
    reference_text: &str,
    raw_allowed_text: &str,
    oriented_allowed_text: &str,
    raw_forbidden_text: &str,
    oriented_forbidden_text: &str,
    forbidden_index: usize,
) {
    let semantics = strand_semantics(strand);
    let reference = seq(reference_text);
    let raw_allowed = seq(raw_allowed_text);
    let raw_forbidden = seq(raw_forbidden_text);
    let original_inputs = (
        reference.clone(),
        raw_allowed.clone(),
        raw_forbidden.clone(),
    );

    let oriented_allowed = orient_fixture_query(&raw_allowed, semantics.orientation());
    let oriented_forbidden = orient_fixture_query(&raw_forbidden, semantics.orientation());
    assert_eq!(oriented_allowed, seq(oriented_allowed_text), "{strand}");
    assert_eq!(oriented_forbidden, seq(oriented_forbidden_text), "{strand}");

    let reference_projection = three_letter_convert(&reference, semantics.search_conversion());
    assert_eq!(reference_projection, oriented_allowed, "{strand}");
    assert_eq!(
        three_letter_convert(&oriented_forbidden, semantics.search_conversion()),
        reference_projection,
        "{strand}: search projection deliberately collapses the forbidden direction"
    );

    let allowed = fixture_relations(&reference, &oriented_allowed, semantics.cytosine_strand());
    assert_eq!(
        allowed
            .iter()
            .filter(|&&relation| relation == BaseRelation::BisulfiteCompatible)
            .count(),
        2,
        "{strand}"
    );
    assert!(
        allowed.iter().all(|relation| relation.cost() == 0),
        "{strand}"
    );

    let forbidden = fixture_relations(&reference, &oriented_forbidden, semantics.cytosine_strand());
    assert_eq!(
        forbidden[forbidden_index],
        BaseRelation::Mismatch,
        "{strand}"
    );
    assert_eq!(
        forbidden
            .iter()
            .map(|relation| relation.cost())
            .sum::<u64>(),
        1,
        "{strand}"
    );

    assert_eq!(reference, original_inputs.0, "{strand}");
    assert_eq!(raw_allowed, original_inputs.1, "{strand}");
    assert_eq!(raw_forbidden, original_inputs.2, "{strand}");
}

#[test]
fn non_palindromic_four_strand_fixtures_compose_orientation_and_relation() {
    let fixtures = [
        (
            BisulfiteStrand::OT,
            "ACCGTA",
            "ATTGTA",
            "ATTGTA",
            "ACCGCA",
            "ACCGCA",
            4,
        ),
        (
            BisulfiteStrand::OB,
            "AGGCTA",
            "TAGTTT",
            "AAACTA",
            "TAGCCC",
            "GGGCTA",
            0,
        ),
        (
            BisulfiteStrand::CTOT,
            "ACCGTA",
            "TACAAT",
            "ATTGTA",
            "TGCGGT",
            "ACCGCA",
            4,
        ),
        (
            BisulfiteStrand::CTOB,
            "AGGCTA",
            "AAACTA",
            "AAACTA",
            "GGGCTA",
            "GGGCTA",
            0,
        ),
    ];

    for (
        strand,
        reference_text,
        raw_allowed_text,
        oriented_allowed_text,
        raw_forbidden_text,
        oriented_forbidden_text,
        forbidden_index,
    ) in fixtures
    {
        assert_four_strand_fixture(
            strand,
            reference_text,
            raw_allowed_text,
            oriented_allowed_text,
            raw_forbidden_text,
            oriented_forbidden_text,
            forbidden_index,
        );
    }
}
