#![allow(missing_docs)]

use std::fmt::Write as _;
use std::sync::Arc;
use std::thread;

use bsbit_align::score::EditDistance;
use bsbit_align::verification::cigar::{
    CigarEvaluationError, CigarEvaluationField, evaluate_cigar,
};
use bsbit_align::verification::distance::{
    AlignmentInvariant, DistanceError, DpCellLimit, TraceScoreField, global_bs_alignment,
    global_bs_distance,
};
use bsbit_core::alphabet::Base;
use bsbit_core::bisulfite::CytosineStrand;
use bsbit_core::cigar::{CigarError, CoreCigar, CoreCigarOp, parse_core_cigar, validate_cigar};
use bsbit_core::sequence::NormalizedSequence;

const STRANDS: [CytosineStrand; 2] = [CytosineStrand::Top, CytosineStrand::Bottom];

#[derive(Clone, Debug, Eq, Ord, PartialEq, PartialOrd)]
struct OraclePathScore {
    distance: u64,
    gap_bases: u64,
    gap_runs: u64,
    operations: Vec<CoreCigarOp>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
struct OracleAlignment {
    distance: u64,
    cigar: String,
    multiple_optimal_paths: bool,
}

#[derive(Default)]
struct EnumerationState {
    minimum_distance: Option<u64>,
    minimum_distance_paths: u8,
    best: Option<OraclePathScore>,
}

#[test]
fn named_distance_and_cigar_ground_truth_fixtures() {
    assert_named_fixture(
        "L1-DIST-001",
        CytosineStrand::Top,
        "ACCGT",
        "ATTGT",
        0,
        "5M",
        None,
    );
    assert_named_fixture(
        "L1-DIST-002",
        CytosineStrand::Top,
        "ATTGT",
        "ACCGT",
        2,
        "5M",
        None,
    );
    assert_named_fixture(
        "L1-DIST-003",
        CytosineStrand::Bottom,
        "ACGGT",
        "ACAAT",
        0,
        "5M",
        None,
    );
    for strand in STRANDS {
        assert_named_fixture(
            "L1-DIST-004",
            strand,
            "ACGT",
            "ACT",
            1,
            "2M1D1M",
            Some(false),
        );
        assert_named_fixture(
            "L1-DIST-005",
            strand,
            "ACT",
            "ACGT",
            1,
            "2M1I1M",
            Some(false),
        );
        assert_named_fixture("L1-DIST-006", strand, "", "AC", 2, "2I", Some(false));
        assert_named_fixture("L1-DIST-007", strand, "AC", "", 2, "2D", Some(false));
    }
    assert_named_fixture(
        "L1-DIST-008",
        CytosineStrand::Top,
        "ACCGT",
        "ACCGT",
        0,
        "5M",
        Some(false),
    );
    assert_named_fixture(
        "L1-DIST-009",
        CytosineStrand::Bottom,
        "ACGGT",
        "ACGGT",
        0,
        "5M",
        Some(false),
    );
    assert_named_fixture(
        "L1-CIGAR-001",
        CytosineStrand::Top,
        "AAA",
        "AA",
        1,
        "1D2M",
        Some(true),
    );
    assert_named_fixture(
        "L1-CIGAR-002",
        CytosineStrand::Top,
        "AA",
        "AAA",
        1,
        "1I2M",
        Some(true),
    );
    assert_named_fixture(
        "L1-CIGAR-003",
        CytosineStrand::Top,
        "AC",
        "CA",
        2,
        "2M",
        Some(true),
    );
}

#[test]
fn exhaustive_distance_matches_independent_full_matrix_through_length_four() {
    let sequences = all_sequences(&Base::CANONICAL, 4);

    for strand in STRANDS {
        for reference in &sequences {
            for query in &sequences {
                let expected = full_matrix_oracle(reference, query, strand);
                let actual = global_bs_distance(reference, query, strand, DpCellLimit::MAX)
                    .expect("small exact distance must fit");
                assert_eq!(
                    actual.get(),
                    expected,
                    "strand={strand:?}, reference={reference}, query={query}"
                );
                assert!(actual.get() <= reference.len().max(query.len()));
            }
        }
    }
}

#[test]
fn exhaustive_traceback_matches_complete_path_enumerator_through_length_three() {
    let sequences = all_sequences(&Base::CANONICAL, 3);

    for strand in STRANDS {
        for reference in &sequences {
            for query in &sequences {
                let oracle =
                    brute_force_alignment(reference, query, strand, CoreCigarOp::LEXICOGRAPHIC);
                assert_eq!(
                    oracle.distance,
                    full_matrix_oracle(reference, query, strand),
                    "the two independent test oracles disagree for strand={strand:?}, reference={reference}, query={query}"
                );

                let distance = global_bs_distance(reference, query, strand, DpCellLimit::MAX)
                    .expect("short exact distance must fit");
                let alignment = global_bs_alignment(reference, query, strand, DpCellLimit::MAX)
                    .expect("short exact traceback must fit");

                assert_eq!(
                    distance.get(),
                    oracle.distance,
                    "distance: strand={strand:?}, reference={reference}, query={query}"
                );
                assert_eq!(
                    alignment.distance().get(),
                    oracle.distance,
                    "traceback distance: strand={strand:?}, reference={reference}, query={query}"
                );
                assert_eq!(
                    alignment.cigar().to_string(),
                    oracle.cigar,
                    "canonical CIGAR: strand={strand:?}, reference={reference}, query={query}"
                );
                assert_eq!(
                    alignment.multiple_optimal_paths(),
                    oracle.multiple_optimal_paths,
                    "ambiguity: strand={strand:?}, reference={reference}, query={query}"
                );

                assert_canonical_runs(alignment.cigar());
                let consumption = validate_cigar(alignment.cigar(), reference.len(), query.len())
                    .expect("traceback CIGAR must consume both inputs exactly");
                assert_eq!(consumption.reference(), reference.len());
                assert_eq!(consumption.query(), query.len());
                let replay = evaluate_cigar(alignment.cigar(), reference, query, strand)
                    .expect("traceback CIGAR must replay");
                assert_eq!(replay.distance(), alignment.distance());
            }
        }
    }
}

#[test]
fn complete_path_oracle_is_independent_of_enumeration_order() {
    let reference = sequence("ACAC");
    let query = sequence("CAAC");
    let orders = [
        [CoreCigarOp::D, CoreCigarOp::I, CoreCigarOp::M],
        [CoreCigarOp::D, CoreCigarOp::M, CoreCigarOp::I],
        [CoreCigarOp::I, CoreCigarOp::D, CoreCigarOp::M],
        [CoreCigarOp::I, CoreCigarOp::M, CoreCigarOp::D],
        [CoreCigarOp::M, CoreCigarOp::D, CoreCigarOp::I],
        [CoreCigarOp::M, CoreCigarOp::I, CoreCigarOp::D],
    ];

    for strand in STRANDS {
        let expected = brute_force_alignment(&reference, &query, strand, orders[0]);
        for order in &orders[1..] {
            assert_eq!(
                brute_force_alignment(&reference, &query, strand, *order),
                expected
            );
        }
    }
}

#[test]
fn reverse_complement_duality_holds_with_unknown_bases() {
    let sequences = all_sequences(&Base::ALL, 3);

    for reference in &sequences {
        let reverse_reference = reference.reverse_complement();
        for query in &sequences {
            let reverse_query = query.reverse_complement();
            let top = global_bs_distance(reference, query, CytosineStrand::Top, DpCellLimit::MAX)
                .expect("short top-strand distance must fit");
            let bottom = global_bs_distance(
                &reverse_reference,
                &reverse_query,
                CytosineStrand::Bottom,
                DpCellLimit::MAX,
            )
            .expect("short bottom-strand distance must fit");
            assert_eq!(
                top, bottom,
                "reference={reference}, query={query}, rc(reference)={reverse_reference}, rc(query)={reverse_query}"
            );
        }
    }
}

#[test]
fn unknown_columns_have_unit_cost_and_match_the_independent_oracle() {
    let sequences = all_sequences(&Base::ALL, 2);
    for strand in STRANDS {
        for reference in &sequences {
            for query in &sequences {
                let expected = full_matrix_oracle(reference, query, strand);
                let actual = global_bs_distance(reference, query, strand, DpCellLimit::MAX)
                    .expect("tiny N-containing distance must fit");
                assert_eq!(
                    actual.get(),
                    expected,
                    "strand={strand:?}, reference={reference}, query={query}"
                );
            }
        }
    }

    for strand in STRANDS {
        assert_named_fixture(
            "N/N is not identity",
            strand,
            "N",
            "N",
            1,
            "1M",
            Some(false),
        );
        assert_eq!(
            global_bs_distance(&sequence("N"), &sequence("A"), strand, DpCellLimit::MAX,)
                .expect("single-column distance")
                .get(),
            1
        );
    }
}

#[test]
fn every_retained_or_converted_methylation_mask_remains_zero_cost() {
    assert_conversion_masks("ACCGCC", Base::C, Base::T, CytosineStrand::Top);
    assert_conversion_masks("AGGCGG", Base::G, Base::A, CytosineStrand::Bottom);
}

#[test]
fn cigar_replay_counts_literal_conversion_mismatch_unknown_and_gaps() {
    let reference = sequence("ACGTN");
    let query = sequence("ATCGN");
    let reference_before = reference.clone();
    let query_before = query.clone();
    let cigar = parse_core_cigar("2M1D1M1I1M").expect("canonical fixture CIGAR");
    let evaluation = evaluate_cigar(&cigar, &reference, &query, CytosineStrand::Top)
        .expect("fixture CIGAR consumes both sequences");

    assert_eq!(evaluation.distance(), EditDistance::new(4));
    assert_eq!(evaluation.literal_matches(), 1);
    assert_eq!(evaluation.bisulfite_compatible(), 1);
    assert_eq!(evaluation.mismatches(), 1);
    assert_eq!(evaluation.unknown_columns(), 1);
    assert_eq!(evaluation.inserted_bases(), 1);
    assert_eq!(evaluation.deleted_bases(), 1);
    assert_eq!(evaluation.gap_runs(), 2);
    assert_eq!(reference, reference_before);
    assert_eq!(query, query_before);
}

#[test]
fn computation_limits_are_checked_at_the_exact_logical_cell_boundary() {
    let reference = sequence("AC");
    let query = sequence("AGT");
    let requested_cells = 12;

    for error in [
        global_bs_distance(
            &reference,
            &query,
            CytosineStrand::Top,
            DpCellLimit::new(requested_cells - 1),
        )
        .expect_err("distance must reject an undersized cell limit"),
        global_bs_alignment(
            &reference,
            &query,
            CytosineStrand::Top,
            DpCellLimit::new(requested_cells - 1),
        )
        .expect_err("traceback must reject an undersized cell limit"),
    ] {
        assert_eq!(
            error,
            DistanceError::ComputationLimitExceeded {
                requested_cells,
                limit: requested_cells - 1,
            }
        );
    }

    let at_limit = global_bs_alignment(
        &reference,
        &query,
        CytosineStrand::Top,
        DpCellLimit::new(requested_cells),
    )
    .expect("exact logical cell limit is inclusive");
    assert_eq!(
        at_limit.distance().get(),
        full_matrix_oracle(&reference, &query, CytosineStrand::Top)
    );

    let empty = sequence("");
    assert_eq!(
        global_bs_distance(&empty, &empty, CytosineStrand::Top, DpCellLimit::new(0),),
        Err(DistanceError::ComputationLimitExceeded {
            requested_cells: 1,
            limit: 0,
        })
    );
    assert_eq!(DpCellLimit::new(7).get(), 7);
    assert_eq!(DpCellLimit::MAX.get(), u64::MAX);
}

#[test]
fn alignment_is_deterministic_concurrent_and_does_not_mutate_shared_inputs() {
    let reference = Arc::new(sequence("ACCGTNNACGTAC"));
    let query = Arc::new(sequence("ATTGTNAACGAC"));
    let reference_before = reference.to_ascii();
    let query_before = query.to_ascii();
    let expected = global_bs_alignment(&reference, &query, CytosineStrand::Top, DpCellLimit::MAX)
        .expect("fixture traceback must fit");

    for _ in 0..64 {
        assert_eq!(
            global_bs_alignment(&reference, &query, CytosineStrand::Top, DpCellLimit::MAX,)
                .expect("repeated fixture traceback must fit"),
            expected
        );
    }

    let mut workers = Vec::new();
    for _ in 0..8 {
        let reference = Arc::clone(&reference);
        let query = Arc::clone(&query);
        let expected = expected.clone();
        workers.push(thread::spawn(move || {
            for _ in 0..32 {
                let actual =
                    global_bs_alignment(&reference, &query, CytosineStrand::Top, DpCellLimit::MAX)
                        .expect("concurrent fixture traceback must fit");
                assert_eq!(actual, expected);
            }
        }));
    }
    for worker in workers {
        worker.join().expect("alignment worker must not panic");
    }

    assert_eq!(reference.to_ascii(), reference_before);
    assert_eq!(query.to_ascii(), query_before);
}

fn assert_named_fixture(
    id: &str,
    strand: CytosineStrand,
    reference_ascii: &str,
    query_ascii: &str,
    expected_distance: u64,
    expected_cigar: &str,
    expected_multiple: Option<bool>,
) {
    let reference = sequence(reference_ascii);
    let query = sequence(query_ascii);
    let reference_before = reference.clone();
    let query_before = query.clone();
    let oracle = brute_force_alignment(&reference, &query, strand, CoreCigarOp::LEXICOGRAPHIC);
    let distance = global_bs_distance(&reference, &query, strand, DpCellLimit::MAX)
        .unwrap_or_else(|error| panic!("{id}: distance failed: {error}"));
    let alignment = global_bs_alignment(&reference, &query, strand, DpCellLimit::MAX)
        .unwrap_or_else(|error| panic!("{id}: traceback failed: {error}"));

    assert_eq!(distance.get(), expected_distance, "{id}: distance");
    assert_eq!(alignment.distance(), distance, "{id}: distance APIs");
    assert_eq!(alignment.cigar().to_string(), expected_cigar, "{id}: CIGAR");
    assert_eq!(oracle.distance, expected_distance, "{id}: oracle distance");
    assert_eq!(oracle.cigar, expected_cigar, "{id}: oracle CIGAR");
    assert_eq!(
        alignment.multiple_optimal_paths(),
        oracle.multiple_optimal_paths,
        "{id}: implementation/oracle ambiguity"
    );
    if let Some(expected_multiple) = expected_multiple {
        assert_eq!(
            alignment.multiple_optimal_paths(),
            expected_multiple,
            "{id}: expected ambiguity"
        );
    }
    let replay = evaluate_cigar(alignment.cigar(), &reference, &query, strand)
        .unwrap_or_else(|error| panic!("{id}: replay failed: {error}"));
    assert_eq!(replay.distance(), distance, "{id}: replay distance");
    assert_eq!(reference, reference_before, "{id}: reference mutated");
    assert_eq!(query, query_before, "{id}: query mutated");
}

fn full_matrix_oracle(
    reference: &NormalizedSequence,
    query: &NormalizedSequence,
    strand: CytosineStrand,
) -> u64 {
    let rows = reference.bases().len() + 1;
    let columns = query.bases().len() + 1;
    let mut matrix = vec![vec![0_u64; columns]; rows];

    for (row, cells) in matrix.iter_mut().enumerate().skip(1) {
        cells[0] = len_u64(row);
    }
    for (column, cell) in matrix[0].iter_mut().enumerate().skip(1) {
        *cell = len_u64(column);
    }
    for row in 1..rows {
        for column in 1..columns {
            let deletion = matrix[row - 1][column] + 1;
            let insertion = matrix[row][column - 1] + 1;
            let diagonal = matrix[row - 1][column - 1]
                + oracle_substitution_cost(
                    reference.bases()[row - 1],
                    query.bases()[column - 1],
                    strand,
                );
            matrix[row][column] = deletion.min(insertion).min(diagonal);
        }
    }
    matrix[rows - 1][columns - 1]
}

fn brute_force_alignment(
    reference: &NormalizedSequence,
    query: &NormalizedSequence,
    strand: CytosineStrand,
    operation_order: [CoreCigarOp; 3],
) -> OracleAlignment {
    let mut state = EnumerationState::default();
    let mut operations = Vec::with_capacity(reference.bases().len() + query.bases().len());
    enumerate_paths(
        reference,
        query,
        strand,
        operation_order,
        0,
        0,
        0,
        0,
        0,
        &mut operations,
        &mut state,
    );

    let best = state
        .best
        .expect("every finite sequence pair has a global path");
    OracleAlignment {
        distance: best.distance,
        cigar: oracle_cigar(&best.operations),
        multiple_optimal_paths: state.minimum_distance_paths > 1,
    }
}

#[allow(clippy::too_many_arguments)]
fn enumerate_paths(
    reference: &NormalizedSequence,
    query: &NormalizedSequence,
    strand: CytosineStrand,
    operation_order: [CoreCigarOp; 3],
    reference_index: usize,
    query_index: usize,
    distance: u64,
    gap_bases: u64,
    gap_runs: u64,
    operations: &mut Vec<CoreCigarOp>,
    state: &mut EnumerationState,
) {
    if reference_index == reference.bases().len() && query_index == query.bases().len() {
        match state.minimum_distance {
            None => {
                state.minimum_distance = Some(distance);
                state.minimum_distance_paths = 1;
            }
            Some(current) if distance < current => {
                state.minimum_distance = Some(distance);
                state.minimum_distance_paths = 1;
            }
            Some(current) if distance == current => {
                state.minimum_distance_paths = (state.minimum_distance_paths + 1).min(2);
            }
            Some(_) => {}
        }

        let candidate = OraclePathScore {
            distance,
            gap_bases,
            gap_runs,
            operations: operations.clone(),
        };
        if state.best.as_ref().is_none_or(|best| candidate < *best) {
            state.best = Some(candidate);
        }
        return;
    }

    for &operation in &operation_order {
        let (next_reference, next_query, cost) = match operation {
            CoreCigarOp::D if reference_index < reference.bases().len() => {
                (reference_index + 1, query_index, 1)
            }
            CoreCigarOp::I if query_index < query.bases().len() => {
                (reference_index, query_index + 1, 1)
            }
            CoreCigarOp::M
                if reference_index < reference.bases().len()
                    && query_index < query.bases().len() =>
            {
                (
                    reference_index + 1,
                    query_index + 1,
                    oracle_substitution_cost(
                        reference.bases()[reference_index],
                        query.bases()[query_index],
                        strand,
                    ),
                )
            }
            CoreCigarOp::D | CoreCigarOp::I | CoreCigarOp::M => continue,
        };
        let new_gap_run =
            u64::from(operation.is_gap() && operations.last().copied() != Some(operation));
        operations.push(operation);
        enumerate_paths(
            reference,
            query,
            strand,
            operation_order,
            next_reference,
            next_query,
            distance + cost,
            gap_bases + u64::from(operation.is_gap()),
            gap_runs + new_gap_run,
            operations,
            state,
        );
        operations.pop();
    }
}

fn oracle_substitution_cost(reference: Base, query: Base, strand: CytosineStrand) -> u64 {
    if reference == Base::N || query == Base::N {
        1
    } else {
        u64::from(
            !(reference == query
                || matches!(
                    (strand, reference, query),
                    (CytosineStrand::Top, Base::C, Base::T)
                        | (CytosineStrand::Bottom, Base::G, Base::A)
                )),
        )
    }
}

fn oracle_cigar(operations: &[CoreCigarOp]) -> String {
    let mut cigar = String::new();
    let Some(mut operation) = operations.first().copied() else {
        return cigar;
    };
    let mut run_length = 0_u64;
    for &next in operations {
        if next == operation {
            run_length += 1;
        } else {
            write!(cigar, "{run_length}{operation}").expect("writing to String cannot fail");
            operation = next;
            run_length = 1;
        }
    }
    write!(cigar, "{run_length}{operation}").expect("writing to String cannot fail");
    cigar
}

fn assert_canonical_runs(cigar: &CoreCigar) {
    let mut previous = None;
    for run in cigar.runs() {
        assert!(run.length() > 0);
        assert_ne!(previous, Some(run.operation()));
        previous = Some(run.operation());
    }
}

fn assert_conversion_masks(
    reference_ascii: &str,
    retained: Base,
    converted: Base,
    strand: CytosineStrand,
) {
    let reference = sequence(reference_ascii);
    let positions: Vec<usize> = reference
        .bases()
        .iter()
        .enumerate()
        .filter_map(|(index, &base)| (base == retained).then_some(index))
        .collect();
    for mask in 0..(1_usize << positions.len()) {
        let mut query_bases = reference.bases().to_vec();
        for (bit, &position) in positions.iter().enumerate() {
            if mask & (1_usize << bit) != 0 {
                query_bases[position] = converted;
            }
        }
        let query = NormalizedSequence::from_bases(query_bases);
        let alignment = global_bs_alignment(&reference, &query, strand, DpCellLimit::MAX)
            .expect("small methylation mask must fit");
        assert_eq!(alignment.distance(), EditDistance::new(0), "mask={mask:#b}");
        assert_eq!(
            alignment.cigar().to_string(),
            format!("{}M", reference.len()),
            "mask={mask:#b}"
        );
    }
}

fn all_sequences(alphabet: &[Base], maximum_length: usize) -> Vec<NormalizedSequence> {
    let mut output = Vec::new();
    let mut prefix = Vec::new();
    append_sequences(alphabet, maximum_length, &mut prefix, &mut output);
    output
}

fn append_sequences(
    alphabet: &[Base],
    maximum_length: usize,
    prefix: &mut Vec<Base>,
    output: &mut Vec<NormalizedSequence>,
) {
    output.push(NormalizedSequence::from_bases(prefix.iter().copied()));
    if prefix.len() == maximum_length {
        return;
    }
    for &base in alphabet {
        prefix.push(base);
        append_sequences(alphabet, maximum_length, prefix, output);
        prefix.pop();
    }
}

fn sequence(ascii: &str) -> NormalizedSequence {
    let bases = ascii.bytes().map(|byte| match byte {
        b'A' => Base::A,
        b'C' => Base::C,
        b'G' => Base::G,
        b'T' => Base::T,
        b'N' => Base::N,
        _ => panic!("test helper received non-normalized byte 0x{byte:02X}"),
    });
    NormalizedSequence::from_bases(bases)
}

fn len_u64(length: usize) -> u64 {
    u64::try_from(length).expect("test corpus length fits u64")
}

#[test]
fn invariant_and_counter_error_variants_preserve_fields_and_diagnostics() {
    for field in [
        CigarEvaluationField::LiteralMatches,
        CigarEvaluationField::BisulfiteCompatible,
        CigarEvaluationField::Mismatches,
        CigarEvaluationField::UnknownColumns,
        CigarEvaluationField::InsertedBases,
        CigarEvaluationField::DeletedBases,
        CigarEvaluationField::GapRuns,
        CigarEvaluationField::Distance,
    ] {
        let error = CigarEvaluationError::CounterOverflow { field };
        assert!(matches!(
            &error,
            CigarEvaluationError::CounterOverflow { field: observed } if *observed == field
        ));
        assert_eq!(
            error.to_string(),
            format!("CIGAR evaluation counter {field:?} overflowed")
        );
        let standard_error: &dyn std::error::Error = &error;
        assert!(!standard_error.to_string().is_empty());
    }

    for field in [TraceScoreField::GapBases, TraceScoreField::GapRuns] {
        let error = DistanceError::TraceScoreOverflow {
            field,
            accumulated: u64::MAX,
            increment: 1,
        };
        assert!(matches!(
            &error,
            DistanceError::TraceScoreOverflow {
                field: observed,
                accumulated: u64::MAX,
                increment: 1,
            } if *observed == field
        ));
        assert_eq!(
            error.to_string(),
            format!("trace score {field:?} addition {} + 1 overflowed", u64::MAX)
        );
    }

    let nested = CigarError::ZeroLengthCigarRun {
        run_index: 7,
        operation: CoreCigarOp::I,
    };
    let cigar_invariant = DistanceError::from(nested.clone());
    assert_eq!(
        cigar_invariant,
        DistanceError::CigarInvariant { error: nested }
    );
    assert_eq!(
        cigar_invariant.to_string(),
        "constructed CIGAR failed validation: CIGAR run 7 has zero length for operation I"
    );

    for invariant in [
        AlignmentInvariant::MissingSuffixPath,
        AlignmentInvariant::MissingTraceStep,
        AlignmentInvariant::PrimaryDistanceMismatch,
        AlignmentInvariant::CigarDistanceMismatch,
    ] {
        let error = DistanceError::AlignmentInvariant {
            invariant,
            reference_index: 11,
            query_index: 13,
            expected: Some(17),
            observed: Some(19),
        };
        assert!(matches!(
            &error,
            DistanceError::AlignmentInvariant {
                invariant: observed_invariant,
                reference_index: 11,
                query_index: 13,
                expected: Some(17),
                observed: Some(19),
            } if *observed_invariant == invariant
        ));
        assert_eq!(
            error.to_string(),
            format!(
                "alignment invariant {invariant:?} failed at prefixes 11/13; expected=Some(17), observed=Some(19)"
            )
        );
    }
}
