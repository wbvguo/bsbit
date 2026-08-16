//! Independent scientific and API tests for Level 3 scalar extension.

#[path = "support/extension_oracle.rs"]
mod extension_oracle;

use extension_oracle::{
    OraclePlacement, OracleStrand, candidate_window_best, enumerate_strings, to_u64,
    whole_contig_best,
};

use std::error::Error as _;

use bsbit_align::extension::{
    ExtensionCounter, ExtensionError, ExtensionLimits, VerifiedAlignment, extend_candidate_window,
};
use bsbit_align::score::EditDistance;
use bsbit_align::search::candidate::{
    CandidateLimits, CandidateSet, FixedSeedPlan, FixedSeedRequest, SeedPlanLimits,
    candidates_for_fixed_seeds,
};
use bsbit_core::bisulfite::BisulfiteStrand;
use bsbit_core::coordinate::{CoordinateShift, QueryInterval, QueryLength};
use bsbit_core::sequence::{NormalizedSequence, normalize_dna};
use bsbit_index::reference::{
    ContigInput, ReferenceBuildLimits, ReferenceIndex, ReferenceQueryLimits,
};

const ALPHABET: [u8; 3] = [b'A', b'C', b'T'];

fn normalized(raw: &[u8]) -> NormalizedSequence {
    normalize_dna(raw).expect("test input uses only A/C/G/T/N")
}

fn build_reference(catalog: &[(&[u8], &[u8])]) -> ReferenceIndex {
    let contigs = catalog
        .iter()
        .map(|(name, sequence)| ContigInput::new(name.to_vec(), normalized(sequence)))
        .collect();
    ReferenceIndex::build(contigs, ReferenceBuildLimits::MAX).expect("bounded reference builds")
}

fn fixed_candidates(
    reference: &ReferenceIndex,
    query: &[u8],
    strand: BisulfiteStrand,
    seed_start: usize,
    seed_end: usize,
) -> CandidateSet {
    let query = normalized(query);
    let interval = QueryInterval::new(
        to_u64(seed_start),
        to_u64(seed_end),
        QueryLength::new(query.len()),
    )
    .expect("seed interval fits");
    let plan = FixedSeedPlan::new(
        query,
        &[FixedSeedRequest::new(strand, interval)],
        SeedPlanLimits::MAX,
    )
    .expect("plan builds");
    candidates_for_fixed_seeds(
        reference,
        &plan,
        ReferenceQueryLimits::MAX,
        CandidateLimits::MAX,
    )
    .expect("candidate generation succeeds")
}

fn diagonal(candidates: &CandidateSet, ordinal: usize) -> i128 {
    match candidates.anchors()[ordinal].diagonal().shift() {
        CoordinateShift::Zero => 0,
        CoordinateShift::Forward(value) => i128::from(value.get()),
        CoordinateShift::Backward(value) => -i128::from(value.get()),
    }
}

const fn oracle_strand(strand: BisulfiteStrand) -> OracleStrand {
    match strand {
        BisulfiteStrand::OT => OracleStrand::Ot,
        BisulfiteStrand::OB => OracleStrand::Ob,
        BisulfiteStrand::CTOT => OracleStrand::Ctot,
        BisulfiteStrand::CTOB => OracleStrand::Ctob,
    }
}

fn placement(alignment: &VerifiedAlignment) -> OraclePlacement {
    OraclePlacement {
        start: usize::try_from(alignment.interval().start()).expect("bounded start"),
        end: usize::try_from(alignment.interval().end()).expect("bounded end"),
        distance: alignment.distance().get(),
    }
}

#[test]
#[allow(clippy::too_many_lines, clippy::type_complexity)]
fn named_exact_conversion_reverse_indel_and_boundary_ground_truth() {
    let cases: &[(
        &[u8],
        &[u8],
        BisulfiteStrand,
        usize,
        usize,
        u64,
        usize,
        usize,
        u64,
        &str,
    )] = &[
        (
            b"GGACCTAA",
            b"ACCT",
            BisulfiteStrand::OT,
            0,
            4,
            0,
            2,
            6,
            0,
            "4M",
        ),
        (
            b"GGACCGAA",
            b"ATTG",
            BisulfiteStrand::OT,
            0,
            4,
            0,
            2,
            6,
            0,
            "4M",
        ),
        (
            b"TTAACGAA",
            b"CGTT",
            BisulfiteStrand::OB,
            0,
            4,
            0,
            2,
            6,
            0,
            "4M",
        ),
        (
            b"AGGT",
            b"ATTT",
            BisulfiteStrand::OB,
            0,
            4,
            0,
            0,
            4,
            0,
            "4M",
        ),
        (
            b"ACGT",
            b"TACGT",
            BisulfiteStrand::OT,
            1,
            5,
            1,
            0,
            4,
            1,
            "1I4M",
        ),
        (
            b"GACTGA",
            b"GACGA",
            BisulfiteStrand::OT,
            3,
            5,
            1,
            0,
            6,
            1,
            "3M1D2M",
        ),
        (
            b"ACGTGG",
            b"TACGT",
            BisulfiteStrand::OT,
            1,
            5,
            1,
            0,
            4,
            1,
            "1I4M",
        ),
    ];
    for &(
        reference_raw,
        query,
        strand,
        seed_start,
        seed_end,
        budget,
        start,
        end,
        distance,
        cigar,
    ) in cases
    {
        let reference = build_reference(&[(b"chr", reference_raw)]);
        let candidates = fixed_candidates(&reference, query, strand, seed_start, seed_end);
        let ordinal = candidates
            .anchors()
            .iter()
            .enumerate()
            .find_map(|(ordinal, anchor)| {
                let oracle = candidate_window_best(
                    reference_raw,
                    query,
                    oracle_strand(anchor.strand()),
                    diagonal(&candidates, ordinal),
                    budget,
                );
                oracle
                    .placements
                    .iter()
                    .any(|placement| placement.start == start && placement.end == end)
                    .then_some(ordinal)
            })
            .expect("a true-origin candidate window retains the named placement");
        let result = extend_candidate_window(
            &reference,
            &candidates,
            to_u64(ordinal),
            EditDistance::new(budget),
            ExtensionLimits::MAX,
        )
        .expect("extension succeeds");
        let expected = candidate_window_best(
            reference_raw,
            query,
            oracle_strand(strand),
            diagonal(&candidates, ordinal),
            budget,
        );
        assert_eq!(result.window().start(), to_u64(expected.start));
        assert_eq!(result.window().end(), to_u64(expected.end));
        assert_eq!(
            result
                .alignments()
                .iter()
                .map(placement)
                .collect::<Vec<_>>(),
            expected.placements
        );
        let selected = result
            .alignments()
            .iter()
            .find(|alignment| {
                alignment.interval().start() == to_u64(start)
                    && alignment.interval().end() == to_u64(end)
            })
            .expect("named placement is retained");
        assert_eq!(selected.distance().get(), distance);
        assert_eq!(selected.cigar().to_string(), cigar);
    }
}

#[test]
fn exhaustive_canonical_exact_windows_equal_independent_oracle() {
    let references = enumerate_strings(&ALPHABET, 4);
    let queries = enumerate_strings(&ALPHABET, 3);
    let mut windows = 0_u64;
    let mut evaluated_intervals = 0_u64;
    for reference_raw in &references {
        let reference = build_reference(&[(b"chr", reference_raw)]);
        for query in &queries {
            for strand in [BisulfiteStrand::OT, BisulfiteStrand::OB] {
                let candidates = fixed_candidates(&reference, query, strand, 0, query.len());
                for (ordinal, anchor) in candidates.anchors().iter().enumerate() {
                    let result = extend_candidate_window(
                        &reference,
                        &candidates,
                        to_u64(ordinal),
                        EditDistance::new(0),
                        ExtensionLimits::MAX,
                    )
                    .expect("bounded extension succeeds");
                    let oracle = candidate_window_best(
                        reference_raw,
                        query,
                        oracle_strand(anchor.strand()),
                        diagonal(&candidates, ordinal),
                        0,
                    );
                    assert_eq!(
                        (result.window().start(), result.window().end()),
                        (to_u64(oracle.start), to_u64(oracle.end))
                    );
                    assert_eq!(result.metrics().interval_alignments(), oracle.intervals);
                    assert_eq!(result.metrics().aggregate_dp_cells(), oracle.dp_cells);
                    assert_eq!(
                        result
                            .alignments()
                            .iter()
                            .map(placement)
                            .collect::<Vec<_>>(),
                        oracle.placements,
                        "reference={}, query={}, strand={strand:?}, diagonal={}",
                        String::from_utf8_lossy(reference_raw),
                        String::from_utf8_lossy(query),
                        diagonal(&candidates, ordinal)
                    );
                    windows += 1;
                    evaluated_intervals += oracle.intervals;
                }
            }
        }
    }
    assert_eq!(windows, 3_356);
    assert_eq!(evaluated_intervals, 3_356);
}

#[test]
fn direct_endpoint_sweeps_reduce_short_and_long_filter_work_exactly() {
    let short_reference = build_reference(&[(b"chr", b"GGACGTACGTTT")]);
    let short_query = b"ACGTACGT";
    let short_candidates = fixed_candidates(
        &short_reference,
        short_query,
        BisulfiteStrand::OT,
        0,
        short_query.len(),
    );
    let short_ordinal = short_candidates
        .anchors()
        .iter()
        .enumerate()
        .find_map(|(ordinal, _)| (diagonal(&short_candidates, ordinal) == 2).then_some(ordinal))
        .expect("exact short origin is a candidate");
    let short = extend_candidate_window(
        &short_reference,
        &short_candidates,
        to_u64(short_ordinal),
        EditDistance::new(2),
        ExtensionLimits::MAX,
    )
    .expect("short direct sweep succeeds");
    let short_expected = expected_sweep_metrics(short.window().len(), to_u64(short_query.len()), 2);
    assert_eq!(short.metrics().distance_sweeps(), short_expected.0);
    assert_eq!(
        short.metrics().distance_filter_updates(),
        short_expected.1,
        "one-word Myers performs one update per scanned reference base"
    );
    assert!(short.metrics().interval_alignments() > short.metrics().distance_sweeps());
    assert_eq!(
        short.metrics().traceback_alignments(),
        short.metrics().passing_alignments()
    );

    let mut long_reference_raw = vec![b'T'; 3];
    long_reference_raw.extend(std::iter::repeat_n(b'G', 80));
    long_reference_raw.extend(std::iter::repeat_n(b'A', 3));
    let mut long_query = vec![b'G'; 80];
    long_query[40] = b'C';
    let long_reference = build_reference(&[(b"chr", &long_reference_raw)]);
    let long_candidates =
        fixed_candidates(&long_reference, &long_query, BisulfiteStrand::OT, 0, 20);
    let long_ordinal = long_candidates
        .anchors()
        .iter()
        .enumerate()
        .find_map(|(ordinal, _)| (diagonal(&long_candidates, ordinal) == 3).then_some(ordinal))
        .expect("long true origin is a candidate");
    let long = extend_candidate_window(
        &long_reference,
        &long_candidates,
        to_u64(long_ordinal),
        EditDistance::new(2),
        ExtensionLimits::MAX,
    )
    .expect("long scalar direct sweep succeeds");
    let long_expected = expected_sweep_metrics(long.window().len(), to_u64(long_query.len()), 2);
    assert_eq!(long.metrics().distance_sweeps(), long_expected.0);
    assert_eq!(
        long.metrics().distance_filter_updates(),
        long_expected.1,
        "two-word Myers performs one update per scanned reference base"
    );
    assert!(long.metrics().interval_alignments() > long.metrics().distance_sweeps());
    assert!(long.metrics().traceback_alignments() < long.metrics().interval_alignments());
    assert_eq!(
        long.metrics().traceback_alignments(),
        long.metrics().passing_alignments()
    );
    let true_origin = long
        .alignments()
        .iter()
        .find(|alignment| alignment.interval().start() == 3 && alignment.interval().end() == 83)
        .expect("true long origin remains in the complete best set");
    assert_eq!(true_origin.distance().get(), 1);
    assert_eq!(true_origin.cigar().to_string(), "80M");
}

fn expected_sweep_metrics(window_bases: u64, query_bases: u64, budget: u64) -> (u64, u64) {
    let minimum_length = query_bases.saturating_sub(budget).max(1);
    let maximum_length = query_bases.saturating_add(budget).min(window_bases);
    let mut sweeps = 0_u64;
    let mut reference_updates = 0_u64;
    for start in 0..window_bases {
        let remaining = window_bases - start;
        if remaining < minimum_length {
            break;
        }
        sweeps += 1;
        reference_updates += maximum_length.min(remaining);
    }
    (sweeps, reference_updates)
}

#[test]
fn whole_contig_oracle_proves_seedless_nonempty_placement_policy() {
    let with_match = whole_contig_best(b"TACG", b"AC", OracleStrand::Ot, 2);
    assert_eq!(with_match.best_distance, Some(0));
    assert!(with_match.placements.contains(&OraclePlacement {
        start: 1,
        end: 3,
        distance: 0,
    }));
    let empty = whole_contig_best(b"", b"AC", OracleStrand::Ot, 2);
    assert_eq!(empty.best_distance, None);
    assert!(empty.placements.is_empty());
    assert_eq!(empty.intervals, 0);
}

#[test]
fn tied_placements_and_best_result_cap_are_explicit() {
    let reference = build_reference(&[(b"chr", b"ACAC")]);
    let candidates = fixed_candidates(&reference, b"AC", BisulfiteStrand::OT, 0, 2);
    let first = extend_candidate_window(
        &reference,
        &candidates,
        0,
        EditDistance::new(0),
        ExtensionLimits::MAX,
    )
    .expect("extension succeeds");
    assert_eq!(
        first.alignments().len(),
        1,
        "one exact seed anchor owns one local window"
    );

    let error = extend_candidate_window(
        &reference,
        &candidates,
        0,
        EditDistance::new(0),
        ExtensionLimits::new(u64::MAX, u64::MAX, u64::MAX, 0),
    )
    .expect_err("zero result cap rejects the complete result");
    assert!(matches!(
        error,
        ExtensionError::LimitExceeded {
            counter: ExtensionCounter::BestAlignments,
            requested: 1,
            maximum: 0,
        }
    ));
}

#[test]
fn public_error_sources_and_limit_getters_are_stable() {
    let limits = ExtensionLimits::new(1, 2, 3, 4);
    assert_eq!(limits.max_window_bases(), 1);
    assert_eq!(limits.max_interval_alignments(), 2);
    assert_eq!(limits.max_aggregate_dp_cells(), 3);
    assert_eq!(limits.max_best_alignments(), 4);
    let error = ExtensionError::CandidateOrdinalOutOfBounds {
        ordinal: 7,
        candidate_count: 2,
    };
    assert_eq!(
        error.to_string(),
        "candidate ordinal 7 is outside candidate count 2"
    );
    assert!(error.source().is_none());
}
