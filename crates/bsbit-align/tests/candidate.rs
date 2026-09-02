//! Independent scientific and API tests for Level 2C fixed-seed candidates.

#[path = "support/candidate_oracle.rs"]
mod candidate_oracle;

use candidate_oracle::{
    CANONICAL, OracleAnchor, OracleContig, OracleEvidence, OracleMetrics, OracleRequest,
    OracleSnapshot, OracleStrand, REFERENCE_ALPHABET, candidate_snapshot, enumerate_strings,
    group_evidence, reverse_complement, to_u64,
};

use std::collections::HashSet;
use std::error::Error as _;
use std::num::NonZeroU64;
use std::thread;

use bsbit_align::search::candidate::{
    CandidateAnchor, CandidateDiagonal, CandidateError, CandidateLimits, CandidateSet,
    candidates_for_fixed_seeds,
};
use bsbit_align::search::fixed_seed::{
    FixedSeedPlan, FixedSeedRequest, SeedPlanError, SeedPlanLimits,
};
use bsbit_core::bisulfite::BisulfiteStrand;
use bsbit_core::coordinate::{
    CoordinateDomain, CoordinateError, CoordinateOperation, CoordinateShift, QueryInterval,
    QueryLength,
};
use bsbit_core::sequence::{NormalizedSequence, normalize_dna};
use bsbit_index::reference::{
    ContigInput, ReferenceBuildLimits, ReferenceIndex, ReferenceQueryError, ReferenceQueryLimits,
};

const UNBOUNDED_PLAN: SeedPlanLimits = SeedPlanLimits::new(u64::MAX, u64::MAX);
const UNBOUNDED_CANDIDATES: CandidateLimits = CandidateLimits::new(u64::MAX, u64::MAX);

fn normalized(raw: &[u8]) -> NormalizedSequence {
    normalize_dna(raw).expect("test DNA contains only A/C/G/T/N")
}

fn build_catalog(catalog: &[(&[u8], &[u8])]) -> ReferenceIndex {
    let inputs = catalog
        .iter()
        .map(|(name, sequence)| ContigInput::new(name.to_vec(), normalized(sequence)))
        .collect();
    ReferenceIndex::build(inputs, ReferenceBuildLimits::MAX)
        .expect("bounded test reference should build")
}

fn oracle_catalog<'a>(catalog: &'a [(&'a [u8], &'a [u8])]) -> Vec<OracleContig<'a>> {
    catalog
        .iter()
        .map(|(name, sequence)| OracleContig { name, sequence })
        .collect()
}

const fn implementation_strand(strand: OracleStrand) -> BisulfiteStrand {
    match strand {
        OracleStrand::Ot => BisulfiteStrand::OT,
        OracleStrand::Ob => BisulfiteStrand::OB,
        OracleStrand::Ctot => BisulfiteStrand::CTOT,
        OracleStrand::Ctob => BisulfiteStrand::CTOB,
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

fn implementation_requests(
    query_length: usize,
    requests: &[OracleRequest],
) -> Vec<FixedSeedRequest> {
    let length = QueryLength::new(to_u64(query_length));
    requests
        .iter()
        .map(|request| {
            let interval = QueryInterval::new(to_u64(request.start), to_u64(request.end), length)
                .expect("oracle request fits its query");
            FixedSeedRequest::new(implementation_strand(request.strand), interval)
        })
        .collect()
}

fn build_plan(query: &[u8], requests: &[OracleRequest]) -> FixedSeedPlan {
    let requests = implementation_requests(query.len(), requests);
    FixedSeedPlan::new(normalized(query), &requests, UNBOUNDED_PLAN)
        .expect("bounded valid seed plan should build")
}

fn produce(
    index: &ReferenceIndex,
    plan: &FixedSeedPlan,
    query_limits: ReferenceQueryLimits,
    candidate_limits: CandidateLimits,
) -> CandidateSet {
    candidates_for_fixed_seeds(index, plan, query_limits, candidate_limits)
        .expect("bounded candidate generation should succeed")
}

fn signed_diagonal(anchor: &CandidateAnchor) -> i128 {
    match anchor.diagonal().shift() {
        CoordinateShift::Zero => 0,
        CoordinateShift::Forward(value) => i128::from(value.get()),
        CoordinateShift::Backward(value) => -i128::from(value.get()),
    }
}

fn anchor_snapshot(anchor: &CandidateAnchor) -> OracleAnchor {
    OracleAnchor {
        contig_ordinal: anchor.contig().ordinal(),
        diagonal: signed_diagonal(anchor),
        strand: oracle_strand(anchor.strand()),
        support: anchor.support().get(),
    }
}

fn implementation_snapshot(candidates: &CandidateSet) -> OracleSnapshot {
    let metrics = candidates.metrics();
    OracleSnapshot {
        anchors: candidates.anchors().iter().map(anchor_snapshot).collect(),
        metrics: OracleMetrics {
            request_count: metrics.request_count(),
            total_seed_bases: metrics.total_seed_bases(),
            total_exact_hits: metrics.total_exact_hits(),
            matched_intervals: metrics.matched_intervals(),
            unique_candidates: metrics.unique_candidates(),
            duplicate_evidence: metrics.duplicate_evidence(),
            maximum_support: metrics.maximum_support(),
            zero_hit_requests: metrics.zero_hit_requests(),
        },
    }
}

fn assert_implementation_equals_oracle(
    catalog: &[(&[u8], &[u8])],
    query: &[u8],
    requests: &[OracleRequest],
) -> OracleSnapshot {
    let index = build_catalog(catalog);
    let plan = build_plan(query, requests);
    let actual = implementation_snapshot(&produce(
        &index,
        &plan,
        ReferenceQueryLimits::MAX,
        UNBOUNDED_CANDIDATES,
    ));
    let expected = candidate_snapshot(&oracle_catalog(catalog), query, requests);
    assert_eq!(actual, expected);
    actual
}

fn all_requests(query_length: usize) -> Vec<OracleRequest> {
    let mut requests = Vec::new();
    for strand in OracleStrand::ALL {
        for start in 0..query_length {
            for end in start + 1..=query_length {
                requests.push(OracleRequest::new(strand, start, end));
            }
        }
    }
    requests
}

#[test]
fn named_negative_boundary_candidate_retains_later_1i4m_evidence() {
    let catalog = [(&b"contig"[..], &b"AGGT"[..])];
    let query = b"TAGGT";
    let requests = [OracleRequest::new(OracleStrand::Ot, 1, 5)];
    let snapshot = assert_implementation_equals_oracle(&catalog, query, &requests);
    assert_eq!(
        snapshot.anchors,
        vec![OracleAnchor {
            contig_ordinal: 0,
            diagonal: -1,
            strand: OracleStrand::Ot,
            support: 1,
        }]
    );
    assert_eq!(snapshot.metrics.total_exact_hits, 1);
}

#[test]
fn all_four_nonpalindromic_strands_use_a_or_l_minus_b() {
    let catalog = [(&b"nonpal"[..], &b"ACGTACTGCA"[..])];
    let query = b"AGTAATACGG";
    let requests = [
        OracleRequest::new(OracleStrand::Ot, 1, 4),
        OracleRequest::new(OracleStrand::Ob, 5, 8),
        OracleRequest::new(OracleStrand::Ctot, 5, 8),
        OracleRequest::new(OracleStrand::Ctob, 1, 4),
    ];
    let snapshot = assert_implementation_equals_oracle(&catalog, query, &requests);
    for strand in OracleStrand::ALL {
        let expected_diagonal = i128::from(!strand.is_reverse());
        assert!(snapshot.anchors.contains(&OracleAnchor {
            contig_ordinal: 0,
            diagonal: expected_diagonal,
            strand,
            support: 1,
        }));
    }
}

#[test]
fn projection_only_false_equalities_remain_unverified_candidate_evidence() {
    let catalog = [(&b"ct"[..], &b"C"[..]), (&b"ga"[..], &b"G"[..])];
    let query = b"TA";
    let requests = [
        OracleRequest::new(OracleStrand::Ot, 0, 1),
        OracleRequest::new(OracleStrand::Ctob, 1, 2),
    ];
    assert_ne!(catalog[0].1[0], query[0]);
    assert_ne!(catalog[1].1[0], query[1]);
    let snapshot = assert_implementation_equals_oracle(&catalog, query, &requests);
    assert!(snapshot.anchors.contains(&OracleAnchor {
        contig_ordinal: 0,
        diagonal: 0,
        strand: OracleStrand::Ot,
        support: 1,
    }));
    assert!(snapshot.anchors.contains(&OracleAnchor {
        contig_ordinal: 1,
        diagonal: -1,
        strand: OracleStrand::Ctob,
        support: 1,
    }));
}

#[test]
fn distinct_seed_support_deduplicates_by_contig_diagonal_and_strand() {
    let catalog = [(&b"repeat"[..], &b"ACGTACGT"[..])];
    let query = b"ACGT";
    let requests = [
        OracleRequest::new(OracleStrand::Ot, 0, 2),
        OracleRequest::new(OracleStrand::Ot, 2, 4),
    ];
    let snapshot = assert_implementation_equals_oracle(&catalog, query, &requests);
    assert_eq!(
        snapshot.anchors,
        vec![
            OracleAnchor {
                contig_ordinal: 0,
                diagonal: 0,
                strand: OracleStrand::Ot,
                support: 2,
            },
            OracleAnchor {
                contig_ordinal: 0,
                diagonal: 4,
                strand: OracleStrand::Ot,
                support: 2,
            },
        ]
    );
    assert_eq!(snapshot.metrics.total_exact_hits, 4);
    assert_eq!(snapshot.metrics.unique_candidates, 2);
    assert_eq!(snapshot.metrics.duplicate_evidence, 2);
    assert_eq!(snapshot.metrics.maximum_support, 2);
}

#[test]
fn equal_local_coordinates_on_different_contigs_and_strands_never_merge() {
    let catalog = [(&b"left"[..], &b"AT"[..]), (&b"right"[..], &b"AT"[..])];
    let query = b"AT";
    let requests = OracleStrand::ALL.map(|strand| OracleRequest::new(strand, 0, 2));
    let snapshot = assert_implementation_equals_oracle(&catalog, query, &requests);
    assert_eq!(snapshot.anchors.len(), 8);
    for contig_ordinal in 0..=1 {
        for strand in OracleStrand::ALL {
            assert!(snapshot.anchors.contains(&OracleAnchor {
                contig_ordinal,
                diagonal: 0,
                strand,
                support: 1,
            }));
        }
    }
}

#[test]
fn all_n_and_fragmented_reference_barriers_preserve_run_boundaries() {
    let catalog = [
        (&b"all-n"[..], &b"NNN"[..]),
        (&b"fragmented"[..], &b"ACNNTGN"[..]),
        (&b"single"[..], &b"A"[..]),
    ];
    let query = b"AC";
    let requests = all_requests(query.len());
    let snapshot = assert_implementation_equals_oracle(&catalog, query, &requests);
    assert!(
        snapshot
            .anchors
            .iter()
            .all(|anchor| anchor.contig_ordinal != 0)
    );
}

#[test]
fn oracle_grouping_is_permutation_idempotent_and_rejects_duplicate_request_evidence() {
    let evidence = vec![
        OracleEvidence {
            contig_ordinal: 1,
            diagonal: -2,
            strand: OracleStrand::Ob,
            request_ordinal: 2,
        },
        OracleEvidence {
            contig_ordinal: 0,
            diagonal: 3,
            strand: OracleStrand::Ot,
            request_ordinal: 1,
        },
        OracleEvidence {
            contig_ordinal: 1,
            diagonal: -2,
            strand: OracleStrand::Ob,
            request_ordinal: 0,
        },
    ];
    let expected = group_evidence(evidence.clone());
    let mut reversed = evidence.clone();
    reversed.reverse();
    assert_eq!(group_evidence(reversed), expected);
    let mut rotated = evidence;
    rotated.rotate_left(1);
    assert_eq!(group_evidence(rotated), expected);
}

#[test]
#[should_panic(expected = "one request produced duplicate evidence")]
fn oracle_rejects_duplicate_evidence_from_one_request_for_one_key() {
    let duplicate = OracleEvidence {
        contig_ordinal: 0,
        diagonal: 0,
        strand: OracleStrand::Ot,
        request_ordinal: 0,
    };
    let _ = group_evidence(vec![duplicate, duplicate]);
}

#[test]
fn signed_diagonal_helpers_cover_full_magnitude_and_mathematical_order() {
    let ordered = [
        CandidateDiagonal::before_contig(
            NonZeroU64::new(u64::MAX).expect("maximum magnitude is nonzero"),
        ),
        CandidateDiagonal::before_contig(NonZeroU64::new(1).expect("one is nonzero")),
        CandidateDiagonal::at_or_after_contig(0),
        CandidateDiagonal::at_or_after_contig(u64::MAX),
    ];
    assert!(ordered.windows(2).all(|pair| pair[0] < pair[1]));
    assert!(ordered[0].is_before_contig());
    assert!(ordered[1].is_before_contig());
    assert!(!ordered[2].is_before_contig());
    assert!(!ordered[3].is_before_contig());
    assert_eq!(ordered[0].magnitude(), u64::MAX);
    assert_eq!(ordered[1].magnitude(), 1);
    assert_eq!(ordered[2].magnitude(), 0);
    assert_eq!(ordered[3].magnitude(), u64::MAX);
    assert_eq!(ordered[0].shift().magnitude(), u64::MAX);
    assert_eq!(ordered[1].shift().magnitude(), 1);
    assert_eq!(ordered[2].shift(), CoordinateShift::Zero);
    assert_eq!(ordered[3].shift().magnitude(), u64::MAX);
    let hashed = ordered.into_iter().collect::<HashSet<_>>();
    assert_eq!(hashed.len(), 4);
    assert!(hashed.contains(&CandidateDiagonal::at_or_after_contig(0)));
}

#[test]
fn zero_requests_empty_query_and_empty_result_retain_exact_owners() {
    let index = build_catalog(&[(&b"unknown"[..], &b"N"[..])]);
    let plan = build_plan(b"", &[]);
    let query_id = plan.query_instance_id();
    let reference_id = index.instance_id();
    let candidates = produce(
        &index,
        &plan,
        ReferenceQueryLimits::MAX,
        CandidateLimits::new(0, 0),
    );
    assert!(candidates.anchors().is_empty());
    assert!(candidates.query().is_empty());
    assert!(candidates.belongs_to_query(&query_id));
    assert!(candidates.belongs_to_reference(&reference_id));
    assert!(
        candidates
            .query_instance_id()
            .is_same_instance(&plan.query_instance_id())
    );
    assert!(
        candidates
            .reference_instance_id()
            .is_same_instance(&index.instance_id())
    );
    assert_eq!(
        implementation_snapshot(&candidates)
            .metrics
            .total_exact_hits,
        0
    );
}

#[test]
fn copied_plans_share_query_owner_but_equal_rebuilds_do_not() {
    let requests = [OracleRequest::new(OracleStrand::Ot, 0, 2)];
    let original = build_plan(b"AC", &requests);
    let copied = original
        .try_clone()
        .expect("small plan copy should allocate");
    let rebuilt = build_plan(b"AC", &requests);
    assert!(
        original
            .query_instance_id()
            .is_same_instance(&copied.query_instance_id())
    );
    assert!(
        !original
            .query_instance_id()
            .is_same_instance(&rebuilt.query_instance_id())
    );
    assert_eq!(original.query().to_ascii(), rebuilt.query().to_ascii());
    assert_eq!(original.requests().len(), copied.requests().len());
}

#[test]
fn equal_content_reference_and_query_rebuilds_have_distinct_owners_equal_semantics() {
    let catalog = [(&b"same"[..], &b"ACGTACGT"[..])];
    let requests = all_requests(4);
    let left_index = build_catalog(&catalog);
    let right_index = build_catalog(&catalog);
    let left_plan = build_plan(b"ACGT", &requests);
    let right_plan = build_plan(b"ACGT", &requests);
    let left = produce(
        &left_index,
        &left_plan,
        ReferenceQueryLimits::MAX,
        UNBOUNDED_CANDIDATES,
    );
    let right = produce(
        &right_index,
        &right_plan,
        ReferenceQueryLimits::MAX,
        UNBOUNDED_CANDIDATES,
    );
    assert!(!left.belongs_to_reference(&right_index.instance_id()));
    assert!(!left.belongs_to_query(&right_plan.query_instance_id()));
    assert!(
        !left
            .reference_instance_id()
            .is_same_instance(&right.reference_instance_id())
    );
    assert!(
        !left
            .query_instance_id()
            .is_same_instance(&right.query_instance_id())
    );
    assert_eq!(
        implementation_snapshot(&left),
        implementation_snapshot(&right)
    );
}

#[test]
fn plan_validation_priority_and_boundary_lengths_are_explicit() {
    let outside_interval =
        QueryInterval::new(2, 2, QueryLength::new(2)).expect("foreign boundary is valid");
    let outside = FixedSeedRequest::new(BisulfiteStrand::OT, outside_interval);
    assert_eq!(
        FixedSeedPlan::new(normalized(b"A"), &[outside], UNBOUNDED_PLAN)
            .expect_err("actual-query bounds precede empty-seed validation"),
        SeedPlanError::InvalidInterval {
            request_ordinal: 0,
            source: CoordinateError::OutOfBounds {
                domain: CoordinateDomain::Query,
                operation: CoordinateOperation::IntervalConstruction,
                start: 2,
                end: 2,
                length: 1,
            },
        }
    );

    let empty_interval =
        QueryInterval::new(0, 0, QueryLength::new(1)).expect("empty interval is typed");
    let empty = FixedSeedRequest::new(BisulfiteStrand::OT, empty_interval);
    assert_eq!(
        FixedSeedPlan::new(normalized(b"A"), &[empty], UNBOUNDED_PLAN)
            .expect_err("empty seed is rejected after actual-query validation"),
        SeedPlanError::EmptySeed {
            request_ordinal: 0,
            interval: empty_interval,
        }
    );
    assert_eq!(
        FixedSeedPlan::new(normalized(b"A"), &[empty], SeedPlanLimits::new(0, u64::MAX))
            .expect_err("request cap precedes per-request validation"),
        SeedPlanError::RequestLimitExceeded {
            requested: 1,
            maximum: 0,
        }
    );

    let n_seed = FixedSeedRequest::new(
        BisulfiteStrand::OT,
        QueryInterval::new(1, 3, QueryLength::new(3)).expect("interval fits"),
    );
    assert_eq!(
        FixedSeedPlan::new(normalized(b"ACN"), &[n_seed], SeedPlanLimits::new(1, 0))
            .expect_err("N validation precedes total-seed limit"),
        SeedPlanError::UnsearchableBase {
            request_ordinal: 0,
            query_offset: 2,
        }
    );

    for (query, offset) in [(b"NAA".as_slice(), 0_u64), (b"ANA", 1), (b"AAN", 2)] {
        let request = FixedSeedRequest::new(
            BisulfiteStrand::OT,
            QueryInterval::new(0, 3, QueryLength::new(3)).expect("interval fits"),
        );
        assert_eq!(
            FixedSeedPlan::new(normalized(query), &[request], UNBOUNDED_PLAN)
                .expect_err("N at every seed position is rejected"),
            SeedPlanError::UnsearchableBase {
                request_ordinal: 0,
                query_offset: offset,
            }
        );
    }

    let two_base_interval = QueryInterval::new(0, 2, QueryLength::new(2)).expect("interval fits");
    let two_base = FixedSeedRequest::new(BisulfiteStrand::OT, two_base_interval);
    assert_eq!(
        FixedSeedPlan::new(normalized(b"AA"), &[two_base], SeedPlanLimits::new(1, 1))
            .expect_err("first exact prefix exceeds the seed-base cap"),
        SeedPlanError::TotalSeedBasesLimitExceeded {
            request_ordinal: 0,
            requested: 2,
            maximum: 1,
        }
    );

    let duplicate_interval = QueryInterval::new(0, 1, QueryLength::new(1)).expect("interval fits");
    let duplicate = FixedSeedRequest::new(BisulfiteStrand::OT, duplicate_interval);
    assert_eq!(
        FixedSeedPlan::new(normalized(b"A"), &[duplicate, duplicate], UNBOUNDED_PLAN)
            .expect_err("an exact duplicate request must be rejected"),
        SeedPlanError::DuplicateRequest {
            strand: BisulfiteStrand::OT,
            interval: duplicate_interval,
        }
    );

    for length in [1_usize, 17, 18, 29, 30] {
        let query = vec![b'A'; length];
        let plan = build_plan(&query, &[OracleRequest::new(OracleStrand::Ot, 0, length)]);
        assert_eq!(plan.query_length().get(), to_u64(length));
        assert_eq!(plan.metrics().request_count(), 1);
        assert_eq!(plan.metrics().total_seed_bases(), to_u64(length));
    }
}

#[test]
fn reachable_seed_plan_error_displays_and_sources_are_exact() {
    let foreign = QueryInterval::new(2, 2, QueryLength::new(2)).expect("interval fits");
    let empty = QueryInterval::new(0, 0, QueryLength::new(1)).expect("interval fits");
    let full_one = QueryInterval::new(0, 1, QueryLength::new(1)).expect("interval fits");
    let full_two = QueryInterval::new(0, 2, QueryLength::new(2)).expect("interval fits");
    let full_three = QueryInterval::new(0, 3, QueryLength::new(3)).expect("interval fits");
    let errors = vec![
        (
            FixedSeedPlan::new(
                normalized(b"A"),
                &[FixedSeedRequest::new(BisulfiteStrand::OT, empty)],
                SeedPlanLimits::new(0, u64::MAX),
            )
            .expect_err("request limit must reject"),
            "seed request count 1 exceeds configured maximum 0",
        ),
        (
            FixedSeedPlan::new(
                normalized(b"A"),
                &[FixedSeedRequest::new(BisulfiteStrand::OT, foreign)],
                UNBOUNDED_PLAN,
            )
            .expect_err("foreign interval must reject"),
            "seed request 0 is invalid for the actual query: Query interval [2, 2) is outside length 1 during IntervalConstruction",
        ),
        (
            FixedSeedPlan::new(
                normalized(b"A"),
                &[FixedSeedRequest::new(BisulfiteStrand::OT, empty)],
                UNBOUNDED_PLAN,
            )
            .expect_err("empty seed must reject"),
            "seed request 0 has empty interval query:[0,0)",
        ),
        (
            FixedSeedPlan::new(
                normalized(b"AAN"),
                &[FixedSeedRequest::new(BisulfiteStrand::OT, full_three)],
                UNBOUNDED_PLAN,
            )
            .expect_err("N seed must reject"),
            "seed request 0 contains unsearchable N at absolute query offset 2",
        ),
        (
            FixedSeedPlan::new(
                normalized(b"AA"),
                &[FixedSeedRequest::new(BisulfiteStrand::OT, full_two)],
                SeedPlanLimits::new(1, 1),
            )
            .expect_err("seed total must reject"),
            "seed request 0 raises prefix total to 2, exceeding 1",
        ),
        (
            FixedSeedPlan::new(
                normalized(b"A"),
                &[
                    FixedSeedRequest::new(BisulfiteStrand::OT, full_one),
                    FixedSeedRequest::new(BisulfiteStrand::OT, full_one),
                ],
                UNBOUNDED_PLAN,
            )
            .expect_err("duplicate request must reject"),
            "duplicate seed request OT query:[0,1)",
        ),
    ];
    for (error, display) in &errors {
        assert_eq!(error.to_string(), *display);
    }
    assert_eq!(
        errors[1].0.source().map(ToString::to_string),
        Some("Query interval [2, 2) is outside length 1 during IntervalConstruction".to_owned())
    );
    for index in [0, 2, 3, 4, 5] {
        assert!(errors[index].0.source().is_none());
    }
}

#[test]
fn seed_start_middle_end_and_full_query_are_distinct_requests() {
    let requests = [
        OracleRequest::new(OracleStrand::Ot, 0, 1),
        OracleRequest::new(OracleStrand::Ot, 2, 3),
        OracleRequest::new(OracleStrand::Ot, 4, 5),
        OracleRequest::new(OracleStrand::Ot, 0, 5),
    ];
    let plan = build_plan(b"ACGTA", &requests);
    assert_eq!(plan.metrics().request_count(), 4);
    assert_eq!(plan.metrics().total_seed_bases(), 8);
    assert_eq!(plan.requests().len(), 4);
}

#[test]
fn public_limit_request_and_metric_getters_report_exact_values() {
    let plan_limits = SeedPlanLimits::new(7, 11);
    assert_eq!(plan_limits.max_requests(), 7);
    assert_eq!(plan_limits.max_total_seed_bases(), 11);
    let candidate_limits = CandidateLimits::new(13, 17);
    assert_eq!(candidate_limits.max_total_exact_hits(), 13);
    assert_eq!(candidate_limits.max_unique_candidates(), 17);

    let interval = QueryInterval::new(0, 2, QueryLength::new(2)).expect("interval fits");
    let request = FixedSeedRequest::new(BisulfiteStrand::OT, interval);
    assert_eq!(request.strand(), BisulfiteStrand::OT);
    assert_eq!(request.interval(), interval);
    let plan = FixedSeedPlan::new(normalized(b"AC"), &[request], plan_limits)
        .expect("small plan fits explicit limits");
    assert_eq!(plan.metrics().query_bases(), 2);
    assert_eq!(plan.metrics().request_count(), 1);
    assert_eq!(plan.metrics().total_seed_bases(), 2);

    let index = build_catalog(&[(&b"one"[..], &b"AC"[..])]);
    let candidates = produce(&index, &plan, ReferenceQueryLimits::MAX, candidate_limits);
    let expected = candidate_snapshot(
        &[OracleContig {
            name: b"one",
            sequence: b"AC",
        }],
        b"AC",
        &[OracleRequest::new(OracleStrand::Ot, 0, 2)],
    );
    assert_eq!(implementation_snapshot(&candidates), expected);
    let metrics = candidates.metrics();
    assert_eq!(metrics.search_rank_operations(), 8);
    assert_eq!(metrics.locate_calls(), 1);
    assert_eq!(metrics.located_coordinates(), 1);
    assert_eq!(metrics.locate_lf_steps(), 0);
    assert_eq!(metrics.locate_rank_operations(), 0);
    assert_eq!(metrics.locate_interval_nodes(), 0);
    assert_eq!(metrics.candidate_key_materializations(), 1);
    assert_eq!(metrics.peak_request_candidate_keys(), 1);
}

#[test]
fn caller_hit_cap_wins_a_tie_and_aggregate_cap_wins_when_smaller() {
    let index = build_catalog(&[(&b"repeat"[..], &b"AAAA"[..])]);
    let plan = build_plan(b"A", &[OracleRequest::new(OracleStrand::Ot, 0, 1)]);

    let tie = candidates_for_fixed_seeds(
        &index,
        &plan,
        ReferenceQueryLimits::new(1, 3),
        CandidateLimits::new(3, u64::MAX),
    )
    .expect_err("four hits exceed both equal remaining caps");
    let seed_interval = QueryInterval::new(0, 1, QueryLength::new(1)).expect("interval fits");
    let expected_search_source = ReferenceQueryError::HitLimitExceeded {
        requested: 4,
        maximum: 3,
    };
    assert_eq!(
        tie,
        CandidateError::Search {
            request_ordinal: 0,
            strand: BisulfiteStrand::OT,
            interval: seed_interval,
            source: expected_search_source.clone(),
        }
    );
    assert_eq!(
        tie.source()
            .and_then(|source| source.downcast_ref::<ReferenceQueryError>()),
        Some(&expected_search_source)
    );
    assert_eq!(
        tie.to_string(),
        "candidate search failed for request 0 OT query:[0,1): exact hit count 4 exceeds configured maximum 3"
    );

    let aggregate = candidates_for_fixed_seeds(
        &index,
        &plan,
        ReferenceQueryLimits::new(1, 4),
        CandidateLimits::new(3, u64::MAX),
    )
    .expect_err("aggregate remaining capacity is smaller");
    assert_eq!(
        aggregate,
        CandidateError::AggregateHitLimitExceeded {
            accumulated: 0,
            request_hits: 4,
            requested: 4,
            maximum: 3,
        }
    );
    assert!(aggregate.source().is_none());
    assert_eq!(
        aggregate.to_string(),
        "aggregate exact hits 0 plus request count 4 is 4, exceeding 3"
    );
}

#[test]
fn zero_and_unique_candidate_limits_are_complete_result_gates() {
    let index = build_catalog(&[(&b"repeat"[..], &b"ACAC"[..])]);
    let empty_plan = build_plan(b"", &[]);
    assert!(
        candidates_for_fixed_seeds(
            &index,
            &empty_plan,
            ReferenceQueryLimits::MAX,
            CandidateLimits::new(0, 0),
        )
        .expect("zero caps admit a complete empty result")
        .anchors()
        .is_empty()
    );

    let plan = build_plan(b"AC", &[OracleRequest::new(OracleStrand::Ot, 0, 2)]);
    let error = candidates_for_fixed_seeds(
        &index,
        &plan,
        ReferenceQueryLimits::MAX,
        CandidateLimits::new(2, 1),
    )
    .expect_err("two complete unique candidates exceed a limit of one");
    assert_eq!(
        error,
        CandidateError::UniqueCandidateLimitExceeded {
            requested: 2,
            maximum: 1,
        }
    );
    assert!(error.source().is_none());
    assert_eq!(
        error.to_string(),
        "unique candidate count 2 exceeds configured maximum 1"
    );
}

#[test]
fn request_permutations_and_repeated_execution_preserve_exact_order() {
    let catalog = [
        (&b"alpha"[..], &b"ACGTCNTA"[..]),
        (&b"beta"[..], &b"GCA"[..]),
    ];
    let query = b"ACGT";
    let canonical = all_requests(query.len());
    let mut reversed = canonical.clone();
    reversed.reverse();
    let mut rotated = canonical.clone();
    rotated.rotate_left(7);

    let index = build_catalog(&catalog);
    let expected_oracle = candidate_snapshot(&oracle_catalog(&catalog), query, &canonical);
    let mut prior = None;
    for requests in [&canonical, &reversed, &rotated] {
        let plan = build_plan(query, requests);
        for _ in 0..3 {
            let observed = implementation_snapshot(&produce(
                &index,
                &plan,
                ReferenceQueryLimits::MAX,
                UNBOUNDED_CANDIDATES,
            ));
            assert_eq!(observed, expected_oracle);
            if let Some(previous) = &prior {
                assert_eq!(&observed, previous);
            }
            prior = Some(observed);
        }
    }
}

#[test]
fn single_request_partitions_recombine_to_the_whole_semantic_multiset() {
    let catalog = [(&b"repeat"[..], &b"ACGTACGT"[..])];
    let query = b"ACGT";
    let requests = all_requests(query.len());
    let index = build_catalog(&catalog);
    let whole_plan = build_plan(query, &requests);
    let whole = implementation_snapshot(&produce(
        &index,
        &whole_plan,
        ReferenceQueryLimits::MAX,
        UNBOUNDED_CANDIDATES,
    ));

    let mut partition_evidence = Vec::new();
    for (request_ordinal, request) in requests.iter().enumerate() {
        let plan = build_plan(query, &[*request]);
        let part = produce(
            &index,
            &plan,
            ReferenceQueryLimits::MAX,
            UNBOUNDED_CANDIDATES,
        );
        for anchor in part.anchors() {
            assert_eq!(anchor.support().get(), 1);
            partition_evidence.push(OracleEvidence {
                contig_ordinal: anchor.contig().ordinal(),
                diagonal: signed_diagonal(anchor),
                strand: oracle_strand(anchor.strand()),
                request_ordinal: to_u64(request_ordinal),
            });
        }
    }
    assert_eq!(group_evidence(partition_evidence), whole.anchors);
}

#[test]
fn appending_a_contig_preserves_prior_ordinal_keys_and_adds_only_the_new_ordinal() {
    let base_catalog = [(&b"base"[..], &b"ACGT"[..])];
    let appended_catalog = [(&b"base"[..], &b"ACGT"[..]), (&b"extra"[..], &b"TTTT"[..])];
    let query = b"ACGT";
    let requests = all_requests(query.len());
    let base = assert_implementation_equals_oracle(&base_catalog, query, &requests);
    let appended = assert_implementation_equals_oracle(&appended_catalog, query, &requests);
    let retained = appended
        .anchors
        .iter()
        .copied()
        .filter(|anchor| anchor.contig_ordinal == 0)
        .collect::<Vec<_>>();
    assert_eq!(retained, base.anchors);
    assert!(
        appended
            .anchors
            .iter()
            .filter(|anchor| anchor.contig_ordinal != 0)
            .all(|anchor| anchor.contig_ordinal == 1)
    );
}

#[test]
fn bounded_one_and_two_contig_catalogs_equal_the_independent_direct_scan() {
    let references = enumerate_strings(&REFERENCE_ALPHABET, 2)
        .into_iter()
        .filter(|sequence| !sequence.is_empty())
        .collect::<Vec<_>>();
    let queries = enumerate_strings(&CANONICAL, 2)
        .into_iter()
        .filter(|query| !query.is_empty())
        .collect::<Vec<_>>();
    let short_references = enumerate_strings(&REFERENCE_ALPHABET, 1)
        .into_iter()
        .filter(|sequence| !sequence.is_empty())
        .collect::<Vec<_>>();
    let mut catalog_count = 0_u64;
    let mut snapshot_count = 0_u64;
    let mut request_count = 0_u64;

    for sequence in &references {
        let catalog = [(&b"one"[..], sequence.as_slice())];
        compare_catalog_queries(&catalog, &queries, &mut snapshot_count, &mut request_count);
        catalog_count += 1;
    }
    for left in &short_references {
        for right in &short_references {
            let catalog = [
                (&b"left"[..], left.as_slice()),
                (&b"right"[..], right.as_slice()),
            ];
            compare_catalog_queries(&catalog, &queries, &mut snapshot_count, &mut request_count);
            catalog_count += 1;
        }
    }

    assert_eq!(catalog_count, 55);
    assert_eq!(snapshot_count, 1_100);
    assert_eq!(request_count, 11_440);
}

fn compare_catalog_queries(
    catalog: &[(&[u8], &[u8])],
    queries: &[Vec<u8>],
    snapshot_count: &mut u64,
    request_count: &mut u64,
) {
    let index = build_catalog(catalog);
    let views = oracle_catalog(catalog);
    for query in queries {
        let requests = all_requests(query.len());
        let plan = build_plan(query, &requests);
        let actual = implementation_snapshot(&produce(
            &index,
            &plan,
            ReferenceQueryLimits::MAX,
            UNBOUNDED_CANDIDATES,
        ));
        assert_eq!(actual, candidate_snapshot(&views, query, &requests));
        *snapshot_count += 1;
        *request_count += to_u64(requests.len());
    }
}

#[test]
fn reverse_interval_content_identity_is_independently_exhaustive() {
    let queries = enumerate_strings(&CANONICAL, 4);
    let mut cases = 0_u64;
    for query in queries {
        let reverse = reverse_complement(&query);
        for start in 0..query.len() {
            for end in start + 1..=query.len() {
                let reverse_start = query.len() - end;
                let reverse_end = query.len() - start;
                assert_eq!(
                    &reverse[reverse_start..reverse_end],
                    reverse_complement(&query[start..end])
                );
                cases += 1;
            }
        }
    }
    assert_eq!(cases, 2_996);
}

#[test]
fn eight_workers_with_rotated_and_reversed_arrival_orders_equal_serial() {
    let catalog = [
        (&b"alpha"[..], &b"ACGTCNTA"[..]),
        (&b"beta"[..], &b"GCAACGT"[..]),
    ];
    let query = b"ACGT";
    let requests = all_requests(query.len());
    let index = build_catalog(&catalog);
    let serial_plan = build_plan(query, &requests);
    let serial = implementation_snapshot(&produce(
        &index,
        &serial_plan,
        ReferenceQueryLimits::MAX,
        UNBOUNDED_CANDIDATES,
    ));

    thread::scope(|scope| {
        let mut shared_plan_workers = Vec::new();
        for _ in 0..8 {
            let index = &index;
            let serial_plan = &serial_plan;
            let serial = &serial;
            shared_plan_workers.push(scope.spawn(move || {
                let observed = implementation_snapshot(&produce(
                    index,
                    serial_plan,
                    ReferenceQueryLimits::MAX,
                    UNBOUNDED_CANDIDATES,
                ));
                assert_eq!(&observed, serial);
            }));
        }
        for worker in shared_plan_workers {
            worker.join().expect("shared-plan worker should not panic");
        }
    });

    thread::scope(|scope| {
        let mut workers = Vec::new();
        for worker_id in 0..8 {
            let index = &index;
            let serial = &serial;
            let mut arrival = requests.clone();
            if worker_id % 2 == 1 {
                arrival.reverse();
            }
            let arrival_len = arrival.len();
            arrival.rotate_left(worker_id % arrival_len);
            workers.push(scope.spawn(move || {
                let plan = build_plan(query, &arrival);
                let observed = implementation_snapshot(&produce(
                    index,
                    &plan,
                    ReferenceQueryLimits::MAX,
                    UNBOUNDED_CANDIDATES,
                ));
                assert_eq!(&observed, serial);
            }));
        }
        for worker in workers {
            worker.join().expect("candidate worker should not panic");
        }
    });
}
