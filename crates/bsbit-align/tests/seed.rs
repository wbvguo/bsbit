//! Independent scientific and API tests for Level 2D proof-seed scheduling.

#[path = "support/seed_oracle.rs"]
mod seed_oracle;

use seed_oracle::{
    OracleCertificate, OracleOutcome, OracleRequest, exhaustive_edit_footprint_stats,
    schedule as oracle_schedule,
};

use std::thread;

use bsbit_align::score::EditDistance;
use bsbit_align::search::candidate::{
    CandidateAnchor, CandidateError, CandidateLimits, SeedPlanLimits, candidates_for_fixed_seeds,
};
use bsbit_align::search::seed::{
    AdmissibleStrands, ProofSeedLimits, ProofSeedOutcome, schedule_proof_seeds,
};
use bsbit_core::bisulfite::BisulfiteStrand;
use bsbit_core::coordinate::CoordinateShift;
use bsbit_core::sequence::{NormalizedSequence, normalize_dna};
use bsbit_index::reference::{
    ContigInput, ReferenceBuildLimits, ReferenceIndex, ReferenceQueryError, ReferenceQueryLimits,
};

const UNBOUNDED_CANDIDATES: CandidateLimits = CandidateLimits::new(u64::MAX, u64::MAX);

#[derive(Clone, Debug, Eq, PartialEq)]
enum ImplementationOutcome {
    Certified {
        certificate: OracleCertificate,
        requests: Vec<OracleRequest>,
    },
    SeedlessFallbackRequired {
        query: Vec<u8>,
        strand_ranks: Vec<u8>,
        max_edit_distance: usize,
    },
    NoAlignmentWithinBudget {
        query: Vec<u8>,
        strand_ranks: Vec<u8>,
        max_edit_distance: usize,
        unknown_bases: usize,
    },
}

fn normalized(raw: &[u8]) -> NormalizedSequence {
    normalize_dna(raw).expect("test sequence is normalized A/C/G/T/N")
}

const fn strand_rank(strand: BisulfiteStrand) -> u8 {
    match strand {
        BisulfiteStrand::OT => 0,
        BisulfiteStrand::OB => 1,
        BisulfiteStrand::CTOT => 2,
        BisulfiteStrand::CTOB => 3,
    }
}

const fn strand_from_rank(rank: u8) -> BisulfiteStrand {
    match rank {
        0 => BisulfiteStrand::OT,
        1 => BisulfiteStrand::OB,
        2 => BisulfiteStrand::CTOT,
        3 => BisulfiteStrand::CTOB,
        _ => panic!("oracle rank is outside the four-strand table"),
    }
}

fn to_usize(value: u64) -> usize {
    usize::try_from(value).expect("bounded test value fits usize")
}

fn implementation_snapshot(
    query: &[u8],
    supplied_ranks: &[u8],
    budget: usize,
) -> ImplementationOutcome {
    let supplied = supplied_ranks
        .iter()
        .rev()
        .copied()
        .map(strand_from_rank)
        .collect::<Vec<_>>();
    let strands = AdmissibleStrands::new(&supplied).expect("oracle subset is distinct");
    match schedule_proof_seeds(
        normalized(query),
        strands,
        EditDistance::new(u64::try_from(budget).expect("budget fits u64")),
        ProofSeedLimits::MAX,
    )
    .expect("bounded schedule succeeds")
    {
        ProofSeedOutcome::Certified(plan) => {
            let certificate = plan.certificate();
            let block_count = to_usize(certificate.block_count());
            let blocks = seed_oracle::balanced_blocks(query.len(), block_count);
            ImplementationOutcome::Certified {
                certificate: OracleCertificate {
                    query_bases: to_usize(certificate.query_bases()),
                    max_edit_distance: to_usize(certificate.max_edit_distance().get()),
                    strand_count: to_usize(certificate.strand_count()),
                    blocks,
                    emitted_blocks: to_usize(certificate.emitted_blocks()),
                    omitted_unknown_blocks: to_usize(certificate.omitted_unknown_blocks()),
                    unknown_bases: to_usize(certificate.unknown_bases()),
                    total_seed_bases: to_usize(certificate.total_seed_bases()),
                    minimum_block_bases: to_usize(certificate.minimum_block_bases()),
                    maximum_block_bases: to_usize(certificate.maximum_block_bases()),
                },
                requests: plan
                    .fixed_plan()
                    .requests()
                    .iter()
                    .map(|request| OracleRequest {
                        strand_rank: strand_rank(request.strand()),
                        start: to_usize(request.interval().start()),
                        end: to_usize(request.interval().end()),
                    })
                    .collect(),
            }
        }
        ProofSeedOutcome::SeedlessFallbackRequired(fallback) => {
            ImplementationOutcome::SeedlessFallbackRequired {
                query: fallback.query().to_ascii(),
                strand_ranks: fallback.strands().iter().map(strand_rank).collect(),
                max_edit_distance: to_usize(fallback.max_edit_distance().get()),
            }
        }
        ProofSeedOutcome::NoAlignmentWithinBudget(proof) => {
            ImplementationOutcome::NoAlignmentWithinBudget {
                query: proof.query().to_ascii(),
                strand_ranks: proof.strands().iter().map(strand_rank).collect(),
                max_edit_distance: to_usize(proof.max_edit_distance().get()),
                unknown_bases: to_usize(proof.unknown_bases()),
            }
        }
    }
}

fn oracle_snapshot(query: &[u8], supplied_ranks: &[u8], budget: usize) -> ImplementationOutcome {
    match oracle_schedule(query, supplied_ranks, budget) {
        OracleOutcome::Certified {
            certificate,
            requests,
        } => ImplementationOutcome::Certified {
            certificate,
            requests,
        },
        OracleOutcome::SeedlessFallbackRequired {
            query,
            strand_ranks,
            max_edit_distance,
        } => ImplementationOutcome::SeedlessFallbackRequired {
            query,
            strand_ranks,
            max_edit_distance,
        },
        OracleOutcome::NoAlignmentWithinBudget {
            query,
            strand_ranks,
            max_edit_distance,
            unknown_bases,
        } => ImplementationOutcome::NoAlignmentWithinBudget {
            query,
            strand_ranks,
            max_edit_distance,
            unknown_bases,
        },
    }
}

fn rank_subset(mask: u8) -> Vec<u8> {
    (0_u8..4).filter(|rank| mask & (1 << rank) != 0).collect()
}

fn permutations(values: &[u8]) -> Vec<Vec<u8>> {
    fn visit(values: &mut [u8], start: usize, output: &mut Vec<Vec<u8>>) {
        if start == values.len() {
            output.push(values.to_vec());
            return;
        }
        for index in start..values.len() {
            values.swap(start, index);
            visit(values, start + 1, output);
            values.swap(start, index);
        }
    }

    let mut owned = values.to_vec();
    let mut output = Vec::new();
    visit(&mut owned, 0, &mut output);
    output
}

fn query_from_n_mask(length: usize, mask: usize) -> Vec<u8> {
    (0..length)
        .map(|position| {
            if mask & (1 << position) == 0 {
                b'A'
            } else {
                b'N'
            }
        })
        .collect()
}

fn build_index(sequence: &[u8]) -> ReferenceIndex {
    ReferenceIndex::build(
        vec![ContigInput::new(b"origin".to_vec(), normalized(sequence))],
        ReferenceBuildLimits::MAX,
    )
    .expect("bounded reference builds")
}

fn ot_certified(query: &[u8], budget: u64) -> bsbit_align::search::seed::CertifiedSeedPlan {
    let strands = AdmissibleStrands::new(&[BisulfiteStrand::OT]).expect("one strand");
    match schedule_proof_seeds(
        normalized(query),
        strands,
        EditDistance::new(budget),
        ProofSeedLimits::MAX,
    )
    .expect("bounded schedule succeeds")
    {
        ProofSeedOutcome::Certified(plan) => plan,
        other => panic!("expected certified schedule, observed {other:?}"),
    }
}

fn signed_diagonal(anchor: &CandidateAnchor) -> i128 {
    match anchor.diagonal().shift() {
        CoordinateShift::Zero => 0,
        CoordinateShift::Forward(value) => i128::from(value.get()),
        CoordinateShift::Backward(value) => -i128::from(value.get()),
    }
}

#[test]
fn implementation_equals_independent_raw_byte_oracle_exhaustively() {
    let mut comparisons = 0_u64;
    for length in 0_usize..=7 {
        for n_mask in 0_usize..(1_usize << length) {
            let query = query_from_n_mask(length, n_mask);
            for budget in 0..=length + 1 {
                for strand_mask in 1_u8..16 {
                    let ranks = rank_subset(strand_mask);
                    assert_eq!(
                        implementation_snapshot(&query, &ranks, budget),
                        oracle_snapshot(&query, &ranks, budget),
                        "length={length} n_mask={n_mask:#x} budget={budget} strands={strand_mask:#x}"
                    );
                    comparisons += 1;
                }
            }
        }
    }

    for length in [9_usize, 10, 17, 18, 29, 30] {
        for budget in [0, 1, length - 1, length, length + 1] {
            let mut query = vec![b'A'; length];
            if length > 1 {
                query[length / 2] = b'N';
            }
            for strand_mask in 1_u8..16 {
                let ranks = rank_subset(strand_mask);
                assert_eq!(
                    implementation_snapshot(&query, &ranks, budget),
                    oracle_snapshot(&query, &ranks, budget)
                );
                comparisons += 1;
            }
        }
    }
    assert_eq!(comparisons, 31_170);
}

#[test]
fn every_bounded_independent_edit_footprint_leaves_an_emitted_block() {
    let mut footprints = 0_usize;
    let mut displacement_assignments = 0_usize;
    for length in 1_usize..=10 {
        for budget in 0..=3.min(length - 1) {
            let blocks = budget + 1;
            let n_patterns = [
                Vec::new(),
                vec![0],
                vec![length / 2],
                vec![length - 1],
                vec![0, length - 1],
            ];
            for n_positions in n_patterns {
                let mut query = vec![b'A'; length];
                for position in n_positions.into_iter().take(budget) {
                    query[position] = b'N';
                }
                let observed = exhaustive_edit_footprint_stats(&query, budget);
                assert!(
                    observed.footprints > 0,
                    "length={length} budget={budget} blocks={blocks}"
                );
                assert!(observed.maximum_absolute_displacement <= budget);
                footprints += observed.footprints;
                displacement_assignments += observed.displacement_assignments;
            }
        }
    }
    assert_eq!(footprints, 6_491);
    assert_eq!(displacement_assignments, 133_115);
}

#[test]
fn every_strand_subset_is_invariant_under_every_supplied_order() {
    let query = b"ACGTNACGTTGCA";
    for strand_mask in 1_u8..16 {
        let ranks = rank_subset(strand_mask);
        let expected = implementation_snapshot(query, &ranks, 2);
        for permutation in permutations(&ranks) {
            assert_eq!(implementation_snapshot(query, &permutation, 2), expected);
        }
    }
}

#[test]
fn discrepancy_0010_exact17_composes_to_the_true_origin_and_caps_fail_whole() {
    let sequence = b"ACGTTGCAACGATTCGA";
    let index = build_index(sequence);
    let plan = ot_certified(sequence, 0);
    let candidates = candidates_for_fixed_seeds(
        &index,
        plan.fixed_plan(),
        ReferenceQueryLimits::MAX,
        UNBOUNDED_CANDIDATES,
    )
    .expect("exact17 candidate generation succeeds");
    assert!(
        candidates
            .anchors()
            .iter()
            .any(|anchor| anchor.strand() == BisulfiteStrand::OT && signed_diagonal(anchor) == 0)
    );

    assert!(matches!(
        candidates_for_fixed_seeds(
            &index,
            plan.fixed_plan(),
            ReferenceQueryLimits::new(u64::MAX, 0),
            UNBOUNDED_CANDIDATES,
        ),
        Err(CandidateError::Search {
            source: ReferenceQueryError::HitLimitExceeded {
                requested: 1,
                maximum: 0,
            },
            ..
        })
    ));
}

#[test]
fn insertion_and_deletion_origins_have_an_anchor_within_the_certified_budget() {
    let reference = b"ACGTTGCAACGATTCGAGCTTAGC";
    let index = build_index(reference);

    let mut insertion = vec![b'T'];
    insertion.extend_from_slice(reference);
    let insertion_plan = ot_certified(&insertion, 1);
    let insertion_candidates = candidates_for_fixed_seeds(
        &index,
        insertion_plan.fixed_plan(),
        ReferenceQueryLimits::MAX,
        UNBOUNDED_CANDIDATES,
    )
    .expect("insertion-origin candidates succeed");
    assert!(
        insertion_candidates
            .anchors()
            .iter()
            .any(|anchor| anchor.strand() == BisulfiteStrand::OT && signed_diagonal(anchor) == -1)
    );

    let mut deletion = reference[..5].to_vec();
    deletion.extend_from_slice(&reference[6..]);
    let deletion_plan = ot_certified(&deletion, 1);
    let deletion_candidates = candidates_for_fixed_seeds(
        &index,
        deletion_plan.fixed_plan(),
        ReferenceQueryLimits::MAX,
        UNBOUNDED_CANDIDATES,
    )
    .expect("deletion-origin candidates succeed");
    assert!(
        deletion_candidates
            .anchors()
            .iter()
            .any(|anchor| anchor.strand() == BisulfiteStrand::OT && signed_diagonal(anchor) == 1)
    );

    assert_eq!(
        insertion_plan.certificate().maximum_diagonal_displacement(),
        1
    );
    assert_eq!(
        deletion_plan.certificate().maximum_diagonal_displacement(),
        1
    );
}

#[test]
fn repeated_and_eight_worker_schedules_are_semantically_identical() {
    let query = b"ACGTNACGTTGCAACGATTCGAGCTTAGC";
    let ranks = [0_u8, 1, 2, 3];
    let expected = implementation_snapshot(query, &ranks, 3);
    assert_eq!(implementation_snapshot(query, &ranks, 3), expected);

    thread::scope(|scope| {
        let handles = (0..8)
            .map(|worker| {
                scope.spawn(move || {
                    let mut supplied = ranks;
                    let shift = worker % supplied.len();
                    supplied.rotate_left(shift);
                    implementation_snapshot(query, &supplied, 3)
                })
            })
            .collect::<Vec<_>>();
        for handle in handles {
            assert_eq!(handle.join().expect("worker did not panic"), expected);
        }
    });
}

#[test]
fn public_plan_limits_reject_complete_results_without_truncation() {
    let strands = AdmissibleStrands::new(&BisulfiteStrand::ALL).expect("all strands");
    let query = normalized(b"ACGTACGTAA");
    let error = schedule_proof_seeds(
        query,
        strands,
        EditDistance::new(2),
        ProofSeedLimits::new(3, SeedPlanLimits::new(12, 39)),
    )
    .expect_err("complete forty-base plan exceeds limit");
    assert_eq!(
        error,
        bsbit_align::search::seed::ProofSeedError::TotalSeedBasesLimitExceeded {
            requested: 40,
            maximum: 39,
        }
    );
}
