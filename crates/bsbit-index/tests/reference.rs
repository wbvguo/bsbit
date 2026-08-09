//! Independent scientific and API tests for the Level 2B projected reference.

#[path = "support/reference_oracle.rs"]
mod reference_oracle;

use reference_oracle::{
    CANONICAL, OracleContig, OracleHit, OracleStrand, canonical_runs, direct_search,
    enumerate_strings, origin_cases, project, reverse_complement,
};

use std::io::Write;
use std::mem::size_of;
use std::process::{Command, Stdio};
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};
use std::thread;

use bsbit_core::bisulfite::BisulfiteStrand;
use bsbit_core::coordinate::{CoordinateDomain, CoordinateError};
use bsbit_core::sequence::{NormalizedSequence, normalize_dna};
use bsbit_index::reference::{
    ContigId, ContigInput, ProjectedMatches, ReferenceAccessError, ReferenceAllocation,
    ReferenceArithmetic, ReferenceBuildError, ReferenceBuildLimits, ReferenceIndex,
    ReferenceInstanceId, ReferenceLocateError, ReferenceLocateInvariant, ReferenceQueryCounter,
    ReferenceQueryError, ReferenceQueryLimits, ReferenceResource,
};
use bsbit_index::storage::fm::{FmAllocation, FmError, FmIndex};

const NAMED_CATALOG: [OracleContig<'static>; 3] = [
    OracleContig {
        name: b"alpha",
        sequence: b"ACGTCNTA",
    },
    OracleContig {
        name: b"beta",
        sequence: b"GCA",
    },
    OracleContig {
        name: b"unknown",
        sequence: b"NN",
    },
];

fn oracle_hit(contig_ordinal: u64, start: u64, end: u64, strand: OracleStrand) -> OracleHit {
    OracleHit {
        contig_ordinal,
        start,
        end,
        strand,
    }
}

fn normalized(raw: &[u8]) -> NormalizedSequence {
    normalize_dna(raw).expect("test sequence is normalized A/C/G/T/N")
}

fn build_catalog(catalog: &[(&[u8], &[u8])]) -> ReferenceIndex {
    let inputs = catalog
        .iter()
        .map(|(name, sequence)| ContigInput::new(name.to_vec(), normalized(sequence)))
        .collect();
    ReferenceIndex::build(inputs, ReferenceBuildLimits::MAX)
        .expect("small test reference should build")
}

fn expected_strand(strand: OracleStrand) -> BisulfiteStrand {
    match strand {
        OracleStrand::Ot => BisulfiteStrand::OT,
        OracleStrand::Ob => BisulfiteStrand::OB,
        OracleStrand::Ctot => BisulfiteStrand::CTOT,
        OracleStrand::Ctob => BisulfiteStrand::CTOB,
    }
}

fn snapshot(index: &ReferenceIndex, strand: OracleStrand, raw_pattern: &[u8]) -> Vec<OracleHit> {
    let matches = index
        .exact_search(
            expected_strand(strand),
            &normalized(raw_pattern),
            ReferenceQueryLimits::MAX,
        )
        .expect("small canonical query should search");
    assert_eq!(matches.strand(), expected_strand(strand));
    assert_eq!(
        matches.pattern_len(),
        u64::try_from(raw_pattern.len()).expect("test pattern length fits u64")
    );
    assert!(matches.matched_interval_count() <= matches.exact_hit_count());
    let hits = index
        .locate(&matches)
        .expect("matches belong to the receiving reference");
    assert_eq!(
        matches.exact_hit_count(),
        u64::try_from(hits.len()).expect("test hit count fits u64")
    );
    hits.iter()
        .map(|hit| OracleHit {
            contig_ordinal: hit.contig().ordinal(),
            start: hit.interval().start(),
            end: hit.interval().end(),
            strand,
        })
        .collect()
}

fn named_build() -> ReferenceIndex {
    build_catalog(&[
        (b"alpha", b"ACGTCNTA"),
        (b"beta", b"GCA"),
        (b"unknown", b"NN"),
    ])
}

#[test]
fn borrowed_pattern_slices_preserve_complete_search_semantics_and_offsets() {
    let index = named_build();
    let owned = normalized(b"GT");
    let via_deref = index
        .exact_search(BisulfiteStrand::OT, &owned, ReferenceQueryLimits::MAX)
        .expect("owned sequence deref-coerces to the borrowed API");
    let via_direct_slice = index
        .exact_search(
            BisulfiteStrand::OT,
            owned.bases(),
            ReferenceQueryLimits::MAX,
        )
        .expect("direct normalized slice is accepted");
    let padded = normalized(b"NNGTNN");
    let via_nonzero_subslice = index
        .exact_search(
            BisulfiteStrand::OT,
            &padded.bases()[2..4],
            ReferenceQueryLimits::MAX,
        )
        .expect("nonzero normalized subslice is accepted without copying");

    let semantic_hits = |matches: &ProjectedMatches| {
        index
            .locate(matches)
            .expect("all compared artifacts belong to this index")
            .into_iter()
            .map(|hit| {
                (
                    hit.contig().ordinal(),
                    hit.interval().start(),
                    hit.interval().end(),
                    hit.strand(),
                )
            })
            .collect::<Vec<_>>()
    };
    assert_eq!(via_deref.pattern_len(), 2);
    assert_eq!(
        via_deref.exact_hit_count(),
        via_direct_slice.exact_hit_count()
    );
    assert_eq!(
        via_deref.exact_hit_count(),
        via_nonzero_subslice.exact_hit_count()
    );
    assert_eq!(semantic_hits(&via_deref), semantic_hits(&via_direct_slice));
    assert_eq!(
        semantic_hits(&via_deref),
        semantic_hits(&via_nonzero_subslice)
    );

    assert_eq!(
        index
            .exact_search(
                BisulfiteStrand::OT,
                &owned.bases()[0..0],
                ReferenceQueryLimits::MAX,
            )
            .unwrap_err(),
        ReferenceQueryError::EmptyPattern
    );
    let n_with_nonzero_origin = normalized(b"CCNGTN");
    assert_eq!(
        index
            .exact_search(
                BisulfiteStrand::OT,
                &n_with_nonzero_origin.bases()[1..5],
                ReferenceQueryLimits::MAX,
            )
            .unwrap_err(),
        ReferenceQueryError::UnsearchableBase { offset: 1 }
    );

    let limits = ReferenceQueryLimits::new(11, 13).with_max_exact_hits(17);
    assert_eq!(limits.max_exact_hits(), 17);
}

#[test]
fn independent_oracle_reproduces_accepted_scientific_evidence() {
    assert_eq!(canonical_runs(b"ACNNTGN"), vec![0..2, 4..6]);

    let mut duality_cases = 0;
    for sequence in enumerate_strings(&CANONICAL, 5) {
        assert_eq!(
            project(&reverse_complement(&sequence), b'C', b'T'),
            reverse_complement(&project(&sequence, b'G', b'A'))
        );
        assert_eq!(
            project(&reverse_complement(&sequence), b'G', b'A'),
            reverse_complement(&project(&sequence, b'C', b'T'))
        );
        duality_cases += 2;
    }
    assert_eq!(duality_cases, 2_730);

    let named = [
        (
            OracleStrand::Ot,
            &b"GT"[..],
            vec![
                oracle_hit(0, 2, 4, OracleStrand::Ot),
                oracle_hit(1, 0, 2, OracleStrand::Ot),
            ],
        ),
        (
            OracleStrand::Ctob,
            &b"ACA"[..],
            vec![
                oracle_hit(0, 0, 3, OracleStrand::Ctob),
                oracle_hit(1, 0, 3, OracleStrand::Ctob),
            ],
        ),
        (
            OracleStrand::Ob,
            &b"GAT"[..],
            vec![oracle_hit(0, 2, 5, OracleStrand::Ob)],
        ),
        (
            OracleStrand::Ob,
            &b"TGT"[..],
            vec![
                oracle_hit(0, 0, 3, OracleStrand::Ob),
                oracle_hit(1, 0, 3, OracleStrand::Ob),
            ],
        ),
        (
            OracleStrand::Ctot,
            &b"TAC"[..],
            vec![oracle_hit(1, 0, 3, OracleStrand::Ctot)],
        ),
        (OracleStrand::Ot, &b"TTT"[..], vec![]),
        (OracleStrand::Ot, &b"TAG"[..], vec![]),
    ];
    for (strand, pattern, expected) in named {
        assert_eq!(
            direct_search(&NAMED_CATALOG, strand, pattern),
            expected,
            "{} {}",
            strand.label(),
            String::from_utf8_lossy(pattern)
        );
    }

    let cases = origin_cases(&NAMED_CATALOG);
    for strand in OracleStrand::ALL {
        let strand_cases = cases
            .iter()
            .filter(|case| case.expected_origin.strand == strand)
            .collect::<Vec<_>>();
        let expected_count = match strand {
            OracleStrand::Ot | OracleStrand::Ctot => 43,
            OracleStrand::Ob | OracleStrand::Ctob => 36,
        };
        assert_eq!(strand_cases.len(), expected_count);
        for case in strand_cases {
            assert!(
                direct_search(&NAMED_CATALOG, strand, &case.raw_pattern)
                    .contains(&case.expected_origin),
                "origin recall failed for {} raw {} at {:?}",
                strand.label(),
                String::from_utf8_lossy(&case.raw_pattern),
                case.expected_origin
            );
        }
    }
    assert_eq!(cases.len(), 158);
}

fn build_owned_catalog(catalog: &[(Vec<u8>, Vec<u8>)]) -> ReferenceIndex {
    let inputs = catalog
        .iter()
        .map(|(name, sequence)| ContigInput::new(name.clone(), normalized(sequence)))
        .collect();
    ReferenceIndex::build(inputs, ReferenceBuildLimits::MAX)
        .expect("small owned test reference should build")
}

fn oracle_views(catalog: &[(Vec<u8>, Vec<u8>)]) -> Vec<OracleContig<'_>> {
    catalog
        .iter()
        .map(|(name, sequence)| OracleContig { name, sequence })
        .collect()
}

fn assert_run_metrics_equal_oracle(index: &ReferenceIndex, catalog: &[(Vec<u8>, Vec<u8>)]) {
    let expected_original_bases = catalog
        .iter()
        .map(|(_, sequence)| u64::try_from(sequence.len()).expect("test length fits u64"))
        .sum::<u64>();
    let expected_canonical_bases = catalog
        .iter()
        .flat_map(|(_, sequence)| sequence.iter())
        .filter(|base| **base != b'N')
        .count();
    let expected_canonical_bases =
        u64::try_from(expected_canonical_bases).expect("test canonical count fits u64");
    let expected_runs = catalog
        .iter()
        .map(|(_, sequence)| canonical_runs(sequence).len())
        .sum::<usize>();
    let expected_runs = u64::try_from(expected_runs).expect("test run count fits u64");

    let metrics = index.metrics();
    assert_eq!(
        metrics.contig_count(),
        u64::try_from(catalog.len()).expect("test contig count fits u64")
    );
    assert_eq!(metrics.total_reference_bases(), expected_original_bases);
    assert_eq!(metrics.canonical_bases(), expected_canonical_bases);
    assert_eq!(metrics.canonical_run_count(), expected_runs);
    assert_eq!(metrics.lane_count(), 4 * expected_runs);
    assert_eq!(metrics.projected_bases(), 4 * expected_canonical_bases);
    assert_eq!(
        metrics.projected_suffix_rows(),
        4 * (expected_canonical_bases + expected_runs)
    );
}

#[test]
fn named_hits_equal_the_independent_oracle() {
    let index = named_build();
    assert_eq!(index.contig_count(), 3);

    for (ordinal, name, sequence) in [
        (0, &b"alpha"[..], &b"ACGTCNTA"[..]),
        (1, &b"beta"[..], &b"GCA"[..]),
        (2, &b"unknown"[..], &b"NN"[..]),
    ] {
        let id = index.contig_id(ordinal).expect("named ordinal exists");
        let view = index
            .resolve_contig(&id)
            .expect("ID belongs to the named reference");
        assert_eq!(view.ordinal(), ordinal);
        assert_eq!(view.name(), name);
        assert_eq!(view.sequence().to_ascii(), sequence);
    }

    let cases = [
        (OracleStrand::Ot, &b"GT"[..]),
        (OracleStrand::Ctob, &b"ACA"[..]),
        (OracleStrand::Ob, &b"GAT"[..]),
        (OracleStrand::Ob, &b"TGT"[..]),
        (OracleStrand::Ctot, &b"TAC"[..]),
        (OracleStrand::Ot, &b"TTT"[..]),
        (OracleStrand::Ot, &b"TAG"[..]),
    ];
    for (strand, pattern) in cases {
        assert_eq!(
            snapshot(&index, strand, pattern),
            direct_search(&NAMED_CATALOG, strand, pattern),
            "{} {}",
            strand.label(),
            String::from_utf8_lossy(pattern)
        );
    }
}

#[test]
fn every_origin_aware_case_is_recalled_by_its_intended_lane() {
    let index = named_build();
    let cases = origin_cases(&NAMED_CATALOG);
    let expected_counts = [
        (OracleStrand::Ot, 43),
        (OracleStrand::Ob, 36),
        (OracleStrand::Ctot, 43),
        (OracleStrand::Ctob, 36),
    ];

    for (strand, expected_count) in expected_counts {
        let strand_cases = cases
            .iter()
            .filter(|case| case.expected_origin.strand == strand)
            .collect::<Vec<_>>();
        assert_eq!(strand_cases.len(), expected_count);
        for case in strand_cases {
            let observed = snapshot(&index, strand, &case.raw_pattern);
            assert!(
                observed.contains(&case.expected_origin),
                "projected origin recall failed for {} raw {} at {:?}",
                strand.label(),
                String::from_utf8_lossy(&case.raw_pattern),
                case.expected_origin
            );
        }
    }
}

#[test]
fn exhaustive_bounded_catalogs_equal_direct_raw_byte_scanning() {
    let reference_strings = enumerate_strings(&reference_oracle::REFERENCE_ALPHABET, 3);
    let one_contig = reference_strings
        .iter()
        .filter(|sequence| !sequence.is_empty())
        .collect::<Vec<_>>();
    let patterns = enumerate_strings(&CANONICAL, 2)
        .into_iter()
        .filter(|pattern| !pattern.is_empty())
        .collect::<Vec<_>>();

    assert_eq!(one_contig.len(), 155);
    assert_eq!(patterns.len(), 20);

    let mut searches = 0_u64;
    for sequence in one_contig {
        let catalog = vec![(b"x".to_vec(), sequence.clone())];
        let index = build_owned_catalog(&catalog);
        let oracle = oracle_views(&catalog);
        assert_run_metrics_equal_oracle(&index, &catalog);
        for strand in OracleStrand::ALL {
            for pattern in &patterns {
                assert_eq!(
                    snapshot(&index, strand, pattern),
                    direct_search(&oracle, strand, pattern),
                    "one-contig sequence={} strand={} pattern={}",
                    String::from_utf8_lossy(sequence),
                    strand.label(),
                    String::from_utf8_lossy(pattern)
                );
                searches += 1;
            }
        }
    }

    let singleton_references = reference_oracle::REFERENCE_ALPHABET;
    for left in singleton_references {
        for right in singleton_references {
            let catalog = vec![(b"x".to_vec(), vec![left]), (b"y".to_vec(), vec![right])];
            let index = build_owned_catalog(&catalog);
            let oracle = oracle_views(&catalog);
            assert_run_metrics_equal_oracle(&index, &catalog);
            for strand in OracleStrand::ALL {
                for pattern in &patterns {
                    assert_eq!(
                        snapshot(&index, strand, pattern),
                        direct_search(&oracle, strand, pattern),
                        "two-contig sequence={}/{} strand={} pattern={}",
                        char::from(left),
                        char::from(right),
                        strand.label(),
                        String::from_utf8_lossy(pattern)
                    );
                    searches += 1;
                }
            }
        }
    }

    assert_eq!(searches, 14_400);
}

#[test]
fn barriers_duplicate_content_and_original_bases_are_preserved() {
    let duplicates = build_catalog(&[(b"x", b"AC"), (b"y", b"AC")]);
    assert_eq!(
        snapshot(&duplicates, OracleStrand::Ot, b"AT"),
        vec![
            oracle_hit(0, 0, 2, OracleStrand::Ot),
            oracle_hit(1, 0, 2, OracleStrand::Ot),
        ]
    );
    assert!(snapshot(&duplicates, OracleStrand::Ot, b"ATAT").is_empty());

    let collapsed = build_catalog(&[(b"collapse", b"TC")]);
    let collapsed_hits = snapshot(&collapsed, OracleStrand::Ot, b"T");
    assert_eq!(
        collapsed_hits,
        vec![
            oracle_hit(0, 0, 1, OracleStrand::Ot),
            oracle_hit(0, 1, 2, OracleStrand::Ot),
        ]
    );
    let collapse_id = collapsed.contig_id(0).expect("contig exists");
    assert_eq!(
        collapsed
            .resolve_contig(&collapse_id)
            .expect("ID is local")
            .sequence()
            .to_ascii(),
        b"TC"
    );

    let ambiguity = [b"ACGT".as_slice(), &[b'N'; 64], b"TGCA".as_slice()].concat();
    let n_index = build_catalog(&[(b"ambiguous", &ambiguity)]);
    let n_id = n_index.contig_id(0).expect("ambiguous contig exists");
    let retained = n_index
        .resolve_contig(&n_id)
        .expect("ID belongs to ambiguity index")
        .sequence()
        .to_ascii();
    assert_eq!(retained, ambiguity);
    assert_eq!(&retained[4..68], &[b'N'; 64]);
    assert!(snapshot(&n_index, OracleStrand::Ot, b"ACGTTGCA").is_empty());

    let mut appended_catalog = vec![
        (b"alpha".to_vec(), b"ACGTCNTA".to_vec()),
        (b"beta".to_vec(), b"GCA".to_vec()),
        (b"unknown".to_vec(), b"NN".to_vec()),
    ];
    let baseline = build_owned_catalog(&appended_catalog);
    appended_catalog.push((b"unrelated".to_vec(), b"TTAC".to_vec()));
    let appended = build_owned_catalog(&appended_catalog);
    for strand in OracleStrand::ALL {
        for pattern in [b"A".as_slice(), b"GT", b"TAC"] {
            let existing = snapshot(&appended, strand, pattern)
                .into_iter()
                .filter(|hit| hit.contig_ordinal < 3)
                .collect::<Vec<_>>();
            assert_eq!(existing, snapshot(&baseline, strand, pattern));
        }
    }
}

#[test]
fn discrepancy_0004_cross_contig_read_is_not_an_exact_ot_hit() {
    let reference = include_str!("../../../tests/fixtures/contig-boundary-exact/reference.fa");
    let mut catalog = Vec::<(Vec<u8>, Vec<u8>)>::new();
    for record in reference.split('>').filter(|record| !record.is_empty()) {
        let mut lines = record.lines();
        let name = lines
            .next()
            .expect("fixture record has a name")
            .as_bytes()
            .to_vec();
        let sequence = lines.flat_map(str::bytes).collect();
        catalog.push((name, sequence));
    }
    assert_eq!(catalog.len(), 2);

    let fastq = include_str!("../../../tests/fixtures/contig-boundary-exact/cross-boundary.fastq");
    let pattern = fastq.lines().nth(1).expect("fixture FASTQ has a sequence");
    assert_eq!(pattern.len(), 75);

    let index = build_owned_catalog(&catalog);
    assert!(snapshot(&index, OracleStrand::Ot, pattern.as_bytes()).is_empty());
    assert!(
        direct_search(
            &oracle_views(&catalog),
            OracleStrand::Ot,
            pattern.as_bytes()
        )
        .is_empty()
    );
}

#[test]
fn owner_identity_is_exact_and_checked_before_empty_artifacts() {
    let index = named_build();
    let clone = index.clone();
    let matches = index
        .exact_search(
            BisulfiteStrand::OT,
            &normalized(b"GT"),
            ReferenceQueryLimits::MAX,
        )
        .expect("local query succeeds");
    assert!(clone.locate(&matches).is_ok());

    let initial_instance = index.instance_id();
    let clone_instance = clone.instance_id();
    assert!(initial_instance.is_same_instance(&clone_instance));
    assert!(!format!("{initial_instance:?}").contains("0x"));
    let clone_matches = clone
        .exact_search(
            BisulfiteStrand::OT,
            &normalized(b"GT"),
            ReferenceQueryLimits::MAX,
        )
        .expect("clone query succeeds");
    assert!(matches.is_same_instance(&clone_matches));
    assert!(matches.belongs_to_instance(&initial_instance));
    assert!(clone_matches.belongs_to_instance(&clone_instance));

    let initial_contig = index.contig_id(0).expect("contig exists");
    let clone_contig = clone.contig_id(0).expect("clone has same contig");
    assert!(initial_contig.is_same_instance(&clone_contig));
    assert!(initial_contig.is_same_contig(&clone_contig));

    let rebuild = named_build();
    let rebuilt_instance = rebuild.instance_id();
    let rebuilt_contig = rebuild.contig_id(0).expect("rebuilt contig exists");
    let rebuilt_matches = rebuild
        .exact_search(
            BisulfiteStrand::OT,
            &normalized(b"GT"),
            ReferenceQueryLimits::MAX,
        )
        .expect("rebuilt query succeeds");
    assert!(!matches.is_same_instance(&rebuilt_matches));
    assert!(!matches.belongs_to_instance(&rebuilt_instance));
    assert!(!initial_instance.is_same_instance(&rebuilt_instance));
    assert!(!initial_contig.is_same_instance(&rebuilt_contig));
    assert!(!initial_contig.is_same_contig(&rebuilt_contig));
    for strand in OracleStrand::ALL {
        for pattern in [b"A".as_slice(), b"GT", b"GAT", b"TAC", b"ACA"] {
            assert_eq!(
                snapshot(&rebuild, strand, pattern),
                snapshot(&index, strand, pattern),
                "independent rebuild differs for {} {}",
                strand.label(),
                String::from_utf8_lossy(pattern)
            );
        }
    }
    assert!(rebuild.locate(&matches).is_err());
    assert!(rebuild.resolve_contig(&initial_contig).is_err());

    let different = build_catalog(&[
        (b"alpha", b"TTTTTNTA"),
        (b"beta", b"AAA"),
        (b"unknown", b"NN"),
    ]);
    assert!(different.locate(&matches).is_err());

    let empty_matches = index
        .exact_search(
            BisulfiteStrand::OT,
            &normalized(b"CCCCCCCC"),
            ReferenceQueryLimits::MAX,
        )
        .expect("zero-hit query still returns an owner-bound artifact");
    assert!(empty_matches.is_empty());
    assert!(rebuild.locate(&empty_matches).is_err());
    assert!(index.contig_id(index.contig_count()).is_err());
}

#[test]
fn shared_immutable_queries_are_deterministic_across_eight_workers() {
    fn assert_send_sync<T: Send + Sync>() {}
    assert_send_sync::<ReferenceIndex>();

    let index = Arc::new(named_build());
    let workload = [
        (OracleStrand::Ot, &b"GT"[..]),
        (OracleStrand::Ob, &b"GAT"[..]),
        (OracleStrand::Ctot, &b"TAC"[..]),
        (OracleStrand::Ctob, &b"ACA"[..]),
    ];
    let serial = workload
        .iter()
        .map(|(strand, pattern)| snapshot(&index, *strand, pattern))
        .collect::<Vec<_>>();

    thread::scope(|scope| {
        let mut workers = Vec::new();
        for worker_index in 0..8 {
            let shared = Arc::clone(&index);
            workers.push(scope.spawn(move || {
                let mut observed = vec![Vec::new(); workload.len()];
                for step in 0..workload.len() {
                    let forward = (step + worker_index) % workload.len();
                    let query_index = if worker_index % 2 == 0 {
                        forward
                    } else {
                        workload.len() - 1 - forward
                    };
                    let (strand, pattern) = workload[query_index];
                    observed[query_index] = snapshot(&shared, strand, pattern);
                }
                observed
            }));
        }
        for worker in workers {
            assert_eq!(worker.join().expect("worker must not panic"), serial);
        }
    });
}

#[test]
fn metrics_and_every_configured_resource_limit_are_exact() {
    let index = named_build();
    let metrics = index.metrics();
    assert_eq!(metrics.contig_count(), 3);
    assert_eq!(metrics.total_name_bytes(), 16);
    assert_eq!(metrics.total_reference_bases(), 13);
    assert_eq!(metrics.canonical_bases(), 10);
    assert_eq!(metrics.canonical_run_count(), 3);
    assert_eq!(metrics.lane_count(), 12);
    assert_eq!(metrics.projected_bases(), 40);
    assert_eq!(metrics.projected_suffix_rows(), 52);

    let expected_retained = 12_u64 * u64::try_from(size_of::<FmIndex>()).unwrap()
        + 4 * (10 + 3)
            * (u64::try_from(size_of::<usize>()).unwrap()
                + u64::try_from(size_of::<u8>()).unwrap())
        + 4 * (10 + 2 * 3) * u64::try_from(size_of::<[u64; 4]>()).unwrap();
    assert_eq!(metrics.estimated_retained_fm_bytes(), expected_retained);

    let exact = ReferenceBuildLimits::MAX
        .with_max_contigs(3)
        .with_max_total_name_bytes(16)
        .with_max_total_reference_bases(13)
        .with_max_canonical_runs(3)
        .with_max_suffix_rows_per_lane(6)
        .with_max_lanes(12)
        .with_max_projected_bases(40)
        .with_max_projected_suffix_rows(52)
        .with_max_estimated_retained_fm_bytes(expected_retained);
    let above = ReferenceBuildLimits::MAX
        .with_max_contigs(4)
        .with_max_total_name_bytes(17)
        .with_max_total_reference_bases(14)
        .with_max_canonical_runs(4)
        .with_max_suffix_rows_per_lane(7)
        .with_max_lanes(13)
        .with_max_projected_bases(41)
        .with_max_projected_suffix_rows(53)
        .with_max_estimated_retained_fm_bytes(expected_retained + 1);
    for limits in [exact, above] {
        let inputs = vec![
            ContigInput::new(b"alpha".to_vec(), normalized(b"ACGTCNTA")),
            ContigInput::new(b"beta".to_vec(), normalized(b"GCA")),
            ContigInput::new(b"unknown".to_vec(), normalized(b"NN")),
        ];
        ReferenceIndex::build(inputs, limits).expect("at/above every limit must build");
    }

    let below = [
        ReferenceBuildLimits::MAX.with_max_contigs(2),
        ReferenceBuildLimits::MAX.with_max_total_name_bytes(15),
        ReferenceBuildLimits::MAX.with_max_total_reference_bases(12),
        ReferenceBuildLimits::MAX.with_max_canonical_runs(2),
        ReferenceBuildLimits::MAX.with_max_suffix_rows_per_lane(5),
        ReferenceBuildLimits::MAX.with_max_lanes(11),
        ReferenceBuildLimits::MAX.with_max_projected_bases(39),
        ReferenceBuildLimits::MAX.with_max_projected_suffix_rows(51),
        ReferenceBuildLimits::MAX.with_max_estimated_retained_fm_bytes(expected_retained - 1),
    ];
    for limits in below {
        let inputs = vec![
            ContigInput::new(b"alpha".to_vec(), normalized(b"ACGTCNTA")),
            ContigInput::new(b"beta".to_vec(), normalized(b"GCA")),
            ContigInput::new(b"unknown".to_vec(), normalized(b"NN")),
        ];
        assert!(ReferenceIndex::build(inputs, limits).is_err());
    }
}

#[test]
fn validation_and_query_limits_fail_whole_operations() {
    assert!(ReferenceIndex::build(Vec::new(), ReferenceBuildLimits::MAX).is_err());
    assert!(
        ReferenceIndex::build(
            vec![ContigInput::new(Vec::new(), normalized(b"A"))],
            ReferenceBuildLimits::MAX,
        )
        .is_err()
    );
    assert!(
        ReferenceIndex::build(
            vec![ContigInput::new(b"empty".to_vec(), normalized(b""))],
            ReferenceBuildLimits::MAX,
        )
        .is_err()
    );
    assert!(
        ReferenceIndex::build(
            vec![
                ContigInput::new(b"same".to_vec(), normalized(b"A")),
                ContigInput::new(b"same".to_vec(), normalized(b"T")),
            ],
            ReferenceBuildLimits::MAX,
        )
        .is_err()
    );

    let index = named_build();
    assert!(
        index
            .exact_search(
                BisulfiteStrand::OT,
                &normalized(b""),
                ReferenceQueryLimits::MAX,
            )
            .is_err()
    );
    assert!(
        index
            .exact_search(
                BisulfiteStrand::OT,
                &normalized(b"AN"),
                ReferenceQueryLimits::new(1, u64::MAX),
            )
            .is_err()
    );
    assert!(
        index
            .exact_search(
                BisulfiteStrand::OT,
                &normalized(b"AN"),
                ReferenceQueryLimits::new(2, u64::MAX),
            )
            .is_err()
    );

    let exact = index
        .exact_search(
            BisulfiteStrand::OT,
            &normalized(b"GT"),
            ReferenceQueryLimits::new(2, 2),
        )
        .expect("at-limit query returns its complete result");
    assert_eq!(exact.exact_hit_count(), 2);
    assert_eq!(index.locate(&exact).expect("local matches locate").len(), 2);
    index
        .exact_search(
            BisulfiteStrand::OT,
            &normalized(b"GT"),
            ReferenceQueryLimits::new(3, 3),
        )
        .expect("above-limit query succeeds");
    assert!(
        index
            .exact_search(
                BisulfiteStrand::OT,
                &normalized(b"GT"),
                ReferenceQueryLimits::new(1, u64::MAX),
            )
            .is_err()
    );
    assert!(
        index
            .exact_search(
                BisulfiteStrand::OT,
                &normalized(b"GT"),
                ReferenceQueryLimits::new(2, 1),
            )
            .is_err()
    );

    let zero = index
        .exact_search(
            BisulfiteStrand::OT,
            &normalized(b"CCCCCCCC"),
            ReferenceQueryLimits::new(8, 0),
        )
        .expect("zero-hit query is admitted by a zero hit limit");
    assert!(zero.is_empty());
}

fn build_named_with_limits(
    limits: ReferenceBuildLimits,
) -> Result<ReferenceIndex, ReferenceBuildError> {
    ReferenceIndex::build(
        vec![
            ContigInput::new(b"alpha".to_vec(), normalized(b"ACGTCNTA")),
            ContigInput::new(b"beta".to_vec(), normalized(b"GCA")),
            ContigInput::new(b"unknown".to_vec(), normalized(b"NN")),
        ],
        limits,
    )
}

#[allow(
    clippy::too_many_lines,
    reason = "the public total error order is intentionally asserted in one readable sequence"
)]
#[test]
fn reachable_error_fields_and_priority_are_exact() {
    assert_eq!(
        ReferenceIndex::build(Vec::new(), ReferenceBuildLimits::MAX).unwrap_err(),
        ReferenceBuildError::EmptyReference
    );
    assert_eq!(
        ReferenceIndex::build(
            vec![ContigInput::new(Vec::new(), normalized(b""))],
            ReferenceBuildLimits::MAX.with_max_contigs(0),
        )
        .unwrap_err(),
        ReferenceBuildError::LimitExceeded {
            resource: ReferenceResource::Contigs,
            requested: 1,
            maximum: 0,
        }
    );
    assert_eq!(
        ReferenceIndex::build(
            vec![ContigInput::new(b"x".to_vec(), normalized(b""))],
            ReferenceBuildLimits::MAX.with_max_total_name_bytes(0),
        )
        .unwrap_err(),
        ReferenceBuildError::LimitExceeded {
            resource: ReferenceResource::TotalNameBytes,
            requested: 1,
            maximum: 0,
        }
    );
    assert_eq!(
        ReferenceIndex::build(
            vec![ContigInput::new(Vec::new(), normalized(b""))],
            ReferenceBuildLimits::MAX,
        )
        .unwrap_err(),
        ReferenceBuildError::EmptyContigName { contig_ordinal: 0 }
    );
    assert_eq!(
        ReferenceIndex::build(
            vec![
                ContigInput::new(b"same".to_vec(), normalized(b"A")),
                ContigInput::new(b"same".to_vec(), normalized(b"")),
            ],
            ReferenceBuildLimits::MAX,
        )
        .unwrap_err(),
        ReferenceBuildError::DuplicateContigName {
            first_ordinal: 0,
            duplicate_ordinal: 1,
        }
    );
    assert_eq!(
        ReferenceIndex::build(
            vec![
                ContigInput::new(b"first".to_vec(), normalized(b"")),
                ContigInput::new(Vec::new(), normalized(b"A")),
            ],
            ReferenceBuildLimits::MAX,
        )
        .unwrap_err(),
        ReferenceBuildError::EmptyContigSequence { contig_ordinal: 0 }
    );
    assert_eq!(
        ReferenceIndex::build(
            vec![
                ContigInput::new(b"a".to_vec(), normalized(b"A")),
                ContigInput::new(b"b".to_vec(), normalized(b"C")),
                ContigInput::new(b"a".to_vec(), normalized(b"G")),
                ContigInput::new(b"a".to_vec(), normalized(b"T")),
            ],
            ReferenceBuildLimits::MAX,
        )
        .unwrap_err(),
        ReferenceBuildError::DuplicateContigName {
            first_ordinal: 0,
            duplicate_ordinal: 2,
        }
    );
    assert_eq!(
        ReferenceIndex::build(
            vec![
                ContigInput::new(b"a".to_vec(), normalized(b"A")),
                ContigInput::new(b"b".to_vec(), normalized(b"T")),
            ],
            ReferenceBuildLimits::MAX.with_max_total_reference_bases(1),
        )
        .unwrap_err(),
        ReferenceBuildError::LimitExceeded {
            resource: ReferenceResource::TotalReferenceBases,
            requested: 2,
            maximum: 1,
        }
    );
    assert_eq!(
        ReferenceIndex::build(
            vec![ContigInput::new(b"runs".to_vec(), normalized(b"ANANA"))],
            ReferenceBuildLimits::MAX.with_max_canonical_runs(2),
        )
        .unwrap_err(),
        ReferenceBuildError::LimitExceeded {
            resource: ReferenceResource::CanonicalRuns,
            requested: 3,
            maximum: 2,
        }
    );

    let retained = named_build().metrics().estimated_retained_fm_bytes();
    let limit_cases = [
        (
            ReferenceBuildLimits::MAX.with_max_contigs(2),
            ReferenceBuildError::LimitExceeded {
                resource: ReferenceResource::Contigs,
                requested: 3,
                maximum: 2,
            },
        ),
        (
            ReferenceBuildLimits::MAX.with_max_total_name_bytes(15),
            ReferenceBuildError::LimitExceeded {
                resource: ReferenceResource::TotalNameBytes,
                requested: 16,
                maximum: 15,
            },
        ),
        (
            ReferenceBuildLimits::MAX.with_max_total_reference_bases(12),
            ReferenceBuildError::LimitExceeded {
                resource: ReferenceResource::TotalReferenceBases,
                requested: 13,
                maximum: 12,
            },
        ),
        (
            ReferenceBuildLimits::MAX.with_max_canonical_runs(2),
            ReferenceBuildError::LimitExceeded {
                resource: ReferenceResource::CanonicalRuns,
                requested: 3,
                maximum: 2,
            },
        ),
        (
            ReferenceBuildLimits::MAX.with_max_suffix_rows_per_lane(5),
            ReferenceBuildError::SuffixRowsPerLaneLimitExceeded {
                requested: 6,
                maximum: 5,
                contig_ordinal: 0,
                run_start: 0,
            },
        ),
        (
            ReferenceBuildLimits::MAX.with_max_lanes(11),
            ReferenceBuildError::LimitExceeded {
                resource: ReferenceResource::Lanes,
                requested: 12,
                maximum: 11,
            },
        ),
        (
            ReferenceBuildLimits::MAX.with_max_projected_bases(39),
            ReferenceBuildError::LimitExceeded {
                resource: ReferenceResource::ProjectedBases,
                requested: 40,
                maximum: 39,
            },
        ),
        (
            ReferenceBuildLimits::MAX.with_max_projected_suffix_rows(51),
            ReferenceBuildError::LimitExceeded {
                resource: ReferenceResource::ProjectedSuffixRows,
                requested: 52,
                maximum: 51,
            },
        ),
        (
            ReferenceBuildLimits::MAX.with_max_estimated_retained_fm_bytes(retained - 1),
            ReferenceBuildError::LimitExceeded {
                resource: ReferenceResource::EstimatedRetainedFmBytes,
                requested: retained,
                maximum: retained - 1,
            },
        ),
    ];
    for (limits, expected) in limit_cases {
        assert_eq!(build_named_with_limits(limits).unwrap_err(), expected);
    }

    let index = named_build();
    assert_eq!(
        index
            .exact_search(
                BisulfiteStrand::OT,
                &normalized(b""),
                ReferenceQueryLimits::MAX,
            )
            .unwrap_err(),
        ReferenceQueryError::EmptyPattern
    );
    assert_eq!(
        index
            .exact_search(
                BisulfiteStrand::OT,
                &normalized(b"AN"),
                ReferenceQueryLimits::new(1, u64::MAX),
            )
            .unwrap_err(),
        ReferenceQueryError::PatternLimitExceeded {
            requested: 2,
            maximum: 1,
        }
    );
    assert_eq!(
        index
            .exact_search(
                BisulfiteStrand::OT,
                &normalized(b"AN"),
                ReferenceQueryLimits::new(2, u64::MAX),
            )
            .unwrap_err(),
        ReferenceQueryError::UnsearchableBase { offset: 1 }
    );
    assert_eq!(
        index
            .exact_search(
                BisulfiteStrand::OT,
                &normalized(b"GT"),
                ReferenceQueryLimits::new(2, 1),
            )
            .unwrap_err(),
        ReferenceQueryError::HitLimitExceeded {
            requested: 2,
            maximum: 1,
        }
    );

    assert_eq!(
        index.contig_id(3).unwrap_err(),
        ReferenceAccessError::ContigOrdinalOutOfBounds {
            ordinal: 3,
            contig_count: 3,
        }
    );
    let foreign_owner = build_catalog(&[(b"one", b"A")]);
    let foreign_id = index.contig_id(2).expect("source ordinal exists");
    assert_eq!(
        foreign_owner.resolve_contig(&foreign_id).unwrap_err(),
        ReferenceAccessError::ForeignContigId
    );
    let empty_matches = index
        .exact_search(
            BisulfiteStrand::OT,
            &normalized(b"CCCCCCCC"),
            ReferenceQueryLimits::MAX,
        )
        .expect("zero-hit artifact is still owner-bound");
    assert_eq!(
        foreign_owner.locate(&empty_matches).unwrap_err(),
        ReferenceLocateError::ForeignMatches
    );

    let display = ReferenceQueryError::UnsearchableBase { offset: 1 };
    assert_eq!(
        display.to_string(),
        "query contains unsearchable N at offset 1"
    );
    assert!(std::error::Error::source(&display).is_none());
}

fn assert_error_display_and_source<E>(
    error: &E,
    expected_display: &str,
    expected_source: Option<&str>,
) where
    E: std::error::Error,
{
    assert_eq!(error.to_string(), expected_display);
    assert_eq!(
        std::error::Error::source(error).map(ToString::to_string),
        expected_source.map(str::to_owned)
    );
}

#[allow(
    clippy::too_many_lines,
    reason = "the public diagnostic contract is exhaustively enumerated"
)]
#[test]
fn every_public_error_variant_formats_fields_and_sources_exactly() {
    let build_cases = [
        (
            ReferenceBuildError::EmptyReference,
            "reference contains no contigs".to_owned(),
            None,
        ),
        (
            ReferenceBuildError::CountNotRepresentable {
                resource: ReferenceResource::CanonicalBases,
                value: 17,
            },
            "CanonicalBases count 17 is not representable as u64".to_owned(),
            None,
        ),
        (
            ReferenceBuildError::ArithmeticOverflow {
                resource: ReferenceResource::ProjectedBases,
                operation: ReferenceArithmetic::Multiply,
                lhs: u64::MAX,
                rhs: 4,
            },
            format!(
                "ProjectedBases arithmetic {} Multiply 4 overflowed",
                u64::MAX
            ),
            None,
        ),
        (
            ReferenceBuildError::LimitExceeded {
                resource: ReferenceResource::CanonicalRuns,
                requested: 3,
                maximum: 2,
            },
            "CanonicalRuns value 3 exceeds configured maximum 2".to_owned(),
            None,
        ),
        (
            ReferenceBuildError::EmptyContigName { contig_ordinal: 5 },
            "contig 5 has an empty name".to_owned(),
            None,
        ),
        (
            ReferenceBuildError::DuplicateContigName {
                first_ordinal: 2,
                duplicate_ordinal: 7,
            },
            "contig 7 duplicates the exact name of contig 2".to_owned(),
            None,
        ),
        (
            ReferenceBuildError::EmptyContigSequence { contig_ordinal: 4 },
            "contig 4 has an empty sequence".to_owned(),
            None,
        ),
        (
            ReferenceBuildError::SuffixRowsPerLaneLimitExceeded {
                requested: 9,
                maximum: 8,
                contig_ordinal: 3,
                run_start: 11,
            },
            "run at contig 3:11 needs 9 suffix rows, exceeding 8".to_owned(),
            None,
        ),
        (
            ReferenceBuildError::AllocationSizeOverflow {
                allocation: ReferenceAllocation::RunMetadata,
                elements: 17,
                element_size: 8,
            },
            "cannot size RunMetadata: 17 elements of 8 bytes".to_owned(),
            None,
        ),
        (
            ReferenceBuildError::AllocationSizeOverflow {
                allocation: ReferenceAllocation::ProjectionScratch,
                elements: 19,
                element_size: 1,
            },
            "cannot size ProjectionScratch: 19 elements of 1 bytes".to_owned(),
            None,
        ),
        (
            ReferenceBuildError::AllocationFailed {
                allocation: ReferenceAllocation::RunMetadata,
                elements: 23,
            },
            "failed to reserve 23 elements for RunMetadata".to_owned(),
            None,
        ),
        (
            ReferenceBuildError::AllocationFailed {
                allocation: ReferenceAllocation::ProjectionScratch,
                elements: 29,
            },
            "failed to reserve 29 elements for ProjectionScratch".to_owned(),
            None,
        ),
        {
            let source = FmError::AllocationFailed {
                component: FmAllocation::SuffixArray,
                elements: 31,
            };
            (
                ReferenceBuildError::FmBuild {
                    contig_ordinal: 2,
                    run_start: 5,
                    strand: BisulfiteStrand::CTOT,
                    source,
                },
                "FM build failed for contig 2 run 5 lane CTOT: failed to reserve 31 elements for SuffixArray".to_owned(),
                Some("failed to reserve 31 elements for SuffixArray".to_owned()),
            )
        },
        (
            ReferenceBuildError::InternalInvariant {
                expected: 37,
                observed: 41,
            },
            "reference build invariant expected 37, observed 41".to_owned(),
            None,
        ),
    ];
    for (error, display, source) in &build_cases {
        assert_error_display_and_source(error, display, source.as_deref());
    }

    let access_cases = [
        (
            ReferenceAccessError::ForeignContigId,
            "contig identifier belongs to another reference instance",
        ),
        (
            ReferenceAccessError::ContigOrdinalOutOfBounds {
                ordinal: 7,
                contig_count: 3,
            },
            "contig ordinal 7 is outside catalog count 3",
        ),
    ];
    for (error, display) in &access_cases {
        assert_error_display_and_source(error, display, None);
    }

    let query_cases = [
        (
            ReferenceQueryError::PatternLengthNotRepresentable { pattern_len: 11 },
            "physical pattern length 11 is not representable as u64".to_owned(),
        ),
        (
            ReferenceQueryError::EmptyPattern,
            "exact-search pattern is empty".to_owned(),
        ),
        (
            ReferenceQueryError::PatternLimitExceeded {
                requested: 13,
                maximum: 8,
            },
            "pattern length 13 exceeds configured maximum 8".to_owned(),
        ),
        (
            ReferenceQueryError::UnsearchableBase { offset: 5 },
            "query contains unsearchable N at offset 5".to_owned(),
        ),
        (
            ReferenceQueryError::CountOverflow {
                counter: ReferenceQueryCounter::ExactHits,
                accumulated: u64::MAX,
                next: 1,
            },
            format!("ExactHits count {} plus 1 overflowed", u64::MAX),
        ),
        (
            ReferenceQueryError::HitLimitExceeded {
                requested: 17,
                maximum: 16,
            },
            "exact hit count 17 exceeds configured maximum 16".to_owned(),
        ),
        (
            ReferenceQueryError::AllocationSizeOverflow {
                allocation: ReferenceAllocation::ProjectedPattern,
                elements: 43,
                element_size: 1,
            },
            "cannot size ProjectedPattern: 43 elements of 1 bytes".to_owned(),
        ),
        (
            ReferenceQueryError::AllocationSizeOverflow {
                allocation: ReferenceAllocation::OpaqueMatches,
                elements: 47,
                element_size: 24,
            },
            "cannot size OpaqueMatches: 47 elements of 24 bytes".to_owned(),
        ),
        (
            ReferenceQueryError::AllocationFailed {
                allocation: ReferenceAllocation::ProjectedPattern,
                elements: 53,
            },
            "failed to reserve 53 elements for ProjectedPattern".to_owned(),
        ),
        (
            ReferenceQueryError::AllocationFailed {
                allocation: ReferenceAllocation::OpaqueMatches,
                elements: 59,
            },
            "failed to reserve 59 elements for OpaqueMatches".to_owned(),
        ),
        (
            ReferenceQueryError::InvariantMismatch {
                counter: ReferenceQueryCounter::NonemptyIntervals,
                expected: 61,
                observed: 67,
            },
            "NonemptyIntervals count pass produced 61, materialization produced 67".to_owned(),
        ),
        (
            ReferenceQueryError::CapacityInvariant {
                reserved: 71,
                materialized: 73,
            },
            "opaque-match reservation 71 cannot accept entry 73".to_owned(),
        ),
    ];
    for (error, display) in &query_cases {
        assert_error_display_and_source(error, display, None);
    }

    let mut locate_cases = vec![
        (
            ReferenceLocateError::ForeignMatches,
            "projected matches belong to another reference instance".to_owned(),
            None,
        ),
        (
            ReferenceLocateError::AllocationSizeOverflow {
                allocation: ReferenceAllocation::FinalHits,
                elements: 79,
                element_size: 40,
            },
            "cannot size FinalHits: 79 elements of 40 bytes".to_owned(),
            None,
        ),
        (
            ReferenceLocateError::AllocationFailed {
                allocation: ReferenceAllocation::FinalHits,
                elements: 83,
            },
            "failed to reserve 83 elements for FinalHits".to_owned(),
            None,
        ),
        {
            let source = FmError::AllocationFailed {
                component: FmAllocation::LocateResults,
                elements: 89,
            };
            (
                ReferenceLocateError::FmLocate {
                    contig_ordinal: 1,
                    run_start: 2,
                    strand: BisulfiteStrand::OB,
                    source,
                },
                "FM locate failed for contig 1 run 2 lane OB: failed to reserve 89 elements for LocateResults".to_owned(),
                Some("failed to reserve 89 elements for LocateResults".to_owned()),
            )
        },
        {
            let source = CoordinateError::InvertedInterval {
                domain: CoordinateDomain::Reference,
                start: 5,
                end: 4,
            };
            (
                ReferenceLocateError::CoordinateRecovery {
                    contig_ordinal: 1,
                    run_start: 2,
                    strand: BisulfiteStrand::CTOB,
                    source,
                },
                "coordinate recovery failed for contig 1 run 2 lane CTOB: Reference half-open interval [5, 4) is inverted".to_owned(),
                Some("Reference half-open interval [5, 4) is inverted".to_owned()),
            )
        },
        (
            ReferenceLocateError::CoordinateArithmetic {
                contig_ordinal: 3,
                run_start: 7,
                offset: 11,
                pattern_len: 13,
            },
            "coordinate arithmetic failed at contig 3 run 7, offset 11, pattern 13".to_owned(),
            None,
        ),
    ];
    for invariant in [
        ReferenceLocateInvariant::MissingRun,
        ReferenceLocateInvariant::OffsetCount,
        ReferenceLocateInvariant::TerminalSuffix,
        ReferenceLocateInvariant::RunBounds,
        ReferenceLocateInvariant::FinalHitCapacity,
        ReferenceLocateInvariant::FinalHitCount,
    ] {
        locate_cases.push((
            ReferenceLocateError::Invariant {
                invariant,
                expected: 97,
                observed: 101,
            },
            format!("{invariant:?} invariant expected 97, observed 101"),
            None,
        ));
    }
    for (error, display, source) in &locate_cases {
        assert_error_display_and_source(error, display, source.as_deref());
    }
}

macro_rules! assert_not_impl {
    ($type:ty, $forbidden:path) => {{
        trait AmbiguousIfImplemented<Disambiguator> {
            fn marker() {}
        }
        impl<T: ?Sized> AmbiguousIfImplemented<()> for T {}
        struct ForbiddenImplementation;
        impl<T: ?Sized + $forbidden> AmbiguousIfImplemented<ForbiddenImplementation> for T {}
        let _ = <$type as AmbiguousIfImplemented<_>>::marker;
    }};
}

fn linked_bsbit_index_rlibs() -> (std::path::PathBuf, Vec<std::path::PathBuf>) {
    let executable = std::env::current_exe().expect("test executable path is available");
    let dependency_dir = executable
        .parent()
        .expect("integration test executable is under target dependencies")
        .to_path_buf();
    let mut candidates = std::fs::read_dir(&dependency_dir)
        .expect("target dependency directory is readable")
        .filter_map(Result::ok)
        .map(|entry| entry.path())
        .filter(|path| {
            path.file_name().is_some_and(|name| {
                let name = name.to_string_lossy();
                name.starts_with("libbsbit_index-") && name.ends_with(".rlib")
            })
        })
        .collect::<Vec<_>>();
    // Prefer the artifact produced by the current Cargo invocation. A shared
    // target directory can contain many older feature-hash variants; probing
    // every compatible historical rlib makes this four-case API gate grow
    // without bound over the lifetime of a checkout.
    candidates.sort_unstable_by(|left, right| {
        let modified = |path: &std::path::Path| {
            std::fs::metadata(path)
                .and_then(|metadata| metadata.modified())
                .unwrap_or(std::time::SystemTime::UNIX_EPOCH)
        };
        modified(right)
            .cmp(&modified(left))
            .then_with(|| right.cmp(left))
    });
    assert!(
        !candidates.is_empty(),
        "compile-fail gate requires at least one linked bsbit_index rlib"
    );
    (dependency_dir, candidates)
}

fn assert_rust_source_fails(source: &str, required_diagnostics: &[&str]) {
    static NEXT_METADATA_ID: AtomicU64 = AtomicU64::new(0);

    let (dependency_dir, libraries) = linked_bsbit_index_rlibs();
    let compiler = std::env::var_os("RUSTC").unwrap_or_else(|| "rustc".into());
    let metadata_path = std::env::temp_dir().join(format!(
        "bsbit-reference-api-opacity-{}-{}.rmeta",
        std::process::id(),
        NEXT_METADATA_ID.fetch_add(1, Ordering::Relaxed)
    ));
    for library in libraries {
        let mut emit_metadata = std::ffi::OsString::from("metadata=");
        emit_metadata.push(&metadata_path);
        let mut child = Command::new(&compiler)
            .arg("--edition=2024")
            .arg("--crate-name=bsbit_reference_api_opacity_probe")
            .arg("--crate-type=lib")
            .arg("--emit")
            .arg(&emit_metadata)
            // rustc 1.94's annotated-snippet renderer can ICE while drawing a
            // private-field diagnostic for stdin.  The short renderer carries
            // the same error code and privacy text without that presentation-
            // only failure.
            .arg("--error-format=short")
            .arg("-")
            .arg("--extern")
            .arg(format!("bsbit_index={}", library.display()))
            .arg("-L")
            .arg(format!("dependency={}", dependency_dir.display()))
            .stdin(Stdio::piped())
            .stdout(Stdio::null())
            .stderr(Stdio::piped())
            .spawn()
            .expect("rustc must be available to the compile-fail gate");
        {
            let mut stdin = child.stdin.take().expect("rustc stdin is piped");
            stdin
                .write_all(source.as_bytes())
                .expect("compile-fail source fits rustc stdin");
        }
        let output = child
            .wait_with_output()
            .expect("compile-fail rustc process completes");
        match std::fs::remove_file(&metadata_path) {
            Ok(()) => {}
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => {}
            Err(error) => panic!(
                "compile-fail gate could not remove temporary metadata {}: {error}",
                metadata_path.display()
            ),
        }
        let stderr = String::from_utf8_lossy(&output.stderr);
        // A shared target directory can retain rlibs built with a different
        // feature/dependency hash.  rustc cannot link those candidates against
        // the dependency set of this integration-test process; they are not
        // valid probes of the current API and must not turn the opacity gate
        // into an order-dependent false failure.
        if stderr.contains("can't find crate for `bsbit_index`")
            || stderr.contains("could not find `reference` in `bsbit_index`")
            || stderr.contains("could not find `index` in `bsbit_index`")
        {
            continue;
        }
        assert!(
            !output.status.success(),
            "source unexpectedly compiled against {}; API opacity regressed",
            library.display()
        );
        for &diagnostic in required_diagnostics {
            assert!(
                stderr.contains(diagnostic),
                "compile failure against {} did not contain {diagnostic:?}:\n{stderr}",
                library.display()
            );
        }
        return;
    }
    panic!("compile-fail gate found no bsbit_index rlib compatible with the current test process");
}

#[test]
fn owner_bound_artifacts_have_no_value_equality_hash_or_order() {
    assert_not_impl!(ReferenceInstanceId, PartialEq);
    assert_not_impl!(ReferenceInstanceId, Eq);
    assert_not_impl!(ReferenceInstanceId, std::hash::Hash);
    assert_not_impl!(ReferenceInstanceId, Ord);
    assert_not_impl!(ContigId, PartialEq);
    assert_not_impl!(ContigId, Eq);
    assert_not_impl!(ContigId, std::hash::Hash);
    assert_not_impl!(ContigId, Ord);
    assert_not_impl!(ProjectedMatches, PartialEq);
    assert_not_impl!(ProjectedMatches, Eq);
    assert_not_impl!(ProjectedMatches, std::hash::Hash);
    assert_not_impl!(ProjectedMatches, Ord);
}

#[test]
fn compile_fail_gate_keeps_owner_bound_fields_private() {
    assert_rust_source_fails(
        r"
use std::hash::Hash;
use bsbit_index::reference::{ContigId, ProjectedMatches, ReferenceInstanceId};

fn requires_eq<T: Eq>() {}
fn requires_hash<T: Hash>() {}
fn unavailable<T>() -> T { panic!() }

pub fn trait_gate() {
    requires_eq::<ReferenceInstanceId>();
    requires_eq::<ContigId>();
    requires_eq::<ProjectedMatches>();
    requires_hash::<ReferenceInstanceId>();
    requires_hash::<ContigId>();
    requires_hash::<ProjectedMatches>();
}

pub fn forge_instance() -> ReferenceInstanceId {
    ReferenceInstanceId { owner: unavailable() }
}

pub fn forge_contig() -> ContigId {
    ContigId { owner: unavailable(), ordinal: 0 }
}

pub fn forge_matches() -> ProjectedMatches {
    ProjectedMatches {
        owner: unavailable(),
        strand: unavailable(),
        pattern_len: 1,
        exact_hit_count: 0,
        entries: unavailable(),
    }
}
",
        &[
            "ReferenceInstanceId",
            "ContigId",
            "ProjectedMatches",
            "Eq",
            "Hash",
            "private",
        ],
    );
}
