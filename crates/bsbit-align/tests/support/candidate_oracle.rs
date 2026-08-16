//! Independent raw-byte oracle for Level 2C fixed-seed candidates.
//!
//! This module deliberately does not import implementation sequence, bisulfite,
//! reference, coordinate, candidate, ordering, or deduplication code.

use core::cmp::Ordering;
use core::ops::Range;

pub(crate) const CANONICAL: [u8; 4] = [b'A', b'C', b'G', b'T'];
pub(crate) const REFERENCE_ALPHABET: [u8; 5] = [b'A', b'C', b'G', b'T', b'N'];

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum OracleStrand {
    Ot,
    Ob,
    Ctot,
    Ctob,
}

impl OracleStrand {
    pub(crate) const ALL: [Self; 4] = [Self::Ot, Self::Ob, Self::Ctot, Self::Ctob];

    pub(crate) const fn rank(self) -> u8 {
        match self {
            Self::Ot => 0,
            Self::Ob => 1,
            Self::Ctot => 2,
            Self::Ctob => 3,
        }
    }

    pub(crate) const fn is_reverse(self) -> bool {
        matches!(self, Self::Ob | Self::Ctot)
    }
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct OracleContig<'a> {
    pub(crate) name: &'a [u8],
    pub(crate) sequence: &'a [u8],
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct OracleRequest {
    pub(crate) strand: OracleStrand,
    pub(crate) start: usize,
    pub(crate) end: usize,
}

impl OracleRequest {
    pub(crate) const fn new(strand: OracleStrand, start: usize, end: usize) -> Self {
        Self { strand, start, end }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct OracleEvidence {
    pub(crate) contig_ordinal: u64,
    pub(crate) diagonal: i128,
    pub(crate) strand: OracleStrand,
    pub(crate) request_ordinal: u64,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct OracleAnchor {
    pub(crate) contig_ordinal: u64,
    pub(crate) diagonal: i128,
    pub(crate) strand: OracleStrand,
    pub(crate) support: u64,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct OracleMetrics {
    pub(crate) request_count: u64,
    pub(crate) total_seed_bases: u64,
    pub(crate) total_exact_hits: u64,
    pub(crate) matched_intervals: u64,
    pub(crate) unique_candidates: u64,
    pub(crate) duplicate_evidence: u64,
    pub(crate) maximum_support: u64,
    pub(crate) zero_hit_requests: u64,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct OracleSnapshot {
    pub(crate) anchors: Vec<OracleAnchor>,
    pub(crate) metrics: OracleMetrics,
}

pub(crate) fn candidate_snapshot(
    catalog: &[OracleContig<'_>],
    query: &[u8],
    requests: &[OracleRequest],
) -> OracleSnapshot {
    assert_canonical(query);
    let mut evidence = Vec::new();
    let mut matched_intervals = 0_u64;
    let mut zero_hit_requests = 0_u64;
    let mut total_seed_bases = 0_u64;

    for (request_ordinal, request) in requests.iter().enumerate() {
        assert!(request.start < request.end, "oracle seed must be nonempty");
        assert!(request.end <= query.len(), "oracle seed must fit the query");
        total_seed_bases = total_seed_bases
            .checked_add(to_u64(request.end - request.start))
            .expect("bounded oracle seed total fits u64");
        let (mut request_evidence, request_intervals) =
            direct_request_evidence(catalog, query, *request, to_u64(request_ordinal));
        if request_evidence.is_empty() {
            zero_hit_requests += 1;
        }
        matched_intervals = matched_intervals
            .checked_add(request_intervals)
            .expect("bounded oracle interval count fits u64");
        evidence.append(&mut request_evidence);
    }

    let anchors = group_evidence(evidence.clone());
    let total_exact_hits = to_u64(evidence.len());
    let unique_candidates = to_u64(anchors.len());
    let maximum_support = anchors
        .iter()
        .map(|anchor| anchor.support)
        .max()
        .unwrap_or(0);
    OracleSnapshot {
        anchors,
        metrics: OracleMetrics {
            request_count: to_u64(requests.len()),
            total_seed_bases,
            total_exact_hits,
            matched_intervals,
            unique_candidates,
            duplicate_evidence: total_exact_hits - unique_candidates,
            maximum_support,
            zero_hit_requests,
        },
    }
}

pub(crate) fn direct_request_evidence(
    catalog: &[OracleContig<'_>],
    query: &[u8],
    request: OracleRequest,
    request_ordinal: u64,
) -> (Vec<OracleEvidence>, u64) {
    assert_canonical(query);
    assert!(request.start < request.end);
    assert!(request.end <= query.len());
    let raw_seed = &query[request.start..request.end];
    let projected_seed = project_raw_view(raw_seed, request.strand);
    let oriented_start = if request.strand.is_reverse() {
        query.len() - request.end
    } else {
        request.start
    };
    let mut evidence = Vec::new();
    let mut matched_intervals = 0_u64;

    for (contig_ordinal, contig) in catalog.iter().enumerate() {
        assert!(
            !contig.name.is_empty(),
            "oracle contig name must be nonempty"
        );
        assert_reference(contig.sequence);
        for run in maximal_canonical_runs(contig.sequence) {
            let original = &contig.sequence[run.clone()];
            let raw_view = if request.strand.is_reverse() {
                reverse_complement(original)
            } else {
                original.to_vec()
            };
            let projected_reference = project_raw_view(&raw_view, request.strand);
            if projected_seed.len() > projected_reference.len() {
                continue;
            }
            let mut run_matched = false;
            for lane_offset in 0..=projected_reference.len() - projected_seed.len() {
                let lane_end = lane_offset + projected_seed.len();
                if projected_reference[lane_offset..lane_end] != projected_seed {
                    continue;
                }
                run_matched = true;
                let forward_start = if request.strand.is_reverse() {
                    run.end - lane_end
                } else {
                    run.start + lane_offset
                };
                evidence.push(OracleEvidence {
                    contig_ordinal: to_u64(contig_ordinal),
                    diagonal: to_i128(forward_start) - to_i128(oriented_start),
                    strand: request.strand,
                    request_ordinal,
                });
            }
            if run_matched {
                matched_intervals += 1;
            }
        }
    }

    (evidence, matched_intervals)
}

pub(crate) fn group_evidence(mut evidence: Vec<OracleEvidence>) -> Vec<OracleAnchor> {
    evidence.sort_unstable_by(compare_evidence);
    let mut anchors = Vec::new();
    let mut cursor = 0;
    while cursor < evidence.len() {
        let first = evidence[cursor];
        let mut support = 0_u64;
        let mut prior_request = None;
        while cursor < evidence.len() && same_key(evidence[cursor], first) {
            let current_request = evidence[cursor].request_ordinal;
            assert_ne!(
                prior_request,
                Some(current_request),
                "one request produced duplicate evidence for one candidate key"
            );
            prior_request = Some(current_request);
            support += 1;
            cursor += 1;
        }
        anchors.push(OracleAnchor {
            contig_ordinal: first.contig_ordinal,
            diagonal: first.diagonal,
            strand: first.strand,
            support,
        });
    }
    anchors
}

pub(crate) fn maximal_canonical_runs(sequence: &[u8]) -> Vec<Range<usize>> {
    assert_reference(sequence);
    let mut runs = Vec::new();
    let mut start = 0;
    for (offset, &base) in sequence.iter().enumerate() {
        if base == b'N' {
            if start < offset {
                runs.push(start..offset);
            }
            start = offset + 1;
        }
    }
    if start < sequence.len() {
        runs.push(start..sequence.len());
    }
    runs
}

pub(crate) fn project_raw_view(sequence: &[u8], strand: OracleStrand) -> Vec<u8> {
    assert_canonical(sequence);
    sequence
        .iter()
        .map(|&base| match strand {
            OracleStrand::Ot | OracleStrand::Ob if base == b'C' => b'T',
            OracleStrand::Ctot | OracleStrand::Ctob if base == b'G' => b'A',
            _ => base,
        })
        .collect()
}

pub(crate) fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    assert_canonical(sequence);
    sequence
        .iter()
        .rev()
        .map(|base| match base {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => unreachable!("canonical assertion rejected this byte"),
        })
        .collect()
}

pub(crate) fn enumerate_strings(alphabet: &[u8], maximum_length: usize) -> Vec<Vec<u8>> {
    let mut all = vec![Vec::new()];
    let mut frontier = vec![Vec::new()];
    for _ in 0..maximum_length {
        let mut next = Vec::new();
        for prefix in &frontier {
            for &base in alphabet {
                let mut value = prefix.clone();
                value.push(base);
                next.push(value);
            }
        }
        all.extend(next.iter().cloned());
        frontier = next;
    }
    all
}

fn compare_evidence(left: &OracleEvidence, right: &OracleEvidence) -> Ordering {
    left.contig_ordinal
        .cmp(&right.contig_ordinal)
        .then_with(|| left.diagonal.cmp(&right.diagonal))
        .then_with(|| left.strand.rank().cmp(&right.strand.rank()))
        .then_with(|| left.request_ordinal.cmp(&right.request_ordinal))
}

fn same_key(left: OracleEvidence, right: OracleEvidence) -> bool {
    left.contig_ordinal == right.contig_ordinal
        && left.diagonal == right.diagonal
        && left.strand == right.strand
}

fn assert_canonical(sequence: &[u8]) {
    assert!(
        sequence.iter().all(|base| CANONICAL.contains(base)),
        "oracle sequence contains a noncanonical byte"
    );
}

fn assert_reference(sequence: &[u8]) {
    assert!(
        sequence
            .iter()
            .all(|base| REFERENCE_ALPHABET.contains(base)),
        "oracle reference contains an unsupported byte"
    );
}

fn to_i128(value: usize) -> i128 {
    i128::try_from(value).expect("bounded oracle coordinate fits i128")
}

pub(crate) fn to_u64(value: usize) -> u64 {
    u64::try_from(value).expect("bounded oracle count fits u64")
}
