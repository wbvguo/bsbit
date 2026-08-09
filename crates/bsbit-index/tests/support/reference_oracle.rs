//! Raw-byte scientific oracle for projected-reference integration tests.
//!
//! This module intentionally does not import production reference, FM,
//! projection, run-splitting, strand-table, or coordinate-recovery code.

use core::ops::Range;

pub(crate) const CANONICAL: [u8; 4] = [b'A', b'C', b'G', b'T'];
pub(crate) const REFERENCE_ALPHABET: [u8; 5] = [b'A', b'C', b'G', b'T', b'N'];

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) enum OracleStrand {
    Ot,
    Ob,
    Ctot,
    Ctob,
}

impl OracleStrand {
    pub(crate) const ALL: [Self; 4] = [Self::Ot, Self::Ob, Self::Ctot, Self::Ctob];

    pub(crate) const fn label(self) -> &'static str {
        match self {
            Self::Ot => "OT",
            Self::Ob => "OB",
            Self::Ctot => "CTOT",
            Self::Ctob => "CTOB",
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct OracleContig<'a> {
    pub(crate) name: &'a [u8],
    pub(crate) sequence: &'a [u8],
}

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct OracleHit {
    pub(crate) contig_ordinal: u64,
    pub(crate) start: u64,
    pub(crate) end: u64,
    pub(crate) strand: OracleStrand,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct OriginCase {
    pub(crate) raw_pattern: Vec<u8>,
    pub(crate) expected_origin: OracleHit,
}

pub(crate) fn canonical_runs(sequence: &[u8]) -> Vec<Range<usize>> {
    assert_reference_bytes(sequence);
    let mut runs = Vec::new();
    let mut run_start = 0;
    for (offset, &base) in sequence.iter().enumerate() {
        if base == b'N' {
            if run_start < offset {
                runs.push(run_start..offset);
            }
            run_start = offset + 1;
        }
    }
    if run_start < sequence.len() {
        runs.push(run_start..sequence.len());
    }
    runs
}

pub(crate) fn direct_search(
    catalog: &[OracleContig<'_>],
    strand: OracleStrand,
    raw_pattern: &[u8],
) -> Vec<OracleHit> {
    assert!(!raw_pattern.is_empty(), "oracle query must be nonempty");
    assert_canonical(raw_pattern);

    let (reverse, source, observed) = raw_view_rule(strand);
    let projected_pattern = project(raw_pattern, source, observed);
    let mut hits = Vec::new();

    for (contig_ordinal, contig) in catalog.iter().enumerate() {
        assert!(
            !contig.name.is_empty(),
            "oracle contig name must be nonempty"
        );
        for run in canonical_runs(contig.sequence) {
            let original_run = &contig.sequence[run.clone()];
            let reference_view = if reverse {
                reverse_complement(original_run)
            } else {
                original_run.to_vec()
            };
            let projected_run = project(&reference_view, source, observed);
            let width = projected_pattern.len();
            if width > projected_run.len() {
                continue;
            }
            for lane_offset in 0..=projected_run.len() - width {
                if projected_run[lane_offset..lane_offset + width] != projected_pattern {
                    continue;
                }
                let (start, end) = if reverse {
                    (run.end - (lane_offset + width), run.end - lane_offset)
                } else {
                    (run.start + lane_offset, run.start + lane_offset + width)
                };
                hits.push(OracleHit {
                    contig_ordinal: to_u64(contig_ordinal),
                    start: to_u64(start),
                    end: to_u64(end),
                    strand,
                });
            }
        }
    }

    hits.sort_unstable();
    hits
}

/// Generates known-origin raw reads from the independent four-letter relation.
pub(crate) fn origin_cases(catalog: &[OracleContig<'_>]) -> Vec<OriginCase> {
    let mut cases = Vec::new();
    for strand in OracleStrand::ALL {
        let (reverse, source, observed) = scientific_axes(strand);
        for (contig_ordinal, contig) in catalog.iter().enumerate() {
            for run in canonical_runs(contig.sequence) {
                let original_run = &contig.sequence[run.clone()];
                for left in 0..original_run.len() {
                    for right in left + 1..=original_run.len() {
                        let reference = &original_run[left..right];
                        for oriented in legal_oriented_reads(reference, source, observed) {
                            let raw_pattern = if reverse {
                                reverse_complement(&oriented)
                            } else {
                                oriented
                            };
                            cases.push(OriginCase {
                                raw_pattern,
                                expected_origin: OracleHit {
                                    contig_ordinal: to_u64(contig_ordinal),
                                    start: to_u64(run.start + left),
                                    end: to_u64(run.start + right),
                                    strand,
                                },
                            });
                        }
                    }
                }
            }
        }
    }
    cases
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
            _ => unreachable!("canonical assertion already rejected this byte"),
        })
        .collect()
}

pub(crate) fn project(sequence: &[u8], source: u8, observed: u8) -> Vec<u8> {
    assert_canonical(sequence);
    assert!(CANONICAL.contains(&source));
    assert!(CANONICAL.contains(&observed));
    sequence
        .iter()
        .map(|&base| if base == source { observed } else { base })
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

pub(crate) fn to_u64(value: usize) -> u64 {
    u64::try_from(value).expect("test sizes fit u64")
}

// These scientific axes describe an oriented query. For a raw reverse-lane
// view, complementing source and observed derives the required dual conversion.
fn scientific_axes(strand: OracleStrand) -> (bool, u8, u8) {
    match strand {
        OracleStrand::Ot => (false, b'C', b'T'),
        OracleStrand::Ob => (true, b'G', b'A'),
        OracleStrand::Ctot => (true, b'C', b'T'),
        OracleStrand::Ctob => (false, b'G', b'A'),
    }
}

fn raw_view_rule(strand: OracleStrand) -> (bool, u8, u8) {
    let (reverse, mut source, mut observed) = scientific_axes(strand);
    if reverse {
        source = complement(source);
        observed = complement(observed);
    }
    (reverse, source, observed)
}

fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        _ => panic!("oracle complement received a noncanonical byte"),
    }
}

fn legal_oriented_reads(reference: &[u8], source: u8, observed: u8) -> Vec<Vec<u8>> {
    let mut reads = vec![Vec::new()];
    for &base in reference {
        let prior = reads.len();
        for index in 0..prior {
            reads[index].push(base);
            if base == source {
                let mut converted = reads[index].clone();
                *converted
                    .last_mut()
                    .expect("one base was just appended to the oracle read") = observed;
                reads.push(converted);
            }
        }
    }
    reads
}

fn assert_reference_bytes(sequence: &[u8]) {
    assert!(
        sequence
            .iter()
            .all(|base| REFERENCE_ALPHABET.contains(base)),
        "oracle reference contains a byte outside A/C/G/T/N"
    );
}

fn assert_canonical(sequence: &[u8]) {
    assert!(
        sequence.iter().all(|base| CANONICAL.contains(base)),
        "oracle canonical sequence contains a byte outside A/C/G/T"
    );
}
