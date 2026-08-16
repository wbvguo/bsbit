//! Independent raw-byte oracle for Level 3 scalar extension.
//!
//! This module imports no implementation sequence, strand, coordinate, candidate,
//! distance, CIGAR, extension, ordering, or result type.

#![allow(dead_code)]

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum OracleStrand {
    Ot,
    Ob,
    Ctot,
    Ctob,
}

impl OracleStrand {
    const fn is_reverse(self) -> bool {
        matches!(self, Self::Ob | Self::Ctot)
    }

    const fn is_top(self) -> bool {
        matches!(self, Self::Ot | Self::Ctot)
    }
}

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct OraclePlacement {
    pub(crate) start: usize,
    pub(crate) end: usize,
    pub(crate) distance: u64,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct OracleWindow {
    pub(crate) start: usize,
    pub(crate) end: usize,
    pub(crate) intervals: u64,
    pub(crate) dp_cells: u64,
    pub(crate) best_distance: Option<u64>,
    pub(crate) placements: Vec<OraclePlacement>,
}

pub(crate) fn candidate_window_best(
    reference: &[u8],
    query: &[u8],
    strand: OracleStrand,
    diagonal: i128,
    budget: u64,
) -> OracleWindow {
    assert!(!query.is_empty());
    let reference_length = i128::try_from(reference.len()).expect("bounded reference length");
    let query_length = i128::try_from(query.len()).expect("bounded query length");
    let budget_signed = i128::from(budget);
    let start = (diagonal - budget_signed).clamp(0, reference_length);
    let end = (diagonal + query_length + budget_signed).clamp(0, reference_length);
    let start = usize::try_from(start).expect("clipped start fits usize");
    let end = usize::try_from(end.max(i128::try_from(start).expect("start fits i128")))
        .expect("clipped end fits usize");
    interval_best(reference, query, strand, start, end, budget)
}

pub(crate) fn whole_contig_best(
    reference: &[u8],
    query: &[u8],
    strand: OracleStrand,
    budget: u64,
) -> OracleWindow {
    assert!(!query.is_empty());
    interval_best(reference, query, strand, 0, reference.len(), budget)
}

pub(crate) fn whole_contig_passing(
    reference: &[u8],
    query: &[u8],
    strand: OracleStrand,
    budget: u64,
) -> Vec<OraclePlacement> {
    assert!(!query.is_empty());
    let oriented_query = if strand.is_reverse() {
        reverse_complement(query)
    } else {
        query.to_vec()
    };
    let minimum_length = query.len().saturating_sub(to_usize(budget)).max(1);
    let maximum_length = query
        .len()
        .saturating_add(to_usize(budget))
        .min(reference.len());
    let mut placements = Vec::new();
    if minimum_length > maximum_length {
        return placements;
    }
    for start in 0..reference.len() {
        let remaining = reference.len() - start;
        if remaining < minimum_length {
            break;
        }
        for length in minimum_length..=maximum_length.min(remaining) {
            let end = start + length;
            let distance = full_matrix_distance(&reference[start..end], &oriented_query, strand);
            if distance <= budget {
                placements.push(OraclePlacement {
                    start,
                    end,
                    distance,
                });
            }
        }
    }
    placements
}

fn interval_best(
    reference: &[u8],
    query: &[u8],
    strand: OracleStrand,
    start: usize,
    end: usize,
    budget: u64,
) -> OracleWindow {
    assert!(start <= end);
    assert!(end <= reference.len());
    let oriented_query = if strand.is_reverse() {
        reverse_complement(query)
    } else {
        query.to_vec()
    };
    let minimum_length = query.len().saturating_sub(to_usize(budget)).max(1);
    let maximum_length = query
        .len()
        .saturating_add(to_usize(budget))
        .min(end - start);
    let mut intervals = 0_u64;
    let mut dp_cells = 0_u64;
    let mut best_distance = None;
    let mut placements = Vec::new();
    if minimum_length <= maximum_length {
        for local_start in 0..end - start {
            let remaining = end - start - local_start;
            if remaining < minimum_length {
                break;
            }
            for length in minimum_length..=maximum_length.min(remaining) {
                intervals += 1;
                dp_cells += to_u64((length + 1) * (query.len() + 1));
                let absolute_start = start + local_start;
                let absolute_end = absolute_start + length;
                let distance = full_matrix_distance(
                    &reference[absolute_start..absolute_end],
                    &oriented_query,
                    strand,
                );
                if distance > budget {
                    continue;
                }
                match best_distance {
                    None => {
                        best_distance = Some(distance);
                        placements.clear();
                    }
                    Some(current) if distance < current => {
                        best_distance = Some(distance);
                        placements.clear();
                    }
                    Some(current) if distance > current => continue,
                    Some(_) => {}
                }
                placements.push(OraclePlacement {
                    start: absolute_start,
                    end: absolute_end,
                    distance,
                });
            }
        }
    }
    OracleWindow {
        start,
        end,
        intervals,
        dp_cells,
        best_distance,
        placements,
    }
}

fn full_matrix_distance(reference: &[u8], query: &[u8], strand: OracleStrand) -> u64 {
    let mut matrix = vec![vec![0_u64; query.len() + 1]; reference.len() + 1];
    for (index, cell) in matrix[0].iter_mut().enumerate() {
        *cell = to_u64(index);
    }
    for (index, row) in matrix.iter_mut().enumerate() {
        row[0] = to_u64(index);
    }
    for reference_index in 1..=reference.len() {
        for query_index in 1..=query.len() {
            let deletion = matrix[reference_index - 1][query_index] + 1;
            let insertion = matrix[reference_index][query_index - 1] + 1;
            let diagonal = matrix[reference_index - 1][query_index - 1]
                + relation_cost(
                    reference[reference_index - 1],
                    query[query_index - 1],
                    strand,
                );
            matrix[reference_index][query_index] = deletion.min(insertion).min(diagonal);
        }
    }
    matrix[reference.len()][query.len()]
}

const fn relation_cost(reference: u8, query: u8, strand: OracleStrand) -> u64 {
    if reference == b'N' || query == b'N' {
        1
    } else if reference == query
        || (strand.is_top() && reference == b'C' && query == b'T')
        || (!strand.is_top() && reference == b'G' && query == b'A')
    {
        0
    } else {
        1
    }
}

pub(crate) fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|base| complement(*base))
        .collect()
}

const fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'N' => b'N',
        _ => panic!("oracle input must be A/C/G/T/N"),
    }
}

pub(crate) fn enumerate_strings(alphabet: &[u8], maximum_length: usize) -> Vec<Vec<u8>> {
    let mut output = Vec::new();
    for length in 1..=maximum_length {
        let count = alphabet
            .len()
            .pow(u32::try_from(length).expect("bounded exponent"));
        for mut code in 0..count {
            let mut sequence = vec![alphabet[0]; length];
            for slot in (0..length).rev() {
                sequence[slot] = alphabet[code % alphabet.len()];
                code /= alphabet.len();
            }
            output.push(sequence);
        }
    }
    output
}

pub(crate) fn to_u64(value: usize) -> u64 {
    u64::try_from(value).expect("bounded oracle value fits u64")
}

fn to_usize(value: u64) -> usize {
    usize::try_from(value).unwrap_or(usize::MAX)
}
