//! Independent raw-byte oracle for Level 2D proof-seed scheduling.
//!
//! This module deliberately imports no implementation sequence, coordinate,
//! distance, strand, seed, candidate, or reference code.

use core::ops::Range;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct OracleRequest {
    pub(crate) strand_rank: u8,
    pub(crate) start: usize,
    pub(crate) end: usize,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct OracleCertificate {
    pub(crate) query_bases: usize,
    pub(crate) max_edit_distance: usize,
    pub(crate) strand_count: usize,
    pub(crate) blocks: Vec<Range<usize>>,
    pub(crate) emitted_blocks: usize,
    pub(crate) omitted_unknown_blocks: usize,
    pub(crate) unknown_bases: usize,
    pub(crate) total_seed_bases: usize,
    pub(crate) minimum_block_bases: usize,
    pub(crate) maximum_block_bases: usize,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) enum OracleOutcome {
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

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub(crate) struct OracleFootprintStats {
    pub(crate) footprints: usize,
    pub(crate) displacement_assignments: usize,
    pub(crate) maximum_absolute_displacement: usize,
}

pub(crate) fn schedule(
    query: &[u8],
    supplied_strand_ranks: &[u8],
    max_edit_distance: usize,
) -> OracleOutcome {
    assert!(!supplied_strand_ranks.is_empty());
    assert!(supplied_strand_ranks.iter().all(|rank| *rank < 4));
    let mut strand_ranks = supplied_strand_ranks.to_vec();
    strand_ranks.sort_unstable();
    assert!(strand_ranks.windows(2).all(|pair| pair[0] != pair[1]));
    assert!(query.iter().all(|base| b"ACGTN".contains(base)));

    let mut unknown_bases = 0_usize;
    for base in query {
        if *base == b'N' {
            unknown_bases += 1;
        }
    }
    if unknown_bases > max_edit_distance {
        return OracleOutcome::NoAlignmentWithinBudget {
            query: query.to_vec(),
            strand_ranks,
            max_edit_distance,
            unknown_bases,
        };
    }
    if query.len() <= max_edit_distance {
        return OracleOutcome::SeedlessFallbackRequired {
            query: query.to_vec(),
            strand_ranks,
            max_edit_distance,
        };
    }

    let blocks = balanced_blocks(query.len(), max_edit_distance + 1);
    let mut emitted = Vec::new();
    let mut omitted_unknown_blocks = 0_usize;
    for block in &blocks {
        if query[block.clone()].contains(&b'N') {
            omitted_unknown_blocks += 1;
        } else {
            emitted.push(block.clone());
        }
    }
    assert!(!emitted.is_empty());

    let mut requests = Vec::new();
    for strand_rank in &strand_ranks {
        for block in &emitted {
            requests.push(OracleRequest {
                strand_rank: *strand_rank,
                start: block.start,
                end: block.end,
            });
        }
    }
    let seed_bases_per_strand = emitted.iter().map(Range::len).sum::<usize>();
    OracleOutcome::Certified {
        certificate: OracleCertificate {
            query_bases: query.len(),
            max_edit_distance,
            strand_count: strand_ranks.len(),
            emitted_blocks: emitted.len(),
            omitted_unknown_blocks,
            unknown_bases,
            total_seed_bases: seed_bases_per_strand * strand_ranks.len(),
            minimum_block_bases: blocks.iter().map(Range::len).min().expect("nonempty"),
            maximum_block_bases: blocks.iter().map(Range::len).max().expect("nonempty"),
            blocks,
        },
        requests,
    }
}

pub(crate) fn balanced_blocks(query_bases: usize, block_count: usize) -> Vec<Range<usize>> {
    assert!(block_count > 0);
    assert!(query_bases >= block_count);
    let mut blocks = Vec::with_capacity(block_count);
    for ordinal in 0..block_count {
        let start = ordinal * query_bases / block_count;
        let end = (ordinal + 1) * query_bases / block_count;
        assert!(start < end);
        blocks.push(start..end);
    }
    assert_eq!(blocks.first().expect("nonempty").start, 0);
    assert_eq!(blocks.last().expect("nonempty").end, query_bases);
    assert!(blocks.windows(2).all(|pair| pair[0].end == pair[1].start));
    blocks
}

/// Exhaustively checks abstract unit-edit footprints for the pigeonhole proof.
///
/// Query-base sites model substitutions or query-consuming gaps. Internal
/// boundary sites model reference-consuming gaps. An internal gap disrupts the
/// one block whose interior it splits and no block when it lies on a partition
/// boundary. Every N base is a mandatory unit edit.
pub(crate) fn exhaustive_edit_footprint_stats(
    query: &[u8],
    max_edit_distance: usize,
) -> OracleFootprintStats {
    assert!(query.len() > max_edit_distance);
    let blocks = balanced_blocks(query.len(), max_edit_distance + 1);
    let emitted = blocks
        .iter()
        .enumerate()
        .filter_map(|(ordinal, block)| (!query[block.clone()].contains(&b'N')).then_some(ordinal))
        .collect::<Vec<_>>();
    let mandatory_sites = query
        .iter()
        .enumerate()
        .filter_map(|(position, base)| (*base == b'N').then_some(position))
        .collect::<Vec<_>>();
    assert!(mandatory_sites.len() <= max_edit_distance);

    let site_count = query.len() + query.len().saturating_sub(1);
    let optional_sites = (0..site_count)
        .filter(|site| !mandatory_sites.contains(site))
        .collect::<Vec<_>>();
    let mut stats = OracleFootprintStats::default();
    let optional_budget = max_edit_distance - mandatory_sites.len();
    for extra_count in 0..=optional_budget {
        enumerate_combinations(
            &optional_sites,
            extra_count,
            0,
            &mut Vec::new(),
            &mut |chosen| {
                let edit_sites = mandatory_sites
                    .iter()
                    .chain(chosen.iter())
                    .copied()
                    .collect::<Vec<_>>();
                let mut affected = vec![false; blocks.len()];
                for site in &edit_sites {
                    if let Some(block) = affected_block(*site, query.len(), &blocks) {
                        affected[block] = true;
                    }
                }
                let clean_block = emitted
                    .iter()
                    .copied()
                    .find(|block| !affected[*block])
                    .expect("bounded edit footprint leaves an emitted proof block");
                enumerate_displacements(
                    &edit_sites,
                    0,
                    blocks[clean_block].start,
                    query.len(),
                    0,
                    &mut |displacement| {
                        let absolute = displacement.unsigned_abs();
                        assert!(absolute <= max_edit_distance);
                        stats.maximum_absolute_displacement =
                            stats.maximum_absolute_displacement.max(absolute);
                        stats.displacement_assignments += 1;
                    },
                );
                stats.footprints += 1;
            },
        );
    }
    stats
}

fn affected_block(site: usize, query_bases: usize, blocks: &[Range<usize>]) -> Option<usize> {
    if site < query_bases {
        return blocks
            .iter()
            .position(|block| block.start <= site && site < block.end);
    }
    let boundary = site - query_bases + 1;
    assert!(boundary < query_bases);
    blocks
        .iter()
        .position(|block| block.start < boundary && boundary < block.end)
}

fn enumerate_combinations(
    values: &[usize],
    remaining: usize,
    start: usize,
    selected: &mut Vec<usize>,
    visit: &mut impl FnMut(&[usize]),
) {
    if remaining == 0 {
        visit(selected);
        return;
    }
    if values.len().saturating_sub(start) < remaining {
        return;
    }
    for index in start..=values.len() - remaining {
        selected.push(values[index]);
        enumerate_combinations(values, remaining - 1, index + 1, selected, visit);
        selected.pop();
    }
}

fn enumerate_displacements(
    edit_sites: &[usize],
    index: usize,
    clean_block_start: usize,
    query_bases: usize,
    displacement: isize,
    visit: &mut impl FnMut(isize),
) {
    if index == edit_sites.len() {
        visit(displacement);
        return;
    }
    let site = edit_sites[index];
    let precedes_seed = if site < query_bases {
        site < clean_block_start
    } else {
        let boundary = site - query_bases + 1;
        boundary <= clean_block_start
    };
    for delta in [-1_isize, 0, 1] {
        let next = if precedes_seed {
            displacement + delta
        } else {
            displacement
        };
        enumerate_displacements(
            edit_sites,
            index + 1,
            clean_block_start,
            query_bases,
            next,
            visit,
        );
    }
}
