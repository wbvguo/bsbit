//! Candidate-frontier construction and ranked exact-block proofs.

use super::{
    AlignmentError, AlignmentOrientation, Base, BisulfiteStrand, CombinedSearchReferenceExt,
    FLEXIBLE_NOMINAL_PROOF, MAX_EDIT_DISTANCE, MAX_READ_BASES, MateRescueWindow,
    PairAlignmentMetrics, ProjectedBase, ProofBlock, RESCUE_BLOCKS, RankedBlockPartition,
    RankedBlockSeed, RankedBlockSelection, ReadCandidate, ReferenceIndex,
    SENSITIVE_ADAPTIVE_BOUNDARY_SHIFTS, SENSITIVE_ADAPTIVE_MIN_BLOCK_BASES,
    SENSITIVE_BALANCED_BOUNDARY_SHIFTS, SENSITIVE_PROOF_BLOCKS,
    SENSITIVE_SELECTIVE_UNMAPPED_MAX_RETAINED_HITS, SENSITIVE_SELECTIVE_UNMAPPED_MIN_RETAINED_HITS,
    SearchBase, strand_semantics,
};

pub(super) fn selective_unmapped_frontier_deepening_required(
    selections: [Option<RankedBlockSelection>; 2],
    ordinary_frontier_metrics: Option<PairAlignmentMetrics>,
) -> bool {
    let [Some(first), Some(second)] = selections else {
        return false;
    };
    let retained_hits = first.retained_hits.saturating_add(second.retained_hits);
    let both_incomplete = !first.complete && !second.complete;
    let inside_bounded_window = both_incomplete
        && (SENSITIVE_SELECTIVE_UNMAPPED_MIN_RETAINED_HITS
            ..SENSITIVE_SELECTIVE_UNMAPPED_MAX_RETAINED_HITS)
            .contains(&retained_hits);
    {
        let _ = ordinary_frontier_metrics;
        inside_bounded_window
    }
}

/// Emits a complete flexible candidate proof frontier inside one mate window.
///
/// `maximum_edit_distance + 1` disjoint query blocks guarantee that every
/// in-budget gapped placement leaves at least one exact projected block.  The
/// exact block can be displaced from the true alignment origin by at most the
/// edit budget, which is precisely the start domain covered by the flexible
/// verifier.  This is used only after a whole-genome rescue interval exceeds
/// its locate cap, so work is proportional to the paired genomic window rather
/// than to the global repeat count.
pub(super) fn append_local_flexible_proof_candidates(
    read: &[Base],
    contig: &[Base],
    window: MateRescueWindow,
    maximum_edit_distance: u8,
    candidates: &mut Vec<ReadCandidate>,
) {
    append_scalar_local_flexible_proof_candidates(
        read,
        contig,
        window,
        maximum_edit_distance,
        candidates,
    );
}

pub(super) fn append_scalar_local_flexible_proof_candidates(
    read: &[Base],
    contig: &[Base],
    window: MateRescueWindow,
    maximum_edit_distance: u8,
    candidates: &mut Vec<ReadCandidate>,
) {
    let block_count = usize::from(maximum_edit_distance) + 1;
    debug_assert!(block_count <= usize::from(MAX_EDIT_DISTANCE) + 1);
    let short = read.len() / block_count;
    let long_count = read.len() % block_count;
    let semantics = strand_semantics(window.strand);
    let edit_budget = u64::from(maximum_edit_distance);
    let mut query_start = 0_usize;

    for block_ordinal in 0..block_count {
        let block_len = short + usize::from(block_ordinal < long_count);
        let query_end = query_start + block_len;
        let mut query_code = 0_u128;
        let mut canonical = true;
        for oriented_position in query_start..query_end {
            let query_base = match semantics.orientation() {
                AlignmentOrientation::Forward => read[oriented_position],
                AlignmentOrientation::Reverse => {
                    read[read.len() - oriented_position - 1].complement()
                }
            };
            let Some(code) = projected_code(query_base, semantics.cytosine_strand()) else {
                canonical = false;
                break;
            };
            query_code = (query_code << 2) | u128::from(code);
        }
        if !canonical {
            query_start = query_end;
            continue;
        }

        let query_offset = u64::try_from(query_start).expect("bounded query offset fits u64");
        let scan_start = window
            .start
            .saturating_add(query_offset)
            .saturating_sub(edit_budget);
        let scan_end = window
            .end
            .saturating_add(query_offset)
            .saturating_add(edit_budget);
        let Some(mut position) = usize::try_from(scan_start).ok() else {
            query_start = query_end;
            continue;
        };
        let Some(last) = usize::try_from(scan_end).ok() else {
            query_start = query_end;
            continue;
        };
        let last = last.min(contig.len().saturating_sub(block_len));
        if position > last || position.saturating_add(block_len) > contig.len() {
            query_start = query_end;
            continue;
        }
        let bits = block_len * 2;
        let mask = if bits == u128::BITS as usize {
            u128::MAX
        } else {
            (1_u128 << bits) - 1
        };
        let mut reference_code = pack_projected(
            &contig[position..position + block_len],
            semantics.cytosine_strand(),
        );
        loop {
            if reference_code == query_code {
                let observed = u64::try_from(position).expect("reference position fits u64");
                if let Some(nominal) = observed.checked_sub(query_offset) {
                    let extended_start = window.start.saturating_sub(edit_budget);
                    let extended_end = window.end.saturating_add(edit_budget);
                    if (extended_start..=extended_end).contains(&nominal) {
                        candidates.push(ReadCandidate {
                            contig_ordinal: window.contig_ordinal,
                            start: nominal,
                            strand: window.strand,
                            proof_mask: FLEXIBLE_NOMINAL_PROOF | (1_u8 << block_ordinal),
                        });
                    }
                }
            }
            if position == last {
                break;
            }
            position += 1;
            let incoming = projected_code(
                contig[position + block_len - 1],
                semantics.cytosine_strand(),
            )
            .unwrap_or(3);
            reference_code = ((reference_code << 2) & mask) | u128::from(incoming);
        }
        query_start = query_end;
    }
}

/// Builds a rarity-first disjoint-block seed set for one mate.
///
/// The balanced `d + 1` partition is inspected before any suffix-array rows are
/// located. If it is empty or over limit, a small boundary-state dynamic
/// program chooses the complete shifted partition with the smallest
/// FM-interval sum. Every path remains a disjoint cover and independently
/// retains the pigeonhole completeness guarantee.
///
/// A repetitive block must not make a different, informative block unusable.
/// Retaining only the rarest prefix also keeps the global locate and verifier
/// work bounded by `SENSITIVE_RANKED_BLOCK_HITS` per attempted anchor.
pub(super) fn shifted_ranked_boundary(nominal: usize, shift: i8) -> Option<usize> {
    if shift.is_negative() {
        nominal.checked_sub(usize::from(shift.unsigned_abs()))
    } else {
        nominal.checked_add(usize::from(shift.unsigned_abs()))
    }
}

pub(super) fn ranked_block_boundaries(
    read_len: usize,
    block_count: usize,
    boundary_shifts: [i8; SENSITIVE_PROOF_BLOCKS - 1],
    minimum_block_bases: usize,
) -> Option<[usize; SENSITIVE_PROOF_BLOCKS + 1]> {
    if block_count == 0 || block_count > SENSITIVE_PROOF_BLOCKS {
        return None;
    }
    let short = read_len / block_count;
    let long_count = read_len % block_count;
    let mut boundaries = [0_usize; SENSITIVE_PROOF_BLOCKS + 1];
    for ordinal in 1..block_count {
        let nominal = short
            .checked_mul(ordinal)?
            .checked_add(long_count.min(ordinal))?;
        boundaries[ordinal] = shifted_ranked_boundary(nominal, boundary_shifts[ordinal - 1])?;
    }
    boundaries[block_count] = read_len;
    for ordinal in 0..block_count {
        if boundaries[ordinal + 1].checked_sub(boundaries[ordinal])? < minimum_block_bases {
            return None;
        }
    }
    Some(boundaries)
}

#[allow(clippy::too_many_arguments)]
pub(super) fn ranked_block_seed_for_range(
    reference: &ReferenceIndex,
    read: &[Base],
    projected_search: &[SearchBase; MAX_READ_BASES],
    mate1: bool,
    block_ordinal: usize,
    query_start: usize,
    query_end: usize,
) -> Result<Option<RankedBlockSeed>, AlignmentError> {
    let block_len = query_end - query_start;
    let source = if mate1 {
        &read[query_start..query_end]
    } else {
        &read[read.len() - query_end..read.len() - query_start]
    };
    if source.contains(&Base::N) {
        return Ok(None);
    }
    let reversed_start = read.len() - query_end;
    let reversed_end = read.len() - query_start;
    let pattern = &projected_search[reversed_start..reversed_end];
    let Some(matches) = reference
        .combined_exact_seed(pattern)
        .map_err(|_| AlignmentError::CombinedIndex)?
    else {
        return Ok(None);
    };
    if matches.matched_bases() != u64::try_from(block_len).expect("bounded block length fits u64") {
        return Err(AlignmentError::CombinedIndex);
    }
    Ok(Some(RankedBlockSeed {
        matches,
        query_offset: u64::try_from(query_start).expect("bounded query offset fits u64"),
        proof_mask: 1_u8 << block_ordinal,
    }))
}

#[allow(clippy::too_many_arguments)]
pub(super) fn fill_ranked_block_seed_partition(
    reference: &ReferenceIndex,
    read: &[Base],
    projected_search: &[SearchBase; MAX_READ_BASES],
    mate1: bool,
    block_count: usize,
    boundaries: &[usize; SENSITIVE_PROOF_BLOCKS + 1],
    output: &mut [Option<RankedBlockSeed>; SENSITIVE_PROOF_BLOCKS],
) -> Result<u64, AlignmentError> {
    output.fill(None);
    for (block_ordinal, slot) in output[..block_count].iter_mut().enumerate() {
        let query_start = boundaries[block_ordinal];
        let query_end = boundaries[block_ordinal + 1];
        *slot = ranked_block_seed_for_range(
            reference,
            read,
            projected_search,
            mate1,
            block_ordinal,
            query_start,
            query_end,
        )?;
    }
    Ok(output[..block_count]
        .iter()
        .flatten()
        .fold(0_u64, |total, seed| {
            total.saturating_add(seed.matches.exact_hit_count())
        }))
}

// This small dynamic program deliberately keeps edge evaluation and
// predecessor reconstruction together so completeness is easy to audit.
#[allow(clippy::too_many_lines)]
pub(super) fn optimal_adaptive_ranked_block_partition(
    reference: &ReferenceIndex,
    read: &[Base],
    projected_search: &[SearchBase; MAX_READ_BASES],
    mate1: bool,
    block_count: usize,
    maximum_ranked_block_hits: u64,
) -> Result<RankedBlockPartition, AlignmentError> {
    const STATES: usize = SENSITIVE_ADAPTIVE_BOUNDARY_SHIFTS.len();
    if !(2..=SENSITIVE_PROOF_BLOCKS).contains(&block_count)
        || read.len() < SENSITIVE_ADAPTIVE_MIN_BLOCK_BASES * block_count
    {
        return Ok(None);
    }

    let short = read.len() / block_count;
    let long_count = read.len() % block_count;
    let mut positions = [[0_usize; STATES]; SENSITIVE_PROOF_BLOCKS + 1];
    positions[block_count][0] = read.len();
    for (ordinal, layer) in positions.iter_mut().enumerate().take(block_count).skip(1) {
        let nominal = short * ordinal + long_count.min(ordinal);
        for (position, &shift) in layer.iter_mut().zip(&SENSITIVE_ADAPTIVE_BOUNDARY_SHIFTS) {
            *position = shifted_ranked_boundary(nominal, shift)
                .expect("qualified adaptive boundary remains inside the read");
        }
    }

    let mut edges = [[[None::<Option<RankedBlockSeed>>; STATES]; STATES]; SENSITIVE_PROOF_BLOCKS];
    let mut costs = [[u64::MAX; 2]; STATES];
    costs[0][0] = 0;
    let mut predecessors = [[[None::<(u8, u8)>; 2]; STATES]; SENSITIVE_PROOF_BLOCKS + 1];

    for block_ordinal in 0..block_count {
        let from_states = if block_ordinal == 0 { 1 } else { STATES };
        let to_states = if block_ordinal + 1 == block_count {
            1
        } else {
            STATES
        };
        let mut next_costs = [[u64::MAX; 2]; STATES];
        for from_state in 0..from_states {
            // Edge weights are non-negative and the caller only accepts a
            // complete partition within the locate cap. Once both hit states
            // exceed that cap, no continuation through this boundary can
            // become admissible, so avoid its FM interval queries entirely.
            if costs[from_state]
                .iter()
                .all(|&cost| cost > maximum_ranked_block_hits)
            {
                continue;
            }
            for to_state in 0..to_states {
                let query_start = positions[block_ordinal][from_state];
                let query_end = positions[block_ordinal + 1][to_state];
                if query_end <= query_start
                    || query_end - query_start < SENSITIVE_ADAPTIVE_MIN_BLOCK_BASES
                {
                    continue;
                }
                let seed = ranked_block_seed_for_range(
                    reference,
                    read,
                    projected_search,
                    mate1,
                    block_ordinal,
                    query_start,
                    query_end,
                )?;
                let hits = seed.map_or(0, |seed| seed.matches.exact_hit_count());
                edges[block_ordinal][from_state][to_state] = Some(seed);
                for (had_hits, &cost) in costs[from_state].iter().enumerate() {
                    if cost == u64::MAX {
                        continue;
                    }
                    let has_hits = usize::from(had_hits != 0 || hits != 0);
                    let candidate_cost = cost.saturating_add(hits);
                    if candidate_cost <= maximum_ranked_block_hits
                        && candidate_cost < next_costs[to_state][has_hits]
                    {
                        next_costs[to_state][has_hits] = candidate_cost;
                        predecessors[block_ordinal + 1][to_state][has_hits] = Some((
                            u8::try_from(from_state).expect("adaptive state fits u8"),
                            u8::try_from(had_hits).expect("hit state fits u8"),
                        ));
                    }
                }
            }
        }
        if next_costs
            .iter()
            .flatten()
            .all(|&cost| cost > maximum_ranked_block_hits)
        {
            return Ok(None);
        }
        costs = next_costs;
    }

    let all_hits = costs[0][1];
    if all_hits == u64::MAX {
        return Ok(None);
    }
    let mut seeds = [None; SENSITIVE_PROOF_BLOCKS];
    let mut state = 0_usize;
    let mut had_hits = 1_usize;
    for block_ordinal in (0..block_count).rev() {
        let (from_state, previous_had_hits) = predecessors[block_ordinal + 1][state][had_hits]
            .expect("reachable adaptive state has a predecessor");
        let from_state = usize::from(from_state);
        seeds[block_ordinal] =
            edges[block_ordinal][from_state][state].expect("selected adaptive edge was evaluated");
        state = from_state;
        had_hits = usize::from(previous_had_hits);
    }
    debug_assert_eq!(state, 0);
    debug_assert_eq!(had_hits, 0);
    Ok(Some((all_hits, seeds)))
}

pub(super) fn collect_ranked_block_seeds(
    reference: &ReferenceIndex,
    read: &[Base],
    reversed_projected: &[ProjectedBase],
    mate1: bool,
    maximum_edit_distance: u8,
    maximum_ranked_block_hits: u64,
    output: &mut [Option<RankedBlockSeed>; SENSITIVE_PROOF_BLOCKS],
) -> Result<Option<RankedBlockSelection>, AlignmentError> {
    output.fill(None);
    let block_count = usize::from(maximum_edit_distance) + 1;
    let mut projected_search = [SearchBase::A; MAX_READ_BASES];
    for (destination, &source) in projected_search.iter_mut().zip(reversed_projected) {
        *destination = match source {
            ProjectedBase::A => SearchBase::A,
            ProjectedBase::G => SearchBase::G,
            ProjectedBase::T => SearchBase::T,
        };
    }
    let balanced_boundaries = ranked_block_boundaries(
        read.len(),
        block_count,
        SENSITIVE_BALANCED_BOUNDARY_SHIFTS,
        0,
    )
    .expect("a bounded read has a balanced d+1 partition");
    let mut best_seeds = [None; SENSITIVE_PROOF_BLOCKS];
    let mut all_hits = fill_ranked_block_seed_partition(
        reference,
        read,
        &projected_search,
        mate1,
        block_count,
        &balanced_boundaries,
        &mut best_seeds,
    )?;
    if (all_hits == 0 || all_hits > maximum_ranked_block_hits)
        && let Some((candidate_hits, candidate_seeds)) = optimal_adaptive_ranked_block_partition(
            reference,
            read,
            &projected_search,
            mate1,
            block_count,
            maximum_ranked_block_hits,
        )?
        && candidate_hits <= maximum_ranked_block_hits
    {
        best_seeds = candidate_seeds;
        all_hits = candidate_hits;
    }
    *output = best_seeds;
    output[..block_count]
        .sort_by_key(|seed| seed.map_or(u64::MAX, |seed| seed.matches.exact_hit_count()));
    let mut retained_hits = 0_u64;
    for slot in &mut output[..block_count] {
        let Some(seed) = *slot else {
            break;
        };
        let next = retained_hits.saturating_add(seed.matches.exact_hit_count());
        if next > maximum_ranked_block_hits {
            *slot = None;
            continue;
        }
        retained_hits = next;
    }
    Ok((retained_hits != 0).then_some(RankedBlockSelection {
        retained_hits,
        complete: retained_hits == all_hits,
    }))
}

pub(super) fn append_ranked_block_candidates(
    reference: &ReferenceIndex,
    read_len: usize,
    mate1: bool,
    seeds: &[Option<RankedBlockSeed>; SENSITIVE_PROOF_BLOCKS],
    candidates: &mut Vec<ReadCandidate>,
) -> Result<u64, AlignmentError> {
    let query_len = u64::try_from(read_len).expect("bounded read length fits u64");
    let mut located_rows = 0_u64;
    for seed in seeds.iter().flatten() {
        let metrics = reference
            .visit_combined_seed(seed.matches, seed.query_offset, query_len, &mut |hit| {
                let strand = if mate1 {
                    hit.strand()
                } else {
                    match hit.strand() {
                        BisulfiteStrand::OT => BisulfiteStrand::CTOT,
                        BisulfiteStrand::OB => BisulfiteStrand::CTOB,
                        BisulfiteStrand::CTOT | BisulfiteStrand::CTOB => return true,
                    }
                };
                candidates.push(ReadCandidate {
                    contig_ordinal: hit.contig_ordinal(),
                    start: hit.start(),
                    strand,
                    proof_mask: FLEXIBLE_NOMINAL_PROOF | seed.proof_mask,
                });
                true
            })
            .map_err(|_| AlignmentError::CombinedIndex)?;
        located_rows = located_rows.saturating_add(metrics.located_coordinates());
    }
    Ok(located_rows)
}

pub(super) fn balanced_rescue_blocks(read_len: usize) -> [ProofBlock; RESCUE_BLOCKS] {
    let short = read_len / RESCUE_BLOCKS;
    let long_count = read_len % RESCUE_BLOCKS;
    let mut cursor = 0_usize;
    core::array::from_fn(|ordinal| {
        let length = short + usize::from(ordinal < long_count);
        let start = cursor;
        cursor += length;
        ProofBlock {
            query_start: u16::try_from(start).expect("bounded rescue block start fits u16"),
            query_end: u16::try_from(cursor).expect("bounded rescue block end fits u16"),
        }
    })
}

pub(super) fn pack_projected(
    bases: &[Base],
    strand: bsbit_core::bisulfite::CytosineStrand,
) -> u128 {
    bases.iter().fold(0_u128, |packed, &base| {
        (packed << 2) | u128::from(projected_code(base, strand).unwrap_or(3))
    })
}

pub(super) const fn projected_code(
    base: Base,
    strand: bsbit_core::bisulfite::CytosineStrand,
) -> Option<u8> {
    use bsbit_core::bisulfite::CytosineStrand::{Bottom, Top};
    match (strand, base) {
        (_, Base::A) | (Bottom, Base::G) => Some(0),
        (Top, Base::C | Base::T) | (Bottom, Base::C) => Some(1),
        (Top, Base::G) | (Bottom, Base::T) => Some(2),
        _ => None,
    }
}
