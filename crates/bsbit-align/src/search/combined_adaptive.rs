//! Shared combined-index search scheduling for one or two independent read lanes.
//!
//! A lane may be a paired-end mate or an unrelated single-end read. Pair
//! geometry and result classification deliberately remain outside this module.

use bsbit_core::alphabet::Base;
use bsbit_index::reference::ReferenceIndex;
use bsbit_index::storage::fm::{ProjectedBase, SearchBase};

use crate::AlignmentError;
use crate::alignment_policy::{
    CombinedSearchLimits, DEFAULT_MAXIMUM_SEED_ROUNDS, DEFAULT_SEARCH_LIMITS, EMPTY_SEED_STEP,
    INITIAL_MAXIMUM_SEED_ROUNDS, SEED_PROOF_ROUNDS,
};
use crate::library::ConversionPass;
use crate::read_mapping::{ReadCandidate, ungapped_distance};
use crate::read_mapping_limits::{MAX_READ_BASES, MIN_READ_BASES, MIN_SUFFIX_BASES};
use crate::search::combined_query::{CombinedSearchReferenceExt, CombinedSeedMatches};

/// Exact offset-seed rounds occupy bits zero through nine. Proof-kind flags
/// and the direct candidate's cached edit distance live in disjoint fields so
/// candidate deduplication can safely union independently observed evidence.
pub(crate) const SEED_ROUND_PROOF_MASK: u16 = (1 << SEED_PROOF_ROUNDS) - 1;
pub(crate) const FLEXIBLE_NOMINAL_PROOF: u16 = 1 << SEED_PROOF_ROUNDS;
pub(crate) const DIRECT_SINGLETON_PROOF: u16 = 1 << (SEED_PROOF_ROUNDS + 1);
const DIRECT_SINGLETON_DISTANCE_SHIFT: usize = SEED_PROOF_ROUNDS + 2;
const DIRECT_SINGLETON_DISTANCE_MASK: u16 = 0b111 << DIRECT_SINGLETON_DISTANCE_SHIFT;

const fn seed_round_proof(round: usize) -> u16 {
    debug_assert!(round < SEED_PROOF_ROUNDS);
    1_u16 << round
}

pub(crate) const fn direct_singleton_proof(distance: u8) -> u16 {
    DIRECT_SINGLETON_PROOF | ((distance as u16) << DIRECT_SINGLETON_DISTANCE_SHIFT)
}

pub(crate) const fn direct_singleton_distance(proof: u16) -> u8 {
    ((proof & DIRECT_SINGLETON_DISTANCE_MASK) >> DIRECT_SINGLETON_DISTANCE_SHIFT) as u8
}

#[allow(clippy::cast_possible_truncation)]
pub(crate) const fn seed_round_support(proof: u16) -> u8 {
    let rounds = (proof & SEED_ROUND_PROOF_MASK).count_ones();
    // The mask contains exactly ten bits, so this conversion is bounded.
    debug_assert!(rounds <= 10);
    rounds as u8
}

#[derive(Clone, Copy)]
pub(crate) struct DeferredCombinedSeed {
    pub(crate) matches: CombinedSeedMatches,
    pub(crate) offset: usize,
    pub(crate) round: usize,
}

type CombinedTwoLaneRoundSummary = ([u64; 2], [usize; 2], [bool; 2]);

#[derive(Clone, Copy)]
pub(crate) struct CombinedTwoLaneSearchState {
    pub(crate) located: [u64; 2],
    pub(crate) offsets: [usize; 2],
    pub(crate) active: [bool; 2],
    pub(crate) direct: [bool; 2],
    pub(crate) completed_rounds: usize,
    pub(crate) deferred: [[Option<DeferredCombinedSeed>; INITIAL_MAXIMUM_SEED_ROUNDS]; 2],
    pub(crate) deferred_len: [usize; 2],
    pub(crate) initialized: bool,
}

impl CombinedTwoLaneSearchState {
    pub(crate) const fn new() -> Self {
        Self {
            located: [0; 2],
            offsets: [0; 2],
            active: [true; 2],
            direct: [false; 2],
            completed_rounds: 0,
            deferred: [[None; INITIAL_MAXIMUM_SEED_ROUNDS]; 2],
            deferred_len: [0; 2],
            initialized: false,
        }
    }

    pub(crate) fn defer(&mut self, lane: usize, seed: DeferredCombinedSeed) {
        let offset = self.deferred_len[lane];
        if offset < self.deferred[lane].len() {
            self.deferred[lane][offset] = Some(seed);
            self.deferred_len[lane] += 1;
        }
    }
}

// Every argument is one explicit component of the bounded locate transaction.
#[allow(clippy::too_many_arguments)]
fn visit_combined_seed_round(
    reference: &ReferenceIndex,
    read: &[Base],
    conversion_pass: ConversionPass,
    round: usize,
    offset: usize,
    seed_matches: CombinedSeedMatches,
    limits: CombinedSearchLimits,
    candidates: &mut Vec<ReadCandidate>,
) -> Result<(u64, usize, bool), AlignmentError> {
    let matched_bases = usize::try_from(seed_matches.matched_bases())
        .map_err(|_| AlignmentError::LocatedCountOverflow)?;
    let hits = seed_matches.exact_hit_count();
    if hits != 1
        && (matched_bases < limits.minimum_multi_hit_seed_bases || hits > limits.maximum_seed_hits)
    {
        return Ok((0, matched_bases, false));
    }
    let before = candidates.len();
    let mut direct = false;
    let metrics = reference
        .visit_combined_seed(
            seed_matches,
            u64::try_from(offset).unwrap_or(u64::MAX),
            u64::try_from(read.len()).unwrap_or(u64::MAX),
            &mut |hit| {
                let Some(strand) = conversion_pass.relabel_combined_hit(hit.strand()) else {
                    return true;
                };
                let mut candidate = ReadCandidate {
                    contig_ordinal: hit.contig_ordinal(),
                    start: hit.start(),
                    strand,
                    proof_mask: FLEXIBLE_NOMINAL_PROOF | seed_round_proof(round),
                };
                if round == 0
                    && hits == 1
                    && let Some(distance) = ungapped_distance(reference, read, candidate)
                {
                    candidate.proof_mask = direct_singleton_proof(distance);
                    direct = true;
                }
                candidates.push(candidate);
                true
            },
        )
        .map_err(|_| AlignmentError::CombinedIndex)?;
    if direct {
        candidates.copy_within(before.., 0);
        candidates.truncate(candidates.len() - before);
    }
    Ok((metrics.located_coordinates(), matched_bases, direct))
}

pub(crate) fn combined_seed_round_is_locatable(
    seed_matches: CombinedSeedMatches,
    limits: CombinedSearchLimits,
) -> bool {
    let Ok(matched_bases) = usize::try_from(seed_matches.matched_bases()) else {
        return false;
    };
    let hits = seed_matches.exact_hit_count();
    hits == 1
        || (matched_bases >= limits.minimum_multi_hit_seed_bases
            && hits <= limits.maximum_seed_hits)
}

#[allow(clippy::too_many_arguments)]
fn visit_combined_seed_round_two_lanes(
    reference: &ReferenceIndex,
    reads: [&[Base]; 2],
    conversion_passes: [ConversionPass; 2],
    round: usize,
    offsets: [usize; 2],
    seed_matches: [CombinedSeedMatches; 2],
    mate1_candidates: &mut Vec<ReadCandidate>,
    mate2_candidates: &mut Vec<ReadCandidate>,
) -> Result<CombinedTwoLaneRoundSummary, AlignmentError> {
    let matched_bases = [
        usize::try_from(seed_matches[0].matched_bases())
            .map_err(|_| AlignmentError::LocatedCountOverflow)?,
        usize::try_from(seed_matches[1].matched_bases())
            .map_err(|_| AlignmentError::LocatedCountOverflow)?,
    ];
    let hits = [
        seed_matches[0].exact_hit_count(),
        seed_matches[1].exact_hit_count(),
    ];
    let before = [mate1_candidates.len(), mate2_candidates.len()];
    let mut direct = [false; 2];
    let metrics = reference
        .visit_combined_seed_two_lanes_complete(
            seed_matches,
            offsets.map(|offset| u64::try_from(offset).unwrap_or(u64::MAX)),
            reads.map(|read| u64::try_from(read.len()).unwrap_or(u64::MAX)),
            &mut |lane, hit| {
                let Some(strand) = conversion_passes[lane].relabel_combined_hit(hit.strand())
                else {
                    return;
                };
                let mut candidate = ReadCandidate {
                    contig_ordinal: hit.contig_ordinal(),
                    start: hit.start(),
                    strand,
                    proof_mask: FLEXIBLE_NOMINAL_PROOF | seed_round_proof(round),
                };
                if round == 0
                    && hits[lane] == 1
                    && let Some(distance) = ungapped_distance(reference, reads[lane], candidate)
                {
                    candidate.proof_mask = direct_singleton_proof(distance);
                    direct[lane] = true;
                }
                if lane == 0 {
                    mate1_candidates.push(candidate);
                } else {
                    mate2_candidates.push(candidate);
                }
            },
        )
        .map_err(|_| AlignmentError::CombinedIndex)?;
    if direct[0] {
        mate1_candidates.copy_within(before[0].., 0);
        mate1_candidates.truncate(mate1_candidates.len() - before[0]);
    }
    if direct[1] {
        mate2_candidates.copy_within(before[1].., 0);
        mate2_candidates.truncate(mate2_candidates.len() - before[1]);
    }
    Ok((
        metrics.map(bsbit_index::reference::ReferenceLocateMetrics::located_coordinates),
        matched_bases,
        direct,
    ))
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn start_combined_two_lane_search(
    reference: &ReferenceIndex,
    reads: [&[Base]; 2],
    reversed_projected: [&[ProjectedBase]; 2],
    first_seeds: [Option<CombinedSeedMatches>; 2],
    conversion_passes: [ConversionPass; 2],
    limits: CombinedSearchLimits,
    mate1_candidates: &mut Vec<ReadCandidate>,
    mate2_candidates: &mut Vec<ReadCandidate>,
) -> Result<CombinedTwoLaneSearchState, AlignmentError> {
    let mut state = CombinedTwoLaneSearchState::new();
    state.initialized = true;
    visit_combined_two_lane_search_rounds(
        reference,
        reads,
        reversed_projected,
        first_seeds,
        conversion_passes,
        limits,
        &mut state,
        mate1_candidates,
        mate2_candidates,
    )?;
    Ok(state)
}

// The two-lane wavefront intentionally advances both lanes in one loop so
// index queries can stay batched and their completion evidence stays aligned.
#[allow(
    clippy::needless_range_loop,
    clippy::too_many_arguments,
    clippy::too_many_lines
)]
fn visit_combined_two_lane_search_rounds(
    reference: &ReferenceIndex,
    reads: [&[Base]; 2],
    reversed_projected: [&[ProjectedBase]; 2],
    first_seeds: [Option<CombinedSeedMatches>; 2],
    conversion_passes: [ConversionPass; 2],
    limits: CombinedSearchLimits,
    state: &mut CombinedTwoLaneSearchState,
    mate1_candidates: &mut Vec<ReadCandidate>,
    mate2_candidates: &mut Vec<ReadCandidate>,
) -> Result<(), AlignmentError> {
    for round in state.completed_rounds..limits.maximum_seed_rounds {
        let available = [
            reads[0].len().saturating_sub(state.offsets[0]),
            reads[1].len().saturating_sub(state.offsets[1]),
        ];
        for (active, &available_bases) in state.active.iter_mut().zip(&available) {
            *active &= available_bases >= MIN_SUFFIX_BASES;
        }
        if !state.active[0] && !state.active[1] {
            break;
        }
        state.completed_rounds = round + 1;
        let matches = if round == 0 {
            first_seeds
        } else if state.active[0] && state.active[1] {
            reference
                .combined_maximal_suffix_projected_two_lanes(
                    [
                        &reversed_projected[0][..available[0]],
                        &reversed_projected[1][..available[1]],
                    ],
                    MIN_SUFFIX_BASES,
                )
                .map_err(|_| AlignmentError::CombinedIndex)?
        } else {
            let mut scalar = [None, None];
            let lane = usize::from(!state.active[0]);
            scalar[lane] = reference
                .combined_maximal_suffix_projected(
                    &reversed_projected[lane][..available[lane]],
                    MIN_SUFFIX_BASES,
                )
                .map_err(|_| AlignmentError::CombinedIndex)?;
            scalar
        };

        if limits.maximum_seed_rounds < DEFAULT_MAXIMUM_SEED_ROUNDS {
            let default_limits = DEFAULT_SEARCH_LIMITS;
            let active_lanes = state.active;
            for (lane, (&active, &seed)) in active_lanes.iter().zip(&matches).enumerate() {
                if active
                    && let Some(seed) = seed
                    && !combined_seed_round_is_locatable(seed, limits)
                    && combined_seed_round_is_locatable(seed, default_limits)
                {
                    state.defer(
                        lane,
                        DeferredCombinedSeed {
                            matches: seed,
                            offset: state.offsets[lane],
                            round,
                        },
                    );
                }
            }
        }

        let mut consumed = [false; 2];
        if state.active[0]
            && state.active[1]
            && let [Some(first), Some(second)] = matches
            && combined_seed_round_is_locatable(first, limits)
            && combined_seed_round_is_locatable(second, limits)
        {
            let (rows, matched_bases, direct) = visit_combined_seed_round_two_lanes(
                reference,
                reads,
                conversion_passes,
                round,
                state.offsets,
                [first, second],
                mate1_candidates,
                mate2_candidates,
            )?;
            for lane in 0..2 {
                state.located[lane] = state.located[lane]
                    .checked_add(rows[lane])
                    .ok_or(AlignmentError::LocatedCountOverflow)?;
                state.offsets[lane] = state.offsets[lane]
                    .saturating_add((matched_bases[lane].saturating_mul(3) / 4).max(1));
                state.active[lane] &= !direct[lane];
                state.direct[lane] |= direct[lane];
                consumed[lane] = true;
            }
        }

        macro_rules! consume_lane {
            ($lane:literal, $candidates:expr, $relabel:expr) => {
                if state.active[$lane] && !consumed[$lane] {
                    if let Some(seed) = matches[$lane] {
                        let (rows, matched, direct) = visit_combined_seed_round(
                            reference,
                            reads[$lane],
                            $relabel,
                            round,
                            state.offsets[$lane],
                            seed,
                            limits,
                            $candidates,
                        )?;
                        state.located[$lane] = state.located[$lane]
                            .checked_add(rows)
                            .ok_or(AlignmentError::LocatedCountOverflow)?;
                        state.offsets[$lane] = state.offsets[$lane]
                            .saturating_add((matched.saturating_mul(3) / 4).max(1));
                        state.active[$lane] &= !direct;
                        state.direct[$lane] |= direct;
                    } else {
                        state.offsets[$lane] = state.offsets[$lane].saturating_add(EMPTY_SEED_STEP);
                    }
                }
            };
        }
        consume_lane!(0, mate1_candidates, conversion_passes[0]);
        consume_lane!(1, mate2_candidates, conversion_passes[1]);
    }
    Ok(())
}

pub(crate) fn continue_combined_two_lane_search(
    reference: &ReferenceIndex,
    reads: [&[Base]; 2],
    reversed_projected: [&[ProjectedBase]; 2],
    conversion_passes: [ConversionPass; 2],
    state: &mut CombinedTwoLaneSearchState,
    mate1_candidates: &mut Vec<ReadCandidate>,
    mate2_candidates: &mut Vec<ReadCandidate>,
) -> Result<[u64; 2], AlignmentError> {
    continue_combined_two_lane_search_with_limits(
        reference,
        reads,
        reversed_projected,
        conversion_passes,
        DEFAULT_SEARCH_LIMITS,
        false,
        state,
        mate1_candidates,
        mate2_candidates,
    )
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn continue_combined_two_lane_search_with_limits(
    reference: &ReferenceIndex,
    reads: [&[Base]; 2],
    reversed_projected: [&[ProjectedBase]; 2],
    conversion_passes: [ConversionPass; 2],
    limits: CombinedSearchLimits,
    complete_direct_frontier: bool,
    state: &mut CombinedTwoLaneSearchState,
    mate1_candidates: &mut Vec<ReadCandidate>,
    mate2_candidates: &mut Vec<ReadCandidate>,
) -> Result<[u64; 2], AlignmentError> {
    let before = state.located;
    if complete_direct_frontier {
        for lane in 0..2 {
            state.active[lane] |= state.direct[lane];
        }
    }
    for (lane, read) in reads.into_iter().enumerate() {
        if state.direct[lane] && !complete_direct_frontier {
            continue;
        }
        let candidates = if lane == 0 {
            &mut *mate1_candidates
        } else {
            &mut *mate2_candidates
        };
        for deferred in state.deferred[lane][..state.deferred_len[lane]]
            .iter()
            .flatten()
        {
            let (rows, _, direct) = visit_combined_seed_round(
                reference,
                read,
                conversion_passes[lane],
                deferred.round,
                deferred.offset,
                deferred.matches,
                limits,
                candidates,
            )?;
            state.located[lane] = state.located[lane]
                .checked_add(rows)
                .ok_or(AlignmentError::LocatedCountOverflow)?;
            state.direct[lane] |= direct;
            if direct {
                state.active[lane] = false;
                break;
            }
        }
    }
    state.deferred = [[None; INITIAL_MAXIMUM_SEED_ROUNDS]; 2];
    state.deferred_len = [0; 2];
    visit_combined_two_lane_search_rounds(
        reference,
        reads,
        reversed_projected,
        [None, None],
        conversion_passes,
        limits,
        state,
        mate1_candidates,
        mate2_candidates,
    )?;
    Ok([
        state.located[0].saturating_sub(before[0]),
        state.located[1].saturating_sub(before[1]),
    ])
}

pub(crate) fn prepare_combined_projection(
    read: &[Base],
    conversion_pass: ConversionPass,
    output: &mut [ProjectedBase; MAX_READ_BASES],
) -> Result<(), AlignmentError> {
    if !(MIN_READ_BASES..=MAX_READ_BASES).contains(&read.len()) {
        return Err(AlignmentError::UnsupportedReadLength { length: read.len() });
    }
    if conversion_pass.reverse_complement_query() {
        for (destination, &base) in output.iter_mut().zip(read) {
            *destination = combined_projected_base(base.complement());
        }
    } else {
        for (destination, &base) in output.iter_mut().zip(read.iter().rev()) {
            *destination = combined_projected_base(base);
        }
    }
    Ok(())
}

pub(crate) fn prepare_combined_search_projection(
    read: &[Base],
    conversion_pass: ConversionPass,
    output: &mut [SearchBase; MAX_READ_BASES],
) -> Result<(), AlignmentError> {
    if !(MIN_READ_BASES..=MAX_READ_BASES).contains(&read.len()) {
        return Err(AlignmentError::UnsupportedReadLength { length: read.len() });
    }
    if conversion_pass.reverse_complement_query() {
        for (destination, &base) in output.iter_mut().zip(read) {
            *destination = combined_search_base(base.complement());
        }
    } else {
        for (destination, &base) in output.iter_mut().zip(read.iter().rev()) {
            *destination = combined_search_base(base);
        }
    }
    Ok(())
}

const fn combined_projected_base(base: Base) -> ProjectedBase {
    match base {
        Base::C | Base::T => ProjectedBase::T,
        Base::G => ProjectedBase::G,
        _ => ProjectedBase::A,
    }
}

const fn combined_search_base(base: Base) -> SearchBase {
    match combined_projected_base(base) {
        ProjectedBase::A => SearchBase::A,
        ProjectedBase::G => SearchBase::G,
        ProjectedBase::T => SearchBase::T,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn direct_distance_and_seed_rounds_survive_candidate_evidence_union() {
        let direct = direct_singleton_proof(2);
        let rediscovered = FLEXIBLE_NOMINAL_PROOF | seed_round_proof(0) | seed_round_proof(7);
        let combined = direct | rediscovered;

        assert_ne!(combined & DIRECT_SINGLETON_PROOF, 0);
        assert_ne!(combined & FLEXIBLE_NOMINAL_PROOF, 0);
        assert_eq!(direct_singleton_distance(combined), 2);
        assert_eq!(seed_round_support(combined), 2);
    }
}
