//! Shared combined-index search scheduling for one or two independent read lanes.
//!
//! A lane may be a paired-end mate or an unrelated single-end read. Pair
//! geometry and result classification deliberately remain outside this module.

use bsbit_core::alphabet::Base;
use bsbit_core::bisulfite::BisulfiteStrand;
use bsbit_index::reference::ReferenceIndex;
use bsbit_index::storage::fm::{ProjectedBase, SearchBase};

use crate::AlignmentError;
use crate::read_mapping::{ReadCandidate, ungapped_distance};
use crate::read_mapping_limits::{MAX_READ_BASES, MIN_SUFFIX_BASES};
use crate::search::combined_query::{CombinedSearchReferenceExt, CombinedSeedMatches};

const MINIMUM_READ_BASES: usize = 3;
const INITIAL_MINIMUM_MULTI_HIT_SEED_BASES: usize = 17;
const INITIAL_MAXIMUM_SEED_HITS: u64 = 1_000;
const INITIAL_MAXIMUM_COMBINED_RESCUE_HITS: u64 = 4_096;
const INITIAL_MAXIMUM_SEED_ROUNDS: usize = 5;
pub(crate) const DEFAULT_MINIMUM_MULTI_HIT_SEED_BASES: usize = 16;
pub(crate) const DEFAULT_MAXIMUM_SEED_HITS: u64 = 1_000;
pub(crate) const DEFAULT_MAXIMUM_COMBINED_RESCUE_HITS: u64 = 4_096;
pub(crate) const DEFAULT_MAXIMUM_SEED_ROUNDS: usize = 6;
const SENSITIVE_MAXIMUM_SEED_HITS: u64 = 4_096;
const SENSITIVE_MAXIMUM_COMBINED_RESCUE_HITS: u64 = 4_096;
const SENSITIVE_MAXIMUM_SEED_ROUNDS: usize = 6;
pub(crate) const EMPTY_SEED_STEP: usize = 8;
pub(crate) const DIRECT_SINGLETON_PROOF: u8 = 1 << 7;
pub(crate) const FLEXIBLE_NOMINAL_PROOF: u8 = 1 << 6;

#[derive(Clone, Copy)]
pub(crate) struct CombinedSearchLimits {
    pub(crate) minimum_multi_hit_seed_bases: usize,
    pub(crate) maximum_seed_hits: u64,
    pub(crate) maximum_combined_rescue_hits: u64,
    pub(crate) maximum_seed_rounds: usize,
}

pub(crate) const INITIAL_SEARCH_LIMITS: CombinedSearchLimits = CombinedSearchLimits {
    minimum_multi_hit_seed_bases: INITIAL_MINIMUM_MULTI_HIT_SEED_BASES,
    maximum_seed_hits: INITIAL_MAXIMUM_SEED_HITS,
    maximum_combined_rescue_hits: INITIAL_MAXIMUM_COMBINED_RESCUE_HITS,
    maximum_seed_rounds: INITIAL_MAXIMUM_SEED_ROUNDS,
};

pub(crate) const DEFAULT_SEARCH_LIMITS: CombinedSearchLimits = CombinedSearchLimits {
    minimum_multi_hit_seed_bases: DEFAULT_MINIMUM_MULTI_HIT_SEED_BASES,
    maximum_seed_hits: DEFAULT_MAXIMUM_SEED_HITS,
    maximum_combined_rescue_hits: DEFAULT_MAXIMUM_COMBINED_RESCUE_HITS,
    maximum_seed_rounds: DEFAULT_MAXIMUM_SEED_ROUNDS,
};

pub(crate) const SENSITIVE_SEARCH_LIMITS: CombinedSearchLimits = CombinedSearchLimits {
    minimum_multi_hit_seed_bases: DEFAULT_MINIMUM_MULTI_HIT_SEED_BASES,
    maximum_seed_hits: SENSITIVE_MAXIMUM_SEED_HITS,
    maximum_combined_rescue_hits: SENSITIVE_MAXIMUM_COMBINED_RESCUE_HITS,
    maximum_seed_rounds: SENSITIVE_MAXIMUM_SEED_ROUNDS,
};

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
    relabel_mate2: bool,
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
                let strand = if relabel_mate2 {
                    match hit.strand() {
                        BisulfiteStrand::OT => BisulfiteStrand::CTOT,
                        BisulfiteStrand::OB => BisulfiteStrand::CTOB,
                        BisulfiteStrand::CTOT | BisulfiteStrand::CTOB => return true,
                    }
                } else {
                    hit.strand()
                };
                let mut candidate = ReadCandidate {
                    contig_ordinal: hit.contig_ordinal(),
                    start: hit.start(),
                    strand,
                    proof_mask: FLEXIBLE_NOMINAL_PROOF | (1_u8 << round),
                };
                if round == 0
                    && hits == 1
                    && let Some(distance) = ungapped_distance(reference, read, candidate)
                {
                    candidate.proof_mask = DIRECT_SINGLETON_PROOF | distance;
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
    reverse_second_lane_hits: bool,
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
                let strand = if lane == 1 && reverse_second_lane_hits {
                    match hit.strand() {
                        BisulfiteStrand::OT => BisulfiteStrand::CTOT,
                        BisulfiteStrand::OB => BisulfiteStrand::CTOB,
                        BisulfiteStrand::CTOT | BisulfiteStrand::CTOB => return,
                    }
                } else {
                    hit.strand()
                };
                let mut candidate = ReadCandidate {
                    contig_ordinal: hit.contig_ordinal(),
                    start: hit.start(),
                    strand,
                    proof_mask: FLEXIBLE_NOMINAL_PROOF | (1_u8 << round),
                };
                if round == 0
                    && hits[lane] == 1
                    && let Some(distance) = ungapped_distance(reference, reads[lane], candidate)
                {
                    candidate.proof_mask = DIRECT_SINGLETON_PROOF | distance;
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
    reverse_second_lane_hits: bool,
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
        reverse_second_lane_hits,
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
    reverse_second_lane_hits: bool,
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
                reverse_second_lane_hits,
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
        consume_lane!(0, mate1_candidates, false);
        consume_lane!(1, mate2_candidates, reverse_second_lane_hits);
    }
    Ok(())
}

pub(crate) fn continue_combined_two_lane_search(
    reference: &ReferenceIndex,
    reads: [&[Base]; 2],
    reversed_projected: [&[ProjectedBase]; 2],
    reverse_second_lane_hits: bool,
    state: &mut CombinedTwoLaneSearchState,
    mate1_candidates: &mut Vec<ReadCandidate>,
    mate2_candidates: &mut Vec<ReadCandidate>,
) -> Result<[u64; 2], AlignmentError> {
    continue_combined_two_lane_search_with_limits(
        reference,
        reads,
        reversed_projected,
        reverse_second_lane_hits,
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
    reverse_second_lane_hits: bool,
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
                lane == 1 && reverse_second_lane_hits,
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
        reverse_second_lane_hits,
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
    reverse_complement_query: bool,
    output: &mut [ProjectedBase; MAX_READ_BASES],
) -> Result<(), AlignmentError> {
    if !(MINIMUM_READ_BASES..=MAX_READ_BASES).contains(&read.len()) {
        return Err(AlignmentError::UnsupportedReadLength { length: read.len() });
    }
    if reverse_complement_query {
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
    reverse_complement_query: bool,
    output: &mut [SearchBase; MAX_READ_BASES],
) -> Result<(), AlignmentError> {
    if !(MINIMUM_READ_BASES..=MAX_READ_BASES).contains(&read.len()) {
        return Err(AlignmentError::UnsupportedReadLength { length: read.len() });
    }
    if reverse_complement_query {
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
