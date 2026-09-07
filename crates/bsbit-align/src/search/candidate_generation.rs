//! Deterministic fixed-seed candidate generation.
//!
//! Public contracts and owner-bound result types remain in the parent module;
//! this module owns transient evidence, vote merging, and allocation preflight.

use core::cmp::Ordering;
use core::mem::size_of;
use core::num::NonZeroU64;

use bsbit_core::alphabet::Base;
use bsbit_core::bisulfite::{AlignmentOrientation, BisulfiteStrand, strand_semantics};
use bsbit_core::coordinate::{CoordinateError, QueryInterval, QueryLength};
use bsbit_core::sequence::NormalizedSequence;
#[cfg(test)]
use bsbit_index::reference::ContigId;
use bsbit_index::reference::{
    ProjectedMatches, ReferenceIndex, ReferenceQueryError, ReferenceQueryLimits,
};

use super::candidate::{
    CandidateAllocation, CandidateAnchor, CandidateCounter, CandidateDiagonal, CandidateError,
    CandidateInvariant, CandidateLimits, CandidateMetrics, CandidateSet, FixedSeedPlan,
    FixedSeedRequest, QueryBoundary, SeedPlanError, SeedPlanLimits,
};

struct RetainedMatches {
    request_ordinal: u64,
    request: FixedSeedRequest,
    matches: ProjectedMatches,
}

#[cfg(test)]
pub(super) struct RawEvidence {
    pub(super) contig: ContigId,
    pub(super) strand: BisulfiteStrand,
    pub(super) diagonal: CandidateDiagonal,
    pub(super) request_ordinal: u64,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct CandidateVoteKey {
    contig_ordinal: u64,
    strand: BisulfiteStrand,
    diagonal: CandidateDiagonal,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct CandidateVote {
    key: CandidateVoteKey,
    support: u64,
}

/// Generates a complete deterministic fixed-seed candidate set.
///
/// Every exact occurrence streams into a one-request lightweight key buffer.
/// Sorted request keys are merged into global `(contig, diagonal, strand)`
/// votes without constructing owner-bearing raw-evidence records. Distinct
/// requests contribute support to one final anchor. Candidate limits reject
/// the whole result and never truncate evidence.
///
/// # Errors
///
/// Returns [`CandidateError`] for Level 2B search/locate failures, aggregate
/// resource limits, allocation failures, coordinate failures, or defensive
/// invariant violations.
#[allow(clippy::too_many_lines)]
pub fn candidates_for_fixed_seeds(
    reference: &ReferenceIndex,
    plan: &FixedSeedPlan,
    query_limits: ReferenceQueryLimits,
    candidate_limits: CandidateLimits,
) -> Result<CandidateSet, CandidateError> {
    let request_count = plan.metrics.request_count;
    let retained_storage = preflight_candidate_allocation::<RetainedMatches>(
        request_count,
        CandidateAllocation::RetainedMatches,
    )?;
    let mut retained = Vec::new();
    retained
        .try_reserve_exact(retained_storage)
        .map_err(|_| CandidateError::AllocationFailed {
            allocation: CandidateAllocation::RetainedMatches,
            elements: request_count,
        })?;

    let caller_hit_limit = query_limits.max_exact_hits();
    let mut total_exact_hits = 0_u64;
    let mut matched_intervals = 0_u64;
    let mut zero_hit_requests = 0_u64;
    let mut search_rank_operations = 0_u64;

    for (storage, request) in plan.requests.iter().copied().enumerate() {
        let request_ordinal = candidate_physical_to_logical(CandidateCounter::Requests, storage)?;
        let remaining = candidate_limits
            .max_total_exact_hits
            .checked_sub(total_exact_hits)
            .ok_or(CandidateError::Invariant {
                invariant: CandidateInvariant::AggregateHitLimit,
                expected: candidate_limits.max_total_exact_hits,
                observed: total_exact_hits,
            })?;
        let effective_hit_limit = caller_hit_limit.min(remaining);
        let seed = candidate_seed_slice(plan, request, request_ordinal)?;
        let effective_limits = query_limits.with_max_exact_hits(effective_hit_limit);
        let matches = match reference.exact_search(request.strand, seed, effective_limits) {
            Ok(matches) => matches,
            Err(source @ ReferenceQueryError::HitLimitExceeded { requested, .. }) => {
                let requested_total = total_exact_hits.checked_add(requested).ok_or(
                    CandidateError::AggregateHitCountOverflow {
                        accumulated: total_exact_hits,
                        request_hits: requested,
                    },
                )?;
                if caller_hit_limit <= remaining {
                    return Err(CandidateError::Search {
                        request_ordinal,
                        strand: request.strand,
                        interval: request.interval,
                        source,
                    });
                }
                return Err(CandidateError::AggregateHitLimitExceeded {
                    accumulated: total_exact_hits,
                    request_hits: requested,
                    requested: requested_total,
                    maximum: candidate_limits.max_total_exact_hits,
                });
            }
            Err(source) => {
                return Err(CandidateError::Search {
                    request_ordinal,
                    strand: request.strand,
                    interval: request.interval,
                    source,
                });
            }
        };

        total_exact_hits = checked_candidate_add(
            CandidateCounter::TotalExactHits,
            total_exact_hits,
            matches.exact_hit_count(),
        )?;
        search_rank_operations = checked_candidate_add(
            CandidateCounter::SearchRankOperations,
            search_rank_operations,
            matches.search_rank_operations(),
        )?;
        if total_exact_hits > candidate_limits.max_total_exact_hits {
            return Err(CandidateError::Invariant {
                invariant: CandidateInvariant::AggregateHitLimit,
                expected: candidate_limits.max_total_exact_hits,
                observed: total_exact_hits,
            });
        }
        matched_intervals = checked_candidate_add(
            CandidateCounter::MatchedIntervals,
            matched_intervals,
            matches.matched_interval_count(),
        )?;
        if matches.is_empty() {
            zero_hit_requests =
                checked_candidate_add(CandidateCounter::ZeroHitRequests, zero_hit_requests, 1)?;
        }
        ensure_candidate_capacity(
            CandidateInvariant::RetainedMatchCapacity,
            request_count,
            candidate_physical_to_logical(CandidateCounter::Requests, retained.len())?,
        )?;
        retained.push(RetainedMatches {
            request_ordinal,
            request,
            matches,
        });
    }

    if matched_intervals > total_exact_hits {
        return Err(CandidateError::Invariant {
            invariant: CandidateInvariant::MatchedIntervalsWithinHits,
            expected: total_exact_hits,
            observed: matched_intervals,
        });
    }

    let mut votes = Vec::<CandidateVote>::new();
    let mut locate_calls = 0_u64;
    let mut located_coordinates = 0_u64;
    let mut locate_lf_steps = 0_u64;
    let mut locate_rank_operations = 0_u64;
    let mut locate_interval_nodes = 0_u64;
    let mut candidate_key_materializations = 0_u64;
    let mut peak_request_candidate_keys = 0_u64;
    for retained_request in &retained {
        let request_hits = retained_request.matches.exact_hit_count();
        let key_storage = preflight_candidate_allocation::<CandidateVoteKey>(
            request_hits,
            CandidateAllocation::RequestCandidateKeys,
        )?;
        let mut request_keys = Vec::new();
        request_keys.try_reserve_exact(key_storage).map_err(|_| {
            CandidateError::AllocationFailed {
                allocation: CandidateAllocation::RequestCandidateKeys,
                elements: request_hits,
            }
        })?;
        let oriented_start = oriented_seed_start(
            retained_request.request,
            plan.query_length(),
            retained_request.request_ordinal,
        )?;
        let mut visitor_error = None;
        let locate_metrics = reference
            .visit_located_matches(&retained_request.matches, &mut |hit| {
                if let Err(error) = validate_hit_semantics(
                    retained_request.request_ordinal,
                    retained_request.request,
                    hit.strand(),
                    hit.interval().len(),
                ) {
                    visitor_error = Some(error);
                    return false;
                }
                let materialized = match candidate_physical_to_logical(
                    CandidateCounter::CandidateKeyMaterializations,
                    request_keys.len(),
                ) {
                    Ok(materialized) => materialized,
                    Err(error) => {
                        visitor_error = Some(error);
                        return false;
                    }
                };
                if let Err(error) = ensure_candidate_capacity(
                    CandidateInvariant::CandidateKeyCapacity,
                    request_hits,
                    materialized,
                ) {
                    visitor_error = Some(error);
                    return false;
                }
                request_keys.push(CandidateVoteKey {
                    contig_ordinal: hit.contig().ordinal(),
                    strand: hit.strand(),
                    diagonal: CandidateDiagonal::from_difference(
                        hit.interval().start(),
                        oriented_start,
                    ),
                });
                true
            })
            .map_err(|source| CandidateError::Locate {
                request_ordinal: retained_request.request_ordinal,
                strand: retained_request.request.strand,
                interval: retained_request.request.interval,
                source,
            })?;
        if let Some(error) = visitor_error {
            return Err(error);
        }
        locate_calls = checked_candidate_add(CandidateCounter::LocateCalls, locate_calls, 1)?;
        located_coordinates = checked_candidate_add(
            CandidateCounter::LocatedCoordinates,
            located_coordinates,
            locate_metrics.located_coordinates(),
        )?;
        locate_lf_steps = checked_candidate_add(
            CandidateCounter::LocateLfSteps,
            locate_lf_steps,
            locate_metrics.lf_steps(),
        )?;
        locate_rank_operations = checked_candidate_add(
            CandidateCounter::LocateRankOperations,
            locate_rank_operations,
            locate_metrics.rank_operations(),
        )?;
        locate_interval_nodes = checked_candidate_add(
            CandidateCounter::LocateIntervalNodes,
            locate_interval_nodes,
            locate_metrics.interval_nodes(),
        )?;
        let request_key_count = candidate_physical_to_logical(
            CandidateCounter::CandidateKeyMaterializations,
            request_keys.len(),
        )?;
        candidate_key_materializations = checked_candidate_add(
            CandidateCounter::CandidateKeyMaterializations,
            candidate_key_materializations,
            request_key_count,
        )?;
        peak_request_candidate_keys = peak_request_candidate_keys.max(request_key_count);
        if request_key_count != request_hits {
            return Err(CandidateError::Invariant {
                invariant: CandidateInvariant::LocatedHitCount,
                expected: request_hits,
                observed: request_key_count,
            });
        }
        request_keys.sort_unstable_by(compare_candidate_vote_keys);
        if let Some(key) = request_keys.windows(2).find_map(|pair| {
            (compare_candidate_vote_keys(&pair[0], &pair[1]) == Ordering::Equal).then_some(pair[0])
        }) {
            return Err(CandidateError::DuplicateRequestEvidence {
                request_ordinal: retained_request.request_ordinal,
                contig_ordinal: key.contig_ordinal,
                strand: key.strand,
                diagonal: key.diagonal,
            });
        }
        merge_candidate_votes(&mut votes, &request_keys)?;
    }

    if located_coordinates != total_exact_hits {
        return Err(CandidateError::Invariant {
            invariant: CandidateInvariant::LocatedHitCount,
            expected: total_exact_hits,
            observed: located_coordinates,
        });
    }
    if candidate_key_materializations != total_exact_hits {
        return Err(CandidateError::Invariant {
            invariant: CandidateInvariant::LocatedHitCount,
            expected: total_exact_hits,
            observed: candidate_key_materializations,
        });
    }

    let unique_candidates =
        candidate_physical_to_logical(CandidateCounter::UniqueCandidates, votes.len())?;
    if unique_candidates > candidate_limits.max_unique_candidates {
        return Err(CandidateError::UniqueCandidateLimitExceeded {
            requested: unique_candidates,
            maximum: candidate_limits.max_unique_candidates,
        });
    }

    let anchor_storage = preflight_candidate_allocation::<CandidateAnchor>(
        unique_candidates,
        CandidateAllocation::FinalAnchors,
    )?;
    let mut anchors = Vec::new();
    anchors
        .try_reserve_exact(anchor_storage)
        .map_err(|_| CandidateError::AllocationFailed {
            allocation: CandidateAllocation::FinalAnchors,
            elements: unique_candidates,
        })?;

    let mut support_sum = 0_u64;
    let mut maximum_support = 0_u64;
    for vote in votes {
        let support = NonZeroU64::new(vote.support).ok_or(CandidateError::Invariant {
            invariant: CandidateInvariant::SupportSum,
            expected: 1,
            observed: 0,
        })?;
        support_sum =
            checked_candidate_add(CandidateCounter::SupportSum, support_sum, support.get())?;
        maximum_support = maximum_support.max(support.get());
        ensure_candidate_capacity(
            CandidateInvariant::FinalAnchorCapacity,
            unique_candidates,
            candidate_physical_to_logical(CandidateCounter::UniqueCandidates, anchors.len())?,
        )?;
        let contig = reference.contig_id(vote.key.contig_ordinal).map_err(|_| {
            CandidateError::Invariant {
                invariant: CandidateInvariant::LocatedContigOrdinal,
                expected: reference.contig_count(),
                observed: vote.key.contig_ordinal,
            }
        })?;
        anchors.push(CandidateAnchor {
            contig,
            strand: vote.key.strand,
            diagonal: vote.key.diagonal,
            support,
        });
    }

    let final_count =
        candidate_physical_to_logical(CandidateCounter::UniqueCandidates, anchors.len())?;
    let output_ordered = anchors
        .windows(2)
        .all(|pair| compare_anchors(&pair[0], &pair[1]).is_lt());
    let duplicate_evidence = validate_final_candidate_invariants(
        total_exact_hits,
        unique_candidates,
        support_sum,
        final_count,
        output_ordered,
    )?;

    Ok(CandidateSet {
        reference: reference.instance_id(),
        query: plan.query_instance_id(),
        anchors,
        metrics: CandidateMetrics {
            request_count,
            total_seed_bases: plan.metrics.total_seed_bases,
            total_exact_hits,
            matched_intervals,
            unique_candidates,
            duplicate_evidence,
            maximum_support,
            zero_hit_requests,
            search_rank_operations,
            locate_calls,
            located_coordinates,
            locate_lf_steps,
            locate_rank_operations,
            locate_interval_nodes,
            candidate_key_materializations,
            peak_request_candidate_keys,
        },
    })
}

pub(super) fn validate_supplied_requests(
    query: &NormalizedSequence,
    supplied: &[FixedSeedRequest],
    query_length: QueryLength,
    limits: SeedPlanLimits,
) -> Result<u64, SeedPlanError> {
    let mut total_seed_bases = 0_u64;
    for (storage, request) in supplied.iter().copied().enumerate() {
        let request_ordinal = request_count_to_u64(storage)?;
        let interval = QueryInterval::new(
            request.interval.start(),
            request.interval.end(),
            query_length,
        )
        .map_err(|source| SeedPlanError::InvalidInterval {
            request_ordinal,
            source,
        })?;
        if interval.is_empty() {
            return Err(SeedPlanError::EmptySeed {
                request_ordinal,
                interval,
            });
        }
        let seed = query_slice(query, interval, request_ordinal)?;
        if let Some(local) = seed.iter().position(|base| *base == Base::N) {
            let local =
                u64::try_from(local).map_err(|_| SeedPlanError::SeedOffsetNotRepresentable {
                    request_ordinal,
                    value: local,
                })?;
            let query_offset =
                interval
                    .start()
                    .checked_add(local)
                    .ok_or(SeedPlanError::QueryOffsetOverflow {
                        request_ordinal,
                        start: interval.start(),
                        local_offset: local,
                    })?;
            return Err(SeedPlanError::UnsearchableBase {
                request_ordinal,
                query_offset,
            });
        }
        total_seed_bases = total_seed_bases.checked_add(interval.len()).ok_or(
            SeedPlanError::TotalSeedBasesOverflow {
                accumulated: total_seed_bases,
                next: interval.len(),
            },
        )?;
        if total_seed_bases > limits.max_total_seed_bases {
            return Err(SeedPlanError::TotalSeedBasesLimitExceeded {
                request_ordinal,
                requested: total_seed_bases,
                maximum: limits.max_total_seed_bases,
            });
        }
    }
    Ok(total_seed_bases)
}

fn query_slice(
    query: &NormalizedSequence,
    interval: QueryInterval,
    request_ordinal: u64,
) -> Result<&[Base], SeedPlanError> {
    let start =
        usize::try_from(interval.start()).map_err(|_| SeedPlanError::BoundaryNotRepresentable {
            request_ordinal,
            boundary: QueryBoundary::Start,
            value: interval.start(),
        })?;
    let end =
        usize::try_from(interval.end()).map_err(|_| SeedPlanError::BoundaryNotRepresentable {
            request_ordinal,
            boundary: QueryBoundary::End,
            value: interval.end(),
        })?;
    query
        .bases()
        .get(start..end)
        .ok_or(SeedPlanError::InvalidInterval {
            request_ordinal,
            source: CoordinateError::OutOfBounds {
                domain: bsbit_core::coordinate::CoordinateDomain::Query,
                operation: bsbit_core::coordinate::CoordinateOperation::IntervalConstruction,
                start: interval.start(),
                end: interval.end(),
                length: query.len(),
            },
        })
}

fn candidate_seed_slice(
    plan: &FixedSeedPlan,
    request: FixedSeedRequest,
    request_ordinal: u64,
) -> Result<&[Base], CandidateError> {
    let Ok(start) = usize::try_from(request.interval.start()) else {
        return Err(CandidateError::PlanIntervalStorage {
            request_ordinal,
            start: request.interval.start(),
            end: request.interval.end(),
            query_bases: plan.metrics.query_bases,
        });
    };
    let Ok(end) = usize::try_from(request.interval.end()) else {
        return Err(CandidateError::PlanIntervalStorage {
            request_ordinal,
            start: request.interval.start(),
            end: request.interval.end(),
            query_bases: plan.metrics.query_bases,
        });
    };
    plan.query()
        .bases()
        .get(start..end)
        .ok_or(CandidateError::PlanIntervalStorage {
            request_ordinal,
            start: request.interval.start(),
            end: request.interval.end(),
            query_bases: plan.metrics.query_bases,
        })
}

fn oriented_seed_start(
    request: FixedSeedRequest,
    query_length: QueryLength,
    request_ordinal: u64,
) -> Result<u64, CandidateError> {
    match strand_semantics(request.strand).orientation() {
        AlignmentOrientation::Forward => Ok(request.interval.start()),
        AlignmentOrientation::Reverse => request
            .interval
            .reverse(query_length)
            .map(QueryInterval::start)
            .map_err(|source| CandidateError::OrientedInterval {
                request_ordinal,
                source,
            }),
    }
}

pub(super) fn validate_final_candidate_invariants(
    total_exact_hits: u64,
    unique_candidates: u64,
    support_sum: u64,
    final_count: u64,
    output_ordered: bool,
) -> Result<u64, CandidateError> {
    if support_sum != total_exact_hits {
        return Err(CandidateError::Invariant {
            invariant: CandidateInvariant::SupportSum,
            expected: total_exact_hits,
            observed: support_sum,
        });
    }
    if final_count != unique_candidates {
        return Err(CandidateError::Invariant {
            invariant: CandidateInvariant::CandidateCount,
            expected: unique_candidates,
            observed: final_count,
        });
    }
    if !output_ordered {
        return Err(CandidateError::Invariant {
            invariant: CandidateInvariant::OutputOrder,
            expected: 1,
            observed: 0,
        });
    }
    total_exact_hits
        .checked_sub(unique_candidates)
        .ok_or(CandidateError::Invariant {
            invariant: CandidateInvariant::DuplicateEvidence,
            expected: total_exact_hits,
            observed: unique_candidates,
        })
}

pub(super) fn validate_hit_semantics(
    request_ordinal: u64,
    request: FixedSeedRequest,
    observed_strand: BisulfiteStrand,
    observed_length: u64,
) -> Result<(), CandidateError> {
    if observed_strand != request.strand {
        return Err(CandidateError::HitStrandMismatch {
            request_ordinal,
            expected: request.strand,
            observed: observed_strand,
        });
    }
    if observed_length != request.interval.len() {
        return Err(CandidateError::HitLengthMismatch {
            request_ordinal,
            expected: request.interval.len(),
            observed: observed_length,
        });
    }
    Ok(())
}

fn compare_candidate_vote_keys(lhs: &CandidateVoteKey, rhs: &CandidateVoteKey) -> Ordering {
    lhs.contig_ordinal
        .cmp(&rhs.contig_ordinal)
        .then_with(|| lhs.diagonal.cmp(&rhs.diagonal))
        .then_with(|| strand_rank(lhs.strand).cmp(&strand_rank(rhs.strand)))
}

fn merge_candidate_votes(
    votes: &mut Vec<CandidateVote>,
    request_keys: &[CandidateVoteKey],
) -> Result<(), CandidateError> {
    if request_keys.is_empty() {
        return Ok(());
    }
    let old_count = candidate_physical_to_logical(CandidateCounter::UniqueCandidates, votes.len())?;
    let request_count = candidate_physical_to_logical(
        CandidateCounter::CandidateKeyMaterializations,
        request_keys.len(),
    )?;
    let expanded_count =
        checked_candidate_add(CandidateCounter::UniqueCandidates, old_count, request_count)?;
    let expanded_storage = preflight_candidate_allocation::<CandidateVote>(
        expanded_count,
        CandidateAllocation::CandidateVotes,
    )?;
    let additional =
        expanded_storage
            .checked_sub(votes.len())
            .ok_or(CandidateError::Invariant {
                invariant: CandidateInvariant::CandidateCount,
                expected: old_count,
                observed: expanded_count,
            })?;
    votes
        .try_reserve_exact(additional)
        .map_err(|_| CandidateError::AllocationFailed {
            allocation: CandidateAllocation::CandidateVotes,
            elements: expanded_count,
        })?;

    let placeholder = CandidateVote {
        key: request_keys[0],
        support: 1,
    };
    votes.resize(expanded_storage, placeholder);
    let mut left = usize::try_from(old_count).expect("existing vote count fits usize");
    let mut right = request_keys.len();
    let mut write = expanded_storage;
    while left != 0 || right != 0 {
        let ordering = match (left.checked_sub(1), right.checked_sub(1)) {
            (Some(left_index), Some(right_index)) => {
                compare_candidate_vote_keys(&votes[left_index].key, &request_keys[right_index])
            }
            (Some(_), None) => Ordering::Greater,
            (None, Some(_)) => Ordering::Less,
            (None, None) => break,
        };
        write -= 1;
        match ordering {
            Ordering::Greater => {
                left -= 1;
                votes[write] = votes[left];
            }
            Ordering::Less => {
                right -= 1;
                votes[write] = CandidateVote {
                    key: request_keys[right],
                    support: 1,
                };
            }
            Ordering::Equal => {
                left -= 1;
                right -= 1;
                votes[write] = CandidateVote {
                    key: votes[left].key,
                    support: checked_candidate_add(
                        CandidateCounter::Support,
                        votes[left].support,
                        1,
                    )?,
                };
            }
        }
    }
    votes.copy_within(write..expanded_storage, 0);
    votes.truncate(expanded_storage - write);
    Ok(())
}

#[cfg(test)]
pub(super) fn count_unique_candidates(raw: &[RawEvidence]) -> Result<u64, CandidateError> {
    let mut unique = 0_u64;
    let mut prior: Option<&RawEvidence> = None;
    for evidence in raw {
        if let Some(previous) = prior {
            if same_candidate(previous, evidence)
                && previous.request_ordinal == evidence.request_ordinal
            {
                return Err(CandidateError::DuplicateRequestEvidence {
                    request_ordinal: evidence.request_ordinal,
                    contig_ordinal: evidence.contig.ordinal(),
                    strand: evidence.strand,
                    diagonal: evidence.diagonal,
                });
            }
            if !same_candidate(previous, evidence) {
                unique = checked_candidate_add(CandidateCounter::UniqueCandidates, unique, 1)?;
            }
        } else {
            unique = 1;
        }
        prior = Some(evidence);
    }
    Ok(unique)
}

#[cfg(test)]
fn same_candidate(lhs: &RawEvidence, rhs: &RawEvidence) -> bool {
    lhs.contig.ordinal() == rhs.contig.ordinal()
        && lhs.diagonal == rhs.diagonal
        && lhs.strand == rhs.strand
}

fn compare_anchors(lhs: &CandidateAnchor, rhs: &CandidateAnchor) -> Ordering {
    lhs.contig
        .ordinal()
        .cmp(&rhs.contig.ordinal())
        .then_with(|| lhs.diagonal.cmp(&rhs.diagonal))
        .then_with(|| strand_rank(lhs.strand).cmp(&strand_rank(rhs.strand)))
}

pub(super) const fn strand_rank(strand: BisulfiteStrand) -> u8 {
    match strand {
        BisulfiteStrand::OT => 0,
        BisulfiteStrand::OB => 1,
        BisulfiteStrand::CTOT => 2,
        BisulfiteStrand::CTOB => 3,
    }
}

pub(super) fn request_count_to_u64(value: usize) -> Result<u64, SeedPlanError> {
    u64::try_from(value).map_err(|_| SeedPlanError::RequestCountNotRepresentable { value })
}

fn candidate_physical_to_logical(
    counter: CandidateCounter,
    value: usize,
) -> Result<u64, CandidateError> {
    u64::try_from(value).map_err(|_| CandidateError::CountNotRepresentable { counter, value })
}

pub(super) fn checked_candidate_add(
    counter: CandidateCounter,
    accumulated: u64,
    next: u64,
) -> Result<u64, CandidateError> {
    accumulated
        .checked_add(next)
        .ok_or(CandidateError::CounterOverflow {
            counter,
            accumulated,
            next,
        })
}

pub(super) fn ensure_candidate_capacity(
    invariant: CandidateInvariant,
    reserved: u64,
    materialized: u64,
) -> Result<(), CandidateError> {
    if materialized >= reserved {
        Err(CandidateError::Invariant {
            invariant,
            expected: reserved,
            observed: materialized,
        })
    } else {
        Ok(())
    }
}

pub(super) fn preflight_seed_allocation<T>(
    elements: u64,
    allocation: CandidateAllocation,
) -> Result<usize, SeedPlanError> {
    preflight_storage::<T>(elements).map_err(|(elements, element_size)| {
        SeedPlanError::AllocationSizeOverflow {
            allocation,
            elements,
            element_size,
        }
    })
}

pub(super) fn preflight_candidate_allocation<T>(
    elements: u64,
    allocation: CandidateAllocation,
) -> Result<usize, CandidateError> {
    preflight_storage::<T>(elements).map_err(|(elements, element_size)| {
        CandidateError::AllocationSizeOverflow {
            allocation,
            elements,
            element_size,
        }
    })
}

fn preflight_storage<T>(elements: u64) -> Result<usize, (u64, u64)> {
    let element_size = u64::try_from(size_of::<T>()).map_err(|_| (elements, u64::MAX))?;
    elements
        .checked_mul(element_size)
        .ok_or((elements, element_size))?;
    let storage = usize::try_from(elements).map_err(|_| (elements, element_size))?;
    if size_of::<T>() != 0 && storage > isize::MAX.unsigned_abs() / size_of::<T>() {
        return Err((elements, element_size));
    }
    Ok(storage)
}
