//! Deterministic grouping of fixed-seed reference evidence into candidates.
//!
//! Candidate anchors are pre-extension evidence, not verified alignments,
//! placements, ambiguity decisions, or MAPQ values.
//!
//! Anchors cannot be cloned away from their owner-bound candidate set:
//!
//! ```compile_fail
//! use bsbit_align::search::candidate::CandidateAnchor;
//!
//! fn requires_clone<T: Clone>() {}
//! requires_clone::<CandidateAnchor>();
//! ```
//!
//! Candidate sets use shared borrowing or a caller-owned `Arc`; ordinary
//! cloning would make their variable-sized allocation infallible:
//!
//! ```compile_fail
//! use bsbit_align::search::candidate::CandidateSet;
//!
//! fn requires_clone<T: Clone>() {}
//! requires_clone::<CandidateSet>();
//! ```

use core::cmp::Ordering;
use core::fmt;
use core::hash::{Hash, Hasher};
use core::num::NonZeroU64;

#[cfg(test)]
use super::fixed_seed::preflight_seed_allocation;
pub use super::fixed_seed::{
    CandidateAllocation, FixedSeedPlan, FixedSeedRequest, QueryBoundary, QueryInstanceId,
    SeedPlanError, SeedPlanLimits, SeedPlanMetrics,
};
use super::fixed_seed::{preflight_storage, strand_rank};
use bsbit_core::alphabet::Base;
use bsbit_core::bisulfite::{AlignmentOrientation, BisulfiteStrand, strand_semantics};
use bsbit_core::coordinate::{CoordinateError, CoordinateShift, QueryInterval, QueryLength};
use bsbit_core::sequence::NormalizedSequence;
use bsbit_index::reference::{
    ContigId, ProjectedMatches, ReferenceAccessError, ReferenceIndex, ReferenceInstanceId,
    ReferenceLocateError, ReferenceQueryError, ReferenceQueryLimits,
};

/// A full-range signed candidate diagonal relative to contig coordinate zero.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct CandidateDiagonal(CoordinateShift);

impl CandidateDiagonal {
    /// Creates a negative diagonal with a nonzero magnitude.
    #[must_use]
    pub const fn before_contig(magnitude: NonZeroU64) -> Self {
        Self(CoordinateShift::Backward(magnitude))
    }

    /// Creates zero or a positive contig-relative diagonal.
    #[must_use]
    pub const fn at_or_after_contig(value: u64) -> Self {
        Self(CoordinateShift::forward(value))
    }

    /// Forms the exact mathematical difference `reference_start - query_start`.
    #[must_use]
    pub const fn from_difference(reference_start: u64, query_start: u64) -> Self {
        if reference_start >= query_start {
            Self::at_or_after_contig(reference_start - query_start)
        } else {
            match NonZeroU64::new(query_start - reference_start) {
                Some(magnitude) => Self::before_contig(magnitude),
                None => Self::at_or_after_contig(0),
            }
        }
    }

    /// Returns the underlying accepted signed-magnitude primitive.
    #[must_use]
    pub const fn shift(self) -> CoordinateShift {
        self.0
    }

    /// Returns whether the diagonal lies before contig coordinate zero.
    #[must_use]
    pub const fn is_before_contig(self) -> bool {
        matches!(self.0, CoordinateShift::Backward(_))
    }

    /// Returns the unsigned magnitude.
    #[must_use]
    pub const fn magnitude(self) -> u64 {
        self.0.magnitude()
    }
}

impl Hash for CandidateDiagonal {
    fn hash<H: Hasher>(&self, state: &mut H) {
        match self.0 {
            CoordinateShift::Backward(value) => {
                0_u8.hash(state);
                value.get().hash(state);
            }
            CoordinateShift::Zero => {
                1_u8.hash(state);
                0_u64.hash(state);
            }
            CoordinateShift::Forward(value) => {
                2_u8.hash(state);
                value.get().hash(state);
            }
        }
    }
}

impl Ord for CandidateDiagonal {
    fn cmp(&self, other: &Self) -> Ordering {
        match (self.0, other.0) {
            (CoordinateShift::Backward(lhs), CoordinateShift::Backward(rhs)) => {
                rhs.get().cmp(&lhs.get())
            }
            (CoordinateShift::Backward(_), _)
            | (CoordinateShift::Zero, CoordinateShift::Forward(_)) => Ordering::Less,
            (_, CoordinateShift::Backward(_))
            | (CoordinateShift::Forward(_), CoordinateShift::Zero) => Ordering::Greater,
            (CoordinateShift::Zero, CoordinateShift::Zero) => Ordering::Equal,
            (CoordinateShift::Forward(lhs), CoordinateShift::Forward(rhs)) => {
                lhs.get().cmp(&rhs.get())
            }
        }
    }
}

impl PartialOrd for CandidateDiagonal {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl fmt::Display for CandidateDiagonal {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.fmt(formatter)
    }
}

/// Complete limits for one fixed-seed candidate generation call.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct CandidateLimits {
    max_total_exact_hits: u64,
    max_unique_candidates: u64,
}

impl CandidateLimits {
    /// Limits admitting every representable complete result.
    pub const MAX: Self = Self {
        max_total_exact_hits: u64::MAX,
        max_unique_candidates: u64::MAX,
    };

    /// Creates aggregate exact-hit and unique-candidate limits.
    #[must_use]
    pub const fn new(max_total_exact_hits: u64, max_unique_candidates: u64) -> Self {
        Self {
            max_total_exact_hits,
            max_unique_candidates,
        }
    }

    /// Returns the aggregate exact-hit limit.
    #[must_use]
    pub const fn max_total_exact_hits(self) -> u64 {
        self.max_total_exact_hits
    }

    /// Returns the unique-candidate limit.
    #[must_use]
    pub const fn max_unique_candidates(self) -> u64 {
        self.max_unique_candidates
    }
}

impl Default for CandidateLimits {
    fn default() -> Self {
        Self::MAX
    }
}

/// A logical candidate counter used in structured diagnostics.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum CandidateCounter {
    /// Canonical request ordinal or materialized request count.
    Requests,
    /// Total exact occurrences.
    TotalExactHits,
    /// Sum of nonempty Level 2B match intervals.
    MatchedIntervals,
    /// Number of zero-hit requests.
    ZeroHitRequests,
    /// Number of unique candidate keys.
    UniqueCandidates,
    /// One candidate support value.
    Support,
    /// Sum of all final support values.
    SupportSum,
    /// Rank-boundary operations performed by exact searches.
    SearchRankOperations,
    /// Completed locate API calls.
    LocateCalls,
    /// Coordinates streamed by locate.
    LocatedCoordinates,
    /// Logical LF steps represented by locate traversal.
    LocateLfSteps,
    /// Physical locate rank-boundary operations.
    LocateRankOperations,
    /// Shared locate interval-tree nodes.
    LocateIntervalNodes,
    /// Lightweight candidate keys materialized from locate.
    CandidateKeyMaterializations,
    /// Candidate seed starts tested inside constrained reference windows.
    RegionalSeedStarts,
    /// Base comparisons performed by constrained-window seed scanning.
    RegionalBaseComparisons,
    /// Constrained windows containing at least one hit for one request.
    RegionalMatchedWindows,
}

/// A defensive private candidate invariant.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum CandidateInvariant {
    /// Retained match materialization exceeded its reservation.
    RetainedMatchCapacity,
    /// Aggregate hits exceeded the configured maximum after a successful search.
    AggregateHitLimit,
    /// Matched intervals exceeded exact hits.
    MatchedIntervalsWithinHits,
    /// Raw evidence materialization exceeded its reservation.
    RawEvidenceCapacity,
    /// Per-request candidate-key materialization exceeded its reservation.
    CandidateKeyCapacity,
    /// Locate returned a different hit count from exact search.
    LocatedHitCount,
    /// Final anchor materialization exceeded its reservation.
    FinalAnchorCapacity,
    /// Final support did not sum to exact raw evidence.
    SupportSum,
    /// Final anchor count differed from the unique-key count.
    CandidateCount,
    /// Duplicate evidence arithmetic disagreed.
    DuplicateEvidence,
    /// Final anchors were not in normative total order.
    OutputOrder,
    /// A streamed hit referenced a missing contig ordinal.
    LocatedContigOrdinal,
}

/// A complete candidate generation failure.
#[non_exhaustive]
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum CandidateError {
    /// A candidate-owned allocation cannot fit this architecture.
    AllocationSizeOverflow {
        /// Allocation site.
        allocation: CandidateAllocation,
        /// Requested elements.
        elements: u64,
        /// Element width.
        element_size: u64,
    },
    /// A fallible candidate-owned reservation failed.
    AllocationFailed {
        /// Allocation site.
        allocation: CandidateAllocation,
        /// Requested elements.
        elements: u64,
    },
    /// A physical count cannot fit the logical width.
    CountNotRepresentable {
        /// Counter being converted.
        counter: CandidateCounter,
        /// Physical value.
        value: usize,
    },
    /// A logical counter overflowed.
    CounterOverflow {
        /// Counter being accumulated.
        counter: CandidateCounter,
        /// Accumulated value.
        accumulated: u64,
        /// Next increment.
        next: u64,
    },
    /// Level 2B exact search failed for one canonical request.
    Search {
        /// Canonical request ordinal.
        request_ordinal: u64,
        /// Request strand.
        strand: BisulfiteStrand,
        /// Raw query interval.
        interval: QueryInterval,
        /// Underlying complete-search failure.
        source: ReferenceQueryError,
    },
    /// Aggregate exact-hit prefix arithmetic overflowed.
    AggregateHitCountOverflow {
        /// Prior complete request hits.
        accumulated: u64,
        /// Current request's complete hit count.
        request_hits: u64,
    },
    /// The first exact prefix exceeds the aggregate hit limit.
    AggregateHitLimitExceeded {
        /// Prior complete request hits.
        accumulated: u64,
        /// Current request's complete hit count.
        request_hits: u64,
        /// Exact first exceeding prefix total.
        requested: u64,
        /// Configured aggregate maximum.
        maximum: u64,
    },
    /// Level 2B locate failed for one canonical request.
    Locate {
        /// Canonical request ordinal.
        request_ordinal: u64,
        /// Request strand.
        strand: BisulfiteStrand,
        /// Raw query interval.
        interval: QueryInterval,
        /// Underlying locate failure.
        source: ReferenceLocateError,
    },
    /// A constrained candidate window could not resolve its owner-bound contig.
    RegionalReferenceAccess {
        /// Zero-based constrained-window ordinal.
        window_ordinal: u64,
        /// Underlying owner/access failure.
        source: ReferenceAccessError,
    },
    /// A located hit has the wrong strand.
    HitStrandMismatch {
        /// Canonical request ordinal.
        request_ordinal: u64,
        /// Expected request strand.
        expected: BisulfiteStrand,
        /// Observed hit strand.
        observed: BisulfiteStrand,
    },
    /// A located hit has the wrong interval length.
    HitLengthMismatch {
        /// Canonical request ordinal.
        request_ordinal: u64,
        /// Expected seed length.
        expected: u64,
        /// Observed hit length.
        observed: u64,
    },
    /// Reverse-oriented seed coordinate recovery failed.
    OrientedInterval {
        /// Canonical request ordinal.
        request_ordinal: u64,
        /// Underlying coordinate failure.
        source: CoordinateError,
    },
    /// A valid plan interval was not physically addressable.
    PlanIntervalStorage {
        /// Canonical request ordinal.
        request_ordinal: u64,
        /// Raw query start.
        start: u64,
        /// Raw query end.
        end: u64,
        /// Actual query bases.
        query_bases: u64,
    },
    /// One request produced duplicate evidence for one candidate key.
    DuplicateRequestEvidence {
        /// Canonical request ordinal.
        request_ordinal: u64,
        /// Owner-bound contig ordinal.
        contig_ordinal: u64,
        /// Bisulfite strand.
        strand: BisulfiteStrand,
        /// Signed diagonal.
        diagonal: CandidateDiagonal,
    },
    /// Unique candidate keys exceed the configured complete-result limit.
    UniqueCandidateLimitExceeded {
        /// Complete unique-key count.
        requested: u64,
        /// Configured maximum.
        maximum: u64,
    },
    /// A defensive candidate invariant failed.
    Invariant {
        /// Invariant category.
        invariant: CandidateInvariant,
        /// Expected value.
        expected: u64,
        /// Observed value.
        observed: u64,
    },
}

impl fmt::Display for CandidateError {
    #[allow(clippy::too_many_lines)]
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::AllocationSizeOverflow {
                allocation,
                elements,
                element_size,
            } => write!(
                formatter,
                "cannot size {allocation:?}: {elements} elements of {element_size} bytes"
            ),
            Self::AllocationFailed {
                allocation,
                elements,
            } => write!(
                formatter,
                "failed to reserve {elements} elements for {allocation:?}"
            ),
            Self::CountNotRepresentable { counter, value } => {
                write!(
                    formatter,
                    "{counter:?} physical count {value} is not representable as u64"
                )
            }
            Self::CounterOverflow {
                counter,
                accumulated,
                next,
            } => write!(
                formatter,
                "{counter:?} count {accumulated} plus {next} overflowed"
            ),
            Self::Search {
                request_ordinal,
                strand,
                interval,
                source,
            } => write!(
                formatter,
                "candidate search failed for request {request_ordinal} {strand:?} {interval}: {source}"
            ),
            Self::AggregateHitCountOverflow {
                accumulated,
                request_hits,
            } => write!(
                formatter,
                "aggregate exact hits {accumulated} plus request count {request_hits} overflowed"
            ),
            Self::AggregateHitLimitExceeded {
                accumulated,
                request_hits,
                requested,
                maximum,
            } => write!(
                formatter,
                "aggregate exact hits {accumulated} plus request count {request_hits} is {requested}, exceeding {maximum}"
            ),
            Self::Locate {
                request_ordinal,
                strand,
                interval,
                source,
            } => write!(
                formatter,
                "candidate locate failed for request {request_ordinal} {strand:?} {interval}: {source}"
            ),
            Self::RegionalReferenceAccess {
                window_ordinal,
                source,
            } => write!(
                formatter,
                "regional candidate window {window_ordinal} could not resolve its contig: {source}"
            ),
            Self::HitStrandMismatch {
                request_ordinal,
                expected,
                observed,
            } => write!(
                formatter,
                "request {request_ordinal} expected hit strand {expected:?}, observed {observed:?}"
            ),
            Self::HitLengthMismatch {
                request_ordinal,
                expected,
                observed,
            } => write!(
                formatter,
                "request {request_ordinal} expected hit length {expected}, observed {observed}"
            ),
            Self::OrientedInterval {
                request_ordinal,
                source,
            } => write!(
                formatter,
                "request {request_ordinal} oriented interval recovery failed: {source}"
            ),
            Self::PlanIntervalStorage {
                request_ordinal,
                start,
                end,
                query_bases,
            } => write!(
                formatter,
                "request {request_ordinal} interval [{start},{end}) is not physically addressable in query length {query_bases}"
            ),
            Self::DuplicateRequestEvidence {
                request_ordinal,
                contig_ordinal,
                strand,
                diagonal,
            } => write!(
                formatter,
                "request {request_ordinal} produced duplicate evidence for contig {contig_ordinal}, strand {strand:?}, diagonal {diagonal}"
            ),
            Self::UniqueCandidateLimitExceeded { requested, maximum } => write!(
                formatter,
                "unique candidate count {requested} exceeds configured maximum {maximum}"
            ),
            Self::Invariant {
                invariant,
                expected,
                observed,
            } => write!(
                formatter,
                "{invariant:?} invariant expected {expected}, observed {observed}"
            ),
        }
    }
}

impl std::error::Error for CandidateError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Search { source, .. } => Some(source),
            Self::Locate { source, .. } => Some(source),
            Self::RegionalReferenceAccess { source, .. } => Some(source),
            Self::OrientedInterval { source, .. } => Some(source),
            _ => None,
        }
    }
}

/// Complete deterministic metrics for one candidate set.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct CandidateMetrics {
    request_count: u64,
    total_seed_bases: u64,
    total_exact_hits: u64,
    matched_intervals: u64,
    unique_candidates: u64,
    duplicate_evidence: u64,
    maximum_support: u64,
    zero_hit_requests: u64,
    search_rank_operations: u64,
    locate_calls: u64,
    located_coordinates: u64,
    locate_lf_steps: u64,
    locate_rank_operations: u64,
    locate_interval_nodes: u64,
    candidate_key_materializations: u64,
    peak_request_candidate_keys: u64,
}

impl CandidateMetrics {
    /// Returns the number of fixed seed requests.
    #[must_use]
    pub const fn request_count(self) -> u64 {
        self.request_count
    }

    /// Returns aggregate seed bases, counting overlap per request.
    #[must_use]
    pub const fn total_seed_bases(self) -> u64 {
        self.total_seed_bases
    }

    /// Returns the complete exact occurrence count.
    #[must_use]
    pub const fn total_exact_hits(self) -> u64 {
        self.total_exact_hits
    }

    /// Returns summed nonempty Level 2B match intervals.
    #[must_use]
    pub const fn matched_intervals(self) -> u64 {
        self.matched_intervals
    }

    /// Returns the number of unique pre-extension candidate keys.
    #[must_use]
    pub const fn unique_candidates(self) -> u64 {
        self.unique_candidates
    }

    /// Returns `total_exact_hits - unique_candidates`.
    #[must_use]
    pub const fn duplicate_evidence(self) -> u64 {
        self.duplicate_evidence
    }

    /// Returns maximum distinct-request support, or zero for an empty set.
    #[must_use]
    pub const fn maximum_support(self) -> u64 {
        self.maximum_support
    }

    /// Returns the number of accepted requests with zero complete hits.
    #[must_use]
    pub const fn zero_hit_requests(self) -> u64 {
        self.zero_hit_requests
    }

    /// Returns rank-boundary operations performed by exact searches.
    #[must_use]
    pub const fn search_rank_operations(self) -> u64 {
        self.search_rank_operations
    }

    /// Returns the number of complete locate calls.
    #[must_use]
    pub const fn locate_calls(self) -> u64 {
        self.locate_calls
    }

    /// Returns coordinates streamed through the candidate visitor seam.
    #[must_use]
    pub const fn located_coordinates(self) -> u64 {
        self.located_coordinates
    }

    /// Returns logical LF transitions represented by locate traversal.
    #[must_use]
    pub const fn locate_lf_steps(self) -> u64 {
        self.locate_lf_steps
    }

    /// Returns physical rank-boundary operations performed by locate.
    #[must_use]
    pub const fn locate_rank_operations(self) -> u64 {
        self.locate_rank_operations
    }

    /// Returns shared locate interval-tree nodes processed.
    #[must_use]
    pub const fn locate_interval_nodes(self) -> u64 {
        self.locate_interval_nodes
    }

    /// Returns lightweight candidate keys materialized before merging.
    #[must_use]
    pub const fn candidate_key_materializations(self) -> u64 {
        self.candidate_key_materializations
    }

    /// Returns the largest one-request candidate-key buffer.
    #[must_use]
    pub const fn peak_request_candidate_keys(self) -> u64 {
        self.peak_request_candidate_keys
    }
}

/// One borrowed pre-extension candidate and its distinct-request support.
///
/// Values cannot be constructed or cloned outside this module. Consume anchors
/// only while borrowing their owner-bound [`CandidateSet`].
pub struct CandidateAnchor {
    contig: ContigId,
    strand: BisulfiteStrand,
    diagonal: CandidateDiagonal,
    support: NonZeroU64,
}

impl CandidateAnchor {
    /// Returns the owner-bound contig identifier.
    #[must_use]
    pub const fn contig(&self) -> &ContigId {
        &self.contig
    }

    /// Returns the candidate bisulfite strand.
    #[must_use]
    pub const fn strand(&self) -> BisulfiteStrand {
        self.strand
    }

    /// Returns the signed contig-relative diagonal.
    #[must_use]
    pub const fn diagonal(&self) -> CandidateDiagonal {
        self.diagonal
    }

    /// Returns the nonzero number of distinct supporting requests.
    #[must_use]
    pub const fn support(&self) -> NonZeroU64 {
        self.support
    }
}

impl fmt::Debug for CandidateAnchor {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("CandidateAnchor")
            .field("contig_ordinal", &self.contig.ordinal())
            .field("strand", &self.strand)
            .field("diagonal", &self.diagonal)
            .field("support", &self.support)
            .finish()
    }
}

/// A complete immutable owner-bound candidate result.
pub struct CandidateSet {
    reference: ReferenceInstanceId,
    query: QueryInstanceId,
    anchors: Vec<CandidateAnchor>,
    metrics: CandidateMetrics,
}

impl CandidateSet {
    /// Returns the exact normalized query.
    #[must_use]
    pub fn query(&self) -> &NormalizedSequence {
        &self.query.owner.query
    }

    /// Returns an opaque exact query-instance identifier.
    #[must_use]
    pub fn query_instance_id(&self) -> QueryInstanceId {
        self.query.clone()
    }

    /// Returns an opaque exact reference-instance identifier.
    #[must_use]
    pub fn reference_instance_id(&self) -> ReferenceInstanceId {
        self.reference.clone()
    }

    /// Reports whether this set belongs to the supplied query instance.
    #[must_use]
    pub fn belongs_to_query(&self, query: &QueryInstanceId) -> bool {
        self.query.is_same_instance(query)
    }

    /// Reports whether this set belongs to the supplied reference instance.
    #[must_use]
    pub fn belongs_to_reference(&self, reference: &ReferenceInstanceId) -> bool {
        self.reference.is_same_instance(reference)
    }

    /// Returns deterministic candidate anchors.
    #[must_use]
    pub fn anchors(&self) -> &[CandidateAnchor] {
        &self.anchors
    }

    /// Returns complete deterministic result metrics.
    #[must_use]
    pub const fn metrics(&self) -> CandidateMetrics {
        self.metrics
    }
}

impl fmt::Debug for CandidateSet {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("CandidateSet")
            .field("metrics", &self.metrics)
            .field("anchors", &self.anchors)
            .finish_non_exhaustive()
    }
}

struct RetainedMatches {
    request_ordinal: u64,
    request: FixedSeedRequest,
    matches: ProjectedMatches,
}

#[cfg(test)]
struct RawEvidence {
    contig: ContigId,
    strand: BisulfiteStrand,
    diagonal: CandidateDiagonal,
    request_ordinal: u64,
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

fn validate_final_candidate_invariants(
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

fn validate_hit_semantics(
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
fn count_unique_candidates(raw: &[RawEvidence]) -> Result<u64, CandidateError> {
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

fn candidate_physical_to_logical(
    counter: CandidateCounter,
    value: usize,
) -> Result<u64, CandidateError> {
    u64::try_from(value).map_err(|_| CandidateError::CountNotRepresentable { counter, value })
}

fn checked_candidate_add(
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

fn ensure_candidate_capacity(
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

fn preflight_candidate_allocation<T>(
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

#[cfg(test)]
mod tests {
    use super::*;
    use bsbit_core::sequence::normalize_dna;
    use bsbit_index::reference::{ContigInput, ReferenceBuildLimits};

    fn normalized(raw: &[u8]) -> NormalizedSequence {
        normalize_dna(raw).expect("test sequence is normalized")
    }

    fn interval(start: u64, end: u64, length: u64) -> QueryInterval {
        QueryInterval::new(start, end, QueryLength::new(length)).expect("test interval is valid")
    }

    #[test]
    fn diagonal_has_full_range_unique_zero_and_mathematical_order() {
        let values = [
            CandidateDiagonal::before_contig(NonZeroU64::new(u64::MAX).unwrap()),
            CandidateDiagonal::before_contig(NonZeroU64::new(1).unwrap()),
            CandidateDiagonal::at_or_after_contig(0),
            CandidateDiagonal::at_or_after_contig(u64::MAX),
        ];
        assert!(values.windows(2).all(|pair| pair[0] < pair[1]));
        assert_eq!(CandidateDiagonal::from_difference(0, u64::MAX), values[0]);
        assert_eq!(CandidateDiagonal::from_difference(u64::MAX, 0), values[3]);
        assert_eq!(
            CandidateDiagonal::from_difference(4, 4).shift(),
            CoordinateShift::Zero
        );
    }

    #[test]
    fn plan_validation_is_supplied_ordered_and_canonicalizes_only_valid_inputs() {
        let query = normalized(b"ACNGT");
        let empty = FixedSeedRequest::new(BisulfiteStrand::OT, interval(1, 1, 5));
        let n_seed = FixedSeedRequest::new(BisulfiteStrand::OT, interval(1, 4, 5));
        assert!(matches!(
            FixedSeedPlan::new(query.clone(), &[empty, n_seed], SeedPlanLimits::MAX),
            Err(SeedPlanError::EmptySeed {
                request_ordinal: 0,
                ..
            })
        ));
        assert_eq!(
            FixedSeedPlan::new(query, &[n_seed], SeedPlanLimits::MAX).unwrap_err(),
            SeedPlanError::UnsearchableBase {
                request_ordinal: 0,
                query_offset: 2,
            }
        );

        let valid_query = normalized(b"ACGT");
        let later = FixedSeedRequest::new(BisulfiteStrand::CTOB, interval(2, 4, 4));
        let earlier = FixedSeedRequest::new(BisulfiteStrand::OT, interval(0, 2, 4));
        let plan = FixedSeedPlan::new(valid_query, &[later, earlier], SeedPlanLimits::new(2, 4))
            .expect("valid plan");
        assert_eq!(plan.requests(), &[earlier, later]);
        assert_eq!(
            plan.metrics(),
            SeedPlanMetrics {
                query_bases: 4,
                request_count: 2,
                total_seed_bases: 4,
            }
        );
    }

    #[test]
    fn duplicate_plan_and_prefix_limit_fail_without_publication() {
        let query = normalized(b"ACGT");
        let request = FixedSeedRequest::new(BisulfiteStrand::OT, interval(0, 2, 4));
        assert_eq!(
            FixedSeedPlan::new(query.clone(), &[request, request], SeedPlanLimits::MAX)
                .unwrap_err(),
            SeedPlanError::DuplicateRequest {
                strand: BisulfiteStrand::OT,
                interval: request.interval(),
            }
        );
        let later = FixedSeedRequest::new(BisulfiteStrand::OB, interval(2, 4, 4));
        assert_eq!(
            FixedSeedPlan::new(query, &[request, later], SeedPlanLimits::new(2, 3)).unwrap_err(),
            SeedPlanError::TotalSeedBasesLimitExceeded {
                request_ordinal: 1,
                requested: 4,
                maximum: 3,
            }
        );
    }

    #[test]
    fn owner_identity_clone_and_empty_candidate_set_are_exact() {
        let query = normalized(b"ACGT");
        let plan = FixedSeedPlan::new(query, &[], SeedPlanLimits::MAX).expect("empty plan");
        let copy = plan.try_clone().expect("fallible clone");
        assert!(
            plan.query_instance_id()
                .is_same_instance(&copy.query_instance_id())
        );

        let reference = ReferenceIndex::build(
            vec![ContigInput::new(b"c".to_vec(), normalized(b"ACGT"))],
            ReferenceBuildLimits::MAX,
        )
        .expect("reference");
        let set = candidates_for_fixed_seeds(
            &reference,
            &plan,
            ReferenceQueryLimits::MAX,
            CandidateLimits::new(0, 0),
        )
        .expect("empty result");
        assert!(set.anchors().is_empty());
        assert!(set.belongs_to_query(&plan.query_instance_id()));
        assert!(set.belongs_to_reference(&reference.instance_id()));
        assert_eq!(set.metrics().maximum_support(), 0);
    }

    #[test]
    fn allocation_preflight_preserves_context_and_rejects_extreme_counts() {
        let error = preflight_candidate_allocation::<RawEvidence>(
            u64::MAX,
            CandidateAllocation::RawEvidence,
        )
        .unwrap_err();
        assert!(matches!(
            error,
            CandidateError::AllocationSizeOverflow {
                allocation: CandidateAllocation::RawEvidence,
                elements: u64::MAX,
                ..
            }
        ));
        let error = preflight_seed_allocation::<FixedSeedRequest>(
            u64::MAX,
            CandidateAllocation::CanonicalRequests,
        )
        .unwrap_err();
        assert!(matches!(
            error,
            SeedPlanError::AllocationSizeOverflow {
                allocation: CandidateAllocation::CanonicalRequests,
                elements: u64::MAX,
                ..
            }
        ));
    }

    #[test]
    fn defensive_hit_semantics_report_exact_request_context() {
        let request = FixedSeedRequest::new(BisulfiteStrand::OT, interval(1, 3, 4));
        assert_eq!(
            validate_hit_semantics(7, request, BisulfiteStrand::OT, 2),
            Ok(())
        );
        assert_eq!(
            validate_hit_semantics(7, request, BisulfiteStrand::OB, 2),
            Err(CandidateError::HitStrandMismatch {
                request_ordinal: 7,
                expected: BisulfiteStrand::OT,
                observed: BisulfiteStrand::OB,
            })
        );
        assert_eq!(
            validate_hit_semantics(7, request, BisulfiteStrand::OT, 3),
            Err(CandidateError::HitLengthMismatch {
                request_ordinal: 7,
                expected: 2,
                observed: 3,
            })
        );
    }

    #[test]
    fn seed_plan_validation_errors_have_stable_fields_display_and_source_policy() {
        let empty = interval(0, 0, 1);
        let coordinate = CoordinateError::OutOfBounds {
            domain: bsbit_core::coordinate::CoordinateDomain::Query,
            operation: bsbit_core::coordinate::CoordinateOperation::IntervalConstruction,
            start: 2,
            end: 3,
            length: 1,
        };
        let errors = vec![
            (
                SeedPlanError::RequestCountNotRepresentable { value: 9 },
                "physical seed request count 9 is not representable as u64",
                false,
            ),
            (
                SeedPlanError::RequestLimitExceeded {
                    requested: 3,
                    maximum: 2,
                },
                "seed request count 3 exceeds configured maximum 2",
                false,
            ),
            (
                SeedPlanError::InvalidInterval {
                    request_ordinal: 4,
                    source: coordinate,
                },
                "seed request 4 is invalid for the actual query: Query interval [2, 3) is outside length 1 during IntervalConstruction",
                true,
            ),
            (
                SeedPlanError::EmptySeed {
                    request_ordinal: 5,
                    interval: empty,
                },
                "seed request 5 has empty interval query:[0,0)",
                false,
            ),
            (
                SeedPlanError::UnsearchableBase {
                    request_ordinal: 6,
                    query_offset: 7,
                },
                "seed request 6 contains unsearchable N at absolute query offset 7",
                false,
            ),
            (
                SeedPlanError::SeedOffsetNotRepresentable {
                    request_ordinal: 8,
                    value: 9,
                },
                "seed request 8 local offset 9 is not representable as u64",
                false,
            ),
            (
                SeedPlanError::QueryOffsetOverflow {
                    request_ordinal: 10,
                    start: u64::MAX,
                    local_offset: 1,
                },
                "seed request 10 absolute query offset 18446744073709551615 plus local offset 1 overflowed",
                false,
            ),
            (
                SeedPlanError::BoundaryNotRepresentable {
                    request_ordinal: 11,
                    boundary: QueryBoundary::Start,
                    value: 12,
                },
                "seed request 11 Start boundary 12 does not fit this architecture",
                false,
            ),
        ];
        for (error, display, has_source) in errors {
            assert_eq!(error.to_string(), display);
            assert_eq!(std::error::Error::source(&error).is_some(), has_source);
        }
    }

    #[test]
    fn seed_plan_resource_errors_have_stable_fields_display_and_source_policy() {
        let empty = interval(0, 0, 1);
        let errors = vec![
            (
                SeedPlanError::TotalSeedBasesOverflow {
                    accumulated: u64::MAX,
                    next: 1,
                },
                "total seed bases 18446744073709551615 plus 1 overflowed",
                false,
            ),
            (
                SeedPlanError::TotalSeedBasesLimitExceeded {
                    request_ordinal: 13,
                    requested: 15,
                    maximum: 14,
                },
                "seed request 13 raises prefix total to 15, exceeding 14",
                false,
            ),
            (
                SeedPlanError::AllocationSizeOverflow {
                    allocation: CandidateAllocation::CanonicalRequests,
                    elements: 16,
                    element_size: 24,
                },
                "cannot size CanonicalRequests: 16 elements of 24 bytes",
                false,
            ),
            (
                SeedPlanError::AllocationFailed {
                    allocation: CandidateAllocation::CanonicalRequests,
                    elements: 17,
                },
                "failed to reserve 17 elements for CanonicalRequests",
                false,
            ),
            (
                SeedPlanError::DuplicateRequest {
                    strand: BisulfiteStrand::OT,
                    interval: empty,
                },
                "duplicate seed request OT query:[0,0)",
                false,
            ),
            (
                SeedPlanError::CapacityInvariant {
                    reserved: 18,
                    materialized: 19,
                },
                "canonical request reservation 18 cannot accept entry 19",
                false,
            ),
        ];
        for (error, display, has_source) in errors {
            assert_eq!(error.to_string(), display);
            assert_eq!(std::error::Error::source(&error).is_some(), has_source);
        }
    }

    #[test]
    fn candidate_search_errors_have_stable_fields_display_and_source_policy() {
        let seed_interval = interval(1, 2, 3);
        let errors = vec![
            (
                CandidateError::AllocationSizeOverflow {
                    allocation: CandidateAllocation::RawEvidence,
                    elements: 2,
                    element_size: 32,
                },
                "cannot size RawEvidence: 2 elements of 32 bytes",
                false,
            ),
            (
                CandidateError::AllocationFailed {
                    allocation: CandidateAllocation::FinalAnchors,
                    elements: 3,
                },
                "failed to reserve 3 elements for FinalAnchors",
                false,
            ),
            (
                CandidateError::CountNotRepresentable {
                    counter: CandidateCounter::Requests,
                    value: 4,
                },
                "Requests physical count 4 is not representable as u64",
                false,
            ),
            (
                CandidateError::CounterOverflow {
                    counter: CandidateCounter::SupportSum,
                    accumulated: u64::MAX,
                    next: 1,
                },
                "SupportSum count 18446744073709551615 plus 1 overflowed",
                false,
            ),
            (
                CandidateError::Search {
                    request_ordinal: 5,
                    strand: BisulfiteStrand::OT,
                    interval: seed_interval,
                    source: ReferenceQueryError::EmptyPattern,
                },
                "candidate search failed for request 5 OT query:[1,2): exact-search pattern is empty",
                true,
            ),
            (
                CandidateError::AggregateHitCountOverflow {
                    accumulated: u64::MAX,
                    request_hits: 1,
                },
                "aggregate exact hits 18446744073709551615 plus request count 1 overflowed",
                false,
            ),
            (
                CandidateError::AggregateHitLimitExceeded {
                    accumulated: 6,
                    request_hits: 7,
                    requested: 13,
                    maximum: 12,
                },
                "aggregate exact hits 6 plus request count 7 is 13, exceeding 12",
                false,
            ),
            (
                CandidateError::Locate {
                    request_ordinal: 8,
                    strand: BisulfiteStrand::OB,
                    interval: seed_interval,
                    source: ReferenceLocateError::ForeignMatches,
                },
                "candidate locate failed for request 8 OB query:[1,2): projected matches belong to another reference instance",
                true,
            ),
        ];
        for (error, display, has_source) in errors {
            assert_eq!(error.to_string(), display);
            assert_eq!(std::error::Error::source(&error).is_some(), has_source);
        }
    }

    #[test]
    fn candidate_evidence_errors_have_stable_fields_display_and_source_policy() {
        let reverse_error = CoordinateError::CoordinateUnderflow {
            domain: bsbit_core::coordinate::CoordinateDomain::Query,
            operation: bsbit_core::coordinate::CoordinateOperation::ReverseTransform,
            lhs: 0,
            rhs: 1,
        };
        let negative = CandidateDiagonal::before_contig(NonZeroU64::new(2).unwrap());
        let errors = vec![
            (
                CandidateError::HitStrandMismatch {
                    request_ordinal: 9,
                    expected: BisulfiteStrand::OT,
                    observed: BisulfiteStrand::CTOB,
                },
                "request 9 expected hit strand OT, observed CTOB",
                false,
            ),
            (
                CandidateError::HitLengthMismatch {
                    request_ordinal: 10,
                    expected: 11,
                    observed: 12,
                },
                "request 10 expected hit length 11, observed 12",
                false,
            ),
            (
                CandidateError::OrientedInterval {
                    request_ordinal: 13,
                    source: reverse_error,
                },
                "request 13 oriented interval recovery failed: Query coordinate subtraction 0 - 1 underflowed during ReverseTransform",
                true,
            ),
            (
                CandidateError::PlanIntervalStorage {
                    request_ordinal: 14,
                    start: 15,
                    end: 16,
                    query_bases: 17,
                },
                "request 14 interval [15,16) is not physically addressable in query length 17",
                false,
            ),
            (
                CandidateError::DuplicateRequestEvidence {
                    request_ordinal: 18,
                    contig_ordinal: 19,
                    strand: BisulfiteStrand::CTOT,
                    diagonal: negative,
                },
                "request 18 produced duplicate evidence for contig 19, strand CTOT, diagonal -2",
                false,
            ),
            (
                CandidateError::UniqueCandidateLimitExceeded {
                    requested: 21,
                    maximum: 20,
                },
                "unique candidate count 21 exceeds configured maximum 20",
                false,
            ),
            (
                CandidateError::Invariant {
                    invariant: CandidateInvariant::OutputOrder,
                    expected: 1,
                    observed: 0,
                },
                "OutputOrder invariant expected 1, observed 0",
                false,
            ),
        ];
        for (error, display, has_source) in errors {
            assert_eq!(error.to_string(), display);
            assert_eq!(std::error::Error::source(&error).is_some(), has_source);
        }
    }

    #[test]
    fn private_final_candidate_invariants_fail_in_normative_order() {
        assert_eq!(validate_final_candidate_invariants(5, 3, 5, 3, true), Ok(2));
        assert_eq!(
            validate_final_candidate_invariants(5, 3, 4, 2, false),
            Err(CandidateError::Invariant {
                invariant: CandidateInvariant::SupportSum,
                expected: 5,
                observed: 4,
            })
        );
        assert_eq!(
            validate_final_candidate_invariants(5, 3, 5, 2, false),
            Err(CandidateError::Invariant {
                invariant: CandidateInvariant::CandidateCount,
                expected: 3,
                observed: 2,
            })
        );
        assert_eq!(
            validate_final_candidate_invariants(5, 3, 5, 3, false),
            Err(CandidateError::Invariant {
                invariant: CandidateInvariant::OutputOrder,
                expected: 1,
                observed: 0,
            })
        );
        assert_eq!(
            validate_final_candidate_invariants(2, 3, 2, 3, true),
            Err(CandidateError::Invariant {
                invariant: CandidateInvariant::DuplicateEvidence,
                expected: 2,
                observed: 3,
            })
        );
    }

    #[test]
    fn duplicate_request_evidence_and_counter_guards_are_structured() {
        let reference = ReferenceIndex::build(
            vec![ContigInput::new(b"c".to_vec(), normalized(b"ACGT"))],
            ReferenceBuildLimits::MAX,
        )
        .expect("reference");
        let contig = reference.contig_id(0).expect("contig");
        let diagonal = CandidateDiagonal::at_or_after_contig(1);
        let raw = [
            RawEvidence {
                contig: contig.clone(),
                strand: BisulfiteStrand::OT,
                diagonal,
                request_ordinal: 4,
            },
            RawEvidence {
                contig,
                strand: BisulfiteStrand::OT,
                diagonal,
                request_ordinal: 4,
            },
        ];
        assert_eq!(
            count_unique_candidates(&raw),
            Err(CandidateError::DuplicateRequestEvidence {
                request_ordinal: 4,
                contig_ordinal: 0,
                strand: BisulfiteStrand::OT,
                diagonal,
            })
        );
        assert_eq!(
            checked_candidate_add(CandidateCounter::SupportSum, u64::MAX, 1),
            Err(CandidateError::CounterOverflow {
                counter: CandidateCounter::SupportSum,
                accumulated: u64::MAX,
                next: 1,
            })
        );
        assert_eq!(
            ensure_candidate_capacity(CandidateInvariant::RawEvidenceCapacity, 1, 1),
            Err(CandidateError::Invariant {
                invariant: CandidateInvariant::RawEvidenceCapacity,
                expected: 1,
                observed: 1,
            })
        );
    }
}
