//! Fixed-seed bisulfite candidate anchoring and deterministic support grouping.
//!
//! This module composes the owner-bound projected-reference API with explicit
//! raw-query seed intervals. Candidate anchors are pre-extension evidence, not
//! verified alignments, placements, ambiguity decisions, or MAPQ values.
//!
//! Query instance identifiers deliberately have no generic equality:
//!
//! ```compile_fail
//! use bsbit_align::search::candidate::QueryInstanceId;
//!
//! fn requires_eq<T: Eq>() {}
//! requires_eq::<QueryInstanceId>();
//! ```
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
use std::sync::Arc;

use bsbit_core::bisulfite::BisulfiteStrand;
use bsbit_core::coordinate::{CoordinateError, CoordinateShift, QueryInterval, QueryLength};
use bsbit_core::sequence::NormalizedSequence;
use bsbit_index::reference::{
    ContigId, ReferenceAccessError, ReferenceInstanceId, ReferenceLocateError, ReferenceQueryError,
};
#[cfg(test)]
use bsbit_index::reference::{ReferenceIndex, ReferenceQueryLimits};

pub use super::candidate_generation::candidates_for_fixed_seeds;
#[cfg(test)]
use super::candidate_generation::{
    RawEvidence, checked_candidate_add, count_unique_candidates, ensure_candidate_capacity,
    preflight_candidate_allocation, validate_final_candidate_invariants, validate_hit_semantics,
};
use super::candidate_generation::{
    preflight_seed_allocation, request_count_to_u64, strand_rank, validate_supplied_requests,
};

/// One explicit fixed seed in raw sequencing-order query coordinates.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct FixedSeedRequest {
    pub(super) strand: BisulfiteStrand,
    pub(super) interval: QueryInterval,
}

impl FixedSeedRequest {
    /// Creates one explicit strand and raw-query interval request.
    #[must_use]
    pub const fn new(strand: BisulfiteStrand, interval: QueryInterval) -> Self {
        Self { strand, interval }
    }

    /// Returns the requested bisulfite strand.
    #[must_use]
    pub const fn strand(self) -> BisulfiteStrand {
        self.strand
    }

    /// Returns the raw sequencing-order query interval.
    #[must_use]
    pub const fn interval(self) -> QueryInterval {
        self.interval
    }
}

/// Complete construction limits for one fixed seed plan.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SeedPlanLimits {
    pub(super) max_requests: u64,
    pub(super) max_total_seed_bases: u64,
}

impl SeedPlanLimits {
    /// Limits admitting every representable logical plan.
    pub const MAX: Self = Self {
        max_requests: u64::MAX,
        max_total_seed_bases: u64::MAX,
    };

    /// Creates explicit request-count and aggregate seed-base limits.
    #[must_use]
    pub const fn new(max_requests: u64, max_total_seed_bases: u64) -> Self {
        Self {
            max_requests,
            max_total_seed_bases,
        }
    }

    /// Returns the maximum request count.
    #[must_use]
    pub const fn max_requests(self) -> u64 {
        self.max_requests
    }

    /// Returns the maximum aggregate seed bases.
    #[must_use]
    pub const fn max_total_seed_bases(self) -> u64 {
        self.max_total_seed_bases
    }
}

impl Default for SeedPlanLimits {
    fn default() -> Self {
        Self::MAX
    }
}

/// Complete deterministic dimensions of a fixed seed plan.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SeedPlanMetrics {
    pub(super) query_bases: u64,
    pub(super) request_count: u64,
    pub(super) total_seed_bases: u64,
}

impl SeedPlanMetrics {
    /// Returns the normalized query length.
    #[must_use]
    pub const fn query_bases(self) -> u64 {
        self.query_bases
    }

    /// Returns the number of canonical distinct requests.
    #[must_use]
    pub const fn request_count(self) -> u64 {
        self.request_count
    }

    /// Returns the aggregate seed bases, counting overlap per request.
    #[must_use]
    pub const fn total_seed_bases(self) -> u64 {
        self.total_seed_bases
    }
}

/// A variable-sized allocation owned by the fixed-seed/candidate layer.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum CandidateAllocation {
    /// Canonical fixed seed requests.
    CanonicalRequests,
    /// Retained owner-bound per-request match artifacts.
    RetainedMatches,
    /// One raw evidence record per exact occurrence.
    RawEvidence,
    /// Lightweight candidate keys for one exact-request locate stream.
    RequestCandidateKeys,
    /// Globally merged unique candidate keys and support counters.
    CandidateVotes,
    /// One final anchor per unique pre-extension key.
    FinalAnchors,
}

/// A query interval boundary that failed physical conversion.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum QueryBoundary {
    /// Inclusive interval start.
    Start,
    /// Exclusive interval end.
    End,
}

/// A fixed seed plan construction failure.
#[non_exhaustive]
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum SeedPlanError {
    /// The physical request count cannot fit the logical width.
    RequestCountNotRepresentable {
        /// Physical request count.
        value: usize,
    },
    /// The request count exceeds its configured limit.
    RequestLimitExceeded {
        /// Requested count.
        requested: u64,
        /// Configured maximum.
        maximum: u64,
    },
    /// A request interval is invalid for the actual query.
    InvalidInterval {
        /// Supplied-order request ordinal.
        request_ordinal: u64,
        /// Underlying coordinate failure.
        source: CoordinateError,
    },
    /// A request interval is empty.
    EmptySeed {
        /// Supplied-order request ordinal.
        request_ordinal: u64,
        /// Empty interval.
        interval: QueryInterval,
    },
    /// A normalized seed contains N.
    UnsearchableBase {
        /// Supplied-order request ordinal.
        request_ordinal: u64,
        /// Absolute zero-based query offset.
        query_offset: u64,
    },
    /// A physical seed offset cannot fit the logical width.
    SeedOffsetNotRepresentable {
        /// Supplied-order request ordinal.
        request_ordinal: u64,
        /// Physical local offset.
        value: usize,
    },
    /// An absolute query offset overflowed while identifying an unsearchable base.
    QueryOffsetOverflow {
        /// Supplied-order request ordinal.
        request_ordinal: u64,
        /// Raw interval start.
        start: u64,
        /// Local offset within the interval.
        local_offset: u64,
    },
    /// A query boundary cannot fit this architecture.
    BoundaryNotRepresentable {
        /// Supplied-order request ordinal.
        request_ordinal: u64,
        /// Failed interval boundary.
        boundary: QueryBoundary,
        /// Logical boundary value.
        value: u64,
    },
    /// Aggregate seed bases overflowed.
    TotalSeedBasesOverflow {
        /// Prefix total before the current request.
        accumulated: u64,
        /// Current seed length.
        next: u64,
    },
    /// The first supplied-order prefix exceeds its aggregate limit.
    TotalSeedBasesLimitExceeded {
        /// Supplied-order request ordinal that first exceeds the limit.
        request_ordinal: u64,
        /// First exceeding prefix total.
        requested: u64,
        /// Configured maximum.
        maximum: u64,
    },
    /// Canonical request storage cannot fit this architecture.
    AllocationSizeOverflow {
        /// Allocation site.
        allocation: CandidateAllocation,
        /// Requested elements.
        elements: u64,
        /// Element width.
        element_size: u64,
    },
    /// Canonical request reservation failed.
    AllocationFailed {
        /// Allocation site.
        allocation: CandidateAllocation,
        /// Requested elements.
        elements: u64,
    },
    /// Two supplied requests have the same exact semantic key.
    DuplicateRequest {
        /// Duplicated strand.
        strand: BisulfiteStrand,
        /// Duplicated raw query interval.
        interval: QueryInterval,
    },
    /// A private plan capacity invariant failed.
    CapacityInvariant {
        /// Exact reserved count.
        reserved: u64,
        /// Number already materialized.
        materialized: u64,
    },
}

impl fmt::Display for SeedPlanError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::RequestCountNotRepresentable { value } => {
                write!(
                    formatter,
                    "physical seed request count {value} is not representable as u64"
                )
            }
            Self::RequestLimitExceeded { requested, maximum } => write!(
                formatter,
                "seed request count {requested} exceeds configured maximum {maximum}"
            ),
            Self::InvalidInterval {
                request_ordinal,
                source,
            } => write!(
                formatter,
                "seed request {request_ordinal} is invalid for the actual query: {source}"
            ),
            Self::EmptySeed {
                request_ordinal,
                interval,
            } => write!(
                formatter,
                "seed request {request_ordinal} has empty interval {interval}"
            ),
            Self::UnsearchableBase {
                request_ordinal,
                query_offset,
            } => write!(
                formatter,
                "seed request {request_ordinal} contains unsearchable N at absolute query offset {query_offset}"
            ),
            Self::SeedOffsetNotRepresentable {
                request_ordinal,
                value,
            } => write!(
                formatter,
                "seed request {request_ordinal} local offset {value} is not representable as u64"
            ),
            Self::QueryOffsetOverflow {
                request_ordinal,
                start,
                local_offset,
            } => write!(
                formatter,
                "seed request {request_ordinal} absolute query offset {start} plus local offset {local_offset} overflowed"
            ),
            Self::BoundaryNotRepresentable {
                request_ordinal,
                boundary,
                value,
            } => write!(
                formatter,
                "seed request {request_ordinal} {boundary:?} boundary {value} does not fit this architecture"
            ),
            Self::TotalSeedBasesOverflow { accumulated, next } => write!(
                formatter,
                "total seed bases {accumulated} plus {next} overflowed"
            ),
            Self::TotalSeedBasesLimitExceeded {
                request_ordinal,
                requested,
                maximum,
            } => write!(
                formatter,
                "seed request {request_ordinal} raises prefix total to {requested}, exceeding {maximum}"
            ),
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
            Self::DuplicateRequest { strand, interval } => {
                write!(formatter, "duplicate seed request {strand:?} {interval}")
            }
            Self::CapacityInvariant {
                reserved,
                materialized,
            } => write!(
                formatter,
                "canonical request reservation {reserved} cannot accept entry {materialized}"
            ),
        }
    }
}

impl std::error::Error for SeedPlanError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::InvalidInterval { source, .. } => Some(source),
            _ => None,
        }
    }
}

struct QueryOwner {
    query: NormalizedSequence,
}

/// An opaque process-local identifier for one exact query instance.
#[derive(Clone)]
pub struct QueryInstanceId {
    owner: Arc<QueryOwner>,
}

impl QueryInstanceId {
    /// Reports exact shared runtime ownership.
    #[must_use]
    pub fn is_same_instance(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.owner, &other.owner)
    }
}

impl fmt::Debug for QueryInstanceId {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("QueryInstanceId")
            .finish_non_exhaustive()
    }
}

/// One immutable normalized query and canonical distinct fixed seed plan.
pub struct FixedSeedPlan {
    owner: Arc<QueryOwner>,
    pub(super) requests: Vec<FixedSeedRequest>,
    pub(super) metrics: SeedPlanMetrics,
}

impl FixedSeedPlan {
    /// Validates, canonicalizes, and owns a complete fixed seed plan.
    ///
    /// Validation is performed in supplied order before canonical request
    /// allocation. The first invalid request or first exceeding seed-base
    /// prefix is returned, and no partial plan is published.
    ///
    /// # Errors
    ///
    /// Returns [`SeedPlanError`] for invalid/duplicate requests, overflow,
    /// configured limits, architecture sizing, or fallible allocation failure.
    pub fn new(
        query: NormalizedSequence,
        supplied: &[FixedSeedRequest],
        limits: SeedPlanLimits,
    ) -> Result<Self, SeedPlanError> {
        let request_count = request_count_to_u64(supplied.len())?;
        if request_count > limits.max_requests {
            return Err(SeedPlanError::RequestLimitExceeded {
                requested: request_count,
                maximum: limits.max_requests,
            });
        }

        let query_length = QueryLength::new(query.len());
        let total_seed_bases = validate_supplied_requests(&query, supplied, query_length, limits)?;

        let storage = preflight_seed_allocation::<FixedSeedRequest>(
            request_count,
            CandidateAllocation::CanonicalRequests,
        )?;
        let mut requests = Vec::new();
        requests
            .try_reserve_exact(storage)
            .map_err(|_| SeedPlanError::AllocationFailed {
                allocation: CandidateAllocation::CanonicalRequests,
                elements: request_count,
            })?;
        for request in supplied {
            let materialized = request_count_to_u64(requests.len())?;
            if materialized >= request_count {
                return Err(SeedPlanError::CapacityInvariant {
                    reserved: request_count,
                    materialized,
                });
            }
            requests.push(*request);
        }
        requests.sort_unstable_by_key(|request| {
            (
                strand_rank(request.strand),
                request.interval.start(),
                request.interval.end(),
            )
        });
        if let Some(pair) = requests.windows(2).find(|pair| pair[0] == pair[1]) {
            return Err(SeedPlanError::DuplicateRequest {
                strand: pair[0].strand,
                interval: pair[0].interval,
            });
        }

        let owner = Arc::new(QueryOwner { query });
        Ok(Self {
            owner,
            requests,
            metrics: SeedPlanMetrics {
                query_bases: query_length.get(),
                request_count,
                total_seed_bases,
            },
        })
    }

    /// Returns the exact immutable normalized query.
    #[must_use]
    pub fn query(&self) -> &NormalizedSequence {
        &self.owner.query
    }

    /// Returns the typed query length.
    #[must_use]
    pub const fn query_length(&self) -> QueryLength {
        QueryLength::new(self.metrics.query_bases)
    }

    /// Returns an opaque exact query-instance identifier.
    #[must_use]
    pub fn query_instance_id(&self) -> QueryInstanceId {
        QueryInstanceId {
            owner: Arc::clone(&self.owner),
        }
    }

    /// Returns canonical requests in explicit strand/start/end order.
    #[must_use]
    pub fn requests(&self) -> &[FixedSeedRequest] {
        &self.requests
    }

    /// Returns complete deterministic plan dimensions.
    #[must_use]
    pub const fn metrics(&self) -> SeedPlanMetrics {
        self.metrics
    }

    /// Fallibly copies canonical requests while sharing the exact query owner.
    ///
    /// # Errors
    ///
    /// Returns [`SeedPlanError`] when request storage cannot fit this
    /// architecture or its reservation fails.
    pub fn try_clone(&self) -> Result<Self, SeedPlanError> {
        let request_count = self.metrics.request_count;
        let storage = preflight_seed_allocation::<FixedSeedRequest>(
            request_count,
            CandidateAllocation::CanonicalRequests,
        )?;
        let mut requests = Vec::new();
        requests
            .try_reserve_exact(storage)
            .map_err(|_| SeedPlanError::AllocationFailed {
                allocation: CandidateAllocation::CanonicalRequests,
                elements: request_count,
            })?;
        requests.extend_from_slice(&self.requests);
        Ok(Self {
            owner: Arc::clone(&self.owner),
            requests,
            metrics: self.metrics,
        })
    }
}

impl fmt::Debug for FixedSeedPlan {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("FixedSeedPlan")
            .field("metrics", &self.metrics)
            .field("requests", &self.requests)
            .finish_non_exhaustive()
    }
}

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
    pub(super) max_total_exact_hits: u64,
    pub(super) max_unique_candidates: u64,
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
    pub(super) request_count: u64,
    pub(super) total_seed_bases: u64,
    pub(super) total_exact_hits: u64,
    pub(super) matched_intervals: u64,
    pub(super) unique_candidates: u64,
    pub(super) duplicate_evidence: u64,
    pub(super) maximum_support: u64,
    pub(super) zero_hit_requests: u64,
    pub(super) search_rank_operations: u64,
    pub(super) locate_calls: u64,
    pub(super) located_coordinates: u64,
    pub(super) locate_lf_steps: u64,
    pub(super) locate_rank_operations: u64,
    pub(super) locate_interval_nodes: u64,
    pub(super) candidate_key_materializations: u64,
    pub(super) peak_request_candidate_keys: u64,
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
    pub(super) contig: ContigId,
    pub(super) strand: BisulfiteStrand,
    pub(super) diagonal: CandidateDiagonal,
    pub(super) support: NonZeroU64,
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
    pub(super) reference: ReferenceInstanceId,
    pub(super) query: QueryInstanceId,
    pub(super) anchors: Vec<CandidateAnchor>,
    pub(super) metrics: CandidateMetrics,
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
