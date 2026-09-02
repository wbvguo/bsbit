//! Fixed-seed request planning and query ownership.
//!
//! Plans validate and canonicalize explicit raw-query intervals before the
//! candidate layer consults a reference. Query instance identifiers preserve
//! exact runtime ownership across the planning and candidate stages.
//!
//! Query instance identifiers deliberately have no generic equality:
//!
//! ```compile_fail
//! use bsbit_align::search::fixed_seed::QueryInstanceId;
//!
//! fn requires_eq<T: Eq>() {}
//! requires_eq::<QueryInstanceId>();
//! ```

use core::fmt;
use core::mem::size_of;
use std::sync::Arc;

use bsbit_core::alphabet::Base;
use bsbit_core::bisulfite::BisulfiteStrand;
use bsbit_core::coordinate::{CoordinateError, QueryInterval, QueryLength};
use bsbit_core::sequence::NormalizedSequence;

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
    max_requests: u64,
    max_total_seed_bases: u64,
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

pub(super) struct QueryOwner {
    pub(super) query: NormalizedSequence,
}

/// An opaque process-local identifier for one exact query instance.
#[derive(Clone)]
pub struct QueryInstanceId {
    pub(super) owner: Arc<QueryOwner>,
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
    pub(super) owner: Arc<QueryOwner>,
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

fn validate_supplied_requests(
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

pub(super) fn preflight_storage<T>(elements: u64) -> Result<usize, (u64, u64)> {
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
