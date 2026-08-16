//! Proof-carrying, reference-independent seed scheduling.
//!
//! For a query of length `m` and full-query unit edit budget `k`, the
//! reference scheduler partitions raw query coordinates into `k + 1`
//! non-empty blocks when `m > k`. Every N-free block is requested on every
//! explicitly admitted bisulfite strand. The resulting certificate guarantees
//! candidate evidence recall within the stated budget when downstream exact
//! search completes without a resource error.
//!
//! This module does not inspect a reference, infer a library profile, verify an
//! alignment, or implement the original 18/30/1000 adaptive heuristic.

use core::fmt;
use core::iter::FusedIterator;
use core::mem::size_of;

use crate::score::EditDistance;
use crate::search::candidate::{FixedSeedPlan, FixedSeedRequest, SeedPlanError, SeedPlanLimits};
use bsbit_core::bisulfite::BisulfiteStrand;
use bsbit_core::coordinate::{CoordinateError, QueryInterval, QueryLength};
use bsbit_core::sequence::NormalizedSequence;

const STRAND_COUNT: u8 = 4;

/// A caller-owned, nonempty set of admissible bisulfite strands.
///
/// The set is allocation-free and iterates in [`BisulfiteStrand::ALL`] order.
/// It contains no library-profile or mate inference.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct AdmissibleStrands {
    mask: u8,
    count: u8,
}

impl AdmissibleStrands {
    /// Validates an explicit supplied-order strand list.
    ///
    /// # Errors
    ///
    /// Returns [`AdmissibleStrandError`] for an empty list, an unrepresentable
    /// supplied count, or the first duplicate in supplied order.
    pub fn new(supplied: &[BisulfiteStrand]) -> Result<Self, AdmissibleStrandError> {
        if supplied.is_empty() {
            return Err(AdmissibleStrandError::Empty);
        }
        let _supplied_count = u64::try_from(supplied.len()).map_err(|_| {
            AdmissibleStrandError::SuppliedCountNotRepresentable {
                value: supplied.len(),
            }
        })?;

        let mut mask = 0_u8;
        let mut count = 0_u8;
        for (physical_ordinal, strand) in supplied.iter().copied().enumerate() {
            let ordinal = u64::try_from(physical_ordinal).map_err(|_| {
                AdmissibleStrandError::SuppliedOrdinalNotRepresentable {
                    value: physical_ordinal,
                }
            })?;
            let bit = strand_bit(strand);
            if mask & bit != 0 {
                return Err(AdmissibleStrandError::Duplicate { ordinal, strand });
            }
            mask |= bit;
            count = count
                .checked_add(1)
                .ok_or(AdmissibleStrandError::CountInvariant {
                    expected_maximum: STRAND_COUNT,
                    observed: u64::from(count) + 1,
                })?;
        }

        if count == 0 || count > STRAND_COUNT || mask.count_ones() != u32::from(count) {
            return Err(AdmissibleStrandError::CountInvariant {
                expected_maximum: STRAND_COUNT,
                observed: u64::from(count),
            });
        }
        Ok(Self { mask, count })
    }

    /// Returns the number of admitted strands.
    #[must_use]
    pub const fn len(self) -> u64 {
        self.count as u64
    }

    /// Reports whether no strand is admitted.
    ///
    /// Valid public values are never empty; this method supports generic set
    /// inspection without weakening construction.
    #[must_use]
    pub const fn is_empty(self) -> bool {
        self.count == 0
    }

    /// Reports whether `strand` is admitted.
    #[must_use]
    pub const fn contains(self, strand: BisulfiteStrand) -> bool {
        self.mask & strand_bit(strand) != 0
    }

    /// Iterates admitted strands in canonical order.
    #[must_use]
    pub const fn iter(
        self,
    ) -> impl DoubleEndedIterator<Item = BisulfiteStrand> + ExactSizeIterator + FusedIterator {
        AdmissibleStrandIter {
            mask: self.mask,
            front: 0,
            back: STRAND_COUNT,
            remaining: self.count,
        }
    }
}

/// Validation failure for an explicit admissible-strand set.
#[non_exhaustive]
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum AdmissibleStrandError {
    /// No strand was supplied.
    Empty,
    /// The physical supplied count does not fit the logical width.
    SuppliedCountNotRepresentable {
        /// Physical count.
        value: usize,
    },
    /// A supplied ordinal does not fit the logical width.
    SuppliedOrdinalNotRepresentable {
        /// Physical ordinal.
        value: usize,
    },
    /// A strand was supplied more than once.
    Duplicate {
        /// Supplied-order ordinal of the first repeated value.
        ordinal: u64,
        /// Repeated strand.
        strand: BisulfiteStrand,
    },
    /// A private fixed-set cardinality invariant failed.
    CountInvariant {
        /// Maximum representable strands.
        expected_maximum: u8,
        /// Observed count.
        observed: u64,
    },
}

impl fmt::Display for AdmissibleStrandError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => formatter.write_str("admissible strand set is empty"),
            Self::SuppliedCountNotRepresentable { value } => write!(
                formatter,
                "physical admissible strand count {value} is not representable as u64"
            ),
            Self::SuppliedOrdinalNotRepresentable { value } => write!(
                formatter,
                "physical admissible strand ordinal {value} is not representable as u64"
            ),
            Self::Duplicate { ordinal, strand } => write!(
                formatter,
                "admissible strand {strand} is duplicated at supplied ordinal {ordinal}"
            ),
            Self::CountInvariant {
                expected_maximum,
                observed,
            } => write!(
                formatter,
                "admissible strand count invariant expected at most {expected_maximum}, observed {observed}"
            ),
        }
    }
}

impl std::error::Error for AdmissibleStrandError {}

#[derive(Debug)]
struct AdmissibleStrandIter {
    mask: u8,
    front: u8,
    back: u8,
    remaining: u8,
}

impl Iterator for AdmissibleStrandIter {
    type Item = BisulfiteStrand;

    fn next(&mut self) -> Option<Self::Item> {
        while self.front < self.back {
            let rank = self.front;
            self.front += 1;
            let strand = BisulfiteStrand::ALL[usize::from(rank)];
            if self.mask & strand_bit(strand) != 0 {
                self.remaining -= 1;
                return Some(strand);
            }
        }
        None
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let remaining = usize::from(self.remaining);
        (remaining, Some(remaining))
    }
}

impl DoubleEndedIterator for AdmissibleStrandIter {
    fn next_back(&mut self) -> Option<Self::Item> {
        while self.front < self.back {
            self.back -= 1;
            let strand = BisulfiteStrand::ALL[usize::from(self.back)];
            if self.mask & strand_bit(strand) != 0 {
                self.remaining -= 1;
                return Some(strand);
            }
        }
        None
    }
}

impl ExactSizeIterator for AdmissibleStrandIter {
    fn len(&self) -> usize {
        usize::from(self.remaining)
    }
}

impl FusedIterator for AdmissibleStrandIter {}

/// Complete construction limits for one proof-seed schedule.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ProofSeedLimits {
    max_blocks: u64,
    plan_limits: SeedPlanLimits,
}

impl ProofSeedLimits {
    /// Limits admitting every representable logical schedule.
    pub const MAX: Self = Self {
        max_blocks: u64::MAX,
        plan_limits: SeedPlanLimits::MAX,
    };

    /// Creates an explicit partition-block limit and Level 2C plan limits.
    #[must_use]
    pub const fn new(max_blocks: u64, plan_limits: SeedPlanLimits) -> Self {
        Self {
            max_blocks,
            plan_limits,
        }
    }

    /// Returns the maximum partition block count.
    #[must_use]
    pub const fn max_blocks(self) -> u64 {
        self.max_blocks
    }

    /// Returns the exact downstream fixed-plan limits.
    #[must_use]
    pub const fn plan_limits(self) -> SeedPlanLimits {
        self.plan_limits
    }
}

impl Default for ProofSeedLimits {
    fn default() -> Self {
        Self::MAX
    }
}

/// Deterministic proof dimensions for one certified schedule.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SeedCompletenessCertificate {
    query_bases: u64,
    max_edit_distance: EditDistance,
    strand_count: u64,
    block_count: u64,
    emitted_blocks: u64,
    omitted_unknown_blocks: u64,
    unknown_bases: u64,
    request_count: u64,
    total_seed_bases: u64,
    minimum_block_bases: u64,
    maximum_block_bases: u64,
}

impl SeedCompletenessCertificate {
    /// Returns the full normalized query length.
    #[must_use]
    pub const fn query_bases(self) -> u64 {
        self.query_bases
    }

    /// Returns the full-query unit edit budget certified by this schedule.
    #[must_use]
    pub const fn max_edit_distance(self) -> EditDistance {
        self.max_edit_distance
    }

    /// Returns the number of explicitly admitted strands.
    #[must_use]
    pub const fn strand_count(self) -> u64 {
        self.strand_count
    }

    /// Returns the balanced partition block count (`k + 1`).
    #[must_use]
    pub const fn block_count(self) -> u64 {
        self.block_count
    }

    /// Returns the number of N-free blocks emitted on each strand.
    #[must_use]
    pub const fn emitted_blocks(self) -> u64 {
        self.emitted_blocks
    }

    /// Returns the number of blocks omitted because they contain N.
    #[must_use]
    pub const fn omitted_unknown_blocks(self) -> u64 {
        self.omitted_unknown_blocks
    }

    /// Returns the total number of query N bases.
    #[must_use]
    pub const fn unknown_bases(self) -> u64 {
        self.unknown_bases
    }

    /// Returns the complete fixed request count across all strands.
    #[must_use]
    pub const fn request_count(self) -> u64 {
        self.request_count
    }

    /// Returns aggregate requested seed bases, counting each strand.
    #[must_use]
    pub const fn total_seed_bases(self) -> u64 {
        self.total_seed_bases
    }

    /// Returns the shortest balanced block before N omission.
    #[must_use]
    pub const fn minimum_block_bases(self) -> u64 {
        self.minimum_block_bases
    }

    /// Returns the longest balanced block before N omission.
    #[must_use]
    pub const fn maximum_block_bases(self) -> u64 {
        self.maximum_block_bases
    }

    /// Returns the Level 3 candidate-diagonal displacement obligation.
    #[must_use]
    pub const fn maximum_diagonal_displacement(self) -> u64 {
        self.max_edit_distance.get()
    }
}

/// One certified proof basis owning the only query instance.
pub struct CertifiedSeedPlan {
    fixed: FixedSeedPlan,
    strands: AdmissibleStrands,
    certificate: SeedCompletenessCertificate,
}

impl CertifiedSeedPlan {
    /// Returns the complete fixed plan for Level 2C candidate generation.
    #[must_use]
    pub const fn fixed_plan(&self) -> &FixedSeedPlan {
        &self.fixed
    }

    /// Returns the canonical admitted strand set.
    #[must_use]
    pub const fn strands(&self) -> AdmissibleStrands {
        self.strands
    }

    /// Returns the complete seed-recall certificate.
    #[must_use]
    pub const fn certificate(&self) -> SeedCompletenessCertificate {
        self.certificate
    }
}

impl fmt::Debug for CertifiedSeedPlan {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("CertifiedSeedPlan")
            .field("strands", &self.strands)
            .field("certificate", &self.certificate)
            .field("fixed", &self.fixed)
            .finish()
    }
}

/// A valid query that cannot have a nonempty `k + 1` partition.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SeedlessFallback {
    query: NormalizedSequence,
    strands: AdmissibleStrands,
    max_edit_distance: EditDistance,
}

impl SeedlessFallback {
    /// Returns the unchanged normalized query.
    #[must_use]
    pub const fn query(&self) -> &NormalizedSequence {
        &self.query
    }

    /// Returns the canonical admitted strand set.
    #[must_use]
    pub const fn strands(&self) -> AdmissibleStrands {
        self.strands
    }

    /// Returns the requested edit budget.
    #[must_use]
    pub const fn max_edit_distance(&self) -> EditDistance {
        self.max_edit_distance
    }

    /// Returns all owned routing state without reparsing the query.
    #[must_use]
    pub fn into_parts(self) -> (NormalizedSequence, AdmissibleStrands, EditDistance) {
        (self.query, self.strands, self.max_edit_distance)
    }
}

/// Proof that query-N count alone exceeds the full-query edit budget.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct NoAlignmentWithinBudget {
    query: NormalizedSequence,
    strands: AdmissibleStrands,
    max_edit_distance: EditDistance,
    unknown_bases: u64,
}

impl NoAlignmentWithinBudget {
    /// Returns the unchanged normalized query.
    #[must_use]
    pub const fn query(&self) -> &NormalizedSequence {
        &self.query
    }

    /// Returns the canonical admitted strand set.
    #[must_use]
    pub const fn strands(&self) -> AdmissibleStrands {
        self.strands
    }

    /// Returns the requested full-query edit budget.
    #[must_use]
    pub const fn max_edit_distance(&self) -> EditDistance {
        self.max_edit_distance
    }

    /// Returns the query-N count, which is strictly greater than the budget.
    #[must_use]
    pub const fn unknown_bases(&self) -> u64 {
        self.unknown_bases
    }

    /// Returns all owned routing state without reparsing the query.
    #[must_use]
    pub fn into_parts(self) -> (NormalizedSequence, AdmissibleStrands, EditDistance, u64) {
        (
            self.query,
            self.strands,
            self.max_edit_distance,
            self.unknown_bases,
        )
    }
}

/// Complete semantic outcome of proof-seed scheduling.
#[derive(Debug)]
pub enum ProofSeedOutcome {
    /// A nonempty proof basis is available for complete Level 2C search.
    Certified(CertifiedSeedPlan),
    /// The query is too short relative to the budget; a seedless path is needed.
    SeedlessFallbackRequired(SeedlessFallback),
    /// Query-N count proves no full-query alignment can fit the budget.
    NoAlignmentWithinBudget(NoAlignmentWithinBudget),
}

/// Boundary of a balanced block that failed physical representation.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ProofSeedBoundary {
    /// Inclusive block start.
    Start,
    /// Exclusive block end.
    End,
}

/// Variable-sized storage owned by the scheduler.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ProofSeedAllocation {
    /// Temporary requests transferred into `FixedSeedPlan`.
    Requests,
}

/// Defensive invariant checked before publishing a certificate.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ProofSeedInvariant {
    /// A block did not begin at the preceding block end.
    ContiguousPartition,
    /// A block was empty despite `m > k`.
    NonemptyBlock,
    /// The final block did not end at query length.
    CompleteCoverage,
    /// Unknown-block classification did not cover every block.
    ClassifiedBlockCount,
    /// A feasible certified query produced no N-free block.
    EmittedBlockExists,
    /// Materialized request count disagreed with the preflight count.
    RequestCount,
    /// `FixedSeedPlan` metrics disagreed with the certificate.
    FixedPlanMetrics,
}

/// A checked proof-seed construction failure.
#[non_exhaustive]
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ProofSeedError {
    /// Counting query N bases overflowed.
    UnknownBaseCountOverflow {
        /// Count before adding the current N.
        accumulated: u64,
    },
    /// `k + 1` overflowed.
    BlockCountOverflow {
        /// Supplied edit budget.
        max_edit_distance: u64,
    },
    /// The proof block count exceeds its explicit limit.
    BlockLimitExceeded {
        /// Requested blocks.
        requested: u64,
        /// Configured maximum.
        maximum: u64,
    },
    /// Emitted block count times strand count overflowed.
    RequestCountOverflow {
        /// N-free blocks.
        emitted_blocks: u64,
        /// Admitted strands.
        strand_count: u64,
    },
    /// N-free seed bases times strand count overflowed.
    TotalSeedBasesOverflow {
        /// Aggregate N-free block bases before strand multiplication.
        seed_bases_per_strand: u64,
        /// Admitted strands.
        strand_count: u64,
    },
    /// The complete request count exceeds its configured plan limit.
    RequestLimitExceeded {
        /// Requested fixed seeds.
        requested: u64,
        /// Configured maximum.
        maximum: u64,
    },
    /// Aggregate requested seed bases exceed their configured plan limit.
    TotalSeedBasesLimitExceeded {
        /// Requested aggregate bases.
        requested: u64,
        /// Configured maximum.
        maximum: u64,
    },
    /// A logical block boundary cannot index physical query storage.
    BoundaryNotRepresentable {
        /// Zero-based partition block ordinal.
        block_ordinal: u64,
        /// Failed boundary.
        boundary: ProofSeedBoundary,
        /// Logical boundary value.
        value: u64,
    },
    /// Request storage cannot fit this architecture.
    AllocationSizeOverflow {
        /// Allocation site.
        allocation: ProofSeedAllocation,
        /// Requested elements.
        elements: u64,
        /// Element width.
        element_size: u64,
    },
    /// Fallible request reservation failed.
    AllocationFailed {
        /// Allocation site.
        allocation: ProofSeedAllocation,
        /// Requested elements.
        elements: u64,
    },
    /// A balanced block could not become a typed query interval.
    Coordinate {
        /// Zero-based partition block ordinal.
        block_ordinal: u64,
        /// Underlying coordinate error.
        source: CoordinateError,
    },
    /// Level 2C rejected the fully preflighted fixed plan.
    PlanConstruction {
        /// Underlying fixed-plan error.
        source: SeedPlanError,
    },
    /// A private postcondition failed.
    Invariant {
        /// Failed invariant.
        invariant: ProofSeedInvariant,
        /// Expected value.
        expected: u64,
        /// Observed value.
        observed: u64,
    },
}

impl fmt::Display for ProofSeedError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnknownBaseCountOverflow { accumulated } => {
                write!(formatter, "query N count {accumulated} plus one overflowed")
            }
            Self::BlockCountOverflow { max_edit_distance } => write!(
                formatter,
                "proof block count for edit budget {max_edit_distance} plus one overflowed"
            ),
            Self::BlockLimitExceeded { requested, maximum } => write!(
                formatter,
                "proof block count {requested} exceeds configured maximum {maximum}"
            ),
            Self::RequestCountOverflow {
                emitted_blocks,
                strand_count,
            } => write!(
                formatter,
                "emitted block count {emitted_blocks} times strand count {strand_count} overflowed"
            ),
            Self::TotalSeedBasesOverflow {
                seed_bases_per_strand,
                strand_count,
            } => write!(
                formatter,
                "seed bases per strand {seed_bases_per_strand} times strand count {strand_count} overflowed"
            ),
            Self::RequestLimitExceeded { requested, maximum } => write!(
                formatter,
                "proof seed request count {requested} exceeds configured maximum {maximum}"
            ),
            Self::TotalSeedBasesLimitExceeded { requested, maximum } => write!(
                formatter,
                "proof seed bases {requested} exceeds configured maximum {maximum}"
            ),
            Self::BoundaryNotRepresentable {
                block_ordinal,
                boundary,
                value,
            } => write!(
                formatter,
                "proof block {block_ordinal} {boundary:?} boundary {value} does not fit this architecture"
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
            Self::Coordinate {
                block_ordinal,
                source,
            } => write!(
                formatter,
                "proof block {block_ordinal} has invalid query coordinates: {source}"
            ),
            Self::PlanConstruction { source } => {
                write!(formatter, "proof fixed-plan construction failed: {source}")
            }
            Self::Invariant {
                invariant,
                expected,
                observed,
            } => write!(
                formatter,
                "proof seed invariant {invariant:?} expected {expected}, observed {observed}"
            ),
        }
    }
}

impl std::error::Error for ProofSeedError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Coordinate { source, .. } => Some(source),
            Self::PlanConstruction { source } => Some(source),
            _ => None,
        }
    }
}

/// Builds the complete reference proof-seed outcome.
///
/// The partition depends only on query content, explicit strands, edit budget,
/// and limits. It is intentionally independent of every reference instance.
///
/// # Errors
///
/// Returns [`ProofSeedError`] for checked arithmetic, resource limits,
/// architecture sizing, allocation failure, coordinate failure, downstream
/// fixed-plan construction failure, or a defensive invariant violation.
#[allow(clippy::too_many_lines)]
pub fn schedule_proof_seeds(
    query: NormalizedSequence,
    strands: AdmissibleStrands,
    max_edit_distance: EditDistance,
    limits: ProofSeedLimits,
) -> Result<ProofSeedOutcome, ProofSeedError> {
    let query_bases = query.len();
    let edit_budget = max_edit_distance.get();
    let unknown_bases = count_unknown_bases(&query)?;

    if unknown_bases > edit_budget {
        return Ok(ProofSeedOutcome::NoAlignmentWithinBudget(
            NoAlignmentWithinBudget {
                query,
                strands,
                max_edit_distance,
                unknown_bases,
            },
        ));
    }

    if query_bases <= edit_budget {
        return Ok(ProofSeedOutcome::SeedlessFallbackRequired(
            SeedlessFallback {
                query,
                strands,
                max_edit_distance,
            },
        ));
    }

    let block_count = edit_budget
        .checked_add(1)
        .ok_or(ProofSeedError::BlockCountOverflow {
            max_edit_distance: edit_budget,
        })?;
    if block_count > limits.max_blocks {
        return Err(ProofSeedError::BlockLimitExceeded {
            requested: block_count,
            maximum: limits.max_blocks,
        });
    }

    let mut emitted_blocks = 0_u64;
    let mut omitted_unknown_blocks = 0_u64;
    let mut seed_bases_per_strand = 0_u64;
    let mut minimum_block_bases = u64::MAX;
    let mut maximum_block_bases = 0_u64;
    let mut previous_end = 0_u64;

    for block_ordinal in 0..block_count {
        let (start, end) = balanced_block_bounds(query_bases, block_count, block_ordinal)?;
        if start != previous_end {
            return Err(ProofSeedError::Invariant {
                invariant: ProofSeedInvariant::ContiguousPartition,
                expected: previous_end,
                observed: start,
            });
        }
        if start >= end {
            return Err(ProofSeedError::Invariant {
                invariant: ProofSeedInvariant::NonemptyBlock,
                expected: start.saturating_add(1),
                observed: end,
            });
        }
        let block_length = end - start;
        minimum_block_bases = minimum_block_bases.min(block_length);
        maximum_block_bases = maximum_block_bases.max(block_length);
        let block = query_block(&query, block_ordinal, start, end)?;
        if block
            .iter()
            .copied()
            .any(bsbit_core::alphabet::Base::is_unknown)
        {
            omitted_unknown_blocks =
                omitted_unknown_blocks
                    .checked_add(1)
                    .ok_or(ProofSeedError::Invariant {
                        invariant: ProofSeedInvariant::ClassifiedBlockCount,
                        expected: block_count,
                        observed: block_count.saturating_add(1),
                    })?;
        } else {
            emitted_blocks = emitted_blocks
                .checked_add(1)
                .ok_or(ProofSeedError::Invariant {
                    invariant: ProofSeedInvariant::ClassifiedBlockCount,
                    expected: block_count,
                    observed: block_count.saturating_add(1),
                })?;
            seed_bases_per_strand = seed_bases_per_strand.checked_add(block_length).ok_or(
                ProofSeedError::TotalSeedBasesOverflow {
                    seed_bases_per_strand,
                    strand_count: 1,
                },
            )?;
        }
        previous_end = end;
    }

    if previous_end != query_bases {
        return Err(ProofSeedError::Invariant {
            invariant: ProofSeedInvariant::CompleteCoverage,
            expected: query_bases,
            observed: previous_end,
        });
    }
    let classified_blocks =
        emitted_blocks
            .checked_add(omitted_unknown_blocks)
            .ok_or(ProofSeedError::Invariant {
                invariant: ProofSeedInvariant::ClassifiedBlockCount,
                expected: block_count,
                observed: u64::MAX,
            })?;
    if classified_blocks != block_count {
        return Err(ProofSeedError::Invariant {
            invariant: ProofSeedInvariant::ClassifiedBlockCount,
            expected: block_count,
            observed: classified_blocks,
        });
    }
    if emitted_blocks == 0 {
        return Err(ProofSeedError::Invariant {
            invariant: ProofSeedInvariant::EmittedBlockExists,
            expected: 1,
            observed: 0,
        });
    }

    let strand_count = strands.len();
    let request_count =
        emitted_blocks
            .checked_mul(strand_count)
            .ok_or(ProofSeedError::RequestCountOverflow {
                emitted_blocks,
                strand_count,
            })?;
    let total_seed_bases = seed_bases_per_strand.checked_mul(strand_count).ok_or(
        ProofSeedError::TotalSeedBasesOverflow {
            seed_bases_per_strand,
            strand_count,
        },
    )?;
    let plan_limits = limits.plan_limits;
    if request_count > plan_limits.max_requests() {
        return Err(ProofSeedError::RequestLimitExceeded {
            requested: request_count,
            maximum: plan_limits.max_requests(),
        });
    }
    if total_seed_bases > plan_limits.max_total_seed_bases() {
        return Err(ProofSeedError::TotalSeedBasesLimitExceeded {
            requested: total_seed_bases,
            maximum: plan_limits.max_total_seed_bases(),
        });
    }

    let request_storage = preflight_request_storage(request_count)?;
    let mut requests = Vec::new();
    requests
        .try_reserve_exact(request_storage)
        .map_err(|_| ProofSeedError::AllocationFailed {
            allocation: ProofSeedAllocation::Requests,
            elements: request_count,
        })?;
    let query_length = QueryLength::new(query_bases);
    for block_ordinal in 0..block_count {
        let (start, end) = balanced_block_bounds(query_bases, block_count, block_ordinal)?;
        let block = query_block(&query, block_ordinal, start, end)?;
        if block
            .iter()
            .copied()
            .any(bsbit_core::alphabet::Base::is_unknown)
        {
            continue;
        }
        let interval = QueryInterval::new(start, end, query_length).map_err(|source| {
            ProofSeedError::Coordinate {
                block_ordinal,
                source,
            }
        })?;
        for strand in strands.iter() {
            let materialized = u64::try_from(requests.len()).map_err(|_| {
                ProofSeedError::AllocationSizeOverflow {
                    allocation: ProofSeedAllocation::Requests,
                    elements: request_count,
                    element_size: element_size_u64::<FixedSeedRequest>(),
                }
            })?;
            if materialized >= request_count {
                return Err(ProofSeedError::Invariant {
                    invariant: ProofSeedInvariant::RequestCount,
                    expected: request_count,
                    observed: materialized.saturating_add(1),
                });
            }
            requests.push(FixedSeedRequest::new(strand, interval));
        }
    }

    let materialized =
        u64::try_from(requests.len()).map_err(|_| ProofSeedError::AllocationSizeOverflow {
            allocation: ProofSeedAllocation::Requests,
            elements: request_count,
            element_size: element_size_u64::<FixedSeedRequest>(),
        })?;
    if materialized != request_count {
        return Err(ProofSeedError::Invariant {
            invariant: ProofSeedInvariant::RequestCount,
            expected: request_count,
            observed: materialized,
        });
    }

    let fixed = FixedSeedPlan::new(query, &requests, plan_limits)
        .map_err(|source| ProofSeedError::PlanConstruction { source })?;
    let plan_metrics = fixed.metrics();
    if plan_metrics.query_bases() != query_bases {
        return Err(ProofSeedError::Invariant {
            invariant: ProofSeedInvariant::FixedPlanMetrics,
            expected: query_bases,
            observed: plan_metrics.query_bases(),
        });
    }
    if plan_metrics.request_count() != request_count {
        return Err(ProofSeedError::Invariant {
            invariant: ProofSeedInvariant::FixedPlanMetrics,
            expected: request_count,
            observed: plan_metrics.request_count(),
        });
    }
    if plan_metrics.total_seed_bases() != total_seed_bases {
        return Err(ProofSeedError::Invariant {
            invariant: ProofSeedInvariant::FixedPlanMetrics,
            expected: total_seed_bases,
            observed: plan_metrics.total_seed_bases(),
        });
    }

    let certificate = SeedCompletenessCertificate {
        query_bases,
        max_edit_distance,
        strand_count,
        block_count,
        emitted_blocks,
        omitted_unknown_blocks,
        unknown_bases,
        request_count,
        total_seed_bases,
        minimum_block_bases,
        maximum_block_bases,
    };
    Ok(ProofSeedOutcome::Certified(CertifiedSeedPlan {
        fixed,
        strands,
        certificate,
    }))
}

const fn strand_bit(strand: BisulfiteStrand) -> u8 {
    match strand {
        BisulfiteStrand::OT => 1,
        BisulfiteStrand::OB => 1 << 1,
        BisulfiteStrand::CTOT => 1 << 2,
        BisulfiteStrand::CTOB => 1 << 3,
    }
}

fn count_unknown_bases(query: &NormalizedSequence) -> Result<u64, ProofSeedError> {
    let mut count = 0_u64;
    for base in query.iter() {
        if base.is_unknown() {
            count = count
                .checked_add(1)
                .ok_or(ProofSeedError::UnknownBaseCountOverflow { accumulated: count })?;
        }
    }
    Ok(count)
}

fn balanced_block_bounds(
    query_bases: u64,
    block_count: u64,
    ordinal: u64,
) -> Result<(u64, u64), ProofSeedError> {
    let denominator = u128::from(block_count);
    let query = u128::from(query_bases);
    let start = (u128::from(ordinal) * query) / denominator;
    let next_ordinal = ordinal.checked_add(1).ok_or(ProofSeedError::Invariant {
        invariant: ProofSeedInvariant::CompleteCoverage,
        expected: block_count,
        observed: ordinal,
    })?;
    let end = (u128::from(next_ordinal) * query) / denominator;
    let start = u64::try_from(start).map_err(|_| ProofSeedError::BoundaryNotRepresentable {
        block_ordinal: ordinal,
        boundary: ProofSeedBoundary::Start,
        value: u64::MAX,
    })?;
    let end = u64::try_from(end).map_err(|_| ProofSeedError::BoundaryNotRepresentable {
        block_ordinal: ordinal,
        boundary: ProofSeedBoundary::End,
        value: u64::MAX,
    })?;
    Ok((start, end))
}

fn query_block(
    query: &NormalizedSequence,
    block_ordinal: u64,
    start: u64,
    end: u64,
) -> Result<&[bsbit_core::alphabet::Base], ProofSeedError> {
    let physical_start =
        usize::try_from(start).map_err(|_| ProofSeedError::BoundaryNotRepresentable {
            block_ordinal,
            boundary: ProofSeedBoundary::Start,
            value: start,
        })?;
    let physical_end =
        usize::try_from(end).map_err(|_| ProofSeedError::BoundaryNotRepresentable {
            block_ordinal,
            boundary: ProofSeedBoundary::End,
            value: end,
        })?;
    query
        .bases()
        .get(physical_start..physical_end)
        .ok_or(ProofSeedError::Invariant {
            invariant: ProofSeedInvariant::CompleteCoverage,
            expected: query.len(),
            observed: end,
        })
}

fn preflight_request_storage(elements: u64) -> Result<usize, ProofSeedError> {
    let element_size = element_size_u64::<FixedSeedRequest>();
    let storage =
        usize::try_from(elements).map_err(|_| ProofSeedError::AllocationSizeOverflow {
            allocation: ProofSeedAllocation::Requests,
            elements,
            element_size,
        })?;
    storage.checked_mul(size_of::<FixedSeedRequest>()).ok_or(
        ProofSeedError::AllocationSizeOverflow {
            allocation: ProofSeedAllocation::Requests,
            elements,
            element_size,
        },
    )?;
    Ok(storage)
}

fn element_size_u64<T>() -> u64 {
    u64::try_from(size_of::<T>()).unwrap_or(u64::MAX)
}

#[cfg(test)]
mod tests {
    use super::*;
    use bsbit_core::sequence::normalize_dna;

    fn normalized(raw: &[u8]) -> NormalizedSequence {
        normalize_dna(raw).expect("test input is normalized DNA")
    }

    fn all_strands() -> AdmissibleStrands {
        AdmissibleStrands::new(&BisulfiteStrand::ALL).expect("all strands are distinct")
    }

    fn certified(raw: &[u8], budget: u64) -> CertifiedSeedPlan {
        match schedule_proof_seeds(
            normalized(raw),
            all_strands(),
            EditDistance::new(budget),
            ProofSeedLimits::MAX,
        )
        .expect("schedule succeeds")
        {
            ProofSeedOutcome::Certified(plan) => plan,
            other => panic!("expected certified outcome, observed {other:?}"),
        }
    }

    #[test]
    fn admissible_strands_reject_empty_and_first_duplicate_and_iterate_canonically() {
        assert_eq!(
            AdmissibleStrands::new(&[]),
            Err(AdmissibleStrandError::Empty)
        );
        assert_eq!(
            AdmissibleStrands::new(&[
                BisulfiteStrand::CTOB,
                BisulfiteStrand::OT,
                BisulfiteStrand::CTOB,
                BisulfiteStrand::OT,
            ]),
            Err(AdmissibleStrandError::Duplicate {
                ordinal: 2,
                strand: BisulfiteStrand::CTOB,
            })
        );

        let set = AdmissibleStrands::new(&[
            BisulfiteStrand::CTOB,
            BisulfiteStrand::OT,
            BisulfiteStrand::CTOT,
        ])
        .expect("distinct set");
        assert_eq!(set.len(), 3);
        assert!(!set.is_empty());
        assert!(set.contains(BisulfiteStrand::OT));
        assert!(!set.contains(BisulfiteStrand::OB));
        assert_eq!(
            set.iter().collect::<Vec<_>>(),
            vec![
                BisulfiteStrand::OT,
                BisulfiteStrand::CTOT,
                BisulfiteStrand::CTOB,
            ]
        );
        assert_eq!(
            set.iter().rev().collect::<Vec<_>>(),
            vec![
                BisulfiteStrand::CTOB,
                BisulfiteStrand::CTOT,
                BisulfiteStrand::OT,
            ]
        );
    }

    #[test]
    fn exact_seventeen_has_one_full_block_per_strand() {
        let plan = certified(b"ACGTTGCAACGATTCGA", 0);
        let certificate = plan.certificate();
        assert_eq!(certificate.query_bases(), 17);
        assert_eq!(certificate.max_edit_distance(), EditDistance::new(0));
        assert_eq!(certificate.block_count(), 1);
        assert_eq!(certificate.emitted_blocks(), 1);
        assert_eq!(certificate.omitted_unknown_blocks(), 0);
        assert_eq!(certificate.unknown_bases(), 0);
        assert_eq!(certificate.strand_count(), 4);
        assert_eq!(certificate.request_count(), 4);
        assert_eq!(certificate.total_seed_bases(), 68);
        assert_eq!(certificate.minimum_block_bases(), 17);
        assert_eq!(certificate.maximum_block_bases(), 17);
        assert_eq!(certificate.maximum_diagonal_displacement(), 0);
        for request in plan.fixed_plan().requests() {
            assert_eq!(request.interval().start(), 0);
            assert_eq!(request.interval().end(), 17);
        }
    }

    #[test]
    fn balanced_blocks_cover_exactly_and_omit_unknown_blocks() {
        let strands = AdmissibleStrands::new(&[BisulfiteStrand::OB, BisulfiteStrand::OT])
            .expect("distinct strands");
        let outcome = schedule_proof_seeds(
            normalized(b"ACGNACGTAA"),
            strands,
            EditDistance::new(2),
            ProofSeedLimits::MAX,
        )
        .expect("schedule succeeds");
        let ProofSeedOutcome::Certified(plan) = outcome else {
            panic!("expected certified outcome");
        };
        let certificate = plan.certificate();
        assert_eq!(certificate.block_count(), 3);
        assert_eq!(certificate.emitted_blocks(), 2);
        assert_eq!(certificate.omitted_unknown_blocks(), 1);
        assert_eq!(certificate.unknown_bases(), 1);
        assert_eq!(certificate.minimum_block_bases(), 3);
        assert_eq!(certificate.maximum_block_bases(), 4);
        assert_eq!(certificate.request_count(), 4);
        assert_eq!(certificate.total_seed_bases(), 14);
        assert_eq!(
            plan.fixed_plan()
                .requests()
                .iter()
                .map(|request| (
                    request.strand(),
                    request.interval().start(),
                    request.interval().end()
                ))
                .collect::<Vec<_>>(),
            vec![
                (BisulfiteStrand::OT, 0, 3),
                (BisulfiteStrand::OT, 6, 10),
                (BisulfiteStrand::OB, 0, 3),
                (BisulfiteStrand::OB, 6, 10),
            ]
        );
    }

    #[test]
    fn seedless_and_no_alignment_outcomes_retain_exact_routing_state() {
        let strands =
            AdmissibleStrands::new(&[BisulfiteStrand::CTOT]).expect("one strand is nonempty");
        let fallback = schedule_proof_seeds(
            normalized(b"AC"),
            strands,
            EditDistance::new(2),
            ProofSeedLimits::new(0, SeedPlanLimits::new(0, 0)),
        )
        .expect("fallback ignores construction limits");
        let ProofSeedOutcome::SeedlessFallbackRequired(fallback) = fallback else {
            panic!("expected seedless fallback");
        };
        assert_eq!(fallback.query().to_ascii(), b"AC");
        assert_eq!(fallback.strands(), strands);
        assert_eq!(fallback.max_edit_distance(), EditDistance::new(2));
        let (query, returned_strands, budget) = fallback.into_parts();
        assert_eq!(query.to_ascii(), b"AC");
        assert_eq!(returned_strands, strands);
        assert_eq!(budget, EditDistance::new(2));

        let no_alignment = schedule_proof_seeds(
            normalized(b"ANNA"),
            strands,
            EditDistance::new(1),
            ProofSeedLimits::new(0, SeedPlanLimits::new(0, 0)),
        )
        .expect("N proof ignores construction limits");
        let ProofSeedOutcome::NoAlignmentWithinBudget(no_alignment) = no_alignment else {
            panic!("expected N-count proof");
        };
        assert_eq!(no_alignment.query().to_ascii(), b"ANNA");
        assert_eq!(no_alignment.strands(), strands);
        assert_eq!(no_alignment.max_edit_distance(), EditDistance::new(1));
        assert_eq!(no_alignment.unknown_bases(), 2);
        let (query, returned_strands, budget, unknown_bases) = no_alignment.into_parts();
        assert_eq!(query.to_ascii(), b"ANNA");
        assert_eq!(returned_strands, strands);
        assert_eq!(budget, EditDistance::new(1));
        assert_eq!(unknown_bases, 2);
    }

    #[test]
    fn construction_limits_fail_before_publication_with_exact_counts() {
        let query = normalized(b"ACGTACGTAA");
        let strands = all_strands();
        assert_eq!(
            schedule_proof_seeds(
                query.clone(),
                strands,
                EditDistance::new(2),
                ProofSeedLimits::new(2, SeedPlanLimits::MAX),
            )
            .expect_err("three blocks exceed two"),
            ProofSeedError::BlockLimitExceeded {
                requested: 3,
                maximum: 2,
            }
        );
        assert_eq!(
            schedule_proof_seeds(
                query.clone(),
                strands,
                EditDistance::new(2),
                ProofSeedLimits::new(3, SeedPlanLimits::new(11, u64::MAX)),
            )
            .expect_err("twelve requests exceed eleven"),
            ProofSeedError::RequestLimitExceeded {
                requested: 12,
                maximum: 11,
            }
        );
        assert_eq!(
            schedule_proof_seeds(
                query,
                strands,
                EditDistance::new(2),
                ProofSeedLimits::new(3, SeedPlanLimits::new(12, 39)),
            )
            .expect_err("forty seed bases exceed thirty-nine"),
            ProofSeedError::TotalSeedBasesLimitExceeded {
                requested: 40,
                maximum: 39,
            }
        );
    }

    #[test]
    fn error_displays_sources_and_private_preflight_are_stable() {
        let duplicate = AdmissibleStrandError::Duplicate {
            ordinal: 3,
            strand: BisulfiteStrand::OB,
        };
        assert_eq!(
            duplicate.to_string(),
            "admissible strand OB is duplicated at supplied ordinal 3"
        );
        assert!(std::error::Error::source(&duplicate).is_none());

        let block_limit = ProofSeedError::BlockLimitExceeded {
            requested: 4,
            maximum: 3,
        };
        assert_eq!(
            block_limit.to_string(),
            "proof block count 4 exceeds configured maximum 3"
        );
        assert!(std::error::Error::source(&block_limit).is_none());

        let source = CoordinateError::InvertedInterval {
            domain: bsbit_core::coordinate::CoordinateDomain::Query,
            start: 3,
            end: 2,
        };
        let coordinate = ProofSeedError::Coordinate {
            block_ordinal: 1,
            source,
        };
        assert!(std::error::Error::source(&coordinate).is_some());
        assert!(coordinate.to_string().contains("proof block 1"));

        if usize::BITS < u64::BITS {
            assert!(matches!(
                preflight_request_storage(u64::MAX),
                Err(ProofSeedError::AllocationSizeOverflow {
                    allocation: ProofSeedAllocation::Requests,
                    elements: u64::MAX,
                    ..
                })
            ));
        } else {
            let element_size = element_size_u64::<FixedSeedRequest>();
            let overflowing = u64::try_from(usize::MAX / size_of::<FixedSeedRequest>())
                .expect("usize fits u64")
                + 1;
            assert_eq!(
                preflight_request_storage(overflowing),
                Err(ProofSeedError::AllocationSizeOverflow {
                    allocation: ProofSeedAllocation::Requests,
                    elements: overflowing,
                    element_size,
                })
            );
        }
    }

    #[test]
    fn every_strand_error_and_limit_getter_has_stable_diagnostics() {
        let strand_errors = [
            AdmissibleStrandError::Empty,
            AdmissibleStrandError::SuppliedCountNotRepresentable { value: 7 },
            AdmissibleStrandError::SuppliedOrdinalNotRepresentable { value: 6 },
            AdmissibleStrandError::Duplicate {
                ordinal: 5,
                strand: BisulfiteStrand::CTOT,
            },
            AdmissibleStrandError::CountInvariant {
                expected_maximum: 4,
                observed: 5,
            },
        ];
        for error in strand_errors {
            assert!(!error.to_string().is_empty());
            assert!(std::error::Error::source(&error).is_none());
        }

        let limits = ProofSeedLimits::new(7, SeedPlanLimits::new(8, 99));
        assert_eq!(limits.max_blocks(), 7);
        assert_eq!(limits.plan_limits(), SeedPlanLimits::new(8, 99));
    }

    #[test]
    fn every_proof_error_variant_has_stable_diagnostics_and_source_policy() {
        let coordinate_source = CoordinateError::InvertedInterval {
            domain: bsbit_core::coordinate::CoordinateDomain::Query,
            start: 4,
            end: 3,
        };
        let plan_source = SeedPlanError::RequestLimitExceeded {
            requested: 2,
            maximum: 1,
        };
        let errors = [
            ProofSeedError::UnknownBaseCountOverflow {
                accumulated: u64::MAX,
            },
            ProofSeedError::BlockCountOverflow {
                max_edit_distance: u64::MAX,
            },
            ProofSeedError::BlockLimitExceeded {
                requested: 3,
                maximum: 2,
            },
            ProofSeedError::RequestCountOverflow {
                emitted_blocks: u64::MAX,
                strand_count: 4,
            },
            ProofSeedError::TotalSeedBasesOverflow {
                seed_bases_per_strand: u64::MAX,
                strand_count: 4,
            },
            ProofSeedError::RequestLimitExceeded {
                requested: 9,
                maximum: 8,
            },
            ProofSeedError::TotalSeedBasesLimitExceeded {
                requested: 101,
                maximum: 100,
            },
            ProofSeedError::BoundaryNotRepresentable {
                block_ordinal: 3,
                boundary: ProofSeedBoundary::End,
                value: u64::MAX,
            },
            ProofSeedError::AllocationSizeOverflow {
                allocation: ProofSeedAllocation::Requests,
                elements: u64::MAX,
                element_size: element_size_u64::<FixedSeedRequest>(),
            },
            ProofSeedError::AllocationFailed {
                allocation: ProofSeedAllocation::Requests,
                elements: 17,
            },
            ProofSeedError::Coordinate {
                block_ordinal: 2,
                source: coordinate_source,
            },
            ProofSeedError::PlanConstruction {
                source: plan_source,
            },
            ProofSeedError::Invariant {
                invariant: ProofSeedInvariant::FixedPlanMetrics,
                expected: 4,
                observed: 5,
            },
        ];
        for (ordinal, error) in errors.iter().enumerate() {
            assert!(
                !error.to_string().is_empty(),
                "error ordinal {ordinal} has empty diagnostics"
            );
            let source = std::error::Error::source(error);
            if matches!(
                error,
                ProofSeedError::Coordinate { .. } | ProofSeedError::PlanConstruction { .. }
            ) {
                assert!(source.is_some(), "error ordinal {ordinal} lost its source");
            } else {
                assert!(
                    source.is_none(),
                    "error ordinal {ordinal} invented a source"
                );
            }
        }
    }
}
