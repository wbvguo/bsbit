//! Owner-bound four-lane projected reference index.
//!
//! The implementation is a safe, deterministic correctness backend. Canonical
//! runs are separated by N and contig boundaries, and every run owns one exact
//! FM index for each bisulfite strand.
//!
//! Runtime owner handles deliberately do not implement generic value equality.
//!
//! ```compile_fail
//! use bsbit_index::reference::ReferenceInstanceId;
//!
//! fn requires_eq<T: Eq>() {}
//! requires_eq::<ReferenceInstanceId>();
//! ```
//!
//! Contig identifiers deliberately do not implement generic hashing.
//!
//! ```compile_fail
//! use bsbit_index::reference::ContigId;
//!
//! fn requires_hash<T: std::hash::Hash>() {}
//! requires_hash::<ContigId>();
//! ```
//!
//! Opaque matches deliberately do not implement generic value equality.
//!
//! ```compile_fail
//! use bsbit_index::reference::ProjectedMatches;
//!
//! fn requires_eq<T: Eq>() {}
//! requires_eq::<ProjectedMatches>();
//! ```

use core::fmt;

use bsbit_core::bisulfite::BisulfiteStrand;
use bsbit_core::coordinate::CoordinateError;

use crate::storage::fm::FmError;
#[cfg(feature = "combined-index")]
use crate::storage::fm::{FmInterval, ProjectedBase, SearchBase};

mod catalog;
pub use catalog::*;

mod runtime;
pub use runtime::*;
use runtime::{apply_limit, checked_build_add, physical_to_logical};

/// A logical aggregate resource controlled during reference construction.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ReferenceResource {
    /// Number of contigs.
    Contigs,
    /// Sum of contig-name bytes.
    TotalNameBytes,
    /// Sum of original reference bases, including N.
    TotalReferenceBases,
    /// Sum of canonical A, C, G, or T bases.
    CanonicalBases,
    /// Number of maximal canonical runs.
    CanonicalRuns,
    /// Largest suffix-row count for one run and lane.
    SuffixRowsPerLane,
    /// Aggregate number of run lanes.
    Lanes,
    /// Aggregate projected text bases.
    ProjectedBases,
    /// Aggregate projected suffix rows.
    ProjectedSuffixRows,
    /// Estimated retained FM bytes.
    EstimatedRetainedFmBytes,
}

/// A checked arithmetic operation used in diagnostics.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ReferenceArithmetic {
    /// Checked addition.
    Add,
    /// Checked multiplication.
    Multiply,
}

/// A Level 2B-owned allocation site.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ReferenceAllocation {
    /// Canonical-run metadata and lane handles.
    RunMetadata,
    /// Reusable lane-projection scratch.
    ProjectionScratch,
    /// One projected query pattern.
    ProjectedPattern,
    /// Opaque nonempty FM intervals.
    OpaqueMatches,
    /// Final recovered hits.
    FinalHits,
}

/// Explicit limits for one complete reference build.
#[allow(clippy::struct_field_names)]
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ReferenceBuildLimits {
    max_contigs: u64,
    max_total_name_bytes: u64,
    max_total_reference_bases: u64,
    max_canonical_runs: u64,
    max_suffix_rows_per_lane: u64,
    max_lanes: u64,
    max_projected_bases: u64,
    max_projected_suffix_rows: u64,
    max_estimated_retained_fm_bytes: u64,
}

impl ReferenceBuildLimits {
    /// Limits that admit every representable logical value.
    pub const MAX: Self = Self {
        max_contigs: u64::MAX,
        max_total_name_bytes: u64::MAX,
        max_total_reference_bases: u64::MAX,
        max_canonical_runs: u64::MAX,
        max_suffix_rows_per_lane: u64::MAX,
        max_lanes: u64::MAX,
        max_projected_bases: u64::MAX,
        max_projected_suffix_rows: u64::MAX,
        max_estimated_retained_fm_bytes: u64::MAX,
    };

    /// Sets the maximum contig count.
    #[must_use]
    pub const fn with_max_contigs(mut self, value: u64) -> Self {
        self.max_contigs = value;
        self
    }

    /// Sets the maximum total name bytes.
    #[must_use]
    pub const fn with_max_total_name_bytes(mut self, value: u64) -> Self {
        self.max_total_name_bytes = value;
        self
    }

    /// Sets the maximum original reference bases.
    #[must_use]
    pub const fn with_max_total_reference_bases(mut self, value: u64) -> Self {
        self.max_total_reference_bases = value;
        self
    }

    /// Sets the maximum canonical-run count.
    #[must_use]
    pub const fn with_max_canonical_runs(mut self, value: u64) -> Self {
        self.max_canonical_runs = value;
        self
    }

    /// Sets the maximum suffix rows in one run and lane.
    #[must_use]
    pub const fn with_max_suffix_rows_per_lane(mut self, value: u64) -> Self {
        self.max_suffix_rows_per_lane = value;
        self
    }

    /// Sets the maximum aggregate lane count.
    #[must_use]
    pub const fn with_max_lanes(mut self, value: u64) -> Self {
        self.max_lanes = value;
        self
    }

    /// Sets the maximum aggregate projected bases.
    #[must_use]
    pub const fn with_max_projected_bases(mut self, value: u64) -> Self {
        self.max_projected_bases = value;
        self
    }

    /// Sets the maximum aggregate projected suffix rows.
    #[must_use]
    pub const fn with_max_projected_suffix_rows(mut self, value: u64) -> Self {
        self.max_projected_suffix_rows = value;
        self
    }

    /// Sets the maximum estimated retained FM bytes.
    #[must_use]
    pub const fn with_max_estimated_retained_fm_bytes(mut self, value: u64) -> Self {
        self.max_estimated_retained_fm_bytes = value;
        self
    }
}

impl Default for ReferenceBuildLimits {
    fn default() -> Self {
        Self::MAX
    }
}

/// Explicit limits for one complete exact query.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ReferenceQueryLimits {
    max_pattern_bases: u64,
    max_exact_hits: u64,
}

impl ReferenceQueryLimits {
    /// Limits that admit every representable logical value.
    pub const MAX: Self = Self {
        max_pattern_bases: u64::MAX,
        max_exact_hits: u64::MAX,
    };

    /// Creates explicit pattern and exact-hit limits.
    #[must_use]
    pub const fn new(max_pattern_bases: u64, max_exact_hits: u64) -> Self {
        Self {
            max_pattern_bases,
            max_exact_hits,
        }
    }

    /// Sets the maximum exact-hit count.
    #[must_use]
    pub const fn with_max_exact_hits(mut self, value: u64) -> Self {
        self.max_exact_hits = value;
        self
    }

    /// Returns the maximum admitted exact-search pattern length.
    #[must_use]
    pub const fn max_pattern_bases(self) -> u64 {
        self.max_pattern_bases
    }

    /// Returns the maximum exact-hit count admitted by one complete query.
    #[must_use]
    pub const fn max_exact_hits(self) -> u64 {
        self.max_exact_hits
    }
}

impl Default for ReferenceQueryLimits {
    fn default() -> Self {
        Self::MAX
    }
}

/// Complete deterministic dimensions of a built projected reference.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ReferenceMetrics {
    contig_count: u64,
    total_name_bytes: u64,
    total_reference_bases: u64,
    canonical_bases: u64,
    canonical_run_count: u64,
    lane_count: u64,
    projected_bases: u64,
    projected_suffix_rows: u64,
    estimated_retained_fm_bytes: u64,
}

impl ReferenceMetrics {
    /// Returns the contig count.
    #[must_use]
    pub const fn contig_count(self) -> u64 {
        self.contig_count
    }

    /// Returns the sum of exact name bytes.
    #[must_use]
    pub const fn total_name_bytes(self) -> u64 {
        self.total_name_bytes
    }

    /// Returns the sum of original bases, including N.
    #[must_use]
    pub const fn total_reference_bases(self) -> u64 {
        self.total_reference_bases
    }

    /// Returns the number of canonical A, C, G, or T bases.
    #[must_use]
    pub const fn canonical_bases(self) -> u64 {
        self.canonical_bases
    }

    /// Returns the maximal canonical-run count.
    #[must_use]
    pub const fn canonical_run_count(self) -> u64 {
        self.canonical_run_count
    }

    /// Returns four times the canonical-run count.
    #[must_use]
    pub const fn lane_count(self) -> u64 {
        self.lane_count
    }

    /// Returns four times the canonical-base count.
    #[must_use]
    pub const fn projected_bases(self) -> u64 {
        self.projected_bases
    }

    /// Returns four times canonical bases plus runs.
    #[must_use]
    pub const fn projected_suffix_rows(self) -> u64 {
        self.projected_suffix_rows
    }

    /// Returns the checked retained-FM byte estimate.
    #[must_use]
    pub const fn estimated_retained_fm_bytes(self) -> u64 {
        self.estimated_retained_fm_bytes
    }
}

/// Physical work performed while locating one combined-index interval.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
#[doc(hidden)]
#[cfg(feature = "combined-index")]
pub(crate) struct PrivateCombinedLocateMetrics {
    located_rows: u64,
    lf_steps: u64,
    rank_operations: u64,
    interval_nodes: u64,
}

#[cfg(feature = "combined-index")]
impl PrivateCombinedLocateMetrics {
    /// Creates one checked backend report.
    #[must_use]
    pub(crate) const fn new(
        located_rows: u64,
        lf_steps: u64,
        rank_operations: u64,
        interval_nodes: u64,
    ) -> Self {
        Self {
            located_rows,
            lf_steps,
            rank_operations,
            interval_nodes,
        }
    }

    /// Returns emitted suffix rows.
    #[must_use]
    pub(crate) const fn located_rows(self) -> u64 {
        self.located_rows
    }

    /// Returns logical LF transitions represented by interval expansion.
    #[must_use]
    pub(crate) const fn lf_steps(self) -> u64 {
        self.lf_steps
    }

    /// Returns physical rank-boundary operations.
    #[must_use]
    pub(crate) const fn rank_operations(self) -> u64 {
        self.rank_operations
    }

    /// Returns processed shared interval-tree nodes.
    #[must_use]
    pub(crate) const fn interval_nodes(self) -> u64 {
        self.interval_nodes
    }
}

/// Query failure reported by the validated combined-index backend.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum CombinedIndexBackendError {
    /// A checked FM interval was malformed or outside this image.
    Interval,
    /// A previously validated image violated a query-time invariant.
    Structure,
}

impl fmt::Display for CombinedIndexBackendError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Interval => formatter.write_str("combined-index interval is invalid"),
            Self::Structure => formatter.write_str("combined-index query invariant failed"),
        }
    }
}

impl std::error::Error for CombinedIndexBackendError {}

/// Object-safe boundary for the frozen combined-direction index.
#[doc(hidden)]
#[cfg(feature = "combined-index")]
pub(crate) trait PrivateCombinedIndex: Send + Sync {
    /// Returns the length of one forward-reference half.
    fn reference_len(&self) -> u64;

    /// Returns the complete exact interval for one reversed projected pattern.
    fn exact_interval(
        &self,
        reversed_projected_pattern: &[SearchBase],
    ) -> Result<Option<FmInterval>, CombinedIndexBackendError>;

    /// Projected-digit counterpart for complete exact interval lookup.
    fn exact_projected_interval(
        &self,
        _reversed_projected_pattern: &[ProjectedBase],
    ) -> Result<Option<FmInterval>, CombinedIndexBackendError> {
        Err(CombinedIndexBackendError::Structure)
    }

    /// Looks up one caller-selected exact suffix in the dense table.
    fn lookup_interval(
        &self,
        _projected_suffix: &[SearchBase],
    ) -> Result<Option<FmInterval>, CombinedIndexBackendError> {
        Ok(None)
    }

    /// Projected-digit counterpart for dense exact lookup.
    fn lookup_projected_interval(
        &self,
        _projected_suffix: &[ProjectedBase],
    ) -> Result<Option<FmInterval>, CombinedIndexBackendError> {
        Err(CombinedIndexBackendError::Structure)
    }

    /// Resolves independent projected suffixes under caller-selected stopping rules.
    ///
    /// The backend owns only physical lookup and rank scheduling. The caller
    /// chooses the minimum accepted suffix length and the interval size at
    /// which further extension is unnecessary.
    fn resolve_projected_suffix_intervals(
        &self,
        patterns: &[&[ProjectedBase]],
        minimum_suffix_bases: usize,
        stop_interval_length: u64,
        output: &mut [Option<(FmInterval, u64)>],
    ) -> Result<(), CombinedIndexBackendError> {
        if patterns.len() != output.len() || stop_interval_length == 0 {
            return Err(CombinedIndexBackendError::Structure);
        }
        output.fill(None);
        for (destination, pattern) in output.iter_mut().zip(patterns) {
            if minimum_suffix_bases < COMBINED_EXACT_LOOKUP_BASES
                || pattern.len() < minimum_suffix_bases
            {
                continue;
            }
            let suffix_start = pattern.len() - COMBINED_EXACT_LOOKUP_BASES;
            let Some(mut interval) = self.lookup_projected_interval(&pattern[suffix_start..])?
            else {
                continue;
            };
            let mut matched_bases = COMBINED_EXACT_LOOKUP_BASES;
            let mut remaining_prefix_bases = suffix_start;
            while interval.len() > stop_interval_length && remaining_prefix_bases != 0 {
                let extended = self.backward_extend_projected_interval(
                    interval,
                    pattern[remaining_prefix_bases - 1],
                )?;
                if extended.is_empty() {
                    break;
                }
                interval = extended;
                matched_bases += 1;
                remaining_prefix_bases -= 1;
            }
            *destination = Some((
                interval,
                u64::try_from(matched_bases).map_err(|_| CombinedIndexBackendError::Structure)?,
            ));
        }
        Ok(())
    }

    /// Prepends one projected symbol to an owner-validated interval.
    fn backward_extend_interval(
        &self,
        interval: FmInterval,
        symbol: SearchBase,
    ) -> Result<FmInterval, CombinedIndexBackendError>;

    /// Projected-digit counterpart for one-symbol backward extension.
    fn backward_extend_projected_interval(
        &self,
        interval: FmInterval,
        symbol: ProjectedBase,
    ) -> Result<FmInterval, CombinedIndexBackendError>;

    /// Prepends one symbol to each independent interval in one backend round.
    fn backward_extend_intervals(
        &self,
        intervals: &[FmInterval],
        symbols: &[SearchBase],
        output: &mut [FmInterval],
    ) -> Result<(), CombinedIndexBackendError> {
        if intervals.len() != symbols.len() || intervals.len() != output.len() {
            return Err(CombinedIndexBackendError::Structure);
        }
        for ((output, &interval), &symbol) in output.iter_mut().zip(intervals).zip(symbols) {
            *output = self.backward_extend_interval(interval, symbol)?;
        }
        Ok(())
    }

    /// Projected-digit counterpart for one backend extension round.
    fn backward_extend_projected_intervals(
        &self,
        intervals: &[FmInterval],
        symbols: &[ProjectedBase],
        output: &mut [FmInterval],
    ) -> Result<(), CombinedIndexBackendError> {
        if intervals.len() != symbols.len() || intervals.len() != output.len() {
            return Err(CombinedIndexBackendError::Structure);
        }
        for ((output, &interval), &symbol) in output.iter_mut().zip(intervals).zip(symbols) {
            *output = self.backward_extend_projected_interval(interval, symbol)?;
        }
        Ok(())
    }

    /// Streams raw suffix coordinates from one checked interval.
    fn visit_interval(
        &self,
        interval: FmInterval,
        visitor: &mut dyn FnMut(u64) -> bool,
    ) -> Result<PrivateCombinedLocateMetrics, CombinedIndexBackendError>;

    /// Streams two complete intervals while overlapping locate state when possible.
    fn visit_intervals_two_lanes_complete(
        &self,
        intervals: [FmInterval; 2],
        visitor: &mut dyn FnMut(usize, u64),
    ) -> Result<[PrivateCombinedLocateMetrics; 2], CombinedIndexBackendError> {
        let first = self.visit_interval(intervals[0], &mut |position| {
            visitor(0, position);
            true
        })?;
        let second = self.visit_interval(intervals[1], &mut |position| {
            visitor(1, position);
            true
        })?;
        Ok([first, second])
    }
}

/// One optional combined index supplied to a combined reference owner.
#[doc(hidden)]
#[cfg(feature = "combined-index")]
pub(crate) struct PrivateCombinedReference {
    index: Box<dyn PrivateCombinedIndex>,
}

#[cfg(feature = "combined-index")]
impl PrivateCombinedReference {
    /// Erases a validated combined-index implementation behind the core trait.
    #[must_use]
    pub(crate) fn new(index: Box<dyn PrivateCombinedIndex>) -> Self {
        Self { index }
    }
}

/// Failure to assemble a combined index around one reference catalog.
#[derive(Clone, Debug, Eq, PartialEq)]
#[doc(hidden)]
#[non_exhaustive]
#[cfg(feature = "combined-index")]
pub(crate) enum PrivateCombinedReferenceError {
    /// The retained catalog failed normal deterministic validation.
    Catalog {
        /// Underlying catalog failure.
        source: ReferenceBuildError,
    },
    /// The combined index covers a different reference.
    CombinedDimensions {
        /// Total bases retained by the semantic reference catalog.
        expected_reference_len: u64,
        /// Forward-half length declared by the combined index.
        observed_reference_len: u64,
    },
}

#[cfg(feature = "combined-index")]
impl fmt::Display for PrivateCombinedReferenceError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Catalog { source } => write!(formatter, "combined reference catalog: {source}"),
            Self::CombinedDimensions {
                expected_reference_len,
                observed_reference_len,
            } => write!(
                formatter,
                "combined index expected reference length {expected_reference_len}, observed {observed_reference_len}"
            ),
        }
    }
}

#[cfg(feature = "combined-index")]
impl std::error::Error for PrivateCombinedReferenceError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Catalog { source } => Some(source),
            _ => None,
        }
    }
}

/// A structured construction failure.
#[non_exhaustive]
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ReferenceBuildError {
    /// No contigs were supplied.
    EmptyReference,
    /// A physical count cannot be represented in the logical domain.
    CountNotRepresentable {
        /// Resource whose physical count failed conversion.
        resource: ReferenceResource,
        /// Physical value.
        value: usize,
    },
    /// Checked logical arithmetic overflowed.
    ArithmeticOverflow {
        /// Resource being computed.
        resource: ReferenceResource,
        /// Arithmetic operation.
        operation: ReferenceArithmetic,
        /// Left operand.
        lhs: u64,
        /// Right operand.
        rhs: u64,
    },
    /// A configured limit rejected a complete build.
    LimitExceeded {
        /// Rejected resource.
        resource: ReferenceResource,
        /// Requested value.
        requested: u64,
        /// Configured maximum.
        maximum: u64,
    },
    /// One contig name is empty.
    EmptyContigName {
        /// Zero-based contig ordinal.
        contig_ordinal: u64,
    },
    /// A duplicate exact name was found.
    DuplicateContigName {
        /// Earliest prior exact-name ordinal.
        first_ordinal: u64,
        /// Smallest duplicate ordinal.
        duplicate_ordinal: u64,
    },
    /// One contig sequence is empty.
    EmptyContigSequence {
        /// Zero-based contig ordinal.
        contig_ordinal: u64,
    },
    /// The largest run exceeds the per-lane suffix-row limit.
    SuffixRowsPerLaneLimitExceeded {
        /// Requested suffix rows.
        requested: u64,
        /// Configured maximum.
        maximum: u64,
        /// Earliest contig containing the maximum.
        contig_ordinal: u64,
        /// Run start in the contig.
        run_start: u64,
    },
    /// A Level 2B-owned allocation cannot fit this architecture.
    AllocationSizeOverflow {
        /// Allocation site.
        allocation: ReferenceAllocation,
        /// Requested element count.
        elements: u64,
        /// Element width.
        element_size: u64,
    },
    /// A fallible Level 2B-owned reservation failed.
    AllocationFailed {
        /// Allocation site.
        allocation: ReferenceAllocation,
        /// Requested element count.
        elements: u64,
    },
    /// A private Level 2A lane build failed.
    FmBuild {
        /// Contig ordinal.
        contig_ordinal: u64,
        /// Canonical-run start.
        run_start: u64,
        /// Lane being built.
        strand: BisulfiteStrand,
        /// Underlying FM failure.
        source: FmError,
    },
    /// A checked internal build invariant failed.
    InternalInvariant {
        /// Expected value.
        expected: u64,
        /// Observed value.
        observed: u64,
    },
}

impl fmt::Display for ReferenceBuildError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyReference => formatter.write_str("reference contains no contigs"),
            Self::CountNotRepresentable { resource, value } => {
                write!(
                    formatter,
                    "{resource:?} count {value} is not representable as u64"
                )
            }
            Self::ArithmeticOverflow {
                resource,
                operation,
                lhs,
                rhs,
            } => write!(
                formatter,
                "{resource:?} arithmetic {lhs} {operation:?} {rhs} overflowed"
            ),
            Self::LimitExceeded {
                resource,
                requested,
                maximum,
            } => write!(
                formatter,
                "{resource:?} value {requested} exceeds configured maximum {maximum}"
            ),
            Self::EmptyContigName { contig_ordinal } => {
                write!(formatter, "contig {contig_ordinal} has an empty name")
            }
            Self::DuplicateContigName {
                first_ordinal,
                duplicate_ordinal,
            } => write!(
                formatter,
                "contig {duplicate_ordinal} duplicates the exact name of contig {first_ordinal}"
            ),
            Self::EmptyContigSequence { contig_ordinal } => {
                write!(formatter, "contig {contig_ordinal} has an empty sequence")
            }
            Self::SuffixRowsPerLaneLimitExceeded {
                requested,
                maximum,
                contig_ordinal,
                run_start,
            } => write!(
                formatter,
                "run at contig {contig_ordinal}:{run_start} needs {requested} suffix rows, exceeding {maximum}"
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
            Self::FmBuild {
                contig_ordinal,
                run_start,
                strand,
                source,
            } => write!(
                formatter,
                "FM build failed for contig {contig_ordinal} run {run_start} lane {strand:?}: {source}"
            ),
            Self::InternalInvariant { expected, observed } => write!(
                formatter,
                "reference build invariant expected {expected}, observed {observed}"
            ),
        }
    }
}

impl std::error::Error for ReferenceBuildError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::FmBuild { source, .. } => Some(source),
            _ => None,
        }
    }
}

/// An owner-bound contig access failure.
#[non_exhaustive]
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ReferenceAccessError {
    /// The contig identifier belongs to another index instance.
    ForeignContigId,
    /// The ordinal is outside this catalog.
    ContigOrdinalOutOfBounds {
        /// Requested ordinal.
        ordinal: u64,
        /// Number of contigs.
        contig_count: u64,
    },
}

impl fmt::Display for ReferenceAccessError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ForeignContigId => {
                formatter.write_str("contig identifier belongs to another reference instance")
            }
            Self::ContigOrdinalOutOfBounds {
                ordinal,
                contig_count,
            } => write!(
                formatter,
                "contig ordinal {ordinal} is outside catalog count {contig_count}"
            ),
        }
    }
}

impl std::error::Error for ReferenceAccessError {}

/// A query counter used in exact-search diagnostics.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ReferenceQueryCounter {
    /// Aggregate exact hits.
    ExactHits,
    /// Aggregate nonempty FM intervals.
    NonemptyIntervals,
    /// Physical rank-boundary operations performed by exact search.
    RankOperations,
}

/// A complete-search failure.
#[non_exhaustive]
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ReferenceQueryError {
    /// A physical pattern length cannot be represented by the logical width.
    PatternLengthNotRepresentable {
        /// Physical pattern length.
        pattern_len: usize,
    },
    /// The query pattern is empty.
    EmptyPattern,
    /// The pattern exceeds its configured limit.
    PatternLimitExceeded {
        /// Requested bases.
        requested: u64,
        /// Configured maximum.
        maximum: u64,
    },
    /// The normalized pattern contains N.
    UnsearchableBase {
        /// First zero-based N offset.
        offset: u64,
    },
    /// A query counter overflowed.
    CountOverflow {
        /// Counter being accumulated.
        counter: ReferenceQueryCounter,
        /// Accumulated count.
        accumulated: u64,
        /// Next increment.
        next: u64,
    },
    /// Exact hits exceed the configured complete-result limit.
    HitLimitExceeded {
        /// Requested complete hit count.
        requested: u64,
        /// Configured maximum.
        maximum: u64,
    },
    /// A query allocation cannot fit this architecture.
    AllocationSizeOverflow {
        /// Allocation site.
        allocation: ReferenceAllocation,
        /// Requested elements.
        elements: u64,
        /// Element width.
        element_size: u64,
    },
    /// A fallible query reservation failed.
    AllocationFailed {
        /// Allocation site.
        allocation: ReferenceAllocation,
        /// Requested elements.
        elements: u64,
    },
    /// The count and materialization passes disagreed.
    InvariantMismatch {
        /// Counter that disagreed.
        counter: ReferenceQueryCounter,
        /// First-pass value.
        expected: u64,
        /// Second-pass value.
        observed: u64,
    },
    /// Materialization would exceed the exact reserved entry count.
    CapacityInvariant {
        /// Exact reserved entry count.
        reserved: u64,
        /// Entries already materialized.
        materialized: u64,
    },
    /// The validated combined index rejected the query.
    CombinedIndex {
        /// Underlying combined-index query failure.
        source: CombinedIndexBackendError,
    },
}

impl fmt::Display for ReferenceQueryError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::PatternLengthNotRepresentable { pattern_len } => write!(
                formatter,
                "physical pattern length {pattern_len} is not representable as u64"
            ),
            Self::EmptyPattern => formatter.write_str("exact-search pattern is empty"),
            Self::PatternLimitExceeded { requested, maximum } => write!(
                formatter,
                "pattern length {requested} exceeds configured maximum {maximum}"
            ),
            Self::UnsearchableBase { offset } => {
                write!(
                    formatter,
                    "query contains unsearchable N at offset {offset}"
                )
            }
            Self::CountOverflow {
                counter,
                accumulated,
                next,
            } => write!(
                formatter,
                "{counter:?} count {accumulated} plus {next} overflowed"
            ),
            Self::HitLimitExceeded { requested, maximum } => write!(
                formatter,
                "exact hit count {requested} exceeds configured maximum {maximum}"
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
            Self::InvariantMismatch {
                counter,
                expected,
                observed,
            } => write!(
                formatter,
                "{counter:?} count pass produced {expected}, materialization produced {observed}"
            ),
            Self::CapacityInvariant {
                reserved,
                materialized,
            } => write!(
                formatter,
                "opaque-match reservation {reserved} cannot accept entry {materialized}"
            ),
            Self::CombinedIndex { source } => {
                write!(formatter, "combined-index search failed: {source}")
            }
        }
    }
}

impl std::error::Error for ReferenceQueryError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::CombinedIndex { source } => Some(source),
            _ => None,
        }
    }
}

/// A recovered-hit invariant that failed.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ReferenceLocateInvariant {
    /// A private run index was missing.
    MissingRun,
    /// FM locate returned a different count.
    OffsetCount,
    /// A terminal suffix appeared for a nonempty pattern.
    TerminalSuffix,
    /// A located interval exceeded its canonical run.
    RunBounds,
    /// Final hit materialization exceeded its exact reservation.
    FinalHitCapacity,
    /// Final recovered hit count differed from the search count.
    FinalHitCount,
    /// A physical locate counter overflowed.
    MetricOverflow,
}

/// Physical work performed while locating and recovering one match artifact.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
#[doc(hidden)]
pub struct ReferenceLocateMetrics {
    located_coordinates: u64,
    lf_steps: u64,
    rank_operations: u64,
    interval_nodes: u64,
}

impl ReferenceLocateMetrics {
    #[doc(hidden)]
    #[must_use]
    pub const fn located_coordinates(self) -> u64 {
        self.located_coordinates
    }

    #[doc(hidden)]
    #[must_use]
    pub const fn lf_steps(self) -> u64 {
        self.lf_steps
    }

    #[doc(hidden)]
    #[must_use]
    pub const fn rank_operations(self) -> u64 {
        self.rank_operations
    }

    #[doc(hidden)]
    #[must_use]
    pub const fn interval_nodes(self) -> u64 {
        self.interval_nodes
    }
}

/// A locate and coordinate-recovery failure.
#[non_exhaustive]
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ReferenceLocateError {
    /// The matches belong to another index instance.
    ForeignMatches,
    /// A final-hit allocation cannot fit this architecture.
    AllocationSizeOverflow {
        /// Allocation site.
        allocation: ReferenceAllocation,
        /// Requested elements.
        elements: u64,
        /// Element width.
        element_size: u64,
    },
    /// A fallible final-hit reservation failed.
    AllocationFailed {
        /// Allocation site.
        allocation: ReferenceAllocation,
        /// Requested elements.
        elements: u64,
    },
    /// A private FM locate operation failed.
    FmLocate {
        /// Contig ordinal.
        contig_ordinal: u64,
        /// Run start.
        run_start: u64,
        /// Lane.
        strand: BisulfiteStrand,
        /// Underlying FM failure.
        source: FmError,
    },
    /// The validated combined index rejected interval location.
    CombinedIndex {
        /// Underlying combined-index query failure.
        source: CombinedIndexBackendError,
    },
    /// Coordinate construction rejected a recovered interval.
    CoordinateRecovery {
        /// Contig ordinal.
        contig_ordinal: u64,
        /// Run start.
        run_start: u64,
        /// Lane.
        strand: BisulfiteStrand,
        /// Underlying coordinate failure.
        source: CoordinateError,
    },
    /// Checked coordinate arithmetic overflowed or underflowed.
    CoordinateArithmetic {
        /// Contig ordinal.
        contig_ordinal: u64,
        /// Run start.
        run_start: u64,
        /// Lane offset.
        offset: u64,
        /// Pattern length.
        pattern_len: u64,
    },
    /// A private trust-boundary invariant failed.
    Invariant {
        /// Invariant category.
        invariant: ReferenceLocateInvariant,
        /// Expected value.
        expected: u64,
        /// Observed value.
        observed: u64,
    },
}

impl fmt::Display for ReferenceLocateError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ForeignMatches => {
                formatter.write_str("projected matches belong to another reference instance")
            }
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
            Self::FmLocate {
                contig_ordinal,
                run_start,
                strand,
                source,
            } => write!(
                formatter,
                "FM locate failed for contig {contig_ordinal} run {run_start} lane {strand:?}: {source}"
            ),
            Self::CombinedIndex { source } => {
                write!(formatter, "combined-index locate failed: {source}")
            }
            Self::CoordinateRecovery {
                contig_ordinal,
                run_start,
                strand,
                source,
            } => write!(
                formatter,
                "coordinate recovery failed for contig {contig_ordinal} run {run_start} lane {strand:?}: {source}"
            ),
            Self::CoordinateArithmetic {
                contig_ordinal,
                run_start,
                offset,
                pattern_len,
            } => write!(
                formatter,
                "coordinate arithmetic failed at contig {contig_ordinal} run {run_start}, offset {offset}, pattern {pattern_len}"
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

impl std::error::Error for ReferenceLocateError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::FmLocate { source, .. } => Some(source),
            Self::CombinedIndex { source } => Some(source),
            Self::CoordinateRecovery { source, .. } => Some(source),
            _ => None,
        }
    }
}
