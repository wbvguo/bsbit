//! Public construction, query, locate, resource, and diagnostic contracts.
//!
//! The concrete owner-bound reference index remains in the parent module.
//! Items are re-exported there so the public API path is unchanged.

use core::fmt;

use bsbit_core::bisulfite::BisulfiteStrand;
use bsbit_core::coordinate::CoordinateError;
use bsbit_core::reference::ReferenceSequenceMd5;
use bsbit_core::sequence::NormalizedSequence;

use crate::storage::fm::FmError;

use super::{CombinedIndexBackendError, validate_catalog_and_measure};

/// One owned contig supplied to reference construction.
#[derive(Clone, Debug)]
pub struct ContigInput {
    pub(crate) name: Vec<u8>,
    pub(crate) sequence: NormalizedSequence,
    pub(crate) md5: ReferenceSequenceMd5,
}

impl ContigInput {
    /// Creates one owned contig.
    #[must_use]
    pub fn new(name: Vec<u8>, sequence: NormalizedSequence) -> Self {
        let md5 = ReferenceSequenceMd5::from_normalized(sequence.bases());
        Self {
            name,
            sequence,
            md5,
        }
    }

    /// Returns the exact contig name bytes.
    #[must_use]
    pub fn name(&self) -> &[u8] {
        &self.name
    }

    /// Returns the retained normalized sequence.
    #[must_use]
    pub const fn sequence(&self) -> &NormalizedSequence {
        &self.sequence
    }

    /// Returns the standard SAM reference-sequence checksum.
    #[must_use]
    pub const fn md5(&self) -> ReferenceSequenceMd5 {
        self.md5
    }
}

/// Aggregate dimensions of a validated ordered reference catalog before any
/// search index is constructed.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ReferenceCatalogMetrics {
    pub(super) contig_count: u64,
    pub(super) total_name_bytes: u64,
    pub(super) total_reference_bases: u64,
}

impl ReferenceCatalogMetrics {
    /// Returns the number of ordered contigs.
    #[must_use]
    pub const fn contig_count(self) -> u64 {
        self.contig_count
    }

    /// Returns aggregate exact contig-name bytes.
    #[must_use]
    pub const fn total_name_bytes(self) -> u64 {
        self.total_name_bytes
    }

    /// Returns aggregate normalized bases including `N`.
    #[must_use]
    pub const fn total_reference_bases(self) -> u64 {
        self.total_reference_bases
    }
}

/// Limits for catalog validation that do not model or construct FM lanes.
#[allow(clippy::struct_field_names)]
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ReferenceCatalogLimits {
    pub(super) max_contigs: u64,
    pub(super) max_total_name_bytes: u64,
    pub(super) max_total_reference_bases: u64,
}

impl ReferenceCatalogLimits {
    /// Limits admitting every representable catalog dimension.
    pub const MAX: Self = Self {
        max_contigs: u64::MAX,
        max_total_name_bytes: u64::MAX,
        max_total_reference_bases: u64::MAX,
    };

    /// Sets the maximum ordered contig count.
    #[must_use]
    pub const fn with_max_contigs(mut self, value: u64) -> Self {
        self.max_contigs = value;
        self
    }

    /// Sets the maximum aggregate exact name bytes.
    #[must_use]
    pub const fn with_max_total_name_bytes(mut self, value: u64) -> Self {
        self.max_total_name_bytes = value;
        self
    }

    /// Sets the maximum aggregate normalized bases.
    #[must_use]
    pub const fn with_max_total_reference_bases(mut self, value: u64) -> Self {
        self.max_total_reference_bases = value;
        self
    }
}

impl Default for ReferenceCatalogLimits {
    fn default() -> Self {
        Self::MAX
    }
}

/// Validates ordered catalog semantics and dimensions without constructing FM
/// lanes or allocating a reference owner.
///
/// # Errors
///
/// Returns the same catalog-prefix validation errors and priority used by
/// [`ReferenceIndex::build`].
pub fn validate_reference_catalog(
    contigs: &[ContigInput],
    limits: ReferenceCatalogLimits,
) -> Result<ReferenceCatalogMetrics, ReferenceBuildError> {
    validate_catalog_and_measure(contigs, limits)
}

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
    pub(super) max_contigs: u64,
    pub(super) max_total_name_bytes: u64,
    pub(super) max_total_reference_bases: u64,
    pub(super) max_canonical_runs: u64,
    pub(super) max_suffix_rows_per_lane: u64,
    pub(super) max_lanes: u64,
    pub(super) max_projected_bases: u64,
    pub(super) max_projected_suffix_rows: u64,
    pub(super) max_estimated_retained_fm_bytes: u64,
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
    pub(super) max_pattern_bases: u64,
    pub(super) max_exact_hits: u64,
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
    pub(super) contig_count: u64,
    pub(super) total_name_bytes: u64,
    pub(super) total_reference_bases: u64,
    pub(super) canonical_bases: u64,
    pub(super) canonical_run_count: u64,
    pub(super) lane_count: u64,
    pub(super) projected_bases: u64,
    pub(super) projected_suffix_rows: u64,
    pub(super) estimated_retained_fm_bytes: u64,
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
    pub(super) located_coordinates: u64,
    pub(super) lf_steps: u64,
    pub(super) rank_operations: u64,
    pub(super) interval_nodes: u64,
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
