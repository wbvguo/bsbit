//! Safe scalar four-letter extension of owner-bound candidate windows.
//!
//! A candidate diagonal is search evidence, not a placement. This module
//! constructs the complete checked window implied by the Level 2D edit budget,
//! enumerates every nonempty admissible reference interval, and verifies each
//! with the Level 1 global bisulfite aligner.

use core::fmt;
use core::mem::size_of;

use crate::score::EditDistance;
use crate::search::candidate::{CandidateDiagonal, CandidateSet};
use crate::verification::distance::global_bs_alignment_banded_exact;
use crate::verification::distance::{DistanceError, DpCellLimit, global_bs_alignment};
use crate::verification::prefix_filter::ungapped_traceback_at_most_two_certified_cached_nm;
use crate::verification::prefix_filter::{
    BandedPrefixDistanceWorkspace, MIN_FILTER_QUERY_BASES, WordMyersQuery,
    ungapped_traceback_at_most_one,
};
use bsbit_core::bisulfite::{
    AlignmentOrientation, BisulfiteStrand, CytosineStrand, strand_semantics,
};
use bsbit_core::cigar::CoreCigar;
use bsbit_core::coordinate::{CoordinateShift, ReferenceInterval, ReferenceLength};
use bsbit_core::sequence::NormalizedSequence;
use bsbit_index::reference::{ContigId, ReferenceAccessError, ReferenceIndex};

/// Complete fail-whole limits for one candidate-window extension.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[allow(clippy::struct_field_names)]
pub struct ExtensionLimits {
    max_window_bases: u64,
    max_interval_alignments: u64,
    max_aggregate_dp_cells: u64,
    max_best_alignments: u64,
}

impl ExtensionLimits {
    /// Limits admitting every representable logical request.
    pub const MAX: Self = Self {
        max_window_bases: u64::MAX,
        max_interval_alignments: u64::MAX,
        max_aggregate_dp_cells: u64::MAX,
        max_best_alignments: u64::MAX,
    };

    /// Creates explicit window, interval, DP-cell, and retained-best limits.
    #[must_use]
    pub const fn new(
        max_window_bases: u64,
        max_interval_alignments: u64,
        max_aggregate_dp_cells: u64,
        max_best_alignments: u64,
    ) -> Self {
        Self {
            max_window_bases,
            max_interval_alignments,
            max_aggregate_dp_cells,
            max_best_alignments,
        }
    }

    /// Returns the largest admitted candidate-window length.
    #[must_use]
    pub const fn max_window_bases(self) -> u64 {
        self.max_window_bases
    }

    /// Returns the maximum complete interval-evaluation count.
    #[must_use]
    pub const fn max_interval_alignments(self) -> u64 {
        self.max_interval_alignments
    }

    /// Returns the maximum aggregate logical DP-cell count.
    #[must_use]
    pub const fn max_aggregate_dp_cells(self) -> u64 {
        self.max_aggregate_dp_cells
    }

    /// Returns the maximum retained equal-best alignment count.
    #[must_use]
    pub const fn max_best_alignments(self) -> u64 {
        self.max_best_alignments
    }
}

impl Default for ExtensionLimits {
    fn default() -> Self {
        Self::MAX
    }
}

/// A logical quantity measured by the scalar extension backend.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ExtensionCounter {
    /// Candidate anchors in one owner-bound set.
    CandidateAnchors,
    /// Candidate-window bases after contig clipping.
    WindowBases,
    /// Nonempty reference intervals submitted to global alignment.
    IntervalAlignments,
    /// Sum of logical DP cells across all submitted intervals.
    AggregateDpCells,
    /// Endpoint-distance scratch values for one direct sweep.
    EndpointDistances,
    /// Direct multi-endpoint distance sweeps.
    DistanceSweeps,
    /// Scalar DP cells or Myers text-symbol updates used by distance sweeps.
    DistanceFilterUpdates,
    /// Intervals submitted to canonical traceback after distance filtering.
    TracebackAlignments,
    /// In-budget alignments observed before local-best selection.
    PassingAlignments,
    /// Equal-best alignments retained in the complete result.
    BestAlignments,
}

/// A coordinate boundary that could not address physical sequence storage.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ExtensionBoundary {
    /// Inclusive candidate-window start.
    WindowStart,
    /// Exclusive candidate-window end.
    WindowEnd,
    /// Inclusive evaluated-interval start.
    IntervalStart,
    /// Exclusive evaluated-interval end.
    IntervalEnd,
}

/// A private consistency condition checked before publishing a result.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ExtensionInvariant {
    /// The closed-form interval count differed from enumeration.
    IntervalCount,
    /// The closed-form aggregate DP count differed from enumeration.
    AggregateDpCells,
    /// A retained alignment exceeded the caller's edit budget.
    BestDistanceWithinBudget,
    /// Equal-best output order was not strict interval order.
    ResultOrder,
    /// A direct distance sweep disagreed with canonical traceback.
    FilterDistance,
    #[cfg(test)]
    /// Candidate-start coverage events did not return to zero.
    CandidateStartBalance,
}

/// A checked candidate-window construction or extension failure.
#[non_exhaustive]
#[derive(Debug)]
pub enum ExtensionError {
    /// The candidate set belongs to another reference instance.
    ForeignCandidateSet,
    /// The requested candidate ordinal is outside the complete candidate set.
    CandidateOrdinalOutOfBounds {
        /// Supplied zero-based candidate ordinal.
        ordinal: u64,
        /// Complete candidate count.
        candidate_count: u64,
    },
    /// Full-query mapping does not accept an empty read.
    EmptyQuery,
    /// A logical metric exceeded the portable `u64` result domain.
    MetricNotRepresentable {
        /// Metric that could not be represented.
        counter: ExtensionCounter,
    },
    /// A complete logical metric exceeds its configured maximum.
    LimitExceeded {
        /// Limited metric.
        counter: ExtensionCounter,
        /// Complete requested amount.
        requested: u64,
        /// Configured maximum.
        maximum: u64,
    },
    /// A logical boundary cannot address this architecture's storage.
    BoundaryNotRepresentable {
        /// Failed boundary.
        boundary: ExtensionBoundary,
        /// Logical value.
        value: u64,
    },
    /// Resolving the owner-bound candidate contig failed.
    ReferenceAccess {
        /// Underlying owner or ordinal failure.
        source: ReferenceAccessError,
    },
    /// Exact Level 1 alignment failed for one interval.
    Alignment {
        /// Owner-bound contig ordinal.
        contig_ordinal: u64,
        /// Evaluated reference interval.
        interval: ReferenceInterval,
        /// Candidate strand.
        strand: BisulfiteStrand,
        /// Underlying scalar alignment failure.
        source: DistanceError,
    },
    /// A direct multi-endpoint distance sweep failed.
    DistanceSweep {
        /// Owner-bound contig ordinal.
        contig_ordinal: u64,
        /// Checked reference window.
        window: ReferenceInterval,
        /// Candidate strand.
        strand: BisulfiteStrand,
        /// Underlying scalar distance failure.
        source: DistanceError,
    },
    /// Retained-result storage could not be reserved.
    AllocationFailed {
        /// Allocation site.
        counter: ExtensionCounter,
        /// Requested retained element count.
        elements: u64,
    },
    /// A defensive postcondition failed.
    Invariant {
        /// Failed condition.
        invariant: ExtensionInvariant,
        /// Expected value.
        expected: u64,
        /// Observed value.
        observed: u64,
    },
}

impl fmt::Display for ExtensionError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ForeignCandidateSet => {
                formatter.write_str("candidate set belongs to another reference instance")
            }
            Self::CandidateOrdinalOutOfBounds {
                ordinal,
                candidate_count,
            } => write!(
                formatter,
                "candidate ordinal {ordinal} is outside candidate count {candidate_count}"
            ),
            Self::EmptyQuery => {
                formatter.write_str("single-read extension requires a nonempty query")
            }
            Self::MetricNotRepresentable { counter } => {
                write!(
                    formatter,
                    "extension metric {counter:?} is not representable as u64"
                )
            }
            Self::LimitExceeded {
                counter,
                requested,
                maximum,
            } => write!(
                formatter,
                "extension metric {counter:?} requested {requested}, exceeding maximum {maximum}"
            ),
            Self::BoundaryNotRepresentable { boundary, value } => write!(
                formatter,
                "extension {boundary:?} boundary {value} does not fit this architecture"
            ),
            Self::ReferenceAccess { source } => {
                write!(formatter, "candidate contig resolution failed: {source}")
            }
            Self::Alignment {
                contig_ordinal,
                interval,
                strand,
                source,
            } => write!(
                formatter,
                "scalar alignment failed for contig {contig_ordinal} {interval} strand {strand}: {source}"
            ),
            Self::DistanceSweep {
                contig_ordinal,
                window,
                strand,
                source,
            } => write!(
                formatter,
                "direct distance sweep failed for contig {contig_ordinal} {window} strand {strand}: {source}"
            ),
            Self::AllocationFailed { counter, elements } => write!(
                formatter,
                "failed to reserve {elements} retained elements for {counter:?}"
            ),
            Self::Invariant {
                invariant,
                expected,
                observed,
            } => write!(
                formatter,
                "extension invariant {invariant:?} expected {expected}, observed {observed}"
            ),
        }
    }
}

impl std::error::Error for ExtensionError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::ReferenceAccess { source } => Some(source),
            Self::Alignment { source, .. } | Self::DistanceSweep { source, .. } => Some(source),
            _ => None,
        }
    }
}

/// Exact preflight and result dimensions for one candidate window.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ExtensionMetrics {
    window_bases: u64,
    interval_alignments: u64,
    aggregate_dp_cells: u64,
    distance_sweeps: u64,
    distance_filter_updates: u64,
    traceback_alignments: u64,
    passing_alignments: u64,
    best_alignments: u64,
}

impl ExtensionMetrics {
    /// Returns candidate-window bases after clipping.
    #[must_use]
    pub const fn window_bases(self) -> u64 {
        self.window_bases
    }

    /// Returns the complete interval-evaluation count.
    #[must_use]
    pub const fn interval_alignments(self) -> u64 {
        self.interval_alignments
    }

    /// Returns the aggregate logical DP-cell count.
    #[must_use]
    pub const fn aggregate_dp_cells(self) -> u64 {
        self.aggregate_dp_cells
    }

    /// Returns direct multi-endpoint distance sweeps.
    #[must_use]
    pub const fn distance_sweeps(self) -> u64 {
        self.distance_sweeps
    }

    /// Returns scalar DP cells or Myers text-symbol updates used by sweeps.
    #[must_use]
    pub const fn distance_filter_updates(self) -> u64 {
        self.distance_filter_updates
    }

    /// Returns intervals submitted to canonical traceback after filtering.
    #[must_use]
    pub const fn traceback_alignments(self) -> u64 {
        self.traceback_alignments
    }

    /// Returns the count of evaluated intervals at or below the edit budget.
    #[must_use]
    pub const fn passing_alignments(self) -> u64 {
        self.passing_alignments
    }

    /// Returns the retained equal-best alignment count.
    #[must_use]
    pub const fn best_alignments(self) -> u64 {
        self.best_alignments
    }
}

/// One four-letter-verified full-query alignment.
#[derive(Clone)]
pub struct VerifiedAlignment {
    contig: ContigId,
    interval: ReferenceInterval,
    strand: BisulfiteStrand,
    orientation: AlignmentOrientation,
    cytosine_strand: CytosineStrand,
    distance: EditDistance,
    cigar: CoreCigar,
    cached_literal_nm: Option<u64>,
    multiple_optimal_paths: bool,
    maximum_seed_support: u64,
}

impl VerifiedAlignment {
    /// Returns the exact owner-bound contig.
    #[must_use]
    pub const fn contig(&self) -> &ContigId {
        &self.contig
    }

    /// Returns the verified zero-based half-open nonempty interval.
    #[must_use]
    pub const fn interval(&self) -> ReferenceInterval {
        self.interval
    }

    /// Returns the molecular/derived bisulfite strand.
    #[must_use]
    pub const fn strand(&self) -> BisulfiteStrand {
        self.strand
    }

    /// Returns query orientation relative to forward FASTA coordinates.
    #[must_use]
    pub const fn orientation(&self) -> AlignmentOrientation {
        self.orientation
    }

    /// Returns the genomic cytosine evidence strand used for verification.
    #[must_use]
    pub const fn cytosine_strand(&self) -> CytosineStrand {
        self.cytosine_strand
    }

    /// Returns the exact unit edit distance.
    #[must_use]
    pub const fn distance(&self) -> EditDistance {
        self.distance
    }

    /// Returns the canonical Level 1 core CIGAR.
    #[must_use]
    pub const fn cigar(&self) -> &CoreCigar {
        &self.cigar
    }

    /// Returns literal NM captured by an authoritative ungapped traceback.
    #[must_use]
    pub const fn cached_literal_nm(&self) -> Option<u64> {
        self.cached_literal_nm
    }

    /// Reports primary-distance path ambiguity within this exact interval.
    #[must_use]
    pub const fn multiple_optimal_paths(&self) -> bool {
        self.multiple_optimal_paths
    }

    /// Returns the candidate's distinct-request support for diagnostics.
    #[must_use]
    pub const fn maximum_seed_support(&self) -> u64 {
        self.maximum_seed_support
    }
}

impl fmt::Debug for VerifiedAlignment {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("VerifiedAlignment")
            .field("contig_ordinal", &self.contig.ordinal())
            .field("interval", &self.interval)
            .field("strand", &self.strand)
            .field("orientation", &self.orientation)
            .field("cytosine_strand", &self.cytosine_strand)
            .field("distance", &self.distance)
            .field("cigar", &self.cigar)
            .field("cached_literal_nm", &self.cached_literal_nm)
            .field("multiple_optimal_paths", &self.multiple_optimal_paths)
            .field("maximum_seed_support", &self.maximum_seed_support)
            .finish()
    }
}

/// Complete result for one owner-bound candidate window.
pub struct WindowExtension {
    window: ReferenceInterval,
    best_distance: Option<EditDistance>,
    alignments: Vec<VerifiedAlignment>,
    metrics: ExtensionMetrics,
}

struct RegionExtension {
    window: ReferenceInterval,
    best_distance: Option<EditDistance>,
    alignments: Vec<VerifiedAlignment>,
    metrics: ExtensionMetrics,
}

impl WindowExtension {
    /// Returns the checked contig-local candidate window.
    #[must_use]
    pub const fn window(&self) -> ReferenceInterval {
        self.window
    }

    /// Returns the minimum in-budget distance, or `None` when no interval passed.
    #[must_use]
    pub const fn best_distance(&self) -> Option<EditDistance> {
        self.best_distance
    }

    /// Returns every equal-best verified interval in canonical interval order.
    #[must_use]
    pub fn alignments(&self) -> &[VerifiedAlignment] {
        &self.alignments
    }

    /// Returns complete preflight and observed dimensions.
    #[must_use]
    pub const fn metrics(&self) -> ExtensionMetrics {
        self.metrics
    }

    /// Consumes this complete window and returns its verified alignments.
    #[must_use]
    pub fn into_alignments(self) -> Vec<VerifiedAlignment> {
        self.alignments
    }
}

impl fmt::Debug for WindowExtension {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("WindowExtension")
            .field("window", &self.window)
            .field("best_distance", &self.best_distance)
            .field("alignments", &self.alignments)
            .field("metrics", &self.metrics)
            .finish()
    }
}

/// Recovers the canonical qualified traceback for a retained paired-end placement.
pub(crate) fn traceback_retained_placement_banded(
    reference: &ReferenceIndex,
    raw_query: &NormalizedSequence,
    contig: &ContigId,
    interval: ReferenceInterval,
    strand: BisulfiteStrand,
    distance: EditDistance,
) -> Result<VerifiedAlignment, ExtensionError> {
    if raw_query.is_empty() {
        return Err(ExtensionError::EmptyQuery);
    }
    let resolved = reference
        .resolve_contig(contig)
        .map_err(|source| ExtensionError::ReferenceAccess { source })?;
    let start = boundary_to_storage(ExtensionBoundary::IntervalStart, interval.start())?;
    let end = boundary_to_storage(ExtensionBoundary::IntervalEnd, interval.end())?;
    let reference_interval =
        resolved
            .sequence()
            .bases()
            .get(start..end)
            .ok_or(ExtensionError::Invariant {
                invariant: ExtensionInvariant::ResultOrder,
                expected: resolved.sequence().len(),
                observed: interval.end(),
            })?;
    let semantics = strand_semantics(strand);
    let reversed;
    let oriented_query = match semantics.orientation() {
        AlignmentOrientation::Forward => raw_query,
        AlignmentOrientation::Reverse => {
            reversed = raw_query.reverse_complement();
            &reversed
        }
    };
    let cells = interval_dp_cells(interval.len(), oriented_query.len())?;
    let fast = ungapped_traceback_at_most_two_certified_cached_nm(
        reference_interval,
        oriented_query.bases(),
        semantics.cytosine_strand(),
    )
    .map_err(|source| ExtensionError::Alignment {
        contig_ordinal: contig.ordinal(),
        interval,
        strand,
        source,
    })?;
    let (traceback, cached_literal_nm) = if let Some((traceback, literal_nm)) = fast {
        (traceback, Some(literal_nm))
    } else {
        (
            global_bs_alignment_banded_exact(
                &NormalizedSequence::from_bases(reference_interval.iter().copied()),
                oriented_query,
                semantics.cytosine_strand(),
                distance,
                DpCellLimit::new(cells),
            )
            .map_err(|source| ExtensionError::Alignment {
                contig_ordinal: contig.ordinal(),
                interval,
                strand,
                source,
            })?,
            None,
        )
    };
    if traceback.distance() != distance {
        return Err(ExtensionError::Invariant {
            invariant: ExtensionInvariant::FilterDistance,
            expected: distance.get(),
            observed: traceback.distance().get(),
        });
    }
    Ok(VerifiedAlignment {
        contig: contig.clone(),
        interval,
        strand,
        orientation: semantics.orientation(),
        cytosine_strand: semantics.cytosine_strand(),
        distance,
        cigar: traceback.cigar().clone(),
        cached_literal_nm,
        multiple_optimal_paths: traceback.multiple_optimal_paths(),
        maximum_seed_support: 0,
    })
}

/// Extends one candidate selected by its owner-set ordinal.
///
/// The ordinal boundary prevents a borrowed anchor from one query from being
/// paired with another candidate set. Reference ownership is checked before
/// candidate ordinal or contig access.
///
/// # Errors
///
/// Returns [`ExtensionError`] for owner/ordinal/input, checked resource,
/// coordinate/storage, Level 1 alignment, allocation, or invariant failures.
pub fn extend_candidate_window(
    reference: &ReferenceIndex,
    candidates: &CandidateSet,
    candidate_ordinal: u64,
    max_edit_distance: EditDistance,
    limits: ExtensionLimits,
) -> Result<WindowExtension, ExtensionError> {
    let extension = extend_candidate_window_inner(
        reference,
        candidates,
        candidate_ordinal,
        max_edit_distance,
        limits,
    )?;
    Ok(WindowExtension {
        window: extension.window,
        best_distance: extension.best_distance,
        alignments: extension.alignments,
        metrics: extension.metrics,
    })
}

fn extend_candidate_window_inner(
    reference: &ReferenceIndex,
    candidates: &CandidateSet,
    candidate_ordinal: u64,
    max_edit_distance: EditDistance,
    limits: ExtensionLimits,
) -> Result<RegionExtension, ExtensionError> {
    if !candidates.belongs_to_reference(&reference.instance_id()) {
        return Err(ExtensionError::ForeignCandidateSet);
    }
    if candidates.query().is_empty() {
        return Err(ExtensionError::EmptyQuery);
    }
    let candidate_count = physical_to_logical(
        ExtensionCounter::CandidateAnchors,
        candidates.anchors().len(),
    )?;
    let storage = usize::try_from(candidate_ordinal).map_err(|_| {
        ExtensionError::CandidateOrdinalOutOfBounds {
            ordinal: candidate_ordinal,
            candidate_count,
        }
    })?;
    let Some(anchor) = candidates.anchors().get(storage) else {
        return Err(ExtensionError::CandidateOrdinalOutOfBounds {
            ordinal: candidate_ordinal,
            candidate_count,
        });
    };
    let contig = reference
        .resolve_contig(anchor.contig())
        .map_err(|source| ExtensionError::ReferenceAccess { source })?;
    let window = candidate_window(
        anchor.diagonal(),
        candidates.query().len(),
        max_edit_distance.get(),
        contig.sequence().len(),
    )?;
    extend_region(
        anchor.contig(),
        contig.sequence(),
        candidates.query(),
        anchor.strand(),
        window,
        max_edit_distance,
        limits,
        anchor.support().get(),
    )
}

fn candidate_window(
    diagonal: CandidateDiagonal,
    query_bases: u64,
    max_edit_distance: u64,
    contig_bases: u64,
) -> Result<ReferenceInterval, ExtensionError> {
    let signed_diagonal = match diagonal.shift() {
        CoordinateShift::Zero => 0_i128,
        CoordinateShift::Forward(value) => i128::from(value.get()),
        CoordinateShift::Backward(value) => -i128::from(value.get()),
    };
    let lower_unclipped = signed_diagonal - i128::from(max_edit_distance);
    let upper_unclipped = signed_diagonal + i128::from(query_bases) + i128::from(max_edit_distance);
    let contig = i128::from(contig_bases);
    let lower = u64::try_from(lower_unclipped.clamp(0, contig)).map_err(|_| {
        ExtensionError::MetricNotRepresentable {
            counter: ExtensionCounter::WindowBases,
        }
    })?;
    let upper = u64::try_from(upper_unclipped.clamp(0, contig)).map_err(|_| {
        ExtensionError::MetricNotRepresentable {
            counter: ExtensionCounter::WindowBases,
        }
    })?;
    let (start, end) = if lower <= upper {
        (lower, upper)
    } else {
        (lower, lower)
    };
    ReferenceInterval::new(start, end, ReferenceLength::new(contig_bases)).map_err(|_| {
        ExtensionError::Invariant {
            invariant: ExtensionInvariant::ResultOrder,
            expected: contig_bases,
            observed: end,
        }
    })
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn extend_region(
    contig: &ContigId,
    reference_sequence: &NormalizedSequence,
    raw_query: &NormalizedSequence,
    strand: BisulfiteStrand,
    window: ReferenceInterval,
    max_edit_distance: EditDistance,
    limits: ExtensionLimits,
    maximum_seed_support: u64,
) -> Result<RegionExtension, ExtensionError> {
    let semantics = strand_semantics(strand);
    let reversed;
    let oriented_query = match semantics.orientation() {
        AlignmentOrientation::Forward => raw_query,
        AlignmentOrientation::Reverse => {
            reversed = raw_query.reverse_complement();
            &reversed
        }
    };
    let preflight = preflight_region(window, oriented_query.len(), max_edit_distance, limits)?;
    let start_storage = boundary_to_storage(ExtensionBoundary::WindowStart, window.start())?;
    let end_storage = boundary_to_storage(ExtensionBoundary::WindowEnd, window.end())?;
    let window_bases = reference_sequence
        .bases()
        .get(start_storage..end_storage)
        .ok_or(ExtensionError::Invariant {
            invariant: ExtensionInvariant::ResultOrder,
            expected: reference_sequence.len(),
            observed: window.end(),
        })?;

    let mut alignments = Vec::new();
    let mut best_distance = None;
    let mut local_best_alignments = 0_u64;
    let mut passing_alignments = 0_u64;
    let mut observed_intervals = 0_u64;
    let mut observed_dp_cells = 0_u64;
    let mut distance_sweeps = 0_u64;
    let mut distance_filter_updates = 0_u64;
    let mut traceback_alignments = 0_u64;
    let minimum_length = minimum_interval_length(oriented_query.len(), max_edit_distance.get());
    let maximum_length = oriented_query
        .len()
        .saturating_add(max_edit_distance.get())
        .min(window.len());
    let has_intervals = preflight.interval_alignments > 0;
    let word_query = (has_intervals && oriented_query.len() >= MIN_FILTER_QUERY_BASES)
        .then(|| WordMyersQuery::new(oriented_query.bases(), semantics.cytosine_strand()))
        .flatten();
    let mut banded_workspace = if has_intervals && word_query.is_none() {
        Some(
            BandedPrefixDistanceWorkspace::new(oriented_query.len(), max_edit_distance).map_err(
                |source| ExtensionError::DistanceSweep {
                    contig_ordinal: contig.ordinal(),
                    window,
                    strand,
                    source,
                },
            )?,
        )
    } else {
        None
    };
    let maximum_length_storage =
        boundary_to_storage(ExtensionBoundary::IntervalEnd, maximum_length)?;
    let mut endpoint_distances = Vec::new();
    if maximum_length > 0 {
        preflight_vector_growth::<u64>(
            ExtensionCounter::EndpointDistances,
            &endpoint_distances,
            maximum_length,
        )?;
        endpoint_distances
            .try_reserve_exact(maximum_length_storage)
            .map_err(|_| ExtensionError::AllocationFailed {
                counter: ExtensionCounter::EndpointDistances,
                elements: maximum_length,
            })?;
        endpoint_distances.resize(maximum_length_storage, u64::MAX);
    }

    for local_start in 0..window.len() {
        let remaining = window.len() - local_start;
        if remaining < minimum_length {
            break;
        }
        let final_length = maximum_length.min(remaining);
        let local_start_storage =
            boundary_to_storage(ExtensionBoundary::IntervalStart, local_start)?;
        let final_length_storage =
            boundary_to_storage(ExtensionBoundary::IntervalEnd, final_length)?;
        let local_final_end_storage = local_start_storage
            .checked_add(final_length_storage)
            .ok_or(ExtensionError::MetricNotRepresentable {
                counter: ExtensionCounter::WindowBases,
            })?;
        let distances = &mut endpoint_distances[..final_length_storage];
        if let Some(query) = word_query.as_ref() {
            let reference_prefix = window_bases
                .get(local_start_storage..local_final_end_storage)
                .ok_or(ExtensionError::Invariant {
                    invariant: ExtensionInvariant::ResultOrder,
                    expected: window.len(),
                    observed: local_start.saturating_add(final_length),
                })?;
            if !query.prefix_distances(reference_prefix, distances) {
                return Err(ExtensionError::Invariant {
                    invariant: ExtensionInvariant::FilterDistance,
                    expected: final_length,
                    observed: physical_to_logical(
                        ExtensionCounter::EndpointDistances,
                        distances.len(),
                    )?,
                });
            }
            distance_filter_updates = distance_filter_updates.checked_add(final_length).ok_or(
                ExtensionError::MetricNotRepresentable {
                    counter: ExtensionCounter::DistanceFilterUpdates,
                },
            )?;
        } else {
            let reference_prefix = window_bases
                .get(local_start_storage..local_final_end_storage)
                .ok_or(ExtensionError::Invariant {
                    invariant: ExtensionInvariant::ResultOrder,
                    expected: window.len(),
                    observed: local_start.saturating_add(final_length),
                })?;
            let updates = banded_workspace
                .as_mut()
                .expect("non-Myers path owns banded workspace")
                .prefix_distances(
                    reference_prefix,
                    oriented_query.bases(),
                    semantics.cytosine_strand(),
                    distances,
                )
                .map_err(|source| ExtensionError::DistanceSweep {
                    contig_ordinal: contig.ordinal(),
                    window,
                    strand,
                    source,
                })?;
            distance_filter_updates = distance_filter_updates.checked_add(updates).ok_or(
                ExtensionError::MetricNotRepresentable {
                    counter: ExtensionCounter::DistanceFilterUpdates,
                },
            )?;
        }
        distance_sweeps =
            distance_sweeps
                .checked_add(1)
                .ok_or(ExtensionError::MetricNotRepresentable {
                    counter: ExtensionCounter::DistanceSweeps,
                })?;

        for interval_length in minimum_length..=final_length {
            let distance_index =
                boundary_to_storage(ExtensionBoundary::IntervalEnd, interval_length - 1)?;
            let filter_distance = endpoint_distances[distance_index];
            observed_intervals = observed_intervals.checked_add(1).ok_or(
                ExtensionError::MetricNotRepresentable {
                    counter: ExtensionCounter::IntervalAlignments,
                },
            )?;
            let cells = interval_dp_cells(interval_length, oriented_query.len())?;
            observed_dp_cells = observed_dp_cells.checked_add(cells).ok_or(
                ExtensionError::MetricNotRepresentable {
                    counter: ExtensionCounter::AggregateDpCells,
                },
            )?;
            let local_end = local_start.checked_add(interval_length).ok_or(
                ExtensionError::MetricNotRepresentable {
                    counter: ExtensionCounter::WindowBases,
                },
            )?;
            let local_end_storage = boundary_to_storage(ExtensionBoundary::IntervalEnd, local_end)?;
            let reference_interval = window_bases
                .get(local_start_storage..local_end_storage)
                .ok_or(ExtensionError::Invariant {
                    invariant: ExtensionInvariant::ResultOrder,
                    expected: window.len(),
                    observed: local_end,
                })?;
            let absolute_start = window.start().checked_add(local_start).ok_or(
                ExtensionError::MetricNotRepresentable {
                    counter: ExtensionCounter::WindowBases,
                },
            )?;
            let absolute_end = absolute_start.checked_add(interval_length).ok_or(
                ExtensionError::MetricNotRepresentable {
                    counter: ExtensionCounter::WindowBases,
                },
            )?;
            let interval = ReferenceInterval::new(
                absolute_start,
                absolute_end,
                ReferenceLength::new(reference_sequence.len()),
            )
            .map_err(|_| ExtensionError::Invariant {
                invariant: ExtensionInvariant::ResultOrder,
                expected: reference_sequence.len(),
                observed: absolute_end,
            })?;
            if filter_distance > max_edit_distance.get() {
                continue;
            }
            traceback_alignments = traceback_alignments.checked_add(1).ok_or(
                ExtensionError::MetricNotRepresentable {
                    counter: ExtensionCounter::TracebackAlignments,
                },
            )?;
            let fast_traceback = ungapped_traceback_at_most_one(
                reference_interval,
                oriented_query.bases(),
                semantics.cytosine_strand(),
            )
            .map_err(|source| ExtensionError::Alignment {
                contig_ordinal: contig.ordinal(),
                interval,
                strand,
                source,
            })?;
            let traceback = if let Some(traceback) = fast_traceback {
                traceback
            } else {
                let owned_reference =
                    NormalizedSequence::from_bases(reference_interval.iter().copied());
                global_bs_alignment(
                    &owned_reference,
                    oriented_query,
                    semantics.cytosine_strand(),
                    DpCellLimit::new(cells),
                )
                .map_err(|source| ExtensionError::Alignment {
                    contig_ordinal: contig.ordinal(),
                    interval,
                    strand,
                    source,
                })?
            };
            if traceback.distance().get() != filter_distance {
                return Err(ExtensionError::Invariant {
                    invariant: ExtensionInvariant::FilterDistance,
                    expected: filter_distance,
                    observed: traceback.distance().get(),
                });
            }
            if traceback.distance() > max_edit_distance {
                continue;
            }
            passing_alignments = passing_alignments.checked_add(1).ok_or(
                ExtensionError::MetricNotRepresentable {
                    counter: ExtensionCounter::PassingAlignments,
                },
            )?;
            let retain = match best_distance {
                None => {
                    best_distance = Some(traceback.distance());
                    local_best_alignments = 1;
                    alignments.clear();
                    true
                }
                Some(current) if traceback.distance() < current => {
                    best_distance = Some(traceback.distance());
                    local_best_alignments = 1;
                    alignments.clear();
                    true
                }
                Some(current) if traceback.distance() > current => false,
                Some(_) => {
                    local_best_alignments = local_best_alignments.checked_add(1).ok_or(
                        ExtensionError::MetricNotRepresentable {
                            counter: ExtensionCounter::BestAlignments,
                        },
                    )?;
                    true
                }
            };
            enforce_limit(
                ExtensionCounter::BestAlignments,
                local_best_alignments,
                limits.max_best_alignments,
            )?;
            if !retain {
                continue;
            }
            let retained_counter = ExtensionCounter::BestAlignments;
            let requested = physical_to_logical(retained_counter, alignments.len())?
                .checked_add(1)
                .ok_or(ExtensionError::MetricNotRepresentable {
                    counter: retained_counter,
                })?;
            preflight_vector_growth::<VerifiedAlignment>(retained_counter, &alignments, requested)?;
            alignments
                .try_reserve(1)
                .map_err(|_| ExtensionError::AllocationFailed {
                    counter: retained_counter,
                    elements: requested,
                })?;
            alignments.push(VerifiedAlignment {
                contig: contig.clone(),
                interval,
                strand,
                orientation: semantics.orientation(),
                cytosine_strand: semantics.cytosine_strand(),
                distance: traceback.distance(),
                cigar: traceback.cigar().clone(),
                cached_literal_nm: None,
                multiple_optimal_paths: traceback.multiple_optimal_paths(),
                maximum_seed_support,
            });
        }
    }

    if observed_intervals != preflight.interval_alignments {
        return Err(ExtensionError::Invariant {
            invariant: ExtensionInvariant::IntervalCount,
            expected: preflight.interval_alignments,
            observed: observed_intervals,
        });
    }
    if observed_dp_cells != preflight.aggregate_dp_cells {
        return Err(ExtensionError::Invariant {
            invariant: ExtensionInvariant::AggregateDpCells,
            expected: preflight.aggregate_dp_cells,
            observed: observed_dp_cells,
        });
    }
    validate_alignment_order(&alignments, max_edit_distance)?;
    Ok(RegionExtension {
        window,
        best_distance,
        alignments,
        metrics: ExtensionMetrics {
            window_bases: window.len(),
            interval_alignments: observed_intervals,
            aggregate_dp_cells: observed_dp_cells,
            distance_sweeps,
            distance_filter_updates,
            traceback_alignments,
            passing_alignments,
            best_alignments: local_best_alignments,
        },
    })
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct RegionPreflight {
    interval_alignments: u64,
    aggregate_dp_cells: u64,
}

fn preflight_region(
    window: ReferenceInterval,
    query_bases: u64,
    max_edit_distance: EditDistance,
    limits: ExtensionLimits,
) -> Result<RegionPreflight, ExtensionError> {
    enforce_limit(
        ExtensionCounter::WindowBases,
        window.len(),
        limits.max_window_bases,
    )?;
    let minimum_length = minimum_interval_length(query_bases, max_edit_distance.get());
    let maximum_length = query_bases
        .saturating_add(max_edit_distance.get())
        .min(window.len());
    if minimum_length > maximum_length || window.is_empty() {
        return Ok(RegionPreflight {
            interval_alignments: 0,
            aggregate_dp_cells: 0,
        });
    }
    let interval_alignments = interval_count(window.len(), minimum_length, maximum_length)?;
    enforce_limit(
        ExtensionCounter::IntervalAlignments,
        interval_alignments,
        limits.max_interval_alignments,
    )?;
    let aggregate_dp_cells =
        aggregate_dp_cells(window.len(), minimum_length, maximum_length, query_bases)?;
    enforce_limit(
        ExtensionCounter::AggregateDpCells,
        aggregate_dp_cells,
        limits.max_aggregate_dp_cells,
    )?;
    Ok(RegionPreflight {
        interval_alignments,
        aggregate_dp_cells,
    })
}

const fn minimum_interval_length(query_bases: u64, max_edit_distance: u64) -> u64 {
    let lower = query_bases.saturating_sub(max_edit_distance);
    if lower == 0 { 1 } else { lower }
}

fn interval_count(
    window_bases: u64,
    minimum_length: u64,
    maximum_length: u64,
) -> Result<u64, ExtensionError> {
    let count = u128::from(maximum_length - minimum_length + 1);
    let width_plus_one = u128::from(window_bases) + 1;
    let sum_lengths = range_sum(minimum_length, maximum_length)?;
    let total = count
        .checked_mul(width_plus_one)
        .and_then(|value| value.checked_sub(sum_lengths))
        .ok_or(ExtensionError::MetricNotRepresentable {
            counter: ExtensionCounter::IntervalAlignments,
        })?;
    u64::try_from(total).map_err(|_| ExtensionError::MetricNotRepresentable {
        counter: ExtensionCounter::IntervalAlignments,
    })
}

fn aggregate_dp_cells(
    window_bases: u64,
    minimum_length: u64,
    maximum_length: u64,
    query_bases: u64,
) -> Result<u64, ExtensionError> {
    let count = u128::from(maximum_length - minimum_length + 1);
    let sum_lengths = range_sum(minimum_length, maximum_length)?;
    let sum_lengths_plus_one =
        sum_lengths
            .checked_add(count)
            .ok_or(ExtensionError::MetricNotRepresentable {
                counter: ExtensionCounter::AggregateDpCells,
            })?;
    let weighted_positive = (u128::from(window_bases) + 1)
        .checked_mul(sum_lengths_plus_one)
        .ok_or(ExtensionError::MetricNotRepresentable {
            counter: ExtensionCounter::AggregateDpCells,
        })?;
    let sum_squares = range_square_sum(minimum_length, maximum_length)?;
    let weighted_negative =
        sum_squares
            .checked_add(sum_lengths)
            .ok_or(ExtensionError::MetricNotRepresentable {
                counter: ExtensionCounter::AggregateDpCells,
            })?;
    let interval_cells_without_query = weighted_positive.checked_sub(weighted_negative).ok_or(
        ExtensionError::MetricNotRepresentable {
            counter: ExtensionCounter::AggregateDpCells,
        },
    )?;
    let total = interval_cells_without_query
        .checked_mul(u128::from(query_bases) + 1)
        .ok_or(ExtensionError::MetricNotRepresentable {
            counter: ExtensionCounter::AggregateDpCells,
        })?;
    u64::try_from(total).map_err(|_| ExtensionError::MetricNotRepresentable {
        counter: ExtensionCounter::AggregateDpCells,
    })
}

fn range_sum(first: u64, last: u64) -> Result<u128, ExtensionError> {
    triangular(last)
        .checked_sub(triangular(first.saturating_sub(1)))
        .ok_or(ExtensionError::MetricNotRepresentable {
            counter: ExtensionCounter::IntervalAlignments,
        })
}

fn triangular(value: u64) -> u128 {
    let left = u128::from(value);
    let right = left + 1;
    if value.is_multiple_of(2) {
        (left / 2) * right
    } else {
        left * (right / 2)
    }
}

fn range_square_sum(first: u64, last: u64) -> Result<u128, ExtensionError> {
    square_sum(last)?
        .checked_sub(square_sum(first.saturating_sub(1))?)
        .ok_or(ExtensionError::MetricNotRepresentable {
            counter: ExtensionCounter::AggregateDpCells,
        })
}

fn square_sum(value: u64) -> Result<u128, ExtensionError> {
    if value == 0 {
        return Ok(0);
    }
    let mut factors = [
        u128::from(value),
        u128::from(value) + 1,
        u128::from(value) * 2 + 1,
    ];
    let mut divisor = 6_u128;
    for factor in &mut factors {
        let common = gcd(*factor, divisor);
        *factor /= common;
        divisor /= common;
    }
    debug_assert_eq!(divisor, 1);
    factors
        .into_iter()
        .try_fold(1_u128, u128::checked_mul)
        .ok_or(ExtensionError::MetricNotRepresentable {
            counter: ExtensionCounter::AggregateDpCells,
        })
}

const fn gcd(mut left: u128, mut right: u128) -> u128 {
    while right != 0 {
        let remainder = left % right;
        left = right;
        right = remainder;
    }
    left
}

fn interval_dp_cells(interval_bases: u64, query_bases: u64) -> Result<u64, ExtensionError> {
    interval_bases
        .checked_add(1)
        .and_then(|rows| {
            query_bases
                .checked_add(1)
                .and_then(|columns| rows.checked_mul(columns))
        })
        .ok_or(ExtensionError::MetricNotRepresentable {
            counter: ExtensionCounter::AggregateDpCells,
        })
}

const fn enforce_limit(
    counter: ExtensionCounter,
    requested: u64,
    maximum: u64,
) -> Result<(), ExtensionError> {
    if requested > maximum {
        Err(ExtensionError::LimitExceeded {
            counter,
            requested,
            maximum,
        })
    } else {
        Ok(())
    }
}

fn boundary_to_storage(boundary: ExtensionBoundary, value: u64) -> Result<usize, ExtensionError> {
    usize::try_from(value).map_err(|_| ExtensionError::BoundaryNotRepresentable { boundary, value })
}

fn physical_to_logical(counter: ExtensionCounter, value: usize) -> Result<u64, ExtensionError> {
    u64::try_from(value).map_err(|_| ExtensionError::MetricNotRepresentable { counter })
}

fn preflight_vector_growth<T>(
    counter: ExtensionCounter,
    values: &[T],
    requested_elements: u64,
) -> Result<(), ExtensionError> {
    let requested_storage = usize::try_from(requested_elements)
        .map_err(|_| ExtensionError::MetricNotRepresentable { counter })?;
    let _bytes = requested_storage
        .checked_mul(size_of::<T>())
        .ok_or(ExtensionError::MetricNotRepresentable { counter })?;
    if requested_storage <= values.len() {
        return Err(ExtensionError::Invariant {
            invariant: ExtensionInvariant::ResultOrder,
            expected: requested_elements,
            observed: physical_to_logical(counter, values.len())?,
        });
    }
    Ok(())
}

fn validate_alignment_order(
    alignments: &[VerifiedAlignment],
    max_edit_distance: EditDistance,
) -> Result<(), ExtensionError> {
    for alignment in alignments {
        if alignment.distance > max_edit_distance {
            return Err(ExtensionError::Invariant {
                invariant: ExtensionInvariant::BestDistanceWithinBudget,
                expected: max_edit_distance.get(),
                observed: alignment.distance.get(),
            });
        }
    }
    for pair in alignments.windows(2) {
        let left = pair[0].interval;
        let right = pair[1].interval;
        if (left.start(), left.end()) >= (right.start(), right.end()) {
            return Err(ExtensionError::Invariant {
                invariant: ExtensionInvariant::ResultOrder,
                expected: left.end(),
                observed: right.end(),
            });
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::search::candidate::{
        CandidateLimits, FixedSeedPlan, FixedSeedRequest, SeedPlanLimits,
        candidates_for_fixed_seeds,
    };
    use bsbit_core::coordinate::{QueryInterval, QueryLength};
    use bsbit_core::sequence::normalize_dna;
    use bsbit_index::reference::{ContigInput, ReferenceBuildLimits, ReferenceQueryLimits};

    fn sequence(input: &[u8]) -> NormalizedSequence {
        normalize_dna(input).expect("test sequence is valid")
    }

    fn one_candidate(
        reference_bytes: &[u8],
        query_bytes: &[u8],
        strand: BisulfiteStrand,
        seed_start: u64,
        seed_end: u64,
    ) -> (ReferenceIndex, CandidateSet) {
        let reference = ReferenceIndex::build(
            vec![ContigInput::new(b"chr".to_vec(), sequence(reference_bytes))],
            ReferenceBuildLimits::MAX,
        )
        .expect("reference builds");
        let query = sequence(query_bytes);
        let interval = QueryInterval::new(seed_start, seed_end, QueryLength::new(query.len()))
            .expect("seed is valid");
        let plan = FixedSeedPlan::new(
            query,
            &[FixedSeedRequest::new(strand, interval)],
            SeedPlanLimits::MAX,
        )
        .expect("plan builds");
        let candidates = candidates_for_fixed_seeds(
            &reference,
            &plan,
            ReferenceQueryLimits::MAX,
            CandidateLimits::MAX,
        )
        .expect("candidates build");
        (reference, candidates)
    }

    #[test]
    fn candidate_window_clips_signed_diagonals_without_unsigned_wrap() {
        let length = 20;
        let negative =
            CandidateDiagonal::before_contig(core::num::NonZeroU64::new(3).expect("nonzero"));
        assert_eq!(
            candidate_window(negative, 10, 2, length).expect("window"),
            ReferenceInterval::new(0, 9, ReferenceLength::new(length)).expect("interval")
        );
        let positive = CandidateDiagonal::at_or_after_contig(15);
        assert_eq!(
            candidate_window(positive, 10, 2, length).expect("window"),
            ReferenceInterval::new(13, 20, ReferenceLength::new(length)).expect("interval")
        );
    }

    #[test]
    fn closed_form_counts_equal_direct_enumeration() {
        for window in 0..=20_u64 {
            for query in 1..=12_u64 {
                for budget in 0..=query + 2 {
                    let minimum = minimum_interval_length(query, budget);
                    let maximum = query.saturating_add(budget).min(window);
                    if minimum > maximum || window == 0 {
                        continue;
                    }
                    let mut intervals = 0_u64;
                    let mut cells = 0_u64;
                    for start in 0..window {
                        for length in minimum..=maximum.min(window - start) {
                            intervals += 1;
                            cells += (length + 1) * (query + 1);
                        }
                    }
                    assert_eq!(
                        interval_count(window, minimum, maximum).expect("count is representable"),
                        intervals
                    );
                    assert_eq!(
                        aggregate_dp_cells(window, minimum, maximum, query)
                            .expect("cells are representable"),
                        cells
                    );
                }
            }
        }
    }

    #[test]
    fn exact_and_projection_only_false_equality_use_four_letter_extension() {
        let (reference, candidates) =
            one_candidate(b"GGGACCTGGG", b"ACCT", BisulfiteStrand::OT, 0, 4);
        let exact = extend_candidate_window(
            &reference,
            &candidates,
            0,
            EditDistance::new(0),
            ExtensionLimits::MAX,
        )
        .expect("extension succeeds");
        assert_eq!(exact.best_distance(), Some(EditDistance::new(0)));
        assert_eq!(exact.alignments().len(), 1);
        assert_eq!(exact.alignments()[0].interval().start(), 3);

        let (reference, candidates) =
            one_candidate(b"GGGATTTGGG", b"ACCC", BisulfiteStrand::OT, 0, 4);
        let rejected = extend_candidate_window(
            &reference,
            &candidates,
            0,
            EditDistance::new(0),
            ExtensionLimits::MAX,
        )
        .expect("extension succeeds");
        assert_eq!(rejected.best_distance(), None);
        assert!(rejected.alignments().is_empty());
    }

    #[test]
    fn limits_fail_whole_in_preflight_order() {
        let (reference, candidates) =
            one_candidate(b"AACCGGTT", b"CCGG", BisulfiteStrand::OT, 0, 4);
        let window_error = extend_candidate_window(
            &reference,
            &candidates,
            0,
            EditDistance::new(1),
            ExtensionLimits::new(5, u64::MAX, u64::MAX, u64::MAX),
        )
        .expect_err("window cap fails");
        assert!(matches!(
            window_error,
            ExtensionError::LimitExceeded {
                counter: ExtensionCounter::WindowBases,
                ..
            }
        ));

        let interval_error = extend_candidate_window(
            &reference,
            &candidates,
            0,
            EditDistance::new(1),
            ExtensionLimits::new(u64::MAX, 0, u64::MAX, u64::MAX),
        )
        .expect_err("interval cap fails");
        assert!(matches!(
            interval_error,
            ExtensionError::LimitExceeded {
                counter: ExtensionCounter::IntervalAlignments,
                ..
            }
        ));
    }

    #[test]
    fn foreign_owner_and_ordinal_errors_are_ordered() {
        let (reference, candidates) = one_candidate(b"AACCGG", b"CCG", BisulfiteStrand::OT, 0, 3);
        let foreign = ReferenceIndex::build(
            vec![ContigInput::new(b"chr".to_vec(), sequence(b"AACCGG"))],
            ReferenceBuildLimits::MAX,
        )
        .expect("reference");
        assert!(matches!(
            extend_candidate_window(
                &foreign,
                &candidates,
                u64::MAX,
                EditDistance::new(0),
                ExtensionLimits::MAX,
            ),
            Err(ExtensionError::ForeignCandidateSet)
        ));
        assert!(matches!(
            extend_candidate_window(
                &reference,
                &candidates,
                u64::MAX,
                EditDistance::new(0),
                ExtensionLimits::MAX,
            ),
            Err(ExtensionError::CandidateOrdinalOutOfBounds { .. })
        ));
    }
}
