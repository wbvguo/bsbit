//! Exact scalar bisulfite-aware global distance and traceback.
//!
//! This module is the Level 1 semantic reference backend. Distance-only
//! execution keeps two rows. Traceback uses a checked full suffix table and
//! applies the complete canonical ordering instead of predecessor order.

use core::mem::size_of;

use core::fmt;

use crate::score::{EditDistance, EditDistanceOverflow};
use bsbit_core::bisulfite::{CytosineStrand, classify_bases};
use bsbit_core::cigar::{
    CigarError, CoreCigar, CoreCigarOp, CoreCigarRun, RawCigarRun, RawCoreCigar,
    canonicalize_operations, try_core_cigar, validate_cigar,
};
use bsbit_core::sequence::NormalizedSequence;

use crate::verification::cigar::{CigarEvaluationError, evaluate_cigar};

/// Maximum permitted logical DP cells, including row zero and column zero.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct DpCellLimit(u64);

impl DpCellLimit {
    /// A limit that admits every representable logical cell count.
    pub const MAX: Self = Self(u64::MAX);

    /// Constructs an explicit logical-cell limit. Zero rejects every request.
    #[must_use]
    pub const fn new(value: u64) -> Self {
        Self(value)
    }

    /// Returns the maximum logical cell count.
    #[must_use]
    pub const fn get(self) -> u64 {
        self.0
    }
}

/// The logical or physical allocation whose sizing could not be represented.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum MatrixAllocation {
    /// The reference prefix extent r + 1.
    LogicalReferenceExtent,
    /// The query prefix extent q + 1.
    LogicalQueryExtent,
    /// The logical product (r + 1) times (q + 1).
    LogicalCells,
    /// Two u64 rows used by distance computation.
    DistanceRows,
    /// Two u8 rows used for capped path multiplicity.
    PathCountRows,
    /// Three suffix scores for every logical DP cell.
    SuffixScores,
    /// Expanded traceback operations.
    TraceOperations,
    /// Worst-case canonical CIGAR runs.
    CigarRuns,
}

/// A secondary traceback field whose checked arithmetic overflowed.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum TraceScoreField {
    /// Total insertion and deletion bases.
    GapBases,
    /// Coalesced insertion and deletion runs.
    GapRuns,
}

/// An internal invariant checked at the public result boundary.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum AlignmentInvariant {
    /// No suffix transition existed at a nonterminal cell.
    MissingSuffixPath,
    /// Greedy canonical traceback could not reproduce the stored suffix score.
    MissingTraceStep,
    /// Prefix and suffix implementations disagreed on primary distance.
    PrimaryDistanceMismatch,
    /// A preceding exact filter supplied a distance different from recomputation.
    ExpectedDistanceMismatch,
    /// Independent CIGAR replay disagreed with the primary distance.
    CigarDistanceMismatch,
}

/// A checked scalar-DP, traceback, or result-validation failure.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum DistanceError {
    /// A logical dimension or physical allocation size is not representable.
    MatrixSizeOverflow {
        /// Logical reference length.
        reference_length: u64,
        /// Logical query length.
        query_length: u64,
        /// Allocation or logical quantity being sized.
        allocation: MatrixAllocation,
        /// Requested elements when that count was representable.
        requested_elements: Option<u64>,
        /// Requested bytes when that count was representable.
        requested_bytes: Option<u64>,
    },
    /// The exact request exceeds the caller's explicit logical-cell limit.
    ComputationLimitExceeded {
        /// Requested logical cells.
        requested_cells: u64,
        /// Caller-provided maximum cells.
        limit: u64,
    },
    /// Unit edit-cost addition overflowed.
    DistanceOverflow {
        /// Value before addition.
        accumulated: u64,
        /// Requested nonnegative increment.
        increment: u64,
    },
    /// Secondary canonical-path arithmetic overflowed.
    TraceScoreOverflow {
        /// Secondary field being accumulated.
        field: TraceScoreField,
        /// Value before addition.
        accumulated: u64,
        /// Requested increment.
        increment: u64,
    },
    /// A constructed CIGAR failed its independent validation boundary.
    CigarInvariant {
        /// Structured CIGAR failure.
        error: CigarError,
    },
    /// Replaying a constructed CIGAR against its sequences failed.
    CigarEvaluationInvariant {
        /// Structured replay failure.
        error: CigarEvaluationError,
    },
    /// Two independently computed alignment invariants disagreed.
    AlignmentInvariant {
        /// Invariant that failed.
        invariant: AlignmentInvariant,
        /// Reference prefix at the failure point.
        reference_index: u64,
        /// Query prefix at the failure point.
        query_index: u64,
        /// Expected value when applicable.
        expected: Option<u64>,
        /// Observed value when applicable.
        observed: Option<u64>,
    },
}

impl fmt::Display for DistanceError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MatrixSizeOverflow {
                reference_length,
                query_length,
                allocation,
                requested_elements,
                requested_bytes,
            } => write!(
                formatter,
                "cannot size {allocation:?} for reference/query lengths {reference_length}/{query_length}; elements={requested_elements:?}, bytes={requested_bytes:?}"
            ),
            Self::ComputationLimitExceeded {
                requested_cells,
                limit,
            } => write!(
                formatter,
                "requested {requested_cells} logical DP cells exceeds limit {limit}"
            ),
            Self::DistanceOverflow {
                accumulated,
                increment,
            } => write!(
                formatter,
                "edit distance addition {accumulated} + {increment} overflowed"
            ),
            Self::TraceScoreOverflow {
                field,
                accumulated,
                increment,
            } => write!(
                formatter,
                "trace score {field:?} addition {accumulated} + {increment} overflowed"
            ),
            Self::CigarInvariant { error } => {
                write!(formatter, "constructed CIGAR failed validation: {error}")
            }
            Self::CigarEvaluationInvariant { error } => {
                write!(formatter, "constructed CIGAR failed replay: {error}")
            }
            Self::AlignmentInvariant {
                invariant,
                reference_index,
                query_index,
                expected,
                observed,
            } => write!(
                formatter,
                "alignment invariant {invariant:?} failed at prefixes {reference_index}/{query_index}; expected={expected:?}, observed={observed:?}"
            ),
        }
    }
}

impl std::error::Error for DistanceError {}

impl From<EditDistanceOverflow> for DistanceError {
    fn from(error: EditDistanceOverflow) -> Self {
        Self::DistanceOverflow {
            accumulated: error.accumulated(),
            increment: error.increment(),
        }
    }
}

impl From<CigarError> for DistanceError {
    fn from(error: CigarError) -> Self {
        Self::CigarInvariant { error }
    }
}

impl From<CigarEvaluationError> for DistanceError {
    fn from(error: CigarEvaluationError) -> Self {
        Self::CigarEvaluationInvariant { error }
    }
}

/// Exact scalar traceback output.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TracebackResult {
    distance: EditDistance,
    cigar: CoreCigar,
    multiple_optimal_paths: bool,
}

impl TracebackResult {
    /// Returns the minimum primary unit edit distance.
    #[must_use]
    pub const fn distance(&self) -> EditDistance {
        self.distance
    }

    /// Returns the canonical Level 1 core CIGAR.
    #[must_use]
    pub const fn cigar(&self) -> &CoreCigar {
        &self.cigar
    }

    /// Reports whether at least two paths share the minimum primary distance.
    #[must_use]
    pub const fn multiple_optimal_paths(&self) -> bool {
        self.multiple_optimal_paths
    }

    pub(crate) fn ungapped(distance: EditDistance, length: u64) -> Result<Self, DistanceError> {
        debug_assert!(distance.get() <= length);
        let raw = if length == 0 {
            RawCoreCigar::default()
        } else {
            RawCoreCigar::new([RawCigarRun::new(CoreCigarOp::M, length)])
        };
        let cigar = try_core_cigar(&raw, length, length)?;
        Ok(Self {
            distance,
            cigar,
            multiple_optimal_paths: false,
        })
    }
}

/// Computes exact global bisulfite-aware unit edit distance.
///
/// Reference and oriented-query roles are never exchanged. Time is O(r*q) and
/// auxiliary storage is two query-length rows.
///
/// # Errors
///
/// Returns a structured size, limit, or arithmetic error before publishing any
/// partial result.
pub fn global_bs_distance(
    reference: &NormalizedSequence,
    oriented_query: &NormalizedSequence,
    cytosine_strand: CytosineStrand,
    limit: DpCellLimit,
) -> Result<EditDistance, DistanceError> {
    let dimensions = validate_dimensions(reference, oriented_query, limit)?;
    preflight_distance(&dimensions)?;
    distance_rows(reference, oriented_query, cytosine_strand, &dimensions)
}

/// Computes an exact canonical global alignment.
///
/// The primary objective is unit edit distance. Ties then minimize gap bases,
/// gap runs, and finally the expanded operation stream under D < I < M.
/// Ambiguity is computed from primary-distance paths before those secondary
/// choices. Time and suffix storage are O(r*q).
///
/// # Errors
///
/// Returns a structured size, limit, arithmetic, or invariant error. All
/// planned capacities are checked before the first allocation.
pub fn global_bs_alignment(
    reference: &NormalizedSequence,
    oriented_query: &NormalizedSequence,
    cytosine_strand: CytosineStrand,
    limit: DpCellLimit,
) -> Result<TracebackResult, DistanceError> {
    let dimensions = validate_dimensions(reference, oriented_query, limit)?;
    let capacities = preflight_alignment(&dimensions)?;

    let (primary_distance, path_count) = primary_distance_and_count(
        reference,
        oriented_query,
        cytosine_strand,
        &dimensions,
        capacities.row_elements,
    )?;
    let suffix_scores = build_suffix_scores(
        reference,
        oriented_query,
        cytosine_strand,
        &dimensions,
        capacities.suffix_elements,
    )?;
    let start_score = suffix_score(&suffix_scores, 0, 0, PREVIOUS_OTHER, &dimensions);
    if start_score.distance != primary_distance.get() {
        return Err(invariant_error(
            AlignmentInvariant::PrimaryDistanceMismatch,
            0,
            0,
            Some(primary_distance.get()),
            Some(start_score.distance),
        ));
    }

    let operations = canonical_trace(
        reference,
        oriented_query,
        cytosine_strand,
        &dimensions,
        &suffix_scores,
        capacities.trace_operations,
    )?;
    let cigar = canonicalize_operations(operations.iter().copied())?;
    validate_cigar(&cigar, dimensions.reference_length, dimensions.query_length)?;
    let evaluation = evaluate_cigar(&cigar, reference, oriented_query, cytosine_strand)?;
    if evaluation.distance() != primary_distance {
        return Err(invariant_error(
            AlignmentInvariant::CigarDistanceMismatch,
            dimensions.reference_length,
            dimensions.query_length,
            Some(primary_distance.get()),
            Some(evaluation.distance().get()),
        ));
    }

    Ok(TracebackResult {
        distance: primary_distance,
        cigar,
        multiple_optimal_paths: path_count > 1,
    })
}

/// Reconstructs the exact canonical alignment inside a proven edit-distance band.
///
/// The caller supplies an exact distance already established by an independent
/// filter. Every path having that primary distance stays within
/// `|reference_prefix - query_prefix| <= expected_distance`, because the
/// displacement at any prefix cannot exceed the total number of gap bases.
/// Restricting both the primary recurrence and suffix traceback to that band
/// therefore preserves the minimum distance, primary-path ambiguity, all
/// secondary objectives, and the final `D < I < M` lexicographic tie-break.
///
/// This is the qualified read-mapping backend. [`global_bs_alignment`] remains
/// the scalar semantic oracle, and both implementations retain the same
/// full-matrix resource preflight.
pub(crate) fn global_bs_alignment_banded_exact(
    reference: &NormalizedSequence,
    oriented_query: &NormalizedSequence,
    cytosine_strand: CytosineStrand,
    expected_distance: EditDistance,
    limit: DpCellLimit,
) -> Result<TracebackResult, DistanceError> {
    let dimensions = validate_dimensions(reference, oriented_query, limit)?;
    let capacities = preflight_alignment(&dimensions)?;
    let band = usize::try_from(expected_distance.get()).unwrap_or(usize::MAX);
    if dimensions
        .reference_storage
        .abs_diff(dimensions.query_storage)
        > band
    {
        return Err(invariant_error(
            AlignmentInvariant::ExpectedDistanceMismatch,
            dimensions.reference_length,
            dimensions.query_length,
            Some(expected_distance.get()),
            None,
        ));
    }

    let (primary_distance, path_count) = primary_distance_and_count_banded(
        reference,
        oriented_query,
        cytosine_strand,
        &dimensions,
        capacities.row_elements,
        band,
    )?;
    if primary_distance != expected_distance {
        return Err(invariant_error(
            AlignmentInvariant::ExpectedDistanceMismatch,
            dimensions.reference_length,
            dimensions.query_length,
            Some(expected_distance.get()),
            Some(primary_distance.get()),
        ));
    }

    let suffix_scores = build_suffix_scores_banded(
        reference,
        oriented_query,
        cytosine_strand,
        &dimensions,
        capacities.suffix_elements,
        band,
    )?;
    let start_score = suffix_scores.get(0, 0, PREVIOUS_OTHER);
    if start_score.distance != primary_distance.get() {
        return Err(invariant_error(
            AlignmentInvariant::PrimaryDistanceMismatch,
            0,
            0,
            Some(primary_distance.get()),
            Some(start_score.distance),
        ));
    }

    let operations = canonical_trace_banded(
        reference,
        oriented_query,
        cytosine_strand,
        &dimensions,
        &suffix_scores,
        capacities.trace_operations,
        band,
    )?;
    let cigar = canonicalize_operations(operations.iter().copied())?;
    validate_cigar(&cigar, dimensions.reference_length, dimensions.query_length)?;
    let evaluation = evaluate_cigar(&cigar, reference, oriented_query, cytosine_strand)?;
    if evaluation.distance() != primary_distance {
        return Err(invariant_error(
            AlignmentInvariant::CigarDistanceMismatch,
            dimensions.reference_length,
            dimensions.query_length,
            Some(primary_distance.get()),
            Some(evaluation.distance().get()),
        ));
    }

    Ok(TracebackResult {
        distance: primary_distance,
        cigar,
        multiple_optimal_paths: path_count > 1,
    })
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct Dimensions {
    reference_length: u64,
    query_length: u64,
    rows: u64,
    columns: u64,
    cells: u64,
    reference_storage: usize,
    query_storage: usize,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct AlignmentCapacities {
    row_elements: usize,
    suffix_elements: usize,
    trace_operations: usize,
}

fn validate_dimensions(
    reference: &NormalizedSequence,
    query: &NormalizedSequence,
    limit: DpCellLimit,
) -> Result<Dimensions, DistanceError> {
    let reference_length = reference.len();
    let query_length = query.len();
    let (rows, columns, cells) =
        validate_logical_dimensions(reference_length, query_length, limit)?;

    let reference_storage = usize::try_from(reference_length).map_err(|_| {
        matrix_error(
            reference_length,
            query_length,
            MatrixAllocation::LogicalReferenceExtent,
            Some(reference_length),
            None,
        )
    })?;
    let query_storage = usize::try_from(query_length).map_err(|_| {
        matrix_error(
            reference_length,
            query_length,
            MatrixAllocation::LogicalQueryExtent,
            Some(query_length),
            None,
        )
    })?;
    if reference_storage != reference.bases().len() {
        return Err(matrix_error(
            reference_length,
            query_length,
            MatrixAllocation::LogicalReferenceExtent,
            Some(reference_length),
            None,
        ));
    }
    if query_storage != query.bases().len() {
        return Err(matrix_error(
            reference_length,
            query_length,
            MatrixAllocation::LogicalQueryExtent,
            Some(query_length),
            None,
        ));
    }

    Ok(Dimensions {
        reference_length,
        query_length,
        rows,
        columns,
        cells,
        reference_storage,
        query_storage,
    })
}

fn validate_logical_dimensions(
    reference_length: u64,
    query_length: u64,
    limit: DpCellLimit,
) -> Result<(u64, u64, u64), DistanceError> {
    let rows = reference_length.checked_add(1).ok_or_else(|| {
        matrix_error(
            reference_length,
            query_length,
            MatrixAllocation::LogicalReferenceExtent,
            None,
            None,
        )
    })?;
    let columns = query_length.checked_add(1).ok_or_else(|| {
        matrix_error(
            reference_length,
            query_length,
            MatrixAllocation::LogicalQueryExtent,
            None,
            None,
        )
    })?;
    let cells = rows.checked_mul(columns).ok_or_else(|| {
        matrix_error(
            reference_length,
            query_length,
            MatrixAllocation::LogicalCells,
            None,
            None,
        )
    })?;
    if cells > limit.get() {
        return Err(DistanceError::ComputationLimitExceeded {
            requested_cells: cells,
            limit: limit.get(),
        });
    }
    Ok((rows, columns, cells))
}

fn preflight_distance(dimensions: &Dimensions) -> Result<usize, DistanceError> {
    let elements = checked_product(
        dimensions.columns,
        2,
        MatrixAllocation::DistanceRows,
        dimensions,
    )?;
    checked_allocation::<u64>(elements, MatrixAllocation::DistanceRows, dimensions)
}

fn preflight_alignment(dimensions: &Dimensions) -> Result<AlignmentCapacities, DistanceError> {
    let row_elements_u64 = checked_product(
        dimensions.columns,
        2,
        MatrixAllocation::DistanceRows,
        dimensions,
    )?;
    let row_elements =
        checked_allocation::<u64>(row_elements_u64, MatrixAllocation::DistanceRows, dimensions)?;
    checked_allocation::<u8>(
        row_elements_u64,
        MatrixAllocation::PathCountRows,
        dimensions,
    )?;

    let suffix_elements_u64 = checked_product(
        dimensions.cells,
        3,
        MatrixAllocation::SuffixScores,
        dimensions,
    )?;
    let suffix_elements = checked_allocation::<PathScore>(
        suffix_elements_u64,
        MatrixAllocation::SuffixScores,
        dimensions,
    )?;

    let trace_u64 = dimensions
        .reference_length
        .checked_add(dimensions.query_length)
        .ok_or_else(|| {
            matrix_error(
                dimensions.reference_length,
                dimensions.query_length,
                MatrixAllocation::TraceOperations,
                None,
                None,
            )
        })?;
    let trace_operations = checked_allocation::<CoreCigarOp>(
        trace_u64,
        MatrixAllocation::TraceOperations,
        dimensions,
    )?;
    checked_allocation::<CoreCigarRun>(trace_u64, MatrixAllocation::CigarRuns, dimensions)?;

    Ok(AlignmentCapacities {
        row_elements,
        suffix_elements,
        trace_operations,
    })
}

fn checked_product(
    left: u64,
    right: u64,
    allocation: MatrixAllocation,
    dimensions: &Dimensions,
) -> Result<u64, DistanceError> {
    left.checked_mul(right).ok_or_else(|| {
        matrix_error(
            dimensions.reference_length,
            dimensions.query_length,
            allocation,
            None,
            None,
        )
    })
}

fn checked_allocation<T>(
    elements: u64,
    allocation: MatrixAllocation,
    dimensions: &Dimensions,
) -> Result<usize, DistanceError> {
    let element_size = u64::try_from(size_of::<T>()).map_err(|_| {
        matrix_error(
            dimensions.reference_length,
            dimensions.query_length,
            allocation,
            Some(elements),
            None,
        )
    })?;
    let bytes = elements.checked_mul(element_size).ok_or_else(|| {
        matrix_error(
            dimensions.reference_length,
            dimensions.query_length,
            allocation,
            Some(elements),
            None,
        )
    })?;
    let capacity = usize::try_from(elements).map_err(|_| {
        matrix_error(
            dimensions.reference_length,
            dimensions.query_length,
            allocation,
            Some(elements),
            Some(bytes),
        )
    })?;
    let byte_capacity = usize::try_from(bytes).map_err(|_| {
        matrix_error(
            dimensions.reference_length,
            dimensions.query_length,
            allocation,
            Some(elements),
            Some(bytes),
        )
    })?;
    if byte_capacity > (usize::MAX >> 1) {
        return Err(matrix_error(
            dimensions.reference_length,
            dimensions.query_length,
            allocation,
            Some(elements),
            Some(bytes),
        ));
    }
    Ok(capacity)
}

const fn matrix_error(
    reference_length: u64,
    query_length: u64,
    allocation: MatrixAllocation,
    requested_elements: Option<u64>,
    requested_bytes: Option<u64>,
) -> DistanceError {
    DistanceError::MatrixSizeOverflow {
        reference_length,
        query_length,
        allocation,
        requested_elements,
        requested_bytes,
    }
}

fn distance_rows(
    reference: &NormalizedSequence,
    query: &NormalizedSequence,
    strand: CytosineStrand,
    dimensions: &Dimensions,
) -> Result<EditDistance, DistanceError> {
    let columns = dimensions.query_storage + 1;
    let mut previous = Vec::with_capacity(columns);
    previous.extend(0..=dimensions.query_length);
    let mut current = vec![0_u64; columns];

    for reference_index in 1..=dimensions.reference_storage {
        current[0] = storage_index_to_u64(reference_index, dimensions)?;
        for query_index in 1..=dimensions.query_storage {
            let deletion = checked_distance_add(previous[query_index], 1)?;
            let insertion = checked_distance_add(current[query_index - 1], 1)?;
            let relation_cost = classify_bases(
                reference.bases()[reference_index - 1],
                query.bases()[query_index - 1],
                strand,
            )
            .cost();
            let diagonal = checked_distance_add(previous[query_index - 1], relation_cost)?;
            current[query_index] = deletion.min(insertion).min(diagonal);
        }
        core::mem::swap(&mut previous, &mut current);
    }
    Ok(EditDistance::new(previous[dimensions.query_storage]))
}

fn primary_distance_and_count(
    reference: &NormalizedSequence,
    query: &NormalizedSequence,
    strand: CytosineStrand,
    dimensions: &Dimensions,
    row_elements: usize,
) -> Result<(EditDistance, u8), DistanceError> {
    let columns = dimensions.query_storage + 1;
    let mut previous_distance = Vec::with_capacity(columns);
    previous_distance.extend(0..=dimensions.query_length);
    let mut current_distance = vec![0_u64; columns];

    let mut previous_count = Vec::with_capacity(columns);
    previous_count.resize(columns, 1_u8);
    let mut current_count = Vec::with_capacity(columns);
    current_count.resize(columns, 1_u8);
    debug_assert_eq!(row_elements, columns * 2);

    for reference_index in 1..=dimensions.reference_storage {
        current_distance[0] = storage_index_to_u64(reference_index, dimensions)?;
        current_count[0] = 1;
        for query_index in 1..=dimensions.query_storage {
            let deletion = checked_distance_add(previous_distance[query_index], 1)?;
            let insertion = checked_distance_add(current_distance[query_index - 1], 1)?;
            let relation_cost = classify_bases(
                reference.bases()[reference_index - 1],
                query.bases()[query_index - 1],
                strand,
            )
            .cost();
            let diagonal = checked_distance_add(previous_distance[query_index - 1], relation_cost)?;
            let best = deletion.min(insertion).min(diagonal);
            current_distance[query_index] = best;

            let mut paths = 0_u8;
            if deletion == best {
                paths = capped_path_sum(paths, previous_count[query_index]);
            }
            if insertion == best {
                paths = capped_path_sum(paths, current_count[query_index - 1]);
            }
            if diagonal == best {
                paths = capped_path_sum(paths, previous_count[query_index - 1]);
            }
            current_count[query_index] = paths;
        }
        core::mem::swap(&mut previous_distance, &mut current_distance);
        core::mem::swap(&mut previous_count, &mut current_count);
    }

    Ok((
        EditDistance::new(previous_distance[dimensions.query_storage]),
        previous_count[dimensions.query_storage],
    ))
}

fn primary_distance_and_count_banded(
    reference: &NormalizedSequence,
    query: &NormalizedSequence,
    strand: CytosineStrand,
    dimensions: &Dimensions,
    row_elements: usize,
    band: usize,
) -> Result<(EditDistance, u8), DistanceError> {
    const UNREACHABLE: u64 = u64::MAX;

    let columns = dimensions.query_storage + 1;
    let mut previous_distance = vec![UNREACHABLE; columns];
    let mut current_distance = vec![UNREACHABLE; columns];
    let mut previous_count = vec![0_u8; columns];
    let mut current_count = vec![0_u8; columns];
    debug_assert_eq!(row_elements, columns * 2);

    let (_, initial_end) = band_bounds(0, dimensions.query_storage, band);
    for query_index in 0..=initial_end {
        previous_distance[query_index] = storage_index_to_u64(query_index, dimensions)?;
        previous_count[query_index] = 1;
    }

    for reference_index in 1..=dimensions.reference_storage {
        current_distance.fill(UNREACHABLE);
        current_count.fill(0);
        let (query_start, query_end) = band_bounds(reference_index, dimensions.query_storage, band);
        for query_index in query_start..=query_end {
            let deletion = banded_distance_add(previous_distance[query_index], 1)?;
            let insertion = if query_index == 0 {
                UNREACHABLE
            } else {
                banded_distance_add(current_distance[query_index - 1], 1)?
            };
            let diagonal = if query_index == 0 {
                UNREACHABLE
            } else {
                let relation_cost = classify_bases(
                    reference.bases()[reference_index - 1],
                    query.bases()[query_index - 1],
                    strand,
                )
                .cost();
                banded_distance_add(previous_distance[query_index - 1], relation_cost)?
            };
            let best = deletion.min(insertion).min(diagonal);
            if best == UNREACHABLE {
                return Err(invariant_error(
                    AlignmentInvariant::MissingTraceStep,
                    storage_index_to_u64(reference_index, dimensions)?,
                    storage_index_to_u64(query_index, dimensions)?,
                    None,
                    None,
                ));
            }
            current_distance[query_index] = best;

            let mut paths = 0_u8;
            if deletion == best {
                paths = capped_path_sum(paths, previous_count[query_index]);
            }
            if insertion == best {
                paths = capped_path_sum(paths, current_count[query_index - 1]);
            }
            if diagonal == best {
                paths = capped_path_sum(paths, previous_count[query_index - 1]);
            }
            current_count[query_index] = paths;
        }
        core::mem::swap(&mut previous_distance, &mut current_distance);
        core::mem::swap(&mut previous_count, &mut current_count);
    }

    Ok((
        EditDistance::new(previous_distance[dimensions.query_storage]),
        previous_count[dimensions.query_storage],
    ))
}

const fn banded_distance_add(accumulated: u64, increment: u64) -> Result<u64, DistanceError> {
    if accumulated == u64::MAX {
        Ok(u64::MAX)
    } else {
        checked_distance_add(accumulated, increment)
    }
}

fn band_bounds(reference_index: usize, query_length: usize, band: usize) -> (usize, usize) {
    (
        reference_index.saturating_sub(band),
        reference_index.saturating_add(band).min(query_length),
    )
}

const fn is_inside_band(reference_index: usize, query_index: usize, band: usize) -> bool {
    reference_index.abs_diff(query_index) <= band
}

const fn capped_path_sum(left: u8, right: u8) -> u8 {
    let sum = left.saturating_add(right);
    if sum > 2 { 2 } else { sum }
}

fn storage_index_to_u64(index: usize, dimensions: &Dimensions) -> Result<u64, DistanceError> {
    u64::try_from(index).map_err(|_| {
        matrix_error(
            dimensions.reference_length,
            dimensions.query_length,
            MatrixAllocation::LogicalReferenceExtent,
            None,
            None,
        )
    })
}

const fn checked_distance_add(accumulated: u64, increment: u64) -> Result<u64, DistanceError> {
    match accumulated.checked_add(increment) {
        Some(value) => Ok(value),
        None => Err(DistanceError::DistanceOverflow {
            accumulated,
            increment,
        }),
    }
}

const PREVIOUS_D: usize = 0;
const PREVIOUS_I: usize = 1;
const PREVIOUS_OTHER: usize = 2;
const PREVIOUS_STATE_COUNT: usize = 3;

#[derive(Clone, Copy, Debug, Default, Eq, Ord, PartialEq, PartialOrd)]
struct PathScore {
    distance: u64,
    gap_bases: u64,
    gap_runs: u64,
}

impl PathScore {
    const ZERO: Self = Self {
        distance: 0,
        gap_bases: 0,
        gap_runs: 0,
    };

    fn prepend(
        self,
        operation: CoreCigarOp,
        previous_state: usize,
        relation_cost: u64,
    ) -> Result<Self, DistanceError> {
        let distance_increment = if operation == CoreCigarOp::M {
            relation_cost
        } else {
            1
        };
        let gap_increment = u64::from(operation.is_gap());
        let run_increment =
            u64::from(operation.is_gap() && previous_state != operation_state(operation));
        Ok(Self {
            distance: checked_distance_add(self.distance, distance_increment)?,
            gap_bases: checked_trace_add(TraceScoreField::GapBases, self.gap_bases, gap_increment)?,
            gap_runs: checked_trace_add(TraceScoreField::GapRuns, self.gap_runs, run_increment)?,
        })
    }
}

const fn checked_trace_add(
    field: TraceScoreField,
    accumulated: u64,
    increment: u64,
) -> Result<u64, DistanceError> {
    match accumulated.checked_add(increment) {
        Some(value) => Ok(value),
        None => Err(DistanceError::TraceScoreOverflow {
            field,
            accumulated,
            increment,
        }),
    }
}

const fn operation_state(operation: CoreCigarOp) -> usize {
    match operation {
        CoreCigarOp::D => PREVIOUS_D,
        CoreCigarOp::I => PREVIOUS_I,
        CoreCigarOp::M => PREVIOUS_OTHER,
    }
}

fn build_suffix_scores(
    reference: &NormalizedSequence,
    query: &NormalizedSequence,
    strand: CytosineStrand,
    dimensions: &Dimensions,
    suffix_elements: usize,
) -> Result<Vec<PathScore>, DistanceError> {
    let mut scores = Vec::with_capacity(suffix_elements);
    scores.resize(suffix_elements, PathScore::ZERO);

    for reference_index in (0..=dimensions.reference_storage).rev() {
        for query_index in (0..=dimensions.query_storage).rev() {
            if reference_index == dimensions.reference_storage
                && query_index == dimensions.query_storage
            {
                continue;
            }
            for previous_state in 0..PREVIOUS_STATE_COUNT {
                let mut best = None;
                for operation in CoreCigarOp::LEXICOGRAPHIC {
                    if let Some(candidate) = suffix_candidate(
                        reference,
                        query,
                        strand,
                        dimensions,
                        &scores,
                        reference_index,
                        query_index,
                        previous_state,
                        operation,
                    )? {
                        best = Some(
                            best.map_or(candidate, |current: PathScore| current.min(candidate)),
                        );
                    }
                }
                let Some(best_score) = best else {
                    return Err(invariant_error(
                        AlignmentInvariant::MissingSuffixPath,
                        storage_index_to_u64(reference_index, dimensions)?,
                        storage_index_to_u64(query_index, dimensions)?,
                        None,
                        None,
                    ));
                };
                let index = suffix_index(reference_index, query_index, previous_state, dimensions);
                scores[index] = best_score;
            }
        }
    }
    Ok(scores)
}

#[allow(clippy::too_many_arguments)]
fn suffix_candidate(
    reference: &NormalizedSequence,
    query: &NormalizedSequence,
    strand: CytosineStrand,
    dimensions: &Dimensions,
    scores: &[PathScore],
    reference_index: usize,
    query_index: usize,
    previous_state: usize,
    operation: CoreCigarOp,
) -> Result<Option<PathScore>, DistanceError> {
    let step = match operation {
        CoreCigarOp::D if reference_index < dimensions.reference_storage => {
            Some((reference_index + 1, query_index, 1))
        }
        CoreCigarOp::I if query_index < dimensions.query_storage => {
            Some((reference_index, query_index + 1, 1))
        }
        CoreCigarOp::M
            if reference_index < dimensions.reference_storage
                && query_index < dimensions.query_storage =>
        {
            Some((
                reference_index + 1,
                query_index + 1,
                classify_bases(
                    reference.bases()[reference_index],
                    query.bases()[query_index],
                    strand,
                )
                .cost(),
            ))
        }
        CoreCigarOp::D | CoreCigarOp::I | CoreCigarOp::M => None,
    };
    let Some((next_reference, next_query, relation_cost)) = step else {
        return Ok(None);
    };
    let successor = suffix_score(
        scores,
        next_reference,
        next_query,
        operation_state(operation),
        dimensions,
    );
    successor
        .prepend(operation, previous_state, relation_cost)
        .map(Some)
}

fn build_suffix_scores_banded(
    reference: &NormalizedSequence,
    query: &NormalizedSequence,
    strand: CytosineStrand,
    dimensions: &Dimensions,
    _suffix_elements: usize,
    band: usize,
) -> Result<BandedSuffixScores, DistanceError> {
    let mut scores = BandedSuffixScores::new(dimensions, band);

    for reference_index in (0..=dimensions.reference_storage).rev() {
        let (query_start, query_end) = band_bounds(reference_index, dimensions.query_storage, band);
        for query_index in (query_start..=query_end).rev() {
            if reference_index == dimensions.reference_storage
                && query_index == dimensions.query_storage
            {
                continue;
            }
            for previous_state in 0..PREVIOUS_STATE_COUNT {
                let mut best = None;
                for operation in CoreCigarOp::LEXICOGRAPHIC {
                    if let Some(candidate) = suffix_candidate_banded(
                        reference,
                        query,
                        strand,
                        dimensions,
                        &scores,
                        reference_index,
                        query_index,
                        previous_state,
                        operation,
                        band,
                    )? {
                        best = Some(
                            best.map_or(candidate, |current: PathScore| current.min(candidate)),
                        );
                    }
                }
                let Some(best_score) = best else {
                    return Err(invariant_error(
                        AlignmentInvariant::MissingSuffixPath,
                        storage_index_to_u64(reference_index, dimensions)?,
                        storage_index_to_u64(query_index, dimensions)?,
                        None,
                        None,
                    ));
                };
                scores.set(reference_index, query_index, previous_state, best_score);
            }
        }
    }
    Ok(scores)
}

#[allow(clippy::too_many_arguments)]
fn suffix_candidate_banded(
    reference: &NormalizedSequence,
    query: &NormalizedSequence,
    strand: CytosineStrand,
    dimensions: &Dimensions,
    scores: &BandedSuffixScores,
    reference_index: usize,
    query_index: usize,
    previous_state: usize,
    operation: CoreCigarOp,
    band: usize,
) -> Result<Option<PathScore>, DistanceError> {
    let step = match operation {
        CoreCigarOp::D if reference_index < dimensions.reference_storage => {
            Some((reference_index + 1, query_index, 1))
        }
        CoreCigarOp::I if query_index < dimensions.query_storage => {
            Some((reference_index, query_index + 1, 1))
        }
        CoreCigarOp::M
            if reference_index < dimensions.reference_storage
                && query_index < dimensions.query_storage =>
        {
            Some((
                reference_index + 1,
                query_index + 1,
                classify_bases(
                    reference.bases()[reference_index],
                    query.bases()[query_index],
                    strand,
                )
                .cost(),
            ))
        }
        CoreCigarOp::D | CoreCigarOp::I | CoreCigarOp::M => None,
    };
    let Some((next_reference, next_query, relation_cost)) = step else {
        return Ok(None);
    };
    if !is_inside_band(next_reference, next_query, band) {
        return Ok(None);
    }
    let successor = scores.get(next_reference, next_query, operation_state(operation));
    successor
        .prepend(operation, previous_state, relation_cost)
        .map(Some)
}

fn canonical_trace(
    reference: &NormalizedSequence,
    query: &NormalizedSequence,
    strand: CytosineStrand,
    dimensions: &Dimensions,
    scores: &[PathScore],
    trace_capacity: usize,
) -> Result<Vec<CoreCigarOp>, DistanceError> {
    let mut operations = Vec::with_capacity(trace_capacity);
    let mut reference_index = 0_usize;
    let mut query_index = 0_usize;
    let mut previous_state = PREVIOUS_OTHER;

    while reference_index < dimensions.reference_storage || query_index < dimensions.query_storage {
        let target = suffix_score(
            scores,
            reference_index,
            query_index,
            previous_state,
            dimensions,
        );
        let mut chosen = None;
        for operation in CoreCigarOp::LEXICOGRAPHIC {
            if let Some(candidate) = suffix_candidate(
                reference,
                query,
                strand,
                dimensions,
                scores,
                reference_index,
                query_index,
                previous_state,
                operation,
            )? && candidate == target
            {
                chosen = Some(operation);
                break;
            }
        }
        let Some(operation) = chosen else {
            return Err(invariant_error(
                AlignmentInvariant::MissingTraceStep,
                storage_index_to_u64(reference_index, dimensions)?,
                storage_index_to_u64(query_index, dimensions)?,
                None,
                None,
            ));
        };
        operations.push(operation);
        match operation {
            CoreCigarOp::D => reference_index += 1,
            CoreCigarOp::I => query_index += 1,
            CoreCigarOp::M => {
                reference_index += 1;
                query_index += 1;
            }
        }
        previous_state = operation_state(operation);
    }
    Ok(operations)
}

#[allow(clippy::too_many_arguments)]
fn canonical_trace_banded(
    reference: &NormalizedSequence,
    query: &NormalizedSequence,
    strand: CytosineStrand,
    dimensions: &Dimensions,
    scores: &BandedSuffixScores,
    trace_capacity: usize,
    band: usize,
) -> Result<Vec<CoreCigarOp>, DistanceError> {
    let mut operations = Vec::with_capacity(trace_capacity);
    let mut reference_index = 0_usize;
    let mut query_index = 0_usize;
    let mut previous_state = PREVIOUS_OTHER;

    while reference_index < dimensions.reference_storage || query_index < dimensions.query_storage {
        let target = scores.get(reference_index, query_index, previous_state);
        let mut chosen = None;
        for operation in CoreCigarOp::LEXICOGRAPHIC {
            if let Some(candidate) = suffix_candidate_banded(
                reference,
                query,
                strand,
                dimensions,
                scores,
                reference_index,
                query_index,
                previous_state,
                operation,
                band,
            )? && candidate == target
            {
                chosen = Some(operation);
                break;
            }
        }
        let Some(operation) = chosen else {
            return Err(invariant_error(
                AlignmentInvariant::MissingTraceStep,
                storage_index_to_u64(reference_index, dimensions)?,
                storage_index_to_u64(query_index, dimensions)?,
                None,
                None,
            ));
        };
        operations.push(operation);
        match operation {
            CoreCigarOp::D => reference_index += 1,
            CoreCigarOp::I => query_index += 1,
            CoreCigarOp::M => {
                reference_index += 1;
                query_index += 1;
            }
        }
        previous_state = operation_state(operation);
    }
    Ok(operations)
}

fn suffix_score(
    scores: &[PathScore],
    reference_index: usize,
    query_index: usize,
    previous_state: usize,
    dimensions: &Dimensions,
) -> PathScore {
    scores[suffix_index(reference_index, query_index, previous_state, dimensions)]
}

struct BandedSuffixScores {
    scores: Vec<PathScore>,
    band: usize,
    row_width: usize,
}

impl BandedSuffixScores {
    fn new(dimensions: &Dimensions, requested_band: usize) -> Self {
        let band = requested_band.min(dimensions.reference_storage.max(dimensions.query_storage));
        let row_width = band.saturating_mul(2).saturating_add(1);
        let rows = dimensions.reference_storage.saturating_add(1);
        let elements = rows
            .checked_mul(row_width)
            .and_then(|value| value.checked_mul(PREVIOUS_STATE_COUNT))
            .expect("validated alignment dimensions fit compact band storage");
        Self {
            scores: vec![PathScore::ZERO; elements],
            band,
            row_width,
        }
    }

    #[inline]
    fn index(&self, reference_index: usize, query_index: usize, previous_state: usize) -> usize {
        debug_assert!(previous_state < PREVIOUS_STATE_COUNT);
        debug_assert!(is_inside_band(reference_index, query_index, self.band));
        let diagonal = query_index + self.band - reference_index;
        debug_assert!(diagonal < self.row_width);
        (reference_index * self.row_width + diagonal) * PREVIOUS_STATE_COUNT + previous_state
    }

    #[inline]
    fn get(&self, reference_index: usize, query_index: usize, previous_state: usize) -> PathScore {
        self.scores[self.index(reference_index, query_index, previous_state)]
    }

    #[inline]
    fn set(
        &mut self,
        reference_index: usize,
        query_index: usize,
        previous_state: usize,
        score: PathScore,
    ) {
        let index = self.index(reference_index, query_index, previous_state);
        self.scores[index] = score;
    }
}

fn suffix_index(
    reference_index: usize,
    query_index: usize,
    previous_state: usize,
    dimensions: &Dimensions,
) -> usize {
    debug_assert!(previous_state < PREVIOUS_STATE_COUNT);
    let columns = dimensions.query_storage + 1;
    ((reference_index * columns) + query_index) * PREVIOUS_STATE_COUNT + previous_state
}

const fn invariant_error(
    invariant: AlignmentInvariant,
    reference_index: u64,
    query_index: u64,
    expected: Option<u64>,
    observed: Option<u64>,
) -> DistanceError {
    DistanceError::AlignmentInvariant {
        invariant,
        reference_index,
        query_index,
        expected,
        observed,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bsbit_core::alphabet::Base;

    fn synthetic_dimensions(
        reference_length: u64,
        query_length: u64,
        columns: u64,
        cells: u64,
    ) -> Dimensions {
        Dimensions {
            reference_length,
            query_length,
            rows: 1,
            columns,
            cells,
            reference_storage: 0,
            query_storage: 0,
        }
    }

    #[test]
    fn logical_dimension_failures_are_ordered_and_structured() {
        assert_eq!(
            validate_logical_dimensions(u64::MAX, 0, DpCellLimit::MAX),
            Err(DistanceError::MatrixSizeOverflow {
                reference_length: u64::MAX,
                query_length: 0,
                allocation: MatrixAllocation::LogicalReferenceExtent,
                requested_elements: None,
                requested_bytes: None,
            })
        );
        assert_eq!(
            validate_logical_dimensions(0, u64::MAX, DpCellLimit::MAX),
            Err(DistanceError::MatrixSizeOverflow {
                reference_length: 0,
                query_length: u64::MAX,
                allocation: MatrixAllocation::LogicalQueryExtent,
                requested_elements: None,
                requested_bytes: None,
            })
        );
        assert_eq!(
            validate_logical_dimensions(u64::MAX - 1, 1, DpCellLimit::MAX),
            Err(DistanceError::MatrixSizeOverflow {
                reference_length: u64::MAX - 1,
                query_length: 1,
                allocation: MatrixAllocation::LogicalCells,
                requested_elements: None,
                requested_bytes: None,
            })
        );
        assert_eq!(
            validate_logical_dimensions(2, 3, DpCellLimit::new(11)),
            Err(DistanceError::ComputationLimitExceeded {
                requested_cells: 12,
                limit: 11,
            })
        );
        assert_eq!(
            validate_logical_dimensions(2, 3, DpCellLimit::new(12)),
            Ok((3, 4, 12))
        );
    }

    #[test]
    fn alignment_preflight_reports_the_first_planned_allocation_failure() {
        let distance_rows = synthetic_dimensions(0, 0, u64::MAX, 1);
        assert!(matches!(
            preflight_alignment(&distance_rows),
            Err(DistanceError::MatrixSizeOverflow {
                allocation: MatrixAllocation::DistanceRows,
                ..
            })
        ));

        let suffix_scores = synthetic_dimensions(0, 0, 1, u64::MAX);
        assert!(matches!(
            preflight_alignment(&suffix_scores),
            Err(DistanceError::MatrixSizeOverflow {
                allocation: MatrixAllocation::SuffixScores,
                ..
            })
        ));

        let trace_operations = synthetic_dimensions(u64::MAX, 1, 1, 1);
        assert!(matches!(
            preflight_alignment(&trace_operations),
            Err(DistanceError::MatrixSizeOverflow {
                allocation: MatrixAllocation::TraceOperations,
                ..
            })
        ));

        let maximum_planned_bytes =
            u64::try_from(usize::MAX >> 1).expect("supported pointer widths fit u64");
        let run_size =
            u64::try_from(size_of::<CoreCigarRun>()).expect("supported pointer widths fit u64");
        let operation_size =
            u64::try_from(size_of::<CoreCigarOp>()).expect("supported pointer widths fit u64");
        assert!(run_size > operation_size);
        let trace_elements = maximum_planned_bytes / run_size + 1;
        let operation_bytes = trace_elements
            .checked_mul(operation_size)
            .expect("operation capacity remains representable");
        assert!(operation_bytes <= maximum_planned_bytes);
        let run_bytes = trace_elements
            .checked_mul(run_size)
            .expect("run capacity remains representable in u64");

        let cigar_runs = synthetic_dimensions(trace_elements, 0, 1, 1);
        assert_eq!(
            preflight_alignment(&cigar_runs),
            Err(DistanceError::MatrixSizeOverflow {
                reference_length: trace_elements,
                query_length: 0,
                allocation: MatrixAllocation::CigarRuns,
                requested_elements: Some(trace_elements),
                requested_bytes: Some(run_bytes),
            })
        );
    }

    #[test]
    fn exact_banded_traceback_matches_full_oracle_exhaustively_through_length_three() {
        let sequences = sequences_through(3);
        for strand in [CytosineStrand::Top, CytosineStrand::Bottom] {
            for reference in &sequences {
                for query in &sequences {
                    let oracle = global_bs_alignment(reference, query, strand, DpCellLimit::MAX)
                        .expect("small oracle alignment must fit");
                    let observed = global_bs_alignment_banded_exact(
                        reference,
                        query,
                        strand,
                        oracle.distance(),
                        DpCellLimit::MAX,
                    )
                    .expect("the oracle distance certifies the exact band");
                    assert_eq!(
                        observed,
                        oracle,
                        "reference={:?} query={:?} strand={strand:?}",
                        reference.bases(),
                        query.bases(),
                    );
                }
            }
        }
    }

    #[test]
    fn exact_banded_traceback_rejects_a_wrong_filter_distance() {
        let reference = NormalizedSequence::from_bases([Base::A]);
        let query = NormalizedSequence::from_bases([Base::C]);
        assert_eq!(
            global_bs_alignment_banded_exact(
                &reference,
                &query,
                CytosineStrand::Top,
                EditDistance::new(0),
                DpCellLimit::MAX,
            ),
            Err(DistanceError::AlignmentInvariant {
                invariant: AlignmentInvariant::ExpectedDistanceMismatch,
                reference_index: 1,
                query_index: 1,
                expected: Some(0),
                observed: Some(1),
            })
        );
    }

    fn sequences_through(maximum_length: usize) -> Vec<NormalizedSequence> {
        let mut sequences = Vec::new();
        for length in 0..=maximum_length {
            let count = Base::ALL
                .len()
                .pow(u32::try_from(length).expect("small length"));
            for ordinal in 0..count {
                let mut remaining = ordinal;
                let mut bases = Vec::with_capacity(length);
                for _ in 0..length {
                    bases.push(Base::ALL[remaining % Base::ALL.len()]);
                    remaining /= Base::ALL.len();
                }
                sequences.push(NormalizedSequence::from_bases(bases));
            }
        }
        sequences
    }
}
