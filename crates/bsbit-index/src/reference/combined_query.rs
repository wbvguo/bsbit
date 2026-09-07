//! Combined-image backend boundary and owner-bound query facade.
//!
//! The public items are re-exported by `reference`, preserving the established
//! API path. The storage implementation remains crate-private while the hot
//! query handle retains its concrete type so rank and locate calls can inline.

use core::fmt;

#[cfg(feature = "combined-index")]
use super::{
    CombinedReferenceCatalog, ReferenceBuildError, ReferenceIndex, ReferenceLocateError,
    ReferenceLocateInvariant, ReferenceLocateMetrics,
};
#[cfg(feature = "combined-index")]
use crate::storage::combined::CombinedIndex;
#[cfg(feature = "combined-index")]
use crate::storage::fm::{FmInterval, ProjectedBase, SearchBase};
#[cfg(feature = "combined-index")]
use bsbit_core::bisulfite::BisulfiteStrand;

#[cfg(all(feature = "combined-index", not(test)))]
pub(crate) type PrivateCombinedBackend = CombinedIndex;
#[cfg(all(feature = "combined-index", test))]
pub(crate) type PrivateCombinedBackend = dyn PrivateCombinedIndex;

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

/// Crate-private operation contract for the frozen combined-direction index.
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
    pub(super) index: Box<PrivateCombinedBackend>,
}

#[cfg(feature = "combined-index")]
impl PrivateCombinedReference {
    /// Wraps the validated crate-private combined-index implementation.
    #[must_use]
    pub(crate) fn new(index: CombinedIndex) -> Self {
        Self {
            index: Box::new(index),
        }
    }

    #[cfg(test)]
    pub(crate) fn from_test(index: Box<dyn PrivateCombinedIndex>) -> Self {
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

/// Number of projected symbols represented by the combined image's dense lookup.
#[cfg(feature = "combined-index")]
pub const COMBINED_EXACT_LOOKUP_BASES: usize = 16;

/// A checked interval owned by one combined index instance.
#[cfg(feature = "combined-index")]
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub struct CombinedIndexInterval {
    owner_token: usize,
    interval: FmInterval,
}

#[cfg(feature = "combined-index")]
impl CombinedIndexInterval {
    /// Returns the number of suffix-array rows in this interval.
    #[must_use]
    pub const fn len(self) -> u64 {
        self.interval.len()
    }

    /// Returns whether this interval has no rows.
    #[must_use]
    pub const fn is_empty(self) -> bool {
        self.interval.is_empty()
    }
}

/// One coordinate recovered from a combined-index interval.
#[cfg(feature = "combined-index")]
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct CombinedIndexCoordinate {
    contig_ordinal: u64,
    strand: BisulfiteStrand,
    start: u64,
}

#[cfg(feature = "combined-index")]
impl CombinedIndexCoordinate {
    /// Returns the zero-based contig ordinal.
    #[must_use]
    pub const fn contig_ordinal(self) -> u64 {
        self.contig_ordinal
    }

    /// Returns the bisulfite reference orientation.
    #[must_use]
    pub const fn strand(self) -> BisulfiteStrand {
        self.strand
    }

    /// Returns the zero-based contig-local start.
    #[must_use]
    pub const fn start(self) -> u64 {
        self.start
    }
}

/// Failure of a neutral combined-index interval or coordinate operation.
#[cfg(feature = "combined-index")]
#[derive(Clone, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum CombinedIndexQueryError {
    /// An interval belongs to another reference instance.
    ForeignInterval,
    /// The packed backend rejected an interval operation.
    Backend(CombinedIndexBackendError),
    /// A checked combined coordinate could not be recovered.
    Coordinate(ReferenceLocateError),
}

#[cfg(feature = "combined-index")]
impl fmt::Display for CombinedIndexQueryError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ForeignInterval => {
                formatter.write_str("combined interval belongs to another reference index")
            }
            Self::Backend(error) => write!(formatter, "combined interval query failed: {error}"),
            Self::Coordinate(error) => {
                write!(formatter, "combined coordinate recovery failed: {error}")
            }
        }
    }
}

#[cfg(feature = "combined-index")]
impl std::error::Error for CombinedIndexQueryError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Backend(error) => Some(error),
            Self::Coordinate(error) => Some(error),
            Self::ForeignInterval => None,
        }
    }
}

/// Borrowed owner-bound interval, rank, and locate primitives for the combined image.
///
/// The handle does not choose seed lengths, stop at singleton intervals, or
/// interpret coordinates as alignment candidates.
#[cfg(feature = "combined-index")]
#[derive(Clone, Copy)]
pub struct CombinedIndexQuery<'index> {
    pub(super) reference: &'index ReferenceIndex,
    pub(super) catalog: &'index CombinedReferenceCatalog,
    pub(super) backend: &'index PrivateCombinedBackend,
    pub(super) owner_token: usize,
}

#[cfg(feature = "combined-index")]
impl fmt::Debug for CombinedIndexQuery<'_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("CombinedIndexQuery")
            .finish_non_exhaustive()
    }
}

#[cfg(feature = "combined-index")]
impl CombinedIndexQuery<'_> {
    /// Returns the complete exact interval for a caller-supplied projected pattern.
    ///
    /// # Errors
    ///
    /// Returns a backend error if the validated image rejects the query.
    pub fn exact_interval(
        self,
        reversed_projected_pattern: &[SearchBase],
    ) -> Result<Option<CombinedIndexInterval>, CombinedIndexQueryError> {
        PrivateCombinedIndex::exact_interval(self.backend, reversed_projected_pattern)
            .map(|interval| interval.map(|interval| self.bind(interval)))
            .map_err(CombinedIndexQueryError::Backend)
    }

    /// Projected-digit counterpart for complete exact lookup.
    ///
    /// # Errors
    ///
    /// Returns a backend error if the validated image rejects the query.
    pub fn exact_projected_interval(
        self,
        reversed_projected_pattern: &[ProjectedBase],
    ) -> Result<Option<CombinedIndexInterval>, CombinedIndexQueryError> {
        PrivateCombinedIndex::exact_projected_interval(self.backend, reversed_projected_pattern)
            .map(|interval| interval.map(|interval| self.bind(interval)))
            .map_err(CombinedIndexQueryError::Backend)
    }

    /// Looks up one caller-selected exact suffix in the dense table.
    ///
    /// # Errors
    ///
    /// Returns a backend error if the validated image rejects the lookup.
    pub fn lookup_interval(
        self,
        projected_suffix: &[SearchBase],
    ) -> Result<Option<CombinedIndexInterval>, CombinedIndexQueryError> {
        PrivateCombinedIndex::lookup_interval(self.backend, projected_suffix)
            .map(|interval| interval.map(|interval| self.bind(interval)))
            .map_err(CombinedIndexQueryError::Backend)
    }

    /// Projected-digit counterpart for dense exact lookup.
    ///
    /// # Errors
    ///
    /// Returns a backend error if the validated image rejects the lookup.
    pub fn lookup_projected_interval(
        self,
        projected_suffix: &[ProjectedBase],
    ) -> Result<Option<CombinedIndexInterval>, CombinedIndexQueryError> {
        PrivateCombinedIndex::lookup_projected_interval(self.backend, projected_suffix)
            .map(|interval| interval.map(|interval| self.bind(interval)))
            .map_err(CombinedIndexQueryError::Backend)
    }

    /// Prepends one projected symbol to a checked interval.
    ///
    /// # Errors
    ///
    /// Rejects a foreign interval or a backend extension failure.
    pub fn backward_extend(
        self,
        interval: CombinedIndexInterval,
        symbol: SearchBase,
    ) -> Result<CombinedIndexInterval, CombinedIndexQueryError> {
        let interval = self.validate(interval)?;
        PrivateCombinedIndex::backward_extend_interval(self.backend, interval, symbol)
            .map(|interval| self.bind(interval))
            .map_err(CombinedIndexQueryError::Backend)
    }

    /// Projected-digit counterpart for one-symbol backward extension.
    ///
    /// # Errors
    ///
    /// Rejects a foreign interval or a backend extension failure.
    pub fn backward_extend_projected(
        self,
        interval: CombinedIndexInterval,
        symbol: ProjectedBase,
    ) -> Result<CombinedIndexInterval, CombinedIndexQueryError> {
        let interval = self.validate(interval)?;
        PrivateCombinedIndex::backward_extend_projected_interval(self.backend, interval, symbol)
            .map(|interval| self.bind(interval))
            .map_err(CombinedIndexQueryError::Backend)
    }

    /// Extends independent intervals by one symbol in one physical rank round.
    ///
    /// # Errors
    ///
    /// Rejects foreign intervals, inconsistent slice lengths, or a backend failure.
    pub fn backward_extend_intervals(
        self,
        intervals: &[CombinedIndexInterval],
        symbols: &[SearchBase],
        output: &mut [CombinedIndexInterval],
    ) -> Result<(), CombinedIndexQueryError> {
        self.backward_extend_intervals_inner(intervals, output, |private, private_output| {
            PrivateCombinedIndex::backward_extend_intervals(
                self.backend,
                private,
                symbols,
                private_output,
            )
        })
    }

    /// Projected-digit counterpart for one physical rank round.
    ///
    /// # Errors
    ///
    /// Rejects foreign intervals, inconsistent slice lengths, or a backend failure.
    pub fn backward_extend_projected_intervals(
        self,
        intervals: &[CombinedIndexInterval],
        symbols: &[ProjectedBase],
        output: &mut [CombinedIndexInterval],
    ) -> Result<(), CombinedIndexQueryError> {
        self.backward_extend_intervals_inner(intervals, output, |private, private_output| {
            PrivateCombinedIndex::backward_extend_projected_intervals(
                self.backend,
                private,
                symbols,
                private_output,
            )
        })
    }

    /// Resolves independent projected suffixes in one physical backend call.
    ///
    /// This is a neutral index operation: the caller supplies both stopping
    /// rules, while the backend decides only how to schedule table reads and
    /// rank rounds.
    ///
    /// # Errors
    ///
    /// Rejects inconsistent output dimensions, more than 64 lanes, or a
    /// backend query failure.
    pub fn resolve_projected_suffix_intervals(
        self,
        patterns: &[&[ProjectedBase]],
        minimum_suffix_bases: usize,
        stop_interval_length: u64,
        output: &mut [Option<(CombinedIndexInterval, u64)>],
    ) -> Result<(), CombinedIndexQueryError> {
        const MAX_LANES: usize = 64;
        if patterns.len() != output.len() || patterns.len() > MAX_LANES || stop_interval_length == 0
        {
            return Err(CombinedIndexQueryError::Backend(
                CombinedIndexBackendError::Structure,
            ));
        }
        let mut private = [None; MAX_LANES];
        PrivateCombinedIndex::resolve_projected_suffix_intervals(
            self.backend,
            patterns,
            minimum_suffix_bases,
            stop_interval_length,
            &mut private[..patterns.len()],
        )
        .map_err(CombinedIndexQueryError::Backend)?;
        for (destination, result) in output.iter_mut().zip(private) {
            *destination =
                result.map(|(interval, matched_bases)| (self.bind(interval), matched_bases));
        }
        Ok(())
    }

    /// Streams one interval as combined-reference coordinates.
    ///
    /// # Errors
    ///
    /// Rejects a foreign interval, backend failure, or coordinate invariant failure.
    pub fn visit_interval(
        self,
        interval: CombinedIndexInterval,
        matched_bases: u64,
        pattern_offset: u64,
        pattern_len: u64,
        visitor: &mut dyn FnMut(CombinedIndexCoordinate) -> bool,
    ) -> Result<ReferenceLocateMetrics, CombinedIndexQueryError> {
        let interval = self.validate(interval)?;
        let mut stopped = false;
        let metrics =
            PrivateCombinedIndex::visit_interval(self.backend, interval, &mut |position| {
                let Some(coordinate) =
                    self.recover_coordinate(position, matched_bases, pattern_offset, pattern_len)
                else {
                    return true;
                };
                let keep_going = visitor(coordinate);
                stopped = !keep_going;
                keep_going
            })
            .map_err(CombinedIndexQueryError::Backend)?;
        if !stopped && metrics.located_rows() != interval.len() {
            return Err(CombinedIndexQueryError::Coordinate(
                ReferenceLocateError::Invariant {
                    invariant: ReferenceLocateInvariant::OffsetCount,
                    expected: interval.len(),
                    observed: metrics.located_rows(),
                },
            ));
        }
        Ok(ReferenceLocateMetrics {
            located_coordinates: metrics.located_rows(),
            lf_steps: metrics.lf_steps(),
            rank_operations: metrics.rank_operations(),
            interval_nodes: metrics.interval_nodes(),
        })
    }

    /// Streams two complete checked intervals through the backend's paired
    /// locate primitive and reports physical work per lane.
    ///
    /// # Errors
    ///
    /// Rejects foreign intervals, backend failures, or locate-count mismatches.
    pub fn visit_raw_intervals_two_lanes_complete(
        self,
        intervals: [CombinedIndexInterval; 2],
        visitor: &mut dyn FnMut(usize, u64),
    ) -> Result<[ReferenceLocateMetrics; 2], CombinedIndexQueryError> {
        let private = [self.validate(intervals[0])?, self.validate(intervals[1])?];
        let metrics = PrivateCombinedIndex::visit_intervals_two_lanes_complete(
            self.backend,
            private,
            visitor,
        )
        .map_err(CombinedIndexQueryError::Backend)?;
        for lane in 0..2 {
            if metrics[lane].located_rows() != private[lane].len() {
                return Err(CombinedIndexQueryError::Coordinate(
                    ReferenceLocateError::Invariant {
                        invariant: ReferenceLocateInvariant::OffsetCount,
                        expected: private[lane].len(),
                        observed: metrics[lane].located_rows(),
                    },
                ));
            }
        }
        Ok(metrics.map(private_locate_metrics))
    }

    fn backward_extend_intervals_inner(
        self,
        intervals: &[CombinedIndexInterval],
        output: &mut [CombinedIndexInterval],
        operation: impl FnOnce(
            &[FmInterval],
            &mut [FmInterval],
        ) -> Result<(), CombinedIndexBackendError>,
    ) -> Result<(), CombinedIndexQueryError> {
        const MAX_LANES: usize = 64;
        if intervals.len() != output.len() || intervals.len() > MAX_LANES {
            return Err(CombinedIndexQueryError::Backend(
                CombinedIndexBackendError::Structure,
            ));
        }
        let Some(&first) = intervals.first() else {
            return Ok(());
        };
        let first = self.validate(first)?;
        let mut private = [first; MAX_LANES];
        let mut private_output = [first; MAX_LANES];
        for (destination, &interval) in private.iter_mut().zip(intervals) {
            *destination = self.validate(interval)?;
        }
        operation(
            &private[..intervals.len()],
            &mut private_output[..intervals.len()],
        )
        .map_err(CombinedIndexQueryError::Backend)?;
        for (destination, interval) in output.iter_mut().zip(private_output) {
            *destination = self.bind(interval);
        }
        Ok(())
    }

    /// Recovers a contig-local coordinate from one raw combined-image suffix
    /// position and caller-supplied text placement facts.
    #[must_use]
    pub fn recover_coordinate(
        self,
        position: u64,
        matched_bases: u64,
        pattern_offset: u64,
        pattern_len: u64,
    ) -> Option<CombinedIndexCoordinate> {
        let reference_len = PrivateCombinedIndex::reference_len(self.backend);
        let combined_len = reference_len.checked_mul(2)?;
        let site = combined_len
            .checked_sub(position)?
            .checked_sub(matched_bases)?
            .checked_sub(pattern_offset)?;
        let (strand, global_start) = if site < reference_len {
            (BisulfiteStrand::OT, site)
        } else {
            let start = combined_len.checked_sub(site)?.checked_sub(pattern_len)?;
            (BisulfiteStrand::OB, start)
        };
        let contig_position = self
            .catalog
            .contig_starts
            .partition_point(|&start| start <= global_start)
            .checked_sub(1)?;
        let contig = self.reference.owner.contigs.get(contig_position)?;
        let local_start = global_start - self.catalog.contig_starts[contig_position];
        let local_end = local_start.checked_add(pattern_len)?;
        if local_end > contig.sequence.len() {
            return None;
        }
        Some(CombinedIndexCoordinate {
            contig_ordinal: u64::try_from(contig_position).ok()?,
            strand,
            start: local_start,
        })
    }

    fn bind(self, interval: FmInterval) -> CombinedIndexInterval {
        CombinedIndexInterval {
            owner_token: self.owner_token,
            interval,
        }
    }

    fn validate(
        self,
        interval: CombinedIndexInterval,
    ) -> Result<FmInterval, CombinedIndexQueryError> {
        if interval.owner_token != self.owner_token {
            return Err(CombinedIndexQueryError::ForeignInterval);
        }
        Ok(interval.interval)
    }
}

#[cfg(feature = "combined-index")]
const fn private_locate_metrics(metrics: PrivateCombinedLocateMetrics) -> ReferenceLocateMetrics {
    ReferenceLocateMetrics {
        located_coordinates: metrics.located_rows(),
        lf_steps: metrics.lf_steps(),
        rank_operations: metrics.rank_operations(),
        interval_nodes: metrics.interval_nodes(),
    }
}
