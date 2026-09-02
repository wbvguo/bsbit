//! Owner-bound reference storage, query, and coordinate recovery runtime.

use core::{fmt, mem::size_of};
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::{Arc, OnceLock};

use bsbit_core::alphabet::Base;
use bsbit_core::bisulfite::{
    AlignmentOrientation, BisulfiteStrand, ThreeLetterConversion, strand_semantics,
};
use bsbit_core::coordinate::{ReferenceInterval, ReferenceLength};
use bsbit_core::sequence::NormalizedSequence;

use crate::storage::fm::{FmBuildLimit, FmError, FmIndex, FmInterval, SearchBase};

#[cfg(feature = "combined-index")]
use super::CombinedIndexBackendError;
use super::{
    ContigInput, ReferenceAccessError, ReferenceAllocation, ReferenceArithmetic,
    ReferenceBuildError, ReferenceBuildLimits, ReferenceCatalogLimits, ReferenceLocateError,
    ReferenceLocateInvariant, ReferenceLocateMetrics, ReferenceMetrics, ReferenceQueryCounter,
    ReferenceQueryError, ReferenceQueryLimits, ReferenceResource, validate_catalog_and_measure,
};
#[cfg(feature = "combined-index")]
use super::{
    PrivateCombinedIndex, PrivateCombinedLocateMetrics, PrivateCombinedReference,
    PrivateCombinedReferenceError,
};

#[cfg(feature = "combined-index")]
mod combined_query;
#[cfg(feature = "combined-index")]
pub use combined_query::*;

/// An immutable projected reference index.
#[derive(Clone)]
pub struct ReferenceIndex {
    owner: Arc<ReferenceOwner>,
}

/// Optional runtime work counters for combined-index search and locate calls.
///
/// These counters are disabled by default. They are intended for the aligner's
/// explicit profiling mode and do not change query or locate results.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
#[doc(hidden)]
pub struct ReferenceQueryDiagnostics {
    suffix_search_lanes: u64,
    suffix_search_rank_operations: u64,
    locate_calls: u64,
    singleton_locate_calls: u64,
    multi_hit_locate_calls: u64,
    located_rows: u64,
    locate_lf_steps: u64,
    locate_rank_operations: u64,
    locate_interval_nodes: u64,
}

impl ReferenceQueryDiagnostics {
    /// Returns projected maximal-suffix lanes submitted to the combined index.
    #[must_use]
    pub const fn suffix_search_lanes(self) -> u64 {
        self.suffix_search_lanes
    }

    /// Returns physical rank-boundary operations in maximal-suffix extension.
    #[must_use]
    pub const fn suffix_search_rank_operations(self) -> u64 {
        self.suffix_search_rank_operations
    }

    /// Returns complete sampled-SA interval locate calls.
    #[must_use]
    pub const fn locate_calls(self) -> u64 {
        self.locate_calls
    }

    /// Returns locate calls whose input interval contained one suffix row.
    #[must_use]
    pub const fn singleton_locate_calls(self) -> u64 {
        self.singleton_locate_calls
    }

    /// Returns locate calls whose input interval contained multiple suffix rows.
    #[must_use]
    pub const fn multi_hit_locate_calls(self) -> u64 {
        self.multi_hit_locate_calls
    }

    /// Returns suffix rows completed by sampled-SA locate.
    #[must_use]
    pub const fn located_rows(self) -> u64 {
        self.located_rows
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
}

impl fmt::Debug for ReferenceIndex {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("ReferenceIndex")
            .field("metrics", &self.owner.metrics)
            .finish_non_exhaustive()
    }
}

/// An opaque process-local instance identifier.
#[derive(Clone)]
pub struct ReferenceInstanceId {
    owner: Arc<ReferenceOwner>,
}

impl ReferenceInstanceId {
    /// Reports exact shared runtime ownership.
    #[must_use]
    pub fn is_same_instance(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.owner, &other.owner)
    }
}

impl fmt::Debug for ReferenceInstanceId {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("ReferenceInstanceId")
            .finish_non_exhaustive()
    }
}

/// An opaque owner-bound contig identifier.
#[derive(Clone)]
pub struct ContigId {
    owner: Arc<ReferenceOwner>,
    ordinal: u64,
}

impl ContigId {
    /// Returns the deterministic catalog ordinal.
    #[must_use]
    pub const fn ordinal(&self) -> u64 {
        self.ordinal
    }

    /// Reports exact shared runtime ownership.
    #[must_use]
    pub fn is_same_instance(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.owner, &other.owner)
    }

    /// Reports exact shared ownership and equal local ordinal.
    #[must_use]
    pub fn is_same_contig(&self, other: &Self) -> bool {
        self.is_same_instance(other) && self.ordinal == other.ordinal
    }
}

impl fmt::Debug for ContigId {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("ContigId")
            .field("ordinal", &self.ordinal)
            .finish_non_exhaustive()
    }
}

/// A borrowed view of one retained contig.
#[derive(Clone, Copy, Debug)]
pub struct ContigView<'index> {
    ordinal: u64,
    name: &'index [u8],
    sequence: &'index NormalizedSequence,
}

impl<'index> ContigView<'index> {
    /// Returns the deterministic catalog ordinal.
    #[must_use]
    pub const fn ordinal(self) -> u64 {
        self.ordinal
    }

    /// Returns the exact name bytes.
    #[must_use]
    pub const fn name(self) -> &'index [u8] {
        self.name
    }

    /// Returns the retained original sequence.
    #[must_use]
    pub const fn sequence(self) -> &'index NormalizedSequence {
        self.sequence
    }
}

/// An opaque owner-bound aggregate of exact FM intervals.
pub struct ProjectedMatches {
    owner: Arc<ReferenceOwner>,
    strand: BisulfiteStrand,
    pattern_len: u64,
    exact_hit_count: u64,
    search_rank_operations: u64,
    entries: MatchEntries,
}

enum MatchEntries {
    PerRun(Vec<RunMatch>),
}

impl MatchEntries {
    fn len(&self) -> usize {
        match self {
            Self::PerRun(entries) => entries.len(),
        }
    }
}

impl ProjectedMatches {
    /// Returns the exact aggregate occurrence count.
    #[must_use]
    pub const fn exact_hit_count(&self) -> u64 {
        self.exact_hit_count
    }

    #[doc(hidden)]
    #[must_use]
    pub const fn search_rank_operations(&self) -> u64 {
        self.search_rank_operations
    }

    /// Reports whether another match artifact has the same runtime owner.
    #[must_use]
    pub fn is_same_instance(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.owner, &other.owner)
    }

    /// Reports whether this artifact belongs to an instance identifier.
    #[must_use]
    pub fn belongs_to_instance(&self, instance: &ReferenceInstanceId) -> bool {
        Arc::ptr_eq(&self.owner, &instance.owner)
    }

    /// Returns the number of nonempty run intervals.
    #[must_use]
    pub fn nonempty_interval_count(&self) -> u64 {
        u64::try_from(self.entries.len()).unwrap_or(u64::MAX)
    }

    /// Returns the number of nonempty run intervals.
    #[must_use]
    pub fn matched_interval_count(&self) -> u64 {
        self.nonempty_interval_count()
    }

    /// Returns whether the complete exact result is empty.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.exact_hit_count == 0
    }

    /// Returns the searched bisulfite strand.
    #[must_use]
    pub const fn strand(&self) -> BisulfiteStrand {
        self.strand
    }

    /// Returns the nonzero searched pattern length.
    #[must_use]
    pub const fn pattern_len(&self) -> u64 {
        self.pattern_len
    }
}

impl fmt::Debug for ProjectedMatches {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("ProjectedMatches")
            .field("strand", &self.strand)
            .field("pattern_len", &self.pattern_len)
            .field("exact_hit_count", &self.exact_hit_count)
            .field("search_rank_operations", &self.search_rank_operations)
            .field("nonempty_interval_count", &self.entries.len())
            .finish_non_exhaustive()
    }
}

/// One recovered exact projected-reference hit.
pub struct ProjectedHit {
    contig: ContigId,
    interval: ReferenceInterval,
    strand: BisulfiteStrand,
}

impl ProjectedHit {
    /// Returns the owner-bound contig identifier.
    #[must_use]
    pub const fn contig(&self) -> &ContigId {
        &self.contig
    }

    /// Returns the forward contig-local interval.
    #[must_use]
    pub const fn interval(&self) -> ReferenceInterval {
        self.interval
    }

    /// Returns the bisulfite strand lane.
    #[must_use]
    pub const fn strand(&self) -> BisulfiteStrand {
        self.strand
    }
}

impl fmt::Debug for ProjectedHit {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("ProjectedHit")
            .field("contig_ordinal", &self.contig.ordinal)
            .field("interval", &self.interval)
            .field("strand", &self.strand)
            .finish()
    }
}

/// One 64-base window of the retained reference's neutral two-bit planes.
///
/// Raw bases use `A/N=00`, `C=01`, `G=10`, and `T=11`. Consumers decide how
/// to compare these facts with a query or chemistry projection.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct PackedReferenceWord {
    low: u64,
    high: u64,
}

impl PackedReferenceWord {
    /// Returns the low bit plane for this reference window.
    #[must_use]
    pub const fn low(self) -> u64 {
        self.low
    }

    /// Returns the high bit plane for this reference window.
    #[must_use]
    pub const fn high(self) -> u64 {
        self.high
    }
}

struct ReferenceOwner {
    contigs: Vec<ContigInput>,
    runs: RunCatalog,
    metrics: ReferenceMetrics,
    packed_reference_words: OnceLock<Result<PrivatePackedReferenceWords, String>>,
    query_diagnostics: AtomicReferenceQueryDiagnostics,
}

#[derive(Default)]
struct AtomicReferenceQueryDiagnostics {
    enabled: AtomicBool,
    suffix_search_lanes: AtomicU64,
    suffix_search_rank_operations: AtomicU64,
    locate_calls: AtomicU64,
    singleton_locate_calls: AtomicU64,
    multi_hit_locate_calls: AtomicU64,
    located_rows: AtomicU64,
    locate_lf_steps: AtomicU64,
    locate_rank_operations: AtomicU64,
    locate_interval_nodes: AtomicU64,
}

impl AtomicReferenceQueryDiagnostics {
    fn enable(&self) {
        self.suffix_search_lanes.store(0, Ordering::Relaxed);
        self.suffix_search_rank_operations
            .store(0, Ordering::Relaxed);
        self.locate_calls.store(0, Ordering::Relaxed);
        self.singleton_locate_calls.store(0, Ordering::Relaxed);
        self.multi_hit_locate_calls.store(0, Ordering::Relaxed);
        self.located_rows.store(0, Ordering::Relaxed);
        self.locate_lf_steps.store(0, Ordering::Relaxed);
        self.locate_rank_operations.store(0, Ordering::Relaxed);
        self.locate_interval_nodes.store(0, Ordering::Relaxed);
        self.enabled.store(true, Ordering::Release);
    }

    #[cfg(feature = "combined-index")]
    fn is_enabled(&self) -> bool {
        self.enabled.load(Ordering::Relaxed)
    }

    #[cfg(feature = "combined-index")]
    fn observe_suffix_search(&self, lanes: u64, rank_operations: u64) {
        self.suffix_search_lanes.fetch_add(lanes, Ordering::Relaxed);
        self.suffix_search_rank_operations
            .fetch_add(rank_operations, Ordering::Relaxed);
    }

    #[cfg(feature = "combined-index")]
    fn observe_locate(&self, interval_rows: u64, metrics: ReferenceLocateMetrics) {
        self.locate_calls.fetch_add(1, Ordering::Relaxed);
        if interval_rows == 1 {
            self.singleton_locate_calls.fetch_add(1, Ordering::Relaxed);
        } else {
            self.multi_hit_locate_calls.fetch_add(1, Ordering::Relaxed);
        }
        self.located_rows
            .fetch_add(metrics.located_coordinates(), Ordering::Relaxed);
        self.locate_lf_steps
            .fetch_add(metrics.lf_steps(), Ordering::Relaxed);
        self.locate_rank_operations
            .fetch_add(metrics.rank_operations(), Ordering::Relaxed);
        self.locate_interval_nodes
            .fetch_add(metrics.interval_nodes(), Ordering::Relaxed);
    }

    fn disable_and_snapshot(&self) -> ReferenceQueryDiagnostics {
        self.enabled.store(false, Ordering::Release);
        ReferenceQueryDiagnostics {
            suffix_search_lanes: self.suffix_search_lanes.load(Ordering::Relaxed),
            suffix_search_rank_operations: self
                .suffix_search_rank_operations
                .load(Ordering::Relaxed),
            locate_calls: self.locate_calls.load(Ordering::Relaxed),
            singleton_locate_calls: self.singleton_locate_calls.load(Ordering::Relaxed),
            multi_hit_locate_calls: self.multi_hit_locate_calls.load(Ordering::Relaxed),
            located_rows: self.located_rows.load(Ordering::Relaxed),
            locate_lf_steps: self.locate_lf_steps.load(Ordering::Relaxed),
            locate_rank_operations: self.locate_rank_operations.load(Ordering::Relaxed),
            locate_interval_nodes: self.locate_interval_nodes.load(Ordering::Relaxed),
        }
    }
}

/// Private two-bit positional representation of retained reference bases.
///
/// Raw bases use `A/N=00`, `C=01`, `G=10`, and `T=11`.  Both bisulfite
/// projections can therefore be derived from two words without retaining four
/// strand-specific one-hot planes.
struct PrivatePackedReferenceWords {
    contigs: Vec<PrivatePackedContigWords>,
}

struct PrivatePackedContigWords {
    length: usize,
    low: Vec<u64>,
    high: Vec<u64>,
}

impl PrivatePackedReferenceWords {
    fn build(contigs: &[ContigInput]) -> Result<Self, String> {
        let mut encoded = Vec::new();
        encoded
            .try_reserve_exact(contigs.len())
            .map_err(|_| "packed reference-word catalog allocation failed".to_owned())?;
        for contig in contigs {
            let bases = contig.sequence().bases();
            let words = bases
                .len()
                .checked_add(63)
                .ok_or_else(|| "packed reference-word contig length overflow".to_owned())?
                / 64;
            let mut low = Vec::new();
            let mut high = Vec::new();
            low.try_reserve_exact(words)
                .map_err(|_| "packed reference low-plane allocation failed".to_owned())?;
            high.try_reserve_exact(words)
                .map_err(|_| "packed reference high-plane allocation failed".to_owned())?;
            low.resize(words, 0);
            high.resize(words, 0);
            for (position, &base) in bases.iter().enumerate() {
                let bit = 1_u64 << (position % 64);
                match base {
                    Base::C => low[position / 64] |= bit,
                    Base::G => high[position / 64] |= bit,
                    Base::T => {
                        low[position / 64] |= bit;
                        high[position / 64] |= bit;
                    }
                    _ => {}
                }
            }
            encoded.push(PrivatePackedContigWords {
                length: bases.len(),
                low,
                high,
            });
        }
        Ok(Self { contigs: encoded })
    }
}

impl PrivatePackedContigWords {
    fn shifted_word(plane: &[u64], start: usize) -> u64 {
        let word = start / 64;
        let shift = start % 64;
        let lower = plane.get(word).copied().unwrap_or(0) >> shift;
        if shift == 0 {
            lower
        } else {
            lower | (plane.get(word.saturating_add(1)).copied().unwrap_or(0) << (64 - shift))
        }
    }

    fn word(&self, start: usize) -> Option<PackedReferenceWord> {
        if start >= self.length {
            return None;
        }
        Some(PackedReferenceWord {
            low: Self::shifted_word(&self.low, start),
            high: Self::shifted_word(&self.high, start),
        })
    }
}

struct RunMetadata {
    contig_ordinal: u64,
    start: u64,
    end: u64,
    lanes: RunLanes,
}

struct RunLanes {
    ot: FmIndex,
    ob: FmIndex,
    ctot: FmIndex,
    ctob: FmIndex,
}

impl RunLanes {
    const fn get(&self, strand: BisulfiteStrand) -> &FmIndex {
        match strand {
            BisulfiteStrand::OT => &self.ot,
            BisulfiteStrand::OB => &self.ob,
            BisulfiteStrand::CTOT => &self.ctot,
            BisulfiteStrand::CTOB => &self.ctob,
        }
    }
}

#[cfg(feature = "combined-index")]
struct CombinedReferenceCatalog {
    runs: Vec<RunCoordinates>,
    contig_starts: Vec<u64>,
    combined_index: Box<dyn PrivateCombinedIndex>,
}

enum RunCatalog {
    Scalar(Vec<RunMetadata>),
    #[cfg(feature = "combined-index")]
    Combined(CombinedReferenceCatalog),
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
struct RunLocateMetrics {
    located_rows: u64,
    lf_steps: u64,
    rank_operations: u64,
    interval_nodes: u64,
}

impl RunCatalog {
    fn len(&self) -> usize {
        match self {
            Self::Scalar(runs) => runs.len(),
            #[cfg(feature = "combined-index")]
            Self::Combined(catalog) => catalog.runs.len(),
        }
    }

    fn coordinates(&self, run_index: usize) -> Option<RunCoordinates> {
        match self {
            Self::Scalar(runs) => runs.get(run_index).map(|run| RunCoordinates {
                contig_ordinal: run.contig_ordinal,
                start: run.start,
                end: run.end,
            }),
            #[cfg(feature = "combined-index")]
            Self::Combined(catalog) => catalog.runs.get(run_index).copied(),
        }
    }

    fn exact_search(
        &self,
        run_index: usize,
        strand: BisulfiteStrand,
        pattern: &[SearchBase],
    ) -> Option<FmInterval> {
        match self {
            Self::Scalar(runs) => {
                Some(runs.get(run_index)?.lanes.get(strand).exact_search(pattern))
            }
            #[cfg(feature = "combined-index")]
            Self::Combined(_) => None,
        }
    }

    fn visit_locate(
        &self,
        run_index: usize,
        strand: BisulfiteStrand,
        interval: FmInterval,
        visitor: &mut dyn FnMut(u64) -> bool,
    ) -> Option<Result<RunLocateMetrics, FmError>> {
        match self {
            Self::Scalar(runs) => Some(
                runs.get(run_index)?
                    .lanes
                    .get(strand)
                    .visit_locate(interval, visitor)
                    .map(|located_rows| RunLocateMetrics {
                        located_rows,
                        ..RunLocateMetrics::default()
                    }),
            ),
            #[cfg(feature = "combined-index")]
            Self::Combined(_) => None,
        }
    }
}

#[derive(Clone, Copy)]
struct RunCoordinates {
    contig_ordinal: u64,
    start: u64,
    end: u64,
}

impl RunCoordinates {
    const fn len(self) -> u64 {
        self.end - self.start
    }
}

struct RunMatch {
    run_index: usize,
    interval: FmInterval,
}

impl ReferenceIndex {
    /// Builds one complete immutable projected reference.
    ///
    /// Validation and resource errors are returned before a partial index can
    /// be published.
    ///
    /// # Errors
    ///
    /// Returns `ReferenceBuildError` under the deterministic validation, resource,
    /// allocation, or contextual FM-build failure order.
    #[allow(clippy::too_many_lines)]
    pub fn build(
        contigs: Vec<ContigInput>,
        limits: ReferenceBuildLimits,
    ) -> Result<Self, ReferenceBuildError> {
        let scan = validate_and_measure(&contigs, limits)?;
        let run_storage = preflight_build_allocation::<RunMetadata>(
            ReferenceAllocation::RunMetadata,
            scan.metrics.canonical_run_count,
        )?;
        let scratch_storage = preflight_build_allocation::<SearchBase>(
            ReferenceAllocation::ProjectionScratch,
            scan.max_run_bases,
        )?;

        let mut runs = Vec::new();
        runs.try_reserve_exact(run_storage)
            .map_err(|_| ReferenceBuildError::AllocationFailed {
                allocation: ReferenceAllocation::RunMetadata,
                elements: scan.metrics.canonical_run_count,
            })?;
        let mut scratch = Vec::new();
        scratch.try_reserve_exact(scratch_storage).map_err(|_| {
            ReferenceBuildError::AllocationFailed {
                allocation: ReferenceAllocation::ProjectionScratch,
                elements: scan.max_run_bases,
            }
        })?;

        for (contig_storage, contig) in contigs.iter().enumerate() {
            let contig_ordinal = physical_to_logical(contig_storage, ReferenceResource::Contigs)?;
            for (start_storage, end_storage) in CanonicalRuns::new(contig.sequence.bases()) {
                ensure_build_capacity(scan.metrics.canonical_run_count, run_storage, runs.len())?;
                let start =
                    physical_to_logical(start_storage, ReferenceResource::TotalReferenceBases)?;
                let end = physical_to_logical(end_storage, ReferenceResource::TotalReferenceBases)?;
                let original = &contig.sequence.bases()[start_storage..end_storage];
                let ot = build_lane(
                    original,
                    BisulfiteStrand::OT,
                    contig_ordinal,
                    start,
                    limits.max_suffix_rows_per_lane,
                    &mut scratch,
                    scratch_storage,
                )?;
                let ob = build_lane(
                    original,
                    BisulfiteStrand::OB,
                    contig_ordinal,
                    start,
                    limits.max_suffix_rows_per_lane,
                    &mut scratch,
                    scratch_storage,
                )?;
                let complementary_top_lane = build_lane(
                    original,
                    BisulfiteStrand::CTOT,
                    contig_ordinal,
                    start,
                    limits.max_suffix_rows_per_lane,
                    &mut scratch,
                    scratch_storage,
                )?;
                let complementary_bottom_lane = build_lane(
                    original,
                    BisulfiteStrand::CTOB,
                    contig_ordinal,
                    start,
                    limits.max_suffix_rows_per_lane,
                    &mut scratch,
                    scratch_storage,
                )?;
                runs.push(RunMetadata {
                    contig_ordinal,
                    start,
                    end,
                    lanes: RunLanes {
                        ot,
                        ob,
                        ctot: complementary_top_lane,
                        ctob: complementary_bottom_lane,
                    },
                });
            }
        }

        let observed_runs = physical_to_logical(runs.len(), ReferenceResource::CanonicalRuns)?;
        if observed_runs != scan.metrics.canonical_run_count {
            return Err(ReferenceBuildError::InternalInvariant {
                expected: scan.metrics.canonical_run_count,
                observed: observed_runs,
            });
        }

        let owner = Arc::new(ReferenceOwner {
            contigs,
            runs: RunCatalog::Scalar(runs),
            metrics: scan.metrics,
            packed_reference_words: OnceLock::new(),
            query_diagnostics: AtomicReferenceQueryDiagnostics::default(),
        });
        Ok(Self { owner })
    }

    /// Assembles an owner backed only by the frozen combined index.
    ///
    /// Combined coordinates come from that index while the same owner retains
    /// original contig sequences. Omitting regular lanes keeps their object
    /// validation, page residency, and retained bytes out of this owner.
    ///
    /// # Errors
    ///
    /// Returns catalog or combined-index dimension failures.
    #[cfg(feature = "combined-index")]
    #[doc(hidden)]
    pub(crate) fn from_private_combined(
        contigs: Vec<ContigInput>,
        combined: PrivateCombinedReference,
    ) -> Result<Self, PrivateCombinedReferenceError> {
        let scan = validate_and_measure(&contigs, ReferenceBuildLimits::MAX)
            .map_err(|source| PrivateCombinedReferenceError::Catalog { source })?;
        if combined.index.reference_len() != scan.metrics.total_reference_bases {
            return Err(PrivateCombinedReferenceError::CombinedDimensions {
                expected_reference_len: scan.metrics.total_reference_bases,
                observed_reference_len: combined.index.reference_len(),
            });
        }
        let run_storage = preflight_build_allocation::<RunCoordinates>(
            ReferenceAllocation::RunMetadata,
            scan.metrics.canonical_run_count,
        )
        .map_err(|source| PrivateCombinedReferenceError::Catalog { source })?;
        let mut runs = Vec::new();
        runs.try_reserve_exact(run_storage).map_err(|_| {
            PrivateCombinedReferenceError::Catalog {
                source: ReferenceBuildError::AllocationFailed {
                    allocation: ReferenceAllocation::RunMetadata,
                    elements: scan.metrics.canonical_run_count,
                },
            }
        })?;
        let contig_storage = preflight_build_allocation::<u64>(
            ReferenceAllocation::RunMetadata,
            scan.metrics.contig_count,
        )
        .map_err(|source| PrivateCombinedReferenceError::Catalog { source })?;
        let mut contig_starts = Vec::new();
        contig_starts
            .try_reserve_exact(contig_storage)
            .map_err(|_| PrivateCombinedReferenceError::Catalog {
                source: ReferenceBuildError::AllocationFailed {
                    allocation: ReferenceAllocation::RunMetadata,
                    elements: scan.metrics.contig_count,
                },
            })?;
        let mut concatenated_start = 0_u64;
        for (contig_storage, contig) in contigs.iter().enumerate() {
            contig_starts.push(concatenated_start);
            concatenated_start = concatenated_start
                .checked_add(contig.sequence.len())
                .ok_or(PrivateCombinedReferenceError::Catalog {
                    source: ReferenceBuildError::ArithmeticOverflow {
                        resource: ReferenceResource::TotalReferenceBases,
                        operation: ReferenceArithmetic::Add,
                        lhs: concatenated_start,
                        rhs: contig.sequence.len(),
                    },
                })?;
            let contig_ordinal = physical_to_logical(contig_storage, ReferenceResource::Contigs)
                .map_err(|source| PrivateCombinedReferenceError::Catalog { source })?;
            for (start_storage, end_storage) in CanonicalRuns::new(contig.sequence.bases()) {
                runs.push(RunCoordinates {
                    contig_ordinal,
                    start: physical_to_logical(
                        start_storage,
                        ReferenceResource::TotalReferenceBases,
                    )
                    .map_err(|source| PrivateCombinedReferenceError::Catalog { source })?,
                    end: physical_to_logical(end_storage, ReferenceResource::TotalReferenceBases)
                        .map_err(|source| PrivateCombinedReferenceError::Catalog { source })?,
                });
            }
        }
        let owner = Arc::new(ReferenceOwner {
            contigs,
            runs: RunCatalog::Combined(CombinedReferenceCatalog {
                runs,
                contig_starts,
                combined_index: combined.index,
            }),
            metrics: scan.metrics,
            packed_reference_words: OnceLock::new(),
            query_diagnostics: AtomicReferenceQueryDiagnostics::default(),
        });
        Ok(Self { owner })
    }

    /// Returns an opaque identifier for this exact runtime instance.
    #[must_use]
    pub fn instance_id(&self) -> ReferenceInstanceId {
        ReferenceInstanceId {
            owner: Arc::clone(&self.owner),
        }
    }

    /// Returns the number of retained contigs.
    #[must_use]
    pub fn contig_count(&self) -> u64 {
        self.owner.metrics.contig_count
    }

    /// Returns complete deterministic build metrics.
    #[must_use]
    pub fn metrics(&self) -> ReferenceMetrics {
        self.owner.metrics
    }

    /// Enables and resets optional combined-index work diagnostics.
    #[doc(hidden)]
    pub fn enable_query_diagnostics(&self) {
        self.owner.query_diagnostics.enable();
    }

    /// Disables optional work diagnostics and returns their final snapshot.
    #[doc(hidden)]
    #[must_use]
    pub fn disable_and_take_query_diagnostics(&self) -> ReferenceQueryDiagnostics {
        self.owner.query_diagnostics.disable_and_snapshot()
    }

    /// Borrows the raw interval, rank, and locate surface for the combined
    /// image.
    ///
    /// The returned handle owns no seeding or alignment policy. It is absent
    /// for ordinary in-memory reference owners.
    #[cfg(feature = "combined-index")]
    #[must_use]
    pub fn combined_index_query(&self) -> Option<CombinedIndexQuery<'_>> {
        let RunCatalog::Combined(catalog) = &self.owner.runs else {
            return None;
        };
        Some(CombinedIndexQuery {
            reference: self,
            catalog,
            backend: catalog.combined_index.as_ref(),
            owner_token: Arc::as_ptr(&self.owner) as usize,
        })
    }

    fn packed_reference_words(&self) -> Result<&PrivatePackedReferenceWords, String> {
        match self
            .owner
            .packed_reference_words
            .get_or_init(|| PrivatePackedReferenceWords::build(&self.owner.contigs))
        {
            Ok(planes) => Ok(planes),
            Err(error) => Err(error.clone()),
        }
    }

    /// Returns one caller-selected 64-base window of neutral reference facts.
    #[doc(hidden)]
    #[must_use]
    pub fn packed_reference_word(
        &self,
        contig_ordinal: u64,
        start: usize,
    ) -> Option<PackedReferenceWord> {
        let storage = usize::try_from(contig_ordinal).ok()?;
        self.packed_reference_words()
            .ok()?
            .contigs
            .get(storage)?
            .word(start)
    }

    /// Creates an owner-bound contig identifier after bounds validation.
    ///
    /// # Errors
    ///
    /// Returns `ReferenceAccessError` when the ordinal is outside this catalog.
    pub fn contig_id(&self, ordinal: u64) -> Result<ContigId, ReferenceAccessError> {
        let _ = contig_storage(ordinal, self.owner.metrics.contig_count)?;
        Ok(ContigId {
            owner: Arc::clone(&self.owner),
            ordinal,
        })
    }

    /// Resolves an owner-bound contig identifier to retained content.
    ///
    /// Exact owner identity is checked before ordinal bounds.
    ///
    /// # Errors
    ///
    /// Returns `ReferenceAccessError` for a foreign identifier or invalid ordinal.
    pub fn resolve_contig<'index>(
        &'index self,
        contig: &ContigId,
    ) -> Result<ContigView<'index>, ReferenceAccessError> {
        if !Arc::ptr_eq(&self.owner, &contig.owner) {
            return Err(ReferenceAccessError::ForeignContigId);
        }
        let storage = contig_storage(contig.ordinal, self.owner.metrics.contig_count)?;
        let retained = &self.owner.contigs[storage];
        Ok(ContigView {
            ordinal: contig.ordinal,
            name: retained.name(),
            sequence: retained.sequence(),
        })
    }

    /// Borrows retained contig content by ordinal without creating an
    /// owner-bound identifier.
    ///
    /// This is an integer-only coordinate boundary.
    /// Ordinals reaching it have already been decoded from this index's
    /// validated packed lanes.
    #[doc(hidden)]
    #[must_use]
    pub fn contig_by_ordinal(&self, ordinal: u64) -> Option<ContigView<'_>> {
        let storage = usize::try_from(ordinal).ok()?;
        let retained = self.owner.contigs.get(storage)?;
        Some(ContigView {
            ordinal,
            name: retained.name(),
            sequence: retained.sequence(),
        })
    }

    /// Searches every canonical run in one bisulfite lane.
    ///
    /// The supplied pattern remains in sequencing order. The complete pattern
    /// is projected once, then searched against each canonical run.
    ///
    /// # Errors
    ///
    /// Returns `ReferenceQueryError` for invalid input, limits, allocation failure,
    /// count overflow, or a count/materialization disagreement.
    #[allow(clippy::too_many_lines)]
    pub fn exact_search(
        &self,
        strand: BisulfiteStrand,
        pattern: &[Base],
        limits: ReferenceQueryLimits,
    ) -> Result<ProjectedMatches, ReferenceQueryError> {
        let pattern_len = query_storage_to_logical(pattern.len())?;
        if pattern_len == 0 {
            return Err(ReferenceQueryError::EmptyPattern);
        }
        if pattern_len > limits.max_pattern_bases {
            return Err(ReferenceQueryError::PatternLimitExceeded {
                requested: pattern_len,
                maximum: limits.max_pattern_bases,
            });
        }
        for (storage, &base) in pattern.iter().enumerate() {
            if base == Base::N {
                let offset = u64::try_from(storage).unwrap_or(u64::MAX);
                return Err(ReferenceQueryError::UnsearchableBase { offset });
            }
        }

        let projected_storage = preflight_query_allocation::<SearchBase>(
            ReferenceAllocation::ProjectedPattern,
            pattern_len,
        )?;
        let mut projected = Vec::new();
        projected
            .try_reserve_exact(projected_storage)
            .map_err(|_| ReferenceQueryError::AllocationFailed {
                allocation: ReferenceAllocation::ProjectedPattern,
                elements: pattern_len,
            })?;
        let conversion = raw_view_conversion(strand);
        for &base in pattern {
            ensure_query_capacity(
                pattern_len,
                u64::try_from(projected.len()).unwrap_or(u64::MAX),
            )?;
            let converted = convert_base(base, conversion);
            let Some(search) = SearchBase::from_base(converted) else {
                return Err(ReferenceQueryError::UnsearchableBase {
                    offset: u64::try_from(projected.len()).unwrap_or(u64::MAX),
                });
            };
            projected.push(search);
        }

        #[cfg(feature = "combined-index")]
        if matches!(&self.owner.runs, RunCatalog::Combined(_)) {
            return Ok(ProjectedMatches {
                owner: Arc::clone(&self.owner),
                strand,
                pattern_len,
                exact_hit_count: 0,
                search_rank_operations: 0,
                entries: MatchEntries::PerRun(Vec::new()),
            });
        }

        let (hit_count, entry_count, mut search_rank_operations) =
            self.count_matches(strand, &projected)?;
        if hit_count > limits.max_exact_hits {
            return Err(ReferenceQueryError::HitLimitExceeded {
                requested: hit_count,
                maximum: limits.max_exact_hits,
            });
        }
        let entry_storage = preflight_query_allocation::<RunMatch>(
            ReferenceAllocation::OpaqueMatches,
            entry_count,
        )?;
        let mut entries = Vec::new();
        entries.try_reserve_exact(entry_storage).map_err(|_| {
            ReferenceQueryError::AllocationFailed {
                allocation: ReferenceAllocation::OpaqueMatches,
                elements: entry_count,
            }
        })?;

        let mut observed_hits = 0_u64;
        let mut observed_entries = 0_u64;
        for run_index in 0..self.owner.runs.len() {
            let interval = self
                .owner
                .runs
                .exact_search(run_index, strand, &projected)
                .ok_or(ReferenceQueryError::InvariantMismatch {
                    counter: ReferenceQueryCounter::NonemptyIntervals,
                    expected: entry_count,
                    observed: observed_entries,
                })?;
            search_rank_operations = checked_query_add(
                ReferenceQueryCounter::RankOperations,
                search_rank_operations,
                pattern_len,
            )?;
            search_rank_operations = checked_query_add(
                ReferenceQueryCounter::RankOperations,
                search_rank_operations,
                pattern_len,
            )?;
            let count = interval.len();
            observed_hits =
                checked_query_add(ReferenceQueryCounter::ExactHits, observed_hits, count)?;
            if !interval.is_empty() {
                observed_entries = checked_query_add(
                    ReferenceQueryCounter::NonemptyIntervals,
                    observed_entries,
                    1,
                )?;
                ensure_query_capacity(entry_count, observed_entries - 1)?;
                entries.push(RunMatch {
                    run_index,
                    interval,
                });
            }
        }
        ensure_query_count(ReferenceQueryCounter::ExactHits, hit_count, observed_hits)?;
        ensure_query_count(
            ReferenceQueryCounter::NonemptyIntervals,
            entry_count,
            observed_entries,
        )?;

        Ok(ProjectedMatches {
            owner: Arc::clone(&self.owner),
            strand,
            pattern_len,
            exact_hit_count: hit_count,
            search_rank_operations,
            entries: MatchEntries::PerRun(entries),
        })
    }

    /// Locates and recovers a complete owner-bound match set.
    ///
    /// Foreign ownership is rejected before bounds checks or allocation.
    ///
    /// # Errors
    ///
    /// Returns `ReferenceLocateError` for foreign ownership, allocation or private
    /// FM failure, or a checked coordinate-recovery invariant failure.
    pub fn locate(
        &self,
        matches: &ProjectedMatches,
    ) -> Result<Vec<ProjectedHit>, ReferenceLocateError> {
        if !Arc::ptr_eq(&self.owner, &matches.owner) {
            return Err(ReferenceLocateError::ForeignMatches);
        }
        let hit_storage = preflight_locate_allocation::<ProjectedHit>(
            ReferenceAllocation::FinalHits,
            matches.exact_hit_count,
        )?;
        let mut hits = Vec::new();
        hits.try_reserve_exact(hit_storage).map_err(|_| {
            ReferenceLocateError::AllocationFailed {
                allocation: ReferenceAllocation::FinalHits,
                elements: matches.exact_hit_count,
            }
        })?;

        let metrics = self.visit_located_matches(matches, &mut |hit| {
            hits.push(hit);
            true
        })?;
        ensure_locate_equal(
            ReferenceLocateInvariant::FinalHitCount,
            matches.exact_hit_count,
            metrics.located_coordinates,
        )?;
        hits.sort_unstable_by_key(|hit| {
            (
                hit.contig.ordinal,
                hit.interval.start(),
                hit.interval.end(),
                strand_rank(hit.strand),
            )
        });
        Ok(hits)
    }

    #[allow(clippy::too_many_lines)]
    #[doc(hidden)]
    pub fn visit_located_matches(
        &self,
        matches: &ProjectedMatches,
        visitor: &mut dyn FnMut(ProjectedHit) -> bool,
    ) -> Result<ReferenceLocateMetrics, ReferenceLocateError> {
        if !Arc::ptr_eq(&self.owner, &matches.owner) {
            return Err(ReferenceLocateError::ForeignMatches);
        }
        let mut totals = ReferenceLocateMetrics::default();
        let mut recovered = 0_u64;
        let mut stopped = false;

        match &matches.entries {
            MatchEntries::PerRun(entries) => {
                for entry in entries {
                    let Some(run) = self.owner.runs.coordinates(entry.run_index) else {
                        return Err(ReferenceLocateError::Invariant {
                            invariant: ReferenceLocateInvariant::MissingRun,
                            expected: u64::try_from(self.owner.runs.len()).unwrap_or(u64::MAX),
                            observed: u64::try_from(entry.run_index).unwrap_or(u64::MAX),
                        });
                    };
                    let mut visitor_error = None;
                    let mut visit_offset = |offset| {
                        if let Err(error) =
                            ensure_locate_capacity(matches.exact_hit_count, recovered)
                        {
                            visitor_error = Some(error);
                            return false;
                        }
                        let hit = match recover_projected_hit(self, matches, run, offset) {
                            Ok(hit) => hit,
                            Err(error) => {
                                visitor_error = Some(error);
                                return false;
                            }
                        };
                        recovered += 1;
                        if !visitor(hit) {
                            stopped = true;
                            return false;
                        }
                        true
                    };
                    let metrics = self
                        .owner
                        .runs
                        .visit_locate(
                            entry.run_index,
                            matches.strand,
                            entry.interval,
                            &mut visit_offset,
                        )
                        .ok_or(ReferenceLocateError::Invariant {
                            invariant: ReferenceLocateInvariant::MissingRun,
                            expected: u64::try_from(self.owner.runs.len()).unwrap_or(u64::MAX),
                            observed: u64::try_from(entry.run_index).unwrap_or(u64::MAX),
                        })?
                        .map_err(|source| ReferenceLocateError::FmLocate {
                            contig_ordinal: run.contig_ordinal,
                            run_start: run.start,
                            strand: matches.strand,
                            source,
                        })?;
                    if let Some(error) = visitor_error {
                        return Err(error);
                    }
                    if !stopped {
                        ensure_locate_equal(
                            ReferenceLocateInvariant::OffsetCount,
                            entry.interval.len(),
                            metrics.located_rows,
                        )?;
                    }
                    accumulate_locate_metrics(
                        &mut totals,
                        metrics.located_rows,
                        metrics.lf_steps,
                        metrics.rank_operations,
                        metrics.interval_nodes,
                    )?;
                    if stopped {
                        break;
                    }
                }
            }
        }

        if !stopped {
            ensure_locate_equal(
                ReferenceLocateInvariant::FinalHitCount,
                matches.exact_hit_count,
                recovered,
            )?;
        }
        ensure_locate_equal(
            ReferenceLocateInvariant::OffsetCount,
            recovered,
            totals.located_coordinates,
        )?;
        Ok(totals)
    }

    fn count_matches(
        &self,
        strand: BisulfiteStrand,
        projected: &[SearchBase],
    ) -> Result<(u64, u64, u64), ReferenceQueryError> {
        let mut hits = 0_u64;
        let mut entries = 0_u64;
        let pattern_len = u64::try_from(projected.len()).map_err(|_| {
            ReferenceQueryError::PatternLengthNotRepresentable {
                pattern_len: projected.len(),
            }
        })?;
        let mut rank_operations = 0_u64;
        for run_index in 0..self.owner.runs.len() {
            let interval = self
                .owner
                .runs
                .exact_search(run_index, strand, projected)
                .ok_or(ReferenceQueryError::InvariantMismatch {
                    counter: ReferenceQueryCounter::NonemptyIntervals,
                    expected: entries,
                    observed: entries.saturating_add(1),
                })?;
            rank_operations = checked_query_add(
                ReferenceQueryCounter::RankOperations,
                rank_operations,
                pattern_len,
            )?;
            rank_operations = checked_query_add(
                ReferenceQueryCounter::RankOperations,
                rank_operations,
                pattern_len,
            )?;
            hits = checked_query_add(ReferenceQueryCounter::ExactHits, hits, interval.len())?;
            if !interval.is_empty() {
                entries = checked_query_add(ReferenceQueryCounter::NonemptyIntervals, entries, 1)?;
            }
        }
        Ok((hits, entries, rank_operations))
    }
}

fn recover_projected_hit(
    index: &ReferenceIndex,
    matches: &ProjectedMatches,
    run: RunCoordinates,
    offset: u64,
) -> Result<ProjectedHit, ReferenceLocateError> {
    let interval = recover_interval(
        run,
        offset,
        matches.pattern_len,
        matches.strand,
        &index.owner.contigs,
    )?;
    Ok(ProjectedHit {
        contig: ContigId {
            owner: Arc::clone(&index.owner),
            ordinal: run.contig_ordinal,
        },
        interval,
        strand: matches.strand,
    })
}

fn accumulate_locate_metrics(
    totals: &mut ReferenceLocateMetrics,
    located_coordinates: u64,
    lf_steps: u64,
    rank_operations: u64,
    interval_nodes: u64,
) -> Result<(), ReferenceLocateError> {
    totals.located_coordinates =
        checked_locate_metric_add(totals.located_coordinates, located_coordinates)?;
    totals.lf_steps = checked_locate_metric_add(totals.lf_steps, lf_steps)?;
    totals.rank_operations = checked_locate_metric_add(totals.rank_operations, rank_operations)?;
    totals.interval_nodes = checked_locate_metric_add(totals.interval_nodes, interval_nodes)?;
    Ok(())
}

const fn checked_locate_metric_add(
    accumulated: u64,
    next: u64,
) -> Result<u64, ReferenceLocateError> {
    match accumulated.checked_add(next) {
        Some(value) => Ok(value),
        None => Err(ReferenceLocateError::Invariant {
            invariant: ReferenceLocateInvariant::MetricOverflow,
            expected: accumulated,
            observed: next,
        }),
    }
}

struct ValidationScan {
    metrics: ReferenceMetrics,
    max_run_bases: u64,
}

#[allow(clippy::too_many_lines)]
fn validate_and_measure(
    contigs: &[ContigInput],
    limits: ReferenceBuildLimits,
) -> Result<ValidationScan, ReferenceBuildError> {
    let catalog = validate_catalog_and_measure(
        contigs,
        ReferenceCatalogLimits::MAX
            .with_max_contigs(limits.max_contigs)
            .with_max_total_name_bytes(limits.max_total_name_bytes)
            .with_max_total_reference_bases(limits.max_total_reference_bases),
    )?;
    let contig_count = catalog.contig_count();
    let total_name_bytes = catalog.total_name_bytes();
    let total_reference_bases = catalog.total_reference_bases();

    let mut canonical_bases = 0_u64;
    let mut canonical_run_count = 0_u64;
    let mut max_suffix_rows = 0_u64;
    let mut max_run_bases = 0_u64;
    let mut max_context = (0_u64, 0_u64);
    for (contig_storage, contig) in contigs.iter().enumerate() {
        let contig_ordinal = physical_to_logical(contig_storage, ReferenceResource::Contigs)?;
        for (start_storage, end_storage) in CanonicalRuns::new(contig.sequence.bases()) {
            let start = physical_to_logical(start_storage, ReferenceResource::TotalReferenceBases)?;
            let run_bases = physical_to_logical(
                end_storage - start_storage,
                ReferenceResource::CanonicalBases,
            )?;
            let suffix_rows =
                checked_build_add(ReferenceResource::SuffixRowsPerLane, run_bases, 1)?;
            canonical_bases = checked_build_add(
                ReferenceResource::CanonicalBases,
                canonical_bases,
                run_bases,
            )?;
            canonical_run_count =
                checked_build_add(ReferenceResource::CanonicalRuns, canonical_run_count, 1)?;
            if suffix_rows > max_suffix_rows {
                max_suffix_rows = suffix_rows;
                max_run_bases = run_bases;
                max_context = (contig_ordinal, start);
            }
        }
    }

    apply_limit(
        ReferenceResource::CanonicalRuns,
        canonical_run_count,
        limits.max_canonical_runs,
    )?;
    if max_suffix_rows > limits.max_suffix_rows_per_lane {
        return Err(ReferenceBuildError::SuffixRowsPerLaneLimitExceeded {
            requested: max_suffix_rows,
            maximum: limits.max_suffix_rows_per_lane,
            contig_ordinal: max_context.0,
            run_start: max_context.1,
        });
    }

    let lane_count = checked_build_mul(ReferenceResource::Lanes, canonical_run_count, 4)?;
    apply_limit(ReferenceResource::Lanes, lane_count, limits.max_lanes)?;

    let projected_bases = checked_build_mul(ReferenceResource::ProjectedBases, canonical_bases, 4)?;
    apply_limit(
        ReferenceResource::ProjectedBases,
        projected_bases,
        limits.max_projected_bases,
    )?;

    let bases_plus_runs = checked_build_add(
        ReferenceResource::ProjectedSuffixRows,
        canonical_bases,
        canonical_run_count,
    )?;
    let projected_suffix_rows =
        checked_build_mul(ReferenceResource::ProjectedSuffixRows, bases_plus_runs, 4)?;
    apply_limit(
        ReferenceResource::ProjectedSuffixRows,
        projected_suffix_rows,
        limits.max_projected_suffix_rows,
    )?;

    let estimated_retained_fm_bytes =
        retained_fm_bytes(canonical_bases, canonical_run_count, lane_count)?;
    apply_limit(
        ReferenceResource::EstimatedRetainedFmBytes,
        estimated_retained_fm_bytes,
        limits.max_estimated_retained_fm_bytes,
    )?;

    Ok(ValidationScan {
        metrics: ReferenceMetrics {
            contig_count,
            total_name_bytes,
            total_reference_bases,
            canonical_bases,
            canonical_run_count,
            lane_count,
            projected_bases,
            projected_suffix_rows,
            estimated_retained_fm_bytes,
        },
        max_run_bases,
    })
}

fn retained_fm_bytes(
    canonical_bases: u64,
    canonical_runs: u64,
    lane_count: u64,
) -> Result<u64, ReferenceBuildError> {
    let fm_size = physical_to_logical(
        size_of::<FmIndex>(),
        ReferenceResource::EstimatedRetainedFmBytes,
    )?;
    let term_one = checked_build_mul(
        ReferenceResource::EstimatedRetainedFmBytes,
        lane_count,
        fm_size,
    )?;

    let bases_plus_runs = checked_build_add(
        ReferenceResource::EstimatedRetainedFmBytes,
        canonical_bases,
        canonical_runs,
    )?;
    let four_bases_plus_runs = checked_build_mul(
        ReferenceResource::EstimatedRetainedFmBytes,
        bases_plus_runs,
        4,
    )?;
    let usize_size = physical_to_logical(
        size_of::<usize>(),
        ReferenceResource::EstimatedRetainedFmBytes,
    )?;
    let u8_size =
        physical_to_logical(size_of::<u8>(), ReferenceResource::EstimatedRetainedFmBytes)?;
    let suffix_and_bwt_width = checked_build_add(
        ReferenceResource::EstimatedRetainedFmBytes,
        usize_size,
        u8_size,
    )?;
    let term_two = checked_build_mul(
        ReferenceResource::EstimatedRetainedFmBytes,
        four_bases_plus_runs,
        suffix_and_bwt_width,
    )?;

    let two_runs = checked_build_mul(
        ReferenceResource::EstimatedRetainedFmBytes,
        canonical_runs,
        2,
    )?;
    let bases_plus_two_runs = checked_build_add(
        ReferenceResource::EstimatedRetainedFmBytes,
        canonical_bases,
        two_runs,
    )?;
    let four_rank_rows = checked_build_mul(
        ReferenceResource::EstimatedRetainedFmBytes,
        bases_plus_two_runs,
        4,
    )?;
    let rank_width = physical_to_logical(
        size_of::<[u64; 4]>(),
        ReferenceResource::EstimatedRetainedFmBytes,
    )?;
    let term_three = checked_build_mul(
        ReferenceResource::EstimatedRetainedFmBytes,
        four_rank_rows,
        rank_width,
    )?;
    let first_two = checked_build_add(
        ReferenceResource::EstimatedRetainedFmBytes,
        term_one,
        term_two,
    )?;
    checked_build_add(
        ReferenceResource::EstimatedRetainedFmBytes,
        first_two,
        term_three,
    )
}

fn build_lane(
    original: &[Base],
    strand: BisulfiteStrand,
    contig_ordinal: u64,
    run_start: u64,
    max_suffix_rows: u64,
    scratch: &mut Vec<SearchBase>,
    reserved: usize,
) -> Result<FmIndex, ReferenceBuildError> {
    scratch.clear();
    let reverse = matches!(
        strand_semantics(strand).orientation(),
        AlignmentOrientation::Reverse
    );
    let conversion = lane_reference_conversion(strand);
    if reverse {
        for base in original.iter().rev() {
            push_lane_base(scratch, reserved, base.complement(), conversion)?;
        }
    } else {
        for base in original {
            push_lane_base(scratch, reserved, *base, conversion)?;
        }
    }
    FmIndex::build_reference(scratch, FmBuildLimit::new(max_suffix_rows)).map_err(|source| {
        ReferenceBuildError::FmBuild {
            contig_ordinal,
            run_start,
            strand,
            source,
        }
    })
}

fn push_lane_base(
    scratch: &mut Vec<SearchBase>,
    reserved: usize,
    base: Base,
    conversion: ThreeLetterConversion,
) -> Result<(), ReferenceBuildError> {
    if scratch.len() >= reserved {
        return Err(ReferenceBuildError::InternalInvariant {
            expected: u64::try_from(reserved).unwrap_or(u64::MAX),
            observed: u64::try_from(scratch.len()).unwrap_or(u64::MAX),
        });
    }
    let converted = convert_base(base, conversion);
    let Some(search) = SearchBase::from_base(converted) else {
        return Err(ReferenceBuildError::InternalInvariant {
            expected: 1,
            observed: 0,
        });
    };
    scratch.push(search);
    Ok(())
}

const fn lane_reference_conversion(strand: BisulfiteStrand) -> ThreeLetterConversion {
    match strand {
        BisulfiteStrand::OT | BisulfiteStrand::OB => ThreeLetterConversion::CToT,
        BisulfiteStrand::CTOT | BisulfiteStrand::CTOB => ThreeLetterConversion::GToA,
    }
}

fn raw_view_conversion(strand: BisulfiteStrand) -> ThreeLetterConversion {
    let semantics = strand_semantics(strand);
    match semantics.orientation() {
        AlignmentOrientation::Forward => semantics.search_conversion(),
        AlignmentOrientation::Reverse => dual_conversion(semantics.search_conversion()),
    }
}

const fn dual_conversion(conversion: ThreeLetterConversion) -> ThreeLetterConversion {
    match conversion {
        ThreeLetterConversion::CToT => ThreeLetterConversion::GToA,
        ThreeLetterConversion::GToA => ThreeLetterConversion::CToT,
    }
}

const fn convert_base(base: Base, conversion: ThreeLetterConversion) -> Base {
    match (base, conversion) {
        (Base::C, ThreeLetterConversion::CToT) => Base::T,
        (Base::G, ThreeLetterConversion::GToA) => Base::A,
        _ => base,
    }
}

fn recover_interval(
    run: RunCoordinates,
    offset: u64,
    pattern_len: u64,
    strand: BisulfiteStrand,
    contigs: &[ContigInput],
) -> Result<ReferenceInterval, ReferenceLocateError> {
    let run_len = run.len();
    if offset == run_len {
        return Err(ReferenceLocateError::Invariant {
            invariant: ReferenceLocateInvariant::TerminalSuffix,
            expected: run_len.saturating_sub(1),
            observed: offset,
        });
    }
    let Some(offset_end) = offset.checked_add(pattern_len) else {
        return Err(ReferenceLocateError::CoordinateArithmetic {
            contig_ordinal: run.contig_ordinal,
            run_start: run.start,
            offset,
            pattern_len,
        });
    };
    if offset_end > run_len {
        return Err(ReferenceLocateError::Invariant {
            invariant: ReferenceLocateInvariant::RunBounds,
            expected: run_len,
            observed: offset_end,
        });
    }

    let (start, end) = match strand_semantics(strand).orientation() {
        AlignmentOrientation::Forward => {
            let Some(start) = run.start.checked_add(offset) else {
                return Err(ReferenceLocateError::CoordinateArithmetic {
                    contig_ordinal: run.contig_ordinal,
                    run_start: run.start,
                    offset,
                    pattern_len,
                });
            };
            let Some(end) = start.checked_add(pattern_len) else {
                return Err(ReferenceLocateError::CoordinateArithmetic {
                    contig_ordinal: run.contig_ordinal,
                    run_start: run.start,
                    offset,
                    pattern_len,
                });
            };
            (start, end)
        }
        AlignmentOrientation::Reverse => {
            let Some(start) = run.end.checked_sub(offset_end) else {
                return Err(ReferenceLocateError::CoordinateArithmetic {
                    contig_ordinal: run.contig_ordinal,
                    run_start: run.start,
                    offset,
                    pattern_len,
                });
            };
            let Some(end) = run.end.checked_sub(offset) else {
                return Err(ReferenceLocateError::CoordinateArithmetic {
                    contig_ordinal: run.contig_ordinal,
                    run_start: run.start,
                    offset,
                    pattern_len,
                });
            };
            (start, end)
        }
    };

    let storage =
        usize::try_from(run.contig_ordinal).map_err(|_| ReferenceLocateError::Invariant {
            invariant: ReferenceLocateInvariant::MissingRun,
            expected: u64::try_from(contigs.len()).unwrap_or(u64::MAX),
            observed: run.contig_ordinal,
        })?;
    let Some(contig) = contigs.get(storage) else {
        return Err(ReferenceLocateError::Invariant {
            invariant: ReferenceLocateInvariant::MissingRun,
            expected: u64::try_from(contigs.len()).unwrap_or(u64::MAX),
            observed: run.contig_ordinal,
        });
    };
    ReferenceInterval::new(start, end, ReferenceLength::new(contig.sequence.len())).map_err(
        |source| ReferenceLocateError::CoordinateRecovery {
            contig_ordinal: run.contig_ordinal,
            run_start: run.start,
            strand,
            source,
        },
    )
}

struct CanonicalRuns<'a> {
    bases: &'a [Base],
    next: usize,
}

impl<'a> CanonicalRuns<'a> {
    const fn new(bases: &'a [Base]) -> Self {
        Self { bases, next: 0 }
    }
}

impl Iterator for CanonicalRuns<'_> {
    type Item = (usize, usize);

    fn next(&mut self) -> Option<Self::Item> {
        while self.next < self.bases.len() && self.bases[self.next] == Base::N {
            self.next += 1;
        }
        if self.next == self.bases.len() {
            return None;
        }
        let start = self.next;
        while self.next < self.bases.len() && self.bases[self.next] != Base::N {
            self.next += 1;
        }
        Some((start, self.next))
    }
}

fn ensure_build_capacity(
    expected: u64,
    reserved: usize,
    materialized: usize,
) -> Result<(), ReferenceBuildError> {
    if materialized >= reserved {
        return Err(ReferenceBuildError::InternalInvariant {
            expected,
            observed: physical_to_logical(materialized, ReferenceResource::CanonicalRuns)?,
        });
    }
    Ok(())
}

const fn ensure_query_capacity(
    reserved: u64,
    materialized: u64,
) -> Result<(), ReferenceQueryError> {
    if materialized >= reserved {
        return Err(ReferenceQueryError::CapacityInvariant {
            reserved,
            materialized,
        });
    }
    Ok(())
}

const fn ensure_query_count(
    counter: ReferenceQueryCounter,
    expected: u64,
    observed: u64,
) -> Result<(), ReferenceQueryError> {
    if observed != expected {
        return Err(ReferenceQueryError::InvariantMismatch {
            counter,
            expected,
            observed,
        });
    }
    Ok(())
}

const fn ensure_locate_capacity(
    reserved: u64,
    materialized: u64,
) -> Result<(), ReferenceLocateError> {
    if materialized >= reserved {
        return Err(ReferenceLocateError::Invariant {
            invariant: ReferenceLocateInvariant::FinalHitCapacity,
            expected: reserved,
            observed: materialized,
        });
    }
    Ok(())
}

const fn ensure_locate_equal(
    invariant: ReferenceLocateInvariant,
    expected: u64,
    observed: u64,
) -> Result<(), ReferenceLocateError> {
    if observed != expected {
        return Err(ReferenceLocateError::Invariant {
            invariant,
            expected,
            observed,
        });
    }
    Ok(())
}

fn query_storage_to_logical(value: usize) -> Result<u64, ReferenceQueryError> {
    u64::try_from(value)
        .map_err(|_| ReferenceQueryError::PatternLengthNotRepresentable { pattern_len: value })
}

pub(super) fn physical_to_logical(
    value: usize,
    resource: ReferenceResource,
) -> Result<u64, ReferenceBuildError> {
    u64::try_from(value).map_err(|_| ReferenceBuildError::CountNotRepresentable { resource, value })
}

pub(super) fn apply_limit(
    resource: ReferenceResource,
    requested: u64,
    maximum: u64,
) -> Result<(), ReferenceBuildError> {
    if requested > maximum {
        Err(ReferenceBuildError::LimitExceeded {
            resource,
            requested,
            maximum,
        })
    } else {
        Ok(())
    }
}

pub(super) fn checked_build_add(
    resource: ReferenceResource,
    lhs: u64,
    rhs: u64,
) -> Result<u64, ReferenceBuildError> {
    lhs.checked_add(rhs)
        .ok_or(ReferenceBuildError::ArithmeticOverflow {
            resource,
            operation: ReferenceArithmetic::Add,
            lhs,
            rhs,
        })
}

fn checked_build_mul(
    resource: ReferenceResource,
    lhs: u64,
    rhs: u64,
) -> Result<u64, ReferenceBuildError> {
    lhs.checked_mul(rhs)
        .ok_or(ReferenceBuildError::ArithmeticOverflow {
            resource,
            operation: ReferenceArithmetic::Multiply,
            lhs,
            rhs,
        })
}

fn checked_query_add(
    counter: ReferenceQueryCounter,
    accumulated: u64,
    next: u64,
) -> Result<u64, ReferenceQueryError> {
    accumulated
        .checked_add(next)
        .ok_or(ReferenceQueryError::CountOverflow {
            counter,
            accumulated,
            next,
        })
}

fn preflight_build_allocation<T>(
    allocation: ReferenceAllocation,
    elements: u64,
) -> Result<usize, ReferenceBuildError> {
    preflight_storage::<T>(elements).map_err(|(elements, element_size)| {
        ReferenceBuildError::AllocationSizeOverflow {
            allocation,
            elements,
            element_size,
        }
    })
}

fn preflight_query_allocation<T>(
    allocation: ReferenceAllocation,
    elements: u64,
) -> Result<usize, ReferenceQueryError> {
    preflight_storage::<T>(elements).map_err(|(elements, element_size)| {
        ReferenceQueryError::AllocationSizeOverflow {
            allocation,
            elements,
            element_size,
        }
    })
}

fn preflight_locate_allocation<T>(
    allocation: ReferenceAllocation,
    elements: u64,
) -> Result<usize, ReferenceLocateError> {
    preflight_storage::<T>(elements).map_err(|(elements, element_size)| {
        ReferenceLocateError::AllocationSizeOverflow {
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

fn contig_storage(ordinal: u64, contig_count: u64) -> Result<usize, ReferenceAccessError> {
    if ordinal >= contig_count {
        return Err(ReferenceAccessError::ContigOrdinalOutOfBounds {
            ordinal,
            contig_count,
        });
    }
    usize::try_from(ordinal).map_err(|_| ReferenceAccessError::ContigOrdinalOutOfBounds {
        ordinal,
        contig_count,
    })
}

const fn strand_rank(strand: BisulfiteStrand) -> u8 {
    match strand {
        BisulfiteStrand::OT => 0,
        BisulfiteStrand::OB => 1,
        BisulfiteStrand::CTOT => 2,
        BisulfiteStrand::CTOB => 3,
    }
}

#[cfg(test)]
#[path = "../../../tests/whitebox/reference.rs"]
mod tests;
