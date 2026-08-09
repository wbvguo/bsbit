//! Read-only adapter for the current combined-directional FM index.
//!
//! This module is intentionally feature-gated and format-specific. It reads the
//! current 16-mer table, Occ64/Occ65536 bit-plane rank, and SA16 sampled-row
//! representation without exposing that layout through bsbit's opaque index
//! command contract.

use core::fmt;
use core::mem::size_of;
use core::ptr::NonNull;
use std::cell::RefCell;
use std::ffi::OsString;
use std::fs::File;
use std::io::{self, Read};
use std::os::fd::AsRawFd;
#[cfg(target_os = "linux")]
use std::path::{Path, PathBuf};

use crate::reference::{
    CombinedIndexBackendError, PrivateCombinedIndex, PrivateCombinedLocateMetrics,
    PrivateCombinedReference, PrivateCombinedReferenceError, ReferenceIndex,
};
use crate::storage::fm::{FmInterval, ProjectedBase, SearchBase};
use crate::storage::reference_catalog::{
    ReferenceCatalogError, ReferenceCatalogSummary, load_reference_catalog,
};
use bsbit_core::reference::ReferenceSemanticDigest;

pub(crate) const BWT_WORDS_PER_128_ROWS: u64 = 5;
pub(crate) const SA_FLAG_WORDS_PER_256_ROWS: u64 = 5;

/// Returns all three LF boundaries from one validated packed-rank boundary.
///
/// Storage validation remains with the combined image reader; the builder
/// shares this arithmetic so encoding and decoding cannot silently diverge.
#[inline]
#[cfg(feature = "index-construction")]
pub(crate) fn lf_all_boundaries(
    boundary: u64,
    suffix_count: u64,
    sentinel_row: u64,
    first_occurrence: [u64; 4],
    mut bwt_word: impl FnMut(u64) -> u64,
    mut high_occ: impl FnMut(u64) -> u64,
) -> Option<[u64; 3]> {
    if boundary > suffix_count {
        return None;
    }
    let line = boundary - u64::from(boundary > sentinel_row);
    let high_word = (line >> 7).checked_mul(BWT_WORDS_PER_128_ROWS)?;
    let low_block = (line & 127) >> 6;
    let plane_start = high_word.checked_add(1 + (low_block << 1))?;
    let high_occ_block = (line >> 16).checked_mul(2)?;
    let counter_word = bwt_word(high_word);
    let first_plane = bwt_word(plane_start);
    let second_plane = bwt_word(plane_start + 1);
    let first_absolute = high_occ(high_occ_block);
    let second_absolute = high_occ(high_occ_block + 1);
    let counter_shift = low_block << 5;
    let packed = counter_word >> (32 - counter_shift);
    let nonzero = ((packed >> 16) & 0xffff) + (packed & 0xffff);
    let at_block = [
        ((line >> 6) << 6).checked_sub(
            first_absolute
                .checked_add(second_absolute)?
                .checked_add(nonzero)?,
        )?,
        first_absolute + ((counter_word >> (48 - counter_shift)) & 0xffff),
        second_absolute + ((counter_word >> (32 - counter_shift)) & 0xffff),
    ];
    let need = u32::try_from(line & 63).expect("six bits fit u32");
    let within = if need == 0 {
        [0_u64; 3]
    } else {
        let shift = 64 - need;
        [
            u64::from(((!(first_plane | second_plane)) >> shift).count_ones()),
            u64::from((first_plane >> shift).count_ones()),
            u64::from((second_plane >> shift).count_ones()),
        ]
    };
    Some([
        first_occurrence[0]
            .checked_add(at_block[0])?
            .checked_add(within[0])?,
        first_occurrence[1]
            .checked_add(at_block[1])?
            .checked_add(within[1])?,
        first_occurrence[2]
            .checked_add(at_block[2])?
            .checked_add(within[2])?,
    ])
}

pub(crate) const META_BYTES: usize = 120;
pub(crate) const META_BYTES_U32: u32 = 120;
pub(crate) const META_EXTENSION_MAGIC: &[u8; 8] = b"BSBICMB1";
pub(crate) const META_EXTENSION_MAJOR: u16 = 1;
pub(crate) const META_EXTENSION_MINOR: u16 = 0;
pub(crate) const META_EXTENSION_OFFSET: usize = 68;
pub(crate) const META_DIGEST_OFFSET: usize = 84;
const LOOKUP_BASES: usize = 16;
const LOOKUP_ENTRIES: u64 = 43_046_722;
const SA_STRIDE: u64 = 16;
const SA_STRIDE_U32: u32 = 16;
const SA_VALUE_MASK: u32 = 0x3fff_ffff;
const MAX_WAVEFRONT_LANES: usize = 64;
const MAX_WAVEFRONT_BOUNDARIES: usize = MAX_WAVEFRONT_LANES * 2;

#[derive(Clone, Copy, Debug)]
struct SameLowBlockRankPlan {
    lines: [u64; 2],
    high_word: u64,
    digit: u8,
}

#[derive(Clone, Copy, Debug)]
struct BoundaryRankPlan {
    line: u64,
    high_word: u64,
    digit: u8,
}

#[derive(Clone, Copy, Debug)]
struct ProjectedSuffixState {
    interval: FmInterval,
    matched_bases: usize,
    remaining_prefix_bases: usize,
    finished: bool,
}

impl ProjectedSuffixState {
    const fn new(
        interval: FmInterval,
        remaining_prefix_bases: usize,
        stop_interval_length: u64,
    ) -> Self {
        Self {
            interval,
            matched_bases: LOOKUP_BASES,
            remaining_prefix_bases,
            finished: interval.len() <= stop_interval_length || remaining_prefix_bases == 0,
        }
    }

    fn accept(&mut self, extended: FmInterval, stop_interval_length: u64) {
        if extended.is_empty() {
            self.finished = true;
            return;
        }
        self.interval = extended;
        self.matched_bases += 1;
        self.remaining_prefix_bases -= 1;
        self.finished = extended.len() <= stop_interval_length || self.remaining_prefix_bases == 0;
    }
}

struct BackwardExtendWavefrontWorkspace {
    plans: [Option<SameLowBlockRankPlan>; MAX_WAVEFRONT_LANES],
    boundary_plans: [Option<BoundaryRankPlan>; MAX_WAVEFRONT_BOUNDARIES],
    counters: [u64; MAX_WAVEFRONT_LANES],
    first_planes: [u64; MAX_WAVEFRONT_LANES],
    second_planes: [u64; MAX_WAVEFRONT_LANES],
    first_high_occ: [u64; MAX_WAVEFRONT_LANES],
    second_high_occ: [u64; MAX_WAVEFRONT_LANES],
    boundary_counters: [u64; MAX_WAVEFRONT_BOUNDARIES],
    boundary_first_planes: [u64; MAX_WAVEFRONT_BOUNDARIES],
    boundary_second_planes: [u64; MAX_WAVEFRONT_BOUNDARIES],
    boundary_first_high_occ: [u64; MAX_WAVEFRONT_BOUNDARIES],
    boundary_second_high_occ: [u64; MAX_WAVEFRONT_BOUNDARIES],
}

impl Default for BackwardExtendWavefrontWorkspace {
    fn default() -> Self {
        Self {
            plans: [None; MAX_WAVEFRONT_LANES],
            boundary_plans: [None; MAX_WAVEFRONT_BOUNDARIES],
            counters: [0; MAX_WAVEFRONT_LANES],
            first_planes: [0; MAX_WAVEFRONT_LANES],
            second_planes: [0; MAX_WAVEFRONT_LANES],
            first_high_occ: [0; MAX_WAVEFRONT_LANES],
            second_high_occ: [0; MAX_WAVEFRONT_LANES],
            boundary_counters: [0; MAX_WAVEFRONT_BOUNDARIES],
            boundary_first_planes: [0; MAX_WAVEFRONT_BOUNDARIES],
            boundary_second_planes: [0; MAX_WAVEFRONT_BOUNDARIES],
            boundary_first_high_occ: [0; MAX_WAVEFRONT_BOUNDARIES],
            boundary_second_high_occ: [0; MAX_WAVEFRONT_BOUNDARIES],
        }
    }
}

struct BackwardExtendRoundWorkspace {
    intervals: [Option<FmInterval>; MAX_WAVEFRONT_LANES],
    digits: [u8; MAX_WAVEFRONT_LANES],
    extended: [Option<FmInterval>; MAX_WAVEFRONT_LANES],
    rank: BackwardExtendWavefrontWorkspace,
}

impl Default for BackwardExtendRoundWorkspace {
    fn default() -> Self {
        Self {
            intervals: [None; MAX_WAVEFRONT_LANES],
            digits: [0; MAX_WAVEFRONT_LANES],
            extended: [None; MAX_WAVEFRONT_LANES],
            rank: BackwardExtendWavefrontWorkspace::default(),
        }
    }
}

thread_local! {
    /// Batched prefix extension performs adjacent rank rounds on the same
    /// mapping thread. Keeping the fixed wavefront buffers here avoids
    /// repeatedly zeroing several KiB of scratch state while retaining one
    /// independent workspace per mapping worker.
    static BACKWARD_EXTEND_ROUND_WORKSPACE: RefCell<BackwardExtendRoundWorkspace> =
        RefCell::new(BackwardExtendRoundWorkspace::default());
}

#[derive(Clone, Copy, Debug)]
struct SampleOrdinalPlan {
    block_word: u64,
    within: u64,
    flag_word: u64,
    bit: u64,
}

fn minimum_bwt_words_for_suffix_count(suffix_count: u64) -> Option<u64> {
    let text_rows = suffix_count.checked_sub(1)?;
    let last_text_line = text_rows.checked_sub(1)?;
    let last_text_low_block = (last_text_line & 127) >> 6;
    let minimum_text_words = (last_text_line >> 7)
        .checked_mul(BWT_WORDS_PER_128_ROWS)?
        .checked_add(3 + (last_text_low_block << 1))?;
    let boundary_low_block = (text_rows & 127) >> 6;
    let minimum_boundary_words = (text_rows >> 7)
        .checked_mul(BWT_WORDS_PER_128_ROWS)?
        .checked_add(if text_rows.is_multiple_of(64) {
            1
        } else {
            3 + (boundary_low_block << 1)
        })?;
    Some(minimum_text_words.max(minimum_boundary_words))
}

fn minimum_sa_flag_entries_for_suffix_count(suffix_count: u64) -> Option<u64> {
    let last_row = suffix_count.checked_sub(1)?;
    (last_row >> 8)
        .checked_mul(SA_FLAG_WORDS_PER_256_ROWS)?
        .checked_add(2 + ((last_row & 255) >> 6))
}

fn minimum_high_occ_entries_for_suffix_count(suffix_count: u64) -> Option<u64> {
    let text_rows = suffix_count.checked_sub(1)?;
    (text_rows >> 16).checked_add(1)?.checked_mul(2)
}

/// A validated format or I/O failure while opening the current combined index.
#[derive(Debug)]
pub enum CombinedIndexError {
    /// A file could not be opened, mapped, or read.
    Io(io::Error),
    /// One combined-index header, dimension, or offset violated the format.
    Structure(&'static str),
    /// The index was built from a different semantic reference catalog.
    ReferenceDigestMismatch {
        /// Digest required by the reference owner.
        expected: ReferenceSemanticDigest,
        /// Digest embedded in the combined-index metadata.
        observed: ReferenceSemanticDigest,
    },
}

impl fmt::Display for CombinedIndexError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(error) => write!(formatter, "combined index I/O: {error}"),
            Self::Structure(message) => {
                write!(formatter, "combined index structure: {message}")
            }
            Self::ReferenceDigestMismatch { expected, observed } => write!(
                formatter,
                "combined index reference digest differs: expected {expected}, observed {observed}"
            ),
        }
    }
}

impl std::error::Error for CombinedIndexError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Io(error) => Some(error),
            Self::Structure(_) | Self::ReferenceDigestMismatch { .. } => None,
        }
    }
}

#[derive(Debug)]
enum CombinedReferenceLoadFailure {
    Catalog(ReferenceCatalogError),
    Index(CombinedIndexError),
    Reference(PrivateCombinedReferenceError),
}

/// Opaque failure while binding a validated reference catalog to its search image.
#[derive(Debug)]
#[doc(hidden)]
pub struct CombinedReferenceLoadError {
    inner: CombinedReferenceLoadFailure,
}

impl fmt::Display for CombinedReferenceLoadError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match &self.inner {
            CombinedReferenceLoadFailure::Catalog(source) => {
                write!(formatter, "reference catalog validation failed: {source}")
            }
            CombinedReferenceLoadFailure::Index(source) => {
                write!(formatter, "combined-index validation failed: {source}")
            }
            CombinedReferenceLoadFailure::Reference(source) => {
                write!(formatter, "combined reference assembly failed: {source}")
            }
        }
    }
}

impl std::error::Error for CombinedReferenceLoadError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match &self.inner {
            CombinedReferenceLoadFailure::Catalog(source) => Some(source),
            CombinedReferenceLoadFailure::Index(source) => Some(source),
            CombinedReferenceLoadFailure::Reference(source) => Some(source),
        }
    }
}

/// A reference catalog bound to the current combined-index image used by alignment.
#[derive(Debug)]
#[doc(hidden)]
pub struct LoadedCombinedReference {
    index: ReferenceIndex,
    summary: ReferenceCatalogSummary,
    mapped_index_bytes: u64,
}

impl LoadedCombinedReference {
    /// Consumes the loaded bundle and returns the reference owner.
    #[must_use]
    pub fn into_index(self) -> ReferenceIndex {
        self.index
    }

    /// Returns the verified catalog summary.
    #[must_use]
    pub const fn summary(&self) -> ReferenceCatalogSummary {
        self.summary
    }

    /// Returns the logical bytes mapped by the combined-index components.
    #[must_use]
    pub const fn mapped_index_bytes(&self) -> u64 {
        self.mapped_index_bytes
    }
}

/// Loads and binds the two pieces created by one `bsbit index` run.
///
/// # Errors
///
/// Returns the first catalog, combined-index, semantic-identity, or reference
/// assembly failure.
#[doc(hidden)]
pub fn load_combined_reference_catalog(
    catalog_path: &Path,
    expected_reference_digest: Option<ReferenceSemanticDigest>,
    combined_index_prefix: &Path,
    threads: usize,
) -> Result<LoadedCombinedReference, CombinedReferenceLoadError> {
    let catalog = load_reference_catalog(catalog_path, expected_reference_digest, threads)
        .map_err(|source| CombinedReferenceLoadError {
            inner: CombinedReferenceLoadFailure::Catalog(source),
        })?;
    let summary = catalog.summary();
    let combined = CombinedIndex::open(combined_index_prefix).map_err(|source| {
        CombinedReferenceLoadError {
            inner: CombinedReferenceLoadFailure::Index(source),
        }
    })?;
    combined
        .verify_reference_semantic_digest(summary.semantic_digest())
        .map_err(|source| CombinedReferenceLoadError {
            inner: CombinedReferenceLoadFailure::Index(source),
        })?;
    let mapped_index_bytes = combined.mapped_bytes();
    let index = ReferenceIndex::from_private_combined(
        catalog.into_contigs(),
        PrivateCombinedReference::new(Box::new(combined)),
    )
    .map_err(|source| CombinedReferenceLoadError {
        inner: CombinedReferenceLoadFailure::Reference(source),
    })?;
    Ok(LoadedCombinedReference {
        index,
        summary,
        mapped_index_bytes,
    })
}

impl From<io::Error> for CombinedIndexError {
    fn from(error: io::Error) -> Self {
        Self::Io(error)
    }
}

#[derive(Debug)]
pub(crate) struct ReadOnlyMapping {
    pointer: NonNull<u8>,
    length: usize,
}

// SAFETY: the mapping is immutable for its complete lifetime and owns no Rust
// references. `munmap` runs only after the last shared owner is dropped.
unsafe impl Send for ReadOnlyMapping {}
// SAFETY: the mapping is immutable; all reads are bounded before dereference.
unsafe impl Sync for ReadOnlyMapping {}

impl ReadOnlyMapping {
    pub(crate) fn map(file: &File) -> Result<Self, CombinedIndexError> {
        let length = usize::try_from(file.metadata()?.len())
            .map_err(|_| CombinedIndexError::Structure("file length exceeds usize"))?;
        if length == 0 {
            return Err(CombinedIndexError::Structure("mapped file is empty"));
        }
        // SAFETY: the descriptor is live, `length` came from that descriptor,
        // and no writable mapping is requested. MAP_FAILED is checked below.
        let mapped = unsafe {
            libc::mmap(
                core::ptr::null_mut(),
                length,
                libc::PROT_READ,
                libc::MAP_PRIVATE,
                file.as_raw_fd(),
                0,
            )
        };
        if mapped == libc::MAP_FAILED {
            return Err(io::Error::last_os_error().into());
        }
        #[cfg(target_os = "linux")]
        // Advisory only.  Modern kernels may back aligned portions of a
        // read-only file mapping with larger folios; unsupported filesystems
        // retain ordinary demand paging without changing index semantics.
        unsafe {
            let _ = libc::madvise(mapped, length, libc::MADV_HUGEPAGE);
        }
        let pointer = NonNull::new(mapped.cast::<u8>()).ok_or(CombinedIndexError::Structure(
            "mmap returned a null pointer",
        ))?;
        Ok(Self { pointer, length })
    }

    #[inline]
    fn read_u8(&self, offset: usize) -> u8 {
        debug_assert!(offset < self.length);
        // SAFETY: every caller uses a validated component range and this
        // one-byte read is therefore inside the live read-only mapping.
        unsafe { self.pointer.as_ptr().add(offset).read() }
    }

    #[inline]
    fn read_u32(&self, offset: usize) -> u32 {
        debug_assert!(offset + size_of::<u32>() <= self.length);
        // SAFETY: the four-byte range was validated at open time. Index files
        // are little-endian; read_unaligned avoids imposing pointer alignment.
        u32::from_le(unsafe {
            self.pointer
                .as_ptr()
                .add(offset)
                .cast::<u32>()
                .read_unaligned()
        })
    }

    #[inline]
    pub(crate) fn read_u64(&self, offset: usize) -> u64 {
        debug_assert!(offset + size_of::<u64>() <= self.length);
        // SAFETY: the eight-byte range was validated at open time.
        u64::from_le(unsafe {
            self.pointer
                .as_ptr()
                .add(offset)
                .cast::<u64>()
                .read_unaligned()
        })
    }
}

impl Drop for ReadOnlyMapping {
    fn drop(&mut self) {
        // SAFETY: this is the exact pointer and nonzero length returned by the
        // successful mmap call, released exactly once by Drop.
        let _ = unsafe { libc::munmap(self.pointer.as_ptr().cast(), self.length) };
    }
}

/// Physical work used by sparse-SA location in the combined-index layout.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct CombinedLocateMetrics {
    /// Number of emitted suffix rows.
    pub located_rows: u64,
    /// Canonical row-wise LF depth represented by the located rows.
    pub lf_steps: u64,
    /// Logical rank boundaries evaluated by direct or shared traversal.
    pub rank_operations: u64,
    /// Number of direct rows or shared interval-tree nodes processed.
    pub interval_nodes: u64,
}

/// Read-only current combined-directional FM index.
#[derive(Debug)]
pub struct CombinedIndex {
    bwt: ReadOnlyMapping,
    sa: ReadOnlyMapping,
    occ: ReadOnlyMapping,
    suffix_count: u64,
    sentinel_row: u64,
    first_occurrence: [u64; 4],
    bwt_words: u64,
    lookup_entries: u64,
    lookup_high_offset: usize,
    lookup_low_offset: usize,
    sparse_sa_entries: u64,
    sa_values_offset: usize,
    sa_flag_entries: u64,
    sa_flags_offset: usize,
    high_occ_entries: u64,
    high_occ_offset: usize,
    reference_semantic_digest: ReferenceSemanticDigest,
}

impl CombinedIndex {
    /// Opens `<prefix>`, `<prefix>.bwt`, `<prefix>.sa`, and `<prefix>.occ`.
    ///
    /// # Errors
    ///
    /// Rejects missing files and every unsupported or inconsistent combined-index
    /// dimension before publishing a mapped index.
    #[allow(clippy::too_many_lines)]
    pub fn open(prefix: &Path) -> Result<Self, CombinedIndexError> {
        let mut meta_file = File::open(prefix)?;
        let meta_length = usize::try_from(meta_file.metadata()?.len())
            .map_err(|_| CombinedIndexError::Structure("metadata file length exceeds usize"))?;
        if meta_length != META_BYTES {
            return Err(CombinedIndexError::Structure(
                "metadata file is not the current 120-byte bound format",
            ));
        }
        let mut meta = [0_u8; META_BYTES];
        meta_file.read_exact(&mut meta)?;
        let suffix_count = slice_u64(&meta, 0);
        let sentinel_row = slice_u64(&meta, 8);
        let first_occurrence = core::array::from_fn(|index| slice_u64(&meta, 16 + index * 8));
        let terminal = slice_u64(&meta, 48);
        let sa_stride = slice_u32(&meta, 56);
        let occ_stride = slice_u32(&meta, 60);
        let high_occ_stride = slice_u32(&meta, 64);
        if &meta[META_EXTENSION_OFFSET..META_EXTENSION_OFFSET + 8] != META_EXTENSION_MAGIC
            || slice_u16(&meta, 76) != META_EXTENSION_MAJOR
            || slice_u16(&meta, 78) != META_EXTENSION_MINOR
            || slice_u32(&meta, 80) != META_BYTES_U32
            || meta[116..120] != [0; 4]
        {
            return Err(CombinedIndexError::Structure(
                "metadata extension version, length, or reserved bytes are invalid",
            ));
        }
        let reference_semantic_digest = ReferenceSemanticDigest::from_bytes(metadata_digest(&meta));
        if suffix_count < 3
            || suffix_count.is_multiple_of(2)
            || sentinel_row >= suffix_count
            || terminal != suffix_count
            || first_occurrence.windows(2).any(|pair| pair[0] > pair[1])
            || first_occurrence.iter().any(|&value| value > terminal)
        {
            return Err(CombinedIndexError::Structure(
                "metadata suffix or cumulative-count domain is invalid",
            ));
        }
        if sa_stride != SA_STRIDE_U32 || occ_stride != 64 || high_occ_stride != 128 {
            return Err(CombinedIndexError::Structure(
                "only the current SA16/Occ64/Occ128 layout is supported",
            ));
        }

        let bwt = map_suffix(prefix, ".bwt")?;
        let bwt_words = bwt.read_u64(0);
        let lookup_header_offset = checked_component_end(8, bwt_words, 8, bwt.length)?;
        let lookup_entries = bwt.read_u64(lookup_header_offset);
        if lookup_entries != LOOKUP_ENTRIES {
            return Err(CombinedIndexError::Structure(
                "combined lookup is not the dense three-letter 16-mer table",
            ));
        }
        let lookup_high_offset = lookup_header_offset
            .checked_add(8)
            .ok_or(CombinedIndexError::Structure("BWT offset overflow"))?;
        let lookup_low_offset =
            checked_component_end(lookup_high_offset, lookup_entries, 4, bwt.length)?;
        let bwt_end = checked_component_end(lookup_low_offset, lookup_entries, 1, bwt.length)?;
        if bwt_end != bwt.length {
            return Err(CombinedIndexError::Structure(
                "BWT file has trailing or missing lookup bytes",
            ));
        }
        let sa = map_suffix(prefix, ".sa")?;
        let sparse_sa_entries = sa.read_u64(0);
        let sa_values_offset = 8;
        let sa_flag_header_offset =
            checked_component_end(sa_values_offset, sparse_sa_entries, 4, sa.length)?;
        let sa_flag_entries = sa.read_u64(sa_flag_header_offset);
        let sa_flags_offset = sa_flag_header_offset
            .checked_add(8)
            .ok_or(CombinedIndexError::Structure("SA offset overflow"))?;
        let flags_end = checked_component_end(sa_flags_offset, sa_flag_entries, 8, sa.length)?;
        if flags_end != sa.length || sparse_sa_entries == 0 || sa_flag_entries == 0 {
            return Err(CombinedIndexError::Structure(
                "SA file dimensions are invalid",
            ));
        }
        let occ = map_suffix(prefix, ".occ")?;
        let high_occ_entries = occ.read_u64(0);
        let high_occ_offset = 8;
        let occ_end = checked_component_end(high_occ_offset, high_occ_entries, 8, occ.length)?;
        if occ_end != occ.length || high_occ_entries < 2 || high_occ_entries % 2 != 0 {
            return Err(CombinedIndexError::Structure(
                "high-occurrence file dimensions are invalid",
            ));
        }
        let high_occ_entries_usize = usize::try_from(high_occ_entries).map_err(|_| {
            CombinedIndexError::Structure("high-occurrence entry count exceeds usize")
        })?;
        for ordinal in 0..high_occ_entries_usize {
            if occ.read_u64(high_occ_offset + ordinal * 8) > suffix_count {
                return Err(CombinedIndexError::Structure(
                    "high-occurrence value exceeds suffix domain",
                ));
            }
        }

        let index = Self {
            bwt,
            sa,
            occ,
            suffix_count,
            sentinel_row,
            first_occurrence,
            bwt_words,
            lookup_entries,
            lookup_high_offset,
            lookup_low_offset,
            sparse_sa_entries,
            sa_values_offset,
            sa_flag_entries,
            sa_flags_offset,
            high_occ_entries,
            high_occ_offset,
            reference_semantic_digest,
        };
        index.validate_runtime_dimensions()?;
        Ok(index)
    }

    /// Returns the combined suffix-row count, including the terminal suffix.
    #[must_use]
    pub const fn suffix_count(&self) -> u64 {
        self.suffix_count
    }

    /// Returns the forward reference length represented by each physical half.
    #[must_use]
    pub const fn reference_length(&self) -> u64 {
        (self.suffix_count - 1) / 2
    }

    /// Returns the semantic reference digest embedded by the index builder.
    #[must_use]
    pub const fn reference_semantic_digest(&self) -> ReferenceSemanticDigest {
        self.reference_semantic_digest
    }

    /// Requires this image to be bound to the supplied reference catalog.
    ///
    /// # Errors
    ///
    /// Rejects every digest mismatch.
    pub fn verify_reference_semantic_digest(
        &self,
        expected: ReferenceSemanticDigest,
    ) -> Result<(), CombinedIndexError> {
        let observed = self.reference_semantic_digest;
        if observed != expected {
            return Err(CombinedIndexError::ReferenceDigestMismatch { expected, observed });
        }
        Ok(())
    }

    /// Returns bytes mapped from the combined-index BWT, SA, Occ, and optional packed
    /// reference files.
    #[must_use]
    pub fn mapped_bytes(&self) -> u64 {
        u64::try_from(self.bwt.length)
            .unwrap_or(u64::MAX)
            .saturating_add(u64::try_from(self.sa.length).unwrap_or(u64::MAX))
            .saturating_add(u64::try_from(self.occ.length).unwrap_or(u64::MAX))
    }

    /// Looks up exactly 16 projected bases in the current dense table.
    #[must_use]
    pub fn lookup16(&self, pattern: &[SearchBase]) -> Option<FmInterval> {
        self.lookup16_symbols(pattern)
    }

    #[inline]
    fn lookup16_symbols<T: ProjectedQuerySymbol>(&self, pattern: &[T]) -> Option<FmInterval> {
        if pattern.len() != LOOKUP_BASES {
            return None;
        }
        let mut key = 0_u64;
        for &base in pattern {
            key = key * 3 + u64::from(base.projected_digit()?);
        }
        self.lookup16_key(key)
    }

    #[inline]
    pub(crate) fn lookup16_key(&self, key: u64) -> Option<FmInterval> {
        if key >= LOOKUP_ENTRIES - 1 {
            return None;
        }
        let next = key + 1;
        debug_assert!(next < self.lookup_entries);
        let lower_high = self.lookup_high(key);
        let upper_high = self.lookup_high(next);
        let key = usize::try_from(key).expect("validated lookup key fits usize");
        let next = usize::try_from(next).expect("validated lookup key fits usize");
        let lower = (u64::from(lower_high & 0x0fff_ffff) << 8)
            | u64::from(self.bwt.read_u8(self.lookup_low_offset + key));
        let upper = ((u64::from(upper_high & 0x0fff_ffff) << 8)
            | u64::from(self.bwt.read_u8(self.lookup_low_offset + next)))
        .checked_sub(u64::from(upper_high >> 28))?;
        FmInterval::private_checked(lower, upper, self.suffix_count).ok()
    }

    /// Prepends one projected symbol to an interval in this combined domain.
    #[must_use]
    pub fn backward_extend(&self, interval: FmInterval, base: SearchBase) -> Option<FmInterval> {
        if interval.private_suffix_count() != self.suffix_count
            || interval.upper() > self.suffix_count
        {
            return None;
        }
        let digit = projected_digit(base)?;
        self.backward_extend_validated(interval, digit)
    }

    #[inline(never)]
    fn backward_extend_validated(&self, interval: FmInterval, digit: u8) -> Option<FmInterval> {
        debug_assert_eq!(interval.private_suffix_count(), self.suffix_count);
        debug_assert!(interval.upper() <= self.suffix_count);
        debug_assert!(digit <= 2);
        let [lower, upper] = self.lf_boundary_pair(interval.lower(), interval.upper(), digit)?;
        FmInterval::private_checked(lower, upper, self.suffix_count).ok()
    }

    fn backward_extend_interval_round(
        &self,
        intervals: &[FmInterval],
        digits: &[u8],
        output: &mut [FmInterval],
    ) -> Result<(), CombinedIndexBackendError> {
        if intervals.len() != digits.len() || intervals.len() != output.len() {
            return Err(CombinedIndexBackendError::Structure);
        }
        if intervals.is_empty() {
            return Ok(());
        }
        if intervals.len() > MAX_WAVEFRONT_LANES {
            for ((output, &interval), &digit) in output.iter_mut().zip(intervals).zip(digits) {
                *output = self
                    .backward_extend_validated(interval, digit)
                    .ok_or(CombinedIndexBackendError::Interval)?;
            }
            return Ok(());
        }
        BACKWARD_EXTEND_ROUND_WORKSPACE.with(|workspace| {
            let mut workspace = workspace.borrow_mut();
            let BackwardExtendRoundWorkspace {
                intervals: workspace_intervals,
                digits: workspace_digits,
                extended,
                rank,
            } = &mut *workspace;
            for (destination, &interval) in workspace_intervals.iter_mut().zip(intervals) {
                *destination = Some(interval);
            }
            workspace_digits[..digits.len()].copy_from_slice(digits);
            self.backward_extend_wavefront_validated_with_workspace(
                workspace_intervals,
                workspace_digits,
                intervals.len(),
                extended,
                rank,
            );
            for (output, interval) in output.iter_mut().zip(&extended[..intervals.len()]) {
                *output = interval.ok_or(CombinedIndexBackendError::Interval)?;
            }
            Ok(())
        })
    }

    // Dense lookup and every dependent rank round form one physical memory
    // transaction; splitting them would recreate the hot-path regression this
    // backend primitive exists to avoid.
    #[allow(clippy::too_many_lines)]
    fn resolve_projected_suffix_intervals_impl(
        &self,
        patterns: &[&[ProjectedBase]],
        minimum_suffix_bases: usize,
        stop_interval_length: u64,
        output: &mut [Option<(FmInterval, u64)>],
    ) -> Result<(), CombinedIndexBackendError> {
        if patterns.len() != output.len()
            || patterns.len() > MAX_WAVEFRONT_LANES
            || stop_interval_length == 0
        {
            return Err(CombinedIndexBackendError::Structure);
        }
        output.fill(None);
        let lane_count = patterns.len();
        let mut states = [None; MAX_WAVEFRONT_LANES];
        let mut lookup_keys = [None; MAX_WAVEFRONT_LANES];
        for lane in 0..lane_count {
            let pattern = patterns[lane];
            if minimum_suffix_bases < LOOKUP_BASES || pattern.len() < minimum_suffix_bases {
                continue;
            }
            let suffix = &pattern[pattern.len() - LOOKUP_BASES..];
            let mut key = 0_u64;
            for &base in suffix {
                key = key * 3 + u64::from(base.digit());
            }
            if key < LOOKUP_ENTRIES - 1 {
                lookup_keys[lane] = Some(key);
            }
        }

        // Issue each lookup component as a separate lane loop. The dense
        // table is much larger than cache, so independent addresses let the
        // CPU overlap memory latency without owning caller stopping policy.
        let mut lower_high = [0_u32; MAX_WAVEFRONT_LANES];
        let mut upper_high = [0_u32; MAX_WAVEFRONT_LANES];
        for lane in 0..lane_count {
            let Some(key) = lookup_keys[lane] else {
                continue;
            };
            lower_high[lane] = self.lookup_high(key);
            upper_high[lane] = self.lookup_high(key + 1);
        }
        let mut lower_low = [0_u8; MAX_WAVEFRONT_LANES];
        let mut upper_low = [0_u8; MAX_WAVEFRONT_LANES];
        for lane in 0..lane_count {
            let Some(key) = lookup_keys[lane] else {
                continue;
            };
            let key = usize::try_from(key).expect("validated lookup key fits usize");
            lower_low[lane] = self.bwt.read_u8(self.lookup_low_offset + key);
            upper_low[lane] = self.bwt.read_u8(self.lookup_low_offset + key + 1);
        }
        for lane in 0..lane_count {
            let Some(_) = lookup_keys[lane] else {
                continue;
            };
            let lower =
                (u64::from(lower_high[lane] & 0x0fff_ffff) << 8) | u64::from(lower_low[lane]);
            let Some(upper) = ((u64::from(upper_high[lane] & 0x0fff_ffff) << 8)
                | u64::from(upper_low[lane]))
            .checked_sub(u64::from(upper_high[lane] >> 28)) else {
                continue;
            };
            let Ok(interval) = FmInterval::private_checked(lower, upper, self.suffix_count) else {
                continue;
            };
            if !interval.is_empty() {
                states[lane] = Some(ProjectedSuffixState::new(
                    interval,
                    patterns[lane].len() - LOOKUP_BASES,
                    stop_interval_length,
                ));
            }
        }

        let mut active_lanes = [0_usize; MAX_WAVEFRONT_LANES];
        BACKWARD_EXTEND_ROUND_WORKSPACE.with(|workspace| {
            let mut workspace = workspace.borrow_mut();
            let BackwardExtendRoundWorkspace {
                intervals,
                digits,
                extended,
                rank,
            } = &mut *workspace;
            loop {
                let mut active_count = 0_usize;
                for lane in 0..lane_count {
                    let Some(state) = states[lane].filter(|state| !state.finished) else {
                        continue;
                    };
                    intervals[active_count] = Some(state.interval);
                    digits[active_count] = patterns[lane][state.remaining_prefix_bases - 1].digit();
                    active_lanes[active_count] = lane;
                    active_count += 1;
                }
                if active_count == 0 {
                    break;
                }
                self.backward_extend_wavefront_validated_with_workspace(
                    intervals,
                    digits,
                    active_count,
                    extended,
                    rank,
                );
                for active in 0..active_count {
                    let lane = active_lanes[active];
                    let interval = extended[active].ok_or(CombinedIndexBackendError::Interval)?;
                    states[lane]
                        .as_mut()
                        .expect("active projected-suffix state exists")
                        .accept(interval, stop_interval_length);
                }
            }
            Ok::<(), CombinedIndexBackendError>(())
        })?;

        for lane in 0..lane_count {
            output[lane] = states[lane]
                .map(|state| {
                    u64::try_from(state.matched_bases)
                        .map(|matched_bases| (state.interval, matched_bases))
                        .map_err(|_| CombinedIndexBackendError::Structure)
                })
                .transpose()?;
        }
        Ok(())
    }

    /// Returns the exact interval for a projected pattern of at least 16 bases.
    #[must_use]
    pub fn exact_search(&self, pattern: &[SearchBase]) -> Option<FmInterval> {
        if pattern.len() < LOOKUP_BASES {
            return None;
        }
        let split = pattern.len() - LOOKUP_BASES;
        let mut interval = self.lookup16(&pattern[split..])?;
        for &base in pattern[..split].iter().rev() {
            interval = self.backward_extend(interval, base)?;
            if interval.is_empty() {
                break;
            }
        }
        Some(interval)
    }

    /// Locates every row in an interval using direct SA completion.
    ///
    /// # Errors
    ///
    /// Rejects a foreign interval, invalid sampled-row rank, or an LF path
    /// longer than the declared sampling distance.
    pub fn visit_interval(
        &self,
        interval: FmInterval,
        visitor: &mut dyn FnMut(u64) -> bool,
    ) -> Result<CombinedLocateMetrics, CombinedIndexError> {
        if interval.private_suffix_count() != self.suffix_count
            || interval.upper() > self.suffix_count
        {
            return Err(CombinedIndexError::Structure(
                "locate interval belongs to another FM domain",
            ));
        }
        self.visit_interval_direct(interval, visitor)
    }

    #[inline(never)]
    fn visit_interval_direct(
        &self,
        interval: FmInterval,
        visitor: &mut dyn FnMut(u64) -> bool,
    ) -> Result<CombinedLocateMetrics, CombinedIndexError> {
        let mut metrics = CombinedLocateMetrics::default();
        let mut row = interval.lower();
        while row.saturating_add(1) < interval.upper() {
            let input = [row, row + 1];
            let located = self.locate_rows_two_lanes(input)?;
            for (position, steps) in located {
                metrics.located_rows += 1;
                metrics.lf_steps += steps;
                metrics.rank_operations += steps;
                if !visitor(position) {
                    metrics.interval_nodes = metrics.located_rows;
                    return Ok(metrics);
                }
            }
            row += 2;
        }
        if row < interval.upper() {
            let (position, steps) = self.locate_row(row)?;
            metrics.located_rows += 1;
            metrics.lf_steps += steps;
            metrics.rank_operations += steps;
            let _ = visitor(position);
        }
        metrics.interval_nodes = metrics.located_rows;
        Ok(metrics)
    }

    fn visit_interval_two_lanes_complete(
        &self,
        intervals: [FmInterval; 2],
        visitor: &mut dyn FnMut(usize, u64),
    ) -> Result<[CombinedLocateMetrics; 2], CombinedIndexError> {
        for interval in intervals {
            if interval.private_suffix_count() != self.suffix_count
                || interval.upper() > self.suffix_count
            {
                return Err(CombinedIndexError::Structure(
                    "locate interval belongs to another FM domain",
                ));
            }
        }
        let mut metrics = [CombinedLocateMetrics::default(); 2];
        let mut rows = [intervals[0].lower(), intervals[1].lower()];
        let upper = [intervals[0].upper(), intervals[1].upper()];
        while rows[0] < upper[0] || rows[1] < upper[1] {
            let lanes = if rows[0] < upper[0] && rows[1] < upper[1] {
                [0_usize, 1_usize]
            } else if rows[0].saturating_add(1) < upper[0] {
                [0_usize, 0_usize]
            } else if rows[1].saturating_add(1) < upper[1] {
                [1_usize, 1_usize]
            } else {
                let lane = usize::from(rows[0] >= upper[0]);
                let (position, steps) = self.locate_row(rows[lane])?;
                rows[lane] += 1;
                metrics[lane].located_rows += 1;
                metrics[lane].lf_steps += steps;
                metrics[lane].rank_operations += steps;
                visitor(lane, position);
                continue;
            };
            let input = if lanes[0] == lanes[1] {
                [rows[lanes[0]], rows[lanes[0]] + 1]
            } else {
                [rows[lanes[0]], rows[lanes[1]]]
            };
            rows[lanes[0]] += 1;
            rows[lanes[1]] += 1;
            let located = self.locate_rows_two_lanes(input)?;
            for ordinal in 0..2 {
                let lane = lanes[ordinal];
                let (position, steps) = located[ordinal];
                metrics[lane].located_rows += 1;
                metrics[lane].lf_steps += steps;
                metrics[lane].rank_operations += steps;
                visitor(lane, position);
            }
        }
        for metric in &mut metrics {
            metric.interval_nodes = metric.located_rows;
        }
        Ok(metrics)
    }

    fn validate_runtime_dimensions(&self) -> Result<(), CombinedIndexError> {
        let minimum_bwt_words = minimum_bwt_words_for_suffix_count(self.suffix_count).ok_or(
            CombinedIndexError::Structure("BWT word dimensions overflow"),
        )?;
        let minimum_sa_flags = minimum_sa_flag_entries_for_suffix_count(self.suffix_count)
            .ok_or(CombinedIndexError::Structure("SA-flag dimensions overflow"))?;
        let minimum_high_occ = minimum_high_occ_entries_for_suffix_count(self.suffix_count).ok_or(
            CombinedIndexError::Structure("high-occurrence dimensions overflow"),
        )?;
        if self.bwt_words < minimum_bwt_words
            || self.sa_flag_entries < minimum_sa_flags
            || self.high_occ_entries < minimum_high_occ
        {
            return Err(CombinedIndexError::Structure(
                "rank or sampled-SA arrays are shorter than their row domain",
            ));
        }
        if self.sample_rank(self.suffix_count)? != self.sparse_sa_entries {
            return Err(CombinedIndexError::Structure(
                "sampled-SA flags and sparse values disagree",
            ));
        }
        Ok(())
    }

    #[inline]
    fn lookup_high(&self, key: u64) -> u32 {
        debug_assert!(key < self.lookup_entries);
        let key = usize::try_from(key).expect("validated lookup key fits usize");
        self.bwt.read_u32(self.lookup_high_offset + key * 4)
    }

    #[inline]
    fn bwt_word(&self, ordinal: u64) -> u64 {
        debug_assert!(ordinal < self.bwt_words);
        let ordinal = usize::try_from(ordinal).expect("validated BWT ordinal fits usize");
        self.bwt.read_u64(8 + ordinal * 8)
    }

    #[inline]
    fn high_occ(&self, ordinal: u64) -> u64 {
        debug_assert!(ordinal < self.high_occ_entries);
        let ordinal = usize::try_from(ordinal).expect("validated Occ ordinal fits usize");
        self.occ.read_u64(self.high_occ_offset + ordinal * 8)
    }

    #[inline]
    fn sparse_sa(&self, ordinal: u64) -> Option<u64> {
        (ordinal < self.sparse_sa_entries).then(|| {
            let ordinal = usize::try_from(ordinal).expect("validated SA ordinal fits usize");
            u64::from(self.sa.read_u32(self.sa_values_offset + ordinal * 4))
        })
    }

    #[inline]
    fn sa_flag_word(&self, ordinal: u64) -> u64 {
        debug_assert!(ordinal < self.sa_flag_entries);
        let ordinal = usize::try_from(ordinal).expect("validated flag ordinal fits usize");
        self.sa.read_u64(self.sa_flags_offset + ordinal * 8)
    }

    #[inline]
    fn lf_boundary(&self, boundary: u64, digit: u8) -> Option<u64> {
        debug_assert!(boundary <= self.suffix_count);
        debug_assert!(digit <= 2);
        let mut line = boundary;
        if line > self.sentinel_row {
            line -= 1;
        }
        let high_block = line >> 7;
        let high_word = high_block.checked_mul(BWT_WORDS_PER_128_ROWS)?;
        let low_block = (line & 127) >> 6;
        let plane_start = high_word.checked_add(1 + (low_block << 1))?;
        let high_occ_block = (line >> 16).checked_mul(2)?;
        let counter_word = self.bwt_word(high_word);
        let at_block = if digit == 0 {
            let g = self.high_occ(high_occ_block);
            let t = self.high_occ(high_occ_block + 1);
            let packed = counter_word >> (32 - (low_block << 5));
            let non_a = ((packed >> 16) & 0xffff) + (packed & 0xffff);
            ((line >> 6) << 6).checked_sub(g.checked_add(t)?.checked_add(non_a)?)?
        } else {
            let plane = u64::from(digit >> 1);
            let absolute = self.high_occ(high_occ_block + plane);
            let shift = (48 - (plane << 4)).checked_sub(low_block << 5)?;
            absolute + ((counter_word >> shift) & 0xffff)
        };
        let need = u32::try_from(line & 63).expect("six bits fit u32");
        let within = if need == 0 {
            0
        } else if digit == 0 {
            let non_a = self.bwt_word(plane_start) | self.bwt_word(plane_start + 1);
            ((!non_a) >> (64 - need)).count_ones()
        } else {
            let plane = u64::from(digit >> 1);
            (self.bwt_word(plane_start + plane) >> (64 - need)).count_ones()
        };
        self.first_occurrence[usize::from(digit)]
            .checked_add(at_block)?
            .checked_add(u64::from(within))
    }

    #[inline]
    fn boundary_rank_plan(&self, boundary: u64, digit: u8) -> Option<BoundaryRankPlan> {
        debug_assert!(boundary <= self.suffix_count);
        debug_assert!(digit <= 2);
        let line = boundary - u64::from(boundary > self.sentinel_row);
        let high_word = (line >> 7).checked_mul(BWT_WORDS_PER_128_ROWS)?;
        let low_block = (line & 127) >> 6;
        high_word.checked_add(1 + (low_block << 1))?;
        (line >> 16).checked_mul(2)?;
        Some(BoundaryRankPlan {
            line,
            high_word,
            digit,
        })
    }

    #[inline]
    fn finish_boundary_rank_plan(
        &self,
        plan: BoundaryRankPlan,
        counter_word: u64,
        planes: [u64; 2],
        high_occ: [u64; 2],
    ) -> Option<u64> {
        let low_block = (plan.line & 127) >> 6;
        let at_block = if plan.digit == 0 {
            let packed = counter_word >> (32 - (low_block << 5));
            let non_a = ((packed >> 16) & 0xffff) + (packed & 0xffff);
            ((plan.line >> 6) << 6)
                .checked_sub(high_occ[0].checked_add(high_occ[1])?.checked_add(non_a)?)?
        } else {
            let plane = usize::from(plan.digit >> 1);
            let shift = (48 - (u64::try_from(plane).ok()? << 4)).checked_sub(low_block << 5)?;
            high_occ[plane] + ((counter_word >> shift) & 0xffff)
        };
        let need = u32::try_from(plan.line & 63).expect("six bits fit u32");
        let within = if need == 0 {
            0
        } else if plan.digit == 0 {
            u64::from(((!(planes[0] | planes[1])) >> (64 - need)).count_ones())
        } else {
            u64::from((planes[usize::from(plan.digit >> 1)] >> (64 - need)).count_ones())
        };
        self.first_occurrence[usize::from(plan.digit)]
            .checked_add(at_block)?
            .checked_add(within)
    }

    #[inline]
    fn same_low_block_rank_plan(
        &self,
        boundaries: [u64; 2],
        digit: u8,
    ) -> Option<SameLowBlockRankPlan> {
        debug_assert!(boundaries[0] <= boundaries[1]);
        debug_assert!(boundaries[1] <= self.suffix_count);
        debug_assert!(digit <= 2);
        let lines = boundaries.map(|boundary| boundary - u64::from(boundary > self.sentinel_row));
        if lines[0] >> 6 != lines[1] >> 6 {
            return None;
        }
        let high_word = (lines[0] >> 7).checked_mul(BWT_WORDS_PER_128_ROWS)?;
        let low_block = (lines[0] & 127) >> 6;
        high_word.checked_add(1 + (low_block << 1))?;
        (lines[0] >> 16).checked_mul(2)?;
        Some(SameLowBlockRankPlan {
            lines,
            high_word,
            digit,
        })
    }

    #[inline]
    fn finish_same_low_block_rank_plan(
        &self,
        plan: SameLowBlockRankPlan,
        counter_word: u64,
        planes: [u64; 2],
        high_occ: [u64; 2],
    ) -> Option<[u64; 2]> {
        let low_block = (plan.lines[0] & 127) >> 6;
        let at_block = if plan.digit == 0 {
            let packed = counter_word >> (32 - (low_block << 5));
            let non_a = ((packed >> 16) & 0xffff) + (packed & 0xffff);
            ((plan.lines[0] >> 6) << 6)
                .checked_sub(high_occ[0].checked_add(high_occ[1])?.checked_add(non_a)?)?
        } else {
            let plane = u64::from(plan.digit >> 1);
            let shift = (48 - (plane << 4)).checked_sub(low_block << 5)?;
            high_occ[usize::try_from(plane).expect("combined plane fits usize")]
                + ((counter_word >> shift) & 0xffff)
        };
        let plane_bits = if plan.digit == 0 {
            !(planes[0] | planes[1])
        } else {
            planes[usize::from(plan.digit >> 1)]
        };
        let first = self.first_occurrence[usize::from(plan.digit)].checked_add(at_block)?;
        let rank = |line| {
            let need = u32::try_from(line & 63).expect("six bits fit u32");
            let within = if need == 0 {
                0
            } else {
                u64::from((plane_bits >> (64 - need)).count_ones())
            };
            first.checked_add(within)
        };
        Some([rank(plan.lines[0])?, rank(plan.lines[1])?])
    }

    #[inline]
    #[allow(clippy::too_many_lines)]
    fn backward_extend_wavefront_validated_with_workspace(
        &self,
        intervals: &[Option<FmInterval>; MAX_WAVEFRONT_LANES],
        digits: &[u8; MAX_WAVEFRONT_LANES],
        lane_count: usize,
        output: &mut [Option<FmInterval>; MAX_WAVEFRONT_LANES],
        workspace: &mut BackwardExtendWavefrontWorkspace,
    ) {
        debug_assert!(lane_count <= MAX_WAVEFRONT_LANES);
        let BackwardExtendWavefrontWorkspace {
            plans,
            boundary_plans,
            counters,
            first_planes,
            second_planes,
            first_high_occ,
            second_high_occ,
            boundary_counters,
            boundary_first_planes,
            boundary_second_planes,
            boundary_first_high_occ,
            boundary_second_high_occ,
        } = workspace;
        for lane in 0..lane_count {
            output[lane] = None;
            boundary_plans[lane * 2] = None;
            boundary_plans[lane * 2 + 1] = None;
            let Some(interval) = intervals[lane] else {
                plans[lane] = None;
                continue;
            };
            debug_assert_eq!(interval.private_suffix_count(), self.suffix_count);
            debug_assert!(interval.upper() <= self.suffix_count);
            debug_assert!(digits[lane] <= 2);
            plans[lane] =
                self.same_low_block_rank_plan([interval.lower(), interval.upper()], digits[lane]);
            if plans[lane].is_none() {
                boundary_plans[lane * 2] = self.boundary_rank_plan(interval.lower(), digits[lane]);
                boundary_plans[lane * 2 + 1] =
                    self.boundary_rank_plan(interval.upper(), digits[lane]);
            }
        }

        // Keep address planning, independent memory issue, and arithmetic in
        // separate loops. This exposes a small exact wavefront to an ordinary
        // out-of-order CPU without adding a second genome-scale index image.
        for lane in 0..lane_count {
            if let Some(plan) = plans[lane] {
                counters[lane] = self.bwt_word(plan.high_word);
            }
        }
        for boundary in 0..lane_count * 2 {
            if let Some(plan) = boundary_plans[boundary] {
                boundary_counters[boundary] = self.bwt_word(plan.high_word);
            }
        }
        for lane in 0..lane_count {
            let Some(plan) = plans[lane] else {
                continue;
            };
            if plan.lines.iter().any(|line| line & 63 != 0) && plan.digit != 2 {
                let low_block = (plan.lines[0] & 127) >> 6;
                first_planes[lane] = self.bwt_word(plan.high_word + 1 + (low_block << 1));
            }
        }
        for boundary in 0..lane_count * 2 {
            let Some(plan) = boundary_plans[boundary] else {
                continue;
            };
            if plan.line & 63 != 0 && plan.digit != 2 {
                let low_block = (plan.line & 127) >> 6;
                boundary_first_planes[boundary] =
                    self.bwt_word(plan.high_word + 1 + (low_block << 1));
            }
        }
        for lane in 0..lane_count {
            let Some(plan) = plans[lane] else {
                continue;
            };
            if plan.lines.iter().any(|line| line & 63 != 0) && plan.digit != 1 {
                let low_block = (plan.lines[0] & 127) >> 6;
                second_planes[lane] = self.bwt_word(plan.high_word + 2 + (low_block << 1));
            }
        }
        for boundary in 0..lane_count * 2 {
            let Some(plan) = boundary_plans[boundary] else {
                continue;
            };
            if plan.line & 63 != 0 && plan.digit != 1 {
                let low_block = (plan.line & 127) >> 6;
                boundary_second_planes[boundary] =
                    self.bwt_word(plan.high_word + 2 + (low_block << 1));
            }
        }
        for lane in 0..lane_count {
            let Some(plan) = plans[lane] else {
                continue;
            };
            if plan.digit != 2 {
                first_high_occ[lane] = self.high_occ((plan.lines[0] >> 16) << 1);
            }
        }
        for boundary in 0..lane_count * 2 {
            let Some(plan) = boundary_plans[boundary] else {
                continue;
            };
            if plan.digit != 2 {
                boundary_first_high_occ[boundary] = self.high_occ((plan.line >> 16) << 1);
            }
        }
        for lane in 0..lane_count {
            let Some(plan) = plans[lane] else {
                continue;
            };
            if plan.digit != 1 {
                second_high_occ[lane] = self.high_occ(((plan.lines[0] >> 16) << 1) + 1);
            }
        }
        for boundary in 0..lane_count * 2 {
            let Some(plan) = boundary_plans[boundary] else {
                continue;
            };
            if plan.digit != 1 {
                boundary_second_high_occ[boundary] = self.high_occ(((plan.line >> 16) << 1) + 1);
            }
        }
        for lane in 0..lane_count {
            if intervals[lane].is_none() {
                continue;
            }
            let boundaries = if let Some(plan) = plans[lane] {
                self.finish_same_low_block_rank_plan(
                    plan,
                    counters[lane],
                    [first_planes[lane], second_planes[lane]],
                    [first_high_occ[lane], second_high_occ[lane]],
                )
            } else {
                let lower = lane * 2;
                let upper = lower + 1;
                boundary_plans[lower].zip(boundary_plans[upper]).and_then(
                    |(lower_plan, upper_plan)| {
                        Some([
                            self.finish_boundary_rank_plan(
                                lower_plan,
                                boundary_counters[lower],
                                [boundary_first_planes[lower], boundary_second_planes[lower]],
                                [
                                    boundary_first_high_occ[lower],
                                    boundary_second_high_occ[lower],
                                ],
                            )?,
                            self.finish_boundary_rank_plan(
                                upper_plan,
                                boundary_counters[upper],
                                [boundary_first_planes[upper], boundary_second_planes[upper]],
                                [
                                    boundary_first_high_occ[upper],
                                    boundary_second_high_occ[upper],
                                ],
                            )?,
                        ])
                    },
                )
            };
            output[lane] = boundaries.and_then(|[lower, upper]| {
                FmInterval::private_checked(lower, upper, self.suffix_count).ok()
            });
        }
    }

    #[inline]
    fn lf_boundary_pair(&self, lower: u64, upper: u64, digit: u8) -> Option<[u64; 2]> {
        // The symbol is fixed for both boundaries. Monomorphizing this small
        // dispatch removes repeated A/G/T branches from the hot rank kernel.
        match digit {
            0 => self.lf_boundary_pair_digit::<0>(lower, upper),
            1 => self.lf_boundary_pair_digit::<1>(lower, upper),
            2 => self.lf_boundary_pair_digit::<2>(lower, upper),
            _ => None,
        }
    }

    #[inline]
    fn lf_boundary_pair_digit<const DIGIT: u8>(&self, lower: u64, upper: u64) -> Option<[u64; 2]> {
        let digit = DIGIT;
        debug_assert!(lower <= upper);
        debug_assert!(upper <= self.suffix_count);
        debug_assert!(digit <= 2);
        let lower_line = lower - u64::from(lower > self.sentinel_row);
        let upper_line = upper - u64::from(upper > self.sentinel_row);
        if lower_line >> 6 != upper_line >> 6 {
            if lower_line >> 7 == upper_line >> 7 {
                let high_word = (lower_line >> 7).checked_mul(BWT_WORDS_PER_128_ROWS)?;
                let high_occ_block = (lower_line >> 16).checked_mul(2)?;
                let counter_word = self.bwt_word(high_word);
                let absolute = if digit == 0 {
                    self.high_occ(high_occ_block)
                        .checked_add(self.high_occ(high_occ_block + 1))?
                } else {
                    self.high_occ(high_occ_block + u64::from(digit >> 1))
                };
                let rank = |line: u64| {
                    let low_block = (line & 127) >> 6;
                    let at_block = if digit == 0 {
                        let packed = counter_word >> (32 - (low_block << 5));
                        let non_a = ((packed >> 16) & 0xffff) + (packed & 0xffff);
                        ((line >> 6) << 6).checked_sub(absolute.checked_add(non_a)?)?
                    } else {
                        let plane = u64::from(digit >> 1);
                        let shift = (48 - (plane << 4)).checked_sub(low_block << 5)?;
                        absolute + ((counter_word >> shift) & 0xffff)
                    };
                    let plane_start = high_word.checked_add(1 + (low_block << 1))?;
                    let plane_bits = if digit == 0 {
                        !(self.bwt_word(plane_start) | self.bwt_word(plane_start + 1))
                    } else {
                        self.bwt_word(plane_start + u64::from(digit >> 1))
                    };
                    let need = u32::try_from(line & 63).expect("six bits fit u32");
                    let within = if need == 0 {
                        0
                    } else {
                        u64::from((plane_bits >> (64 - need)).count_ones())
                    };
                    self.first_occurrence[usize::from(digit)]
                        .checked_add(at_block)?
                        .checked_add(within)
                };
                return Some([rank(lower_line)?, rank(upper_line)?]);
            }
            return Some([
                self.lf_boundary(lower, digit)?,
                self.lf_boundary(upper, digit)?,
            ]);
        }

        let high_word = (lower_line >> 7).checked_mul(BWT_WORDS_PER_128_ROWS)?;
        let low_block = (lower_line & 127) >> 6;
        let plane_start = high_word.checked_add(1 + (low_block << 1))?;
        let high_occ_block = (lower_line >> 16).checked_mul(2)?;
        let counter_word = self.bwt_word(high_word);
        let at_block = if digit == 0 {
            let g = self.high_occ(high_occ_block);
            let t = self.high_occ(high_occ_block + 1);
            let packed = counter_word >> (32 - (low_block << 5));
            let non_a = ((packed >> 16) & 0xffff) + (packed & 0xffff);
            ((lower_line >> 6) << 6).checked_sub(g.checked_add(t)?.checked_add(non_a)?)?
        } else {
            let plane = u64::from(digit >> 1);
            let absolute = self.high_occ(high_occ_block + plane);
            let shift = (48 - (plane << 4)).checked_sub(low_block << 5)?;
            absolute + ((counter_word >> shift) & 0xffff)
        };
        let plane_bits = if digit == 0 {
            !(self.bwt_word(plane_start) | self.bwt_word(plane_start + 1))
        } else {
            self.bwt_word(plane_start + u64::from(digit >> 1))
        };
        let rank = |line: u64| {
            let need = u32::try_from(line & 63).expect("six bits fit u32");
            let within = if need == 0 {
                0
            } else {
                u64::from((plane_bits >> (64 - need)).count_ones())
            };
            self.first_occurrence[usize::from(digit)]
                .checked_add(at_block)?
                .checked_add(within)
        };
        Some([rank(lower_line)?, rank(upper_line)?])
    }

    #[inline]
    fn lf_row(&self, row: u64) -> Option<u64> {
        if row == self.sentinel_row || row >= self.suffix_count {
            return None;
        }
        let line = if row > self.sentinel_row {
            row - 1
        } else {
            row
        };
        let high_word = (line >> 7).checked_mul(BWT_WORDS_PER_128_ROWS)?;
        let low_block = (line & 127) >> 6;
        let plane_start = high_word.checked_add(1 + (low_block << 1))?;
        let plane0 = self.bwt_word(plane_start);
        let plane1 = self.bwt_word(plane_start + 1);
        let bit = 1_u64 << (63 - (line & 63));
        let digit = if plane0 & bit != 0 {
            1_u8
        } else if plane1 & bit != 0 {
            2_u8
        } else {
            0_u8
        };
        let high_occ_block = (line >> 16).checked_mul(2)?;
        let counter_word = self.bwt_word(high_word);
        let at_block = if digit == 0 {
            let g = self.high_occ(high_occ_block);
            let t = self.high_occ(high_occ_block + 1);
            let packed = counter_word >> (32 - (low_block << 5));
            let non_a = ((packed >> 16) & 0xffff) + (packed & 0xffff);
            ((line >> 6) << 6).checked_sub(g.checked_add(t)?.checked_add(non_a)?)?
        } else {
            let plane = u64::from(digit >> 1);
            let absolute = self.high_occ(high_occ_block + plane);
            let shift = (48 - (plane << 4)).checked_sub(low_block << 5)?;
            absolute + ((counter_word >> shift) & 0xffff)
        };
        let need = u32::try_from(line & 63).expect("six bits fit u32");
        let within = if need == 0 {
            0
        } else if digit == 0 {
            ((!(plane0 | plane1)) >> (64 - need)).count_ones()
        } else if digit == 1 {
            (plane0 >> (64 - need)).count_ones()
        } else {
            (plane1 >> (64 - need)).count_ones()
        };
        self.first_occurrence[usize::from(digit)]
            .checked_add(at_block)?
            .checked_add(u64::from(within))
    }

    #[inline]
    fn sample_ordinal_plan(&self, row: u64) -> Result<SampleOrdinalPlan, CombinedIndexError> {
        if row >= self.suffix_count {
            return Err(CombinedIndexError::Structure(
                "sampled-SA row exceeds suffix domain",
            ));
        }
        let block = row >> 8;
        let within = row & 255;
        let block_word =
            block
                .checked_mul(SA_FLAG_WORDS_PER_256_ROWS)
                .ok_or(CombinedIndexError::Structure(
                    "SA-flag block offset overflow",
                ))?;
        let flag_word =
            block_word
                .checked_add(1 + (within >> 6))
                .ok_or(CombinedIndexError::Structure(
                    "SA-flag word offset overflow",
                ))?;
        let bit = 1_u64 << (63 - (within & 63));
        Ok(SampleOrdinalPlan {
            block_word,
            within,
            flag_word,
            bit,
        })
    }

    #[inline]
    fn sample_ordinal_sparse(&self, row: u64) -> Result<Option<u64>, CombinedIndexError> {
        let plan = self.sample_ordinal_plan(row)?;
        let flags = self.sa_flag_word(plan.flag_word);
        let bit = plan.bit;
        if flags & bit == 0 {
            return Ok(None);
        }
        let mut ordinal = self.sa_flag_word(plan.block_word);
        let full_words = plan.within >> 6;
        for word in 0..full_words {
            ordinal = ordinal
                .checked_add(u64::from(
                    self.sa_flag_word(plan.block_word + 1 + word).count_ones(),
                ))
                .ok_or(CombinedIndexError::Structure("sampled-SA ordinal overflow"))?;
        }
        let prefix = plan.within & 63;
        if prefix != 0 {
            ordinal = ordinal
                .checked_add(u64::from(
                    (flags & (u64::MAX << (64 - prefix))).count_ones(),
                ))
                .ok_or(CombinedIndexError::Structure("sampled-SA ordinal overflow"))?;
        }
        Ok(Some(ordinal))
    }

    #[inline]
    fn sample_rank(&self, boundary: u64) -> Result<u64, CombinedIndexError> {
        if boundary > self.suffix_count {
            return Err(CombinedIndexError::Structure(
                "sampled-SA rank boundary exceeds suffix domain",
            ));
        }
        if boundary == self.suffix_count {
            let last_row = boundary - 1;
            let before_last = self.sample_rank(last_row)?;
            let block_word = (last_row >> 8) * SA_FLAG_WORDS_PER_256_ROWS;
            let within = last_row & 255;
            let flags = self.sa_flag_word(block_word + 1 + (within >> 6));
            let sampled = u64::from(flags & (1_u64 << (63 - (within & 63))) != 0);
            return before_last
                .checked_add(sampled)
                .ok_or(CombinedIndexError::Structure("sampled-SA rank overflow"));
        }
        let block_word = (boundary >> 8) * SA_FLAG_WORDS_PER_256_ROWS;
        let within = boundary & 255;
        let mut ordinal = self.sa_flag_word(block_word);
        let full_words = within >> 6;
        for word in 0..full_words {
            ordinal = ordinal
                .checked_add(u64::from(
                    self.sa_flag_word(block_word + 1 + word).count_ones(),
                ))
                .ok_or(CombinedIndexError::Structure("sampled-SA rank overflow"))?;
        }
        let prefix = within & 63;
        if prefix != 0 {
            let flags = self.sa_flag_word(block_word + 1 + full_words);
            ordinal = ordinal
                .checked_add(u64::from(
                    (flags & (u64::MAX << (64 - prefix))).count_ones(),
                ))
                .ok_or(CombinedIndexError::Structure("sampled-SA rank overflow"))?;
        }
        Ok(ordinal)
    }

    fn locate_row(&self, row: u64) -> Result<(u64, u64), CombinedIndexError> {
        let mut row = row;
        let mut steps = 0_u64;
        loop {
            if row == self.sentinel_row {
                return Ok((steps, steps));
            }
            if let Some(ordinal) = self.sample_ordinal_sparse(row)? {
                let packed = self
                    .sparse_sa(ordinal)
                    .ok_or(CombinedIndexError::Structure(
                        "sampled-SA ordinal exceeds sparse array",
                    ))?;
                let sampled = (packed & u64::from(SA_VALUE_MASK))
                    .checked_mul(SA_STRIDE)
                    .ok_or(CombinedIndexError::Structure(
                        "sampled suffix coordinate overflow",
                    ))?;
                let position = sampled
                    .checked_add(steps)
                    .ok_or(CombinedIndexError::Structure("located coordinate overflow"))?;
                return Ok((position, steps));
            }
            if steps + 1 >= SA_STRIDE {
                return Err(CombinedIndexError::Structure(
                    "LF completion exceeded the declared sampling distance",
                ));
            }
            row = self.lf_row(row).ok_or(CombinedIndexError::Structure(
                "LF boundary calculation failed",
            ))?;
            steps += 1;
        }
    }

    #[inline]
    fn locate_rows_two_lanes(&self, rows: [u64; 2]) -> Result<[(u64, u64); 2], CombinedIndexError> {
        let mut rows = rows;
        let mut steps = [0_u64; 2];
        let mut located = [None, None];
        loop {
            for lane in 0..2 {
                if located[lane].is_some() {
                    continue;
                }
                if rows[lane] == self.sentinel_row {
                    located[lane] = Some((steps[lane], steps[lane]));
                    continue;
                }
                let Some(ordinal) = self.sample_ordinal_sparse(rows[lane])? else {
                    continue;
                };
                let packed = self
                    .sparse_sa(ordinal)
                    .ok_or(CombinedIndexError::Structure(
                        "sampled-SA ordinal exceeds sparse array",
                    ))?;
                let sampled = (packed & u64::from(SA_VALUE_MASK))
                    .checked_mul(SA_STRIDE)
                    .ok_or(CombinedIndexError::Structure(
                        "sampled suffix coordinate overflow",
                    ))?;
                let position = sampled
                    .checked_add(steps[lane])
                    .ok_or(CombinedIndexError::Structure("located coordinate overflow"))?;
                located[lane] = Some((position, steps[lane]));
            }
            if let [Some(first), Some(second)] = located {
                return Ok([first, second]);
            }
            for lane in 0..2 {
                if located[lane].is_some() {
                    continue;
                }
                if steps[lane] + 1 >= SA_STRIDE {
                    return Err(CombinedIndexError::Structure(
                        "LF completion exceeded the declared sampling distance",
                    ));
                }
                rows[lane] = self
                    .lf_row(rows[lane])
                    .ok_or(CombinedIndexError::Structure(
                        "LF boundary calculation failed",
                    ))?;
                steps[lane] += 1;
            }
        }
    }
}

impl PrivateCombinedIndex for CombinedIndex {
    fn reference_len(&self) -> u64 {
        self.reference_length()
    }

    fn exact_interval(
        &self,
        reversed_projected_pattern: &[SearchBase],
    ) -> Result<Option<FmInterval>, CombinedIndexBackendError> {
        self.exact_search(reversed_projected_pattern)
            .ok_or(CombinedIndexBackendError::Interval)
            .map(|interval| (!interval.is_empty()).then_some(interval))
    }

    fn exact_projected_interval(
        &self,
        reversed_projected_pattern: &[ProjectedBase],
    ) -> Result<Option<FmInterval>, CombinedIndexBackendError> {
        let Some((suffix, prefix)) = reversed_projected_pattern
            .len()
            .checked_sub(LOOKUP_BASES)
            .map(|split| reversed_projected_pattern.split_at(split))
        else {
            return Ok(None);
        };
        let Some(mut interval) = self.lookup16_symbols(prefix) else {
            return Ok(None);
        };
        for &symbol in suffix.iter().rev() {
            interval = self
                .backward_extend_validated(interval, symbol.digit())
                .ok_or(CombinedIndexBackendError::Interval)?;
            if interval.is_empty() {
                break;
            }
        }
        Ok((!interval.is_empty()).then_some(interval))
    }

    fn lookup_interval(
        &self,
        projected_suffix: &[SearchBase],
    ) -> Result<Option<FmInterval>, CombinedIndexBackendError> {
        Ok(self
            .lookup16_symbols(projected_suffix)
            .filter(|interval| !interval.is_empty()))
    }

    fn lookup_projected_interval(
        &self,
        projected_suffix: &[ProjectedBase],
    ) -> Result<Option<FmInterval>, CombinedIndexBackendError> {
        Ok(self
            .lookup16_symbols(projected_suffix)
            .filter(|interval| !interval.is_empty()))
    }

    fn backward_extend_interval(
        &self,
        interval: FmInterval,
        symbol: SearchBase,
    ) -> Result<FmInterval, CombinedIndexBackendError> {
        let digit = projected_digit(symbol).ok_or(CombinedIndexBackendError::Interval)?;
        self.backward_extend_validated(interval, digit)
            .ok_or(CombinedIndexBackendError::Interval)
    }

    fn backward_extend_projected_interval(
        &self,
        interval: FmInterval,
        symbol: ProjectedBase,
    ) -> Result<FmInterval, CombinedIndexBackendError> {
        self.backward_extend_validated(interval, symbol.digit())
            .ok_or(CombinedIndexBackendError::Interval)
    }

    fn backward_extend_intervals(
        &self,
        intervals: &[FmInterval],
        symbols: &[SearchBase],
        output: &mut [FmInterval],
    ) -> Result<(), CombinedIndexBackendError> {
        if intervals.len() != symbols.len() || intervals.len() != output.len() {
            return Err(CombinedIndexBackendError::Structure);
        }
        if symbols.len() > MAX_WAVEFRONT_LANES {
            for ((output, &interval), &symbol) in output.iter_mut().zip(intervals).zip(symbols) {
                *output = self.backward_extend_interval(interval, symbol)?;
            }
            return Ok(());
        }
        let mut digits = [0_u8; MAX_WAVEFRONT_LANES];
        for (destination, &symbol) in digits.iter_mut().zip(symbols) {
            *destination = projected_digit(symbol).ok_or(CombinedIndexBackendError::Interval)?;
        }
        self.backward_extend_interval_round(intervals, &digits[..symbols.len()], output)
    }

    fn backward_extend_projected_intervals(
        &self,
        intervals: &[FmInterval],
        symbols: &[ProjectedBase],
        output: &mut [FmInterval],
    ) -> Result<(), CombinedIndexBackendError> {
        if intervals.len() != symbols.len() || intervals.len() != output.len() {
            return Err(CombinedIndexBackendError::Structure);
        }
        if symbols.len() > MAX_WAVEFRONT_LANES {
            for ((output, &interval), &symbol) in output.iter_mut().zip(intervals).zip(symbols) {
                *output = self.backward_extend_projected_interval(interval, symbol)?;
            }
            return Ok(());
        }
        let mut digits = [0_u8; MAX_WAVEFRONT_LANES];
        for (destination, &symbol) in digits.iter_mut().zip(symbols) {
            *destination = symbol.digit();
        }
        self.backward_extend_interval_round(intervals, &digits[..symbols.len()], output)
    }

    fn resolve_projected_suffix_intervals(
        &self,
        patterns: &[&[ProjectedBase]],
        minimum_suffix_bases: usize,
        stop_interval_length: u64,
        output: &mut [Option<(FmInterval, u64)>],
    ) -> Result<(), CombinedIndexBackendError> {
        self.resolve_projected_suffix_intervals_impl(
            patterns,
            minimum_suffix_bases,
            stop_interval_length,
            output,
        )
    }

    fn visit_interval(
        &self,
        interval: FmInterval,
        visitor: &mut dyn FnMut(u64) -> bool,
    ) -> Result<PrivateCombinedLocateMetrics, CombinedIndexBackendError> {
        let metrics = CombinedIndex::visit_interval(self, interval, visitor)
            .map_err(|_| CombinedIndexBackendError::Structure)?;
        Ok(PrivateCombinedLocateMetrics::new(
            metrics.located_rows,
            metrics.lf_steps,
            metrics.rank_operations,
            metrics.interval_nodes,
        ))
    }

    fn visit_intervals_two_lanes_complete(
        &self,
        intervals: [FmInterval; 2],
        visitor: &mut dyn FnMut(usize, u64),
    ) -> Result<[PrivateCombinedLocateMetrics; 2], CombinedIndexBackendError> {
        let metrics = self
            .visit_interval_two_lanes_complete(intervals, visitor)
            .map_err(|_| CombinedIndexBackendError::Structure)?;
        Ok(metrics.map(|metric| {
            PrivateCombinedLocateMetrics::new(
                metric.located_rows,
                metric.lf_steps,
                metric.rank_operations,
                metric.interval_nodes,
            )
        }))
    }
}

fn map_suffix(prefix: &Path, suffix: &str) -> Result<ReadOnlyMapping, CombinedIndexError> {
    let file = File::open(suffixed_path(prefix, suffix))?;
    ReadOnlyMapping::map(&file)
}

fn suffixed_path(prefix: &Path, suffix: &str) -> PathBuf {
    let mut name: OsString = prefix.as_os_str().to_owned();
    name.push(suffix);
    PathBuf::from(name)
}

fn checked_component_end(
    offset: usize,
    elements: u64,
    width: usize,
    file_len: usize,
) -> Result<usize, CombinedIndexError> {
    let elements = usize::try_from(elements)
        .map_err(|_| CombinedIndexError::Structure("component length exceeds usize"))?;
    let bytes = elements
        .checked_mul(width)
        .ok_or(CombinedIndexError::Structure(
            "component byte length overflow",
        ))?;
    let end = offset
        .checked_add(bytes)
        .ok_or(CombinedIndexError::Structure(
            "component end offset overflow",
        ))?;
    if end > file_len {
        return Err(CombinedIndexError::Structure(
            "component extends beyond its file",
        ));
    }
    Ok(end)
}

#[inline]
fn projected_digit(base: SearchBase) -> Option<u8> {
    match base {
        SearchBase::G => Some(0),
        SearchBase::T => Some(1),
        SearchBase::A => Some(2),
        SearchBase::C => None,
    }
}

trait ProjectedQuerySymbol: Copy {
    fn projected_digit(self) -> Option<u8>;
}

impl ProjectedQuerySymbol for SearchBase {
    #[inline]
    fn projected_digit(self) -> Option<u8> {
        projected_digit(self)
    }
}

impl ProjectedQuerySymbol for ProjectedBase {
    #[inline]
    fn projected_digit(self) -> Option<u8> {
        Some(self.digit())
    }
}

fn slice_u32(bytes: &[u8], offset: usize) -> u32 {
    u32::from_le_bytes(
        bytes[offset..offset + 4]
            .try_into()
            .expect("fixed metadata field"),
    )
}

fn metadata_digest(meta: &[u8; META_BYTES]) -> [u8; 32] {
    let mut digest = [0_u8; 32];
    digest.copy_from_slice(&meta[META_DIGEST_OFFSET..META_DIGEST_OFFSET + 32]);
    digest
}

fn slice_u16(bytes: &[u8], offset: usize) -> u16 {
    u16::from_le_bytes(
        bytes[offset..offset + 2]
            .try_into()
            .expect("fixed metadata field"),
    )
}

fn slice_u64(bytes: &[u8], offset: usize) -> u64 {
    u64::from_le_bytes(
        bytes[offset..offset + 8]
            .try_into()
            .expect("fixed metadata field"),
    )
}

#[cfg(test)]
#[path = "../../tests/whitebox/storage_combined.rs"]
mod whitebox_tests;

#[cfg(test)]
#[path = "../../tests/qualification/combined_index.rs"]
mod qualification_tests;
