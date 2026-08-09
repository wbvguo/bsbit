//! Rust-owned publisher for the current combined three-letter index image.
//!
//! libsais supplies only the cyclic BWT and inverse-SA samples. Rust owns the
//! reference projection, Occ64/Occ65536 packing, dense 16-mer lookup, SA16
//! row inversion, validation, and create-only multi-file publication.

use core::ffi::c_int;
#[cfg(test)]
use core::ffi::c_longlong;
use core::fmt;
use core::mem::size_of;
use std::ffi::OsString;
use std::fs::File;
#[cfg(test)]
use std::fs::OpenOptions;
use std::io::{self, BufWriter, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};

use crate::reference::ContigInput;
#[cfg(test)]
use bsbit_core::alphabet::Base;
use bsbit_core::reference::{ReferenceSemanticDigest, ReferenceSemanticDigestBuilder};
use bsbit_io::{CompletedFile, PublicationError, PublishedFile, StagedFile};

use super::combined_blocks::{
    BoundedBwt, BoundedBwtConfig, BoundedBwtError, PackedProjectedText, build_bounded_bwt,
    project_combined_packed_text,
};
#[cfg(test)]
use crate::build::libsais::{libsais_bwt_aux_omp, libsais64_bwt_aux_omp};
use crate::storage::combined::{
    BWT_WORDS_PER_128_ROWS, META_BYTES, META_BYTES_U32, META_DIGEST_OFFSET, META_EXTENSION_MAGIC,
    META_EXTENSION_MAJOR, META_EXTENSION_MINOR, META_EXTENSION_OFFSET, SA_FLAG_WORDS_PER_256_ROWS,
    lf_all_boundaries,
};
use crate::storage::combined::{CombinedIndex, CombinedIndexError, ReadOnlyMapping};

/// Default sparse suffix-array stride used by the qualified low-memory index.
pub const DEFAULT_COMBINED_INDEX_SA_STRIDE: u32 = 16;
const OCC_STRIDE: u32 = 64;
const HIGH_OCC_STRIDE: u32 = 128;
const SA_VALUE_BITS: u32 = 30;
const SA_VALUE_MASK: u64 = (1_u64 << SA_VALUE_BITS) - 1;
const LOOKUP_BASES: usize = 16;
const LOOKUP_KEYS: u64 = 43_046_721;
const LOOKUP_KEYS_USIZE: usize = 43_046_721;
const LOOKUP_ENTRIES: u64 = LOOKUP_KEYS + 1;
const LOOKUP_GAP_BITS: u32 = 4;
const LOOKUP_BOUNDARY_HIGH_MASK: u64 = 0x0fff_ffff;
const IO_BUFFER_BYTES: usize = 8 * 1024 * 1024;
#[cfg(test)]
const RADIX_BITS: u32 = 12;
#[cfg(test)]
const RADIX_BUCKETS: usize = 1 << RADIX_BITS;
#[cfg(test)]
const RADIX_MASK: u64 = (1 << RADIX_BITS) - 1;
#[cfg(test)]
const PARALLEL_RADIX_MIN_ENTRIES: usize = 1_000_000;
const MAX_LOOKUP_TASK_DIGITS: usize = 4;
#[cfg(test)]
const PROJECTION_SEGMENT_BASES: usize = 8 * 1024 * 1024;
/// Qualified default working-memory budget for bounded combined-index SA16 builds.
pub const DEFAULT_COMBINED_INDEX_MEMORY_MIB: u64 = 9_300;
const POW3: [u64; LOOKUP_BASES + 1] = [
    1, 3, 9, 27, 81, 243, 729, 2_187, 6_561, 19_683, 59_049, 177_147, 531_441, 1_594_323,
    4_782_969, 14_348_907, 43_046_721,
];

/// Checked options for one combined-index SA16 build.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct CombinedIndexBuildOptions {
    threads: u32,
    memory_mib: u64,
}

impl CombinedIndexBuildOptions {
    /// Creates options with mandatory structural validation and no expensive
    /// post-build query audit.
    ///
    /// # Errors
    ///
    /// Rejects zero or a count outside libsais' signed configuration domain.
    pub fn new(threads: u32) -> Result<Self, CombinedIndexBuildError> {
        if threads == 0 || c_int::try_from(threads).is_err() {
            return Err(CombinedIndexBuildError::Argument(
                "thread count must fit a positive signed 32-bit integer",
            ));
        }
        Ok(Self {
            threads,
            memory_mib: DEFAULT_COMBINED_INDEX_MEMORY_MIB,
        })
    }

    /// Returns the configured libsais and Rust worker count.
    #[must_use]
    pub const fn threads(self) -> u32 {
        self.threads
    }

    /// Returns the bounded constructor's memory budget in MiB.
    #[must_use]
    pub const fn memory_mib(self) -> u64 {
        self.memory_mib
    }
}

/// Owned deterministic frozen-layout combined projected text.
#[derive(Debug)]
#[cfg(test)]
pub(crate) struct CombinedProjectedText {
    digits: Vec<u8>,
    reference_bases: u64,
    replaced_unknown_bases: u64,
}

/// Projection, combined-index construction, validation, or publication failure.
#[derive(Debug)]
pub enum CombinedIndexBuildError {
    /// Invalid caller configuration or projected-text dimensions.
    Argument(&'static str),
    /// A checked allocation could not be reserved.
    Allocation(&'static str),
    /// A native libsais call returned a nonzero status.
    Libsais(i64),
    /// A generated component violated the frozen format contract.
    Structure(&'static str),
    /// A generated component violated a dynamically identified invariant.
    Detail(String),
    /// The supplied snapshot identity did not describe the source catalog.
    ReferenceDigestMismatch {
        /// Snapshot digest supplied by the caller.
        expected: ReferenceSemanticDigest,
        /// Digest derived directly from the builder's source catalog.
        observed: ReferenceSemanticDigest,
    },
    /// Filesystem I/O failed.
    Io(io::Error),
    /// The ordinary combined-index reader rejected the staged image.
    CombinedIndex(CombinedIndexError),
}

impl fmt::Display for CombinedIndexBuildError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Argument(message) => {
                write!(formatter, "combined-index SA16 build argument: {message}")
            }
            Self::Allocation(label) => {
                write!(formatter, "combined-index SA16 allocation failed: {label}")
            }
            Self::Libsais(status) => write!(formatter, "libsais BWT returned status {status}"),
            Self::Structure(message) => {
                write!(formatter, "combined-index SA16 structure: {message}")
            }
            Self::Detail(message) => write!(formatter, "combined-index SA16 structure: {message}"),
            Self::ReferenceDigestMismatch { expected, observed } => write!(
                formatter,
                "combined-index source catalog digest differs: expected {expected}, observed {observed}"
            ),
            Self::Io(error) => write!(formatter, "combined-index SA16 I/O: {error}"),
            Self::CombinedIndex(error) => error.fmt(formatter),
        }
    }
}

impl std::error::Error for CombinedIndexBuildError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Io(error) => Some(error),
            Self::CombinedIndex(error) => Some(error),
            _ => None,
        }
    }
}

impl From<io::Error> for CombinedIndexBuildError {
    fn from(error: io::Error) -> Self {
        Self::Io(error)
    }
}

impl From<PublicationError> for CombinedIndexBuildError {
    fn from(error: PublicationError) -> Self {
        Self::Io(error.into_io_error())
    }
}

impl From<CombinedIndexError> for CombinedIndexBuildError {
    fn from(error: CombinedIndexError) -> Self {
        Self::CombinedIndex(error)
    }
}

impl From<BoundedBwtError> for CombinedIndexBuildError {
    fn from(error: BoundedBwtError) -> Self {
        Self::Detail(format!("bounded BWT: {error}"))
    }
}

#[derive(Clone, Copy)]
#[cfg(test)]
struct ProjectionSegment<'a> {
    coordinate: usize,
    bases: &'a [Base],
}

/// Projects an ordered reference catalog into the frozen combined G/T/A text.
///
/// The first half is complement-then-C-to-T in forward coordinate order. The
/// second half is C-to-T in reverse order over the complete concatenated
/// catalog. `N` substitutions are deterministic functions of reference salt
/// and global coordinate; the runtime's exact N mask remains authoritative.
///
/// # Errors
///
/// Rejects an empty/oversized catalog and checked allocation failure.
#[cfg(test)]
pub(crate) fn project_combined_text(
    contigs: &[ContigInput],
    projection_salt: u64,
) -> Result<CombinedProjectedText, CombinedIndexBuildError> {
    project_combined_text_with_threads(contigs, projection_salt, 1)
}

/// Projects an ordered reference catalog in parallel without changing the
/// deterministic combined-text bytes produced by the scalar projection.
///
/// # Errors
///
/// Rejects a zero worker count, an empty/oversized catalog, worker panic, or
/// checked allocation failure.
#[allow(clippy::too_many_lines)]
#[cfg(test)]
pub(crate) fn project_combined_text_with_threads(
    contigs: &[ContigInput],
    projection_salt: u64,
    threads: u32,
) -> Result<CombinedProjectedText, CombinedIndexBuildError> {
    if threads == 0 {
        return Err(CombinedIndexBuildError::Argument(
            "projection thread count must be positive",
        ));
    }
    let reference_bases = contigs.iter().try_fold(0_u64, |total, contig| {
        total
            .checked_add(contig.sequence().len())
            .ok_or(CombinedIndexBuildError::Argument(
                "reference length overflows u64",
            ))
    })?;
    if reference_bases == 0 {
        return Err(CombinedIndexBuildError::Argument(
            "reference catalog is empty",
        ));
    }
    let projected_symbols =
        reference_bases
            .checked_mul(2)
            .ok_or(CombinedIndexBuildError::Argument(
                "combined reference length overflows u64",
            ))?;
    validate_combined_text_length(projected_symbols, DEFAULT_COMBINED_INDEX_SA_STRIDE)?;
    let capacity = usize::try_from(projected_symbols)
        .map_err(|_| CombinedIndexBuildError::Argument("combined reference exceeds usize"))?;
    let mut digits = reserved_ffi_vec::<u8>(capacity, "combined projected text")?;
    let reference_length = capacity / 2;
    let estimated_segments = reference_length
        .div_ceil(PROJECTION_SEGMENT_BASES)
        .saturating_add(contigs.len());
    let mut segments = Vec::new();
    segments
        .try_reserve_exact(estimated_segments)
        .map_err(|_| CombinedIndexBuildError::Allocation("projection segment catalog"))?;
    let mut coordinate = 0_usize;
    for contig in contigs {
        for bases in contig.sequence().bases().chunks(PROJECTION_SEGMENT_BASES) {
            segments.push(ProjectionSegment { coordinate, bases });
            coordinate =
                coordinate
                    .checked_add(bases.len())
                    .ok_or(CombinedIndexBuildError::Structure(
                        "projection coordinate overflow",
                    ))?;
        }
    }
    if coordinate != reference_length || segments.is_empty() {
        return Err(CombinedIndexBuildError::Structure(
            "projection segment dimensions disagree with the reference",
        ));
    }

    let next_segment = AtomicUsize::new(0);
    let output_address = digits.as_mut_ptr() as usize;
    let worker_limit = usize::try_from(threads)
        .map_err(|_| CombinedIndexBuildError::Argument("projection threads exceed usize"))?;
    let workers = worker_limit.min(segments.len());
    let replaced_unknown_bases = std::thread::scope(|scope| {
        let mut handles = Vec::with_capacity(workers);
        for _ in 0..workers {
            handles.push(scope.spawn(|| {
                let output = output_address as *mut u8;
                let mut replaced = 0_u64;
                loop {
                    let ordinal = next_segment.fetch_add(1, Ordering::Relaxed);
                    let Some(segment) = segments.get(ordinal) else {
                        break;
                    };
                    for (offset, &base) in segment.bases.iter().enumerate() {
                        let coordinate = segment.coordinate + offset;
                        let canonical =
                            deterministic_canonical_code(base, coordinate, projection_salt);
                        replaced += u64::from(matches!(base, Base::N));
                        // SAFETY: segments partition [0, reference_length)
                        // exactly once. Forward and reverse destinations are
                        // disjoint halves and are unique functions of that
                        // coordinate, so scoped workers never alias writes.
                        unsafe {
                            output
                                .add(coordinate)
                                .write([1, 0, 1, 2][usize::from(canonical)]);
                            output
                                .add(capacity - 1 - coordinate)
                                .write([2, 1, 0, 1][usize::from(canonical)]);
                        }
                    }
                }
                replaced
            }));
        }
        let mut total = 0_u64;
        for handle in handles {
            total = total
                .checked_add(handle.join().map_err(|_| {
                    CombinedIndexBuildError::Structure("projection worker panicked")
                })?)
                .ok_or(CombinedIndexBuildError::Structure(
                    "unknown-base count overflow",
                ))?;
        }
        Ok::<_, CombinedIndexBuildError>(total)
    })?;
    // SAFETY: successful worker completion means every slot in both halves was
    // initialized exactly once, as established by the segment partition above.
    unsafe { digits.set_len(capacity) };
    Ok(CombinedProjectedText {
        digits,
        reference_bases,
        replaced_unknown_bases,
    })
}

#[cfg(test)]
fn deterministic_canonical_code(base: Base, coordinate: usize, projection_salt: u64) -> u8 {
    match base {
        Base::A => 0,
        Base::C => 1,
        Base::G => 2,
        Base::T => 3,
        _ => {
            let mut value = u64::try_from(coordinate)
                .unwrap_or(u64::MAX)
                .wrapping_add(projection_salt)
                .wrapping_add(0x9e37_79b9_7f4a_7c15);
            value = (value ^ (value >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
            value = (value ^ (value >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
            u8::try_from((value ^ (value >> 31)) & 3).expect("two bits fit u8")
        }
    }
}

fn validate_combined_text_length(
    length: u64,
    sa_stride: u32,
) -> Result<(), CombinedIndexBuildError> {
    if length < 2 || !length.is_multiple_of(2) {
        return Err(CombinedIndexBuildError::Argument(
            "combined projected text must have positive even length",
        ));
    }
    if length / u64::from(sa_stride) > SA_VALUE_MASK {
        return Err(CombinedIndexBuildError::Argument(
            "combined projected text exceeds the sparse-SA 30-bit quotient domain",
        ));
    }
    Ok(())
}

#[cfg(test)]
struct DirectBwtBuild {
    transformed: Vec<u8>,
    sentinel_row: u64,
    sampled_rows: Vec<u64>,
}

#[cfg(test)]
fn build_direct_bwt(
    text: Vec<u8>,
    threads: u32,
    sa_stride: u32,
) -> Result<DirectBwtBuild, CombinedIndexBuildError> {
    if text.iter().any(|&digit| digit > 2) {
        return Err(CombinedIndexBuildError::Argument(
            "projected text contains a digit outside G/T/A",
        ));
    }
    if c_int::try_from(text.len()).is_ok() {
        build_direct_bwt32(text, threads, sa_stride)
    } else {
        build_direct_bwt64(text, threads, sa_stride)
    }
}

#[cfg(test)]
fn build_direct_bwt32(
    mut text: Vec<u8>,
    threads: u32,
    sa_stride: u32,
) -> Result<DirectBwtBuild, CombinedIndexBuildError> {
    let length = c_int::try_from(text.len())
        .map_err(|_| CombinedIndexBuildError::Argument("text exceeds libsais32 domain"))?;
    let threads = c_int::try_from(threads)
        .map_err(|_| CombinedIndexBuildError::Argument("thread count exceeds c_int"))?;
    let mut temporary = reserved_ffi_vec::<c_int>(text.len(), "libsais32 temporary array")?;
    let sa_stride_usize = usize::try_from(sa_stride).expect("validated SA stride fits usize");
    let sample_count = text.len().div_ceil(sa_stride_usize);
    let mut raw_sampled_rows = reserved_ffi_vec::<c_int>(sample_count, "libsais32 samples")?;
    let text_pointer = text.as_mut_ptr();
    // SAFETY: the pinned libsais API explicitly permits U == T for an in-place
    // BWT. `text_pointer`, `temporary`, and `raw_sampled_rows` name live arrays
    // with the exact required capacities and do not otherwise overlap.
    let status = unsafe {
        libsais_bwt_aux_omp(
            text_pointer.cast_const(),
            text_pointer,
            temporary.as_mut_ptr(),
            length,
            0,
            core::ptr::null_mut(),
            c_int::try_from(sa_stride).expect("validated SA stride fits c_int"),
            raw_sampled_rows.as_mut_ptr(),
            threads,
        )
    };
    if status != 0 {
        return Err(CombinedIndexBuildError::Libsais(i64::from(status)));
    }
    // SAFETY: a successful libsais BWT auxiliary call initializes every
    // requested inverse-SA sample. The in-place text allocation already has
    // its full initialized length and now contains the transformed bytes.
    unsafe { raw_sampled_rows.set_len(sample_count) };
    drop(temporary);
    let expected = text.len() / sa_stride_usize + 1;
    let mut sampled_rows = Vec::new();
    sampled_rows
        .try_reserve_exact(expected)
        .map_err(|_| CombinedIndexBuildError::Allocation("decoded libsais32 samples"))?;
    for row in raw_sampled_rows {
        sampled_rows.push(u64::try_from(row).map_err(|_| {
            CombinedIndexBuildError::Structure("libsais32 returned a negative sample row")
        })?);
    }
    if text.len().is_multiple_of(sa_stride_usize) {
        sampled_rows.push(0);
    }
    finish_direct_bwt(text, sampled_rows, expected)
}

#[cfg(test)]
fn build_direct_bwt64(
    mut text: Vec<u8>,
    threads: u32,
    sa_stride: u32,
) -> Result<DirectBwtBuild, CombinedIndexBuildError> {
    let length = c_longlong::try_from(text.len())
        .map_err(|_| CombinedIndexBuildError::Argument("text exceeds libsais64 domain"))?;
    let mut temporary = reserved_ffi_vec::<c_longlong>(text.len(), "libsais64 temporary array")?;
    let sa_stride_usize = usize::try_from(sa_stride).expect("validated SA stride fits usize");
    let sample_count = text.len().div_ceil(sa_stride_usize);
    let mut sampled_rows = reserved_ffi_vec::<u64>(sample_count, "libsais64 samples")?;
    let text_pointer = text.as_mut_ptr();
    // SAFETY: u64 and c_longlong have equal size/alignment on the supported
    // Linux target. libsais writes nonnegative row values, validated below;
    // every u64 bit pattern is valid while the call is in progress.
    let status = unsafe {
        libsais64_bwt_aux_omp(
            text_pointer.cast_const(),
            text_pointer,
            temporary.as_mut_ptr(),
            length,
            0,
            core::ptr::null_mut(),
            c_longlong::from(sa_stride),
            sampled_rows.as_mut_ptr().cast::<c_longlong>(),
            c_longlong::from(threads),
        )
    };
    if status != 0 {
        return Err(CombinedIndexBuildError::Libsais(status));
    }
    // SAFETY: on success libsais initializes all requested sample slots. The
    // input allocation remains fully initialized and now contains the in-place
    // BWT. `u64` has no invalid bit patterns; the signed sample-domain contract
    // is checked below before values are used.
    unsafe { sampled_rows.set_len(sample_count) };
    drop(temporary);
    let suffix_count = u64::try_from(text.len())
        .ok()
        .and_then(|length| length.checked_add(1))
        .ok_or(CombinedIndexBuildError::Structure(
            "suffix count overflows u64",
        ))?;
    if sampled_rows.iter().any(|&row| row >= suffix_count) {
        return Err(CombinedIndexBuildError::Structure(
            "libsais64 returned an invalid sample row",
        ));
    }
    let expected = text.len() / sa_stride_usize + 1;
    if text.len().is_multiple_of(sa_stride_usize) {
        sampled_rows
            .try_reserve_exact(1)
            .map_err(|_| CombinedIndexBuildError::Allocation("terminal SA16 sample"))?;
        sampled_rows.push(0);
    }
    finish_direct_bwt(text, sampled_rows, expected)
}

#[cfg(test)]
fn finish_direct_bwt(
    transformed: Vec<u8>,
    sampled_rows: Vec<u64>,
    expected_samples: usize,
) -> Result<DirectBwtBuild, CombinedIndexBuildError> {
    if sampled_rows.len() != expected_samples {
        return Err(CombinedIndexBuildError::Structure(
            "libsais sample count differs from the SA16 contract",
        ));
    }
    let sentinel_row = *sampled_rows
        .first()
        .ok_or(CombinedIndexBuildError::Structure(
            "libsais omitted the sentinel sample",
        ))?;
    Ok(DirectBwtBuild {
        transformed,
        sentinel_row,
        sampled_rows,
    })
}

#[cfg(test)]
fn zeroed_vec<T: Default + Clone>(
    length: usize,
    label: &'static str,
) -> Result<Vec<T>, CombinedIndexBuildError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(length)
        .map_err(|_| CombinedIndexBuildError::Allocation(label))?;
    values.resize(length, T::default());
    Ok(values)
}

#[cfg(test)]
fn reserved_ffi_vec<T>(
    length: usize,
    label: &'static str,
) -> Result<Vec<T>, CombinedIndexBuildError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(length)
        .map_err(|_| CombinedIndexBuildError::Allocation(label))?;
    debug_assert!(values.capacity() >= length);
    Ok(values)
}

#[derive(Debug)]
struct IndexComponentPaths {
    meta: PathBuf,
    bwt: PathBuf,
    sa: PathBuf,
    occ: PathBuf,
}

impl IndexComponentPaths {
    fn from_prefix(prefix: &Path) -> Self {
        Self {
            meta: prefix.to_path_buf(),
            bwt: suffixed_path(prefix, ".bwt"),
            sa: suffixed_path(prefix, ".sa"),
            occ: suffixed_path(prefix, ".occ"),
        }
    }

    fn all(&self) -> [&Path; 4] {
        [&self.meta, &self.bwt, &self.sa, &self.occ]
    }

    fn any_exists(&self) -> io::Result<bool> {
        for path in self.all() {
            match bsbit_io::validate_absent(path) {
                Ok(()) => {}
                Err(error) if error.kind() == io::ErrorKind::AlreadyExists => return Ok(true),
                Err(error) => return Err(error),
            }
        }
        Ok(false)
    }
}

#[derive(Debug)]
struct StagedComponent {
    staging: StagedFile,
    file: File,
}

impl StagedComponent {
    fn create(path: &Path) -> Result<Self, CombinedIndexBuildError> {
        let mut staging = StagedFile::create_new(path)?;
        let file = staging.take_file()?;
        Ok(Self { staging, file })
    }

    fn complete(self) -> Result<CompletedFile, CombinedIndexBuildError> {
        self.staging.complete(self.file).map_err(Into::into)
    }
}

#[derive(Debug)]
struct StagedCombinedIndex {
    target: IndexComponentPaths,
    stage: IndexComponentPaths,
    meta: StagedComponent,
    bwt: StagedComponent,
    sa: StagedComponent,
    occ: StagedComponent,
}

impl StagedCombinedIndex {
    fn create(target_prefix: &Path) -> Result<Self, CombinedIndexBuildError> {
        let target_prefix = bsbit_io::absolute_path(target_prefix)?;
        if target_prefix.file_name().is_none() {
            return Err(CombinedIndexBuildError::Argument(
                "combined index prefix must name a file",
            ));
        }
        let target = IndexComponentPaths::from_prefix(&target_prefix);
        if target.any_exists()? {
            return Err(CombinedIndexBuildError::Argument(
                "combined index target or one of its components already exists",
            ));
        }
        let stage_prefix = bsbit_io::select_sibling_staging_path(&target.meta, "combined-index")?;
        let stage = IndexComponentPaths::from_prefix(&stage_prefix);
        let meta = StagedComponent::create(&stage.meta)?;
        let bwt = StagedComponent::create(&stage.bwt)?;
        let sa = StagedComponent::create(&stage.sa)?;
        let occ = StagedComponent::create(&stage.occ)?;
        Ok(Self {
            target,
            stage,
            meta,
            bwt,
            sa,
            occ,
        })
    }

    fn seal(self) -> Result<CompletedCombinedIndex, CombinedIndexBuildError> {
        let Self {
            target,
            stage,
            meta,
            bwt,
            sa,
            occ,
        } = self;
        let meta = meta.complete()?;
        let bwt = bwt.complete()?;
        let sa = sa.complete()?;
        let occ = occ.complete()?;
        Ok(CompletedCombinedIndex {
            target,
            stage,
            meta,
            bwt,
            sa,
            occ,
        })
    }
}

#[derive(Debug)]
struct CompletedCombinedIndex {
    target: IndexComponentPaths,
    stage: IndexComponentPaths,
    meta: CompletedFile,
    bwt: CompletedFile,
    sa: CompletedFile,
    occ: CompletedFile,
}

impl CompletedCombinedIndex {
    fn staging_identities_match(&self) -> Result<bool, CombinedIndexBuildError> {
        for component in [&self.meta, &self.bwt, &self.sa, &self.occ] {
            if !component.staging_identity_matches()? {
                return Ok(false);
            }
        }
        Ok(true)
    }

    fn publish(self) -> Result<(), CombinedIndexBuildError> {
        // Components become visible first; the versioned metadata file is the
        // commit marker opened first by every combined-index reader.
        let Self {
            target,
            stage: _,
            meta,
            bwt,
            sa,
            occ,
        } = self;
        let ordered = [
            (bwt, target.bwt),
            (sa, target.sa),
            (occ, target.occ),
            (meta, target.meta),
        ];
        let mut published = Vec::<PublishedFile>::with_capacity(ordered.len());
        for (completed, target) in ordered {
            match completed.publish_create_new_at(&target) {
                Ok(component) => published.push(component),
                Err(error) => {
                    for component in published.into_iter().rev() {
                        let _ = component.rollback();
                    }
                    return Err(error.into());
                }
            }
        }
        Ok(())
    }
}

fn suffixed_path(prefix: &Path, suffix: &str) -> PathBuf {
    let mut name = prefix.as_os_str().to_os_string();
    name.push(OsString::from(suffix));
    PathBuf::from(name)
}

#[derive(Clone, Copy, Debug)]
struct BwtDimensions {
    suffix_count: u64,
    sentinel_row: u64,
    first_occurrence: [u64; 4],
    bwt_words: u64,
    high_occ_entries: u64,
}

#[allow(clippy::too_many_lines)]
#[cfg(test)]
fn write_bwt_and_occ(
    transformed: &[u8],
    sentinel_row: u64,
    bwt_path: &Path,
    occ_path: &Path,
) -> Result<BwtDimensions, CombinedIndexBuildError> {
    let bwt_file = create_new_file(bwt_path)?;
    let occ_file = create_new_file(occ_path)?;
    write_bwt_and_occ_from_symbols(
        transformed.len(),
        sentinel_row,
        |line| transformed[line],
        &bwt_file,
        &occ_file,
    )
}

fn write_bounded_bwt_and_occ(
    state: &BoundedBwt,
    bwt_file: &File,
    occ_file: &File,
) -> Result<BwtDimensions, CombinedIndexBuildError> {
    write_bwt_and_occ_from_symbols(
        state.text_len(),
        u64::try_from(state.sentinel_row())
            .map_err(|_| CombinedIndexBuildError::Structure("bounded sentinel row exceeds u64"))?,
        |line| state.transformed_digit(line),
        bwt_file,
        occ_file,
    )
}

#[allow(clippy::too_many_lines)]
fn write_bwt_and_occ_from_symbols(
    text_length_usize: usize,
    sentinel_row: u64,
    symbol: impl Fn(usize) -> u8,
    bwt_file: &File,
    occ_file: &File,
) -> Result<BwtDimensions, CombinedIndexBuildError> {
    let text_length = u64::try_from(text_length_usize)
        .map_err(|_| CombinedIndexBuildError::Structure("BWT text length exceeds u64"))?;
    let suffix_count = text_length
        .checked_add(1)
        .ok_or(CombinedIndexBuildError::Structure(
            "BWT suffix count overflow",
        ))?;
    if sentinel_row >= suffix_count {
        return Err(CombinedIndexBuildError::Structure(
            "sentinel row exceeds suffix domain",
        ));
    }
    let bwt_words = bwt_word_count(text_length).ok_or(CombinedIndexBuildError::Structure(
        "BWT word count overflow",
    ))?;
    let bwt_file = reset_component_file(bwt_file)?;
    let mut bwt = BufWriter::with_capacity(IO_BUFFER_BYTES, bwt_file);
    write_u64(&mut bwt, bwt_words)?;

    let expected_high_occ = (text_length >> 16)
        .checked_add(1)
        .and_then(|blocks| blocks.checked_mul(2))
        .ok_or(CombinedIndexBuildError::Structure(
            "high occurrence count overflow",
        ))?;
    let expected_high_occ_usize = usize::try_from(expected_high_occ)
        .map_err(|_| CombinedIndexBuildError::Structure("high occurrence count exceeds usize"))?;
    let mut high_occ = Vec::new();
    high_occ
        .try_reserve_exact(expected_high_occ_usize)
        .map_err(|_| CombinedIndexBuildError::Allocation("high occurrence table"))?;

    let mut counts = [0_u64; 3];
    let mut checkpoint = [0_u64; 3];
    let mut line = 0_usize;
    let mut words_written = 0_u64;
    while line < text_length_usize {
        if line.is_multiple_of(65_536) {
            checkpoint = counts;
            high_occ.extend_from_slice(&[counts[1], counts[2]]);
        }
        let block_end = line.saturating_add(128).min(text_length_usize);
        let first_end = line.saturating_add(64).min(block_end);
        let block_start_counts = counts;
        let mut first_planes = [0_u64; 2];
        encode_bwt_half_from_symbols(line, first_end, &symbol, &mut first_planes, &mut counts)?;
        let mut counter = checked_delta16(block_start_counts[1], checkpoint[1])? << 48;
        counter |= checked_delta16(block_start_counts[2], checkpoint[2])? << 32;
        if first_end == line.saturating_add(64) {
            counter |= checked_delta16(counts[1], checkpoint[1])? << 16;
            counter |= checked_delta16(counts[2], checkpoint[2])?;
        }
        write_u64(&mut bwt, counter)?;
        write_u64(&mut bwt, first_planes[0])?;
        write_u64(&mut bwt, first_planes[1])?;
        words_written += 3;

        if first_end < block_end {
            let mut second_planes = [0_u64; 2];
            encode_bwt_half_from_symbols(
                first_end,
                block_end,
                &symbol,
                &mut second_planes,
                &mut counts,
            )?;
            write_u64(&mut bwt, second_planes[0])?;
            write_u64(&mut bwt, second_planes[1])?;
            words_written += 2;
        }
        line = block_end;
    }
    if text_length_usize.is_multiple_of(128) {
        if text_length_usize.is_multiple_of(65_536) {
            checkpoint = counts;
            high_occ.extend_from_slice(&[counts[1], counts[2]]);
        }
        let counter = (checked_delta16(counts[1], checkpoint[1])? << 48)
            | (checked_delta16(counts[2], checkpoint[2])? << 32);
        write_u64(&mut bwt, counter)?;
        write_u64(&mut bwt, 0)?;
        write_u64(&mut bwt, 0)?;
        words_written += 3;
    } else if text_length_usize.is_multiple_of(64) {
        // The optimized all-symbol rank reader loads both zero planes even at
        // an exact 64-row boundary. Preserve the two C++ tail padding words.
        write_u64(&mut bwt, 0)?;
        write_u64(&mut bwt, 0)?;
        words_written += 2;
    }
    if words_written != bwt_words || high_occ.len() != expected_high_occ_usize {
        return Err(CombinedIndexBuildError::Structure(
            "BWT or high-occurrence writer dimensions disagree",
        ));
    }
    bwt.flush()?;
    bwt.get_ref().sync_all()?;

    let occ_file = reset_component_file(occ_file)?;
    let mut occ = BufWriter::with_capacity(IO_BUFFER_BYTES, occ_file);
    write_u64(&mut occ, expected_high_occ)?;
    for value in high_occ {
        write_u64(&mut occ, value)?;
    }
    occ.flush()?;
    occ.get_ref().sync_all()?;

    let first_occurrence = [
        1,
        1_u64
            .checked_add(counts[0])
            .ok_or(CombinedIndexBuildError::Structure(
                "first occurrence overflow",
            ))?,
        1_u64
            .checked_add(counts[0])
            .and_then(|value| value.checked_add(counts[1]))
            .ok_or(CombinedIndexBuildError::Structure(
                "first occurrence overflow",
            ))?,
        suffix_count,
    ];
    Ok(BwtDimensions {
        suffix_count,
        sentinel_row,
        first_occurrence,
        bwt_words,
        high_occ_entries: expected_high_occ,
    })
}

fn encode_bwt_half_from_symbols(
    first: usize,
    end: usize,
    symbol: &impl Fn(usize) -> u8,
    planes: &mut [u64; 2],
    counts: &mut [u64; 3],
) -> Result<(), CombinedIndexBuildError> {
    if end < first || end - first > 64 {
        return Err(CombinedIndexBuildError::Structure(
            "BWT half-block exceeds 64 symbols",
        ));
    }
    for (offset, line) in (first..end).enumerate() {
        let digit = symbol(line);
        if digit > 2 {
            return Err(CombinedIndexBuildError::Structure(
                "BWT output contains a non-G/T/A digit",
            ));
        }
        let bit = 1_u64 << (63 - offset);
        if digit == 1 {
            planes[0] |= bit;
        } else if digit == 2 {
            planes[1] |= bit;
        }
        counts[usize::from(digit)] =
            counts[usize::from(digit)]
                .checked_add(1)
                .ok_or(CombinedIndexBuildError::Structure(
                    "BWT symbol count overflow",
                ))?;
    }
    Ok(())
}

fn checked_delta16(value: u64, checkpoint: u64) -> Result<u64, CombinedIndexBuildError> {
    let delta = value
        .checked_sub(checkpoint)
        .ok_or(CombinedIndexBuildError::Structure(
            "BWT occurrence checkpoint underflow",
        ))?;
    if delta > u64::from(u16::MAX) {
        return Err(CombinedIndexBuildError::Structure(
            "BWT occurrence delta exceeds 16 bits",
        ));
    }
    Ok(delta)
}

fn bwt_word_count(text_length: u64) -> Option<u64> {
    let complete = (text_length >> 7).checked_mul(BWT_WORDS_PER_128_ROWS)?;
    complete.checked_add(if text_length & 127 < 64 { 3 } else { 5 })
}

#[cfg(test)]
fn create_new_file(path: &Path) -> Result<File, CombinedIndexBuildError> {
    bsbit_io::create_new(path)
        .map(|(file, _)| file)
        .map_err(Into::into)
}

fn reset_component_file(file: &File) -> Result<File, CombinedIndexBuildError> {
    let mut file = bsbit_io::reopen_read_write(file)?;
    file.set_len(0)?;
    file.seek(SeekFrom::Start(0))?;
    Ok(file)
}

fn write_u32(writer: &mut impl Write, value: u32) -> io::Result<()> {
    writer.write_all(&value.to_le_bytes())
}

fn write_u64(writer: &mut impl Write, value: u64) -> io::Result<()> {
    writer.write_all(&value.to_le_bytes())
}

struct BuildRank {
    bwt: ReadOnlyMapping,
    occ: ReadOnlyMapping,
    dimensions: BwtDimensions,
}

impl BuildRank {
    fn open(
        bwt_file: &File,
        occ_file: &File,
        dimensions: BwtDimensions,
    ) -> Result<Self, CombinedIndexBuildError> {
        let bwt = ReadOnlyMapping::map(bwt_file)?;
        let occ = ReadOnlyMapping::map(occ_file)?;
        if bwt.read_u64(0) != dimensions.bwt_words || occ.read_u64(0) != dimensions.high_occ_entries
        {
            return Err(CombinedIndexBuildError::Structure(
                "staged rank headers disagree with writer dimensions",
            ));
        }
        Ok(Self {
            bwt,
            occ,
            dimensions,
        })
    }

    #[inline]
    fn bwt_word(&self, ordinal: u64) -> u64 {
        debug_assert!(ordinal < self.dimensions.bwt_words);
        self.bwt.read_u64(
            8 + usize::try_from(ordinal).expect("BWT ordinal fits usize") * size_of::<u64>(),
        )
    }

    #[inline]
    fn high_occ(&self, ordinal: u64) -> u64 {
        debug_assert!(ordinal < self.dimensions.high_occ_entries);
        self.occ.read_u64(
            8 + usize::try_from(ordinal).expect("Occ ordinal fits usize") * size_of::<u64>(),
        )
    }

    #[inline]
    fn all_boundaries(&self, boundary: u64) -> Option<[u64; 3]> {
        lf_all_boundaries(
            boundary,
            self.dimensions.suffix_count,
            self.dimensions.sentinel_row,
            self.dimensions.first_occurrence,
            |ordinal| self.bwt_word(ordinal),
            |ordinal| self.high_occ(ordinal),
        )
    }
}

#[cfg(test)]
#[derive(Clone, Copy, Debug)]
struct SaDimensions {
    sparse_entries: u64,
    flag_entries: u64,
}

#[allow(clippy::too_many_lines)]
#[cfg(test)]
fn write_sa16(
    sampled_rows: &mut Vec<u64>,
    suffix_count: u64,
    threads: u32,
    sa_stride: u32,
    sa_path: &Path,
) -> Result<SaDimensions, CombinedIndexBuildError> {
    let expected_samples = (suffix_count - 1) / u64::from(sa_stride) + 1;
    if u64::try_from(sampled_rows.len()).ok() != Some(expected_samples) {
        return Err(CombinedIndexBuildError::Structure(
            "inverse-SA sample count disagrees with SA16",
        ));
    }
    for (quotient, row) in sampled_rows.iter_mut().enumerate() {
        let quotient = u64::try_from(quotient)
            .map_err(|_| CombinedIndexBuildError::Structure("SA16 quotient exceeds u64"))?;
        if quotient > SA_VALUE_MASK || *row >= suffix_count {
            return Err(CombinedIndexBuildError::Structure(
                "inverse-SA sample exceeds the combined row/value domain",
            ));
        }
        *row = row
            .checked_shl(SA_VALUE_BITS)
            .and_then(|value| value.checked_add(quotient))
            .ok_or(CombinedIndexBuildError::Structure(
                "packed SA16 sample overflow",
            ))?;
    }
    sort_packed_samples_by_row(sampled_rows, suffix_count - 1, threads)?;
    let mut previous = None;
    for &packed in sampled_rows.iter() {
        let row = packed >> SA_VALUE_BITS;
        if previous.is_some_and(|value| value >= row) {
            return Err(CombinedIndexBuildError::Structure(
                "inverse-SA sample rows are not unique",
            ));
        }
        previous = Some(row);
    }

    let full_blocks = suffix_count / 256;
    let tail_rows = suffix_count % 256;
    let tail_flag_words = tail_rows.div_ceil(64);
    // The current layout leaves a cumulative-count word at the boundary following the
    // final complete 256-row block, then one zero guard word.  A partial final
    // block stores only the flag words that contain rows; it is not padded to
    // four words.  Retaining that exact tail layout makes the Rust publisher's
    // SA16 image byte-compatible with the current reader.
    let flag_entries = full_blocks
        .checked_mul(SA_FLAG_WORDS_PER_256_ROWS)
        .and_then(|entries| entries.checked_add(tail_flag_words))
        .and_then(|entries| entries.checked_add(2))
        .ok_or(CombinedIndexBuildError::Structure(
            "SA16 flag count overflow",
        ))?;
    let sparse_bytes = expected_samples
        .checked_mul(size_of::<u32>() as u64)
        .ok_or(CombinedIndexBuildError::Structure(
            "SA16 sparse byte count overflow",
        ))?;
    let flag_bytes = flag_entries.checked_mul(size_of::<u64>() as u64).ok_or(
        CombinedIndexBuildError::Structure("SA16 flag byte count overflow"),
    )?;
    let flag_header_offset = (size_of::<u64>() as u64).checked_add(sparse_bytes).ok_or(
        CombinedIndexBuildError::Structure("SA16 flag header offset overflow"),
    )?;
    let flags_offset = flag_header_offset
        .checked_add(size_of::<u64>() as u64)
        .ok_or(CombinedIndexBuildError::Structure(
            "SA16 flag offset overflow",
        ))?;
    let final_bytes =
        flags_offset
            .checked_add(flag_bytes)
            .ok_or(CombinedIndexBuildError::Structure(
                "SA16 file length overflow",
            ))?;

    let mut header = create_new_file(sa_path)?;
    write_u64(&mut header, expected_samples)?;
    header.seek(SeekFrom::Start(flag_header_offset))?;
    write_u64(&mut header, flag_entries)?;
    header.set_len(final_bytes)?;

    // Values and flags occupy disjoint final-file ranges. Independent opens
    // give each buffered writer its own cursor, avoiding a temporary flags
    // file and the subsequent ~1 GiB copy for a human-sized reference.
    let mut sparse_file = OpenOptions::new().write(true).open(sa_path)?;
    sparse_file.seek(SeekFrom::Start(size_of::<u64>() as u64))?;
    let mut sa = BufWriter::with_capacity(IO_BUFFER_BYTES, sparse_file);
    let mut flags_file = OpenOptions::new().write(true).open(sa_path)?;
    flags_file.seek(SeekFrom::Start(flags_offset))?;
    let mut flags = BufWriter::with_capacity(IO_BUFFER_BYTES, flags_file);
    let mut sample = 0_usize;
    for block in 0..full_blocks {
        write_u64(
            &mut flags,
            u64::try_from(sample).expect("sample index fits u64"),
        )?;
        let row_start = block * 256;
        let row_end = row_start + 256;
        let mut words = [0_u64; 4];
        while let Some(&packed) = sampled_rows.get(sample) {
            let row = packed >> SA_VALUE_BITS;
            if row >= row_end {
                break;
            }
            if row < row_start {
                return Err(CombinedIndexBuildError::Structure(
                    "SA16 radix order regressed across a flag block",
                ));
            }
            let within = row - row_start;
            let word = usize::try_from(within >> 6).expect("SA flag word fits usize");
            words[word] |= 1_u64 << (63 - (within & 63));
            write_u32(
                &mut sa,
                u32::try_from(packed & SA_VALUE_MASK).expect("30-bit SA16 quotient fits u32"),
            )?;
            sample += 1;
        }
        for word in words {
            write_u64(&mut flags, word)?;
        }
    }

    // The cumulative count at the final full-block boundary is present even
    // when there is no partial block.  It is the prefix word for the partial
    // block when tail_rows is non-zero.
    write_u64(
        &mut flags,
        u64::try_from(sample).expect("sample index fits u64"),
    )?;
    if tail_rows != 0 {
        let row_start = full_blocks * 256;
        let mut words = [0_u64; 4];
        while let Some(&packed) = sampled_rows.get(sample) {
            let row = packed >> SA_VALUE_BITS;
            if row < row_start || row >= suffix_count {
                return Err(CombinedIndexBuildError::Structure(
                    "SA16 radix order regressed in the final flag block",
                ));
            }
            let within = row - row_start;
            let word = usize::try_from(within >> 6).expect("SA flag word fits usize");
            words[word] |= 1_u64 << (63 - (within & 63));
            write_u32(
                &mut sa,
                u32::try_from(packed & SA_VALUE_MASK).expect("30-bit SA16 quotient fits u32"),
            )?;
            sample += 1;
        }
        for &word in words
            .iter()
            .take(usize::try_from(tail_flag_words).expect("tail flag word count fits usize"))
        {
            write_u64(&mut flags, word)?;
        }
    }
    write_u64(&mut flags, 0)?;
    if sample != sampled_rows.len() {
        return Err(CombinedIndexBuildError::Structure(
            "SA16 flags did not consume every sample",
        ));
    }
    sa.flush()?;
    flags.flush()?;
    drop(sa);
    drop(flags);
    header.sync_all()?;
    Ok(SaDimensions {
        sparse_entries: expected_samples,
        flag_entries,
    })
}

#[allow(clippy::too_many_lines)]
fn write_bounded_sa16(
    state: &BoundedBwt,
    sa_stride: u32,
    sa_file: &File,
) -> Result<(), CombinedIndexBuildError> {
    if state.sample_stride() != usize::try_from(sa_stride).expect("validated SA stride fits usize")
    {
        return Err(CombinedIndexBuildError::Structure(
            "bounded BWT and serialized SA strides disagree",
        ));
    }
    let suffix_count = u64::try_from(state.rows())
        .map_err(|_| CombinedIndexBuildError::Structure("bounded suffix count exceeds u64"))?;
    let expected_samples = (suffix_count - 1) / u64::from(sa_stride) + 1;
    if u64::try_from(state.sample_quotients().len()).ok() != Some(expected_samples) {
        return Err(CombinedIndexBuildError::Structure(
            "bounded sample count disagrees with SA16",
        ));
    }
    let full_blocks = suffix_count / 256;
    let tail_rows = suffix_count % 256;
    let tail_flag_words = tail_rows.div_ceil(64);
    let flag_entries = full_blocks
        .checked_mul(SA_FLAG_WORDS_PER_256_ROWS)
        .and_then(|entries| entries.checked_add(tail_flag_words))
        .and_then(|entries| entries.checked_add(2))
        .ok_or(CombinedIndexBuildError::Structure(
            "SA16 flag count overflow",
        ))?;
    let sparse_bytes = expected_samples
        .checked_mul(size_of::<u32>() as u64)
        .ok_or(CombinedIndexBuildError::Structure(
            "SA16 sparse byte count overflow",
        ))?;
    let flag_bytes = flag_entries.checked_mul(size_of::<u64>() as u64).ok_or(
        CombinedIndexBuildError::Structure("SA16 flag byte count overflow"),
    )?;
    let flag_header_offset = (size_of::<u64>() as u64).checked_add(sparse_bytes).ok_or(
        CombinedIndexBuildError::Structure("SA16 flag header offset overflow"),
    )?;
    let flags_offset = flag_header_offset
        .checked_add(size_of::<u64>() as u64)
        .ok_or(CombinedIndexBuildError::Structure(
            "SA16 flag offset overflow",
        ))?;
    let final_bytes =
        flags_offset
            .checked_add(flag_bytes)
            .ok_or(CombinedIndexBuildError::Structure(
                "SA16 file length overflow",
            ))?;

    let mut header = reset_component_file(sa_file)?;
    write_u64(&mut header, expected_samples)?;
    header.seek(SeekFrom::Start(flag_header_offset))?;
    write_u64(&mut header, flag_entries)?;
    header.set_len(final_bytes)?;

    let mut sparse_file = bsbit_io::reopen_read_write(sa_file)?;
    sparse_file.seek(SeekFrom::Start(size_of::<u64>() as u64))?;
    let mut sa = BufWriter::with_capacity(IO_BUFFER_BYTES, sparse_file);
    let mut flags_file = bsbit_io::reopen_read_write(sa_file)?;
    flags_file.seek(SeekFrom::Start(flags_offset))?;
    let mut flags = BufWriter::with_capacity(IO_BUFFER_BYTES, flags_file);

    let mut samples = state.row_ordered_samples().peekable();
    let mut sample_ordinal = 0_u64;
    let mut previous_row = None;
    for block in 0..full_blocks {
        write_u64(&mut flags, sample_ordinal)?;
        let row_start = block * 256;
        let row_end = row_start + 256;
        let mut words = [0_u64; 4];
        while let Some(&(row, quotient)) = samples.peek() {
            if row >= row_end {
                break;
            }
            if row < row_start || previous_row.is_some_and(|previous| previous >= row) {
                return Err(CombinedIndexBuildError::Structure(
                    "bounded SA16 row order regressed",
                ));
            }
            samples.next();
            let within = row - row_start;
            words[usize::try_from(within >> 6).expect("SA flag word fits usize")] |=
                1_u64 << (63 - (within & 63));
            if u64::from(quotient) > SA_VALUE_MASK {
                return Err(CombinedIndexBuildError::Structure(
                    "bounded SA16 quotient exceeds 30 bits",
                ));
            }
            write_u32(&mut sa, quotient)?;
            sample_ordinal += 1;
            previous_row = Some(row);
        }
        for word in words {
            write_u64(&mut flags, word)?;
        }
    }

    write_u64(&mut flags, sample_ordinal)?;
    if tail_rows != 0 {
        let row_start = full_blocks * 256;
        let mut words = [0_u64; 4];
        for (row, quotient) in samples.by_ref() {
            if row < row_start
                || row >= suffix_count
                || previous_row.is_some_and(|previous| previous >= row)
            {
                return Err(CombinedIndexBuildError::Structure(
                    "bounded SA16 final row order regressed",
                ));
            }
            let within = row - row_start;
            words[usize::try_from(within >> 6).expect("SA flag word fits usize")] |=
                1_u64 << (63 - (within & 63));
            if u64::from(quotient) > SA_VALUE_MASK {
                return Err(CombinedIndexBuildError::Structure(
                    "bounded SA16 quotient exceeds 30 bits",
                ));
            }
            write_u32(&mut sa, quotient)?;
            sample_ordinal += 1;
            previous_row = Some(row);
        }
        for &word in words
            .iter()
            .take(usize::try_from(tail_flag_words).expect("tail flag word count fits usize"))
        {
            write_u64(&mut flags, word)?;
        }
    }
    write_u64(&mut flags, 0)?;
    if sample_ordinal != expected_samples || samples.next().is_some() {
        return Err(CombinedIndexBuildError::Structure(
            "bounded SA16 writer did not consume every sample",
        ));
    }
    sa.flush()?;
    flags.flush()?;
    drop(sa);
    drop(flags);
    header.sync_all()?;
    Ok(())
}

#[cfg(test)]
fn sort_packed_samples_by_row(
    samples: &mut Vec<u64>,
    maximum_row: u64,
    threads: u32,
) -> Result<(), CombinedIndexBuildError> {
    if samples.len() < PARALLEL_RADIX_MIN_ENTRIES || threads == 1 {
        samples.sort_unstable_by_key(|value| *value >> SA_VALUE_BITS);
        return Ok(());
    }
    let mut scratch = zeroed_vec::<u64>(samples.len(), "SA16 radix scratch")?;
    let row_bits = (u64::BITS - maximum_row.leading_zeros()).max(1);
    let passes = row_bits.div_ceil(RADIX_BITS);
    let workers = usize::try_from(threads)
        .expect("validated thread count fits usize")
        .min(samples.len());
    for pass in 0..passes {
        let shift = SA_VALUE_BITS + pass * RADIX_BITS;
        parallel_radix_pass(samples, &mut scratch, shift, workers);
        core::mem::swap(samples, &mut scratch);
    }
    Ok(())
}

#[cfg(test)]
fn parallel_radix_pass(source: &[u64], destination: &mut [u64], shift: u32, workers: usize) {
    debug_assert_eq!(source.len(), destination.len());
    let chunk = source.len().div_ceil(workers);
    let histograms = std::thread::scope(|scope| {
        let mut handles = Vec::with_capacity(workers);
        for worker in 0..workers {
            let start = worker * chunk;
            let end = start.saturating_add(chunk).min(source.len());
            let input = &source[start..end];
            handles.push(scope.spawn(move || {
                let mut histogram = vec![0_usize; RADIX_BUCKETS];
                for &value in input {
                    let bucket = usize::try_from((value >> shift) & RADIX_MASK)
                        .expect("radix bucket fits usize");
                    histogram[bucket] += 1;
                }
                histogram
            }));
        }
        handles
            .into_iter()
            .map(|handle| handle.join().expect("SA16 radix counter did not panic"))
            .collect::<Vec<_>>()
    });
    let mut offsets = vec![vec![0_usize; RADIX_BUCKETS]; workers];
    let mut bucket_start = 0_usize;
    for bucket in 0..RADIX_BUCKETS {
        let mut cursor = bucket_start;
        for worker in 0..workers {
            offsets[worker][bucket] = cursor;
            cursor += histograms[worker][bucket];
        }
        bucket_start = cursor;
    }
    debug_assert_eq!(bucket_start, source.len());

    let destination_address = destination.as_mut_ptr() as usize;
    std::thread::scope(|scope| {
        let mut handles = Vec::with_capacity(workers);
        for (worker, mut cursor) in offsets.into_iter().enumerate() {
            let start = worker * chunk;
            let end = start.saturating_add(chunk).min(source.len());
            let input = &source[start..end];
            handles.push(scope.spawn(move || {
                let output = destination_address as *mut u64;
                for &value in input {
                    let bucket = usize::try_from((value >> shift) & RADIX_MASK)
                        .expect("radix bucket fits usize");
                    let position = cursor[bucket];
                    cursor[bucket] += 1;
                    // SAFETY: per-worker bucket starts incorporate all counts
                    // from preceding workers. Every source element therefore
                    // owns one distinct in-bounds destination position.
                    unsafe { output.add(position).write(value) };
                }
            }));
        }
        for handle in handles {
            handle.join().expect("SA16 radix scatter did not panic");
        }
    });
}

trait LookupCount: Copy + Default + Send + 'static {
    fn from_width(width: u64) -> Result<Self, CombinedIndexBuildError>;
    fn as_u64(self) -> u64;
}

impl LookupCount for u32 {
    fn from_width(width: u64) -> Result<Self, CombinedIndexBuildError> {
        u32::try_from(width)
            .map_err(|_| CombinedIndexBuildError::Structure("one exact 16-mer count exceeds u32"))
    }

    fn as_u64(self) -> u64 {
        u64::from(self)
    }
}

impl LookupCount for u64 {
    fn from_width(width: u64) -> Result<Self, CombinedIndexBuildError> {
        Ok(width)
    }

    fn as_u64(self) -> u64 {
        self
    }
}

enum LookupCounts {
    Narrow(Vec<Vec<u32>>),
    Wide(Vec<Vec<u64>>),
}

impl LookupCounts {
    fn write(
        &self,
        short_suffix_thresholds: &[u64],
        text_length: u64,
        bwt_words: u64,
        bwt_file: &File,
    ) -> Result<(), CombinedIndexBuildError> {
        match self {
            Self::Narrow(tasks) => write_lookup16(
                tasks,
                short_suffix_thresholds,
                text_length,
                bwt_words,
                bwt_file,
            ),
            Self::Wide(tasks) => write_lookup16(
                tasks,
                short_suffix_thresholds,
                text_length,
                bwt_words,
                bwt_file,
            ),
        }
    }
}

fn build_lookup_counts(
    rank: &BuildRank,
    threads: u32,
) -> Result<LookupCounts, CombinedIndexBuildError> {
    if u32::try_from(rank.dimensions.suffix_count).is_ok() {
        build_lookup_counts_typed::<u32>(rank, threads).map(LookupCounts::Narrow)
    } else {
        build_lookup_counts_typed::<u64>(rank, threads).map(LookupCounts::Wide)
    }
}

fn build_lookup_counts_typed<T: LookupCount>(
    rank: &BuildRank,
    threads: u32,
) -> Result<Vec<Vec<T>>, CombinedIndexBuildError> {
    let worker_limit = usize::try_from(threads).expect("validated thread count fits usize");
    let desired_tasks = worker_limit.saturating_mul(4).max(1);
    let mut task_digits = 0_usize;
    let mut task_count = 1_usize;
    while task_count < desired_tasks && task_digits < MAX_LOOKUP_TASK_DIGITS {
        task_digits += 1;
        task_count *= 3;
    }
    let keys = LOOKUP_KEYS_USIZE;
    let entries_per_task = keys / task_count;
    debug_assert_eq!(entries_per_task * task_count, keys);
    let next_task = AtomicUsize::new(0);
    let workers = worker_limit.min(task_count);
    let worker_results = std::thread::scope(|scope| {
        let mut handles = Vec::with_capacity(workers);
        for _ in 0..workers {
            handles.push(scope.spawn(|| {
                let mut completed = Vec::new();
                loop {
                    let task = next_task.fetch_add(1, Ordering::Relaxed);
                    if task >= task_count {
                        break;
                    }
                    let counts = build_lookup_task::<T>(
                        rank,
                        task,
                        task_digits,
                        task_count,
                        entries_per_task,
                    )?;
                    completed.push((task, counts));
                }
                Ok::<_, CombinedIndexBuildError>(completed)
            }));
        }
        let mut results = Vec::with_capacity(workers);
        for handle in handles {
            results.push(
                handle
                    .join()
                    .map_err(|_| CombinedIndexBuildError::Structure("lookup worker panicked"))??,
            );
        }
        Ok::<_, CombinedIndexBuildError>(results)
    })?;
    let mut ordered = (0..task_count).map(|_| None).collect::<Vec<_>>();
    for worker in worker_results {
        for (task, counts) in worker {
            if ordered[task].replace(counts).is_some() {
                return Err(CombinedIndexBuildError::Structure(
                    "lookup task completed more than once",
                ));
            }
        }
    }
    ordered
        .into_iter()
        .map(|counts| {
            counts.ok_or(CombinedIndexBuildError::Structure(
                "lookup task did not complete",
            ))
        })
        .collect()
}

fn build_lookup_task<T: LookupCount>(
    rank: &BuildRank,
    task: usize,
    task_digits: usize,
    task_count: usize,
    entries_per_task: usize,
) -> Result<Vec<T>, CombinedIndexBuildError> {
    let mut counts = Vec::new();
    counts
        .try_reserve_exact(entries_per_task)
        .map_err(|_| CombinedIndexBuildError::Allocation("dense 16-mer counts"))?;
    counts.resize(entries_per_task, T::default());
    let mut lower = 0_u64;
    let mut upper = rank.dimensions.suffix_count;
    let mut digits = task;
    for _ in 0..task_digits {
        let digit = digits % 3;
        digits /= 3;
        let lower_children =
            rank.all_boundaries(lower)
                .ok_or(CombinedIndexBuildError::Structure(
                    "lookup lower FM boundary failed",
                ))?;
        let upper_children =
            rank.all_boundaries(upper)
                .ok_or(CombinedIndexBuildError::Structure(
                    "lookup upper FM boundary failed",
                ))?;
        lower = lower_children[digit];
        upper = upper_children[digit];
    }
    enumerate_lookup_subtree(
        rank,
        task_digits,
        u64::try_from(task).expect("lookup task fits u64"),
        u64::try_from(task_count).expect("lookup task count fits u64"),
        lower,
        upper,
        &mut counts,
    )?;
    Ok(counts)
}

fn enumerate_lookup_subtree<T: LookupCount>(
    rank: &BuildRank,
    depth: usize,
    key: u64,
    task_count: u64,
    lower: u64,
    upper: u64,
    counts: &mut [T],
) -> Result<(), CombinedIndexBuildError> {
    if lower > upper {
        return Err(CombinedIndexBuildError::Structure(
            "lookup FM interval is reversed",
        ));
    }
    if lower == upper {
        return Ok(());
    }
    if depth == LOOKUP_BASES {
        let slot = usize::try_from(key / task_count)
            .map_err(|_| CombinedIndexBuildError::Structure("lookup output slot exceeds usize"))?;
        *counts
            .get_mut(slot)
            .ok_or(CombinedIndexBuildError::Structure(
                "lookup output slot is out of range",
            ))? = T::from_width(upper - lower)?;
        return Ok(());
    }
    let lower_children = rank
        .all_boundaries(lower)
        .ok_or(CombinedIndexBuildError::Structure(
            "lookup lower FM boundary failed",
        ))?;
    let upper_children = rank
        .all_boundaries(upper)
        .ok_or(CombinedIndexBuildError::Structure(
            "lookup upper FM boundary failed",
        ))?;
    for digit in 0_usize..3 {
        let child_key = key
            .checked_add(
                u64::try_from(digit)
                    .expect("projected digit fits u64")
                    .checked_mul(POW3[depth])
                    .ok_or(CombinedIndexBuildError::Structure(
                        "lookup key multiplication overflow",
                    ))?,
            )
            .ok_or(CombinedIndexBuildError::Structure(
                "lookup key addition overflow",
            ))?;
        enumerate_lookup_subtree(
            rank,
            depth + 1,
            child_key,
            task_count,
            lower_children[digit],
            upper_children[digit],
            counts,
        )?;
    }
    Ok(())
}

fn short_suffix_thresholds(text: &[u8]) -> Result<Vec<u64>, CombinedIndexBuildError> {
    let maximum = text.len().min(LOOKUP_BASES - 1);
    let mut thresholds = Vec::with_capacity(maximum);
    for length in 1..=maximum {
        let mut code = 0_u64;
        for &digit in &text[text.len() - length..] {
            if digit > 2 {
                return Err(CombinedIndexBuildError::Structure(
                    "short suffix contains a non-G/T/A digit",
                ));
            }
            code = code * 3 + u64::from(digit);
        }
        thresholds.push(code * POW3[LOOKUP_BASES - length]);
    }
    thresholds.sort_unstable();
    Ok(thresholds)
}

fn count_at<T: LookupCount>(tasks: &[Vec<T>], key: usize) -> u64 {
    let task = key % tasks.len();
    let prefix = key / tasks.len();
    tasks[task][prefix].as_u64()
}

fn visit_lookup_intervals<T: LookupCount>(
    tasks: &[Vec<T>],
    short_suffix_thresholds: &[u64],
    mut visit: impl FnMut(u64, u64, u64, u64, u8) -> Result<(), CombinedIndexBuildError>,
) -> Result<(u64, u64), CombinedIndexBuildError> {
    let mut cumulative_long = 0_u64;
    let mut previous_end = 1_u64;
    let mut short_cursor = 0_usize;
    for key in 0..LOOKUP_KEYS_USIZE {
        let count = count_at(tasks, key);
        let key_u64 = u64::try_from(key).expect("lookup key fits u64");
        let (lower, upper, gap) = if count == 0 {
            (previous_end, previous_end, 0_u8)
        } else {
            while short_suffix_thresholds
                .get(short_cursor)
                .is_some_and(|&threshold| threshold <= key_u64)
            {
                short_cursor += 1;
            }
            let lower = 1_u64
                .checked_add(cumulative_long)
                .and_then(|value| value.checked_add(short_cursor as u64))
                .ok_or(CombinedIndexBuildError::Structure(
                    "lookup lower boundary overflow",
                ))?;
            let gap = lower
                .checked_sub(previous_end)
                .ok_or(CombinedIndexBuildError::Structure(
                    "lookup boundaries are not monotone",
                ))?;
            if gap >= 1_u64 << LOOKUP_GAP_BITS {
                return Err(CombinedIndexBuildError::Structure(
                    "lookup short-suffix gap exceeds four bits",
                ));
            }
            let upper = lower
                .checked_add(count)
                .ok_or(CombinedIndexBuildError::Structure(
                    "lookup upper boundary overflow",
                ))?;
            previous_end = upper;
            (
                lower,
                upper,
                u8::try_from(gap).expect("four-bit lookup gap fits u8"),
            )
        };
        visit(key_u64, lower, upper, lower, gap)?;
        cumulative_long =
            cumulative_long
                .checked_add(count)
                .ok_or(CombinedIndexBuildError::Structure(
                    "lookup occurrence total overflow",
                ))?;
    }
    Ok((previous_end, cumulative_long))
}

fn write_lookup16<T: LookupCount>(
    tasks: &[Vec<T>],
    short_suffix_thresholds: &[u64],
    text_length: u64,
    bwt_words: u64,
    bwt_file: &File,
) -> Result<(), CombinedIndexBuildError> {
    if tasks.is_empty()
        || tasks
            .iter()
            .any(|task| task.len().checked_mul(tasks.len()) != Some(LOOKUP_KEYS_USIZE))
    {
        return Err(CombinedIndexBuildError::Structure(
            "lookup task dimensions are inconsistent",
        ));
    }
    let expected_bwt_prefix = 8_u64
        .checked_add(
            bwt_words
                .checked_mul(8)
                .ok_or(CombinedIndexBuildError::Structure(
                    "BWT prefix length overflow",
                ))?,
        )
        .ok_or(CombinedIndexBuildError::Structure(
            "BWT prefix length overflow",
        ))?;
    if bwt_file.metadata()?.len() != expected_bwt_prefix {
        return Err(CombinedIndexBuildError::Structure(
            "BWT prefix length changed before lookup append",
        ));
    }
    let lookup_header_bytes = size_of::<u64>() as u64;
    let high_bytes = LOOKUP_ENTRIES.checked_mul(size_of::<u32>() as u64).ok_or(
        CombinedIndexBuildError::Structure("lookup high byte count overflow"),
    )?;
    let low_bytes = LOOKUP_ENTRIES;
    let high_offset = expected_bwt_prefix.checked_add(lookup_header_bytes).ok_or(
        CombinedIndexBuildError::Structure("lookup high offset overflow"),
    )?;
    let low_offset =
        high_offset
            .checked_add(high_bytes)
            .ok_or(CombinedIndexBuildError::Structure(
                "lookup low offset overflow",
            ))?;
    let final_bytes =
        low_offset
            .checked_add(low_bytes)
            .ok_or(CombinedIndexBuildError::Structure(
                "BWT file length overflow",
            ))?;

    let mut header = bsbit_io::reopen_read_write(bwt_file)?;
    header.seek(SeekFrom::Start(expected_bwt_prefix))?;
    write_u64(&mut header, LOOKUP_ENTRIES)?;
    header.set_len(final_bytes)?;

    // Serialize both packed lookup planes directly into their final ranges.
    // The independent cursors preserve the single traversal over 3^16 keys.
    let mut high_file = bsbit_io::reopen_read_write(bwt_file)?;
    high_file.seek(SeekFrom::Start(high_offset))?;
    let mut high = BufWriter::with_capacity(IO_BUFFER_BYTES, high_file);
    let mut low_file = bsbit_io::reopen_read_write(bwt_file)?;
    low_file.seek(SeekFrom::Start(low_offset))?;
    let mut low = BufWriter::with_capacity(IO_BUFFER_BYTES, low_file);
    let (final_boundary, total_windows) = visit_lookup_intervals(
        tasks,
        short_suffix_thresholds,
        |_key, _lower, _upper, boundary, gap| {
            write_lookup_boundary(&mut high, &mut low, boundary, gap)?;
            Ok(())
        },
    )?;
    write_lookup_boundary(&mut high, &mut low, final_boundary, 0)?;
    let expected_windows = text_length.saturating_sub((LOOKUP_BASES - 1) as u64);
    if total_windows != expected_windows {
        return Err(CombinedIndexBuildError::Structure(
            "dense lookup occurrence total differs from projected windows",
        ));
    }
    high.flush()?;
    low.flush()?;
    drop(high);
    drop(low);
    header.sync_all()?;
    Ok(())
}

fn write_lookup_boundary(
    high: &mut impl Write,
    low: &mut impl Write,
    boundary: u64,
    gap: u8,
) -> Result<(), CombinedIndexBuildError> {
    let boundary_high = boundary >> 8;
    if boundary_high > LOOKUP_BOUNDARY_HIGH_MASK || gap >= 1 << LOOKUP_GAP_BITS {
        return Err(CombinedIndexBuildError::Structure(
            "lookup boundary exceeds the packed 36-bit domain",
        ));
    }
    let packed_high = u32::try_from(boundary_high).expect("28-bit boundary high fits u32")
        | (u32::from(gap) << 28);
    write_u32(high, packed_high)?;
    low.write_all(&[u8::try_from(boundary & 0xff).expect("low byte fits u8")])?;
    Ok(())
}

/// Projects a consumed catalog, releases it before libsais, and publishes a
/// complete compatible SA16 image without overwriting an existing component.
///
/// # Errors
///
/// Returns the first projection, allocation, libsais, serialization,
/// validation, or create-only publication failure.
pub fn build_combined_index_from_catalog_create_new(
    contigs: Vec<ContigInput>,
    reference_semantic_digest: ReferenceSemanticDigest,
    prefix: &Path,
    options: CombinedIndexBuildOptions,
) -> Result<(), CombinedIndexBuildError> {
    let observed_reference_digest = compute_reference_semantic_digest(&contigs);
    if observed_reference_digest != reference_semantic_digest {
        return Err(CombinedIndexBuildError::ReferenceDigestMismatch {
            expected: reference_semantic_digest,
            observed: observed_reference_digest,
        });
    }
    let [
        byte_0,
        byte_1,
        byte_2,
        byte_3,
        byte_4,
        byte_5,
        byte_6,
        byte_7,
        ..,
    ] = observed_reference_digest.into_bytes();
    let projection_salt = u64::from_le_bytes([
        byte_0, byte_1, byte_2, byte_3, byte_4, byte_5, byte_6, byte_7,
    ]);
    let projected = project_combined_packed_text(&contigs, projection_salt, options.threads)?;
    let text_length = u64::try_from(projected.len())
        .map_err(|_| CombinedIndexBuildError::Argument("projected text length exceeds u64"))?;
    validate_combined_text_length(text_length, DEFAULT_COMBINED_INDEX_SA_STRIDE)?;
    if projected.reference_bases().checked_mul(2) != Some(text_length) {
        return Err(CombinedIndexBuildError::Argument(
            "projected text length differs from twice the reference",
        ));
    }
    drop(contigs);
    build_combined_index_bounded_create_new(projected, reference_semantic_digest, prefix, options)
}

fn compute_reference_semantic_digest(contigs: &[ContigInput]) -> ReferenceSemanticDigest {
    let contig_count =
        u64::try_from(contigs.len()).expect("supported pointer width fits the semantic u64 domain");
    let mut semantic = ReferenceSemanticDigestBuilder::new(contig_count);
    for contig in contigs {
        semantic
            .push_normalized_contig(contig.name(), contig.sequence().bases())
            .expect("validated catalog dimensions fit the semantic digest contract");
    }
    semantic
        .finish()
        .expect("declared and observed catalog counts are equal")
}

#[allow(clippy::too_many_lines)]
fn build_combined_index_bounded_create_new(
    projected: PackedProjectedText,
    reference_semantic_digest: ReferenceSemanticDigest,
    prefix: &Path,
    options: CombinedIndexBuildOptions,
) -> Result<(), CombinedIndexBuildError> {
    let text_length = u64::try_from(projected.len())
        .map_err(|_| CombinedIndexBuildError::Argument("projected text length exceeds u64"))?;
    validate_combined_text_length(text_length, DEFAULT_COMBINED_INDEX_SA_STRIDE)?;
    if projected.reference_bases().checked_mul(2) != Some(text_length) {
        return Err(CombinedIndexBuildError::Argument(
            "projected text length differs from twice the reference",
        ));
    }
    let tail_length = projected.len().min(LOOKUP_BASES - 1);
    let tail = (projected.len() - tail_length..projected.len())
        .map(|position| projected.get(position))
        .collect::<Vec<_>>();
    let short_thresholds = short_suffix_thresholds(&tail)?;
    let staging = StagedCombinedIndex::create(prefix)?;

    let config = BoundedBwtConfig::new(options.memory_mib, options.threads)?.with_sample_stride(
        usize::try_from(DEFAULT_COMBINED_INDEX_SA_STRIDE).expect("validated SA stride fits usize"),
    )?;
    let state = build_bounded_bwt(projected, config)?;

    let dimensions = write_bounded_bwt_and_occ(&state, &staging.bwt.file, &staging.occ.file)?;

    write_bounded_sa16(&state, DEFAULT_COMBINED_INDEX_SA_STRIDE, &staging.sa.file)?;
    drop(state);

    let rank = BuildRank::open(&staging.bwt.file, &staging.occ.file, dimensions)?;
    let lookup_counts = build_lookup_counts(&rank, options.threads)?;
    drop(rank);
    lookup_counts.write(
        &short_thresholds,
        text_length,
        dimensions.bwt_words,
        &staging.bwt.file,
    )?;

    write_metadata(
        &staging.meta.file,
        dimensions,
        DEFAULT_COMBINED_INDEX_SA_STRIDE,
        reference_semantic_digest,
    )?;
    let completed = staging.seal()?;
    let index = CombinedIndex::open(&completed.stage.meta)?;
    if index.suffix_count() != dimensions.suffix_count {
        return Err(CombinedIndexBuildError::Structure(
            "runtime reader returned a foreign suffix count",
        ));
    }
    index.verify_reference_semantic_digest(reference_semantic_digest)?;
    drop(index);
    if !completed.staging_identities_match()? {
        return Err(CombinedIndexBuildError::Structure(
            "staged combined-index component identity changed during validation",
        ));
    }
    completed.publish()?;
    Ok(())
}

fn write_metadata(
    file: &File,
    dimensions: BwtDimensions,
    sa_stride: u32,
    reference_semantic_digest: ReferenceSemanticDigest,
) -> Result<(), CombinedIndexBuildError> {
    let mut bytes = [0_u8; META_BYTES];
    put_u64(&mut bytes, 0, dimensions.suffix_count);
    put_u64(&mut bytes, 8, dimensions.sentinel_row);
    for (index, &value) in dimensions.first_occurrence.iter().enumerate() {
        put_u64(&mut bytes, 16 + index * 8, value);
    }
    put_u64(&mut bytes, 48, dimensions.suffix_count);
    put_u32(&mut bytes, 56, sa_stride);
    put_u32(&mut bytes, 60, OCC_STRIDE);
    put_u32(&mut bytes, 64, HIGH_OCC_STRIDE);
    bytes[META_EXTENSION_OFFSET..META_EXTENSION_OFFSET + META_EXTENSION_MAGIC.len()]
        .copy_from_slice(META_EXTENSION_MAGIC);
    put_u16(&mut bytes, 76, META_EXTENSION_MAJOR);
    put_u16(&mut bytes, 78, META_EXTENSION_MINOR);
    put_u32(&mut bytes, 80, META_BYTES_U32);
    bytes[META_DIGEST_OFFSET..META_DIGEST_OFFSET + 32]
        .copy_from_slice(reference_semantic_digest.as_bytes());
    let mut file = reset_component_file(file)?;
    file.write_all(&bytes)?;
    file.sync_all()?;
    Ok(())
}

fn put_u16(bytes: &mut [u8], offset: usize, value: u16) {
    bytes[offset..offset + 2].copy_from_slice(&value.to_le_bytes());
}

fn put_u32(bytes: &mut [u8], offset: usize, value: u32) {
    bytes[offset..offset + 4].copy_from_slice(&value.to_le_bytes());
}

fn put_u64(bytes: &mut [u8], offset: usize, value: u64) {
    bytes[offset..offset + 8].copy_from_slice(&value.to_le_bytes());
}

#[cfg(test)]
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct SaAuditSample {
    row: u64,
    coordinate: u64,
}

#[cfg(test)]
#[derive(Debug, Eq, PartialEq)]
enum SaAudit {
    Disabled,
    Sampled(Vec<SaAuditSample>),
    Exhaustive(Vec<u64>),
}

#[cfg(test)]
impl SaAudit {
    fn prepare(
        packed_samples: Vec<u64>,
        requested: u64,
        sa_stride: u32,
    ) -> Result<Self, CombinedIndexBuildError> {
        if requested == 0 {
            return Ok(Self::Disabled);
        }
        if requested == u64::MAX {
            return Ok(Self::Exhaustive(packed_samples));
        }
        let packed_len = u64::try_from(packed_samples.len()).expect("sample count fits u64");
        let checks = requested.min(packed_len);
        let stride = packed_len.div_ceil(checks.max(1));
        let capacity = usize::try_from(checks)
            .map_err(|_| CombinedIndexBuildError::Allocation("sampled SA16 audit expectations"))?;
        let mut expected = Vec::new();
        expected
            .try_reserve_exact(capacity)
            .map_err(|_| CombinedIndexBuildError::Allocation("sampled SA16 audit expectations"))?;
        for ordinal in 0..checks {
            let slot = (ordinal * stride).min(packed_len - 1);
            let packed = packed_samples[usize::try_from(slot).expect("sample slot fits usize")];
            expected.push(SaAuditSample {
                row: packed >> SA_VALUE_BITS,
                coordinate: (packed & SA_VALUE_MASK) * u64::from(sa_stride),
            });
        }
        Ok(Self::Sampled(expected))
    }
}

#[cfg(test)]
#[path = "../../tests/whitebox/build_combined.rs"]
mod tests;
