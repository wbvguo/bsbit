//! Compact reference catalog used by the indexed alignment path.
//!
//! `bsbit index` writes this file as its opaque public handle. Bases use a
//! three-bit representation, while per-contig checksums allow alignment to
//! validate and decode independent contigs in parallel. Search data remains in
//! private sibling files and is bound to this catalog by the semantic digest.

use core::fmt;
use core::ptr::NonNull;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::os::fd::AsRawFd;
use std::path::Path;
use std::thread;

use bsbit_core::alphabet::Base;
use bsbit_core::reference::{
    ReferenceSemanticDigest, ReferenceSemanticDigestBuildError, ReferenceSemanticDigestBuilder,
};
use bsbit_core::sequence::NormalizedSequence;
use bsbit_io::{PublicationError, PublishedFile, StagedFile};
use sha2::{Digest, Sha256};

use crate::reference::{
    ContigInput, ReferenceBuildError, ReferenceCatalogLimits, validate_reference_catalog,
};

const MAGIC: &[u8; 8] = b"BSBITCAT";
const FORMAT_MAJOR: u16 = 1;
const FORMAT_MINOR: u16 = 0;
const ENDIAN_MARKER: u32 = 0x0102_0304;
const HEADER_LEN: usize = 160;
const HEADER_LEN_U32: u32 = 160;
const ENTRY_LEN: usize = 80;
const ENTRY_LEN_U32: u32 = 80;
const ALIGNMENT_BYTES: usize = 64;
const ALIGNMENT: u64 = 64;
const HEADER_DIGEST_OFFSET: usize = 112;
const HEADER_DIGEST_END: usize = 144;
const HEADER_DOMAIN: &[u8] = b"BSBIT-REFERENCE-CATALOG-HEADER-SHA256-V1\0";
const CONTIG_DOMAIN: &[u8] = b"BSBIT-REFERENCE-CATALOG-CONTIG-SHA256-V1\0";

/// Aggregate dimensions and semantic identity of one compact catalog.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ReferenceCatalogSummary {
    file_length: u64,
    contig_count: u64,
    total_name_bytes: u64,
    total_reference_bases: u64,
    semantic_digest: ReferenceSemanticDigest,
}

impl ReferenceCatalogSummary {
    /// Returns the complete encoded byte length.
    #[must_use]
    pub const fn file_length(self) -> u64 {
        self.file_length
    }

    /// Returns the ordered contig count.
    #[must_use]
    pub const fn contig_count(self) -> u64 {
        self.contig_count
    }

    /// Returns the aggregate exact contig-name bytes.
    #[must_use]
    pub const fn total_name_bytes(self) -> u64 {
        self.total_name_bytes
    }

    /// Returns the aggregate normalized reference bases.
    #[must_use]
    pub const fn total_reference_bases(self) -> u64 {
        self.total_reference_bases
    }

    /// Returns the semantic reference digest shared with the search index.
    #[must_use]
    pub const fn semantic_digest(self) -> ReferenceSemanticDigest {
        self.semantic_digest
    }
}

/// A compact-catalog encoding, validation, or publication failure.
#[derive(Debug)]
pub enum ReferenceCatalogError {
    /// A file operation failed.
    Io(io::Error),
    /// Catalog semantics or dimensions were invalid.
    Catalog(ReferenceBuildError),
    /// Semantic digest construction failed.
    Semantic(ReferenceSemanticDigestBuildError),
    /// Staging, publication, or replacement failed.
    Publication(PublicationError),
    /// A fixed field, offset, or packed-base invariant was invalid.
    Structure(&'static str),
    /// One contig payload checksum disagreed.
    ContigIntegrity {
        /// Zero-based catalog ordinal whose payload failed validation.
        ordinal: u64,
    },
    /// A caller or search index required a different reference.
    ReferenceDigestMismatch {
        /// Digest required by the caller or bound search index.
        expected: ReferenceSemanticDigest,
        /// Digest stored by the catalog.
        observed: ReferenceSemanticDigest,
    },
    /// A parallel decoding worker panicked.
    WorkerPanic,
}

impl fmt::Display for ReferenceCatalogError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(source) => write!(formatter, "reference catalog I/O: {source}"),
            Self::Catalog(source) => write!(formatter, "reference catalog is invalid: {source}"),
            Self::Semantic(source) => {
                write!(
                    formatter,
                    "reference catalog semantic digest failed: {source}"
                )
            }
            Self::Publication(source) => {
                write!(formatter, "reference catalog publication failed: {source}")
            }
            Self::Structure(message) => {
                write!(formatter, "reference catalog structure: {message}")
            }
            Self::ContigIntegrity { ordinal } => {
                write!(
                    formatter,
                    "reference catalog contig {ordinal} checksum differs"
                )
            }
            Self::ReferenceDigestMismatch { expected, observed } => write!(
                formatter,
                "reference catalog digest differs: expected {expected}, observed {observed}"
            ),
            Self::WorkerPanic => formatter.write_str("reference catalog decode worker panicked"),
        }
    }
}

impl std::error::Error for ReferenceCatalogError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Io(source) => Some(source),
            Self::Catalog(source) => Some(source),
            Self::Semantic(source) => Some(source),
            Self::Publication(source) => Some(source),
            Self::Structure(_)
            | Self::ContigIntegrity { .. }
            | Self::ReferenceDigestMismatch { .. }
            | Self::WorkerPanic => None,
        }
    }
}

impl From<io::Error> for ReferenceCatalogError {
    fn from(source: io::Error) -> Self {
        Self::Io(source)
    }
}

impl From<PublicationError> for ReferenceCatalogError {
    fn from(source: PublicationError) -> Self {
        Self::Publication(source)
    }
}

/// Successful reference-catalog publication details.
#[derive(Debug)]
pub struct ReferenceCatalogPublication {
    summary: ReferenceCatalogSummary,
    published: PublishedFile,
}

impl ReferenceCatalogPublication {
    /// Returns the encoded catalog summary.
    #[must_use]
    pub const fn summary(&self) -> ReferenceCatalogSummary {
        self.summary
    }

    /// Returns the private staging path used by the publication.
    #[must_use]
    pub fn staging_path(&self) -> &Path {
        self.published.staging_path()
    }

    /// Returns a post-publication cleanup warning.
    #[must_use]
    pub const fn cleanup_error(&self) -> Option<io::ErrorKind> {
        self.published.cleanup_warning()
    }

    /// Retracts the published catalog while it still names this file.
    ///
    /// # Errors
    ///
    /// Returns an identity-safe publication rollback failure.
    pub fn rollback(self) -> Result<(), PublicationError> {
        self.published.rollback()
    }
}

#[derive(Clone, Copy, Debug)]
struct CatalogEntry {
    name_offset: u64,
    name_len: u64,
    packed_offset: u64,
    base_count: u64,
    packed_len: u64,
    digest: [u8; 32],
}

/// Validates and writes one compact reference catalog.
///
/// # Errors
///
/// Returns semantic, checked-arithmetic, or output failures.
pub fn write_reference_catalog<W: Write>(
    contigs: &[ContigInput],
    writer: &mut W,
) -> Result<ReferenceCatalogSummary, ReferenceCatalogError> {
    let layout = build_catalog_layout(contigs)?;
    let header = encode_header(
        layout.summary,
        layout.names_offset,
        layout.packed_section_offset,
    );
    writer.write_all(&header)?;
    for entry in &layout.entries {
        writer.write_all(&encode_entry(*entry))?;
    }
    let mut written = layout.names_offset;
    for contig in contigs {
        writer.write_all(contig.name())?;
        written = written
            .checked_add(u64::try_from(contig.name().len()).unwrap_or(u64::MAX))
            .ok_or(ReferenceCatalogError::Structure(
                "written name offset overflow",
            ))?;
    }
    write_zero_padding(writer, &mut written, layout.packed_section_offset)?;
    for (contig, entry) in contigs.iter().zip(&layout.entries) {
        write_zero_padding(writer, &mut written, entry.packed_offset)?;
        for_each_packed_chunk(contig.sequence().bases(), |chunk| writer.write_all(chunk))?;
        written = written
            .checked_add(entry.packed_len)
            .ok_or(ReferenceCatalogError::Structure(
                "written packed offset overflow",
            ))?;
    }
    if written != layout.summary.file_length {
        return Err(ReferenceCatalogError::Structure(
            "encoded byte count differs from declared file length",
        ));
    }
    Ok(layout.summary)
}

struct CatalogLayout {
    summary: ReferenceCatalogSummary,
    entries: Vec<CatalogEntry>,
    names_offset: u64,
    packed_section_offset: u64,
}

fn build_catalog_layout(contigs: &[ContigInput]) -> Result<CatalogLayout, ReferenceCatalogError> {
    let metrics = validate_reference_catalog(contigs, ReferenceCatalogLimits::MAX)
        .map_err(ReferenceCatalogError::Catalog)?;
    let contig_count = metrics.contig_count();
    let total_name_bytes = metrics.total_name_bytes();
    let total_reference_bases = metrics.total_reference_bases();
    let directory_bytes =
        contig_count
            .checked_mul(ENTRY_LEN as u64)
            .ok_or(ReferenceCatalogError::Structure(
                "directory length overflow",
            ))?;
    let names_offset = (HEADER_LEN as u64)
        .checked_add(directory_bytes)
        .ok_or(ReferenceCatalogError::Structure("name offset overflow"))?;
    let packed_section_offset = align_up(
        names_offset
            .checked_add(total_name_bytes)
            .ok_or(ReferenceCatalogError::Structure("packed offset overflow"))?,
        ALIGNMENT,
    )
    .ok_or(ReferenceCatalogError::Structure(
        "packed alignment overflow",
    ))?;

    let mut semantic = ReferenceSemanticDigestBuilder::new(contig_count);
    let mut entries = Vec::new();
    entries
        .try_reserve_exact(contigs.len())
        .map_err(|_| ReferenceCatalogError::Structure("directory allocation failed"))?;
    let mut next_name = names_offset;
    let mut next_packed = packed_section_offset;
    for contig in contigs {
        semantic
            .push_normalized_contig(contig.name(), contig.sequence().bases())
            .map_err(ReferenceCatalogError::Semantic)?;
        let name_len = u64::try_from(contig.name().len())
            .map_err(|_| ReferenceCatalogError::Structure("contig name length exceeds u64"))?;
        let base_count = contig.sequence().len();
        let packed_len = packed_length(base_count).ok_or(ReferenceCatalogError::Structure(
            "packed base length overflow",
        ))?;
        next_packed = align_up(next_packed, ALIGNMENT).ok_or(ReferenceCatalogError::Structure(
            "packed contig offset overflow",
        ))?;
        let digest = contig_digest(contig.name(), contig.sequence().bases());
        entries.push(CatalogEntry {
            name_offset: next_name,
            name_len,
            packed_offset: next_packed,
            base_count,
            packed_len,
            digest,
        });
        next_name = next_name
            .checked_add(name_len)
            .ok_or(ReferenceCatalogError::Structure(
                "contig name offset overflow",
            ))?;
        next_packed =
            next_packed
                .checked_add(packed_len)
                .ok_or(ReferenceCatalogError::Structure(
                    "packed contig end overflow",
                ))?;
    }
    let semantic_digest = semantic.finish().map_err(ReferenceCatalogError::Semantic)?;
    Ok(CatalogLayout {
        summary: ReferenceCatalogSummary {
            file_length: next_packed,
            contig_count,
            total_name_bytes,
            total_reference_bases,
            semantic_digest,
        },
        entries,
        names_offset,
        packed_section_offset,
    })
}

/// Publishes one compact catalog without replacing an existing target.
///
/// # Errors
///
/// Returns encoding, staging, synchronization, or create-only publication
/// failures.
pub fn publish_reference_catalog_create_new(
    contigs: &[ContigInput],
    target: &Path,
    staging: &Path,
) -> Result<ReferenceCatalogPublication, ReferenceCatalogError> {
    publish_reference_catalog(contigs, target, staging, false)
}

/// Publishes one compact catalog, atomically replacing an existing target.
///
/// # Errors
///
/// Returns encoding, staging, synchronization, backup, or replacement
/// failures.
pub fn publish_reference_catalog_replace(
    contigs: &[ContigInput],
    target: &Path,
    staging: &Path,
) -> Result<ReferenceCatalogPublication, ReferenceCatalogError> {
    publish_reference_catalog(contigs, target, staging, true)
}

fn publish_reference_catalog(
    contigs: &[ContigInput],
    target: &Path,
    staging: &Path,
    replace: bool,
) -> Result<ReferenceCatalogPublication, ReferenceCatalogError> {
    let mut staged = StagedFile::create_new(staging)?;
    let file = staged.take_file()?;
    let mut writer = BufWriter::new(file);
    let summary = write_reference_catalog(contigs, &mut writer)?;
    writer.flush()?;
    let file = writer
        .into_inner()
        .map_err(|error| ReferenceCatalogError::Io(error.into_error()))?;
    let completed = staged.complete(file)?;
    let published = if replace {
        completed.publish_replace_at(target)?
    } else {
        completed.publish_create_new_at(target)?
    };
    Ok(ReferenceCatalogPublication { summary, published })
}

/// A fully validated and decoded compact catalog.
#[derive(Debug)]
pub(crate) struct LoadedReferenceCatalog {
    contigs: Vec<ContigInput>,
    summary: ReferenceCatalogSummary,
}

impl LoadedReferenceCatalog {
    pub(crate) fn into_contigs(self) -> Vec<ContigInput> {
        self.contigs
    }

    pub(crate) const fn summary(&self) -> ReferenceCatalogSummary {
        self.summary
    }
}

/// Opens, validates, and decodes a compact catalog with bounded parallelism.
pub(crate) fn load_reference_catalog(
    path: &Path,
    expected_digest: Option<ReferenceSemanticDigest>,
    threads: usize,
) -> Result<LoadedReferenceCatalog, ReferenceCatalogError> {
    let file = File::open(path)?;
    let mapping = ReadOnlyCatalogMapping::map(&file)?;
    let bytes = mapping.as_slice();
    let (summary, entries) = decode_header_and_entries(bytes)?;
    if let Some(expected) = expected_digest
        && expected != summary.semantic_digest
    {
        return Err(ReferenceCatalogError::ReferenceDigestMismatch {
            expected,
            observed: summary.semantic_digest,
        });
    }
    let worker_count = threads.max(1).min(entries.len().max(1));
    let ranges = balanced_entry_ranges(&entries, worker_count);
    let mut contigs = Vec::new();
    contigs
        .try_reserve_exact(entries.len())
        .map_err(|_| ReferenceCatalogError::Structure("decoded catalog allocation failed"))?;
    thread::scope(|scope| -> Result<(), ReferenceCatalogError> {
        let mut workers = Vec::new();
        for range in ranges {
            let first_ordinal = range.start;
            let chunk = &entries[range];
            workers.push(scope.spawn(move || decode_contig_chunk(bytes, chunk, first_ordinal)));
        }
        for worker in workers {
            let decoded = worker
                .join()
                .map_err(|_| ReferenceCatalogError::WorkerPanic)??;
            contigs.extend(decoded);
        }
        Ok(())
    })?;
    let decoded_metrics = validate_reference_catalog(&contigs, ReferenceCatalogLimits::MAX)
        .map_err(ReferenceCatalogError::Catalog)?;
    if decoded_metrics.contig_count() != summary.contig_count
        || decoded_metrics.total_name_bytes() != summary.total_name_bytes
        || decoded_metrics.total_reference_bases() != summary.total_reference_bases
    {
        return Err(ReferenceCatalogError::Structure(
            "decoded catalog dimensions differ from the header",
        ));
    }
    Ok(LoadedReferenceCatalog { contigs, summary })
}

fn balanced_entry_ranges(entries: &[CatalogEntry], workers: usize) -> Vec<core::ops::Range<usize>> {
    if entries.is_empty() {
        return Vec::new();
    }
    let workers = workers.max(1).min(entries.len());
    let mut ranges = Vec::with_capacity(workers);
    let mut start = 0_usize;
    let mut remaining_bytes = entries
        .iter()
        .fold(0_u64, |total, entry| total.saturating_add(entry.packed_len));
    for worker in 0..workers {
        let remaining_workers = workers - worker;
        let remaining_entries = entries.len() - start;
        if remaining_workers == 1 {
            ranges.push(start..entries.len());
            break;
        }
        let target = remaining_bytes.div_ceil(remaining_workers as u64);
        let maximum_end = entries.len() - (remaining_workers - 1);
        let mut end = start;
        let mut assigned = 0_u64;
        while end < maximum_end && (end == start || assigned < target) {
            assigned = assigned.saturating_add(entries[end].packed_len);
            end += 1;
        }
        ranges.push(start..end);
        remaining_bytes = remaining_bytes.saturating_sub(assigned);
        start = end;
        debug_assert!(remaining_entries >= remaining_workers);
    }
    ranges
}

fn decode_contig_chunk(
    bytes: &[u8],
    entries: &[CatalogEntry],
    first_ordinal: usize,
) -> Result<Vec<ContigInput>, ReferenceCatalogError> {
    let mut output = Vec::new();
    output
        .try_reserve_exact(entries.len())
        .map_err(|_| ReferenceCatalogError::Structure("decoded contig allocation failed"))?;
    for (within, entry) in entries.iter().copied().enumerate() {
        let ordinal = first_ordinal + within;
        let name = checked_slice(bytes, entry.name_offset, entry.name_len)?;
        let packed = checked_slice(bytes, entry.packed_offset, entry.packed_len)?;
        if digest_contig_payload(name, entry.base_count, packed) != entry.digest {
            return Err(ReferenceCatalogError::ContigIntegrity {
                ordinal: u64::try_from(ordinal).unwrap_or(u64::MAX),
            });
        }
        validate_packed_tail(packed, entry.base_count)?;
        let base_count = usize::try_from(entry.base_count)
            .map_err(|_| ReferenceCatalogError::Structure("contig bases exceed usize"))?;
        let mut decoded = Vec::new();
        decoded
            .try_reserve_exact(base_count)
            .map_err(|_| ReferenceCatalogError::Structure("decoded bases allocation failed"))?;
        for position in 0..base_count {
            let bit = position
                .checked_mul(3)
                .ok_or(ReferenceCatalogError::Structure(
                    "packed bit offset overflow",
                ))?;
            let byte = bit / 8;
            let shift = bit % 8;
            let low = u16::from(packed[byte]);
            let high = packed.get(byte + 1).map_or(0, |value| u16::from(*value));
            let code = u8::try_from(((low | (high << 8)) >> shift) & 0b111)
                .expect("three packed bits fit u8");
            decoded.push(match code {
                0 => Base::A,
                1 => Base::C,
                2 => Base::G,
                3 => Base::T,
                4 => Base::N,
                _ => {
                    return Err(ReferenceCatalogError::Structure(
                        "packed catalog contains an invalid base code",
                    ));
                }
            });
        }
        output.push(ContigInput::new(
            name.to_vec(),
            NormalizedSequence::from(decoded),
        ));
    }
    Ok(output)
}

fn decode_header_and_entries(
    bytes: &[u8],
) -> Result<(ReferenceCatalogSummary, Vec<CatalogEntry>), ReferenceCatalogError> {
    let header = decode_catalog_header(bytes)?;
    let entries = decode_catalog_entries(bytes, &header)?;
    Ok((header.summary, entries))
}

struct DecodedCatalogHeader {
    summary: ReferenceCatalogSummary,
    names_offset: u64,
    packed_section_offset: u64,
}

fn decode_catalog_header(bytes: &[u8]) -> Result<DecodedCatalogHeader, ReferenceCatalogError> {
    let header = bytes
        .get(..HEADER_LEN)
        .ok_or(ReferenceCatalogError::Structure(
            "file is shorter than the header",
        ))?;
    if &header[..8] != MAGIC {
        return Err(ReferenceCatalogError::Structure("magic is invalid"));
    }
    if slice_u16(header, 8) != FORMAT_MAJOR
        || slice_u16(header, 10) != FORMAT_MINOR
        || slice_u32(header, 12) != HEADER_LEN_U32
        || slice_u32(header, 16) != ENDIAN_MARKER
        || slice_u32(header, 20) != ENTRY_LEN_U32
        || header[144..160] != [0; 16]
    {
        return Err(ReferenceCatalogError::Structure(
            "version, widths, endian marker, or reserved bytes are invalid",
        ));
    }
    let expected_header_digest = header_digest(header);
    if header[HEADER_DIGEST_OFFSET..HEADER_DIGEST_END] != expected_header_digest {
        return Err(ReferenceCatalogError::Structure("header checksum differs"));
    }
    let file_length = slice_u64(header, 24);
    if usize::try_from(file_length).ok() != Some(bytes.len()) {
        return Err(ReferenceCatalogError::Structure("file length differs"));
    }
    let contig_count = slice_u64(header, 32);
    let total_name_bytes = slice_u64(header, 40);
    let total_reference_bases = slice_u64(header, 48);
    let directory_offset = slice_u64(header, 56);
    let names_offset = slice_u64(header, 64);
    let packed_section_offset = slice_u64(header, 72);
    if directory_offset != HEADER_LEN as u64 {
        return Err(ReferenceCatalogError::Structure("directory offset differs"));
    }
    let directory_bytes =
        contig_count
            .checked_mul(ENTRY_LEN as u64)
            .ok_or(ReferenceCatalogError::Structure(
                "directory length overflow",
            ))?;
    if names_offset
        != directory_offset
            .checked_add(directory_bytes)
            .ok_or(ReferenceCatalogError::Structure("name offset overflow"))?
        || packed_section_offset
            < names_offset
                .checked_add(total_name_bytes)
                .ok_or(ReferenceCatalogError::Structure("packed offset overflow"))?
        || !packed_section_offset.is_multiple_of(ALIGNMENT)
    {
        return Err(ReferenceCatalogError::Structure(
            "directory, names, or packed section is not canonical",
        ));
    }
    let semantic_digest = ReferenceSemanticDigest::from_bytes(
        header[80..112]
            .try_into()
            .expect("fixed semantic digest width"),
    );
    Ok(DecodedCatalogHeader {
        summary: ReferenceCatalogSummary {
            file_length,
            contig_count,
            total_name_bytes,
            total_reference_bases,
            semantic_digest,
        },
        names_offset,
        packed_section_offset,
    })
}

fn decode_catalog_entries(
    bytes: &[u8],
    header: &DecodedCatalogHeader,
) -> Result<Vec<CatalogEntry>, ReferenceCatalogError> {
    let entry_count = usize::try_from(header.summary.contig_count)
        .map_err(|_| ReferenceCatalogError::Structure("contig count exceeds usize"))?;
    let mut entries = Vec::new();
    entries
        .try_reserve_exact(entry_count)
        .map_err(|_| ReferenceCatalogError::Structure("directory allocation failed"))?;
    let mut next_name = header.names_offset;
    let mut next_packed = header.packed_section_offset;
    let mut observed_names = 0_u64;
    let mut observed_bases = 0_u64;
    for ordinal in 0..entry_count {
        let offset =
            HEADER_LEN
                .checked_add(ordinal.checked_mul(ENTRY_LEN).ok_or(
                    ReferenceCatalogError::Structure("directory entry offset overflow"),
                )?)
                .ok_or(ReferenceCatalogError::Structure(
                    "directory entry offset overflow",
                ))?;
        let encoded =
            bytes
                .get(offset..offset + ENTRY_LEN)
                .ok_or(ReferenceCatalogError::Structure(
                    "directory entry is truncated",
                ))?;
        let entry = decode_entry(encoded);
        let expected_packed = align_up(next_packed, ALIGNMENT).ok_or(
            ReferenceCatalogError::Structure("packed contig offset overflow"),
        )?;
        if entry.name_offset != next_name
            || entry.packed_offset != expected_packed
            || entry.packed_len
                != packed_length(entry.base_count).ok_or(ReferenceCatalogError::Structure(
                    "packed base length overflow",
                ))?
            || encoded[72..80] != [0; 8]
        {
            return Err(ReferenceCatalogError::Structure(
                "directory entry offsets, length, or reserved bytes are invalid",
            ));
        }
        checked_slice(bytes, entry.name_offset, entry.name_len)?;
        checked_slice(bytes, entry.packed_offset, entry.packed_len)?;
        next_name = next_name
            .checked_add(entry.name_len)
            .ok_or(ReferenceCatalogError::Structure("name aggregate overflow"))?;
        next_packed = entry.packed_offset.checked_add(entry.packed_len).ok_or(
            ReferenceCatalogError::Structure("packed aggregate overflow"),
        )?;
        observed_names = observed_names
            .checked_add(entry.name_len)
            .ok_or(ReferenceCatalogError::Structure("name aggregate overflow"))?;
        observed_bases = observed_bases
            .checked_add(entry.base_count)
            .ok_or(ReferenceCatalogError::Structure("base aggregate overflow"))?;
        entries.push(entry);
    }
    let expected_names_end = header
        .names_offset
        .checked_add(header.summary.total_name_bytes)
        .ok_or(ReferenceCatalogError::Structure("name aggregate overflow"))?;
    if next_name != expected_names_end
        || observed_names != header.summary.total_name_bytes
        || observed_bases != header.summary.total_reference_bases
        || next_packed != header.summary.file_length
    {
        return Err(ReferenceCatalogError::Structure(
            "catalog aggregate dimensions differ",
        ));
    }
    Ok(entries)
}

fn encode_header(
    summary: ReferenceCatalogSummary,
    names_offset: u64,
    packed_section_offset: u64,
) -> [u8; HEADER_LEN] {
    let mut header = [0_u8; HEADER_LEN];
    header[..8].copy_from_slice(MAGIC);
    header[8..10].copy_from_slice(&FORMAT_MAJOR.to_le_bytes());
    header[10..12].copy_from_slice(&FORMAT_MINOR.to_le_bytes());
    header[12..16].copy_from_slice(&HEADER_LEN_U32.to_le_bytes());
    header[16..20].copy_from_slice(&ENDIAN_MARKER.to_le_bytes());
    header[20..24].copy_from_slice(&ENTRY_LEN_U32.to_le_bytes());
    header[24..32].copy_from_slice(&summary.file_length.to_le_bytes());
    header[32..40].copy_from_slice(&summary.contig_count.to_le_bytes());
    header[40..48].copy_from_slice(&summary.total_name_bytes.to_le_bytes());
    header[48..56].copy_from_slice(&summary.total_reference_bases.to_le_bytes());
    header[56..64].copy_from_slice(&(HEADER_LEN as u64).to_le_bytes());
    header[64..72].copy_from_slice(&names_offset.to_le_bytes());
    header[72..80].copy_from_slice(&packed_section_offset.to_le_bytes());
    header[80..112].copy_from_slice(&summary.semantic_digest.into_bytes());
    let digest = header_digest(&header);
    header[HEADER_DIGEST_OFFSET..HEADER_DIGEST_END].copy_from_slice(&digest);
    header
}

fn encode_entry(entry: CatalogEntry) -> [u8; ENTRY_LEN] {
    let mut encoded = [0_u8; ENTRY_LEN];
    encoded[0..8].copy_from_slice(&entry.name_offset.to_le_bytes());
    encoded[8..16].copy_from_slice(&entry.name_len.to_le_bytes());
    encoded[16..24].copy_from_slice(&entry.packed_offset.to_le_bytes());
    encoded[24..32].copy_from_slice(&entry.base_count.to_le_bytes());
    encoded[32..40].copy_from_slice(&entry.packed_len.to_le_bytes());
    encoded[40..72].copy_from_slice(&entry.digest);
    encoded
}

fn decode_entry(encoded: &[u8]) -> CatalogEntry {
    CatalogEntry {
        name_offset: slice_u64(encoded, 0),
        name_len: slice_u64(encoded, 8),
        packed_offset: slice_u64(encoded, 16),
        base_count: slice_u64(encoded, 24),
        packed_len: slice_u64(encoded, 32),
        digest: encoded[40..72]
            .try_into()
            .expect("fixed contig digest width"),
    }
}

fn header_digest(header: &[u8]) -> [u8; 32] {
    let mut hasher = Sha256::new();
    hasher.update(HEADER_DOMAIN);
    hasher.update(&header[..HEADER_DIGEST_OFFSET]);
    hasher.update(&header[HEADER_DIGEST_END..HEADER_LEN]);
    hasher.finalize().into()
}

fn contig_digest(name: &[u8], bases: &[Base]) -> [u8; 32] {
    let mut hasher = Sha256::new();
    hasher.update(CONTIG_DOMAIN);
    hasher.update(u64::try_from(name.len()).unwrap_or(u64::MAX).to_le_bytes());
    hasher.update(u64::try_from(bases.len()).unwrap_or(u64::MAX).to_le_bytes());
    hasher.update(name);
    for_each_packed_chunk(bases, |chunk| {
        hasher.update(chunk);
        Ok(())
    })
    .expect("hashing a memory buffer cannot fail");
    hasher.finalize().into()
}

fn digest_contig_payload(name: &[u8], base_count: u64, packed: &[u8]) -> [u8; 32] {
    let mut hasher = Sha256::new();
    hasher.update(CONTIG_DOMAIN);
    hasher.update(u64::try_from(name.len()).unwrap_or(u64::MAX).to_le_bytes());
    hasher.update(base_count.to_le_bytes());
    hasher.update(name);
    hasher.update(packed);
    hasher.finalize().into()
}

fn for_each_packed_chunk(
    bases: &[Base],
    mut consume: impl FnMut(&[u8]) -> io::Result<()>,
) -> io::Result<()> {
    let mut packed = [0_u8; 3 * 1024];
    for chunk in bases.chunks(8 * 1024) {
        let mut output = 0_usize;
        for group in chunk.chunks(8) {
            let mut word = 0_u32;
            for (within, base) in group.iter().copied().enumerate() {
                word |= u32::from(base.storage_code()) << (within * 3);
            }
            let bytes = (group.len() * 3).div_ceil(8);
            packed[output..output + bytes].copy_from_slice(&word.to_le_bytes()[..bytes]);
            output += bytes;
        }
        consume(&packed[..output])?;
    }
    Ok(())
}

fn write_zero_padding<W: Write>(
    writer: &mut W,
    written: &mut u64,
    target: u64,
) -> Result<(), ReferenceCatalogError> {
    if *written > target {
        return Err(ReferenceCatalogError::Structure("encoded sections overlap"));
    }
    let zeros = [0_u8; ALIGNMENT_BYTES];
    while *written < target {
        let remaining = usize::try_from(target - *written).unwrap_or(usize::MAX);
        let chunk = remaining.min(zeros.len());
        writer.write_all(&zeros[..chunk])?;
        *written += u64::try_from(chunk).expect("padding chunk fits u64");
    }
    Ok(())
}

fn validate_packed_tail(packed: &[u8], base_count: u64) -> Result<(), ReferenceCatalogError> {
    let used_bits = base_count
        .checked_mul(3)
        .ok_or(ReferenceCatalogError::Structure(
            "packed bit count overflow",
        ))?;
    let remainder = (used_bits % 8) as u8;
    if remainder != 0 {
        let mask = u8::MAX << remainder;
        if packed.last().is_none_or(|value| value & mask != 0) {
            return Err(ReferenceCatalogError::Structure(
                "packed catalog has nonzero trailing bits",
            ));
        }
    }
    Ok(())
}

fn packed_length(base_count: u64) -> Option<u64> {
    base_count
        .checked_mul(3)?
        .checked_add(7)
        .map(|bits| bits / 8)
}

fn align_up(value: u64, alignment: u64) -> Option<u64> {
    value
        .checked_add(alignment.checked_sub(1)?)
        .map(|value| value / alignment * alignment)
}

fn checked_slice(bytes: &[u8], offset: u64, length: u64) -> Result<&[u8], ReferenceCatalogError> {
    let start = usize::try_from(offset)
        .map_err(|_| ReferenceCatalogError::Structure("catalog offset exceeds usize"))?;
    let length = usize::try_from(length)
        .map_err(|_| ReferenceCatalogError::Structure("catalog length exceeds usize"))?;
    let end = start
        .checked_add(length)
        .ok_or(ReferenceCatalogError::Structure("catalog range overflow"))?;
    bytes
        .get(start..end)
        .ok_or(ReferenceCatalogError::Structure(
            "catalog range exceeds file",
        ))
}

fn slice_u16(bytes: &[u8], offset: usize) -> u16 {
    u16::from_le_bytes(bytes[offset..offset + 2].try_into().expect("two bytes"))
}

fn slice_u32(bytes: &[u8], offset: usize) -> u32 {
    u32::from_le_bytes(bytes[offset..offset + 4].try_into().expect("four bytes"))
}

fn slice_u64(bytes: &[u8], offset: usize) -> u64 {
    u64::from_le_bytes(bytes[offset..offset + 8].try_into().expect("eight bytes"))
}

#[derive(Debug)]
struct ReadOnlyCatalogMapping {
    pointer: NonNull<u8>,
    length: usize,
}

// SAFETY: the mapping is immutable for its complete lifetime.
unsafe impl Send for ReadOnlyCatalogMapping {}
// SAFETY: every shared read is bounded against the immutable mapping length.
unsafe impl Sync for ReadOnlyCatalogMapping {}

impl ReadOnlyCatalogMapping {
    fn map(file: &File) -> Result<Self, ReferenceCatalogError> {
        let length = usize::try_from(file.metadata()?.len())
            .map_err(|_| ReferenceCatalogError::Structure("file length exceeds usize"))?;
        if length == 0 {
            return Err(ReferenceCatalogError::Structure("catalog file is empty"));
        }
        // SAFETY: the descriptor remains live for this call, its exact length
        // came from metadata, and MAP_PRIVATE requests no writable alias.
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
        let pointer = NonNull::new(mapped.cast::<u8>())
            .ok_or(ReferenceCatalogError::Structure("mmap returned null"))?;
        Ok(Self { pointer, length })
    }

    fn as_slice(&self) -> &[u8] {
        // SAFETY: the pointer and nonzero length come from the live immutable
        // mapping, which outlives the returned shared slice.
        unsafe { core::slice::from_raw_parts(self.pointer.as_ptr(), self.length) }
    }
}

impl Drop for ReadOnlyCatalogMapping {
    fn drop(&mut self) {
        // SAFETY: this releases exactly the pointer/length returned by mmap.
        let _ = unsafe { libc::munmap(self.pointer.as_ptr().cast(), self.length) };
    }
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use bsbit_core::sequence::normalize_dna;

    use super::*;

    fn fixture() -> Vec<ContigInput> {
        vec![
            ContigInput::new(
                b"chr1".to_vec(),
                normalize_dna(b"ACGTNACGT").expect("fixture normalizes"),
            ),
            ContigInput::new(
                b"chr2".to_vec(),
                normalize_dna(b"NNNTGCA").expect("fixture normalizes"),
            ),
        ]
    }

    #[test]
    fn compact_encoding_is_three_bits_per_base_plus_metadata() {
        let contigs = fixture();
        let mut encoded = Vec::new();
        let summary = write_reference_catalog(&contigs, &mut encoded).expect("catalog writes");
        assert_eq!(&encoded[..8], MAGIC);
        assert_eq!(summary.file_length(), encoded.len() as u64);
        let packed_bytes: u64 = contigs
            .iter()
            .map(|contig| packed_length(contig.sequence().len()).expect("length"))
            .sum();
        assert!(packed_bytes <= summary.total_reference_bases().div_ceil(2));
    }

    #[test]
    fn header_and_entry_layout_round_trip() {
        let contigs = fixture();
        let mut encoded = Cursor::new(Vec::new());
        let written = write_reference_catalog(&contigs, &mut encoded).expect("catalog writes");
        let (summary, entries) =
            decode_header_and_entries(encoded.get_ref()).expect("catalog validates");
        assert_eq!(summary, written);
        let decoded = decode_contig_chunk(encoded.get_ref(), &entries, 0).expect("catalog decodes");
        assert_eq!(decoded.len(), contigs.len());
        for (actual, expected) in decoded.iter().zip(&contigs) {
            assert_eq!(actual.name(), expected.name());
            assert_eq!(actual.sequence(), expected.sequence());
        }
    }

    #[test]
    fn contig_corruption_is_rejected() {
        let contigs = fixture();
        let mut encoded = Vec::new();
        write_reference_catalog(&contigs, &mut encoded).expect("catalog writes");
        let (_, entries) = decode_header_and_entries(&encoded).expect("directory validates");
        let offset = usize::try_from(entries[0].packed_offset).expect("offset fits");
        encoded[offset] ^= 1;
        assert!(matches!(
            decode_contig_chunk(&encoded, &entries, 0),
            Err(ReferenceCatalogError::ContigIntegrity { ordinal: 0 })
        ));
    }
}
