//! Safe wrappers over the pinned private `HTSlib` shim.
//!
//! The public layer contains no unsafe code and exposes no native pointer.
//! Unsafe ownership of the narrow C ABI is confined to the private `sys`
//! module. The interface accepts only concrete local filesystem paths and
//! thread-confined handles.

#![deny(unsafe_code)]

use core::fmt;
use std::ffi::{CString, NulError};
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Write};
use std::os::fd::AsRawFd;
use std::path::{Path, PathBuf};

use crate::AlignmentRecordError;

#[cfg(test)]
use crate::sys;
use crate::sys::{NativeBgzfWriter, NativeCompression, NativeReader};
pub use crate::sys::{NativeError, NativeStatus};

/// One content-derived source compression class.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Compression {
    /// Uncompressed bytes.
    Plain,
    /// Generic RFC 1952 gzip.
    Gzip,
    /// Block gzip.
    Bgzf,
}

impl From<NativeCompression> for Compression {
    fn from(value: NativeCompression) -> Self {
        match value {
            NativeCompression::Plain => Self::Plain,
            NativeCompression::Gzip => Self::Gzip,
            NativeCompression::Bgzf => Self::Bgzf,
        }
    }
}

/// Stable adapter phase associated with a safe HTS failure.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum HtsOperation {
    /// Verify the exact shim ABI and `HTSlib` runtime.
    Health,
    /// Validate a caller path before native work.
    ValidatePath,
    /// Open one decoded input source.
    OpenReader,
    /// Query content-derived compression.
    DetectCompression,
    /// Decode bytes.
    Read,
    /// Open a BAM together with its BAI or CSI index.
    OpenIndexedBam,
    /// Copy and validate an indexed BAM header.
    ReadIndexedBamHeader,
    /// Select one indexed BAM reference interval.
    QueryIndexedBam,
    /// Decode one indexed BAM record.
    ReadIndexedBamRecord,
    /// Close an indexed BAM reader.
    CloseIndexedBam,
    /// Open a FASTA together with its existing FAI/GZI indexes.
    OpenIndexedFasta,
    /// Copy and validate an indexed FASTA dictionary.
    ReadIndexedFastaHeader,
    /// Fetch one indexed FASTA interval.
    FetchIndexedFasta,
    /// Close an indexed FASTA reader.
    CloseIndexedFasta,
    /// Build and synchronize a create-only BAI.
    BuildBamIndex,
    /// Explicitly close a decoded source.
    CloseReader,
    /// Encode the canonical alignment header.
    EncodeHeader,
    /// Reserve an exclusive private BAM staging path.
    CreateStaging,
    /// Verify that a staging path still names the exclusively created file.
    ValidateStaging,
    /// Validate distinct same-directory staging and publication paths.
    ValidatePublicationPaths,
    /// Verify that a create-only publication target is absent.
    ValidateTarget,
    /// Open the native BAM encoder and write its header.
    OpenBam,
    /// Encode one canonical alignment record.
    EncodeRecord,
    /// Write one canonical record through `HTSlib`.
    WriteRecord,
    /// Finalize one BAM stream.
    FinishBam,
    /// Open one plain-text or BGZF text staging stream.
    OpenTextOutput,
    /// Finalize one plain-text or BGZF text staging stream.
    FinishTextOutput,
    /// Synchronize a completed BAM before publication.
    SyncBam,
    /// Synchronize a completed generic output before publication.
    SyncOutput,
    /// Atomically create the final BAM target.
    PublishBam,
    /// Atomically create a final generic output target.
    PublishOutput,
    /// Remove a previously published generic output during rollback.
    RollbackOutput,
    /// Remove an adapter-owned staging path.
    Cleanup,
}

/// Stable high-level class for one adapter error.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum HtsErrorKind {
    /// The path requests stdin, a URL/plugin scheme, or a non-file source.
    UnsupportedPath,
    /// The path cannot be represented by the UTF-8 boundary.
    PathEncoding,
    /// The path contains an embedded NUL and cannot cross the C ABI.
    PathContainsNul,
    /// A direct filesystem operation failed.
    Io(io::ErrorKind),
    /// Canonical SAM encoding failed.
    Encode,
    /// The native shim returned this typed status.
    Native(NativeStatus),
    /// The one-based record counter overflowed.
    RecordCountOverflow,
    /// The safe writer is terminal after an earlier failure.
    Terminal,
    /// A staging path was removed or replaced after exclusive creation.
    StagingIdentityChanged,
    /// Target and staging are not distinct siblings.
    PublicationPathMismatch,
}

/// A path- and record-contextual safe adapter error.
#[derive(Debug)]
pub struct HtsError {
    operation: HtsOperation,
    path: PathBuf,
    record_ordinal: Option<u64>,
    kind: HtsErrorKind,
    system_errno: Option<i32>,
    source: Option<Box<HtsErrorSource>>,
}

#[derive(Debug)]
enum HtsErrorSource {
    Io(io::Error),
    Encode(AlignmentRecordError),
    Native(NativeError),
}

impl HtsError {
    /// Returns the operation that failed.
    #[must_use]
    pub const fn operation(&self) -> HtsOperation {
        self.operation
    }

    /// Returns the exact caller or staging path.
    #[must_use]
    pub fn path(&self) -> &Path {
        &self.path
    }

    /// Returns the one-based record ordinal for a record-local failure.
    #[must_use]
    pub const fn record_ordinal(&self) -> Option<u64> {
        self.record_ordinal
    }

    /// Returns the stable high-level error class.
    #[must_use]
    pub const fn kind(&self) -> HtsErrorKind {
        self.kind
    }

    /// Returns a nonzero native `errno` when one was reported.
    #[must_use]
    pub const fn system_errno(&self) -> Option<i32> {
        self.system_errno
    }

    /// Returns the copied native diagnostic when applicable.
    #[must_use]
    pub fn native_message(&self) -> Option<&str> {
        match self.source.as_deref() {
            Some(HtsErrorSource::Native(source)) => Some(source.message()),
            _ => None,
        }
    }
}

impl fmt::Display for HtsError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "HTS {:?} failed for {}",
            self.operation,
            self.path.display()
        )?;
        if let Some(ordinal) = self.record_ordinal {
            write!(formatter, " at record {ordinal}")?;
        }
        write!(formatter, ": {:?}", self.kind)?;
        if let Some(message) = self.native_message() {
            write!(formatter, " ({message})")?;
        }
        Ok(())
    }
}

impl std::error::Error for HtsError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self.source.as_deref()? {
            HtsErrorSource::Io(source) => Some(source),
            HtsErrorSource::Encode(source) => Some(source),
            HtsErrorSource::Native(source) => Some(source),
        }
    }
}

/// A thread-confined, content-detected decoded byte source.
pub struct DecodedReader {
    path: PathBuf,
    compression: Compression,
    native: NativeReader,
}

impl DecodedReader {
    /// Opens one concrete local regular-file path read-only.
    ///
    /// Missing files are delegated to the native open so `errno` is preserved.
    /// Existing directories, FIFOs, devices, `-`, and URL/plugin schemes are
    /// rejected before native work.
    ///
    /// # Errors
    ///
    /// Returns a typed path, open, or compression-detection failure.
    pub fn open(path: impl AsRef<Path>) -> Result<Self, HtsError> {
        let path = path.as_ref().to_path_buf();
        validate_reader_path(&path)?;
        let c_path = path_cstring(&path)?;
        let native = NativeReader::open(&c_path)
            .map_err(|source| native_error(HtsOperation::OpenReader, &path, None, source))?;
        let compression = native
            .compression()
            .map(Compression::from)
            .map_err(|source| native_error(HtsOperation::DetectCompression, &path, None, source))?;
        Ok(Self {
            path,
            compression,
            native,
        })
    }

    /// Returns the content-derived compression class.
    #[must_use]
    pub const fn compression(&self) -> Compression {
        self.compression
    }

    /// Returns the caller path.
    #[must_use]
    pub fn path(&self) -> &Path {
        &self.path
    }

    /// Decodes at most `buffer.len()` bytes.
    ///
    /// # Errors
    ///
    /// Returns a terminal native decode failure with no bytes from the failing
    /// call published.
    pub fn read_decoded(&mut self, buffer: &mut [u8]) -> Result<usize, HtsError> {
        self.native
            .read(buffer)
            .map_err(|source| native_error(HtsOperation::Read, &self.path, None, source))
    }

    /// Explicitly closes the source and checks the native close result.
    ///
    /// # Errors
    ///
    /// Returns a copied native close failure. The source is closed even on
    /// error and is consumed by this method.
    pub fn close(mut self) -> Result<(), HtsError> {
        self.native
            .close()
            .map_err(|source| native_error(HtsOperation::CloseReader, &self.path, None, source))
    }
}

impl Read for DecodedReader {
    fn read(&mut self, buffer: &mut [u8]) -> io::Result<usize> {
        self.read_decoded(buffer).map_err(io::Error::other)
    }
}

/// A safe BGZF encoder anchored to an already-open output file.
///
/// `HTSlib` opens a private descriptor through `/proc/self/fd`; the supplied
/// `File` remains owned here so callers can synchronize it after finalization.
pub struct BgzfWriter {
    anchor: Option<File>,
    native: Option<NativeBgzfWriter>,
}

impl BgzfWriter {
    /// Starts BGZF encoding on `file` with optional private compression workers.
    ///
    /// `compression_threads == 0` selects synchronous compression. Values above
    /// 64 are rejected by the native shim.
    ///
    /// # Errors
    ///
    /// Returns an I/O error if `HTSlib` cannot duplicate/open the file descriptor
    /// or initialize the requested compression workers.
    pub fn from_file(file: File, compression_threads: u32) -> io::Result<Self> {
        let descriptor_path =
            CString::new(format!("/proc/self/fd/{}", file.as_raw_fd())).map_err(|_| {
                io::Error::other("numeric file-descriptor path unexpectedly contains NUL")
            })?;
        let native = NativeBgzfWriter::open(&descriptor_path, compression_threads)
            .map_err(io::Error::other)?;
        Ok(Self {
            anchor: Some(file),
            native: Some(native),
        })
    }

    /// Finalizes compression, writes the canonical BGZF EOF block, and returns
    /// the still-open original file for synchronization.
    ///
    /// # Errors
    ///
    /// Returns an I/O error if `HTSlib` cannot finish the stream. The original
    /// file is closed on failure.
    pub fn finish(mut self) -> io::Result<File> {
        self.native_mut()?.finish().map_err(io::Error::other)?;
        drop(self.native.take());
        self.anchor
            .take()
            .ok_or_else(|| io::Error::other("BGZF writer lost its file anchor"))
    }

    fn native_mut(&mut self) -> io::Result<&mut NativeBgzfWriter> {
        self.native
            .as_mut()
            .ok_or_else(|| io::Error::other("BGZF writer is already finished"))
    }
}

impl Write for BgzfWriter {
    fn write(&mut self, buffer: &[u8]) -> io::Result<usize> {
        self.native_mut()?.write(buffer).map_err(io::Error::other)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.native_mut()?.flush().map_err(io::Error::other)
    }
}

pub(crate) struct DecodedBufReader {
    inner: BufReader<DecodedReader>,
}

impl DecodedBufReader {
    const DECODE_BUFFER_BYTES: usize = 1 << 20;

    pub(crate) fn open(path: impl AsRef<Path>) -> Result<Self, HtsError> {
        DecodedReader::open(path).map(|reader| Self {
            inner: BufReader::with_capacity(Self::DECODE_BUFFER_BYTES, reader),
        })
    }

    pub(crate) fn path(&self) -> &Path {
        self.inner.get_ref().path()
    }

    pub(crate) fn compression(&self) -> Compression {
        self.inner.get_ref().compression()
    }

    pub(crate) fn close(self) -> Result<(), HtsError> {
        self.inner.into_inner().close()
    }
}

impl Read for DecodedBufReader {
    fn read(&mut self, buffer: &mut [u8]) -> io::Result<usize> {
        self.inner.read(buffer)
    }
}

impl BufRead for DecodedBufReader {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        self.inner.fill_buf()
    }

    fn consume(&mut self, amount: usize) {
        self.inner.consume(amount);
    }
}

/// Returns the exact linked shim ABI and `HTSlib` runtime versions.
///
/// # Errors
///
/// Returns a native invalid-state failure when the linked runtime is not the
/// audited version.
#[cfg(test)]
fn native_versions() -> Result<(u32, String), HtsError> {
    sys::health_check().map_err(|source| {
        native_error(HtsOperation::Health, Path::new("<runtime>"), None, source)
    })?;
    let runtime = sys::runtime_version().map_err(|source| {
        native_error(HtsOperation::Health, Path::new("<runtime>"), None, source)
    })?;
    Ok((sys::shim_abi_version(), runtime))
}

pub(crate) fn absolute_path(path: &Path, operation: HtsOperation) -> Result<PathBuf, HtsError> {
    validate_path_spelling(path)?;
    bsbit_io::absolute_path(path).map_err(|source| io_error(operation, path, None, source))
}

pub(crate) fn validate_reader_path(path: &Path) -> Result<(), HtsError> {
    validate_path_spelling(path)?;
    bsbit_io::validate_regular_file_or_absent(path).map_err(|source| {
        if source.kind() == io::ErrorKind::Unsupported {
            simple_error(
                HtsOperation::ValidatePath,
                path,
                None,
                HtsErrorKind::UnsupportedPath,
            )
        } else {
            io_error(HtsOperation::ValidatePath, path, None, source)
        }
    })
}

pub(crate) fn path_cstring(path: &Path) -> Result<CString, HtsError> {
    validate_path_spelling(path)?;
    let text = path.to_str().ok_or_else(|| {
        simple_error(
            HtsOperation::ValidatePath,
            path,
            None,
            HtsErrorKind::PathEncoding,
        )
    })?;
    CString::new(text).map_err(|source| nul_error(path, source))
}

fn validate_path_spelling(path: &Path) -> Result<(), HtsError> {
    let text = path.to_str().ok_or_else(|| {
        simple_error(
            HtsOperation::ValidatePath,
            path,
            None,
            HtsErrorKind::PathEncoding,
        )
    })?;
    if text == "-" || text.contains("://") {
        return Err(simple_error(
            HtsOperation::ValidatePath,
            path,
            None,
            HtsErrorKind::UnsupportedPath,
        ));
    }
    if text.contains('\0') {
        return Err(simple_error(
            HtsOperation::ValidatePath,
            path,
            None,
            HtsErrorKind::PathContainsNul,
        ));
    }
    Ok(())
}

pub(crate) fn nul_error(path: &Path, _source: NulError) -> HtsError {
    simple_error(
        HtsOperation::ValidatePath,
        path,
        None,
        HtsErrorKind::PathContainsNul,
    )
}

pub(crate) fn simple_error(
    operation: HtsOperation,
    path: &Path,
    record_ordinal: Option<u64>,
    kind: HtsErrorKind,
) -> HtsError {
    HtsError {
        operation,
        path: path.to_path_buf(),
        record_ordinal,
        kind,
        system_errno: None,
        source: None,
    }
}

pub(crate) fn io_error(
    operation: HtsOperation,
    path: &Path,
    record_ordinal: Option<u64>,
    source: io::Error,
) -> HtsError {
    HtsError {
        operation,
        path: path.to_path_buf(),
        record_ordinal,
        kind: HtsErrorKind::Io(source.kind()),
        system_errno: source.raw_os_error(),
        source: Some(Box::new(HtsErrorSource::Io(source))),
    }
}

pub(crate) fn encode_error(
    operation: HtsOperation,
    path: &Path,
    record_ordinal: Option<u64>,
    source: AlignmentRecordError,
) -> HtsError {
    HtsError {
        operation,
        path: path.to_path_buf(),
        record_ordinal,
        kind: HtsErrorKind::Encode,
        system_errno: None,
        source: Some(Box::new(HtsErrorSource::Encode(source))),
    }
}

pub(crate) fn native_error(
    operation: HtsOperation,
    path: &Path,
    record_ordinal: Option<u64>,
    source: NativeError,
) -> HtsError {
    let system_errno = (source.system_errno() != 0).then_some(source.system_errno());
    HtsError {
        operation,
        path: path.to_path_buf(),
        record_ordinal,
        kind: HtsErrorKind::Native(source.status()),
        system_errno,
        source: Some(Box::new(HtsErrorSource::Native(source))),
    }
}

#[cfg(test)]
mod tests {
    use std::fs::{self, OpenOptions};
    use std::io::{Read, Write};
    use std::path::Path;
    use std::time::{SystemTime, UNIX_EPOCH};

    use super::{
        BgzfWriter, Compression, DecodedReader, HtsErrorKind, HtsOperation, NativeStatus,
        native_versions,
    };
    use crate::{
        IndexedBamHeader, IndexedBamReader, IndexedBamRecord, IndexedFastaReader,
        TextOutputCompression, TextStagingWriter,
    };

    #[test]
    fn exact_native_versions_are_linked() {
        assert_eq!(
            native_versions().expect("versions"),
            (3, String::from("1.24"))
        );
    }

    #[test]
    fn bgzf_writer_preserves_anchor_and_emits_detectable_stream() {
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock")
            .as_nanos();
        let path = std::env::temp_dir().join(format!(
            "bsbit-bgzf-writer-{}-{nonce}.data",
            std::process::id()
        ));
        let file = OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&path)
            .expect("exclusive output");
        let mut writer = BgzfWriter::from_file(file, 2).expect("BGZF writer");
        writer.write_all(b"alpha\nbeta\n").expect("BGZF payload");
        let file = writer.finish().expect("BGZF finalization");
        file.sync_all().expect("BGZF synchronization");
        drop(file);

        let mut reader = DecodedReader::open(&path).expect("BGZF reader");
        assert_eq!(reader.compression(), Compression::Bgzf);
        let mut decoded = Vec::new();
        reader.read_to_end(&mut decoded).expect("decoded payload");
        assert_eq!(decoded, b"alpha\nbeta\n");
        reader.close().expect("BGZF close");
        fs::remove_file(path).expect("fixture cleanup");
    }

    #[test]
    fn indexed_reader_copies_header_and_reuses_region_iterator() {
        let bam =
            Path::new(env!("CARGO_MANIFEST_DIR")).join("../../external/htslib/test/range.bam");
        let mut reader = IndexedBamReader::open(&bam).expect("indexed fixture opens");
        assert!(reader.header().text().starts_with(b"@HD\tVN:1.4"));
        assert_eq!(reader.header().references().len(), 7);
        assert_eq!(reader.header().references()[0].name(), b"CHROMOSOME_I");
        assert_eq!(reader.header().references()[0].length(), 1_009_800);

        reader.query(0, 900, 950).expect("first region queries");
        let first = reader
            .next_record()
            .expect("first region reads")
            .expect("first region contains a record");
        assert_eq!(first.reference_id(), 0);
        assert_eq!(first.position(), 913);
        assert_eq!(first.query_name(), b"HS18_09653:4:1315:19857:61712");
        assert_eq!(first.sequence_length(), 100);
        assert_eq!(first.packed_sequence().len(), 50);
        assert_eq!(first.quality().len(), 100);
        assert_eq!(first.sequence_code(0), Some(1));
        assert_eq!(first.sequence_code(100), None);
        assert!(!first.cigar().is_empty());
        assert!(!first.auxiliary().is_empty());
        let mut sequence = Vec::new();
        first
            .decode_sequence_into(&mut sequence)
            .expect("sequence decodes through the format adapter");
        assert_eq!(sequence.len(), 100);
        assert!(sequence.iter().all(|base| b"ACGTN".contains(base)));
        let mut cigar = Vec::new();
        first
            .decode_cigar_into(&mut cigar)
            .expect("CIGAR decodes through the format adapter");
        assert!(!cigar.is_empty());
        assert!(cigar.iter().all(|run| run.length() != 0));

        reader
            .query(1, 1_100, 1_200)
            .expect("second region queries");
        let second = reader
            .next_record()
            .expect("second region reads")
            .expect("second region contains a record");
        assert_eq!(second.reference_id(), 1);
        assert!(second.position() < 1_200);
        reader.close().expect("indexed fixture closes");
    }

    #[test]
    fn indexed_reader_reuses_caller_owned_record_buffers() {
        let bam =
            Path::new(env!("CARGO_MANIFEST_DIR")).join("../../external/htslib/test/range.bam");
        let mut reader = IndexedBamReader::open(&bam).expect("indexed fixture opens");
        let mut record = IndexedBamRecord::default();

        reader.query(0, 900, 950).expect("first region queries");
        assert!(
            reader
                .next_record_into(&mut record)
                .expect("first record reads")
        );
        let query_name_storage = record.query_name.as_ptr();
        let cigar_storage = record.cigar.as_ptr();
        let sequence_storage = record.packed_sequence.as_ptr();
        let quality_storage = record.quality.as_ptr();
        let auxiliary_storage = record.auxiliary.as_ptr();
        let first = record.clone();

        reader.query(0, 900, 950).expect("same region requeries");
        assert!(
            reader
                .next_record_into(&mut record)
                .expect("same first record rereads")
        );
        assert_eq!(record, first);
        assert_eq!(record.query_name.as_ptr(), query_name_storage);
        assert_eq!(record.cigar.as_ptr(), cigar_storage);
        assert_eq!(record.packed_sequence.as_ptr(), sequence_storage);
        assert_eq!(record.quality.as_ptr(), quality_storage);
        assert_eq!(record.auxiliary.as_ptr(), auxiliary_storage);

        let capacities_before_eof = loop {
            let capacities = (
                record.query_name.capacity(),
                record.cigar.capacity(),
                record.packed_sequence.capacity(),
                record.quality.capacity(),
                record.auxiliary.capacity(),
            );
            if !reader
                .next_record_into(&mut record)
                .expect("remaining region records read")
            {
                break capacities;
            }
        };
        assert!(record.query_name.is_empty());
        assert!(record.cigar.is_empty());
        assert!(record.packed_sequence.is_empty());
        assert!(record.quality.is_empty());
        assert!(record.auxiliary.is_empty());
        assert_eq!(
            capacities_before_eof,
            (
                record.query_name.capacity(),
                record.cigar.capacity(),
                record.packed_sequence.capacity(),
                record.quality.capacity(),
                record.auxiliary.capacity(),
            )
        );
        reader.close().expect("indexed fixture closes");
    }

    #[test]
    fn indexed_fasta_reader_requires_existing_index_and_fetches_half_open_regions() {
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock")
            .as_nanos();
        let directory = std::env::temp_dir().join(format!(
            "bsbit-indexed-fasta-{}-{nonce}",
            std::process::id()
        ));
        fs::create_dir(&directory).expect("fixture directory");
        let fasta = directory.join("reference.fa");
        fs::write(&fasta, b">chr1\nACGT\nTGCA\n>chr2\nNNCC\n").expect("FASTA fixture");
        let missing_index = IndexedFastaReader::open(&fasta)
            .err()
            .expect("FAI is required");
        assert_eq!(missing_index.operation(), HtsOperation::OpenIndexedFasta);
        fs::write(
            fasta.with_extension("fa.fai"),
            b"chr1\t8\t6\t4\t5\nchr2\t4\t22\t4\t5\n",
        )
        .expect("FAI fixture");

        let mut reader = IndexedFastaReader::open(&fasta).expect("indexed FASTA opens");
        assert_eq!(reader.references().len(), 2);
        assert_eq!(reader.references()[0].name(), b"chr1");
        assert_eq!(reader.references()[0].length(), 8);
        assert_eq!(reader.fetch(0, 2, 7).expect("interval fetches"), b"GTTGC");
        assert_eq!(
            reader.fetch(1, 0, 4).expect("ambiguous bases fetch"),
            b"NNCC"
        );
        assert!(reader.fetch(0, 7, 9).is_err());
        reader.close().expect("indexed FASTA closes");
        fs::remove_dir_all(directory).expect("fixture cleanup");
    }

    #[test]
    fn indexed_header_queries_sam_sort_and_program_fields() {
        let header = IndexedBamHeader {
            text: b"@HD\tVN:1.6\tSO:coordinate\n@RG\tID:lane-1\tSM:donor-A\n@RG\tSM:donor-A\tID:lane-2\n@RG\tID:unknown\n@PG\tID:other\tPN:x\n@PG\tPN:bsbit\tID:bsbit\n".to_vec(),
            references: Vec::new(),
        };
        assert!(header.is_coordinate_sorted());
        assert!(header.has_program(b"bsbit", b"bsbit"));
        assert!(!header.has_program(b"other", b"bsbit"));
        assert_eq!(
            header.read_group_sample_names().collect::<Vec<_>>(),
            vec![b"donor-A".as_slice(), b"donor-A".as_slice()]
        );

        let unsorted = IndexedBamHeader {
            text: b"@HD\tVN:1.6\tSO:queryname\n".to_vec(),
            references: Vec::new(),
        };
        assert!(!unsorted.is_coordinate_sorted());
    }

    #[test]
    fn streaming_text_output_publishes_bgzf_create_only_and_rolls_back() {
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock")
            .as_nanos();
        let directory =
            std::env::temp_dir().join(format!("bsbit-text-output-{}-{nonce}", std::process::id()));
        fs::create_dir(&directory).expect("directory");
        let target = directory.join("calls.vcf.gz");
        let mut writer = TextStagingWriter::create_sibling(&target, TextOutputCompression::Bgzf, 1)
            .expect("staging writer");
        let staging = writer.staging_path().to_path_buf();
        writer.write_all(b"header\nregion-1\n").expect("first rows");
        writer.write_all(b"region-2\n").expect("second rows");
        let publication = writer
            .finish()
            .expect("final text")
            .publish_create_new()
            .expect("create-only publication");
        assert_eq!(publication.target_path(), target);
        assert!(!staging.exists());

        let mut decoded = DecodedReader::open(&target).expect("published BGZF opens");
        assert_eq!(decoded.compression(), Compression::Bgzf);
        let mut bytes = Vec::new();
        decoded
            .read_to_end(&mut bytes)
            .expect("published BGZF decodes");
        decoded.close().expect("published BGZF closes");
        assert_eq!(bytes, b"header\nregion-1\nregion-2\n");
        assert!(
            TextStagingWriter::create_sibling(&target, TextOutputCompression::Plain, 0).is_err()
        );
        publication.rollback().expect("published target rolls back");
        assert!(!target.exists());

        let abandoned_target = directory.join("abandoned.txt");
        let abandoned_staging = {
            let mut abandoned = TextStagingWriter::create_sibling(
                &abandoned_target,
                TextOutputCompression::Plain,
                0,
            )
            .expect("abandoned staging");
            abandoned.write_all(b"partial").expect("partial bytes");
            abandoned.staging_path().to_path_buf()
        };
        assert!(!abandoned_staging.exists());
        assert!(!abandoned_target.exists());
        fs::remove_dir(directory).expect("directory cleanup");
    }

    #[test]
    fn streaming_bgzf_vcf_builds_an_actual_tabix_index() {
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock")
            .as_nanos();
        let directory =
            std::env::temp_dir().join(format!("bsbit-tabix-vcf-{}-{nonce}", std::process::id()));
        fs::create_dir(&directory).expect("test directory");
        let target = directory.join("calls.vcf.gz");
        let index = directory.join("calls.vcf.gz.tbi");
        let mut writer = TextStagingWriter::create_sibling(&target, TextOutputCompression::Bgzf, 1)
            .expect("VCF staging");
        writer
            .write_all(
                b"##fileformat=VCFv4.3\n##contig=<ID=chr1,length=10>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\nchr1\t1\t.\tA\tC\t30\tPASS\t.\tGT\t0/1\n",
            )
            .expect("VCF rows");
        writer
            .finish()
            .expect("VCF finalization")
            .publish_create_new()
            .expect("VCF publication");

        let target_c = super::path_cstring(&target).expect("VCF path");
        let index_c = super::path_cstring(&index).expect("index path");
        super::sys::build_tabix_index(&target_c, &index_c, 0, 1)
            .expect("HTSlib builds a VCF tabix index");
        assert!(fs::metadata(&index).expect("tabix metadata").len() > 0);

        fs::remove_file(index).expect("index cleanup");
        fs::remove_file(target).expect("VCF cleanup");
        fs::remove_dir(directory).expect("directory cleanup");
    }

    #[test]
    fn indexed_reader_rejects_invalid_intervals_before_native_query() {
        let bam =
            Path::new(env!("CARGO_MANIFEST_DIR")).join("../../external/htslib/test/range.bam");
        let mut reader = IndexedBamReader::open(&bam).expect("indexed fixture opens");
        let error = reader.query(0, 10, 10).expect_err("empty query rejects");
        assert_eq!(error.operation(), HtsOperation::QueryIndexedBam);
        assert_eq!(
            error.kind(),
            HtsErrorKind::Native(NativeStatus::InvalidArgument)
        );
        reader.close().expect("indexed fixture closes");
    }
}
