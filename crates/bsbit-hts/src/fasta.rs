//! FASTA record model and bounded streaming parser.

use std::io::{self, BufRead, Write};
use std::path::{Path, PathBuf};

use bsbit_core::alphabet::Base;
use bsbit_core::sequence::NormalizedSequence;

use crate::htslib::{
    Compression, DecodedBufReader, HtsError, HtsErrorKind, HtsOperation, native_error,
    path_cstring, simple_error, validate_reader_path,
};
use crate::sys::{NativeIndexedFastaReader, NativeStatus};
use crate::text::{
    BoundedLineReader, ErrorContext, PairSourceSide, PhysicalLine, RecordField, RecordName,
    RecordOrdinal, TextRecordAllocation, TextRecordError, TextRecordErrorKind, TextRecordFormat,
    TextRecordLimits, TextRecordResource, check_limit, checked_add, line_error,
    normalize_sequence_line, parse_header, record_limit_context, storage_len, write_record_name,
    write_sequence,
};

/// One owned normalized FASTA record.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct FastaRecord(FastaRecordData);

impl FastaRecord {
    /// Returns the zero-based input ordinal.
    #[must_use]
    pub const fn ordinal(&self) -> RecordOrdinal {
        self.0.ordinal()
    }

    /// Returns the retained header name and description.
    #[must_use]
    pub const fn record_name(&self) -> &RecordName {
        self.0.record_name()
    }

    /// Returns the normalized sequence.
    #[must_use]
    pub const fn sequence(&self) -> &NormalizedSequence {
        self.0.sequence()
    }

    /// Writes one canonical FASTA record using uppercase sequence and LF lines.
    ///
    /// # Errors
    ///
    /// Returns the first error from `writer`.
    pub fn write_canonical<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        self.0.write_canonical(writer)
    }
}

/// Streaming parser for the accepted bounded FASTA subset.
#[derive(Debug)]
pub struct FastaReader<R> {
    core: FastaReaderCore<R>,
}

impl<R: BufRead> FastaReader<R> {
    /// Creates a parser with a complete explicit limit set.
    #[must_use]
    pub const fn new(reader: R, limits: TextRecordLimits) -> Self {
        Self {
            core: FastaReaderCore::new(reader, limits),
        }
    }

    /// Returns the underlying buffered source without changing parser state.
    #[must_use]
    pub const fn get_ref(&self) -> &R {
        self.core.get_ref()
    }

    /// Recovers the underlying source and discards parser state.
    #[must_use]
    pub fn into_inner(self) -> R {
        self.core.into_inner()
    }

    /// Parses and returns the next FASTA record.
    ///
    /// # Errors
    ///
    /// Returns a structured syntax, resource, allocation, arithmetic, or I/O
    /// error. After the first error, later calls return `TerminalState`.
    pub fn next_record(&mut self) -> Result<Option<FastaRecord>, TextRecordError> {
        self.core
            .next_record()
            .map(|record| record.map(FastaRecord))
    }
}

/// One FASTA index dictionary entry.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct IndexedFastaReference {
    name: Vec<u8>,
    length: u64,
}

impl IndexedFastaReference {
    /// Returns the exact FASTA sequence name bytes.
    #[must_use]
    pub fn name(&self) -> &[u8] {
        &self.name
    }

    /// Returns the indexed FASTA sequence length.
    #[must_use]
    pub const fn length(&self) -> u64 {
        self.length
    }
}

/// A thread-confined random-access FASTA reader backed by existing FAI/GZI indexes.
pub struct IndexedFastaReader {
    path: PathBuf,
    references: Vec<IndexedFastaReference>,
    native: NativeIndexedFastaReader,
}

impl IndexedFastaReader {
    /// Opens a local FASTA and its adjacent `.fai`; BGZF FASTA also requires `.gzi`.
    ///
    /// This operation never creates or replaces either index.
    ///
    /// # Errors
    ///
    /// Returns a path, FASTA, FAI, or GZI loading failure.
    pub fn open(path: impl AsRef<Path>) -> Result<Self, HtsError> {
        let path = path.as_ref().to_path_buf();
        validate_reader_path(&path)?;
        let c_path = path_cstring(&path)?;
        let native = NativeIndexedFastaReader::open(&c_path)
            .map_err(|source| native_error(HtsOperation::OpenIndexedFasta, &path, None, source))?;
        let native_references = native.references().map_err(|source| {
            native_error(HtsOperation::ReadIndexedFastaHeader, &path, None, source)
        })?;
        let mut references = Vec::with_capacity(native_references.len());
        for reference in native_references {
            let length = u64::try_from(reference.length).map_err(|_| {
                simple_error(
                    HtsOperation::ReadIndexedFastaHeader,
                    &path,
                    None,
                    HtsErrorKind::Native(NativeStatus::HeaderFailed),
                )
            })?;
            references.push(IndexedFastaReference {
                name: reference.name,
                length,
            });
        }
        Ok(Self {
            path,
            references,
            native,
        })
    }

    /// Returns the copied FASTA dictionary in index order.
    #[must_use]
    pub fn references(&self) -> &[IndexedFastaReference] {
        &self.references
    }

    /// Fetches one zero-based, half-open interval.
    ///
    /// # Errors
    ///
    /// Returns an argument error for an empty or out-of-range interval, or a
    /// terminal native decoding error when the indexed sequence cannot be read.
    pub fn fetch(&mut self, reference_id: u32, start: u64, end: u64) -> Result<Vec<u8>, HtsError> {
        let reference_ordinal = usize::try_from(reference_id).map_err(|_| {
            simple_error(
                HtsOperation::FetchIndexedFasta,
                &self.path,
                None,
                HtsErrorKind::Native(NativeStatus::InvalidArgument),
            )
        })?;
        let Some(reference) = self.references.get(reference_ordinal) else {
            return Err(simple_error(
                HtsOperation::FetchIndexedFasta,
                &self.path,
                None,
                HtsErrorKind::Native(NativeStatus::InvalidArgument),
            ));
        };
        if start >= end || end > reference.length {
            return Err(simple_error(
                HtsOperation::FetchIndexedFasta,
                &self.path,
                None,
                HtsErrorKind::Native(NativeStatus::InvalidArgument),
            ));
        }
        let reference_id = i32::try_from(reference_id).map_err(|_| {
            simple_error(
                HtsOperation::FetchIndexedFasta,
                &self.path,
                None,
                HtsErrorKind::Native(NativeStatus::InvalidArgument),
            )
        })?;
        let start = i64::try_from(start).map_err(|_| {
            simple_error(
                HtsOperation::FetchIndexedFasta,
                &self.path,
                None,
                HtsErrorKind::Native(NativeStatus::InvalidArgument),
            )
        })?;
        let end = i64::try_from(end).map_err(|_| {
            simple_error(
                HtsOperation::FetchIndexedFasta,
                &self.path,
                None,
                HtsErrorKind::Native(NativeStatus::InvalidArgument),
            )
        })?;
        self.native
            .fetch(reference_id, start, end)
            .map_err(|source| {
                native_error(HtsOperation::FetchIndexedFasta, &self.path, None, source)
            })
    }

    /// Explicitly closes the FASTA and its indexes.
    ///
    /// # Errors
    ///
    /// Returns a copied native close failure.
    pub fn close(mut self) -> Result<(), HtsError> {
        self.native.close().map_err(|source| {
            native_error(HtsOperation::CloseIndexedFasta, &self.path, None, source)
        })
    }
}

/// A bounded FASTA parser over one content-detected local decoded source.
pub struct DecodedFastaReader {
    parser: FastaReader<DecodedBufReader>,
}

impl DecodedFastaReader {
    /// Opens one local source and composes it with the validated Rust parser.
    ///
    /// # Errors
    ///
    /// Returns a typed path, native-open, or compression-detection failure.
    pub fn open(path: impl AsRef<Path>, limits: TextRecordLimits) -> Result<Self, HtsError> {
        let source = DecodedBufReader::open(path)?;
        Ok(Self {
            parser: FastaReader::new(source, limits),
        })
    }

    /// Returns the concrete caller path.
    #[must_use]
    pub fn path(&self) -> &Path {
        self.parser.get_ref().path()
    }

    /// Returns the source compression detected from content.
    #[must_use]
    pub fn compression(&self) -> Compression {
        self.parser.get_ref().compression()
    }

    /// Parses the next bounded FASTA record.
    ///
    /// Native decode errors are retained as the source of
    /// `TextRecordErrorKind::Io`; parser terminal-state semantics are unchanged.
    ///
    /// # Errors
    ///
    /// Returns the validated parser's structured record error.
    pub fn next_record(&mut self) -> Result<Option<FastaRecord>, TextRecordError> {
        self.parser.next_record()
    }

    /// Discards parser state and explicitly checks native source close.
    ///
    /// # Errors
    ///
    /// Returns a copied native close failure.
    pub fn close(self) -> Result<(), HtsError> {
        self.parser.into_inner().close()
    }
}

/// One owned normalized FASTA record.
#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct FastaRecordData {
    ordinal: RecordOrdinal,
    name: RecordName,
    sequence: NormalizedSequence,
}

impl FastaRecordData {
    /// Returns the zero-based input ordinal.
    #[must_use]
    pub const fn ordinal(&self) -> RecordOrdinal {
        self.ordinal
    }

    /// Returns the retained header name and description.
    #[must_use]
    pub const fn record_name(&self) -> &RecordName {
        &self.name
    }

    /// Returns the normalized sequence.
    #[must_use]
    pub const fn sequence(&self) -> &NormalizedSequence {
        &self.sequence
    }

    /// Writes one canonical FASTA record using uppercase sequence and LF lines.
    ///
    /// The canonical header uses one ASCII space before a nonempty description,
    /// and the sequence is emitted on one physical line.
    ///
    /// # Errors
    ///
    /// Returns the first error from `writer`.
    pub fn write_canonical<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        writer.write_all(b">")?;
        write_record_name(writer, &self.name)?;
        writer.write_all(b"\n")?;
        write_sequence(writer, &self.sequence)?;
        writer.write_all(b"\n")
    }
}

/// Streaming parser for the accepted bounded FASTA subset.
#[derive(Debug)]
pub(crate) struct FastaReaderCore<R> {
    lines: BoundedLineReader<R>,
    limits: TextRecordLimits,
    pending_header: Option<PhysicalLine>,
    records_emitted: u64,
    total_bases: u64,
    failed: Option<ErrorContext>,
}

impl<R: BufRead> FastaReaderCore<R> {
    /// Creates a parser with a complete explicit limit set.
    #[must_use]
    pub const fn new(reader: R, limits: TextRecordLimits) -> Self {
        Self {
            lines: BoundedLineReader::new(reader),
            limits,
            pending_header: None,
            records_emitted: 0,
            total_bases: 0,
            failed: None,
        }
    }

    /// Returns the underlying buffered source without changing parser state.
    #[must_use]
    pub const fn get_ref(&self) -> &R {
        self.lines.get_ref()
    }

    /// Recovers the underlying source and discards all parser state.
    ///
    /// This is intended for source lifecycle operations such as an explicit
    /// close. A pending FASTA header may already have been consumed from the
    /// source, so the returned source is not guaranteed to be resumable as a
    /// fresh record stream.
    #[must_use]
    pub fn into_inner(self) -> R {
        self.lines.into_inner()
    }

    /// Parses and returns the next FASTA record.
    ///
    /// # Errors
    ///
    /// Returns a structured syntax, resource, allocation, arithmetic, or I/O
    /// error. After the first error, later calls return `TerminalState` without
    /// reading more input.
    pub fn next_record(&mut self) -> Result<Option<FastaRecordData>, TextRecordError> {
        if let Some(context) = self.failed {
            return Err(context.terminal_error());
        }
        match self.next_record_inner() {
            Ok(record) => Ok(record),
            Err(error) => {
                self.failed = Some(error.context());
                Err(error)
            }
        }
    }

    #[allow(clippy::too_many_lines)]
    fn next_record_inner(&mut self) -> Result<Option<FastaRecordData>, TextRecordError> {
        let ordinal = RecordOrdinal::new(self.records_emitted);
        let header = if let Some(header) = self.pending_header.take() {
            header
        } else {
            let context = ErrorContext {
                format: TextRecordFormat::Fasta,
                side: PairSourceSide::Single,
                ordinal,
                line: self.lines.next_line_number(),
                field: RecordField::Header,
            };
            let Some(line) = self
                .lines
                .next_line(self.limits.max_line_bytes)
                .map_err(|error| line_error(context, error))?
            else {
                return Ok(None);
            };
            line
        };

        let record_context = record_limit_context(
            TextRecordFormat::Fasta,
            PairSourceSide::Single,
            ordinal,
            header.number,
        );
        let observed_records = checked_add(
            self.records_emitted,
            1,
            TextRecordResource::Records,
            record_context,
        )?;
        check_limit(
            observed_records,
            self.limits.max_records,
            TextRecordResource::Records,
            record_context,
        )?;

        let header_context = ErrorContext {
            field: RecordField::Header,
            ..record_context
        };
        let name = parse_header(&header, b'>', self.limits, header_context)?;
        let mut bases = Vec::<Base>::new();
        let mut sequence_length = 0_u64;

        loop {
            let context = ErrorContext {
                format: TextRecordFormat::Fasta,
                side: PairSourceSide::Single,
                ordinal,
                line: self.lines.next_line_number(),
                field: RecordField::Sequence,
            };
            let next = self
                .lines
                .next_line(self.limits.max_line_bytes)
                .map_err(|error| line_error(context, error))?;
            let Some(line) = next else {
                break;
            };
            if line.bytes.first() == Some(&b'>') {
                self.pending_header = Some(line);
                break;
            }
            let line_context = ErrorContext {
                line: line.number,
                ..context
            };
            if line.bytes.is_empty() {
                return Err(line_context.error(TextRecordErrorKind::EmptySequence));
            }
            let line_length = storage_len(line.bytes.len());
            let new_length = checked_add(
                sequence_length,
                line_length,
                TextRecordResource::BasesPerRecord,
                line_context,
            )?;
            check_limit(
                new_length,
                self.limits.max_bases_per_record,
                TextRecordResource::BasesPerRecord,
                line_context,
            )?;
            let normalized = normalize_sequence_line(&line, sequence_length, line_context)?;
            bases.try_reserve(normalized.bases().len()).map_err(|_| {
                line_context.error(TextRecordErrorKind::AllocationFailed {
                    allocation: TextRecordAllocation::Sequence,
                    additional: line_length,
                })
            })?;
            bases.extend_from_slice(normalized.bases());
            sequence_length = new_length;
        }

        if sequence_length == 0 {
            return Err(ErrorContext {
                field: RecordField::Sequence,
                ..header_context
            }
            .error(TextRecordErrorKind::EmptySequence));
        }
        let total_bases = checked_add(
            self.total_bases,
            sequence_length,
            TextRecordResource::TotalBases,
            record_context,
        )?;
        check_limit(
            total_bases,
            self.limits.max_total_bases,
            TextRecordResource::TotalBases,
            record_context,
        )?;

        self.records_emitted = observed_records;
        self.total_bases = total_bases;
        Ok(Some(FastaRecordData {
            ordinal,
            name,
            sequence: bases.into(),
        }))
    }
}
