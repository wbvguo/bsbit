//! FASTQ record models and bounded streaming parsers.

use std::io::{self, BufRead, Write};
use std::path::Path;

use bsbit_core::alphabet::Base;
use bsbit_core::sequence::{NormalizationError, NormalizedSequence};

use crate::htslib::{Compression, DecodedBufReader, HtsError};
#[cfg(test)]
use crate::text::LineReadError;
pub use crate::text::PairSourceSide;
use crate::text::{
    BoundedLineReader, ErrorContext, PhysicalLine, RecordField, RecordName, RecordOrdinal,
    TextRecordAllocation, TextRecordError, TextRecordErrorKind, TextRecordFormat, TextRecordLimits,
    TextRecordResource, check_limit, checked_add, line_error, normalize_sequence_line,
    parse_header, sequence_normalization_error, storage_len, validate_header, write_record_name,
    write_sequence,
};

/// One owned normalized strict four-line FASTQ record.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct FastqRecord(FastqRecordData);

impl FastqRecord {
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

    /// Returns exact printable Sanger quality bytes.
    #[must_use]
    pub fn quality(&self) -> &[u8] {
        self.0.quality()
    }

    /// Writes one canonical four-line FASTQ record using LF line endings.
    ///
    /// # Errors
    ///
    /// Returns the first error from `writer`.
    pub fn write_canonical<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        self.0.write_canonical(writer)
    }
}

/// Contiguous storage for one validated batch of FASTQ records.
///
/// Header descriptions are fully validated but not retained; batch consumers
/// receive the alignment-relevant name token, normalized bases, and quality.
#[derive(Debug, Default)]
pub struct FastqRecordBatch(FastqRecordBatchData);

impl FastqRecordBatch {
    /// Returns the number of records in this batch.
    #[must_use]
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Returns whether this batch has no records.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Returns one borrowed record by batch-relative index.
    #[must_use]
    pub fn get(&self, index: usize) -> Option<BorrowedFastqRecord<'_>> {
        self.0.get(index).map(BorrowedFastqRecord)
    }
}

/// One borrowed record inside a [`FastqRecordBatch`].
#[derive(Clone, Copy, Debug)]
pub struct BorrowedFastqRecord<'a>(BorrowedFastqRecordData<'a>);

impl<'a> BorrowedFastqRecord<'a> {
    /// Returns the zero-based input ordinal.
    #[must_use]
    pub const fn ordinal(self) -> RecordOrdinal {
        self.0.ordinal()
    }

    /// Returns the validated header name token.
    #[must_use]
    pub const fn name(self) -> &'a [u8] {
        self.0.name()
    }

    /// Returns the normalized sequence bases.
    #[must_use]
    pub const fn sequence(self) -> &'a [Base] {
        self.0.sequence()
    }

    /// Returns the validated printable Sanger quality bytes.
    #[must_use]
    pub const fn quality(self) -> &'a [u8] {
        self.0.quality()
    }
}

/// One synchronized pair of FASTQ records.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PairedFastqRecord {
    ordinal: RecordOrdinal,
    first: FastqRecord,
    second: FastqRecord,
}

impl PairedFastqRecord {
    /// Combines independently decoded records after checking synchronization.
    #[doc(hidden)]
    #[must_use]
    pub fn from_synchronized_records(first: FastqRecord, second: FastqRecord) -> Option<Self> {
        let pair = PairedFastqRecordData::from_synchronized_records(first.0, second.0)?;
        Some(Self::from_data(pair))
    }

    pub(crate) fn from_data(pair: PairedFastqRecordData) -> Self {
        let ordinal = pair.ordinal();
        let (first, second) = pair.into_records();
        Self {
            ordinal,
            first: FastqRecord(first),
            second: FastqRecord(second),
        }
    }

    /// Returns the shared zero-based pair ordinal.
    #[must_use]
    pub const fn ordinal(&self) -> RecordOrdinal {
        self.ordinal
    }

    /// Returns the source-one record.
    #[must_use]
    pub const fn first(&self) -> &FastqRecord {
        &self.first
    }

    /// Returns the source-two record.
    #[must_use]
    pub const fn second(&self) -> &FastqRecord {
        &self.second
    }

    /// Returns the shared alignment name for the synchronized pair.
    #[must_use]
    pub fn shared_name(&self) -> &[u8] {
        let first = self.first.record_name().name();
        let second = self.second.record_name().name();
        if first == second {
            first
        } else {
            first.strip_suffix(b"/1").unwrap_or(first)
        }
    }

    /// Splits the pair into its two owned records.
    #[must_use]
    pub fn into_records(self) -> (FastqRecord, FastqRecord) {
        (self.first, self.second)
    }
}

/// Streaming parser for strict bounded four-line FASTQ.
#[derive(Debug)]
pub struct FastqReader<R> {
    core: FastqReaderCore<R>,
}

impl<R: BufRead> FastqReader<R> {
    /// Creates a non-paired parser with a complete explicit limit set.
    #[must_use]
    pub const fn new(reader: R, limits: TextRecordLimits) -> Self {
        Self {
            core: FastqReaderCore::new(reader, limits),
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

    /// Parses and returns the next FASTQ record.
    ///
    /// # Errors
    ///
    /// Returns a structured record error and enters terminal state after one
    /// failure.
    pub fn next_record(&mut self) -> Result<Option<FastqRecord>, TextRecordError> {
        self.core
            .next_record()
            .map(|record| record.map(FastqRecord))
    }

    /// Parses at most `maximum_records` fully validated records into contiguous storage.
    ///
    /// When `maximum_records` is positive, an empty batch denotes end of input.
    ///
    /// # Errors
    ///
    /// Returns the same structured validation and terminal-state failures as
    /// [`Self::next_record`].
    pub fn next_batch(
        &mut self,
        maximum_records: usize,
    ) -> Result<FastqRecordBatch, TextRecordError> {
        self.core.next_batch(maximum_records).map(FastqRecordBatch)
    }
}

/// Lockstep parser for two paired FASTQ sources.
#[derive(Debug)]
pub struct PairedFastqReader<R1, R2> {
    core: PairedFastqReaderCore<R1, R2>,
}

impl<R1: BufRead, R2: BufRead> PairedFastqReader<R1, R2> {
    /// Creates a paired parser with the same explicit limits on both sources.
    #[must_use]
    pub const fn new(first: R1, second: R2, limits: TextRecordLimits) -> Self {
        Self {
            core: PairedFastqReaderCore::new(first, second, limits),
        }
    }

    /// Returns both underlying buffered sources.
    #[must_use]
    pub const fn get_ref(&self) -> (&R1, &R2) {
        self.core.get_ref()
    }

    /// Recovers both underlying sources and discards parser state.
    #[must_use]
    pub fn into_inner(self) -> (R1, R2) {
        self.core.into_inner()
    }

    /// Parses the next synchronized pair.
    ///
    /// # Errors
    ///
    /// Returns a source or synchronization error and enters terminal state.
    pub fn next_pair(&mut self) -> Result<Option<PairedFastqRecord>, TextRecordError> {
        self.core
            .next_pair()
            .map(|pair| pair.map(PairedFastqRecord::from_data))
    }
}

/// A bounded FASTQ parser over one content-detected local decoded source.
pub struct DecodedFastqReader {
    parser: FastqReader<DecodedBufReader>,
}

impl DecodedFastqReader {
    /// Opens one local source and composes it with the validated Rust parser.
    ///
    /// # Errors
    ///
    /// Returns a typed path, native-open, or compression-detection failure.
    pub fn open(path: impl AsRef<Path>, limits: TextRecordLimits) -> Result<Self, HtsError> {
        let source = DecodedBufReader::open(path)?;
        Ok(Self {
            parser: FastqReader::new(source, limits),
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

    /// Parses the next bounded strict four-line FASTQ record.
    ///
    /// Native decode errors are retained as the source of
    /// `TextRecordErrorKind::Io`; parser terminal-state semantics are unchanged.
    ///
    /// # Errors
    ///
    /// Returns the validated parser's structured record error.
    pub fn next_record(&mut self) -> Result<Option<FastqRecord>, TextRecordError> {
        self.parser.next_record()
    }

    /// Parses one fully validated contiguous batch without per-record allocations.
    ///
    /// # Errors
    ///
    /// Returns the validated parser's structured record error.
    pub fn next_batch(
        &mut self,
        maximum_records: usize,
    ) -> Result<FastqRecordBatch, TextRecordError> {
        self.parser.next_batch(maximum_records)
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

/// A bounded paired-FASTQ parser over two independently detected sources.
pub struct DecodedPairedFastqReader {
    parser: PairedFastqReader<DecodedBufReader, DecodedBufReader>,
}

impl DecodedPairedFastqReader {
    /// Opens two local sources in first/second order and composes paired parsing.
    ///
    /// # Errors
    ///
    /// Returns the first source-open failure in call order. If opening the
    /// second source fails, the already opened first source is dropped safely.
    pub fn open(
        first_path: impl AsRef<Path>,
        second_path: impl AsRef<Path>,
        limits: TextRecordLimits,
    ) -> Result<Self, HtsError> {
        let first = DecodedBufReader::open(first_path)?;
        let second = DecodedBufReader::open(second_path)?;
        Ok(Self {
            parser: PairedFastqReader::new(first, second, limits),
        })
    }

    /// Returns both concrete caller paths in first/second order.
    #[must_use]
    pub fn paths(&self) -> (&Path, &Path) {
        let (first, second) = self.parser.get_ref();
        (first.path(), second.path())
    }

    /// Returns both content-derived compression classes in first/second order.
    #[must_use]
    pub fn compressions(&self) -> (Compression, Compression) {
        let (first, second) = self.parser.get_ref();
        (first.compression(), second.compression())
    }

    /// Parses and synchronizes the next bounded read pair.
    ///
    /// Native decode errors retain the first/second source side and an
    /// `HtsError` in their I/O source chain.
    ///
    /// # Errors
    ///
    /// Returns the validated parser's structured source, syntax, or pair error.
    pub fn next_pair(&mut self) -> Result<Option<PairedFastqRecord>, TextRecordError> {
        self.parser.next_pair()
    }

    /// Discards parser state and explicitly attempts both native closes.
    ///
    /// If both closes fail, the first-source error has deterministic priority;
    /// the second close is still attempted.
    ///
    /// # Errors
    ///
    /// Returns the first close failure in source order.
    pub fn close(self) -> Result<(), HtsError> {
        let (first, second) = self.parser.into_inner();
        let first_result = first.close();
        let second_result = second.close();
        first_result.and(second_result)
    }
}

/// One owned normalized strict four-line FASTQ record.
#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct FastqRecordData {
    ordinal: RecordOrdinal,
    name: RecordName,
    sequence: NormalizedSequence,
    quality: Box<[u8]>,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct FastqRecordOffsets {
    ordinal: RecordOrdinal,
    name_start: usize,
    name_len: usize,
    sequence_start: usize,
    sequence_len: usize,
    quality_start: usize,
}

/// Contiguous storage for one batch of FASTQ records.
///
/// Names, normalized bases, and qualities each live in one allocation for the
/// complete batch. This avoids the three independent heap allocations retained
/// by every general-purpose [`FastqRecord`].
#[derive(Debug, Default)]
pub(crate) struct FastqRecordBatchData {
    records: Vec<FastqRecordOffsets>,
    names: Vec<u8>,
    sequences: Vec<Base>,
    qualities: Vec<u8>,
}

/// One borrowed record inside a [`FastqRecordBatch`].
#[derive(Clone, Copy, Debug)]
pub(crate) struct BorrowedFastqRecordData<'a> {
    ordinal: RecordOrdinal,
    name: &'a [u8],
    sequence: &'a [Base],
    quality: &'a [u8],
}

impl FastqRecordBatchData {
    fn try_with_record_capacity(
        records: usize,
        context: ErrorContext,
    ) -> Result<Self, TextRecordError> {
        let name_bytes = records.saturating_mul(32);
        let read_bytes = records.saturating_mul(100);
        let mut batch = Self::default();
        batch.records.try_reserve_exact(records).map_err(|_| {
            context.error(TextRecordErrorKind::AllocationFailed {
                allocation: TextRecordAllocation::Record,
                additional: storage_len(records),
            })
        })?;
        batch.names.try_reserve_exact(name_bytes).map_err(|_| {
            ErrorContext {
                field: RecordField::Name,
                ..context
            }
            .error(TextRecordErrorKind::AllocationFailed {
                allocation: TextRecordAllocation::Name,
                additional: storage_len(name_bytes),
            })
        })?;
        batch.sequences.try_reserve_exact(read_bytes).map_err(|_| {
            ErrorContext {
                field: RecordField::Sequence,
                ..context
            }
            .error(TextRecordErrorKind::AllocationFailed {
                allocation: TextRecordAllocation::Sequence,
                additional: storage_len(read_bytes),
            })
        })?;
        batch.qualities.try_reserve_exact(read_bytes).map_err(|_| {
            ErrorContext {
                field: RecordField::Quality,
                ..context
            }
            .error(TextRecordErrorKind::AllocationFailed {
                allocation: TextRecordAllocation::Quality,
                additional: storage_len(read_bytes),
            })
        })?;
        Ok(batch)
    }

    /// Returns the number of records in this batch.
    #[must_use]
    pub fn len(&self) -> usize {
        self.records.len()
    }

    /// Returns whether this batch has no records.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    /// Returns one borrowed record by batch-relative index.
    #[must_use]
    pub fn get(&self, index: usize) -> Option<BorrowedFastqRecordData<'_>> {
        let offsets = *self.records.get(index)?;
        Some(BorrowedFastqRecordData {
            ordinal: offsets.ordinal,
            name: &self.names[offsets.name_start..offsets.name_start + offsets.name_len],
            sequence: &self.sequences
                [offsets.sequence_start..offsets.sequence_start + offsets.sequence_len],
            quality: &self.qualities
                [offsets.quality_start..offsets.quality_start + offsets.sequence_len],
        })
    }
}

impl<'a> BorrowedFastqRecordData<'a> {
    #[must_use]
    pub const fn ordinal(self) -> RecordOrdinal {
        self.ordinal
    }

    #[must_use]
    pub const fn name(self) -> &'a [u8] {
        self.name
    }

    #[must_use]
    pub const fn sequence(self) -> &'a [Base] {
        self.sequence
    }

    #[must_use]
    pub const fn quality(self) -> &'a [u8] {
        self.quality
    }
}

impl FastqRecordData {
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

    /// Returns exact printable Sanger quality bytes.
    #[must_use]
    pub fn quality(&self) -> &[u8] {
        &self.quality
    }

    /// Writes one canonical four-line FASTQ record using LF line endings.
    ///
    /// The sequence is uppercase and the canonical plus line has no suffix.
    /// Exact quality bytes are retained.
    ///
    /// # Errors
    ///
    /// Returns the first error from `writer`.
    pub fn write_canonical<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        writer.write_all(b"@")?;
        write_record_name(writer, &self.name)?;
        writer.write_all(b"\n")?;
        write_sequence(writer, &self.sequence)?;
        writer.write_all(b"\n+\n")?;
        writer.write_all(&self.quality)?;
        writer.write_all(b"\n")
    }
}

/// One synchronized pair of FASTQ records.
#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct PairedFastqRecordData {
    ordinal: RecordOrdinal,
    first: FastqRecordData,
    second: FastqRecordData,
}

impl PairedFastqRecordData {
    /// Combines independently decoded records after checking ordinal and name
    /// synchronization. This supports parallel decoding of paired sources.
    #[doc(hidden)]
    #[must_use]
    pub fn from_synchronized_records(
        first: FastqRecordData,
        second: FastqRecordData,
    ) -> Option<Self> {
        if first.ordinal != second.ordinal
            || !paired_names_compatible(first.name.name(), second.name.name())
        {
            return None;
        }
        Some(Self {
            ordinal: first.ordinal,
            first,
            second,
        })
    }

    /// Returns the shared zero-based pair ordinal.
    #[must_use]
    pub const fn ordinal(&self) -> RecordOrdinal {
        self.ordinal
    }

    /// Splits the pair into its two owned records.
    #[must_use]
    pub fn into_records(self) -> (FastqRecordData, FastqRecordData) {
        (self.first, self.second)
    }
}

/// Streaming parser for strict bounded four-line FASTQ.
#[derive(Debug)]
pub(crate) struct FastqReaderCore<R> {
    lines: BoundedLineReader<R>,
    reusable_lines: [Vec<u8>; 4],
    limits: TextRecordLimits,
    side: PairSourceSide,
    records_emitted: u64,
    total_bases: u64,
    failed: Option<ErrorContext>,
}

impl<R: BufRead> FastqReaderCore<R> {
    /// Creates a non-paired parser with a complete explicit limit set.
    #[must_use]
    pub const fn new(reader: R, limits: TextRecordLimits) -> Self {
        Self::with_side(reader, limits, PairSourceSide::Single)
    }

    /// Returns the underlying buffered source without changing parser state.
    #[must_use]
    pub const fn get_ref(&self) -> &R {
        self.lines.get_ref()
    }

    /// Recovers the underlying source and discards all parser state.
    ///
    /// This is intended for source lifecycle operations such as an explicit
    /// close, not for resuming parsing through a second parser.
    #[must_use]
    pub fn into_inner(self) -> R {
        self.lines.into_inner()
    }

    const fn with_side(reader: R, limits: TextRecordLimits, side: PairSourceSide) -> Self {
        Self {
            lines: BoundedLineReader::new(reader),
            reusable_lines: [Vec::new(), Vec::new(), Vec::new(), Vec::new()],
            limits,
            side,
            records_emitted: 0,
            total_bases: 0,
            failed: None,
        }
    }

    /// Parses and returns the next FASTQ record.
    ///
    /// # Errors
    ///
    /// Returns a structured syntax, resource, allocation, arithmetic, or I/O
    /// error. After the first error, later calls return `TerminalState` without
    /// reading more input.
    pub fn next_record(&mut self) -> Result<Option<FastqRecordData>, TextRecordError> {
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

    /// Parses up to `maximum_records` into contiguous validated storage.
    /// When `maximum_records` is positive, an empty batch denotes EOF.
    pub fn next_batch(
        &mut self,
        maximum_records: usize,
    ) -> Result<FastqRecordBatchData, TextRecordError> {
        if let Some(context) = self.failed {
            return Err(context.terminal_error());
        }
        let context = ErrorContext {
            format: TextRecordFormat::Fastq,
            side: self.side,
            ordinal: RecordOrdinal::new(self.records_emitted),
            line: self.lines.next_line_number(),
            field: RecordField::Record,
        };
        let mut batch =
            match FastqRecordBatchData::try_with_record_capacity(maximum_records, context) {
                Ok(batch) => batch,
                Err(error) => {
                    self.failed = Some(error.context());
                    return Err(error);
                }
            };
        while batch.len() < maximum_records {
            if !self.next_record_into_batch(&mut batch)? {
                break;
            }
        }
        Ok(batch)
    }

    fn next_record_into_batch(
        &mut self,
        batch: &mut FastqRecordBatchData,
    ) -> Result<bool, TextRecordError> {
        if let Some(context) = self.failed {
            return Err(context.terminal_error());
        }
        match self.next_record_into_batch_inner(batch) {
            Ok(present) => Ok(present),
            Err(error) => {
                self.failed = Some(error.context());
                Err(error)
            }
        }
    }

    #[allow(clippy::too_many_lines)]
    fn next_record_into_batch_inner(
        &mut self,
        batch: &mut FastqRecordBatchData,
    ) -> Result<bool, TextRecordError> {
        let ordinal = RecordOrdinal::new(self.records_emitted);
        let initial_context = ErrorContext {
            format: TextRecordFormat::Fastq,
            side: self.side,
            ordinal,
            line: self.lines.next_line_number(),
            field: RecordField::Header,
        };
        let header_buffer = core::mem::take(&mut self.reusable_lines[0]);
        let Some(header) = self
            .lines
            .next_line_reusing(self.limits.max_line_bytes, header_buffer)
            .map_err(|error| line_error(initial_context, error))?
        else {
            return Ok(false);
        };
        let header_context = ErrorContext {
            line: header.number,
            ..initial_context
        };
        let observed_records = checked_add(
            self.records_emitted,
            1,
            TextRecordResource::Records,
            ErrorContext {
                field: RecordField::Record,
                ..header_context
            },
        )?;
        check_limit(
            observed_records,
            self.limits.max_records,
            TextRecordResource::Records,
            ErrorContext {
                field: RecordField::Record,
                ..header_context
            },
        )?;
        let (name, _) = validate_header(&header, b'@', self.limits, header_context)?;
        let name_len = name.len();
        let name_start = batch.names.len();
        batch.names.try_reserve(name_len).map_err(|_| {
            ErrorContext {
                field: RecordField::Name,
                ..header_context
            }
            .error(TextRecordErrorKind::AllocationFailed {
                allocation: TextRecordAllocation::Name,
                additional: storage_len(name_len),
            })
        })?;
        batch.names.extend_from_slice(name);

        let sequence_buffer = core::mem::take(&mut self.reusable_lines[1]);
        let sequence =
            self.read_required_line_reusing(ordinal, RecordField::Sequence, sequence_buffer)?;
        if sequence.bytes.is_empty() {
            return Err(ErrorContext {
                line: sequence.number,
                field: RecordField::Sequence,
                ..header_context
            }
            .error(TextRecordErrorKind::EmptySequence));
        }
        let sequence_len = sequence.bytes.len();
        let sequence_length = storage_len(sequence_len);
        let sequence_context = ErrorContext {
            line: sequence.number,
            field: RecordField::Sequence,
            ..header_context
        };
        check_limit(
            sequence_length,
            self.limits.max_bases_per_record,
            TextRecordResource::BasesPerRecord,
            sequence_context,
        )?;
        let sequence_start = batch.sequences.len();
        batch.sequences.try_reserve(sequence_len).map_err(|_| {
            sequence_context.error(TextRecordErrorKind::AllocationFailed {
                allocation: TextRecordAllocation::Sequence,
                additional: sequence_length,
            })
        })?;
        for (offset, &byte) in sequence.bytes.iter().enumerate() {
            let base = match byte {
                b'A' | b'a' => Base::A,
                b'C' | b'c' => Base::C,
                b'G' | b'g' => Base::G,
                b'T' | b't' => Base::T,
                b'N' | b'n' => Base::N,
                b'R' | b'r' | b'Y' | b'y' | b'S' | b's' | b'W' | b'w' | b'K' | b'k' | b'M'
                | b'm' | b'B' | b'b' | b'D' | b'd' | b'H' | b'h' | b'V' | b'v' => {
                    return Err(sequence_normalization_error(
                        sequence_context,
                        0,
                        NormalizationError::UnsupportedIupac {
                            byte,
                            offset: storage_len(offset),
                        },
                    ));
                }
                _ => {
                    return Err(sequence_normalization_error(
                        sequence_context,
                        0,
                        NormalizationError::InvalidBaseByte {
                            byte,
                            offset: storage_len(offset),
                        },
                    ));
                }
            };
            batch.sequences.push(base);
        }
        self.reusable_lines[1] = sequence.bytes;

        let total_bases = checked_add(
            self.total_bases,
            sequence_length,
            TextRecordResource::TotalBases,
            ErrorContext {
                field: RecordField::Record,
                ..header_context
            },
        )?;
        check_limit(
            total_bases,
            self.limits.max_total_bases,
            TextRecordResource::TotalBases,
            ErrorContext {
                field: RecordField::Record,
                ..header_context
            },
        )?;

        let plus_buffer = core::mem::take(&mut self.reusable_lines[2]);
        let plus = self.read_required_line_reusing(ordinal, RecordField::Plus, plus_buffer)?;
        let plus_context = ErrorContext {
            line: plus.number,
            field: RecordField::Plus,
            ..header_context
        };
        if plus.bytes.first().copied() != Some(b'+') {
            return Err(plus_context.error(TextRecordErrorKind::InvalidMarker {
                expected: b'+',
                found: plus.bytes.first().copied(),
            }));
        }
        if plus.bytes.len() > 1 && plus.bytes[1..] != header.bytes[1..] {
            return Err(plus_context.error(TextRecordErrorKind::PlusHeaderMismatch));
        }
        self.reusable_lines[2] = plus.bytes;

        let quality_buffer = core::mem::take(&mut self.reusable_lines[3]);
        let quality =
            self.read_required_line_reusing(ordinal, RecordField::Quality, quality_buffer)?;
        let quality_context = ErrorContext {
            line: quality.number,
            field: RecordField::Quality,
            ..header_context
        };
        check_limit(
            storage_len(quality.bytes.len()),
            self.limits.max_quality_bytes,
            TextRecordResource::QualityBytes,
            quality_context,
        )?;
        if quality.bytes.len() != sequence_len {
            return Err(
                quality_context.error(TextRecordErrorKind::QualityLengthMismatch {
                    sequence: sequence_length,
                    quality: storage_len(quality.bytes.len()),
                }),
            );
        }
        for (offset, &byte) in quality.bytes.iter().enumerate() {
            if !(33..=126).contains(&byte) {
                return Err(
                    quality_context.error(TextRecordErrorKind::InvalidQualityByte {
                        byte,
                        column: storage_len(offset).saturating_add(1),
                    }),
                );
            }
        }
        let quality_start = batch.qualities.len();
        batch.qualities.try_reserve(sequence_len).map_err(|_| {
            quality_context.error(TextRecordErrorKind::AllocationFailed {
                allocation: TextRecordAllocation::Quality,
                additional: sequence_length,
            })
        })?;
        batch.qualities.extend_from_slice(&quality.bytes);
        self.reusable_lines[3] = quality.bytes;
        self.reusable_lines[0] = header.bytes;
        batch.records.push(FastqRecordOffsets {
            ordinal,
            name_start,
            name_len,
            sequence_start,
            sequence_len,
            quality_start,
        });
        self.records_emitted = observed_records;
        self.total_bases = total_bases;
        Ok(true)
    }

    #[allow(clippy::too_many_lines)]
    fn next_record_inner(&mut self) -> Result<Option<FastqRecordData>, TextRecordError> {
        let ordinal = RecordOrdinal::new(self.records_emitted);
        let header_context = ErrorContext {
            format: TextRecordFormat::Fastq,
            side: self.side,
            ordinal,
            line: self.lines.next_line_number(),
            field: RecordField::Header,
        };
        let header_buffer = core::mem::take(&mut self.reusable_lines[0]);
        let Some(header) = self
            .lines
            .next_line_reusing(self.limits.max_line_bytes, header_buffer)
            .map_err(|error| line_error(header_context, error))?
        else {
            return Ok(None);
        };
        let header_context = ErrorContext {
            line: header.number,
            ..header_context
        };
        let observed_records = checked_add(
            self.records_emitted,
            1,
            TextRecordResource::Records,
            ErrorContext {
                field: RecordField::Record,
                ..header_context
            },
        )?;
        check_limit(
            observed_records,
            self.limits.max_records,
            TextRecordResource::Records,
            ErrorContext {
                field: RecordField::Record,
                ..header_context
            },
        )?;
        let name = parse_header(&header, b'@', self.limits, header_context)?;

        let sequence_buffer = core::mem::take(&mut self.reusable_lines[1]);
        let sequence =
            self.read_required_line_reusing(ordinal, RecordField::Sequence, sequence_buffer)?;
        if sequence.bytes.is_empty() {
            return Err(ErrorContext {
                line: sequence.number,
                field: RecordField::Sequence,
                ..header_context
            }
            .error(TextRecordErrorKind::EmptySequence));
        }
        let sequence_length = storage_len(sequence.bytes.len());
        let sequence_context = ErrorContext {
            line: sequence.number,
            field: RecordField::Sequence,
            ..header_context
        };
        check_limit(
            sequence_length,
            self.limits.max_bases_per_record,
            TextRecordResource::BasesPerRecord,
            sequence_context,
        )?;
        let normalized = normalize_sequence_line(&sequence, 0, sequence_context)?;
        self.reusable_lines[1] = sequence.bytes;

        let total_bases = checked_add(
            self.total_bases,
            sequence_length,
            TextRecordResource::TotalBases,
            ErrorContext {
                field: RecordField::Record,
                ..header_context
            },
        )?;
        check_limit(
            total_bases,
            self.limits.max_total_bases,
            TextRecordResource::TotalBases,
            ErrorContext {
                field: RecordField::Record,
                ..header_context
            },
        )?;

        let plus_buffer = core::mem::take(&mut self.reusable_lines[2]);
        let plus = self.read_required_line_reusing(ordinal, RecordField::Plus, plus_buffer)?;
        let plus_context = ErrorContext {
            line: plus.number,
            field: RecordField::Plus,
            ..header_context
        };
        if plus.bytes.first().copied() != Some(b'+') {
            return Err(plus_context.error(TextRecordErrorKind::InvalidMarker {
                expected: b'+',
                found: plus.bytes.first().copied(),
            }));
        }
        if plus.bytes.len() > 1 && plus.bytes[1..] != header.bytes[1..] {
            return Err(plus_context.error(TextRecordErrorKind::PlusHeaderMismatch));
        }
        self.reusable_lines[2] = plus.bytes;

        let quality = self.read_required_line(ordinal, RecordField::Quality)?;
        let quality_context = ErrorContext {
            line: quality.number,
            field: RecordField::Quality,
            ..header_context
        };
        let quality_length = storage_len(quality.bytes.len());
        check_limit(
            quality_length,
            self.limits.max_quality_bytes,
            TextRecordResource::QualityBytes,
            quality_context,
        )?;
        if quality_length != sequence_length {
            return Err(
                quality_context.error(TextRecordErrorKind::QualityLengthMismatch {
                    sequence: sequence_length,
                    quality: quality_length,
                }),
            );
        }
        for (offset, &byte) in quality.bytes.iter().enumerate() {
            if !(33..=126).contains(&byte) {
                return Err(
                    quality_context.error(TextRecordErrorKind::InvalidQualityByte {
                        byte,
                        column: storage_len(offset).saturating_add(1),
                    }),
                );
            }
        }

        self.records_emitted = observed_records;
        self.total_bases = total_bases;
        self.reusable_lines[0] = header.bytes;
        Ok(Some(FastqRecordData {
            ordinal,
            name,
            sequence: normalized,
            quality: quality.bytes.into_boxed_slice(),
        }))
    }

    fn read_required_line(
        &mut self,
        ordinal: RecordOrdinal,
        field: RecordField,
    ) -> Result<PhysicalLine, TextRecordError> {
        let context = ErrorContext {
            format: TextRecordFormat::Fastq,
            side: self.side,
            ordinal,
            line: self.lines.next_line_number(),
            field,
        };
        self.lines
            .next_line(self.limits.max_line_bytes)
            .map_err(|error| line_error(context, error))?
            .ok_or_else(|| context.error(TextRecordErrorKind::UnexpectedEof))
    }

    fn read_required_line_reusing(
        &mut self,
        ordinal: RecordOrdinal,
        field: RecordField,
        buffer: Vec<u8>,
    ) -> Result<PhysicalLine, TextRecordError> {
        let context = ErrorContext {
            format: TextRecordFormat::Fastq,
            side: self.side,
            ordinal,
            line: self.lines.next_line_number(),
            field,
        };
        self.lines
            .next_line_reusing(self.limits.max_line_bytes, buffer)
            .map_err(|error| line_error(context, error))?
            .ok_or_else(|| context.error(TextRecordErrorKind::UnexpectedEof))
    }

    fn next_line_number(&self) -> u64 {
        self.lines.next_line_number()
    }
}

/// Streaming synchronizer for two strict FASTQ sources.
#[derive(Debug)]
pub(crate) struct PairedFastqReaderCore<R1, R2> {
    first: FastqReaderCore<R1>,
    second: FastqReaderCore<R2>,
    pairs_emitted: u64,
    failed: Option<ErrorContext>,
}

impl<R1: BufRead, R2: BufRead> PairedFastqReaderCore<R1, R2> {
    /// Creates a paired parser using the same explicit limits for each source.
    #[must_use]
    pub const fn new(first: R1, second: R2, limits: TextRecordLimits) -> Self {
        Self {
            first: FastqReaderCore::with_side(first, limits, PairSourceSide::First),
            second: FastqReaderCore::with_side(second, limits, PairSourceSide::Second),
            pairs_emitted: 0,
            failed: None,
        }
    }

    /// Returns both underlying buffered sources in first/second order.
    #[must_use]
    pub const fn get_ref(&self) -> (&R1, &R2) {
        (self.first.get_ref(), self.second.get_ref())
    }

    /// Recovers both sources in first/second order and discards parser state.
    ///
    /// This is intended for source lifecycle operations such as explicit close,
    /// not for resuming either source through a second paired parser.
    #[must_use]
    pub fn into_inner(self) -> (R1, R2) {
        (self.first.into_inner(), self.second.into_inner())
    }

    /// Parses and synchronizes the next read pair.
    ///
    /// # Errors
    ///
    /// Returns source parsing failures, unequal record counts, or incompatible
    /// name tokens. After any error, later calls return `TerminalState` without
    /// reading more input.
    pub fn next_pair(&mut self) -> Result<Option<PairedFastqRecordData>, TextRecordError> {
        if let Some(context) = self.failed {
            return Err(context.terminal_error());
        }
        match self.next_pair_inner() {
            Ok(pair) => Ok(pair),
            Err(error) => {
                self.failed = Some(error.context());
                Err(error)
            }
        }
    }

    fn next_pair_inner(&mut self) -> Result<Option<PairedFastqRecordData>, TextRecordError> {
        let ordinal = RecordOrdinal::new(self.pairs_emitted);
        let first = self.first.next_record()?;
        let second = self.second.next_record()?;
        let (first, second) = match (first, second) {
            (Some(first), Some(second)) => (first, second),
            (None, None) => return Ok(None),
            (Some(_), None) => {
                return Err(ErrorContext {
                    format: TextRecordFormat::Fastq,
                    side: PairSourceSide::Second,
                    ordinal,
                    line: self.second.next_line_number(),
                    field: RecordField::Pair,
                }
                .error(TextRecordErrorKind::PairCountMismatch {
                    missing: PairSourceSide::Second,
                }));
            }
            (None, Some(_)) => {
                return Err(ErrorContext {
                    format: TextRecordFormat::Fastq,
                    side: PairSourceSide::First,
                    ordinal,
                    line: self.first.next_line_number(),
                    field: RecordField::Pair,
                }
                .error(TextRecordErrorKind::PairCountMismatch {
                    missing: PairSourceSide::First,
                }));
            }
        };

        if !paired_names_compatible(first.name.name(), second.name.name()) {
            return Err(ErrorContext {
                format: TextRecordFormat::Fastq,
                side: PairSourceSide::Single,
                ordinal,
                line: self.first.next_line_number().saturating_sub(4),
                field: RecordField::Pair,
            }
            .error(TextRecordErrorKind::PairNameMismatch {
                first: first.name.name,
                second: second.name.name,
            }));
        }

        self.pairs_emitted = self.pairs_emitted.checked_add(1).ok_or_else(|| {
            ErrorContext {
                format: TextRecordFormat::Fastq,
                side: PairSourceSide::Single,
                ordinal,
                line: self.first.next_line_number(),
                field: RecordField::Pair,
            }
            .error(TextRecordErrorKind::ArithmeticOverflow {
                resource: TextRecordResource::Records,
                current: self.pairs_emitted,
                increment: 1,
            })
        })?;
        Ok(Some(PairedFastqRecordData {
            ordinal,
            first,
            second,
        }))
    }
}

fn paired_names_compatible(first: &[u8], second: &[u8]) -> bool {
    if first == second {
        return true;
    }
    let Some(first_stem) = first.strip_suffix(b"/1") else {
        return false;
    };
    let Some(second_stem) = second.strip_suffix(b"/2") else {
        return false;
    };
    first_stem == second_stem
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use super::{
        BoundedLineReader, FastqReaderCore, LineReadError, PairedFastqReaderCore,
        TextRecordErrorKind, TextRecordLimits, TextRecordResource,
    };

    #[test]
    fn physical_line_counter_overflow_is_reported() {
        let mut reader = BoundedLineReader::new(Cursor::new(b"A\n"));
        reader.lines_read = u64::MAX;
        assert!(matches!(
            reader.next_line(1),
            Err(LineReadError::ArithmeticOverflow {
                current: u64::MAX,
                increment: 1
            })
        ));
    }

    #[test]
    fn record_and_total_base_counter_overflows_are_reported() {
        let mut records =
            FastqReaderCore::new(Cursor::new(b"@r\nA\n+\n!\n"), TextRecordLimits::MAX);
        records.records_emitted = u64::MAX;
        assert!(matches!(
            records.next_record().expect_err("record overflow").kind(),
            TextRecordErrorKind::ArithmeticOverflow {
                resource: TextRecordResource::Records,
                current: u64::MAX,
                increment: 1,
            }
        ));

        let mut bases = FastqReaderCore::new(Cursor::new(b"@r\nA\n+\n!\n"), TextRecordLimits::MAX);
        bases.total_bases = u64::MAX;
        assert!(matches!(
            bases.next_record().expect_err("base overflow").kind(),
            TextRecordErrorKind::ArithmeticOverflow {
                resource: TextRecordResource::TotalBases,
                current: u64::MAX,
                increment: 1,
            }
        ));
    }

    #[test]
    fn paired_counter_overflow_is_reported_after_two_valid_records() {
        let mut reader = PairedFastqReaderCore::new(
            Cursor::new(b"@r\nA\n+\n!\n"),
            Cursor::new(b"@r\nT\n+\n!\n"),
            TextRecordLimits::MAX,
        );
        reader.pairs_emitted = u64::MAX;
        assert!(matches!(
            reader.next_pair().expect_err("pair overflow").kind(),
            TextRecordErrorKind::ArithmeticOverflow {
                resource: TextRecordResource::Records,
                current: u64::MAX,
                increment: 1,
            }
        ));
    }

    #[test]
    fn record_batch_retains_contiguous_record_views() {
        let input = b"@pair/1 description\nAcgTN\n+\n!#$%&\n@next\nTTAA\n+next\nIIII\n";
        let mut reader = FastqReaderCore::new(Cursor::new(input), TextRecordLimits::MAX);
        let batch = reader.next_batch(8).expect("validated batch");
        assert_eq!(batch.len(), 2);
        let first = batch.get(0).expect("first record");
        assert_eq!(first.ordinal().get(), 0);
        assert_eq!(first.name(), b"pair/1");
        assert_eq!(
            first
                .sequence()
                .iter()
                .map(|base| base.as_ascii())
                .collect::<Vec<_>>(),
            b"ACGTN"
        );
        assert_eq!(first.quality(), b"!#$%&");
        let second = batch.get(1).expect("second record");
        assert_eq!(second.ordinal().get(), 1);
        assert_eq!(second.name(), b"next");
        assert_eq!(second.quality(), b"IIII");
        assert!(reader.next_batch(8).expect("EOF batch").is_empty());
    }
}
