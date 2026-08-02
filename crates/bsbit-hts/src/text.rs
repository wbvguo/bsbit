//! Shared bounded mechanics for FASTA and FASTQ record parsing.

use core::fmt;
use std::error::Error;
use std::io::{self, BufRead, Write};

use bsbit_core::alphabet::Base;
use bsbit_core::sequence::{NormalizationError, NormalizedSequence};

/// Zero-based identity of one record in its input source.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct RecordOrdinal(u64);

impl RecordOrdinal {
    /// Creates an ordinal from its zero-based value.
    #[must_use]
    pub const fn new(value: u64) -> Self {
        Self(value)
    }

    /// Returns the zero-based value.
    #[must_use]
    pub const fn get(self) -> u64 {
        self.0
    }
}

/// Exact identifier and optional description retained from a record header.
#[derive(Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct RecordName {
    pub(crate) name: Box<[u8]>,
    pub(crate) description: Box<[u8]>,
}

impl RecordName {
    /// Returns the exact nonempty identifier bytes, excluding the format marker.
    #[must_use]
    pub fn name(&self) -> &[u8] {
        &self.name
    }

    /// Returns the description bytes after leading header separators were removed.
    #[must_use]
    pub fn description(&self) -> &[u8] {
        &self.description
    }

    /// Returns whether the record had no description bytes.
    #[must_use]
    pub fn description_is_empty(&self) -> bool {
        self.description.is_empty()
    }
}

pub(crate) fn write_record_name<W: Write>(writer: &mut W, name: &RecordName) -> io::Result<()> {
    writer.write_all(name.name())?;
    if !name.description_is_empty() {
        writer.write_all(b" ")?;
        writer.write_all(name.description())?;
    }
    Ok(())
}

pub(crate) fn write_sequence<W: Write>(
    writer: &mut W,
    sequence: &NormalizedSequence,
) -> io::Result<()> {
    const CHUNK_BASES: usize = 4_096;
    let mut chunk = [0_u8; CHUNK_BASES];
    for bases in sequence.bases().chunks(CHUNK_BASES) {
        for (destination, base) in chunk.iter_mut().zip(bases) {
            *destination = base.as_ascii();
        }
        writer.write_all(&chunk[..bases.len()])?;
    }
    Ok(())
}

/// Explicit caps applied to one text-record source.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[allow(clippy::struct_field_names)]
pub struct TextRecordLimits {
    pub(crate) max_line_bytes: u64,
    pub(crate) max_records: u64,
    pub(crate) max_name_bytes: u64,
    pub(crate) max_description_bytes: u64,
    pub(crate) max_bases_per_record: u64,
    pub(crate) max_total_bases: u64,
    pub(crate) max_quality_bytes: u64,
}

impl TextRecordLimits {
    /// Limits admitting every representable logical count.
    ///
    /// Selecting this constant is explicit; parsers do not choose it by default.
    pub const MAX: Self = Self {
        max_line_bytes: u64::MAX,
        max_records: u64::MAX,
        max_name_bytes: u64::MAX,
        max_description_bytes: u64::MAX,
        max_bases_per_record: u64::MAX,
        max_total_bases: u64::MAX,
        max_quality_bytes: u64::MAX,
    };

    /// Creates a complete explicit limit set.
    #[must_use]
    #[allow(clippy::too_many_arguments)]
    pub const fn new(
        max_line_bytes: u64,
        max_records: u64,
        max_name_bytes: u64,
        max_description_bytes: u64,
        max_bases_per_record: u64,
        max_total_bases: u64,
        max_quality_bytes: u64,
    ) -> Self {
        Self {
            max_line_bytes,
            max_records,
            max_name_bytes,
            max_description_bytes,
            max_bases_per_record,
            max_total_bases,
            max_quality_bytes,
        }
    }

    /// Returns the maximum physical-line content bytes, excluding LF or CRLF.
    #[must_use]
    pub const fn max_line_bytes(self) -> u64 {
        self.max_line_bytes
    }

    /// Returns the maximum number of emitted records.
    #[must_use]
    pub const fn max_records(self) -> u64 {
        self.max_records
    }

    /// Returns the maximum name bytes per record.
    #[must_use]
    pub const fn max_name_bytes(self) -> u64 {
        self.max_name_bytes
    }

    /// Returns the maximum description bytes per record.
    #[must_use]
    pub const fn max_description_bytes(self) -> u64 {
        self.max_description_bytes
    }

    /// Returns the maximum normalized bases per record.
    #[must_use]
    pub const fn max_bases_per_record(self) -> u64 {
        self.max_bases_per_record
    }

    /// Returns the maximum total emitted bases per source.
    #[must_use]
    pub const fn max_total_bases(self) -> u64 {
        self.max_total_bases
    }

    /// Returns the maximum quality bytes per FASTQ record.
    #[must_use]
    pub const fn max_quality_bytes(self) -> u64 {
        self.max_quality_bytes
    }
}

/// Text record syntax associated with an error.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum TextRecordFormat {
    /// FASTA input.
    Fasta,
    /// Strict four-line FASTQ input.
    Fastq,
}

/// Input source associated with an error.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum PairSourceSide {
    /// A non-paired parser source.
    Single,
    /// Source one of a paired FASTQ parser.
    First,
    /// Source two of a paired FASTQ parser.
    Second,
}

/// Logical record field associated with an error.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum RecordField {
    /// Whole-record boundary or record count.
    Record,
    /// Header marker or complete header.
    Header,
    /// Header identifier token.
    Name,
    /// Header description.
    Description,
    /// DNA sequence.
    Sequence,
    /// FASTQ plus line.
    Plus,
    /// FASTQ quality line.
    Quality,
    /// Paired-source synchronization or identity.
    Pair,
}

/// Logical resource controlled by a configured limit.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum TextRecordResource {
    /// Physical line content bytes.
    LineBytes,
    /// Emitted records.
    Records,
    /// Header name bytes.
    NameBytes,
    /// Header description bytes.
    DescriptionBytes,
    /// Normalized bases in one record.
    BasesPerRecord,
    /// Total emitted normalized bases.
    TotalBases,
    /// FASTQ quality bytes in one record.
    QualityBytes,
}

/// Fallible allocation site reported by a parser.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum TextRecordAllocation {
    /// Current bounded physical line.
    Line,
    /// Contiguous record-batch metadata.
    Record,
    /// Header name.
    Name,
    /// Header description.
    Description,
    /// Aggregated normalized sequence.
    Sequence,
    /// Quality bytes.
    Quality,
}

/// Specific cause of a [`TextRecordError`].
#[derive(Debug)]
pub enum TextRecordErrorKind {
    /// The underlying reader failed.
    Io {
        /// Original I/O error.
        source: io::Error,
    },
    /// A logical configured cap was exceeded.
    LimitExceeded {
        /// Controlled resource.
        resource: TextRecordResource,
        /// First observed count known to exceed the cap.
        observed: u64,
        /// Configured cap.
        limit: u64,
    },
    /// A checked logical count could not be represented.
    ArithmeticOverflow {
        /// Resource being counted.
        resource: TextRecordResource,
        /// Value before addition.
        current: u64,
        /// Attempted increment.
        increment: u64,
    },
    /// A fallible storage reservation failed.
    AllocationFailed {
        /// Allocation site.
        allocation: TextRecordAllocation,
        /// Requested additional bytes or elements.
        additional: u64,
    },
    /// Input ended before a required field was read.
    UnexpectedEof,
    /// A required byte-column-one marker was absent.
    InvalidMarker {
        /// Required marker.
        expected: u8,
        /// Actual first byte, or `None` for an empty line.
        found: Option<u8>,
    },
    /// A header has no identifier token.
    EmptyName,
    /// A header contains a forbidden control or non-ASCII byte.
    InvalidHeaderByte {
        /// Offending byte.
        byte: u8,
        /// One-based byte column.
        column: u64,
    },
    /// A record contains no sequence bases.
    EmptySequence,
    /// Core DNA normalization rejected a byte.
    InvalidSequence {
        /// One-based byte column on the physical sequence line.
        column: u64,
        /// Zero-based offset in the concatenated record sequence.
        record_offset: u64,
        /// Core normalization cause with line-local offset.
        source: NormalizationError,
    },
    /// A FASTQ plus-line suffix did not exactly equal the header tail.
    PlusHeaderMismatch,
    /// FASTQ sequence and quality lengths differ.
    QualityLengthMismatch {
        /// Normalized sequence length.
        sequence: u64,
        /// Quality byte length.
        quality: u64,
    },
    /// A quality byte lies outside printable Phred+33 ASCII.
    InvalidQualityByte {
        /// Offending byte.
        byte: u8,
        /// One-based byte column.
        column: u64,
    },
    /// One paired source ended before the other.
    PairCountMismatch {
        /// Source that had no record at this ordinal.
        missing: PairSourceSide,
    },
    /// Paired record name tokens were incompatible.
    PairNameMismatch {
        /// Exact source-one name token.
        first: Box<[u8]>,
        /// Exact source-two name token.
        second: Box<[u8]>,
    },
    /// A previous parser error made the input boundary unusable.
    TerminalState,
}

/// Structured deterministic failure from text-record parsing.
#[derive(Debug)]
pub struct TextRecordError {
    format: TextRecordFormat,
    side: PairSourceSide,
    ordinal: RecordOrdinal,
    line: u64,
    field: RecordField,
    kind: TextRecordErrorKind,
}

impl TextRecordError {
    /// Returns the record syntax being parsed.
    #[must_use]
    pub const fn format(&self) -> TextRecordFormat {
        self.format
    }

    /// Returns the single or paired source side.
    #[must_use]
    pub const fn side(&self) -> PairSourceSide {
        self.side
    }

    /// Returns the zero-based record ordinal being parsed.
    #[must_use]
    pub const fn ordinal(&self) -> RecordOrdinal {
        self.ordinal
    }

    /// Returns the one-based physical line associated with the failure.
    #[must_use]
    pub const fn line(&self) -> u64 {
        self.line
    }

    /// Returns the logical field associated with the failure.
    #[must_use]
    pub const fn field(&self) -> RecordField {
        self.field
    }

    /// Returns the specific cause.
    #[must_use]
    pub const fn kind(&self) -> &TextRecordErrorKind {
        &self.kind
    }

    pub(crate) const fn context(&self) -> ErrorContext {
        ErrorContext {
            format: self.format,
            side: self.side,
            ordinal: self.ordinal,
            line: self.line,
            field: self.field,
        }
    }
}

impl fmt::Display for TextRecordError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "{:?} {:?} record {} line {} {:?}: ",
            self.format,
            self.side,
            self.ordinal.get(),
            self.line,
            self.field
        )?;
        match &self.kind {
            TextRecordErrorKind::Io { source } => source.fmt(formatter),
            TextRecordErrorKind::LimitExceeded {
                resource,
                observed,
                limit,
            } => write!(
                formatter,
                "{resource:?} count {observed} exceeds configured limit {limit}"
            ),
            TextRecordErrorKind::ArithmeticOverflow {
                resource,
                current,
                increment,
            } => write!(
                formatter,
                "{resource:?} count overflow while adding {increment} to {current}"
            ),
            TextRecordErrorKind::AllocationFailed {
                allocation,
                additional,
            } => write!(
                formatter,
                "failed to reserve {additional} additional units for {allocation:?}"
            ),
            TextRecordErrorKind::UnexpectedEof => formatter.write_str("unexpected end of input"),
            TextRecordErrorKind::InvalidMarker { expected, found } => write!(
                formatter,
                "expected marker '{}' but found {found:?}",
                char::from(*expected)
            ),
            TextRecordErrorKind::EmptyName => formatter.write_str("empty record name"),
            TextRecordErrorKind::InvalidHeaderByte { byte, column } => write!(
                formatter,
                "invalid header byte 0x{byte:02X} at column {column}"
            ),
            TextRecordErrorKind::EmptySequence => formatter.write_str("empty sequence"),
            TextRecordErrorKind::InvalidSequence {
                column,
                record_offset,
                source,
            } => write!(
                formatter,
                "{source} at physical column {column}, record offset {record_offset}"
            ),
            TextRecordErrorKind::PlusHeaderMismatch => {
                formatter.write_str("plus-line header suffix differs from the record header")
            }
            TextRecordErrorKind::QualityLengthMismatch { sequence, quality } => write!(
                formatter,
                "quality length {quality} differs from sequence length {sequence}"
            ),
            TextRecordErrorKind::InvalidQualityByte { byte, column } => write!(
                formatter,
                "invalid quality byte 0x{byte:02X} at column {column}"
            ),
            TextRecordErrorKind::PairCountMismatch { missing } => {
                write!(formatter, "paired source {missing:?} ended early")
            }
            TextRecordErrorKind::PairNameMismatch { first, second } => {
                write!(
                    formatter,
                    "paired names {first:?} and {second:?} are incompatible"
                )
            }
            TextRecordErrorKind::TerminalState => {
                formatter.write_str("parser is in a failed terminal state")
            }
        }
    }
}

impl Error for TextRecordError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match &self.kind {
            TextRecordErrorKind::Io { source } => Some(source),
            TextRecordErrorKind::InvalidSequence { source, .. } => Some(source),
            _ => None,
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct ErrorContext {
    pub(crate) format: TextRecordFormat,
    pub(crate) side: PairSourceSide,
    pub(crate) ordinal: RecordOrdinal,
    pub(crate) line: u64,
    pub(crate) field: RecordField,
}

impl ErrorContext {
    pub(crate) const fn error(self, kind: TextRecordErrorKind) -> TextRecordError {
        TextRecordError {
            format: self.format,
            side: self.side,
            ordinal: self.ordinal,
            line: self.line,
            field: self.field,
            kind,
        }
    }

    pub(crate) const fn terminal_error(self) -> TextRecordError {
        self.error(TextRecordErrorKind::TerminalState)
    }
}

#[derive(Debug)]
pub(crate) struct PhysicalLine {
    pub(crate) number: u64,
    pub(crate) bytes: Vec<u8>,
}

#[derive(Debug)]
pub(crate) enum LineReadError {
    Io(io::Error),
    LimitExceeded { observed: u64, limit: u64 },
    ArithmeticOverflow { current: u64, increment: u64 },
    AllocationFailed { additional: u64 },
}

#[derive(Debug)]
pub(crate) struct BoundedLineReader<R> {
    inner: R,
    pub(crate) lines_read: u64,
    eof: bool,
}

impl<R: BufRead> BoundedLineReader<R> {
    pub(crate) const fn new(inner: R) -> Self {
        Self {
            inner,
            lines_read: 0,
            eof: false,
        }
    }

    pub(crate) fn next_line(&mut self, limit: u64) -> Result<Option<PhysicalLine>, LineReadError> {
        self.next_line_reusing(limit, Vec::new())
    }

    pub(crate) fn next_line_reusing(
        &mut self,
        limit: u64,
        mut bytes: Vec<u8>,
    ) -> Result<Option<PhysicalLine>, LineReadError> {
        if self.eof {
            return Ok(None);
        }

        let number = self
            .lines_read
            .checked_add(1)
            .ok_or(LineReadError::ArithmeticOverflow {
                current: self.lines_read,
                increment: 1,
            })?;
        bytes.clear();
        let mut pending_carriage_return = false;
        let mut observed_any = false;

        loop {
            let (consumed, ended_line) = {
                let buffer = self.inner.fill_buf().map_err(LineReadError::Io)?;
                if buffer.is_empty() {
                    self.eof = true;
                    if pending_carriage_return {
                        append_line_bytes(&mut bytes, b"\r", limit)?;
                    }
                    if !observed_any {
                        return Ok(None);
                    }
                    self.lines_read = number;
                    return Ok(Some(PhysicalLine { number, bytes }));
                }

                let newline = buffer.iter().position(|&byte| byte == b'\n');
                let consumed = newline.map_or(buffer.len(), |position| position + 1);
                let segment_end = newline.unwrap_or(buffer.len());
                let segment = &buffer[..segment_end];
                observed_any = true;

                if pending_carriage_return && !segment.is_empty() {
                    append_line_bytes(&mut bytes, b"\r", limit)?;
                    pending_carriage_return = false;
                }

                if let Some((&last, prefix)) = segment.split_last() {
                    if last == b'\r' {
                        append_line_bytes(&mut bytes, prefix, limit)?;
                        pending_carriage_return = true;
                    } else {
                        append_line_bytes(&mut bytes, segment, limit)?;
                    }
                }

                if newline.is_some() {
                    pending_carriage_return = false;
                }
                (consumed, newline.is_some())
            };

            self.inner.consume(consumed);
            if ended_line {
                self.lines_read = number;
                return Ok(Some(PhysicalLine { number, bytes }));
            }
        }
    }

    pub(crate) fn next_line_number(&self) -> u64 {
        self.lines_read.saturating_add(1)
    }

    pub(crate) const fn get_ref(&self) -> &R {
        &self.inner
    }

    pub(crate) fn into_inner(self) -> R {
        self.inner
    }
}

fn append_line_bytes(
    destination: &mut Vec<u8>,
    bytes: &[u8],
    limit: u64,
) -> Result<(), LineReadError> {
    let current = storage_len(destination.len());
    let increment = storage_len(bytes.len());
    let observed = current
        .checked_add(increment)
        .ok_or(LineReadError::ArithmeticOverflow { current, increment })?;
    if observed > limit {
        return Err(LineReadError::LimitExceeded { observed, limit });
    }
    destination
        .try_reserve(bytes.len())
        .map_err(|_| LineReadError::AllocationFailed {
            additional: increment,
        })?;
    destination.extend_from_slice(bytes);
    Ok(())
}

pub(crate) fn storage_len(length: usize) -> u64 {
    u64::try_from(length).expect("supported pointer widths fit in u64")
}

pub(crate) fn line_error(context: ErrorContext, error: LineReadError) -> TextRecordError {
    let kind = match error {
        LineReadError::Io(source) => TextRecordErrorKind::Io { source },
        LineReadError::LimitExceeded { observed, limit } => TextRecordErrorKind::LimitExceeded {
            resource: TextRecordResource::LineBytes,
            observed,
            limit,
        },
        LineReadError::ArithmeticOverflow { current, increment } => {
            TextRecordErrorKind::ArithmeticOverflow {
                resource: TextRecordResource::LineBytes,
                current,
                increment,
            }
        }
        LineReadError::AllocationFailed { additional } => TextRecordErrorKind::AllocationFailed {
            allocation: TextRecordAllocation::Line,
            additional,
        },
    };
    context.error(kind)
}

pub(crate) fn allocate_bytes(
    source: &[u8],
    allocation: TextRecordAllocation,
    context: ErrorContext,
) -> Result<Box<[u8]>, TextRecordError> {
    let mut destination = Vec::new();
    destination.try_reserve_exact(source.len()).map_err(|_| {
        context.error(TextRecordErrorKind::AllocationFailed {
            allocation,
            additional: storage_len(source.len()),
        })
    })?;
    destination.extend_from_slice(source);
    Ok(destination.into_boxed_slice())
}

pub(crate) fn check_limit(
    observed: u64,
    limit: u64,
    resource: TextRecordResource,
    context: ErrorContext,
) -> Result<(), TextRecordError> {
    if observed > limit {
        return Err(context.error(TextRecordErrorKind::LimitExceeded {
            resource,
            observed,
            limit,
        }));
    }
    Ok(())
}

pub(crate) fn checked_add(
    current: u64,
    increment: u64,
    resource: TextRecordResource,
    context: ErrorContext,
) -> Result<u64, TextRecordError> {
    current.checked_add(increment).ok_or_else(|| {
        context.error(TextRecordErrorKind::ArithmeticOverflow {
            resource,
            current,
            increment,
        })
    })
}

pub(crate) fn parse_header(
    line: &PhysicalLine,
    marker: u8,
    limits: TextRecordLimits,
    context: ErrorContext,
) -> Result<RecordName, TextRecordError> {
    let (name, description) = validate_header(line, marker, limits, context)?;
    Ok(RecordName {
        name: allocate_bytes(
            name,
            TextRecordAllocation::Name,
            ErrorContext {
                field: RecordField::Name,
                ..context
            },
        )?,
        description: allocate_bytes(
            description,
            TextRecordAllocation::Description,
            ErrorContext {
                field: RecordField::Description,
                ..context
            },
        )?,
    })
}

pub(crate) fn validate_header(
    line: &PhysicalLine,
    marker: u8,
    limits: TextRecordLimits,
    context: ErrorContext,
) -> Result<(&[u8], &[u8]), TextRecordError> {
    if line.bytes.first().copied() != Some(marker) {
        return Err(context.error(TextRecordErrorKind::InvalidMarker {
            expected: marker,
            found: line.bytes.first().copied(),
        }));
    }
    let tail = &line.bytes[1..];
    let separator = tail
        .iter()
        .position(|byte| matches!(*byte, b' ' | b'\t'))
        .unwrap_or(tail.len());
    let name = &tail[..separator];
    if name.is_empty() {
        return Err(ErrorContext {
            field: RecordField::Name,
            ..context
        }
        .error(TextRecordErrorKind::EmptyName));
    }

    for (offset, &byte) in name.iter().enumerate() {
        if !byte.is_ascii_graphic() || byte.is_ascii_whitespace() {
            let column = storage_len(offset).saturating_add(2);
            return Err(ErrorContext {
                field: RecordField::Name,
                ..context
            }
            .error(TextRecordErrorKind::InvalidHeaderByte { byte, column }));
        }
    }
    check_limit(
        storage_len(name.len()),
        limits.max_name_bytes,
        TextRecordResource::NameBytes,
        ErrorContext {
            field: RecordField::Name,
            ..context
        },
    )?;

    let description_start = tail[separator..]
        .iter()
        .position(|byte| !matches!(*byte, b' ' | b'\t'))
        .map_or(tail.len(), |relative| separator + relative);
    let description = &tail[description_start..];
    for (offset, &byte) in description.iter().enumerate() {
        let valid = matches!(byte, b' ' | b'\t') || byte.is_ascii_graphic();
        if !valid {
            let column = storage_len(description_start)
                .saturating_add(storage_len(offset))
                .saturating_add(2);
            return Err(ErrorContext {
                field: RecordField::Description,
                ..context
            }
            .error(TextRecordErrorKind::InvalidHeaderByte { byte, column }));
        }
    }
    check_limit(
        storage_len(description.len()),
        limits.max_description_bytes,
        TextRecordResource::DescriptionBytes,
        ErrorContext {
            field: RecordField::Description,
            ..context
        },
    )?;

    Ok((name, description))
}

pub(crate) fn normalize_sequence_line(
    line: &PhysicalLine,
    record_offset: u64,
    context: ErrorContext,
) -> Result<NormalizedSequence, TextRecordError> {
    let additional = storage_len(line.bytes.len());
    let mut bases = Vec::new();
    bases.try_reserve_exact(line.bytes.len()).map_err(|_| {
        context.error(TextRecordErrorKind::AllocationFailed {
            allocation: TextRecordAllocation::Sequence,
            additional,
        })
    })?;
    for (storage_offset, &byte) in line.bytes.iter().enumerate() {
        let offset = storage_len(storage_offset);
        let base = match byte {
            b'A' | b'a' => Base::A,
            b'C' | b'c' => Base::C,
            b'G' | b'g' => Base::G,
            b'T' | b't' => Base::T,
            b'N' | b'n' => Base::N,
            b'R' | b'r' | b'Y' | b'y' | b'S' | b's' | b'W' | b'w' | b'K' | b'k' | b'M' | b'm'
            | b'B' | b'b' | b'D' | b'd' | b'H' | b'h' | b'V' | b'v' => {
                return Err(sequence_normalization_error(
                    context,
                    record_offset,
                    NormalizationError::UnsupportedIupac { byte, offset },
                ));
            }
            _ => {
                return Err(sequence_normalization_error(
                    context,
                    record_offset,
                    NormalizationError::InvalidBaseByte { byte, offset },
                ));
            }
        };
        bases.push(base);
    }
    Ok(bases.into())
}

pub(crate) fn sequence_normalization_error(
    context: ErrorContext,
    record_offset: u64,
    source: NormalizationError,
) -> TextRecordError {
    let column = source.offset().saturating_add(1);
    let concatenated = record_offset.saturating_add(source.offset());
    context.error(TextRecordErrorKind::InvalidSequence {
        column,
        record_offset: concatenated,
        source,
    })
}

pub(crate) fn record_limit_context(
    format: TextRecordFormat,
    side: PairSourceSide,
    ordinal: RecordOrdinal,
    line: u64,
) -> ErrorContext {
    ErrorContext {
        format,
        side,
        ordinal,
        line,
        field: RecordField::Record,
    }
}
