//! Canonical SAM 1.6 text serialization over validated alignment records.

use core::fmt;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};

use bsbit_io::{PublicationError, PublishedFile, StagedFile};

use bsbit_core::bisulfite::AlignmentOrientation;
use bsbit_core::cigar::CoreCigar;

use crate::alignment_record::{
    AlignmentAuxiliaryMode, AlignmentCigarOp, AlignmentCigarRun, AlignmentRecord,
    AlignmentRecordAllocation, AlignmentRecordError, AlignmentRecordField, AlignmentRecordLimits,
    AlignmentRecordResource, BorrowedAlignmentRecord, RecordSegment, allocate_bytes_unbounded,
    append_u64, check_limit, checked_add_resource, cigar_text_length, decimal_digits, storage_len,
    validate_reference_length, validate_reference_name,
};

const PROGRAM_PREFIX: &[u8] = b"@PG\tID:bsbit\tPN:bsbit\tVN:";
const PROGRAM_VERSION: &[u8] = env!("CARGO_PKG_VERSION").as_bytes();
const PROGRAM_DESCRIPTION_PREFIX: &[u8] = b"\tDS:reference-semantic-sha256=";
const PROGRAM_MODE_PREFIX: &[u8] = b";alignment-mode=";

/// Alignment contract recorded in the canonical bsbit `@PG` header line.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum BsbitAlignmentMode {
    /// Caller-compatible directional single-end alignment.
    CallerCompatibleDirectionalSingle,
    /// Caller-compatible directional paired-end alignment.
    CallerCompatibleDirectionalPaired,
    /// Caller-compatible non-directional paired-end alignment.
    CallerCompatibleNondirectionalPaired,
}

impl BsbitAlignmentMode {
    const fn header_value(self) -> &'static [u8] {
        match self {
            Self::CallerCompatibleDirectionalSingle => b"caller-compatible-directional-single",
            Self::CallerCompatibleDirectionalPaired => b"caller-compatible-directional-paired",
            Self::CallerCompatibleNondirectionalPaired => {
                b"caller-compatible-nondirectional-paired"
            }
        }
    }

    fn from_header_value(value: &[u8]) -> Option<Self> {
        match value {
            b"caller-compatible-directional-single" => {
                Some(Self::CallerCompatibleDirectionalSingle)
            }
            b"caller-compatible-directional-paired" => {
                Some(Self::CallerCompatibleDirectionalPaired)
            }
            b"caller-compatible-nondirectional-paired" => {
                Some(Self::CallerCompatibleNondirectionalPaired)
            }
            _ => None,
        }
    }

    /// Returns whether the mapping-quality and auxiliary-tag contract is
    /// suitable for the bsbit caller after coordinate sorting and indexing.
    #[must_use]
    pub const fn is_caller_compatible(self) -> bool {
        matches!(
            self,
            Self::CallerCompatibleDirectionalSingle
                | Self::CallerCompatibleDirectionalPaired
                | Self::CallerCompatibleNondirectionalPaired
        )
    }
}

/// Exact provenance embedded in and recovered from a bsbit `@PG` record.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct BsbitProgramProvenance {
    reference_semantic_digest: [u8; 32],
    alignment_mode: BsbitAlignmentMode,
}

impl BsbitProgramProvenance {
    /// Constructs provenance for one verified reference and alignment mode.
    #[must_use]
    pub const fn new(
        reference_semantic_digest: [u8; 32],
        alignment_mode: BsbitAlignmentMode,
    ) -> Self {
        Self {
            reference_semantic_digest,
            alignment_mode,
        }
    }

    /// Returns the exact semantic reference digest.
    #[must_use]
    pub const fn reference_semantic_digest(self) -> [u8; 32] {
        self.reference_semantic_digest
    }

    /// Returns the declared alignment contract.
    #[must_use]
    pub const fn alignment_mode(self) -> BsbitAlignmentMode {
        self.alignment_mode
    }
}

/// Malformed or ambiguous bsbit provenance in a SAM/BAM header.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum BsbitProgramProvenanceError {
    /// More than one matching bsbit program record was present.
    DuplicateProgramRecord,
    /// The matching program record omitted its structured description.
    MissingDescription,
    /// The structured description did not match the versioned grammar.
    MalformedDescription,
}

impl fmt::Display for BsbitProgramProvenanceError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::DuplicateProgramRecord => "BAM header repeats `@PG ID:bsbit PN:bsbit`",
            Self::MissingDescription => "bsbit @PG record lacks reference/alignment provenance",
            Self::MalformedDescription => "bsbit @PG provenance is malformed or unsupported",
        })
    }
}

impl std::error::Error for BsbitProgramProvenanceError {}

/// One validated reference-dictionary entry for a SAM header.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SamHeaderReference {
    name: Box<[u8]>,
    length: u32,
}

impl SamHeaderReference {
    /// Validates and owns one SAM header dictionary entry.
    ///
    /// # Errors
    ///
    /// Returns reference-name, length, or allocation failures.
    pub fn new(ordinal: u64, name: &[u8], length: u64) -> Result<Self, AlignmentRecordError> {
        validate_reference_name(ordinal, name)?;
        validate_reference_length(ordinal, length)?;
        Ok(Self {
            name: allocate_bytes_unbounded(name, AlignmentRecordAllocation::ReferenceName)?,
            length: u32::try_from(length)
                .map_err(|_| AlignmentRecordError::ReferenceLengthOutOfRange { ordinal, length })?,
        })
    }

    /// Returns the exact SAM reference name.
    #[must_use]
    pub fn name(&self) -> &[u8] {
        &self.name
    }
    /// Returns the reference length.
    #[must_use]
    pub const fn length(&self) -> u32 {
        self.length
    }
}

/// Immutable canonical SAM header metadata.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SamHeader {
    references: Vec<SamHeaderReference>,
    sort_order: SamSortOrder,
    bsbit_provenance: Option<BsbitProgramProvenance>,
}

/// Declared alignment ordering for the canonical SAM/BAM header.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum SamSortOrder {
    /// Records have no declared ordering.
    #[default]
    Unsorted,
    /// Records are ordered by reference dictionary and coordinate.
    Coordinate,
    /// Records are grouped by query name.
    QueryName,
}

impl SamSortOrder {
    const fn header_prefix(self) -> &'static [u8] {
        match self {
            Self::Unsorted => b"@HD\tVN:1.6\tSO:unsorted\n",
            Self::Coordinate => b"@HD\tVN:1.6\tSO:coordinate\n",
            Self::QueryName => b"@HD\tVN:1.6\tSO:queryname\n",
        }
    }
}

impl SamHeader {
    /// Builds an ordered SAM dictionary from validated format-level entries.
    ///
    /// # Errors
    ///
    /// Returns [`AlignmentRecordError`] for configured limits, arithmetic, or
    /// allocation failures. Entry grammar and range are checked by
    /// [`SamHeaderReference::new`].
    pub fn new(
        references: Vec<SamHeaderReference>,
        limits: AlignmentRecordLimits,
    ) -> Result<Self, AlignmentRecordError> {
        check_limit(
            storage_len(references.len()),
            limits.max_header_references(),
            AlignmentRecordResource::HeaderReferences,
        )?;
        let mut name_bytes = 0_u64;
        for reference in &references {
            name_bytes = checked_add_resource(
                name_bytes,
                storage_len(reference.name().len()),
                AlignmentRecordResource::HeaderNameBytes,
            )?;
            check_limit(
                name_bytes,
                limits.max_header_name_bytes(),
                AlignmentRecordResource::HeaderNameBytes,
            )?;
        }
        let header = Self {
            references,
            sort_order: SamSortOrder::Unsorted,
            bsbit_provenance: None,
        };
        let length = header_text_length(&header)?;
        check_limit(
            length,
            limits.max_header_bytes(),
            AlignmentRecordResource::HeaderBytes,
        )?;
        Ok(header)
    }

    /// Returns ordered reference dictionary entries.
    #[must_use]
    pub fn references(&self) -> &[SamHeaderReference] {
        &self.references
    }

    /// Returns a copy declaring the ordering the writer will actually use.
    #[must_use]
    pub fn with_sort_order(mut self, sort_order: SamSortOrder) -> Self {
        self.sort_order = sort_order;
        self
    }

    /// Returns the declared alignment ordering.
    #[must_use]
    pub const fn sort_order(&self) -> SamSortOrder {
        self.sort_order
    }

    /// Returns a copy with exact reference and alignment provenance in `@PG`.
    ///
    /// # Errors
    ///
    /// Returns a header-size failure when the configured limit cannot retain
    /// the fixed structured provenance.
    pub fn with_bsbit_provenance(
        mut self,
        provenance: BsbitProgramProvenance,
        limits: AlignmentRecordLimits,
    ) -> Result<Self, AlignmentRecordError> {
        self.bsbit_provenance = Some(provenance);
        check_limit(
            header_text_length(&self)?,
            limits.max_header_bytes(),
            AlignmentRecordResource::HeaderBytes,
        )?;
        Ok(self)
    }

    /// Returns the structured bsbit program provenance, when configured.
    #[must_use]
    pub const fn bsbit_provenance(&self) -> Option<BsbitProgramProvenance> {
        self.bsbit_provenance
    }
}

/// Phase associated with a generic writer failure.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SamWritePhase {
    /// Header bytes.
    Header,
    /// One alignment-record line.
    Record,
}

/// Encoding or generic writer failure.
#[derive(Debug)]
pub enum SamWriteError {
    /// Validation, limit, arithmetic, or allocation failure before writing.
    Encode {
        /// Underlying record-boundary failure.
        source: AlignmentRecordError,
    },
    /// The caller-owned writer failed and may contain a prefix.
    Io {
        /// Failed output phase.
        phase: SamWritePhase,
        /// Underlying writer error.
        source: io::Error,
    },
}

impl fmt::Display for SamWriteError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Encode { source } => write!(formatter, "SAM encoding failed: {source}"),
            Self::Io { phase, source } => {
                write!(formatter, "SAM {phase:?} write failed: {source}")
            }
        }
    }
}

impl std::error::Error for SamWriteError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Encode { source } => Some(source),
            Self::Io { source, .. } => Some(source),
        }
    }
}

/// Returns the exact SAM FLAG for a validated record.
#[must_use]
pub fn sam_flag(record: &AlignmentRecord) -> u16 {
    let paired = !matches!(record.segment(), RecordSegment::Unpaired);
    let mapped = record.mapping();
    let mate = record.mate();
    let mut flag = 0_u16;
    if paired {
        flag |= 0x1;
        if record.is_proper_pair() {
            flag |= 0x2;
        }
        if mate.is_none() {
            flag |= 0x8;
        }
        if mate.is_some_and(|mate| matches!(mate.orientation(), AlignmentOrientation::Reverse)) {
            flag |= 0x20;
        }
        flag |= match record.segment() {
            RecordSegment::First => 0x40,
            RecordSegment::Last => 0x80,
            RecordSegment::Unpaired => 0,
        };
    }
    if mapped.is_none() {
        flag |= 0x4;
    }
    if mapped.is_some_and(|mapping| matches!(mapping.orientation(), AlignmentOrientation::Reverse))
    {
        flag |= 0x10;
    }
    flag
}

/// Encodes one complete canonical SAM alignment line.
///
/// # Errors
///
/// Returns [`AlignmentRecordError`] if the exact line size overflows, exceeds
/// the configured cap, or cannot be reserved.
pub fn sam_record_bytes(
    record: &AlignmentRecord,
    limits: AlignmentRecordLimits,
) -> Result<Vec<u8>, AlignmentRecordError> {
    if let Some(mapping) = record.mapping() {
        check_limit(
            storage_len(mapping.cigar().run_count()),
            limits.max_cigar_runs(),
            AlignmentRecordResource::CigarRuns,
        )?;
        check_limit(
            cigar_text_length(mapping.cigar())?,
            limits.max_cigar_text_bytes(),
            AlignmentRecordResource::CigarTextBytes,
        )?;
    }
    let length = record_text_length(record)?;
    check_limit(
        length,
        limits.max_sam_line_bytes(),
        AlignmentRecordResource::SamLineBytes,
    )?;
    let capacity = usize::try_from(length).map_err(|_| AlignmentRecordError::AllocationFailed {
        allocation: AlignmentRecordAllocation::SamText,
        requested: length,
    })?;
    let mut output = Vec::new();
    output
        .try_reserve_exact(capacity)
        .map_err(|_| AlignmentRecordError::AllocationFailed {
            allocation: AlignmentRecordAllocation::SamText,
            requested: length,
        })?;
    render_record(record, &mut output);
    debug_assert_eq!(output.len(), capacity);
    Ok(output)
}

/// Encodes the compact borrowed alignment-record contract as SAM text.
///
/// The header supplies the reference names for the ordinal-only compact view;
/// neither the conversion nor the rendering allocates an owned alignment
/// record. BAM and SAM can therefore consume the same batch-backed records.
///
/// # Errors
///
/// Returns an error when an ordinal is absent from `header`, when the exact
/// rendered line exceeds `limits`, or when its output cannot be reserved.
#[allow(clippy::too_many_lines)]
pub fn sam_borrowed_record_bytes(
    record: &BorrowedAlignmentRecord<'_>,
    header: &SamHeader,
    limits: AlignmentRecordLimits,
) -> Result<Vec<u8>, AlignmentRecordError> {
    let reference = header_reference(header, record.reference_ordinal())?;
    let mate = header_reference(header, record.mate_reference_ordinal())?;
    let cigar_bytes = borrowed_cigar_text_length(record, limits)?;
    let rname_bytes = reference.map_or(1, |entry| storage_len(entry.name().len()));
    let rnext_bytes = mate.map_or(1, |entry| {
        if record.reference_ordinal() == record.mate_reference_ordinal() {
            1
        } else {
            storage_len(entry.name().len())
        }
    });
    let quality_bytes = record
        .quality()
        .map_or(1, |quality| storage_len(quality.len()));
    let mut length = storage_len(record.query_name().len()) + 1;
    length = add_length(length, decimal_digits(u64::from(record.flag())) + 1)?;
    length = add_length(length, rname_bytes + 1)?;
    let position = reference.map_or(0, |_| record.position());
    let mate_position = mate.map_or(0, |_| record.mate_position());
    length = add_length(length, decimal_digits(u64::from(position)) + 1)?;
    length = add_length(
        length,
        decimal_digits(u64::from(record.mapping_quality())) + 1,
    )?;
    length = add_length(length, cigar_bytes + 1)?;
    length = add_length(length, rnext_bytes + 1)?;
    length = add_length(length, decimal_digits(u64::from(mate_position)) + 1)?;
    length = add_length(length, signed_digits(record.template_length()) + 1)?;
    length = add_length(length, storage_len(record.sequence().len()) + 1)?;
    length = add_length(length, quality_bytes)?;
    if reference.is_some() {
        length = add_length(length, 6 + decimal_digits(u64::from(record.literal_nm())))?;
        length = add_length(length, 8)?;
        if matches!(
            record.auxiliary_mode(),
            crate::AlignmentAuxiliaryMode::Bismark
        ) {
            let md = record
                .md()
                .ok_or(AlignmentRecordError::InvalidCompactRecord {
                    reason: "Bismark output requires MD",
                })?;
            let xm = record
                .bismark_xm()
                .ok_or(AlignmentRecordError::InvalidCompactRecord {
                    reason: "Bismark output requires XM",
                })?;
            length = add_length(length, 6 + storage_len(md.len()))?;
            length = add_length(length, 14 + storage_len(xm.len()))?;
        }
    }
    length = add_length(length, 1)?;
    check_limit(
        length,
        limits.max_sam_line_bytes(),
        AlignmentRecordResource::SamLineBytes,
    )?;
    let capacity = usize::try_from(length).map_err(|_| AlignmentRecordError::AllocationFailed {
        allocation: AlignmentRecordAllocation::SamText,
        requested: length,
    })?;
    let mut output = Vec::new();
    output
        .try_reserve_exact(capacity)
        .map_err(|_| AlignmentRecordError::AllocationFailed {
            allocation: AlignmentRecordAllocation::SamText,
            requested: length,
        })?;
    output.extend_from_slice(record.query_name());
    output.push(b'\t');
    append_u64(&mut output, u64::from(record.flag()));
    output.push(b'\t');
    match reference {
        Some(entry) => output.extend_from_slice(entry.name()),
        None => output.push(b'*'),
    }
    output.push(b'\t');
    append_u64(&mut output, u64::from(position));
    output.push(b'\t');
    append_u64(&mut output, u64::from(record.mapping_quality()));
    output.push(b'\t');
    if reference.is_some() {
        render_borrowed_cigar(record.cigar(), &mut output);
    } else {
        output.push(b'*');
    }
    output.push(b'\t');
    match mate {
        Some(_) if record.reference_ordinal() == record.mate_reference_ordinal() => {
            output.push(b'=');
        }
        Some(entry) => output.extend_from_slice(entry.name()),
        None => output.push(b'*'),
    }
    output.push(b'\t');
    append_u64(&mut output, u64::from(mate_position));
    output.push(b'\t');
    append_i32(&mut output, record.template_length());
    output.push(b'\t');
    output.extend_from_slice(record.sequence());
    output.push(b'\t');
    match record.quality() {
        Some(quality) => output.extend_from_slice(quality),
        None => output.push(b'*'),
    }
    if reference.is_some() {
        output.extend_from_slice(b"\tNM:i:");
        append_u64(&mut output, u64::from(record.literal_nm()));
        output.extend_from_slice(b"\tXG:Z:");
        output.extend_from_slice(record.bismark_xg());
        if matches!(
            record.auxiliary_mode(),
            crate::AlignmentAuxiliaryMode::Bismark
        ) {
            let md = record
                .md()
                .ok_or(AlignmentRecordError::InvalidCompactRecord {
                    reason: "Bismark output requires MD",
                })?;
            let xm = record
                .bismark_xm()
                .ok_or(AlignmentRecordError::InvalidCompactRecord {
                    reason: "Bismark output requires XM",
                })?;
            output.extend_from_slice(b"\tMD:Z:");
            output.extend_from_slice(md);
            output.extend_from_slice(b"\tXM:Z:");
            output.extend_from_slice(xm);
            output.extend_from_slice(b"\tXR:Z:");
            output.extend_from_slice(record.bismark_xr());
        }
    }
    output.push(b'\n');
    debug_assert_eq!(output.len(), capacity);
    Ok(output)
}

fn header_reference(
    header: &SamHeader,
    ordinal: Option<u64>,
) -> Result<Option<&SamHeaderReference>, AlignmentRecordError> {
    let Some(ordinal) = ordinal else {
        return Ok(None);
    };
    let index = usize::try_from(ordinal).map_err(|_| AlignmentRecordError::FieldOutOfRange {
        field: AlignmentRecordField::ReferenceOrdinal,
        value: ordinal,
    })?;
    header
        .references()
        .get(index)
        .map(Some)
        .ok_or(AlignmentRecordError::FieldOutOfRange {
            field: AlignmentRecordField::ReferenceOrdinal,
            value: ordinal,
        })
}

fn borrowed_cigar_text_length(
    record: &BorrowedAlignmentRecord<'_>,
    limits: AlignmentRecordLimits,
) -> Result<u64, AlignmentRecordError> {
    if record.reference_ordinal().is_none() {
        return Ok(1);
    }
    let mut length = 0;
    for run in record.cigar() {
        length = checked_add_resource(
            length,
            decimal_digits(run.length()) + 1,
            AlignmentRecordResource::CigarTextBytes,
        )?;
    }
    check_limit(
        storage_len(record.cigar().len()),
        limits.max_cigar_runs(),
        AlignmentRecordResource::CigarRuns,
    )?;
    check_limit(
        length,
        limits.max_cigar_text_bytes(),
        AlignmentRecordResource::CigarTextBytes,
    )?;
    Ok(length)
}

fn render_borrowed_cigar(cigar: &[AlignmentCigarRun], output: &mut Vec<u8>) {
    for run in cigar {
        append_u64(output, run.length());
        output.push(match run.operation() {
            AlignmentCigarOp::Match => b'M',
            AlignmentCigarOp::Insertion => b'I',
            AlignmentCigarOp::Deletion => b'D',
            AlignmentCigarOp::SoftClip => b'S',
        });
    }
}

/// Writes one complete canonical SAM alignment line.
///
/// # Errors
///
/// Returns an encode error before writing, or the caller writer's first error.
pub fn write_sam_record<W: Write>(
    writer: &mut W,
    record: &AlignmentRecord,
    limits: AlignmentRecordLimits,
) -> Result<(), SamWriteError> {
    let bytes =
        sam_record_bytes(record, limits).map_err(|source| SamWriteError::Encode { source })?;
    writer
        .write_all(&bytes)
        .map_err(|source| SamWriteError::Io {
            phase: SamWritePhase::Record,
            source,
        })
}

/// Encodes the complete canonical SAM header.
///
/// # Errors
///
/// Returns [`AlignmentRecordError`] for size/limit/allocation failure.
pub fn sam_header_bytes(
    header: &SamHeader,
    limits: AlignmentRecordLimits,
) -> Result<Vec<u8>, AlignmentRecordError> {
    let length = header_text_length(header)?;
    check_limit(
        length,
        limits.max_header_bytes(),
        AlignmentRecordResource::HeaderBytes,
    )?;
    let capacity = usize::try_from(length).map_err(|_| AlignmentRecordError::AllocationFailed {
        allocation: AlignmentRecordAllocation::SamText,
        requested: length,
    })?;
    let mut output = Vec::new();
    output
        .try_reserve_exact(capacity)
        .map_err(|_| AlignmentRecordError::AllocationFailed {
            allocation: AlignmentRecordAllocation::SamText,
            requested: length,
        })?;
    output.extend_from_slice(header.sort_order.header_prefix());
    for reference in &header.references {
        output.extend_from_slice(b"@SQ\tSN:");
        output.extend_from_slice(&reference.name);
        output.extend_from_slice(b"\tLN:");
        append_u64(&mut output, u64::from(reference.length));
        output.push(b'\n');
    }
    output.extend_from_slice(PROGRAM_PREFIX);
    output.extend_from_slice(PROGRAM_VERSION);
    if let Some(provenance) = header.bsbit_provenance {
        output.extend_from_slice(PROGRAM_DESCRIPTION_PREFIX);
        append_hex(&mut output, &provenance.reference_semantic_digest);
        output.extend_from_slice(PROGRAM_MODE_PREFIX);
        output.extend_from_slice(provenance.alignment_mode.header_value());
    }
    output.push(b'\n');
    debug_assert_eq!(output.len(), capacity);
    Ok(output)
}

/// Writes one complete canonical SAM header.
///
/// # Errors
///
/// Returns an encode error before writing, or the caller writer's first error.
pub fn write_sam_header<W: Write>(
    writer: &mut W,
    header: &SamHeader,
    limits: AlignmentRecordLimits,
) -> Result<(), SamWriteError> {
    let bytes =
        sam_header_bytes(header, limits).map_err(|source| SamWriteError::Encode { source })?;
    writer
        .write_all(&bytes)
        .map_err(|source| SamWriteError::Io {
            phase: SamWritePhase::Header,
            source,
        })
}

fn record_text_length(record: &AlignmentRecord) -> Result<u64, AlignmentRecordError> {
    let mapping = record.mapping();
    let mate = record.mate();
    let rname_len = mapping.map_or(1, |mapping| storage_len(mapping.reference().name().len()));
    let position = mapping.map_or(0, |mapping| u64::from(mapping.reference().position()));
    let cigar_len = match mapping {
        Some(mapping) => cigar_text_length(mapping.cigar())?,
        None => 1,
    };
    let rnext_len = mate.map_or(1, |mate| {
        if mapping
            .is_some_and(|mapping| mapping.reference().ordinal() == mate.reference().ordinal())
        {
            1
        } else {
            storage_len(mate.reference().name().len())
        }
    });
    let pnext = mate.map_or(0, |mate| u64::from(mate.reference().position()));
    let quality_len = record
        .quality()
        .map_or(1, |quality| storage_len(quality.len()));

    let mut length = storage_len(record.query_name().len()) + 1;
    length = add_length(length, decimal_digits(u64::from(sam_flag(record))) + 1)?;
    length = add_length(length, rname_len + 1)?;
    length = add_length(length, decimal_digits(position) + 1)?;
    length = add_length(
        length,
        decimal_digits(u64::from(record.mapping_quality().sam_value())) + 1,
    )?;
    length = add_length(length, cigar_len + 1)?;
    length = add_length(length, rnext_len + 1)?;
    length = add_length(length, decimal_digits(pnext) + 1)?;
    length = add_length(length, signed_digits(record.template_length()) + 1)?;
    length = add_length(length, storage_len(record.sequence().len()) + 1)?;
    length = add_length(length, quality_len)?;
    if let Some(mapping) = mapping {
        length = add_length(length, 6 + decimal_digits(u64::from(mapping.literal_nm())))?;
        length = add_length(length, 8)?;
        if matches!(mapping.auxiliary_mode(), AlignmentAuxiliaryMode::Bismark) {
            length = add_length(
                length,
                6 + storage_len(mapping.md().expect("mode requires MD").len()),
            )?;
            length = add_length(
                length,
                14 + storage_len(mapping.bismark_xm().expect("mode requires XM").len()),
            )?;
        }
    }
    add_length(length, 1)
}

fn render_record(record: &AlignmentRecord, output: &mut Vec<u8>) {
    let mapping = record.mapping();
    let mate = record.mate();
    output.extend_from_slice(record.query_name());
    output.push(b'\t');
    append_u64(output, u64::from(sam_flag(record)));
    output.push(b'\t');
    if let Some(mapping) = mapping {
        output.extend_from_slice(mapping.reference().name());
    } else {
        output.push(b'*');
    }
    output.push(b'\t');
    append_u64(
        output,
        mapping.map_or(0, |mapping| u64::from(mapping.reference().position())),
    );
    output.push(b'\t');
    append_u64(output, u64::from(record.mapping_quality().sam_value()));
    output.push(b'\t');
    if let Some(mapping) = mapping {
        render_cigar(mapping.cigar(), output);
    } else {
        output.push(b'*');
    }
    output.push(b'\t');
    if let Some(mate) = mate {
        if mapping
            .is_some_and(|mapping| mapping.reference().ordinal() == mate.reference().ordinal())
        {
            output.push(b'=');
        } else {
            output.extend_from_slice(mate.reference().name());
        }
    } else {
        output.push(b'*');
    }
    output.push(b'\t');
    append_u64(
        output,
        mate.map_or(0, |mate| u64::from(mate.reference().position())),
    );
    output.push(b'\t');
    append_i32(output, record.template_length());
    output.push(b'\t');
    output.extend_from_slice(record.sequence());
    output.push(b'\t');
    if let Some(quality) = record.quality() {
        output.extend_from_slice(quality);
    } else {
        output.push(b'*');
    }
    if let Some(mapping) = mapping {
        output.extend_from_slice(b"\tNM:i:");
        append_u64(output, u64::from(mapping.literal_nm()));
        output.extend_from_slice(b"\tXG:Z:");
        output.extend_from_slice(mapping.bismark_xg());
        if matches!(mapping.auxiliary_mode(), AlignmentAuxiliaryMode::Bismark) {
            output.extend_from_slice(b"\tMD:Z:");
            output.extend_from_slice(mapping.md().expect("mode requires MD"));
            output.extend_from_slice(b"\tXM:Z:");
            output.extend_from_slice(mapping.bismark_xm().expect("mode requires XM"));
            output.extend_from_slice(b"\tXR:Z:");
            output.extend_from_slice(mapping.bismark_xr());
        }
    }
    output.push(b'\n');
}

fn render_cigar(cigar: &CoreCigar, output: &mut Vec<u8>) {
    for run in cigar.runs() {
        append_u64(output, run.length());
        output.push(match run.operation() {
            bsbit_core::cigar::CoreCigarOp::M => b'M',
            bsbit_core::cigar::CoreCigarOp::I => b'I',
            bsbit_core::cigar::CoreCigarOp::D => b'D',
        });
    }
}

fn header_text_length(header: &SamHeader) -> Result<u64, AlignmentRecordError> {
    let mut length = storage_len(header.sort_order.header_prefix().len());
    for reference in &header.references {
        length = add_length(length, 7 + storage_len(reference.name.len()))?;
        length = add_length(length, 4 + decimal_digits(u64::from(reference.length)) + 1)?;
    }
    length = add_length(length, storage_len(PROGRAM_PREFIX.len()))?;
    length = add_length(length, storage_len(PROGRAM_VERSION.len()))?;
    if let Some(provenance) = header.bsbit_provenance {
        length = add_length(length, storage_len(PROGRAM_DESCRIPTION_PREFIX.len()))?;
        length = add_length(length, 64)?;
        length = add_length(length, storage_len(PROGRAM_MODE_PREFIX.len()))?;
        length = add_length(
            length,
            storage_len(provenance.alignment_mode.header_value().len()),
        )?;
    }
    length = add_length(length, 1)?;
    Ok(length)
}

fn append_hex(output: &mut Vec<u8>, bytes: &[u8]) {
    const HEX: &[u8; 16] = b"0123456789abcdef";
    for &byte in bytes {
        output.push(HEX[usize::from(byte >> 4)]);
        output.push(HEX[usize::from(byte & 0x0f)]);
    }
}

pub(crate) fn parse_bsbit_program_provenance(
    text: &[u8],
) -> Result<Option<BsbitProgramProvenance>, BsbitProgramProvenanceError> {
    let mut provenance = None;
    for line in text.split(|byte| *byte == b'\n') {
        let mut fields = line.split(|byte| *byte == b'\t');
        if fields.next() != Some(b"@PG".as_slice()) {
            continue;
        }
        let fields = fields.collect::<Vec<_>>();
        let is_bsbit = fields.iter().any(|field| *field == b"ID:bsbit")
            && fields.iter().any(|field| *field == b"PN:bsbit");
        if !is_bsbit {
            continue;
        }
        if provenance.is_some() {
            return Err(BsbitProgramProvenanceError::DuplicateProgramRecord);
        }
        let description = fields
            .iter()
            .find_map(|field| field.strip_prefix(b"DS:"))
            .ok_or(BsbitProgramProvenanceError::MissingDescription)?;
        provenance = Some(parse_program_description(description)?);
    }
    Ok(provenance)
}

fn parse_program_description(
    description: &[u8],
) -> Result<BsbitProgramProvenance, BsbitProgramProvenanceError> {
    let digest_text = description
        .strip_prefix(b"reference-semantic-sha256=")
        .ok_or(BsbitProgramProvenanceError::MalformedDescription)?;
    let delimiter = b";alignment-mode=";
    let delimiter_offset = digest_text
        .windows(delimiter.len())
        .position(|window| window == delimiter)
        .ok_or(BsbitProgramProvenanceError::MalformedDescription)?;
    let (digest_text, mode_text) = digest_text.split_at(delimiter_offset);
    let mode_text = &mode_text[delimiter.len()..];
    if digest_text.len() != 64 {
        return Err(BsbitProgramProvenanceError::MalformedDescription);
    }
    let mut digest = [0_u8; 32];
    for (target, pair) in digest.iter_mut().zip(digest_text.chunks_exact(2)) {
        *target = (hex_digit(pair[0])? << 4) | hex_digit(pair[1])?;
    }
    let alignment_mode = BsbitAlignmentMode::from_header_value(mode_text)
        .ok_or(BsbitProgramProvenanceError::MalformedDescription)?;
    Ok(BsbitProgramProvenance::new(digest, alignment_mode))
}

fn hex_digit(byte: u8) -> Result<u8, BsbitProgramProvenanceError> {
    match byte {
        b'0'..=b'9' => Ok(byte - b'0'),
        b'a'..=b'f' => Ok(byte - b'a' + 10),
        _ => Err(BsbitProgramProvenanceError::MalformedDescription),
    }
}

fn add_length(current: u64, increment: u64) -> Result<u64, AlignmentRecordError> {
    checked_add_resource(current, increment, AlignmentRecordResource::SamLineBytes)
}

fn signed_digits(value: i32) -> u64 {
    let magnitude = u64::from(value.unsigned_abs());
    decimal_digits(magnitude) + u64::from(value < 0)
}

fn append_i32(output: &mut Vec<u8>, value: i32) {
    if value < 0 {
        output.push(b'-');
        append_u64(output, u64::from(value.unsigned_abs()));
    } else {
        append_u64(output, u64::from(value.unsigned_abs()));
    }
}

/// File lifecycle phase associated with a SAM publication failure.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SamFilePhase {
    /// Target and staging paths are not distinct siblings.
    ValidatePaths,
    /// The target already exists or cannot be inspected.
    ValidateTarget,
    /// The complete canonical header could not be encoded.
    EncodeHeader,
    /// Exclusive staging-file creation failed.
    CreateStaging,
    /// Writing the already encoded header failed.
    WriteHeader,
    /// The next complete canonical record could not be encoded.
    EncodeRecord,
    /// The next already encoded record could not be written.
    WriteRecord,
    /// The record ordinal could not be incremented.
    CountRecord,
    /// The writer was already terminal after an earlier failure.
    Closed,
    /// Buffered output could not be flushed or recovered.
    Flush,
    /// Generic completion or publication failed.
    Publish,
    /// Explicit staging cleanup failed during abort.
    Abort,
}

/// A SAM file failure that never reports a newly published target.
#[derive(Debug)]
pub struct SamFileError {
    phase: SamFilePhase,
    record_ordinal: Option<u64>,
    io_error: Option<io::Error>,
    encode_error: Option<AlignmentRecordError>,
    staging_created: bool,
    cleanup_error: Option<io::ErrorKind>,
}

impl SamFileError {
    /// Returns the lifecycle phase that failed.
    #[must_use]
    pub const fn phase(&self) -> SamFilePhase {
        self.phase
    }

    /// Returns the one-based record ordinal for a record-local failure.
    #[must_use]
    pub const fn record_ordinal(&self) -> Option<u64> {
        self.record_ordinal
    }

    /// Returns the direct filesystem error kind, when applicable.
    #[must_use]
    pub fn kind(&self) -> Option<io::ErrorKind> {
        self.io_error.as_ref().map(io::Error::kind)
    }

    /// Returns the canonical record/header encoding failure, when applicable.
    #[must_use]
    pub const fn encode_error(&self) -> Option<&AlignmentRecordError> {
        self.encode_error.as_ref()
    }

    /// Reports whether this invocation created its staging file.
    #[must_use]
    pub const fn staging_created(&self) -> bool {
        self.staging_created
    }

    /// Returns a best-effort cleanup failure after the primary failure.
    #[must_use]
    pub const fn cleanup_error(&self) -> Option<io::ErrorKind> {
        self.cleanup_error
    }
}

impl fmt::Display for SamFileError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "SAM file failed in {:?}", self.phase)?;
        if let Some(ordinal) = self.record_ordinal {
            write!(formatter, " at record {ordinal}")?;
        }
        if let Some(source) = &self.encode_error {
            write!(formatter, ": {source}")?;
        } else if let Some(source) = &self.io_error {
            write!(formatter, ": {source}")?;
        }
        Ok(())
    }
}

impl std::error::Error for SamFileError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        self.encode_error
            .as_ref()
            .map(|source| source as &(dyn std::error::Error + 'static))
            .or_else(|| {
                self.io_error
                    .as_ref()
                    .map(|source| source as &(dyn std::error::Error + 'static))
            })
    }
}

/// Successful create-only SAM publication details.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SamFilePublication {
    published: PublishedFile,
    records_written: u64,
}

impl SamFilePublication {
    /// Returns the final caller-supplied target path.
    #[must_use]
    pub fn target_path(&self) -> &Path {
        self.published.target_path()
    }

    /// Returns the private staging path.
    #[must_use]
    pub fn staging_path(&self) -> &Path {
        self.published.staging_path()
    }

    /// Returns the number of complete records written after the header.
    #[must_use]
    pub const fn records_written(&self) -> u64 {
        self.records_written
    }

    /// Reports whether the valid staging link could not be removed.
    #[must_use]
    pub const fn staging_retained(&self) -> bool {
        self.published.cleanup_warning().is_some()
    }

    /// Returns a post-publication staging cleanup warning.
    #[must_use]
    pub const fn cleanup_error(&self) -> Option<io::ErrorKind> {
        self.published.cleanup_warning()
    }
}

/// Terminal-on-error canonical SAM encoder over generic publication state.
#[derive(Debug)]
pub struct SamFileWriter {
    target: PathBuf,
    header: SamHeader,
    staged: Option<StagedFile>,
    writer: Option<BufWriter<File>>,
    records_written: u64,
}

impl SamFileWriter {
    /// Creates an exclusive staging file and writes one complete SAM header.
    ///
    /// # Errors
    ///
    /// Returns path, header-encoding, staging-creation, or header-write errors.
    pub fn create_new(
        target: impl AsRef<Path>,
        staging: impl AsRef<Path>,
        header: &SamHeader,
        limits: AlignmentRecordLimits,
    ) -> Result<Self, SamFileError> {
        let target = bsbit_io::absolute_path(target.as_ref()).map_err(|source| {
            sam_file_error(SamFilePhase::ValidatePaths, None, Some(source), None, false)
        })?;
        let staging_path = bsbit_io::absolute_path(staging.as_ref()).map_err(|source| {
            sam_file_error(SamFilePhase::ValidatePaths, None, Some(source), None, false)
        })?;
        if target == staging_path || target.parent() != staging_path.parent() {
            return Err(sam_file_error(
                SamFilePhase::ValidatePaths,
                None,
                Some(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "target and staging must be distinct siblings",
                )),
                None,
                false,
            ));
        }
        bsbit_io::validate_absent(&target).map_err(|source| {
            sam_file_error(
                SamFilePhase::ValidateTarget,
                None,
                Some(source),
                None,
                false,
            )
        })?;
        let header_bytes = sam_header_bytes(header, limits).map_err(|source| {
            sam_file_error(SamFilePhase::EncodeHeader, None, None, Some(source), false)
        })?;
        let mut staged = StagedFile::create_new(&staging_path).map_err(|source| {
            publication_as_sam_error(SamFilePhase::CreateStaging, None, source, false)
        })?;
        let file = staged.take_file().map_err(|source| {
            publication_as_sam_error(SamFilePhase::CreateStaging, None, source, true)
        })?;
        let mut writer = BufWriter::new(file);
        if let Err(source) = writer.write_all(&header_bytes) {
            drop(writer);
            drop(staged);
            return Err(sam_file_error(
                SamFilePhase::WriteHeader,
                None,
                Some(source),
                None,
                true,
            ));
        }
        Ok(Self {
            target,
            header: header.clone(),
            staged: Some(staged),
            writer: Some(writer),
            records_written: 0,
        })
    }

    /// Writes one complete canonical record or poisons this writer.
    ///
    /// # Errors
    ///
    /// Returns record encoding, writing, counter, or terminal-state errors.
    pub fn write_record(
        &mut self,
        record: &AlignmentRecord,
        limits: AlignmentRecordLimits,
    ) -> Result<(), SamFileError> {
        let ordinal = self
            .records_written
            .checked_add(1)
            .ok_or_else(|| sam_file_error(SamFilePhase::CountRecord, None, None, None, true))?;
        let bytes = sam_record_bytes(record, limits).map_err(|source| {
            self.poison(
                SamFilePhase::EncodeRecord,
                Some(ordinal),
                None,
                Some(source),
            )
        })?;
        let write_result = self
            .writer
            .as_mut()
            .ok_or_else(|| sam_file_error(SamFilePhase::Closed, Some(ordinal), None, None, true))?
            .write_all(&bytes);
        if let Err(source) = write_result {
            return Err(self.poison(SamFilePhase::WriteRecord, Some(ordinal), Some(source), None));
        }
        self.records_written = ordinal;
        Ok(())
    }

    /// Writes one compact batch-backed record using this writer's dictionary.
    ///
    /// # Errors
    ///
    /// Returns record encoding, writing, counter, or terminal-state errors.
    pub fn write_borrowed_record(
        &mut self,
        record: &BorrowedAlignmentRecord<'_>,
        limits: AlignmentRecordLimits,
    ) -> Result<(), SamFileError> {
        let ordinal = self
            .records_written
            .checked_add(1)
            .ok_or_else(|| sam_file_error(SamFilePhase::CountRecord, None, None, None, true))?;
        let bytes = sam_borrowed_record_bytes(record, &self.header, limits).map_err(|source| {
            self.poison(
                SamFilePhase::EncodeRecord,
                Some(ordinal),
                None,
                Some(source),
            )
        })?;
        let write_result = self
            .writer
            .as_mut()
            .ok_or_else(|| sam_file_error(SamFilePhase::Closed, Some(ordinal), None, None, true))?
            .write_all(&bytes);
        if let Err(source) = write_result {
            return Err(self.poison(SamFilePhase::WriteRecord, Some(ordinal), Some(source), None));
        }
        self.records_written = ordinal;
        Ok(())
    }

    /// Returns the number of complete records accepted by this writer.
    #[must_use]
    pub const fn records_written(&self) -> u64 {
        self.records_written
    }

    /// Flushes, completes, synchronizes, and publishes this SAM create-only.
    ///
    /// # Errors
    ///
    /// Returns terminal, flush, identity, synchronization, or target-race errors.
    pub fn finish(mut self) -> Result<SamFilePublication, SamFileError> {
        let mut writer = self
            .writer
            .take()
            .ok_or_else(|| sam_file_error(SamFilePhase::Closed, None, None, None, true))?;
        writer.flush().map_err(|source| {
            sam_file_error(SamFilePhase::Flush, None, Some(source), None, true)
        })?;
        let file = writer.into_inner().map_err(|error| {
            sam_file_error(
                SamFilePhase::Flush,
                None,
                Some(error.into_error()),
                None,
                true,
            )
        })?;
        let staged = self
            .staged
            .take()
            .ok_or_else(|| sam_file_error(SamFilePhase::Closed, None, None, None, true))?;
        let completed = staged.complete(file).map_err(|source| {
            publication_as_sam_error(SamFilePhase::Publish, None, source, true)
        })?;
        let published = completed
            .publish_create_new_at(&self.target)
            .map_err(|source| {
                publication_as_sam_error(SamFilePhase::Publish, None, source, true)
            })?;
        Ok(SamFilePublication {
            published,
            records_written: self.records_written,
        })
    }

    /// Explicitly closes and removes the unfinished private staging file.
    ///
    /// # Errors
    ///
    /// Returns identity or cleanup errors without touching the target.
    pub fn abort(mut self) -> Result<(), SamFileError> {
        self.writer.take();
        match self.staged.take() {
            Some(staged) => staged.abort().map_err(|source| {
                publication_as_sam_error(SamFilePhase::Abort, None, source, true)
            }),
            None => Ok(()),
        }
    }

    fn poison(
        &mut self,
        phase: SamFilePhase,
        ordinal: Option<u64>,
        io_error: Option<io::Error>,
        encode_error: Option<AlignmentRecordError>,
    ) -> SamFileError {
        self.writer.take();
        self.staged.take();
        sam_file_error(phase, ordinal, io_error, encode_error, true)
    }
}

impl Drop for SamFileWriter {
    fn drop(&mut self) {
        self.writer.take();
        self.staged.take();
    }
}

fn publication_as_sam_error(
    phase: SamFilePhase,
    ordinal: Option<u64>,
    error: PublicationError,
    staging_created: bool,
) -> SamFileError {
    sam_file_error(
        phase,
        ordinal,
        Some(error.into_io_error()),
        None,
        staging_created,
    )
}

fn sam_file_error(
    phase: SamFilePhase,
    record_ordinal: Option<u64>,
    io_error: Option<io::Error>,
    encode_error: Option<AlignmentRecordError>,
    staging_created: bool,
) -> SamFileError {
    SamFileError {
        phase,
        record_ordinal,
        io_error,
        encode_error,
        staging_created,
        cleanup_error: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn signed_i32_boundaries_have_exact_length_and_text() {
        for (value, expected) in [
            (i32::MIN, b"-2147483648".as_slice()),
            (-1, b"-1".as_slice()),
            (0, b"0".as_slice()),
            (1, b"1".as_slice()),
            (i32::MAX, b"2147483647".as_slice()),
        ] {
            let mut output = Vec::new();
            append_i32(&mut output, value);
            assert_eq!(output, expected);
            assert_eq!(storage_len(output.len()), signed_digits(value));
        }
    }

    #[test]
    fn structured_bsbit_provenance_is_exact_unique_and_fail_closed() {
        let limits = AlignmentRecordLimits::default();
        let expected = BsbitProgramProvenance::new(
            [0xab; 32],
            BsbitAlignmentMode::CallerCompatibleNondirectionalPaired,
        );
        let header = SamHeader::new(
            vec![SamHeaderReference::new(0, b"chr1", 7).expect("reference")],
            limits,
        )
        .expect("header")
        .with_bsbit_provenance(expected, limits)
        .expect("provenance header");
        let bytes = sam_header_bytes(&header, limits).expect("SAM header bytes");
        assert!(
            bytes.windows(64).any(|window| window
                == b"abababababababababababababababababababababababababababababababab")
        );
        assert_eq!(
            parse_bsbit_program_provenance(&bytes).expect("generated provenance parses"),
            Some(expected)
        );

        assert_eq!(
            parse_bsbit_program_provenance(b"@PG\tID:bsbit\tPN:bsbit\n"),
            Err(BsbitProgramProvenanceError::MissingDescription)
        );
        let mut duplicate = bytes.clone();
        duplicate.extend_from_slice(
            b"@PG\tID:bsbit\tPN:bsbit\tDS:reference-semantic-sha256=abababababababababababababababababababababababababababababababab;alignment-mode=caller-compatible-directional-single\n",
        );
        assert_eq!(
            parse_bsbit_program_provenance(&duplicate),
            Err(BsbitProgramProvenanceError::DuplicateProgramRecord)
        );
    }

    #[test]
    fn every_alignment_mode_round_trips_and_is_caller_compatible() {
        for mode in [
            BsbitAlignmentMode::CallerCompatibleDirectionalSingle,
            BsbitAlignmentMode::CallerCompatibleDirectionalPaired,
            BsbitAlignmentMode::CallerCompatibleNondirectionalPaired,
        ] {
            assert_eq!(
                BsbitAlignmentMode::from_header_value(mode.header_value()),
                Some(mode)
            );
            assert!(mode.is_caller_compatible());
        }
    }
}
