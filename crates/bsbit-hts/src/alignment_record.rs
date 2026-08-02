//! Validated alignment data shared by SAM and BAM encoders.
//!
//! Constructors accept only format-level values. Mapping-result and reference-
//! index adaptation belongs to the application layer, so this module has no
//! dependency on the aligner or index crates.

use core::fmt;
use core::mem::size_of;

use bsbit_core::alphabet::Base;
use bsbit_core::bisulfite::{AlignmentOrientation, BisulfiteStrand, CytosineStrand};
use bsbit_core::cigar::{CigarError, CoreCigar, validate_cigar};
use bsbit_core::coordinate::ReferenceInterval;
use bsbit_core::sequence::NormalizedSequence;

/// SAM-representable hard maximum for a query name.
pub const SAM_MAX_QUERY_NAME_BYTES: u64 = 254;
/// SAM-representable maximum reference length and position.
pub const SAM_MAX_REFERENCE_LENGTH: u64 = i32::MAX as u64;

/// CIGAR operation shared by compact SAM and BAM output.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum AlignmentCigarOp {
    /// Consume one reference and one query base.
    Match,
    /// Consume query only.
    Insertion,
    /// Consume reference only.
    Deletion,
    /// Consume unaligned query bases at a CIGAR end.
    SoftClip,
}

/// One positive CIGAR run shared by compact SAM and BAM output.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct AlignmentCigarRun {
    operation: AlignmentCigarOp,
    length: u64,
}

impl AlignmentCigarRun {
    /// Creates one positive run.
    ///
    /// # Errors
    ///
    /// Returns a field error for a zero length.
    pub fn new(operation: AlignmentCigarOp, length: u64) -> Result<Self, AlignmentRecordError> {
        if length == 0 {
            return Err(AlignmentRecordError::FieldOutOfRange {
                field: AlignmentRecordField::CigarRunLength,
                value: 0,
            });
        }
        Ok(Self { operation, length })
    }

    /// Returns the operation.
    #[must_use]
    pub const fn operation(self) -> AlignmentCigarOp {
        self.operation
    }

    /// Returns the positive run length.
    #[must_use]
    pub const fn length(self) -> u64 {
        self.length
    }
}

/// One read and optional quality string borrowed during record composition.
#[derive(Clone, Copy, Debug)]
pub struct AlignmentRead<'a> {
    sequence: &'a NormalizedSequence,
    quality: Option<&'a [u8]>,
}

impl<'a> AlignmentRead<'a> {
    /// Creates a borrowed read view.
    #[must_use]
    pub const fn new(sequence: &'a NormalizedSequence, quality: Option<&'a [u8]>) -> Self {
        Self { sequence, quality }
    }

    /// Returns the normalized sequence.
    #[must_use]
    pub const fn sequence(self) -> &'a NormalizedSequence {
        self.sequence
    }

    /// Returns optional quality bytes.
    #[must_use]
    pub const fn quality(self) -> Option<&'a [u8]> {
        self.quality
    }
}

/// Logical resource controlled while constructing a record or header.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum AlignmentRecordResource {
    /// Query-name bytes.
    QueryNameBytes,
    /// Read bases.
    ReadBases,
    /// Quality bytes.
    QualityBytes,
    /// Coalesced CIGAR runs.
    CigarRuns,
    /// Rendered CIGAR bytes.
    CigarTextBytes,
    /// Rendered MD bytes.
    MdBytes,
    /// Aggregate optional-field bytes.
    OptionalFieldBytes,
    /// Complete SAM alignment-line bytes.
    SamLineBytes,
    /// Header reference entries.
    HeaderReferences,
    /// Aggregate header reference-name bytes.
    HeaderNameBytes,
    /// Complete header bytes.
    HeaderBytes,
}

/// Fallible allocation site in the alignment boundary.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum AlignmentRecordAllocation {
    /// Query-name storage.
    QueryName,
    /// Sequence storage.
    Sequence,
    /// Quality storage.
    Quality,
    /// Reference-name storage.
    ReferenceName,
    /// MD tag storage.
    Md,
    /// Methylation-call storage.
    MethylationCall,
    /// Compact-record CIGAR storage.
    Cigar,
    /// Header reference storage.
    HeaderReferences,
    /// SAM text storage.
    SamText,
}

/// Numeric field that could not be represented by SAM/BAM.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum AlignmentRecordField {
    /// Reference-dictionary ordinal.
    ReferenceOrdinal,
    /// Mapped one-based position.
    Position,
    /// Numeric SAM mapping quality.
    MappingQuality,
    /// Literal NM value.
    Nm,
    /// Signed template length.
    TemplateLength,
    /// A direct CIGAR run length.
    CigarRunLength,
}

/// Exact caps for one alignment record and one header.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[allow(clippy::struct_field_names)]
pub struct AlignmentRecordLimits {
    max_query_name_bytes: u64,
    max_read_bases: u64,
    max_quality_bytes: u64,
    max_cigar_runs: u64,
    max_cigar_text_bytes: u64,
    max_md_bytes: u64,
    max_optional_field_bytes: u64,
    max_sam_line_bytes: u64,
    max_header_references: u64,
    max_header_name_bytes: u64,
    max_header_bytes: u64,
}

impl AlignmentRecordLimits {
    /// Constructs the complete record/header limit set.
    #[must_use]
    #[allow(clippy::too_many_arguments)]
    pub const fn new(
        max_query_name_bytes: u64,
        max_read_bases: u64,
        max_quality_bytes: u64,
        max_cigar_runs: u64,
        max_cigar_text_bytes: u64,
        max_md_bytes: u64,
        max_optional_field_bytes: u64,
        max_sam_line_bytes: u64,
        max_header_references: u64,
        max_header_name_bytes: u64,
        max_header_bytes: u64,
    ) -> Self {
        Self {
            max_query_name_bytes,
            max_read_bases,
            max_quality_bytes,
            max_cigar_runs,
            max_cigar_text_bytes,
            max_md_bytes,
            max_optional_field_bytes,
            max_sam_line_bytes,
            max_header_references,
            max_header_name_bytes,
            max_header_bytes,
        }
    }

    /// Returns the query-name cap.
    #[must_use]
    pub const fn max_query_name_bytes(self) -> u64 {
        self.max_query_name_bytes
    }
    /// Returns the read-base cap.
    #[must_use]
    pub const fn max_read_bases(self) -> u64 {
        self.max_read_bases
    }
    /// Returns the quality-byte cap.
    #[must_use]
    pub const fn max_quality_bytes(self) -> u64 {
        self.max_quality_bytes
    }
    /// Returns the CIGAR-run cap.
    #[must_use]
    pub const fn max_cigar_runs(self) -> u64 {
        self.max_cigar_runs
    }
    /// Returns the rendered CIGAR cap.
    #[must_use]
    pub const fn max_cigar_text_bytes(self) -> u64 {
        self.max_cigar_text_bytes
    }
    /// Returns the rendered MD cap.
    #[must_use]
    pub const fn max_md_bytes(self) -> u64 {
        self.max_md_bytes
    }
    /// Returns the aggregate optional-field cap.
    #[must_use]
    pub const fn max_optional_field_bytes(self) -> u64 {
        self.max_optional_field_bytes
    }
    /// Returns the complete SAM line cap.
    #[must_use]
    pub const fn max_sam_line_bytes(self) -> u64 {
        self.max_sam_line_bytes
    }
    /// Returns the header reference-count cap.
    #[must_use]
    pub const fn max_header_references(self) -> u64 {
        self.max_header_references
    }
    /// Returns the aggregate header-name cap.
    #[must_use]
    pub const fn max_header_name_bytes(self) -> u64 {
        self.max_header_name_bytes
    }
    /// Returns the complete header cap.
    #[must_use]
    pub const fn max_header_bytes(self) -> u64 {
        self.max_header_bytes
    }
}

impl Default for AlignmentRecordLimits {
    fn default() -> Self {
        Self::new(
            SAM_MAX_QUERY_NAME_BYTES,
            10_000_000,
            10_000_000,
            1_000_000,
            20_000_000,
            20_000_000,
            32_000_000,
            64_000_000,
            1_000_000,
            64_000_000,
            256_000_000,
        )
    }
}

/// Structured format-model construction failure.
#[non_exhaustive]
#[derive(Debug)]
#[allow(missing_docs)]
pub enum AlignmentRecordError {
    /// A configured cap was exceeded.
    LimitExceeded {
        resource: AlignmentRecordResource,
        observed: u64,
        limit: u64,
    },
    /// A checked logical sum overflowed.
    ArithmeticOverflow {
        resource: AlignmentRecordResource,
        current: u64,
        increment: u64,
    },
    /// A bounded allocation could not be reserved.
    AllocationFailed {
        allocation: AlignmentRecordAllocation,
        requested: u64,
    },
    /// QNAME was empty.
    EmptyQueryName,
    /// QNAME contained a byte outside the SAM grammar.
    InvalidQueryNameByte { offset: u64, byte: u8 },
    /// A record had no sequence bases.
    EmptySequence,
    /// Sequence and quality lengths differed.
    QualityLengthMismatch { sequence: u64, quality: u64 },
    /// Quality contained a byte outside printable Phred+33 SAM bytes.
    InvalidQualityByte { offset: u64, byte: u8 },
    /// A reference name violated the SAM grammar.
    InvalidReferenceNameByte {
        ordinal: u64,
        offset: u64,
        byte: Option<u8>,
    },
    /// A reference length cannot be represented by SAM.
    ReferenceLengthOutOfRange { ordinal: u64, length: u64 },
    /// A numeric field cannot be represented.
    FieldOutOfRange {
        field: AlignmentRecordField,
        value: u64,
    },
    /// Core CIGAR validation failed.
    Cigar { source: CigarError },
    /// A mapped record has no CIGAR operations.
    EmptyMappedCigar,
    /// Optional alignment fields do not form a supported output contract.
    InvalidAuxiliaryFields { reason: &'static str },
    /// Compact output fields describe an impossible record state.
    InvalidCompactRecord { reason: &'static str },
}

impl fmt::Display for AlignmentRecordError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::LimitExceeded {
                resource,
                observed,
                limit,
            } => {
                write!(
                    formatter,
                    "alignment resource {resource:?} observed {observed}, exceeding {limit}"
                )
            }
            Self::ArithmeticOverflow {
                resource,
                current,
                increment,
            } => {
                write!(
                    formatter,
                    "alignment resource {resource:?} overflowed: {current} + {increment}"
                )
            }
            Self::AllocationFailed {
                allocation,
                requested,
            } => {
                write!(
                    formatter,
                    "failed to reserve {requested} bytes/elements for {allocation:?}"
                )
            }
            Self::EmptyQueryName => formatter.write_str("SAM query name is empty"),
            Self::InvalidQueryNameByte { offset, byte } => {
                write!(
                    formatter,
                    "invalid SAM query-name byte 0x{byte:02X} at offset {offset}"
                )
            }
            Self::EmptySequence => {
                formatter.write_str("alignment record requires a nonempty sequence")
            }
            Self::QualityLengthMismatch { sequence, quality } => {
                write!(
                    formatter,
                    "sequence length {sequence} differs from quality length {quality}"
                )
            }
            Self::InvalidQualityByte { offset, byte } => {
                write!(
                    formatter,
                    "invalid SAM quality byte 0x{byte:02X} at offset {offset}"
                )
            }
            Self::InvalidReferenceNameByte {
                ordinal,
                offset,
                byte,
            } => match byte {
                Some(byte) => write!(
                    formatter,
                    "reference {ordinal} has invalid name byte 0x{byte:02X} at {offset}"
                ),
                None => write!(formatter, "reference {ordinal} has an empty name"),
            },
            Self::ReferenceLengthOutOfRange { ordinal, length } => {
                write!(
                    formatter,
                    "reference {ordinal} length {length} is outside SAM range"
                )
            }
            Self::FieldOutOfRange { field, value } => {
                write!(
                    formatter,
                    "alignment field {field:?} cannot represent {value}"
                )
            }
            Self::Cigar { source } => {
                write!(formatter, "alignment CIGAR failed validation: {source}")
            }
            Self::EmptyMappedCigar => {
                formatter.write_str("mapped alignment requires a nonempty CIGAR")
            }
            Self::InvalidAuxiliaryFields { reason } => {
                write!(formatter, "invalid alignment auxiliary fields: {reason}")
            }
            Self::InvalidCompactRecord { reason } => {
                write!(formatter, "invalid compact alignment record: {reason}")
            }
        }
    }
}

impl std::error::Error for AlignmentRecordError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Cigar { source } => Some(source),
            _ => None,
        }
    }
}

/// Borrowed bases and qualities used while composing an alignment record.
#[derive(Clone, Copy, Debug)]
pub struct BorrowedAlignmentRead<'a> {
    sequence: &'a [Base],
    quality: &'a [u8],
}

impl<'a> BorrowedAlignmentRead<'a> {
    /// Creates a borrowed read view.
    #[must_use]
    pub const fn new(sequence: &'a [Base], quality: &'a [u8]) -> Self {
        Self { sequence, quality }
    }

    /// Returns the normalized bases.
    #[must_use]
    pub const fn sequence(self) -> &'a [Base] {
        self.sequence
    }

    /// Returns the quality bytes.
    #[must_use]
    pub const fn quality(self) -> &'a [u8] {
        self.quality
    }
}

/// Borrowed validated fields shared by compact SAM and BAM output paths.
pub struct BorrowedAlignmentRecord<'a> {
    query_name: &'a [u8],
    flag: u16,
    reference_ordinal: Option<u64>,
    position: u32,
    mapping_quality: u8,
    cigar: &'a [AlignmentCigarRun],
    mate_reference_ordinal: Option<u64>,
    mate_position: u32,
    template_length: i32,
    sequence: &'a [u8],
    quality: Option<&'a [u8]>,
    literal_nm: u32,
    auxiliary_mode: AlignmentAuxiliaryMode,
    md: Option<&'a [u8]>,
    strand: BisulfiteStrand,
    bismark_xm: Option<&'a [u8]>,
}

#[allow(missing_docs)]
impl<'a> BorrowedAlignmentRecord<'a> {
    /// Creates a compact record after application-level mapping composition.
    ///
    /// # Errors
    ///
    /// Returns basic name, sequence, quality, and CIGAR validation failures.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        query_name: &'a [u8],
        flag: u16,
        reference_ordinal: Option<u64>,
        position: u32,
        mapping_quality: u8,
        cigar: &'a [AlignmentCigarRun],
        mate_reference_ordinal: Option<u64>,
        mate_position: u32,
        template_length: i32,
        sequence: &'a [u8],
        quality: Option<&'a [u8]>,
        literal_nm: u32,
        auxiliary_mode: AlignmentAuxiliaryMode,
        md: Option<&'a [u8]>,
        strand: BisulfiteStrand,
        bismark_xm: Option<&'a [u8]>,
        limits: AlignmentRecordLimits,
    ) -> Result<Self, AlignmentRecordError> {
        validate_query_name(query_name, limits)?;
        validate_sequence_quality(sequence, quality, limits)?;
        check_limit(
            storage_len(cigar.len()),
            limits.max_cigar_runs(),
            AlignmentRecordResource::CigarRuns,
        )?;
        if reference_ordinal.is_some() {
            if position == 0 || cigar.is_empty() {
                return Err(AlignmentRecordError::FieldOutOfRange {
                    field: AlignmentRecordField::Position,
                    value: u64::from(position),
                });
            }
            match auxiliary_mode {
                AlignmentAuxiliaryMode::Minimal if md.is_some() || bismark_xm.is_some() => {
                    return Err(AlignmentRecordError::InvalidCompactRecord {
                        reason: "minimal output cannot carry MD or XM",
                    });
                }
                AlignmentAuxiliaryMode::Bismark if md.is_none() || bismark_xm.is_none() => {
                    return Err(AlignmentRecordError::InvalidCompactRecord {
                        reason: "Bismark output requires both MD and XM",
                    });
                }
                _ => {}
            }
        } else if !cigar.is_empty() || md.is_some() || bismark_xm.is_some() {
            return Err(AlignmentRecordError::InvalidCompactRecord {
                reason: "unmapped output cannot carry CIGAR, MD, or XM",
            });
        }
        let mut cigar_text_bytes = 0_u64;
        for run in cigar {
            cigar_text_bytes = checked_add_resource(
                cigar_text_bytes,
                decimal_digits(run.length()).saturating_add(1),
                AlignmentRecordResource::CigarTextBytes,
            )?;
        }
        check_limit(
            cigar_text_bytes,
            limits.max_cigar_text_bytes(),
            AlignmentRecordResource::CigarTextBytes,
        )?;
        let md_bytes = md.map_or(0, |value| storage_len(value.len()));
        check_limit(
            md_bytes,
            limits.max_md_bytes(),
            AlignmentRecordResource::MdBytes,
        )?;
        let xm_bytes = bismark_xm.map_or(0, |value| storage_len(value.len()));
        let optional_bytes = if reference_ordinal.is_none() {
            0
        } else {
            let optional_bytes = checked_add_resource(
                14,
                decimal_digits(u64::from(literal_nm)),
                AlignmentRecordResource::OptionalFieldBytes,
            )?;
            match auxiliary_mode {
                AlignmentAuxiliaryMode::Minimal => optional_bytes,
                AlignmentAuxiliaryMode::Bismark => {
                    let optional_bytes = checked_add_resource(
                        optional_bytes,
                        6 + md_bytes,
                        AlignmentRecordResource::OptionalFieldBytes,
                    )?;
                    checked_add_resource(
                        optional_bytes,
                        14 + xm_bytes,
                        AlignmentRecordResource::OptionalFieldBytes,
                    )?
                }
            }
        };
        check_limit(
            optional_bytes,
            limits.max_optional_field_bytes(),
            AlignmentRecordResource::OptionalFieldBytes,
        )?;
        Ok(Self {
            query_name,
            flag,
            reference_ordinal,
            position,
            mapping_quality,
            cigar,
            mate_reference_ordinal,
            mate_position,
            template_length,
            sequence,
            quality,
            literal_nm,
            auxiliary_mode,
            md,
            strand,
            bismark_xm,
        })
    }

    #[must_use]
    pub const fn query_name(&self) -> &[u8] {
        self.query_name
    }
    #[must_use]
    pub const fn flag(&self) -> u16 {
        self.flag
    }
    #[must_use]
    pub const fn reference_ordinal(&self) -> Option<u64> {
        self.reference_ordinal
    }
    #[must_use]
    pub const fn position(&self) -> u32 {
        self.position
    }
    #[must_use]
    pub const fn mapping_quality(&self) -> u8 {
        self.mapping_quality
    }
    #[must_use]
    pub const fn cigar(&self) -> &[AlignmentCigarRun] {
        self.cigar
    }
    #[must_use]
    pub const fn mate_reference_ordinal(&self) -> Option<u64> {
        self.mate_reference_ordinal
    }
    #[must_use]
    pub const fn mate_position(&self) -> u32 {
        self.mate_position
    }
    #[must_use]
    pub const fn template_length(&self) -> i32 {
        self.template_length
    }
    #[must_use]
    pub const fn sequence(&self) -> &[u8] {
        self.sequence
    }
    #[must_use]
    pub const fn quality(&self) -> Option<&[u8]> {
        self.quality
    }
    #[must_use]
    pub const fn literal_nm(&self) -> u32 {
        self.literal_nm
    }
    #[must_use]
    pub const fn auxiliary_mode(&self) -> AlignmentAuxiliaryMode {
        self.auxiliary_mode
    }
    #[must_use]
    pub const fn md(&self) -> Option<&[u8]> {
        self.md
    }
    #[must_use]
    pub const fn strand(&self) -> BisulfiteStrand {
        self.strand
    }
    #[must_use]
    pub const fn bismark_xm(&self) -> Option<&[u8]> {
        self.bismark_xm
    }
    #[must_use]
    pub const fn bismark_xr(&self) -> &'static [u8; 2] {
        bismark_read_conversion(self.strand)
    }
    #[must_use]
    pub const fn bismark_xg(&self) -> &'static [u8; 2] {
        bismark_genome_conversion(self.strand)
    }
}

#[derive(Clone, Copy, Debug)]
struct BatchBytes {
    start: usize,
    len: usize,
}

impl BatchBytes {
    fn slice(self, pool: &[u8]) -> &[u8] {
        &pool[self.start..self.start + self.len]
    }
}

#[derive(Clone, Copy, Debug)]
struct BatchRecord {
    query_name: BatchBytes,
    flag: u16,
    reference_ordinal: Option<u64>,
    position: u32,
    mapping_quality: u8,
    cigar_start: usize,
    cigar_len: usize,
    mate_reference_ordinal: Option<u64>,
    mate_position: u32,
    template_length: i32,
    sequence: BatchBytes,
    quality: Option<BatchBytes>,
    literal_nm: u32,
    auxiliary_mode: AlignmentAuxiliaryMode,
    md: Option<BatchBytes>,
    strand: BisulfiteStrand,
    bismark_xm: Option<BatchBytes>,
}

/// Worker-local owned storage for compact alignment records.
#[derive(Clone, Debug, Default)]
pub struct AlignmentRecordBatch {
    records: Vec<BatchRecord>,
    bytes: Vec<u8>,
    cigar: Vec<AlignmentCigarRun>,
    md: Vec<u8>,
    bismark_xm: Vec<u8>,
}

#[allow(missing_docs)]
impl AlignmentRecordBatch {
    #[must_use]
    pub const fn new() -> Self {
        Self {
            records: Vec::new(),
            bytes: Vec::new(),
            cigar: Vec::new(),
            md: Vec::new(),
            bismark_xm: Vec::new(),
        }
    }
    #[must_use]
    pub fn len(&self) -> usize {
        self.records.len()
    }
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    /// Removes retained records while keeping worker-local allocations reusable.
    pub fn clear(&mut self) {
        self.records.clear();
        self.bytes.clear();
        self.cigar.clear();
        self.md.clear();
        self.bismark_xm.clear();
    }

    /// Copies one already validated borrowed record into this batch.
    ///
    /// # Errors
    ///
    /// Returns a bounded allocation failure.
    pub fn push(
        &mut self,
        record: &BorrowedAlignmentRecord<'_>,
    ) -> Result<(), AlignmentRecordError> {
        let byte_count = checked_add_resource(
            checked_add_resource(
                storage_len(record.query_name.len()),
                storage_len(record.sequence.len()),
                AlignmentRecordResource::ReadBases,
            )?,
            record
                .quality
                .map_or(0, |quality| storage_len(quality.len())),
            AlignmentRecordResource::ReadBases,
        )?;
        reserve_pool(&mut self.records, 1, AlignmentRecordAllocation::Sequence)?;
        reserve_pool(
            &mut self.bytes,
            byte_count,
            AlignmentRecordAllocation::Sequence,
        )?;
        reserve_pool(
            &mut self.cigar,
            storage_len(record.cigar.len()),
            AlignmentRecordAllocation::Cigar,
        )?;
        reserve_pool(
            &mut self.md,
            record.md.map_or(0, |md| storage_len(md.len())),
            AlignmentRecordAllocation::Md,
        )?;
        reserve_pool(
            &mut self.bismark_xm,
            record.bismark_xm.map_or(0, |xm| storage_len(xm.len())),
            AlignmentRecordAllocation::MethylationCall,
        )?;
        let query_name = append_batch_bytes(&mut self.bytes, record.query_name);
        let sequence = append_batch_bytes(&mut self.bytes, record.sequence);
        let quality = record
            .quality
            .map(|quality| append_batch_bytes(&mut self.bytes, quality));
        let cigar_start = self.cigar.len();
        self.cigar.extend_from_slice(record.cigar);
        let md = record.md.map(|md| append_batch_bytes(&mut self.md, md));
        let bismark_xm = record
            .bismark_xm
            .map(|xm| append_batch_bytes(&mut self.bismark_xm, xm));
        self.records.push(BatchRecord {
            query_name,
            flag: record.flag,
            reference_ordinal: record.reference_ordinal,
            position: record.position,
            mapping_quality: record.mapping_quality,
            cigar_start,
            cigar_len: record.cigar.len(),
            mate_reference_ordinal: record.mate_reference_ordinal,
            mate_position: record.mate_position,
            template_length: record.template_length,
            sequence,
            quality,
            literal_nm: record.literal_nm,
            auxiliary_mode: record.auxiliary_mode,
            md,
            strand: record.strand,
            bismark_xm,
        });
        Ok(())
    }

    /// Visits records in deterministic insertion order.
    #[must_use]
    pub fn records(&self) -> impl ExactSizeIterator<Item = BorrowedAlignmentRecord<'_>> {
        self.records.iter().map(|record| BorrowedAlignmentRecord {
            query_name: record.query_name.slice(&self.bytes),
            flag: record.flag,
            reference_ordinal: record.reference_ordinal,
            position: record.position,
            mapping_quality: record.mapping_quality,
            cigar: &self.cigar[record.cigar_start..record.cigar_start + record.cigar_len],
            mate_reference_ordinal: record.mate_reference_ordinal,
            mate_position: record.mate_position,
            template_length: record.template_length,
            sequence: record.sequence.slice(&self.bytes),
            quality: record.quality.map(|quality| quality.slice(&self.bytes)),
            literal_nm: record.literal_nm,
            auxiliary_mode: record.auxiliary_mode,
            md: record.md.map(|md| md.slice(&self.md)),
            strand: record.strand,
            bismark_xm: record.bismark_xm.map(|xm| xm.slice(&self.bismark_xm)),
        })
    }
}

fn reserve_pool<T>(
    pool: &mut Vec<T>,
    additional: u64,
    allocation: AlignmentRecordAllocation,
) -> Result<(), AlignmentRecordError> {
    let additional =
        usize::try_from(additional).map_err(|_| AlignmentRecordError::AllocationFailed {
            allocation,
            requested: additional,
        })?;
    pool.try_reserve(additional)
        .map_err(|_| AlignmentRecordError::AllocationFailed {
            allocation,
            requested: u64::try_from(additional).unwrap_or(u64::MAX),
        })
}

fn append_batch_bytes(pool: &mut Vec<u8>, bytes: &[u8]) -> BatchBytes {
    let start = pool.len();
    pool.extend_from_slice(bytes);
    BatchBytes {
        start,
        len: bytes.len(),
    }
}

/// Integer-only placement retained while composing compact output records.
#[derive(Clone, Copy, Debug)]
pub struct AlignmentPlacement {
    reference_ordinal: u64,
    interval: ReferenceInterval,
    strand: BisulfiteStrand,
    distance: u8,
}

#[allow(missing_docs)]
impl AlignmentPlacement {
    #[must_use]
    pub const fn new(
        reference_ordinal: u64,
        interval: ReferenceInterval,
        strand: BisulfiteStrand,
        distance: u8,
    ) -> Self {
        Self {
            reference_ordinal,
            interval,
            strand,
            distance,
        }
    }
    #[must_use]
    pub const fn reference_ordinal(self) -> u64 {
        self.reference_ordinal
    }
    #[must_use]
    pub const fn interval(self) -> ReferenceInterval {
        self.interval
    }
    #[must_use]
    pub const fn strand(self) -> BisulfiteStrand {
        self.strand
    }
    #[must_use]
    pub const fn distance(self) -> u8 {
        self.distance
    }
}

/// Sequencing-order identity of one emitted record.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[allow(missing_docs)]
pub enum RecordSegment {
    Unpaired,
    First,
    Last,
}

/// Numeric MAPQ policy already resolved by the application layer.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[allow(missing_docs)]
pub enum RecordMappingQuality {
    Calibrated(u8),
    Unavailable,
    Tied,
    Unmapped,
}

impl RecordMappingQuality {
    /// Returns the exact SAM numeric value.
    #[must_use]
    pub const fn sam_value(self) -> u8 {
        match self {
            Self::Calibrated(value) => value,
            Self::Unavailable => 255,
            Self::Tied | Self::Unmapped => 0,
        }
    }
}

/// Optional alignment-field materialization policy shared by SAM and BAM records.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum AlignmentAuxiliaryMode {
    /// Emit literal NM and the bisulfite genome-conversion tag XG.
    Minimal,
    /// Emit literal NM, canonical MD, and Bismark-compatible XM/XR/XG.
    Bismark,
}

pub(crate) const fn bismark_read_conversion(strand: BisulfiteStrand) -> &'static [u8; 2] {
    match strand {
        BisulfiteStrand::OT | BisulfiteStrand::OB => b"CT",
        BisulfiteStrand::CTOT | BisulfiteStrand::CTOB => b"GA",
    }
}

pub(crate) const fn bismark_genome_conversion(strand: BisulfiteStrand) -> &'static [u8; 2] {
    match strand {
        BisulfiteStrand::OT | BisulfiteStrand::CTOT => b"CT",
        BisulfiteStrand::OB | BisulfiteStrand::CTOB => b"GA",
    }
}

/// Validated reference coordinate retained independently of an index owner.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct RecordReference {
    ordinal: u64,
    name: Box<[u8]>,
    interval: ReferenceInterval,
    position: u32,
}

#[allow(missing_docs)]
impl RecordReference {
    /// Validates and owns one reference coordinate.
    ///
    /// # Errors
    ///
    /// Returns name, range, position, or allocation failures.
    pub fn new(
        ordinal: u64,
        name: &[u8],
        reference_length: u64,
        interval: ReferenceInterval,
    ) -> Result<Self, AlignmentRecordError> {
        validate_reference_name(ordinal, name)?;
        validate_reference_length(ordinal, reference_length)?;
        if interval.end() > reference_length {
            return Err(AlignmentRecordError::ReferenceLengthOutOfRange {
                ordinal,
                length: interval.end(),
            });
        }
        let position_u64 =
            interval
                .start()
                .checked_add(1)
                .ok_or(AlignmentRecordError::FieldOutOfRange {
                    field: AlignmentRecordField::Position,
                    value: interval.start(),
                })?;
        let position =
            u32::try_from(position_u64).map_err(|_| AlignmentRecordError::FieldOutOfRange {
                field: AlignmentRecordField::Position,
                value: position_u64,
            })?;
        Ok(Self {
            ordinal,
            name: allocate_bytes_unbounded(name, AlignmentRecordAllocation::ReferenceName)?,
            interval,
            position,
        })
    }

    #[must_use]
    pub const fn ordinal(&self) -> u64 {
        self.ordinal
    }
    #[must_use]
    pub fn name(&self) -> &[u8] {
        &self.name
    }
    #[must_use]
    pub const fn interval(&self) -> ReferenceInterval {
        self.interval
    }
    #[must_use]
    pub const fn position(&self) -> u32 {
        self.position
    }
}

/// Validated mapped portion shared by SAM and BAM encoders.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct MappedAlignmentRecord {
    reference: RecordReference,
    orientation: AlignmentOrientation,
    strand: BisulfiteStrand,
    cytosine_strand: CytosineStrand,
    cigar: CoreCigar,
    literal_nm: u32,
    auxiliary_mode: AlignmentAuxiliaryMode,
    md: Option<Box<[u8]>>,
    bismark_xm: Option<Box<[u8]>>,
}

#[allow(missing_docs)]
impl MappedAlignmentRecord {
    /// Creates mapped alignment semantics from already selected values.
    ///
    /// # Errors
    ///
    /// Returns CIGAR, tag-size, or allocation failures.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        reference: RecordReference,
        orientation: AlignmentOrientation,
        strand: BisulfiteStrand,
        cytosine_strand: CytosineStrand,
        cigar: CoreCigar,
        query_length: u64,
        literal_nm: u32,
        md: Option<&[u8]>,
        bismark_xm: Option<&[u8]>,
        limits: AlignmentRecordLimits,
    ) -> Result<Self, AlignmentRecordError> {
        if cigar.is_empty() {
            return Err(AlignmentRecordError::EmptyMappedCigar);
        }
        validate_cigar(&cigar, reference.interval().len(), query_length)
            .map_err(|source| AlignmentRecordError::Cigar { source })?;
        check_limit(
            storage_len(cigar.run_count()),
            limits.max_cigar_runs(),
            AlignmentRecordResource::CigarRuns,
        )?;
        if let Some(md) = md {
            check_limit(
                storage_len(md.len()),
                limits.max_md_bytes(),
                AlignmentRecordResource::MdBytes,
            )?;
        }
        if let Some(xm) = bismark_xm {
            check_limit(
                storage_len(xm.len()),
                limits.max_optional_field_bytes(),
                AlignmentRecordResource::OptionalFieldBytes,
            )?;
        }
        let auxiliary_mode = match (md.is_some(), bismark_xm.is_some()) {
            (false, false) => AlignmentAuxiliaryMode::Minimal,
            (true, true) => AlignmentAuxiliaryMode::Bismark,
            _ => {
                return Err(AlignmentRecordError::InvalidAuxiliaryFields {
                    reason: "Bismark output requires both MD and XM",
                });
            }
        };
        let mut optional_bytes = checked_add_resource(
            14,
            decimal_digits(u64::from(literal_nm)),
            AlignmentRecordResource::OptionalFieldBytes,
        )?;
        if let (Some(md), Some(bismark_xm)) = (md, bismark_xm) {
            optional_bytes = checked_add_resource(
                optional_bytes,
                6 + storage_len(md.len()),
                AlignmentRecordResource::OptionalFieldBytes,
            )?;
            optional_bytes = checked_add_resource(
                optional_bytes,
                14 + storage_len(bismark_xm.len()),
                AlignmentRecordResource::OptionalFieldBytes,
            )?;
        }
        check_limit(
            optional_bytes,
            limits.max_optional_field_bytes(),
            AlignmentRecordResource::OptionalFieldBytes,
        )?;
        Ok(Self {
            reference,
            orientation,
            strand,
            cytosine_strand,
            cigar,
            literal_nm,
            auxiliary_mode,
            md: md
                .map(|bytes| allocate_bytes_unbounded(bytes, AlignmentRecordAllocation::Md))
                .transpose()?,
            bismark_xm: bismark_xm
                .map(|bytes| {
                    allocate_bytes_unbounded(bytes, AlignmentRecordAllocation::MethylationCall)
                })
                .transpose()?,
        })
    }

    #[must_use]
    pub const fn reference(&self) -> &RecordReference {
        &self.reference
    }
    #[must_use]
    pub const fn orientation(&self) -> AlignmentOrientation {
        self.orientation
    }
    #[must_use]
    pub const fn strand(&self) -> BisulfiteStrand {
        self.strand
    }
    #[must_use]
    pub const fn cytosine_strand(&self) -> CytosineStrand {
        self.cytosine_strand
    }
    #[must_use]
    pub const fn cigar(&self) -> &CoreCigar {
        &self.cigar
    }
    #[must_use]
    pub const fn literal_nm(&self) -> u32 {
        self.literal_nm
    }
    #[must_use]
    pub const fn auxiliary_mode(&self) -> AlignmentAuxiliaryMode {
        self.auxiliary_mode
    }
    #[must_use]
    pub fn md(&self) -> Option<&[u8]> {
        self.md.as_deref()
    }
    #[must_use]
    pub fn bismark_xm(&self) -> Option<&[u8]> {
        self.bismark_xm.as_deref()
    }
    #[must_use]
    pub const fn bismark_xr(&self) -> &'static [u8; 2] {
        bismark_read_conversion(self.strand)
    }
    #[must_use]
    pub const fn bismark_xg(&self) -> &'static [u8; 2] {
        bismark_genome_conversion(self.strand)
    }
}

/// Location and orientation of a mapped mate.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct RecordMateLocation {
    reference: RecordReference,
    orientation: AlignmentOrientation,
}

#[allow(missing_docs)]
impl RecordMateLocation {
    /// Creates a mapped mate location.
    #[must_use]
    pub const fn new(reference: RecordReference, orientation: AlignmentOrientation) -> Self {
        Self {
            reference,
            orientation,
        }
    }
    #[must_use]
    pub const fn reference(&self) -> &RecordReference {
        &self.reference
    }
    #[must_use]
    pub const fn orientation(&self) -> AlignmentOrientation {
        self.orientation
    }
}

/// One immutable validated SAM/BAM alignment record.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct AlignmentRecord {
    query_name: Box<[u8]>,
    segment: RecordSegment,
    proper_pair: bool,
    mapping_quality: RecordMappingQuality,
    mapping: Option<MappedAlignmentRecord>,
    mate: Option<RecordMateLocation>,
    template_length: i32,
    sequence: Box<[u8]>,
    quality: Option<Box<[u8]>>,
}

#[allow(missing_docs)]
impl AlignmentRecord {
    /// Creates one format-level record after mapping policy has been resolved.
    ///
    /// # Errors
    ///
    /// Returns grammar, consistency, size, or allocation failures.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        query_name: &[u8],
        segment: RecordSegment,
        proper_pair: bool,
        mapping_quality: RecordMappingQuality,
        mapping: Option<MappedAlignmentRecord>,
        mate: Option<RecordMateLocation>,
        template_length: i32,
        sequence: &[u8],
        quality: Option<&[u8]>,
        limits: AlignmentRecordLimits,
    ) -> Result<Self, AlignmentRecordError> {
        validate_query_name(query_name, limits)?;
        validate_sequence_quality(sequence, quality, limits)?;
        if let Some(mapping) = &mapping {
            validate_cigar(
                mapping.cigar(),
                mapping.reference().interval().len(),
                storage_len(sequence.len()),
            )
            .map_err(|source| AlignmentRecordError::Cigar { source })?;
        }
        Ok(Self {
            query_name: allocate_bytes_unbounded(query_name, AlignmentRecordAllocation::QueryName)?,
            segment,
            proper_pair,
            mapping_quality,
            mapping,
            mate,
            template_length,
            sequence: allocate_bytes_unbounded(sequence, AlignmentRecordAllocation::Sequence)?,
            quality: quality
                .map(|bytes| allocate_bytes_unbounded(bytes, AlignmentRecordAllocation::Quality))
                .transpose()?,
        })
    }

    #[must_use]
    pub fn query_name(&self) -> &[u8] {
        &self.query_name
    }
    #[must_use]
    pub const fn segment(&self) -> RecordSegment {
        self.segment
    }
    #[must_use]
    pub const fn is_proper_pair(&self) -> bool {
        self.proper_pair
    }
    #[must_use]
    pub const fn mapping_quality(&self) -> RecordMappingQuality {
        self.mapping_quality
    }
    #[must_use]
    pub const fn mapping(&self) -> Option<&MappedAlignmentRecord> {
        self.mapping.as_ref()
    }
    #[must_use]
    pub const fn mate(&self) -> Option<&RecordMateLocation> {
        self.mate.as_ref()
    }
    #[must_use]
    pub const fn template_length(&self) -> i32 {
        self.template_length
    }
    #[must_use]
    pub fn sequence(&self) -> &[u8] {
        &self.sequence
    }
    #[must_use]
    pub fn quality(&self) -> Option<&[u8]> {
        self.quality.as_deref()
    }
    #[must_use]
    pub const fn is_mapped(&self) -> bool {
        self.mapping.is_some()
    }
}

pub(crate) fn validate_reference_name(
    ordinal: u64,
    name: &[u8],
) -> Result<(), AlignmentRecordError> {
    let Some((&first, rest)) = name.split_first() else {
        return Err(AlignmentRecordError::InvalidReferenceNameByte {
            ordinal,
            offset: 0,
            byte: None,
        });
    };
    if !valid_reference_name_byte(first) || matches!(first, b'*' | b'=') {
        return Err(AlignmentRecordError::InvalidReferenceNameByte {
            ordinal,
            offset: 0,
            byte: Some(first),
        });
    }
    for (offset, &byte) in rest.iter().enumerate() {
        if !valid_reference_name_byte(byte) {
            return Err(AlignmentRecordError::InvalidReferenceNameByte {
                ordinal,
                offset: storage_len(offset) + 1,
                byte: Some(byte),
            });
        }
    }
    Ok(())
}

fn valid_reference_name_byte(byte: u8) -> bool {
    byte.is_ascii_graphic()
        && !matches!(
            byte,
            b'\\' | b',' | b'"' | b'\'' | b'(' | b')' | b'[' | b']' | b'{' | b'}' | b'<' | b'>'
        )
}

pub(crate) fn validate_reference_length(
    ordinal: u64,
    length: u64,
) -> Result<(), AlignmentRecordError> {
    if length == 0 || length > SAM_MAX_REFERENCE_LENGTH {
        Err(AlignmentRecordError::ReferenceLengthOutOfRange { ordinal, length })
    } else {
        Ok(())
    }
}

pub(crate) fn validate_query_name(
    name: &[u8],
    limits: AlignmentRecordLimits,
) -> Result<(), AlignmentRecordError> {
    if name.is_empty() {
        return Err(AlignmentRecordError::EmptyQueryName);
    }
    check_limit(
        storage_len(name.len()),
        limits.max_query_name_bytes(),
        AlignmentRecordResource::QueryNameBytes,
    )?;
    for (offset, &byte) in name.iter().enumerate() {
        if !(33..=126).contains(&byte) || byte == b'@' {
            return Err(AlignmentRecordError::InvalidQueryNameByte {
                offset: storage_len(offset),
                byte,
            });
        }
    }
    Ok(())
}

pub(crate) fn validate_sequence_quality(
    sequence: &[u8],
    quality: Option<&[u8]>,
    limits: AlignmentRecordLimits,
) -> Result<(), AlignmentRecordError> {
    if sequence.is_empty() {
        return Err(AlignmentRecordError::EmptySequence);
    }
    check_limit(
        storage_len(sequence.len()),
        limits.max_read_bases(),
        AlignmentRecordResource::ReadBases,
    )?;
    if let Some(quality) = quality {
        if quality.len() != sequence.len() {
            return Err(AlignmentRecordError::QualityLengthMismatch {
                sequence: storage_len(sequence.len()),
                quality: storage_len(quality.len()),
            });
        }
        check_limit(
            storage_len(quality.len()),
            limits.max_quality_bytes(),
            AlignmentRecordResource::QualityBytes,
        )?;
        for (offset, &byte) in quality.iter().enumerate() {
            if !(33..=126).contains(&byte) {
                return Err(AlignmentRecordError::InvalidQualityByte {
                    offset: storage_len(offset),
                    byte,
                });
            }
        }
    }
    Ok(())
}

pub(crate) fn check_limit(
    observed: u64,
    limit: u64,
    resource: AlignmentRecordResource,
) -> Result<(), AlignmentRecordError> {
    if observed > limit {
        Err(AlignmentRecordError::LimitExceeded {
            resource,
            observed,
            limit,
        })
    } else {
        Ok(())
    }
}

pub(crate) fn checked_add_resource(
    current: u64,
    increment: u64,
    resource: AlignmentRecordResource,
) -> Result<u64, AlignmentRecordError> {
    current
        .checked_add(increment)
        .ok_or(AlignmentRecordError::ArithmeticOverflow {
            resource,
            current,
            increment,
        })
}

pub(crate) fn storage_len(length: usize) -> u64 {
    u64::try_from(length).unwrap_or(u64::MAX)
}

pub(crate) fn allocate_bytes_unbounded(
    bytes: &[u8],
    allocation: AlignmentRecordAllocation,
) -> Result<Box<[u8]>, AlignmentRecordError> {
    let mut output = Vec::new();
    output
        .try_reserve_exact(bytes.len())
        .map_err(|_| AlignmentRecordError::AllocationFailed {
            allocation,
            requested: storage_len(bytes.len()),
        })?;
    output.extend_from_slice(bytes);
    Ok(output.into_boxed_slice())
}

pub(crate) fn decimal_digits(mut value: u64) -> u64 {
    let mut digits = 1;
    while value >= 10 {
        value /= 10;
        digits += 1;
    }
    digits
}

pub(crate) fn append_u64(output: &mut Vec<u8>, value: u64) {
    let mut buffer = [0_u8; 20];
    let mut cursor = buffer.len();
    let mut remaining = value;
    loop {
        cursor -= 1;
        buffer[cursor] = b'0' + u8::try_from(remaining % 10).expect("digit");
        remaining /= 10;
        if remaining == 0 {
            break;
        }
    }
    output.extend_from_slice(&buffer[cursor..]);
}

pub(crate) fn cigar_text_length(cigar: &CoreCigar) -> Result<u64, AlignmentRecordError> {
    let mut length = 0;
    for run in cigar.runs() {
        length = checked_add_resource(
            length,
            decimal_digits(run.length()) + 1,
            AlignmentRecordResource::CigarTextBytes,
        )?;
    }
    Ok(length)
}

const _: () = assert!(size_of::<usize>() <= size_of::<u64>());

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn format_model_rejects_empty_names() {
        assert!(matches!(
            RecordReference::new(
                0,
                b"",
                10,
                ReferenceInterval::new(0, 1, bsbit_core::coordinate::ReferenceLength::new(10),)
                    .unwrap(),
            ),
            Err(AlignmentRecordError::InvalidReferenceNameByte { byte: None, .. })
        ));
    }

    #[test]
    fn direct_cigar_requires_positive_runs() {
        assert!(AlignmentCigarRun::new(AlignmentCigarOp::Match, 0).is_err());
    }

    #[test]
    fn mapped_record_derives_one_consistent_auxiliary_mode() {
        fn mapped(
            md: Option<&[u8]>,
            xm: Option<&[u8]>,
        ) -> Result<MappedAlignmentRecord, AlignmentRecordError> {
            let length = bsbit_core::coordinate::ReferenceLength::new(8);
            let interval = ReferenceInterval::new(0, 4, length).expect("interval");
            let reference =
                RecordReference::new(0, b"chr1", 8, interval).expect("reference coordinate");
            MappedAlignmentRecord::new(
                reference,
                AlignmentOrientation::Forward,
                BisulfiteStrand::OT,
                CytosineStrand::Top,
                CoreCigar::all_matches(4),
                4,
                0,
                md,
                xm,
                AlignmentRecordLimits::default(),
            )
        }

        assert_eq!(
            mapped(None, None)
                .expect("minimal mapping")
                .auxiliary_mode(),
            AlignmentAuxiliaryMode::Minimal
        );
        assert_eq!(
            mapped(Some(b"4"), Some(b"...."))
                .expect("Bismark mapping")
                .auxiliary_mode(),
            AlignmentAuxiliaryMode::Bismark
        );
        assert!(matches!(
            mapped(Some(b"4"), None),
            Err(AlignmentRecordError::InvalidAuxiliaryFields { .. })
        ));
        assert!(matches!(
            mapped(None, Some(b"....")),
            Err(AlignmentRecordError::InvalidAuxiliaryFields { .. })
        ));
    }

    #[test]
    fn compact_record_batch_uses_reusable_pools_and_preserves_views() {
        let cigar = [AlignmentCigarRun::new(AlignmentCigarOp::Match, 4).expect("CIGAR")];
        let record = BorrowedAlignmentRecord::new(
            b"read",
            0,
            Some(0),
            1,
            60,
            &cigar,
            None,
            0,
            0,
            b"ACGT",
            Some(b"IIII"),
            0,
            AlignmentAuxiliaryMode::Minimal,
            None,
            BisulfiteStrand::OT,
            None,
            AlignmentRecordLimits::default(),
        )
        .expect("compact record");
        let mut batch = AlignmentRecordBatch::new();
        batch.push(&record).expect("first pooled push");
        let capacities = (
            batch.records.capacity(),
            batch.bytes.capacity(),
            batch.cigar.capacity(),
            batch.md.capacity(),
            batch.bismark_xm.capacity(),
        );
        let retained = batch.records().next().expect("retained view");
        assert_eq!(retained.query_name(), b"read");
        assert_eq!(retained.sequence(), b"ACGT");
        assert_eq!(retained.cigar(), &cigar);
        batch.clear();
        batch.push(&record).expect("reused pooled push");
        assert_eq!(
            capacities,
            (
                batch.records.capacity(),
                batch.bytes.capacity(),
                batch.cigar.capacity(),
                batch.md.capacity(),
                batch.bismark_xm.capacity(),
            )
        );
    }

    #[test]
    #[allow(clippy::items_after_statements)]
    fn compact_record_rejects_cross_codec_auxiliary_states() {
        let cigar = [AlignmentCigarRun::new(AlignmentCigarOp::Match, 4).expect("CIGAR")];
        fn compact<'a>(
            reference_ordinal: Option<u64>,
            cigar: &'a [AlignmentCigarRun],
            mode: AlignmentAuxiliaryMode,
            md: Option<&'a [u8]>,
            xm: Option<&'a [u8]>,
        ) -> Result<BorrowedAlignmentRecord<'a>, AlignmentRecordError> {
            BorrowedAlignmentRecord::new(
                b"read",
                0,
                reference_ordinal,
                1,
                60,
                cigar,
                None,
                0,
                0,
                b"ACGT",
                Some(b"IIII"),
                0,
                mode,
                md,
                BisulfiteStrand::OT,
                xm,
                AlignmentRecordLimits::default(),
            )
        }
        assert!(matches!(
            compact(None, &cigar, AlignmentAuxiliaryMode::Minimal, None, None),
            Err(AlignmentRecordError::InvalidCompactRecord { .. })
        ));
        assert!(matches!(
            compact(
                Some(0),
                &cigar,
                AlignmentAuxiliaryMode::Bismark,
                None,
                Some(b"....")
            ),
            Err(AlignmentRecordError::InvalidCompactRecord { .. })
        ));
        assert!(matches!(
            compact(
                Some(0),
                &cigar,
                AlignmentAuxiliaryMode::Minimal,
                Some(b"4"),
                None
            ),
            Err(AlignmentRecordError::InvalidCompactRecord { .. })
        ));
    }

    #[test]
    fn compact_minimal_auxiliary_budget_is_exactly_nm_and_xg() {
        let cigar = [AlignmentCigarRun::new(AlignmentCigarOp::Match, 4).expect("CIGAR")];
        let limits =
            AlignmentRecordLimits::new(254, 100, 100, 10, 100, 100, 15, 100, 10, 100, 1_000);
        BorrowedAlignmentRecord::new(
            b"read",
            0,
            Some(0),
            1,
            60,
            &cigar,
            None,
            0,
            0,
            b"ACGT",
            Some(b"IIII"),
            0,
            AlignmentAuxiliaryMode::Minimal,
            None,
            BisulfiteStrand::OT,
            None,
            limits,
        )
        .expect("14 fixed bytes plus one NM digit fits exactly");
    }
}
