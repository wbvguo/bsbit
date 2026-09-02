//! BAM field decoding, indexing, and indexed region reading.

use core::fmt;
use std::fs::File;
use std::path::{Path, PathBuf};

use bsbit_io::FileIdentity;

use crate::htslib::{
    HtsError, HtsErrorKind, HtsOperation, absolute_path, io_error, native_error, path_cstring,
    simple_error, validate_reader_path,
};
use crate::sys::{self, NativeIndexedBamReader, NativeIndexedBamRecordView, NativeStatus};

const BAM_BASES: &[u8; 16] = b"=ACMGRSVTWYHKDBN";

/// One semantic BAM CIGAR operation.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum BamCigarOperation {
    /// Alignment match or mismatch (`M`).
    AlignmentMatch,
    /// Query insertion (`I`).
    Insertion,
    /// Reference deletion (`D`).
    Deletion,
    /// Reference skip (`N`).
    ReferenceSkip,
    /// Soft clipping (`S`).
    SoftClip,
    /// Hard clipping (`H`).
    HardClip,
    /// Padding (`P`).
    Padding,
    /// Exact sequence match (`=`).
    SequenceMatch,
    /// Sequence mismatch (`X`).
    SequenceMismatch,
    /// Deprecated back operation (`B`).
    Back,
}

/// One validated BAM CIGAR run.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct BamCigarRun {
    length: u32,
    operation: BamCigarOperation,
}

/// One reference-consuming alignment column projected from BAM CIGAR, SEQ,
/// and QUAL.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct BamAlignmentColumn {
    /// Zero-based reference position.
    pub position: u32,
    /// Forward-reference base supplied by an external reference, or zero in
    /// the raw BAM projection.
    pub reference_base: u8,
    /// Query base, or `None` for a deletion.
    pub query_base: Option<u8>,
    /// Query Phred quality, or `None` for a deletion or unavailable quality.
    pub query_quality: Option<u8>,
}

/// Reusable storage for decoding and reconstructing one BAM record.
#[derive(Debug, Default)]
pub struct BamRecordDecodeWorkspace {
    sequence: Vec<u8>,
    cigar: Vec<BamCigarRun>,
    columns: Vec<BamAlignmentColumn>,
}

impl BamCigarRun {
    /// Constructs a nonempty semantic CIGAR run.
    #[must_use]
    pub const fn new(length: u32, operation: BamCigarOperation) -> Option<Self> {
        if length == 0 {
            None
        } else {
            Some(Self { length, operation })
        }
    }

    /// Returns the nonzero run length.
    #[must_use]
    pub const fn length(self) -> u32 {
        self.length
    }

    /// Returns the semantic operation.
    #[must_use]
    pub const fn operation(self) -> BamCigarOperation {
        self.operation
    }
}

/// A malformed copied BAM record field.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct BamRecordFieldError {
    message: String,
}

impl BamRecordFieldError {
    fn new(message: impl Into<String>) -> Self {
        Self {
            message: message.into(),
        }
    }
}

impl fmt::Display for BamRecordFieldError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(&self.message)
    }
}

impl std::error::Error for BamRecordFieldError {}

impl IndexedBamRecord {
    /// Decodes the four-bit BAM sequence into a caller-owned reusable buffer.
    ///
    /// Only canonical bases and `N` are accepted by the canonical bsbit BAM
    /// contract.
    ///
    /// # Errors
    ///
    /// Returns an allocation or unsupported-base-code error.
    pub fn decode_sequence_into(
        &self,
        destination: &mut Vec<u8>,
    ) -> Result<(), BamRecordFieldError> {
        destination.clear();
        destination.try_reserve(self.sequence_length).map_err(|_| {
            BamRecordFieldError::new(format!(
                "could not reserve {} decoded sequence bases",
                self.sequence_length
            ))
        })?;
        for index in 0..self.sequence_length {
            let code = self.sequence_code(index).ok_or_else(|| {
                BamRecordFieldError::new("packed BAM sequence is shorter than its declared length")
            })?;
            let base = BAM_BASES[usize::from(code)];
            if !matches!(base, b'A' | b'C' | b'G' | b'T' | b'N') {
                return Err(BamRecordFieldError::new(format!(
                    "bsbit BAM sequence contains unsupported base code {code} at offset {index}"
                )));
            }
            destination.push(base);
        }
        Ok(())
    }

    /// Decodes and validates BAM CIGAR words into a reusable buffer.
    ///
    /// # Errors
    ///
    /// Returns an allocation, zero-length, or unknown-operation error.
    pub fn decode_cigar_into(
        &self,
        destination: &mut Vec<BamCigarRun>,
    ) -> Result<(), BamRecordFieldError> {
        destination.clear();
        destination.try_reserve(self.cigar.len()).map_err(|_| {
            BamRecordFieldError::new(format!(
                "could not reserve {} decoded CIGAR runs",
                self.cigar.len()
            ))
        })?;
        for &encoded in &self.cigar {
            let length = encoded >> 4;
            if length == 0 {
                return Err(BamRecordFieldError::new(format!(
                    "invalid zero-length BAM CIGAR word {encoded}"
                )));
            }
            let operation = match encoded & 0x0f {
                0 => BamCigarOperation::AlignmentMatch,
                1 => BamCigarOperation::Insertion,
                2 => BamCigarOperation::Deletion,
                3 => BamCigarOperation::ReferenceSkip,
                4 => BamCigarOperation::SoftClip,
                5 => BamCigarOperation::HardClip,
                6 => BamCigarOperation::Padding,
                7 => BamCigarOperation::SequenceMatch,
                8 => BamCigarOperation::SequenceMismatch,
                9 => BamCigarOperation::Back,
                operation => {
                    return Err(BamRecordFieldError::new(format!(
                        "unknown BAM CIGAR operation {operation}"
                    )));
                }
            };
            destination.push(BamCigarRun { length, operation });
        }
        Ok(())
    }

    /// Returns one string-valued (`Z`) auxiliary field without copying.
    ///
    /// The full auxiliary payload is validated while searching, so malformed
    /// fields cannot be hidden behind an earlier matching tag.
    ///
    /// # Errors
    ///
    /// Returns a truncated, duplicate, mistyped, or unsupported-field error.
    pub fn string_auxiliary(
        &self,
        requested_tag: [u8; 2],
    ) -> Result<Option<&[u8]>, BamRecordFieldError> {
        let [value] = self.string_auxiliaries_impl([requested_tag])?;
        Ok(value)
    }

    fn string_auxiliaries_impl<const N: usize>(
        &self,
        requested_tags: [[u8; 2]; N],
    ) -> Result<[Option<&[u8]>; N], BamRecordFieldError> {
        let mut cursor = ByteCursor::new(&self.auxiliary);
        let mut found = [None; N];
        while !cursor.rest().is_empty() {
            let tag = cursor.take(2, "auxiliary tag")?;
            let physical_type = cursor.byte("auxiliary type")?;
            let requested_index = requested_tags
                .iter()
                .position(|requested_tag| tag == requested_tag);
            if let Some(index) = requested_index
                && physical_type != b'Z'
            {
                let requested_tag = requested_tags[index];
                return Err(BamRecordFieldError::new(format!(
                    "{}{} auxiliary tag is not type Z",
                    char::from(requested_tag[0]),
                    char::from(requested_tag[1])
                )));
            }
            match physical_type {
                b'A' | b'c' | b'C' => {
                    cursor.take(1, "one-byte auxiliary value")?;
                }
                b's' | b'S' => {
                    cursor.take(2, "two-byte auxiliary value")?;
                }
                b'i' | b'I' | b'f' => {
                    cursor.take(4, "four-byte auxiliary value")?;
                }
                b'd' => {
                    cursor.take(8, "eight-byte auxiliary value")?;
                }
                b'Z' | b'H' => {
                    let length = cursor
                        .rest()
                        .iter()
                        .position(|byte| *byte == 0)
                        .ok_or_else(|| {
                            BamRecordFieldError::new("unterminated string auxiliary value")
                        })?;
                    let value = cursor.take(length, "string auxiliary value")?;
                    cursor.take(1, "string auxiliary terminator")?;
                    if let Some(index) = requested_index
                        && physical_type == b'Z'
                        && found[index].replace(value).is_some()
                    {
                        let requested_tag = requested_tags[index];
                        return Err(BamRecordFieldError::new(format!(
                            "record contains duplicate {}{} tags",
                            char::from(requested_tag[0]),
                            char::from(requested_tag[1])
                        )));
                    }
                }
                b'B' => {
                    let subtype = cursor.byte("array auxiliary subtype")?;
                    let element_bytes = match subtype {
                        b'c' | b'C' => 1_usize,
                        b's' | b'S' => 2,
                        b'i' | b'I' | b'f' => 4,
                        _ => {
                            return Err(BamRecordFieldError::new(format!(
                                "unsupported BAM array subtype {subtype}"
                            )));
                        }
                    };
                    let count = cursor.i32("array auxiliary count")?;
                    let count = usize::try_from(count).map_err(|_| {
                        BamRecordFieldError::new(format!(
                            "array auxiliary count is negative ({count})"
                        ))
                    })?;
                    let bytes = count.checked_mul(element_bytes).ok_or_else(|| {
                        BamRecordFieldError::new("array auxiliary byte length overflowed")
                    })?;
                    cursor.take(bytes, "array auxiliary values")?;
                }
                value => {
                    return Err(BamRecordFieldError::new(format!(
                        "unsupported BAM auxiliary type {value}"
                    )));
                }
            }
        }
        Ok(found)
    }

    /// Projects CIGAR, SEQ, and QUAL into reference-coordinate columns without
    /// reading `MD:Z` or any other auxiliary field.
    ///
    /// Every returned `reference_base` is zero so an authoritative external
    /// reference can supply it. Query bases, qualities, and deletion columns
    /// are preserved.
    ///
    /// # Errors
    ///
    /// Returns a malformed-field, allocation, coordinate, or CIGAR error.
    pub fn project_alignment_into<'a>(
        &self,
        reference_length: u32,
        workspace: &'a mut BamRecordDecodeWorkspace,
    ) -> Result<&'a [BamAlignmentColumn], BamRecordFieldError> {
        let start = self.decode_alignment_inputs_into(workspace)?;
        project_decoded_alignment_into(
            start,
            reference_length,
            &workspace.sequence,
            &self.quality,
            &workspace.cigar,
            &mut workspace.columns,
        )?;
        Ok(&workspace.columns)
    }

    fn decode_alignment_inputs_into(
        &self,
        workspace: &mut BamRecordDecodeWorkspace,
    ) -> Result<u32, BamRecordFieldError> {
        let start = u32::try_from(self.position).map_err(|_| {
            BamRecordFieldError::new(format!(
                "mapped reference position {} is negative or exceeds u32",
                self.position
            ))
        })?;
        if self.quality.len() != self.sequence_length {
            return Err(BamRecordFieldError::new(format!(
                "BAM QUAL has {} values but SEQ has {} bases",
                self.quality.len(),
                self.sequence_length
            )));
        }
        self.decode_sequence_into(&mut workspace.sequence)?;
        self.decode_cigar_into(&mut workspace.cigar)?;
        if workspace.cigar.is_empty() {
            return Err(BamRecordFieldError::new(
                "mapped BAM record has no CIGAR runs",
            ));
        }
        Ok(start)
    }
}

fn project_decoded_alignment_into(
    start: u32,
    reference_length: u32,
    sequence: &[u8],
    qualities: &[u8],
    cigar: &[BamCigarRun],
    columns: &mut Vec<BamAlignmentColumn>,
) -> Result<(), BamRecordFieldError> {
    let reference_bases = cigar.iter().try_fold(0_u64, |total, run| {
        if matches!(
            run.operation(),
            BamCigarOperation::AlignmentMatch
                | BamCigarOperation::Deletion
                | BamCigarOperation::SequenceMatch
                | BamCigarOperation::SequenceMismatch
        ) {
            total
                .checked_add(u64::from(run.length()))
                .ok_or_else(|| BamRecordFieldError::new("CIGAR reference length overflowed u64"))
        } else {
            Ok(total)
        }
    })?;
    let capacity = usize::try_from(reference_bases)
        .map_err(|_| BamRecordFieldError::new("CIGAR reference length is not addressable"))?;
    columns.clear();
    columns.try_reserve(capacity).map_err(|_| {
        BamRecordFieldError::new(format!("could not reserve {capacity} alignment columns"))
    })?;
    let mut query_index = 0_usize;
    let mut reference_position = u64::from(start);
    for run in cigar {
        let length = usize::try_from(run.length()).expect("u32 fits usize on supported Linux");
        match run.operation() {
            BamCigarOperation::AlignmentMatch
            | BamCigarOperation::SequenceMatch
            | BamCigarOperation::SequenceMismatch => {
                let query_end = query_index
                    .checked_add(length)
                    .ok_or_else(|| BamRecordFieldError::new("CIGAR query length overflowed"))?;
                let query = sequence
                    .get(query_index..query_end)
                    .ok_or_else(|| BamRecordFieldError::new("CIGAR consumes beyond BAM SEQ"))?;
                let query_qualities = qualities
                    .get(query_index..query_end)
                    .ok_or_else(|| BamRecordFieldError::new("CIGAR consumes beyond BAM QUAL"))?;
                for (base, quality) in query.iter().zip(query_qualities) {
                    columns.push(BamAlignmentColumn {
                        position: checked_reference_position(reference_position, reference_length)?,
                        reference_base: 0,
                        query_base: Some(*base),
                        query_quality: (*quality != u8::MAX).then_some(*quality),
                    });
                    reference_position += 1;
                }
                query_index = query_end;
            }
            BamCigarOperation::Insertion | BamCigarOperation::SoftClip => {
                query_index = query_index
                    .checked_add(length)
                    .filter(|index| *index <= sequence.len())
                    .ok_or_else(|| BamRecordFieldError::new("CIGAR consumes beyond BAM SEQ"))?;
            }
            BamCigarOperation::Deletion => {
                for _ in 0..length {
                    columns.push(BamAlignmentColumn {
                        position: checked_reference_position(reference_position, reference_length)?,
                        reference_base: 0,
                        query_base: None,
                        query_quality: None,
                    });
                    reference_position += 1;
                }
            }
            BamCigarOperation::HardClip | BamCigarOperation::Padding => {}
            BamCigarOperation::ReferenceSkip => {
                return Err(BamRecordFieldError::new(
                    "reference-skip CIGAR operations are outside the bsbit BAM contract",
                ));
            }
            BamCigarOperation::Back => {
                return Err(BamRecordFieldError::new(
                    "back CIGAR operations are unsupported",
                ));
            }
        }
    }
    if query_index != sequence.len() {
        return Err(BamRecordFieldError::new(format!(
            "CIGAR consumes {query_index} query bases but BAM SEQ has {}",
            sequence.len()
        )));
    }
    Ok(())
}

fn checked_reference_position(position: u64, length: u32) -> Result<u32, BamRecordFieldError> {
    if position >= u64::from(length) {
        return Err(BamRecordFieldError::new(format!(
            "alignment reference position {position} reaches beyond contig length {length}"
        )));
    }
    u32::try_from(position)
        .map_err(|_| BamRecordFieldError::new("reference position does not fit u32"))
}

struct ByteCursor<'a> {
    bytes: &'a [u8],
    offset: usize,
}

impl<'a> ByteCursor<'a> {
    const fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, offset: 0 }
    }

    fn rest(&self) -> &'a [u8] {
        &self.bytes[self.offset..]
    }

    fn take(&mut self, length: usize, label: &str) -> Result<&'a [u8], BamRecordFieldError> {
        let end = self
            .offset
            .checked_add(length)
            .ok_or_else(|| BamRecordFieldError::new(format!("{label} offset overflowed")))?;
        let value = self
            .bytes
            .get(self.offset..end)
            .ok_or_else(|| BamRecordFieldError::new(format!("truncated {label}")))?;
        self.offset = end;
        Ok(value)
    }

    fn byte(&mut self, label: &str) -> Result<u8, BamRecordFieldError> {
        Ok(self.take(1, label)?[0])
    }

    fn i32(&mut self, label: &str) -> Result<i32, BamRecordFieldError> {
        Ok(i32::from_le_bytes(
            self.take(4, label)?
                .try_into()
                .expect("four-byte cursor slice"),
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::{BamRecordDecodeWorkspace, IndexedBamRecord};

    fn synthetic_record(auxiliary: Vec<u8>) -> IndexedBamRecord {
        IndexedBamRecord {
            reference_id: 0,
            position: 0,
            mapping_quality: 60,
            flag: 0,
            mate_reference_id: -1,
            mate_position: -1,
            template_length: 0,
            query_name: b"synthetic".to_vec(),
            cigar: vec![8_u32 << 4],
            packed_sequence: vec![0x84, 0x21, 0x48, 0x81],
            sequence_length: 8,
            quality: vec![30; 8],
            auxiliary,
        }
    }

    #[test]
    fn reference_projection_ignores_md_and_preserves_alignment_structure() {
        let record = synthetic_record(b"MDZnot-a-valid-md-value\0".to_vec());
        let mut workspace = BamRecordDecodeWorkspace::default();
        let columns = record
            .project_alignment_into(8, &mut workspace)
            .expect("reference projection ignores MD");
        assert_eq!(
            columns
                .iter()
                .map(|column| column.position)
                .collect::<Vec<_>>(),
            (0_u32..8).collect::<Vec<_>>()
        );
        assert!(columns.iter().all(|column| column.reference_base == 0));
        assert_eq!(
            columns
                .iter()
                .map(|column| column.query_base.expect("aligned query base"))
                .collect::<Vec<_>>(),
            b"TGCAGTTA"
        );
    }

    #[test]
    fn reference_projection_preserves_deletions_for_external_reference_fill() {
        let mut record = synthetic_record(Vec::new());
        record.cigar = vec![2 << 4, (1 << 4) | 2, 2 << 4];
        record.packed_sequence = vec![0x12, 0x48];
        record.sequence_length = 4;
        record.quality = vec![30; 4];

        let mut workspace = BamRecordDecodeWorkspace::default();
        let columns = record
            .project_alignment_into(5, &mut workspace)
            .expect("reference projection accepts an MD-free deletion");
        assert_eq!(columns.len(), 5);
        assert_eq!(columns[2].position, 2);
        assert_eq!(columns[2].query_base, None);
        assert_eq!(columns[2].query_quality, None);
        assert!(columns.iter().all(|column| column.reference_base == 0));
        assert_eq!(
            columns
                .iter()
                .filter_map(|column| column.query_base)
                .collect::<Vec<_>>(),
            b"ACGT"
        );
    }

    #[test]
    fn string_lookup_rejects_a_same_named_non_string_tag() {
        let record = synthetic_record(b"XGi\x01\0\0\0MDZ8\0".to_vec());
        assert!(
            record
                .string_auxiliary(*b"XG")
                .expect_err("strict string lookup rejects another physical type")
                .to_string()
                .contains("not type Z")
        );

        let hex_record = synthetic_record(b"XGH4354\0".to_vec());
        assert!(
            hex_record
                .string_auxiliary(*b"XG")
                .expect_err("hexadecimal XG is not a string")
                .to_string()
                .contains("not type Z")
        );
    }
}

/// Builds a create-only BAI at an explicit path for a coordinate-sorted BAM.
///
/// The index destination is reserved before native work, so an existing file
/// is never overwritten. The adapter removes only its own reserved path when
/// indexing or synchronization fails.
///
/// # Errors
///
/// Returns a path, exclusive-creation, HTS indexing, identity, or
/// synchronization error.
pub fn build_bam_index_create_new(
    bam: impl AsRef<Path>,
    index: impl AsRef<Path>,
    threads: u32,
) -> Result<(), HtsError> {
    let bam = absolute_path(bam.as_ref(), HtsOperation::ValidatePath)?;
    let index = absolute_path(index.as_ref(), HtsOperation::ValidatePath)?;
    let bam_c = path_cstring(&bam)?;
    let index_c = path_cstring(&index)?;
    let (reservation, identity) = bsbit_io::create_new(&index)
        .map_err(|source| io_error(HtsOperation::BuildBamIndex, &index, None, source))?;
    drop(reservation);
    if let Err(source) = sys::build_bam_index(&bam_c, &index_c, threads) {
        let _ = bsbit_io::remove_if_identity_matches(&index, identity);
        return Err(native_error(
            HtsOperation::BuildBamIndex,
            &bam,
            None,
            source,
        ));
    }
    let indexed_identity = match File::open(&index).and_then(|file| {
        let current = FileIdentity::from_file(&file)?;
        file.sync_all()?;
        Ok(current)
    }) {
        Ok(current) if current == identity => current,
        Ok(_) => {
            return Err(simple_error(
                HtsOperation::BuildBamIndex,
                &index,
                None,
                HtsErrorKind::StagingIdentityChanged,
            ));
        }
        Err(source) => {
            let _ = bsbit_io::remove_if_identity_matches(&index, identity);
            return Err(io_error(HtsOperation::BuildBamIndex, &index, None, source));
        }
    };
    debug_assert_eq!(indexed_identity, identity);
    Ok(())
}

/// One BAM reference dictionary entry copied from an indexed reader.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct IndexedBamReference {
    name: Vec<u8>,
    length: u64,
}

impl IndexedBamReference {
    /// Returns the reference name exactly as stored in the BAM header.
    #[must_use]
    pub fn name(&self) -> &[u8] {
        &self.name
    }

    /// Returns the reference length in bases.
    #[must_use]
    pub const fn length(&self) -> u64 {
        self.length
    }
}

/// One fully copied BAM header and reference dictionary.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct IndexedBamHeader {
    pub(crate) text: Vec<u8>,
    pub(crate) references: Vec<IndexedBamReference>,
}

impl IndexedBamHeader {
    /// Returns the SAM header text.
    #[must_use]
    pub fn text(&self) -> &[u8] {
        &self.text
    }

    /// Returns reference entries in BAM numeric-id order.
    #[must_use]
    pub fn references(&self) -> &[IndexedBamReference] {
        &self.references
    }

    /// Iterates over every `SM` value declared by an `@RG` record.
    ///
    /// Values are returned exactly as encoded in the SAM header. Duplicate
    /// values are retained because deciding whether read groups belong to one
    /// or several biological samples is a caller-level policy.
    pub fn read_group_sample_names(&self) -> impl Iterator<Item = &[u8]> {
        self.text
            .split(|byte| *byte == b'\n')
            .filter_map(|line| {
                let mut fields = line.split(|byte| *byte == b'\t');
                (fields.next() == Some(b"@RG".as_slice())).then_some(fields)
            })
            .flatten()
            .filter_map(|field| field.strip_prefix(b"SM:"))
    }

    /// Returns whether the SAM header declares coordinate sort order.
    #[must_use]
    pub fn is_coordinate_sorted(&self) -> bool {
        self.text.split(|byte| *byte == b'\n').any(|line| {
            let mut fields = line.split(|byte| *byte == b'\t');
            fields.next() == Some(b"@HD".as_slice())
                && fields.any(|field| field == b"SO:coordinate")
        })
    }

    /// Returns whether one `@PG` record contains both the requested ID and
    /// program-name fields.
    #[must_use]
    pub fn has_program(&self, id: &[u8], program_name: &[u8]) -> bool {
        self.text.split(|byte| *byte == b'\n').any(|line| {
            let mut fields = line.split(|byte| *byte == b'\t');
            if fields.next() != Some(b"@PG".as_slice()) {
                return false;
            }
            let mut has_id = false;
            let mut has_program_name = false;
            for field in fields {
                has_id |= field.strip_prefix(b"ID:").is_some_and(|value| value == id);
                has_program_name |= field
                    .strip_prefix(b"PN:")
                    .is_some_and(|value| value == program_name);
            }
            has_id && has_program_name
        })
    }

    /// Returns the exact structured provenance from the unique bsbit `@PG`
    /// record.
    ///
    /// # Errors
    ///
    /// Rejects duplicate records and missing or malformed structured fields.
    pub fn bsbit_program_provenance(
        &self,
    ) -> Result<Option<crate::BsbitProgramProvenance>, crate::BsbitProgramProvenanceError> {
        crate::sam::parse_bsbit_program_provenance(&self.text)
    }
}

/// One fully owned BAM record returned by an indexed region query.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct IndexedBamRecord {
    pub(crate) reference_id: i32,
    pub(crate) position: i64,
    mapping_quality: u8,
    flag: u16,
    mate_reference_id: i32,
    mate_position: i64,
    template_length: i64,
    pub(crate) query_name: Vec<u8>,
    pub(crate) cigar: Vec<u32>,
    pub(crate) packed_sequence: Vec<u8>,
    pub(crate) sequence_length: usize,
    pub(crate) quality: Vec<u8>,
    pub(crate) auxiliary: Vec<u8>,
}

impl IndexedBamRecord {
    fn copy_from_native(&mut self, record: &NativeIndexedBamRecordView<'_>) {
        self.reference_id = record.reference_id;
        self.position = record.position;
        self.mapping_quality = record.mapping_quality;
        self.flag = record.flag;
        self.mate_reference_id = record.mate_reference_id;
        self.mate_position = record.mate_position;
        self.template_length = record.template_length;
        copy_slice_into(&mut self.query_name, record.query_name);
        copy_slice_into(&mut self.cigar, record.cigar);
        copy_slice_into(&mut self.packed_sequence, record.packed_sequence);
        self.sequence_length = record.sequence_length;
        copy_slice_into(&mut self.quality, record.quality);
        copy_slice_into(&mut self.auxiliary, record.auxiliary);
    }

    fn clear_preserving_capacity(&mut self) {
        self.reference_id = -1;
        self.position = -1;
        self.mapping_quality = 0;
        self.flag = 0;
        self.mate_reference_id = -1;
        self.mate_position = -1;
        self.template_length = 0;
        self.query_name.clear();
        self.cigar.clear();
        self.packed_sequence.clear();
        self.sequence_length = 0;
        self.quality.clear();
        self.auxiliary.clear();
    }

    /// Returns the zero-based BAM reference id, or -1 when unmapped.
    #[must_use]
    pub const fn reference_id(&self) -> i32 {
        self.reference_id
    }

    /// Returns the zero-based alignment start, or -1 when unmapped.
    #[must_use]
    pub const fn position(&self) -> i64 {
        self.position
    }

    /// Returns the numeric BAM mapping quality.
    #[must_use]
    pub const fn mapping_quality(&self) -> u8 {
        self.mapping_quality
    }

    /// Returns the complete BAM flag word.
    #[must_use]
    pub const fn flag(&self) -> u16 {
        self.flag
    }

    /// Returns the mate reference id, or -1 when absent.
    #[must_use]
    pub const fn mate_reference_id(&self) -> i32 {
        self.mate_reference_id
    }

    /// Returns the zero-based mate position, or -1 when absent.
    #[must_use]
    pub const fn mate_position(&self) -> i64 {
        self.mate_position
    }

    /// Returns the signed observed template length.
    #[must_use]
    pub const fn template_length(&self) -> i64 {
        self.template_length
    }

    /// Returns the query name without a terminating NUL.
    #[must_use]
    pub fn query_name(&self) -> &[u8] {
        &self.query_name
    }

    /// Returns packed BAM CIGAR words.
    #[must_use]
    pub fn cigar(&self) -> &[u32] {
        &self.cigar
    }

    /// Returns the BAM four-bit packed query sequence.
    #[must_use]
    pub fn packed_sequence(&self) -> &[u8] {
        &self.packed_sequence
    }

    /// Returns the number of query bases represented by the packed sequence.
    #[must_use]
    pub const fn sequence_length(&self) -> usize {
        self.sequence_length
    }

    /// Returns raw BAM base qualities; 255 means unavailable.
    #[must_use]
    pub fn quality(&self) -> &[u8] {
        &self.quality
    }

    /// Returns the raw BAM auxiliary field payload.
    #[must_use]
    pub fn auxiliary(&self) -> &[u8] {
        &self.auxiliary
    }

    /// Returns one BAM four-bit base code, if the offset is in range.
    #[must_use]
    pub fn sequence_code(&self, offset: usize) -> Option<u8> {
        if offset >= self.sequence_length {
            return None;
        }
        let packed = *self.packed_sequence.get(offset / 2)?;
        Some(if offset.is_multiple_of(2) {
            packed >> 4
        } else {
            packed & 0x0f
        })
    }
}

impl Default for IndexedBamRecord {
    fn default() -> Self {
        Self {
            reference_id: -1,
            position: -1,
            mapping_quality: 0,
            flag: 0,
            mate_reference_id: -1,
            mate_position: -1,
            template_length: 0,
            query_name: Vec::new(),
            cigar: Vec::new(),
            packed_sequence: Vec::new(),
            sequence_length: 0,
            quality: Vec::new(),
            auxiliary: Vec::new(),
        }
    }
}

fn copy_slice_into<T: Copy>(destination: &mut Vec<T>, source: &[T]) {
    destination.clear();
    destination.extend_from_slice(source);
}

/// A thread-confined BAM reader with reusable indexed region queries.
pub struct IndexedBamReader {
    path: PathBuf,
    header: IndexedBamHeader,
    native: NativeIndexedBamReader,
    record_ordinal: u64,
}

impl IndexedBamReader {
    /// Opens a local BAM and its adjacent index.
    ///
    /// # Errors
    ///
    /// Returns a path, BAM-header, or index-loading failure. The input must be
    /// a concrete local regular file.
    pub fn open(path: impl AsRef<Path>) -> Result<Self, HtsError> {
        let path = path.as_ref().to_path_buf();
        validate_reader_path(&path)?;
        let c_path = path_cstring(&path)?;
        let native = NativeIndexedBamReader::open(&c_path)
            .map_err(|source| native_error(HtsOperation::OpenIndexedBam, &path, None, source))?;
        let text = native.header_text().map_err(|source| {
            native_error(HtsOperation::ReadIndexedBamHeader, &path, None, source)
        })?;
        let native_references = native.references().map_err(|source| {
            native_error(HtsOperation::ReadIndexedBamHeader, &path, None, source)
        })?;
        let mut references = Vec::with_capacity(native_references.len());
        for reference in native_references {
            let length = u64::try_from(reference.length).map_err(|_| {
                simple_error(
                    HtsOperation::ReadIndexedBamHeader,
                    &path,
                    None,
                    HtsErrorKind::Native(NativeStatus::HeaderFailed),
                )
            })?;
            references.push(IndexedBamReference {
                name: reference.name,
                length,
            });
        }
        Ok(Self {
            path,
            header: IndexedBamHeader { text, references },
            native,
            record_ordinal: 0,
        })
    }

    /// Returns the copied BAM header.
    #[must_use]
    pub const fn header(&self) -> &IndexedBamHeader {
        &self.header
    }

    /// Selects one zero-based, half-open interval for subsequent reads.
    ///
    /// # Errors
    ///
    /// Returns an argument error for an unknown reference or an empty/out-of-
    /// range interval, or an HTS index-query failure.
    pub fn query(&mut self, reference_id: u32, start: u64, end: u64) -> Result<(), HtsError> {
        let reference_ordinal = usize::try_from(reference_id).map_err(|_| {
            simple_error(
                HtsOperation::QueryIndexedBam,
                &self.path,
                None,
                HtsErrorKind::Native(NativeStatus::InvalidArgument),
            )
        })?;
        let Some(reference) = self.header.references.get(reference_ordinal) else {
            return Err(simple_error(
                HtsOperation::QueryIndexedBam,
                &self.path,
                None,
                HtsErrorKind::Native(NativeStatus::InvalidArgument),
            ));
        };
        if start >= end || end > reference.length {
            return Err(simple_error(
                HtsOperation::QueryIndexedBam,
                &self.path,
                None,
                HtsErrorKind::Native(NativeStatus::InvalidArgument),
            ));
        }
        let reference_id = i32::try_from(reference_id).map_err(|_| {
            simple_error(
                HtsOperation::QueryIndexedBam,
                &self.path,
                None,
                HtsErrorKind::Native(NativeStatus::InvalidArgument),
            )
        })?;
        let start = i64::try_from(start).map_err(|_| {
            simple_error(
                HtsOperation::QueryIndexedBam,
                &self.path,
                None,
                HtsErrorKind::Native(NativeStatus::InvalidArgument),
            )
        })?;
        let end = i64::try_from(end).map_err(|_| {
            simple_error(
                HtsOperation::QueryIndexedBam,
                &self.path,
                None,
                HtsErrorKind::Native(NativeStatus::InvalidArgument),
            )
        })?;
        self.native
            .query(reference_id, start, end)
            .map_err(|source| {
                native_error(HtsOperation::QueryIndexedBam, &self.path, None, source)
            })?;
        self.record_ordinal = 0;
        Ok(())
    }

    /// Copies the next BAM record overlapping the selected interval.
    ///
    /// # Errors
    ///
    /// Returns a terminal native decoding failure with the region-local record
    /// ordinal attached.
    pub fn next_record(&mut self) -> Result<Option<IndexedBamRecord>, HtsError> {
        let mut record = IndexedBamRecord::default();
        self.next_record_into(&mut record)
            .map(|has_record| has_record.then_some(record))
    }

    /// Copies the next BAM record into caller-owned reusable storage.
    ///
    /// The variable-length QNAME, CIGAR, packed sequence, quality, and
    /// auxiliary buffers retain their capacity across records and region
    /// queries. At end of interval, the record is cleared without releasing
    /// those allocations.
    ///
    /// # Errors
    ///
    /// Returns a terminal native decoding failure with the region-local record
    /// ordinal attached.
    pub fn next_record_into(&mut self, record: &mut IndexedBamRecord) -> Result<bool, HtsError> {
        let next_ordinal = self.record_ordinal.checked_add(1).ok_or_else(|| {
            simple_error(
                HtsOperation::ReadIndexedBamRecord,
                &self.path,
                None,
                HtsErrorKind::RecordCountOverflow,
            )
        })?;
        let next = self.native.next_record().map_err(|source| {
            native_error(
                HtsOperation::ReadIndexedBamRecord,
                &self.path,
                Some(next_ordinal),
                source,
            )
        })?;
        let Some(native_record) = next else {
            record.clear_preserving_capacity();
            return Ok(false);
        };
        record.copy_from_native(&native_record);
        self.record_ordinal = next_ordinal;
        Ok(true)
    }

    /// Explicitly closes all native BAM and index resources.
    ///
    /// # Errors
    ///
    /// Returns a copied native close failure.
    pub fn close(mut self) -> Result<(), HtsError> {
        self.native
            .close()
            .map_err(|source| native_error(HtsOperation::CloseIndexedBam, &self.path, None, source))
    }
}
