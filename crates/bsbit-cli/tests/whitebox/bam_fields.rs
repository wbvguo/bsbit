//! Independent BAM-layout read-back for every canonical alignment field.

use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use super::record_fixture::single_fixture;
use super::{
    PairedRecordComposer, build_sam_header, build_single_alignment_record_with_auxiliary_mode,
};
use bsbit_core::sequence::{NormalizedSequence, normalize_dna};
use bsbit_hts::{
    AlignmentAuxiliaryMode, AlignmentRead, AlignmentRecord, AlignmentRecordBatch,
    AlignmentRecordLimits, BamStagingWriter, BorrowedAlignmentRead, Compression, DecodedReader,
    SamHeader, sam_flag, sam_header_bytes,
};
use bsbit_index::reference::{ContigInput, ReferenceBuildLimits, ReferenceIndex};

#[derive(Debug, Eq, PartialEq)]
struct DecodeError {
    offset: usize,
    context: &'static str,
}

#[derive(Debug, Eq, PartialEq)]
struct DecodedReference {
    name: Vec<u8>,
    length: i32,
}

#[derive(Debug, Eq, PartialEq)]
enum AuxValue {
    Integer(i64),
    String(Vec<u8>),
}

#[derive(Debug, Eq, PartialEq)]
struct DecodedAux {
    tag: [u8; 2],
    value: AuxValue,
}

#[derive(Debug, Eq, PartialEq)]
struct DecodedRecord {
    reference_id: i32,
    position: i32,
    read_name: Vec<u8>,
    mapping_quality: u8,
    bin: u16,
    cigar: Vec<(u32, u8)>,
    flag: u16,
    mate_reference_id: i32,
    mate_position: i32,
    template_length: i32,
    sequence: Vec<u8>,
    quality: Option<Vec<u8>>,
    aux: Vec<DecodedAux>,
}

#[derive(Debug, Eq, PartialEq)]
struct DecodedBam {
    header_text: Vec<u8>,
    references: Vec<DecodedReference>,
    records: Vec<DecodedRecord>,
}

struct BamCursor<'a> {
    bytes: &'a [u8],
    offset: usize,
}

impl<'a> BamCursor<'a> {
    const fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, offset: 0 }
    }

    const fn remaining(&self) -> usize {
        self.bytes.len() - self.offset
    }

    fn fail(&self, context: &'static str) -> DecodeError {
        DecodeError {
            offset: self.offset,
            context,
        }
    }

    fn take(&mut self, length: usize, context: &'static str) -> Result<&'a [u8], DecodeError> {
        let end = self
            .offset
            .checked_add(length)
            .ok_or_else(|| self.fail(context))?;
        let value = self
            .bytes
            .get(self.offset..end)
            .ok_or_else(|| self.fail(context))?;
        self.offset = end;
        Ok(value)
    }

    fn u8(&mut self, context: &'static str) -> Result<u8, DecodeError> {
        Ok(self.take(1, context)?[0])
    }

    fn u16(&mut self, context: &'static str) -> Result<u16, DecodeError> {
        let bytes: [u8; 2] = self
            .take(2, context)?
            .try_into()
            .expect("checked two-byte slice");
        Ok(u16::from_le_bytes(bytes))
    }

    fn i16(&mut self, context: &'static str) -> Result<i16, DecodeError> {
        let bytes: [u8; 2] = self
            .take(2, context)?
            .try_into()
            .expect("checked two-byte slice");
        Ok(i16::from_le_bytes(bytes))
    }

    fn u32(&mut self, context: &'static str) -> Result<u32, DecodeError> {
        let bytes: [u8; 4] = self
            .take(4, context)?
            .try_into()
            .expect("checked four-byte slice");
        Ok(u32::from_le_bytes(bytes))
    }

    fn i32(&mut self, context: &'static str) -> Result<i32, DecodeError> {
        let bytes: [u8; 4] = self
            .take(4, context)?
            .try_into()
            .expect("checked four-byte slice");
        Ok(i32::from_le_bytes(bytes))
    }

    fn nonnegative_length(&mut self, context: &'static str) -> Result<usize, DecodeError> {
        usize::try_from(self.i32(context)?).map_err(|_| self.fail(context))
    }

    fn nul_terminated(&mut self, context: &'static str) -> Result<Vec<u8>, DecodeError> {
        let Some(length) = self.bytes[self.offset..].iter().position(|byte| *byte == 0) else {
            return Err(self.fail(context));
        };
        let value = self.take(length, context)?.to_vec();
        self.take(1, context)?;
        Ok(value)
    }
}

fn decode_bam(bytes: &[u8]) -> Result<DecodedBam, DecodeError> {
    let mut cursor = BamCursor::new(bytes);
    if cursor.take(4, "BAM magic")? != b"BAM\x01" {
        return Err(cursor.fail("BAM magic"));
    }
    let header_length = cursor.nonnegative_length("header length")?;
    let header_text = cursor.take(header_length, "header text")?.to_vec();
    let reference_count = cursor.nonnegative_length("reference count")?;
    if reference_count > cursor.remaining() / 8 {
        return Err(cursor.fail("reference count"));
    }
    let mut references = Vec::with_capacity(reference_count);
    for _ in 0..reference_count {
        let name_length = cursor.nonnegative_length("reference name length")?;
        let name_bytes = cursor.take(name_length, "reference name")?;
        let Some((&0, name)) = name_bytes.split_last() else {
            return Err(cursor.fail("reference name terminator"));
        };
        if name.is_empty() || name.contains(&0) {
            return Err(cursor.fail("reference name"));
        }
        let length = cursor.i32("reference length")?;
        if length < 0 {
            return Err(cursor.fail("reference length"));
        }
        references.push(DecodedReference {
            name: name.to_vec(),
            length,
        });
    }

    let mut records = Vec::new();
    while cursor.remaining() != 0 {
        let block_length = cursor.nonnegative_length("record block length")?;
        let block = cursor.take(block_length, "record block")?;
        records.push(decode_record(block, references.len())?);
    }
    Ok(DecodedBam {
        header_text,
        references,
        records,
    })
}

fn decode_record(bytes: &[u8], reference_count: usize) -> Result<DecodedRecord, DecodeError> {
    let mut cursor = BamCursor::new(bytes);
    let reference_id = cursor.i32("record reference id")?;
    let position = cursor.i32("record position")?;
    let read_name_length = usize::from(cursor.u8("read-name length")?);
    let mapping_quality = cursor.u8("mapping quality")?;
    let bin = cursor.u16("bin")?;
    let cigar_count = usize::from(cursor.u16("CIGAR count")?);
    let flag = cursor.u16("flag")?;
    let sequence_length = cursor.nonnegative_length("sequence length")?;
    let mate_reference_id = cursor.i32("mate reference id")?;
    let mate_position = cursor.i32("mate position")?;
    let template_length = cursor.i32("template length")?;

    validate_reference_position(reference_id, position, reference_count, &cursor)?;
    validate_reference_position(mate_reference_id, mate_position, reference_count, &cursor)?;

    let read_name_bytes = cursor.take(read_name_length, "read name")?;
    let Some((&0, read_name)) = read_name_bytes.split_last() else {
        return Err(cursor.fail("read-name terminator"));
    };
    if read_name.is_empty() || read_name.contains(&0) {
        return Err(cursor.fail("read name"));
    }

    if cigar_count > cursor.remaining() / 4 {
        return Err(cursor.fail("CIGAR count"));
    }
    let mut cigar = Vec::with_capacity(cigar_count);
    for _ in 0..cigar_count {
        let encoded = cursor.u32("CIGAR operation")?;
        let length = encoded >> 4;
        let operation = u8::try_from(encoded & 0xf).expect("low nibble fits u8");
        if length == 0 || operation > 9 {
            return Err(cursor.fail("CIGAR operation"));
        }
        cigar.push((length, operation));
    }

    let packed_length = sequence_length
        .checked_add(1)
        .ok_or_else(|| cursor.fail("packed sequence length"))?
        / 2;
    let packed = cursor.take(packed_length, "packed sequence")?;
    let sequence = decode_sequence(packed, sequence_length, &cursor)?;
    let raw_quality = cursor.take(sequence_length, "quality")?;
    let quality = decode_quality(raw_quality, &cursor)?;

    let mut aux = Vec::new();
    while cursor.remaining() != 0 {
        aux.push(decode_aux(&mut cursor)?);
    }
    Ok(DecodedRecord {
        reference_id,
        position,
        read_name: read_name.to_vec(),
        mapping_quality,
        bin,
        cigar,
        flag,
        mate_reference_id,
        mate_position,
        template_length,
        sequence,
        quality,
        aux,
    })
}

fn validate_reference_position(
    reference_id: i32,
    position: i32,
    reference_count: usize,
    cursor: &BamCursor<'_>,
) -> Result<(), DecodeError> {
    if reference_id == -1 && position == -1 {
        return Ok(());
    }
    let ordinal = usize::try_from(reference_id).map_err(|_| cursor.fail("reference id"))?;
    if ordinal >= reference_count || position < 0 {
        return Err(cursor.fail("reference id/position"));
    }
    Ok(())
}

fn decode_sequence(
    packed: &[u8],
    sequence_length: usize,
    cursor: &BamCursor<'_>,
) -> Result<Vec<u8>, DecodeError> {
    const BASES: &[u8; 16] = b"=ACMGRSVTWYHKDBN";
    let mut sequence = Vec::with_capacity(sequence_length);
    for index in 0..sequence_length {
        let byte = packed[index / 2];
        let code = if index % 2 == 0 {
            byte >> 4
        } else {
            byte & 0xf
        };
        let base = BASES
            .get(usize::from(code))
            .copied()
            .ok_or_else(|| cursor.fail("sequence base"))?;
        sequence.push(base);
    }
    Ok(sequence)
}

fn decode_quality(raw: &[u8], cursor: &BamCursor<'_>) -> Result<Option<Vec<u8>>, DecodeError> {
    if raw.iter().all(|value| *value == u8::MAX) {
        return Ok(None);
    }
    if raw.iter().any(|value| *value == u8::MAX || *value > 93) {
        return Err(cursor.fail("quality"));
    }
    Ok(Some(raw.iter().map(|value| value + 33).collect()))
}

fn decode_aux(cursor: &mut BamCursor<'_>) -> Result<DecodedAux, DecodeError> {
    let tag_bytes = cursor.take(2, "auxiliary tag")?;
    let tag = [tag_bytes[0], tag_bytes[1]];
    let physical_type = cursor.u8("auxiliary type")?;
    let value = match physical_type {
        b'c' => AuxValue::Integer(i64::from(i8::from_le_bytes([cursor.u8("i8 auxiliary")?]))),
        b'C' => AuxValue::Integer(i64::from(cursor.u8("u8 auxiliary")?)),
        b's' => AuxValue::Integer(i64::from(cursor.i16("i16 auxiliary")?)),
        b'S' => AuxValue::Integer(i64::from(cursor.u16("u16 auxiliary")?)),
        b'i' => AuxValue::Integer(i64::from(cursor.i32("i32 auxiliary")?)),
        b'I' => AuxValue::Integer(i64::from(cursor.u32("u32 auxiliary")?)),
        b'Z' => AuxValue::String(cursor.nul_terminated("string auxiliary")?),
        _ => return Err(cursor.fail("unsupported auxiliary type")),
    };
    Ok(DecodedAux { tag, value })
}

fn normalized(raw: &[u8]) -> NormalizedSequence {
    normalize_dna(raw).expect("test fixture is normalized DNA")
}

fn reference(catalog: &[(&[u8], &[u8])]) -> ReferenceIndex {
    ReferenceIndex::build(
        catalog
            .iter()
            .map(|(name, sequence)| ContigInput::new(name.to_vec(), normalized(sequence)))
            .collect(),
        ReferenceBuildLimits::MAX,
    )
    .expect("bounded reference builds")
}

fn single_record(
    reference: &ReferenceIndex,
    query_name: &[u8],
    raw: &[u8],
    quality: Option<&[u8]>,
    budget: u64,
) -> AlignmentRecord {
    single_record_with_mode(
        reference,
        query_name,
        raw,
        quality,
        budget,
        AlignmentAuxiliaryMode::Minimal,
    )
}

fn single_record_with_mode(
    reference: &ReferenceIndex,
    query_name: &[u8],
    raw: &[u8],
    quality: Option<&[u8]>,
    budget: u64,
    auxiliary_mode: AlignmentAuxiliaryMode,
) -> AlignmentRecord {
    let fixture = single_fixture(reference, raw, budget);
    build_single_alignment_record_with_auxiliary_mode(
        reference,
        query_name,
        AlignmentRead::new(&fixture.query, quality),
        fixture.alignment.as_ref(),
        fixture.mapping_quality,
        AlignmentRecordLimits::default(),
        auxiliary_mode,
    )
    .expect("single record builds")
}

fn unique_directory(label: &str) -> PathBuf {
    let nonce = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("clock")
        .as_nanos();
    std::env::temp_dir().join(format!(
        "bsbit-bam-fields-{label}-{}-{nonce}",
        std::process::id()
    ))
}

fn write_direct_bam_payload(
    directory: &Path,
    header: &SamHeader,
    records: &[&AlignmentRecord],
) -> Vec<u8> {
    fs::create_dir(directory).expect("test directory");
    let staging = directory.join("records.bam.tmp");
    let mut writer =
        BamStagingWriter::create_new(&staging, header, AlignmentRecordLimits::default())
            .expect("direct BAM writer opens");
    for record in records {
        writer
            .write_record_as_bam(record)
            .expect("direct BAM record writes");
    }
    let completed = writer.finish().expect("direct BAM finishes");
    let mut reader = DecodedReader::open(completed.path()).expect("direct BAM opens for decoding");
    assert_eq!(reader.compression(), Compression::Bgzf);
    let mut bytes = Vec::new();
    reader.read_to_end(&mut bytes).expect("direct BGZF decodes");
    reader.close().expect("decoded direct BAM closes");
    drop(completed);
    fs::remove_dir(directory).expect("direct test directory cleanup");
    bytes
}

fn write_borrowed_batch_payload(
    directory: &Path,
    header: &SamHeader,
    batch: &AlignmentRecordBatch,
) -> Vec<u8> {
    fs::create_dir(directory).expect("batch test directory");
    let staging = directory.join("records.bam.tmp");
    let mut writer =
        BamStagingWriter::create_new(&staging, header, AlignmentRecordLimits::default())
            .expect("batch BAM writer opens");
    for record in batch.records() {
        writer
            .write_borrowed_alignment_record(&record)
            .expect("borrowed batch BAM record writes");
    }
    let completed = writer.finish().expect("batch BAM finishes");
    let mut reader = DecodedReader::open(completed.path()).expect("batch BAM decodes");
    let mut bytes = Vec::new();
    reader.read_to_end(&mut bytes).expect("batch BGZF decodes");
    reader.close().expect("batch decoded stream closes");
    drop(completed);
    fs::remove_dir(directory).expect("batch test directory cleanup");
    bytes
}

#[test]
fn direct_bam_minimal_mode_emits_nm_and_xg_without_md() {
    let exact_reference = reference(&[(b"chr", b"GGACCTAA")]);
    let exact = single_record(&exact_reference, b"exact", b"ACCT", Some(b"ABCD"), 0);
    let header = build_sam_header(&exact_reference, AlignmentRecordLimits::default())
        .expect("header builds");
    let directory = unique_directory("single-exact-minimal");
    let payload = write_direct_bam_payload(&directory, &header, &[&exact]);
    let decoded = decode_bam(&payload).expect("minimal BAM independently decodes");
    assert_eq!(decoded.records.len(), 1);
    assert_eq!(decoded.records[0].aux.len(), 2);
    assert_eq!(
        decoded.records[0]
            .aux
            .iter()
            .map(|aux| aux.tag)
            .collect::<Vec<_>>(),
        [*b"NM", *b"XG"]
    );
    assert_eq!(
        aux_value(&decoded.records[0], *b"NM"),
        &AuxValue::Integer(i64::from(
            exact.mapping().expect("mapped fixture").literal_nm()
        ))
    );
    assert_eq!(
        aux_value(&decoded.records[0], *b"XG"),
        &AuxValue::String(
            exact
                .mapping()
                .expect("mapped fixture")
                .bismark_xg()
                .to_vec()
        )
    );
    assert!(decoded.records[0].aux.iter().all(|aux| aux.tag != *b"MD"));
}

#[test]
fn borrowed_unmapped_pair_batch_encodes_two_primary_bam_records() {
    let fixture_reference = reference(&[(b"chr", b"AACCGTGATCTAGGCTTACGGAAT")]);
    let first = normalized(b"ACGTN");
    let second = normalized(b"TGCAN");
    let mut batch = AlignmentRecordBatch::new();
    let mut composer = PairedRecordComposer::new();
    composer
        .push_unmapped_pair(
            b"unmapped-pair",
            BorrowedAlignmentRead::new(first.bases(), b"ABCDE"),
            BorrowedAlignmentRead::new(second.bases(), b"12345"),
            AlignmentRecordLimits::default(),
        )
        .expect("unmapped direct pair builds");
    composer
        .flush_into(&mut batch, AlignmentRecordLimits::default())
        .expect("unmapped pair flushes");
    let header = build_sam_header(&fixture_reference, AlignmentRecordLimits::default())
        .expect("header builds");
    let directory = unique_directory("batch-unmapped-pair");
    let payload = write_borrowed_batch_payload(&directory, &header, &batch);
    let decoded = decode_bam(&payload).expect("unmapped pair BAM independently decodes");

    assert_eq!(decoded.records.len(), 2);
    assert_eq!(decoded.records[0].read_name, b"unmapped-pair");
    assert_eq!(decoded.records[1].read_name, b"unmapped-pair");
    assert_eq!(decoded.records[0].flag, 77);
    assert_eq!(decoded.records[1].flag, 141);
    for record in &decoded.records {
        assert_eq!(record.reference_id, -1);
        assert_eq!(record.position, -1);
        assert_eq!(record.mapping_quality, 0);
        assert!(record.cigar.is_empty());
        assert_eq!(record.mate_reference_id, -1);
        assert_eq!(record.mate_position, -1);
        assert_eq!(record.template_length, 0);
        assert!(record.aux.is_empty());
    }
    assert_eq!(decoded.records[0].sequence, b"ACGTN");
    assert_eq!(decoded.records[0].quality, Some(b"ABCDE".to_vec()));
    assert_eq!(decoded.records[1].sequence, b"TGCAN");
    assert_eq!(decoded.records[1].quality, Some(b"12345".to_vec()));
}

fn assert_header_equal(header: &SamHeader, decoded: &DecodedBam) {
    assert_eq!(
        decoded.header_text,
        sam_header_bytes(header, AlignmentRecordLimits::default()).expect("header encodes")
    );
    assert_eq!(decoded.references.len(), header.references().len());
    for (actual, expected) in decoded.references.iter().zip(header.references()) {
        assert_eq!(actual.name, expected.name());
        assert_eq!(
            actual.length,
            i32::try_from(expected.length()).expect("validated length")
        );
    }
}

fn expected_bin(record: &AlignmentRecord) -> u16 {
    let Some(mapping) = record.mapping() else {
        return 4_680;
    };
    let interval = mapping.reference().interval();
    let start = interval.start();
    let end = interval.end() - 1;
    for (shift, offset) in [(14, 4_681_u64), (17, 585), (20, 73), (23, 9), (26, 1)] {
        if start >> shift == end >> shift {
            return u16::try_from(offset + (start >> shift)).expect("validated BAM bin");
        }
    }
    0
}

fn cigar_text(cigar: &[(u32, u8)]) -> Vec<u8> {
    const OPERATIONS: &[u8; 10] = b"MIDNSHP=XB";
    let mut text = Vec::new();
    for (length, operation) in cigar {
        text.extend_from_slice(length.to_string().as_bytes());
        text.push(OPERATIONS[usize::from(*operation)]);
    }
    text
}

fn aux_value(record: &DecodedRecord, tag: [u8; 2]) -> &AuxValue {
    &record
        .aux
        .iter()
        .find(|aux| aux.tag == tag)
        .expect("required auxiliary tag")
        .value
}

fn assert_record_equal(expected: &AlignmentRecord, actual: &DecodedRecord) {
    assert_eq!(actual.read_name, expected.query_name());
    assert_eq!(actual.flag, sam_flag(expected));
    assert_eq!(
        actual.mapping_quality,
        expected.mapping_quality().sam_value()
    );
    assert_eq!(actual.template_length, expected.template_length());
    assert_eq!(actual.sequence, expected.sequence());
    assert_eq!(actual.quality.as_deref(), expected.quality());
    assert_eq!(actual.bin, expected_bin(expected));

    if let Some(mapping) = expected.mapping() {
        assert_eq!(
            actual.reference_id,
            i32::try_from(mapping.reference().ordinal()).expect("validated ordinal")
        );
        assert_eq!(
            actual.position,
            i32::try_from(mapping.reference().position() - 1).expect("validated position")
        );
        assert_eq!(
            cigar_text(&actual.cigar),
            mapping.cigar().to_string().as_bytes()
        );
        assert_eq!(actual.aux.len(), 2);
        assert_eq!(
            actual.aux.iter().map(|aux| aux.tag).collect::<Vec<_>>(),
            [*b"NM", *b"XG"]
        );
        assert_eq!(
            aux_value(actual, *b"NM"),
            &AuxValue::Integer(i64::from(mapping.literal_nm()))
        );
        assert_eq!(
            aux_value(actual, *b"XG"),
            &AuxValue::String(mapping.bismark_xg().to_vec())
        );
    } else {
        assert_eq!(actual.reference_id, -1);
        assert_eq!(actual.position, -1);
        assert!(actual.cigar.is_empty());
        assert!(actual.aux.is_empty());
    }

    if let Some(mate) = expected.mate() {
        assert_eq!(
            actual.mate_reference_id,
            i32::try_from(mate.reference().ordinal()).expect("validated mate ordinal")
        );
        assert_eq!(
            actual.mate_position,
            i32::try_from(mate.reference().position() - 1).expect("validated mate position")
        );
    } else {
        assert_eq!(actual.mate_reference_id, -1);
        assert_eq!(actual.mate_position, -1);
    }
}

fn assert_round_trip(label: &str, reference: &ReferenceIndex, records: &[&AlignmentRecord]) {
    let header =
        build_sam_header(reference, AlignmentRecordLimits::default()).expect("header builds");
    let directory = unique_directory(label);
    let direct = write_direct_bam_payload(&directory, &header, records);
    let decoded = decode_bam(&direct).expect("direct BAM independently decodes");
    assert_header_equal(&header, &decoded);
    assert_eq!(decoded.records.len(), records.len());
    for (actual, expected) in decoded.records.iter().zip(records) {
        assert_record_equal(expected, actual);
    }
}

#[test]
fn independent_bam_decoder_matches_all_single_record_fields() {
    let exact_reference = reference(&[(b"chr", b"GGACCTAA")]);
    let exact = single_record(&exact_reference, b"exact", b"ACCT", Some(b"ABCD"), 0);
    assert_round_trip("single-exact", &exact_reference, &[&exact]);

    let reverse_reference = reference(&[(b"chr", b"TTAACGAA")]);
    let reverse = single_record(&reverse_reference, b"reverse", b"CGTT", Some(b"ABCD"), 0);
    assert_round_trip("single-reverse", &reverse_reference, &[&reverse]);

    let conversion_reference = reference(&[(b"chr", b"GGACCGAA")]);
    let conversion = single_record(&conversion_reference, b"conversion", b"ATTG", None, 0);
    assert_round_trip("single-conversion", &conversion_reference, &[&conversion]);

    let insertion_reference = reference(&[(b"chr", b"AGTC")]);
    let insertion = single_record(&insertion_reference, b"insertion", b"ACGTC", None, 1);
    assert_round_trip("single-insertion", &insertion_reference, &[&insertion]);

    let deletion_reference = reference(&[(b"chr", b"ACGAC")]);
    let deletion = single_record(&deletion_reference, b"deletion", b"ACAC", None, 1);
    assert_round_trip("single-deletion", &deletion_reference, &[&deletion]);

    let unmapped_reference = reference(&[(b"chr", b"AAAA")]);
    let unmapped = single_record(&unmapped_reference, b"unmapped", b"GGGG", Some(b"!!!!"), 0);
    assert_round_trip("single-unmapped", &unmapped_reference, &[&unmapped]);
}
