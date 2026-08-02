//! Public BAM/HTSlib contracts over format-level values only.

use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use bsbit_core::bisulfite::{AlignmentOrientation, BisulfiteStrand, CytosineStrand};
use bsbit_core::cigar::CoreCigar;
use bsbit_core::coordinate::{ReferenceInterval, ReferenceLength};
use bsbit_hts::{
    AlignmentAuxiliaryMode, AlignmentCigarOp, AlignmentCigarRun, AlignmentRecord,
    AlignmentRecordBatch, AlignmentRecordLimits, AlignmentRecordResource, BamStagingWriter,
    BorrowedAlignmentRecord, BsbitAlignmentMode, BsbitProgramProvenance, Compression,
    DecodedReader, HtsErrorKind, HtsOperation, IndexedBamReader, MappedAlignmentRecord,
    RecordMappingQuality, RecordReference, RecordSegment, SamHeader, SamHeaderReference,
    SamSortOrder, build_bam_index_create_new,
};

fn unique_directory(label: &str) -> PathBuf {
    let nonce = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("clock")
        .as_nanos();
    std::env::temp_dir().join(format!(
        "bsbit-hts-bam-{label}-{}-{nonce}",
        std::process::id()
    ))
}

fn header() -> SamHeader {
    SamHeader::new(
        vec![SamHeaderReference::new(0, b"chr1", 8).expect("dictionary entry")],
        AlignmentRecordLimits::default(),
    )
    .expect("header")
    .with_sort_order(SamSortOrder::Coordinate)
}

fn provenance_writer(path: &Path, provenance: BsbitProgramProvenance) -> BamStagingWriter {
    let provenance_header = header()
        .with_bsbit_provenance(provenance, AlignmentRecordLimits::default())
        .expect("provenance header");
    BamStagingWriter::create_new(path, &provenance_header, AlignmentRecordLimits::default())
        .expect("BAM staging")
}

fn assert_bsbit_provenance(reader: &IndexedBamReader, expected: BsbitProgramProvenance) {
    assert_eq!(
        reader
            .header()
            .bsbit_program_provenance()
            .expect("valid BAM provenance"),
        Some(expected)
    );
}

fn shared_record_fixture() -> AlignmentRecord {
    let limits = AlignmentRecordLimits::default();
    let interval = ReferenceInterval::new(0, 4, ReferenceLength::new(8)).expect("interval");
    let reference = RecordReference::new(0, b"chr1", 8, interval).expect("reference");
    let mapping = MappedAlignmentRecord::new(
        reference,
        AlignmentOrientation::Forward,
        BisulfiteStrand::OT,
        CytosineStrand::Top,
        CoreCigar::all_matches(4),
        4,
        0,
        None,
        None,
        limits,
    )
    .expect("mapping");
    AlignmentRecord::new(
        b"read1",
        RecordSegment::Unpaired,
        false,
        RecordMappingQuality::Calibrated(42),
        Some(mapping),
        None,
        0,
        b"ACGT",
        Some(b"IIII"),
        limits,
    )
    .expect("record")
}

fn remove_if_present(path: &Path) {
    match fs::remove_file(path) {
        Ok(()) => {}
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => {}
        Err(error) => panic!("remove {}: {error}", path.display()),
    }
}

#[test]
fn direct_fields_round_trip_through_public_bam_and_index_contracts() {
    let directory = unique_directory("round-trip");
    fs::create_dir(&directory).expect("directory");
    let staging = directory.join("result.tmp");
    let target = directory.join("result.bam");
    let index = directory.join("result.bam.bai");
    let cigar = [AlignmentCigarRun::new(AlignmentCigarOp::Match, 4).expect("CIGAR")];
    let direct = BorrowedAlignmentRecord::new(
        b"read1",
        0,
        Some(0),
        1,
        42,
        &cigar,
        None,
        0,
        0,
        b"ACGT",
        Some(b"IIII"),
        0,
        AlignmentAuxiliaryMode::Bismark,
        Some(b"4"),
        BisulfiteStrand::OT,
        Some(b"...."),
        AlignmentRecordLimits::default(),
    )
    .expect("direct record");
    let mut batch = AlignmentRecordBatch::new();
    batch.push(&direct).expect("batch retention");

    let expected_provenance = BsbitProgramProvenance::new(
        [0x5a; 32],
        BsbitAlignmentMode::CallerCompatibleDirectionalPaired,
    );
    let mut writer = provenance_writer(&staging, expected_provenance);
    writer
        .write_borrowed_alignment_record(&batch.records().next().expect("retained record"))
        .expect("direct BAM fields");
    let completed = writer.finish().expect("BAM finalization");
    assert_eq!(completed.records_written(), 1);
    let publication = completed
        .publish_create_new(&target)
        .expect("create-only publication");
    assert_eq!(publication.records_written(), 1);
    assert_eq!(publication.target_path(), target);

    let mut decoded = DecodedReader::open(&target).expect("decoded BAM source");
    assert_eq!(decoded.compression(), Compression::Bgzf);
    let mut magic = [0_u8; 4];
    decoded.read_exact(&mut magic).expect("BAM magic");
    assert_eq!(&magic, b"BAM\x01");
    decoded.close().expect("decoded source closes");

    build_bam_index_create_new(&target, &index, 1).expect("BAI creation");
    let original_index = fs::read(&index).expect("published BAI bytes");
    build_bam_index_create_new(&target, &index, 1)
        .expect_err("existing BAI path must win create-only reservation");
    assert_eq!(
        fs::read(&index).expect("preserved BAI bytes"),
        original_index
    );
    let mut reader = IndexedBamReader::open(&target).expect("indexed BAM");
    assert!(reader.header().is_coordinate_sorted());
    assert_eq!(reader.header().references()[0].name(), b"chr1");
    assert_bsbit_provenance(&reader, expected_provenance);
    reader.query(0, 0, 8).expect("region query");
    let record = reader
        .next_record()
        .expect("record read")
        .expect("record present");
    assert_eq!(record.query_name(), b"read1");
    assert_eq!(record.reference_id(), 0);
    assert_eq!(record.position(), 0);
    assert_eq!(record.mapping_quality(), 42);
    assert_eq!(record.cigar(), &[4_u32 << 4]);
    let mut sequence = Vec::new();
    record
        .decode_sequence_into(&mut sequence)
        .expect("sequence decode");
    assert_eq!(sequence, b"ACGT");
    assert_eq!(
        record.string_auxiliary(*b"MD").expect("MD"),
        Some(b"4".as_slice())
    );
    assert_eq!(
        record.string_auxiliary(*b"XM").expect("XM"),
        Some(b"....".as_slice())
    );
    assert_eq!(
        record.string_auxiliary(*b"XR").expect("XR"),
        Some(b"CT".as_slice())
    );
    assert_eq!(
        record.string_auxiliary(*b"XG").expect("XG"),
        Some(b"CT".as_slice())
    );
    assert!(reader.next_record().expect("query EOF").is_none());
    reader.close().expect("indexed reader closes");

    remove_if_present(&index);
    remove_if_present(&target);
    remove_if_present(&staging);
    fs::remove_dir(directory).expect("directory cleanup");
}

#[test]
fn bam_writer_is_terminal_after_encode_failure_and_never_replaces_paths() {
    let directory = unique_directory("faults");
    fs::create_dir(&directory).expect("directory");
    let staging = directory.join("result.tmp");
    fs::write(&staging, b"caller-owned").expect("existing staging");
    let collision =
        BamStagingWriter::create_new(&staging, &header(), AlignmentRecordLimits::default())
            .err()
            .expect("staging collision");
    assert_eq!(collision.operation(), HtsOperation::CreateStaging);
    assert_eq!(fs::read(&staging).expect("sentinel"), b"caller-owned");
    fs::remove_file(&staging).expect("sentinel cleanup");

    let mut writer =
        BamStagingWriter::create_new(&staging, &header(), AlignmentRecordLimits::default())
            .expect("writer");
    let tiny = AlignmentRecordLimits::new(254, 100, 100, 100, 100, 100, 100, 1, 10, 100, 1_000);
    let first = writer
        .write_record(&shared_record_fixture(), tiny)
        .expect_err("SAM-size encode cap");
    assert_eq!(first.kind(), HtsErrorKind::Encode);
    let terminal = writer
        .write_record(&shared_record_fixture(), AlignmentRecordLimits::default())
        .expect_err("terminal writer");
    assert_eq!(terminal.kind(), HtsErrorKind::Terminal);
    drop(writer);
    assert!(!staging.exists());

    let target = directory.join("result.bam");
    let mut writer =
        BamStagingWriter::create_new(&staging, &header(), AlignmentRecordLimits::default())
            .expect("second writer");
    writer
        .write_record_as_bam(&shared_record_fixture())
        .expect("valid record");
    let completed = writer.finish().expect("completed BAM");
    fs::write(&target, b"caller target").expect("target collision");
    let collision = completed
        .publish_create_new(&target)
        .expect_err("target collision");
    assert_eq!(
        collision.operation(),
        HtsOperation::ValidatePublicationPaths
    );
    assert_eq!(
        collision.kind(),
        HtsErrorKind::Io(std::io::ErrorKind::AlreadyExists)
    );
    assert_eq!(
        fs::read(&target).expect("target sentinel"),
        b"caller target"
    );
    assert!(!staging.exists());

    fs::remove_file(target).expect("target cleanup");
    fs::remove_dir(directory).expect("directory cleanup");
}

#[test]
fn explicit_bam_compression_preserves_uncompressed_bam_bytes() {
    let directory = unique_directory("compression-level");
    fs::create_dir(&directory).expect("directory");
    let default_path = directory.join("default.bam.tmp");
    let level_one_path = directory.join("level-one.bam.tmp");
    let limits = AlignmentRecordLimits::default();
    let header = header();
    let record = shared_record_fixture();

    let mut default_writer =
        BamStagingWriter::create_new_with_threads(&default_path, &header, limits, 2)
            .expect("default writer");
    default_writer
        .write_record(&record, limits)
        .expect("default record");
    let default_completed = default_writer.finish().expect("default finish");

    let mut level_one_writer = BamStagingWriter::create_new_with_threads_and_compression_level(
        &level_one_path,
        &header,
        limits,
        2,
        1,
    )
    .expect("level-one writer");
    level_one_writer
        .write_record(&record, limits)
        .expect("level-one record");
    let level_one_completed = level_one_writer.finish().expect("level-one finish");

    let mut default_bytes = Vec::new();
    let mut default_reader = DecodedReader::open(&default_path).expect("default opens");
    default_reader
        .read_to_end(&mut default_bytes)
        .expect("default decodes");
    default_reader.close().expect("default closes");
    let mut level_one_bytes = Vec::new();
    let mut level_one_reader = DecodedReader::open(&level_one_path).expect("level one opens");
    level_one_reader
        .read_to_end(&mut level_one_bytes)
        .expect("level one decodes");
    level_one_reader.close().expect("level one closes");
    assert_eq!(level_one_bytes, default_bytes);

    drop(default_completed);
    drop(level_one_completed);
    assert!(!default_path.exists());
    assert!(!level_one_path.exists());
    fs::remove_dir(directory).expect("directory cleanup");
}

#[test]
fn staging_is_create_only_and_unfinished_drop_cleans_up() {
    let directory = unique_directory("unfinished-drop");
    fs::create_dir(&directory).expect("directory");
    let staging = directory.join("output.tmp");
    let writer =
        BamStagingWriter::create_new(&staging, &header(), AlignmentRecordLimits::default())
            .expect("fresh staging");
    assert!(staging.exists());
    drop(writer);
    assert!(!staging.exists());
    fs::remove_dir(directory).expect("directory cleanup");
}

#[test]
fn sibling_staging_is_private_unique_and_publishes_to_the_requested_target() {
    let directory = unique_directory("sibling-staging");
    fs::create_dir(&directory).expect("directory");
    let target = directory.join("output.bam");
    let writer =
        BamStagingWriter::create_sibling(&target, &header(), AlignmentRecordLimits::default())
            .expect("sibling staging");
    let staging = writer.path().to_path_buf();
    assert_ne!(staging, target);
    assert_eq!(staging.parent(), target.parent());
    assert!(staging.exists());

    let publication = writer
        .finish()
        .expect("finish")
        .publish_create_new(&target)
        .expect("publish");
    assert_eq!(publication.target_path(), target);
    assert_eq!(publication.staging_path(), staging);
    assert!(target.exists());
    assert!(!staging.exists());

    fs::remove_file(target).expect("target cleanup");
    fs::remove_dir(directory).expect("directory cleanup");
}

#[test]
fn active_staging_replacement_is_never_truncated_or_removed() {
    let directory = unique_directory("active-replacement");
    fs::create_dir(&directory).expect("directory");
    let staging = directory.join("output.tmp");
    let mut writer =
        BamStagingWriter::create_new(&staging, &header(), AlignmentRecordLimits::default())
            .expect("writer");

    fs::remove_file(&staging).expect("unlink owned staging");
    fs::write(&staging, b"replacement bytes").expect("replacement");
    writer
        .write_record(&shared_record_fixture(), AlignmentRecordLimits::default())
        .expect("native writer retains its descriptor");
    let error = writer.finish().err().expect("identity change fails finish");
    assert_eq!(error.kind(), HtsErrorKind::StagingIdentityChanged);
    assert_eq!(
        fs::read(&staging).expect("replacement survives"),
        b"replacement bytes"
    );

    fs::remove_file(staging).expect("replacement cleanup");
    fs::remove_dir(directory).expect("directory cleanup");
}

#[test]
fn completed_staging_replacement_is_never_removed() {
    let directory = unique_directory("completed-replacement");
    fs::create_dir(&directory).expect("directory");
    let staging = directory.join("output.tmp");
    let writer =
        BamStagingWriter::create_new(&staging, &header(), AlignmentRecordLimits::default())
            .expect("writer");
    let completed = writer.finish().expect("finish");

    fs::remove_file(&staging).expect("unlink completed staging");
    fs::write(&staging, b"replacement bytes").expect("replacement");
    let error = completed
        .remove()
        .expect_err("identity change blocks cleanup");
    assert_eq!(error.kind(), HtsErrorKind::StagingIdentityChanged);
    assert_eq!(
        fs::read(&staging).expect("replacement survives"),
        b"replacement bytes"
    );

    fs::remove_file(staging).expect("replacement cleanup");
    fs::remove_dir(directory).expect("directory cleanup");
}

#[test]
fn concurrent_bam_publishers_have_exactly_one_winner() {
    const WORKERS: usize = 8;

    let directory = unique_directory("publish-race");
    fs::create_dir(&directory).expect("directory");
    let target = directory.join("output.bam");
    let header = header();
    let mut completed = Vec::new();
    let mut staging_paths = Vec::new();
    for worker in 0..WORKERS {
        let staging = directory.join(format!("output-{worker}.tmp"));
        let writer =
            BamStagingWriter::create_new(&staging, &header, AlignmentRecordLimits::default())
                .expect("writer");
        completed.push(writer.finish().expect("finish"));
        staging_paths.push(staging);
    }
    let handles: Vec<_> = completed
        .into_iter()
        .map(|bam| {
            let target = target.clone();
            std::thread::spawn(move || bam.publish_create_new(target))
        })
        .collect();
    let mut successes = 0;
    let mut occupied = 0;
    let mut unsupported = 0;
    for handle in handles {
        match handle.join().expect("publisher thread") {
            Ok(publication) => {
                successes += 1;
                assert_eq!(publication.cleanup_warning(), None);
            }
            Err(error) if error.kind() == HtsErrorKind::Io(std::io::ErrorKind::AlreadyExists) => {
                occupied += 1;
            }
            Err(error) if error.kind() == HtsErrorKind::Io(std::io::ErrorKind::Unsupported) => {
                unsupported += 1;
            }
            Err(error) => panic!("unexpected publication error: {error}"),
        }
    }
    if unsupported == WORKERS {
        assert_eq!(successes, 0);
        assert_eq!(occupied, 0);
    } else {
        assert_eq!(unsupported, 0);
        assert_eq!(successes, 1);
        assert_eq!(occupied, WORKERS - 1);
        fs::remove_file(&target).expect("target cleanup");
    }
    assert!(staging_paths.iter().all(|path| !path.exists()));
    fs::remove_dir(directory).expect("directory cleanup");
}

#[test]
fn direct_record_constructor_enforces_bam_text_and_optional_caps() {
    let cigar = [AlignmentCigarRun::new(AlignmentCigarOp::Match, 10).expect("CIGAR")];
    let limits = AlignmentRecordLimits::new(254, 100, 100, 10, 2, 100, 100, 100, 10, 100, 1_000);
    let error = BorrowedAlignmentRecord::new(
        b"read",
        0,
        Some(0),
        1,
        1,
        &cigar,
        None,
        0,
        0,
        b"ACGTACGTAA",
        Some(b"IIIIIIIIII"),
        0,
        AlignmentAuxiliaryMode::Bismark,
        Some(b"10"),
        BisulfiteStrand::OT,
        Some(b".........."),
        limits,
    )
    .err()
    .expect("CIGAR text cap");
    assert!(matches!(
        error,
        bsbit_hts::AlignmentRecordError::LimitExceeded {
            resource: AlignmentRecordResource::CigarTextBytes,
            ..
        }
    ));
}

#[allow(dead_code)]
mod mutation_oracle {
    use super::*;

    #[derive(Debug, Eq, PartialEq)]
    struct DecodeError {
        offset: usize,
        context: &'static str,
    }

    #[derive(Debug, Eq, PartialEq)]
    struct DecodedBam {
        record_count: usize,
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

        fn nul_terminated(&mut self, context: &'static str) -> Result<(), DecodeError> {
            let Some(length) = self.bytes[self.offset..].iter().position(|byte| *byte == 0) else {
                return Err(self.fail(context));
            };
            self.take(length, context)?;
            self.take(1, context)?;
            Ok(())
        }
    }

    fn decode_bam(bytes: &[u8]) -> Result<DecodedBam, DecodeError> {
        let mut cursor = BamCursor::new(bytes);
        if cursor.take(4, "BAM magic")? != b"BAM\x01" {
            return Err(cursor.fail("BAM magic"));
        }
        let header_length = cursor.nonnegative_length("header length")?;
        cursor.take(header_length, "header text")?;
        let reference_count = cursor.nonnegative_length("reference count")?;
        if reference_count > cursor.remaining() / 8 {
            return Err(cursor.fail("reference count"));
        }
        for _ in 0..reference_count {
            let name_length = cursor.nonnegative_length("reference name length")?;
            let name_bytes = cursor.take(name_length, "reference name")?;
            let Some((&0, name)) = name_bytes.split_last() else {
                return Err(cursor.fail("reference name terminator"));
            };
            if name.is_empty() || name.contains(&0) {
                return Err(cursor.fail("reference name"));
            }
            if cursor.i32("reference length")? < 0 {
                return Err(cursor.fail("reference length"));
            }
        }

        let mut record_count = 0_usize;
        while cursor.remaining() != 0 {
            let block_length = cursor.nonnegative_length("record block length")?;
            let block = cursor.take(block_length, "record block")?;
            decode_record(block, reference_count)?;
            record_count += 1;
        }
        Ok(DecodedBam { record_count })
    }

    fn decode_record(bytes: &[u8], reference_count: usize) -> Result<(), DecodeError> {
        let mut cursor = BamCursor::new(bytes);
        let reference_id = cursor.i32("record reference id")?;
        let position = cursor.i32("record position")?;
        let read_name_length = usize::from(cursor.u8("read-name length")?);
        cursor.u8("mapping quality")?;
        cursor.u16("bin")?;
        let cigar_count = usize::from(cursor.u16("CIGAR count")?);
        cursor.u16("flag")?;
        let sequence_length = cursor.nonnegative_length("sequence length")?;
        let mate_reference_id = cursor.i32("mate reference id")?;
        let mate_position = cursor.i32("mate position")?;
        cursor.i32("template length")?;

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
        for _ in 0..cigar_count {
            let encoded = cursor.u32("CIGAR operation")?;
            let length = encoded >> 4;
            let operation = u8::try_from(encoded & 0xf).expect("low nibble fits u8");
            if length == 0 || operation > 9 {
                return Err(cursor.fail("CIGAR operation"));
            }
        }

        let packed_length = sequence_length
            .checked_add(1)
            .ok_or_else(|| cursor.fail("packed sequence length"))?
            / 2;
        let packed = cursor.take(packed_length, "packed sequence")?;
        decode_sequence(packed, sequence_length, &cursor)?;
        let raw_quality = cursor.take(sequence_length, "quality")?;
        decode_quality(raw_quality, &cursor)?;
        while cursor.remaining() != 0 {
            decode_aux(&mut cursor)?;
        }
        Ok(())
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
    ) -> Result<(), DecodeError> {
        const BASES: &[u8; 16] = b"=ACMGRSVTWYHKDBN";
        for index in 0..sequence_length {
            let byte = packed[index / 2];
            let code = if index % 2 == 0 {
                byte >> 4
            } else {
                byte & 0xf
            };
            BASES
                .get(usize::from(code))
                .ok_or_else(|| cursor.fail("sequence base"))?;
        }
        Ok(())
    }

    fn decode_quality(raw: &[u8], cursor: &BamCursor<'_>) -> Result<(), DecodeError> {
        if raw.iter().all(|value| *value == u8::MAX) {
            return Ok(());
        }
        if raw.iter().any(|value| *value == u8::MAX || *value > 93) {
            return Err(cursor.fail("quality"));
        }
        Ok(())
    }

    fn decode_aux(cursor: &mut BamCursor<'_>) -> Result<(), DecodeError> {
        cursor.take(2, "auxiliary tag")?;
        match cursor.u8("auxiliary type")? {
            b'c' => {
                cursor.u8("i8 auxiliary")?;
            }
            b'C' => {
                cursor.u8("u8 auxiliary")?;
            }
            b's' => {
                cursor.i16("i16 auxiliary")?;
            }
            b'S' => {
                cursor.u16("u16 auxiliary")?;
            }
            b'i' => {
                cursor.i32("i32 auxiliary")?;
            }
            b'I' => {
                cursor.u32("u32 auxiliary")?;
            }
            b'Z' => cursor.nul_terminated("string auxiliary")?,
            _ => return Err(cursor.fail("unsupported auxiliary type")),
        }
        Ok(())
    }

    struct FirstRecordOffsets {
        reference_count: usize,
        reference_count_offset: usize,
        block: usize,
        core: usize,
        read_name: usize,
        read_name_length: usize,
        cigar: usize,
        auxiliary: usize,
    }

    fn first_record_offsets(bytes: &[u8]) -> FirstRecordOffsets {
        let header_length = usize::try_from(i32::from_le_bytes(
            bytes[4..8].try_into().expect("header length bytes"),
        ))
        .expect("nonnegative header length");
        let reference_count_offset = 8 + header_length;
        let reference_count = usize::try_from(i32::from_le_bytes(
            bytes[reference_count_offset..reference_count_offset + 4]
                .try_into()
                .expect("reference count bytes"),
        ))
        .expect("nonnegative reference count");
        let mut block = reference_count_offset + 4;
        for _ in 0..reference_count {
            let name_length = usize::try_from(i32::from_le_bytes(
                bytes[block..block + 4]
                    .try_into()
                    .expect("name length bytes"),
            ))
            .expect("nonnegative name length");
            block += 4 + name_length + 4;
        }
        let core = block + 4;
        let read_name_length = usize::from(bytes[core + 8]);
        let cigar_count = usize::from(u16::from_le_bytes(
            bytes[core + 12..core + 14]
                .try_into()
                .expect("CIGAR count bytes"),
        ));
        let sequence_length = usize::try_from(i32::from_le_bytes(
            bytes[core + 16..core + 20]
                .try_into()
                .expect("sequence length bytes"),
        ))
        .expect("nonnegative sequence length");
        let read_name = core + 32;
        let cigar = read_name + read_name_length;
        let packed_sequence = cigar + cigar_count * 4;
        let auxiliary = packed_sequence + sequence_length.div_ceil(2) + sequence_length;
        FirstRecordOffsets {
            reference_count,
            reference_count_offset,
            block,
            core,
            read_name,
            read_name_length,
            cigar,
            auxiliary,
        }
    }

    fn targeted_structural_mutations(bytes: &[u8]) -> Vec<Vec<u8>> {
        let offsets = first_record_offsets(bytes);
        let mut mutations = Vec::new();
        let mut bad_magic = bytes.to_vec();
        bad_magic[0] ^= 0xff;
        mutations.push(bad_magic);
        let mut negative_header = bytes.to_vec();
        negative_header[4..8].copy_from_slice(&(-1_i32).to_le_bytes());
        mutations.push(negative_header);
        let mut negative_references = bytes.to_vec();
        negative_references[offsets.reference_count_offset..offsets.reference_count_offset + 4]
            .copy_from_slice(&(-1_i32).to_le_bytes());
        mutations.push(negative_references);
        let mut oversized_block = bytes.to_vec();
        oversized_block[offsets.block..offsets.block + 4].copy_from_slice(&i32::MAX.to_le_bytes());
        mutations.push(oversized_block);
        let mut bad_reference = bytes.to_vec();
        bad_reference[offsets.core..offsets.core + 4].copy_from_slice(
            &i32::try_from(offsets.reference_count)
                .expect("tiny reference count")
                .to_le_bytes(),
        );
        mutations.push(bad_reference);
        let mut unterminated_name = bytes.to_vec();
        unterminated_name[offsets.read_name + offsets.read_name_length - 1] = b'X';
        mutations.push(unterminated_name);
        let mut illegal_cigar = bytes.to_vec();
        illegal_cigar[offsets.cigar] = (illegal_cigar[offsets.cigar] & 0xf0) | 0x0f;
        mutations.push(illegal_cigar);
        let mut invalid_auxiliary = bytes.to_vec();
        invalid_auxiliary[offsets.auxiliary + 2] = b'Q';
        mutations.push(invalid_auxiliary);
        mutations
    }

    fn write_bam_payload(directory: &Path) -> Vec<u8> {
        fs::create_dir(directory).expect("test directory");
        let staging = directory.join("records.bam.tmp");
        let mut writer =
            BamStagingWriter::create_new(&staging, &header(), AlignmentRecordLimits::default())
                .expect("BAM writer opens");
        writer
            .write_record(&shared_record_fixture(), AlignmentRecordLimits::default())
            .expect("BAM record writes");
        let completed = writer.finish().expect("BAM finishes");
        let mut reader = DecodedReader::open(completed.path()).expect("BAM opens for decoding");
        assert_eq!(reader.compression(), Compression::Bgzf);
        let mut bytes = Vec::new();
        reader.read_to_end(&mut bytes).expect("BGZF decodes");
        reader.close().expect("decoded BAM closes");
        drop(completed);
        fs::remove_dir(directory).expect("test directory cleanup");
        bytes
    }

    pub(super) fn run_serialized_bam_mutation_contract() {
        let bytes = write_bam_payload(&unique_directory("mutations"));

        let baseline = decode_bam(&bytes).expect("baseline BAM is valid");
        assert_eq!(baseline.record_count, 1);
        for mutation in targeted_structural_mutations(&bytes) {
            assert!(decode_bam(&mutation).is_err());
        }

        let mut accepted_bit_flips = 0_usize;
        let mut rejected_bit_flips = 0_usize;
        for index in 0..bytes.len() {
            let mut mutation = bytes.clone();
            mutation[index] ^= 1_u8 << (index % 8);
            let first = decode_bam(&mutation);
            let second = decode_bam(&mutation);
            assert_eq!(first, second);
            if first.is_ok() {
                accepted_bit_flips += 1;
            } else {
                rejected_bit_flips += 1;
            }
        }
        assert!(
            accepted_bit_flips > 0,
            "semantic byte mutations can remain valid"
        );
        assert!(
            rejected_bit_flips > 0,
            "structural byte mutations fail closed"
        );

        let mut accepted_truncations = 0_usize;
        let mut rejected_truncations = 0_usize;
        for length in 0..bytes.len() {
            let first = decode_bam(&bytes[..length]);
            let second = decode_bam(&bytes[..length]);
            assert_eq!(first, second);
            if first.is_ok() {
                accepted_truncations += 1;
            } else {
                rejected_truncations += 1;
            }
        }
        assert!(
            accepted_truncations > 0,
            "a complete header is a valid zero-record BAM"
        );
        assert!(rejected_truncations > 0, "partial structures fail closed");
    }
}

#[test]
fn serialized_bam_mutations_are_bounded_deterministic_and_structurally_checked() {
    mutation_oracle::run_serialized_bam_mutation_contract();
}
