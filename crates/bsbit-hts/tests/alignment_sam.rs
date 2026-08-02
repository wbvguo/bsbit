//! Public contracts shared by the alignment model and SAM/BAM encoders.

use std::fs;
use std::path::PathBuf;
use std::time::{SystemTime, UNIX_EPOCH};

use bsbit_core::bisulfite::{AlignmentOrientation, BisulfiteStrand, CytosineStrand};
use bsbit_core::cigar::CoreCigar;
use bsbit_core::coordinate::{ReferenceInterval, ReferenceLength};
use bsbit_hts::{
    AlignmentAuxiliaryMode, AlignmentCigarOp, AlignmentCigarRun, AlignmentRecord,
    AlignmentRecordError, AlignmentRecordLimits, AlignmentRecordResource, BamStagingWriter,
    BorrowedAlignmentRecord, MappedAlignmentRecord, RecordMappingQuality, RecordReference,
    RecordSegment, SamFileWriter, SamHeader, SamHeaderReference, sam_borrowed_record_bytes,
    sam_header_bytes, sam_record_bytes,
};

fn unique_path(label: &str) -> PathBuf {
    let nonce = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("clock")
        .as_nanos();
    std::env::temp_dir().join(format!("bsbit-hts-{label}-{}-{nonce}", std::process::id()))
}

fn fixture() -> (SamHeader, AlignmentRecord) {
    let limits = AlignmentRecordLimits::default();
    let header = SamHeader::new(
        vec![SamHeaderReference::new(0, b"chr1", 8).expect("dictionary entry")],
        limits,
    )
    .expect("header");
    let interval = ReferenceInterval::new(0, 4, ReferenceLength::new(8)).expect("interval");
    let reference = RecordReference::new(0, b"chr1", 8, interval).expect("coordinate");
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
    let record = AlignmentRecord::new(
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
    .expect("record");
    (header, record)
}

#[test]
fn sam_and_bam_share_one_validated_alignment_model() {
    let (header, record) = fixture();
    assert_eq!(
        sam_record_bytes(&record, AlignmentRecordLimits::default()).expect("SAM record"),
        b"read1\t0\tchr1\t1\t42\t4M\t*\t0\t0\tACGT\tIIII\tNM:i:0\tXG:Z:CT\n"
    );
    assert!(
        sam_header_bytes(&header, AlignmentRecordLimits::default())
            .expect("SAM header")
            .starts_with(b"@HD\tVN:1.6")
    );

    let staging = unique_path("bam-stage");
    let target = staging.with_extension("bam");
    let mut writer =
        BamStagingWriter::create_new(&staging, &header, AlignmentRecordLimits::default())
            .expect("BAM staging");
    writer
        .write_record_as_bam(&record)
        .expect("direct BAM record");
    match writer
        .finish()
        .expect("complete BAM")
        .publish_create_new(&target)
    {
        Ok(publication) => {
            assert_eq!(publication.records_written(), 1);
            assert!(fs::metadata(&target).expect("BAM target").len() > 0);
            fs::remove_file(target).expect("target cleanup");
        }
        Err(error)
            if error.kind() == bsbit_hts::HtsErrorKind::Io(std::io::ErrorKind::Unsupported) => {}
        Err(error) => panic!("publish BAM: {error}"),
    }
}

#[test]
fn owned_sam_uses_the_records_bismark_auxiliary_mode() {
    let limits = AlignmentRecordLimits::default();
    let interval = ReferenceInterval::new(0, 4, ReferenceLength::new(8)).expect("interval");
    let reference = RecordReference::new(0, b"chr1", 8, interval).expect("coordinate");
    let mapping = MappedAlignmentRecord::new(
        reference,
        AlignmentOrientation::Forward,
        BisulfiteStrand::OT,
        CytosineStrand::Top,
        CoreCigar::all_matches(4),
        4,
        0,
        Some(b"4"),
        Some(b"...."),
        limits,
    )
    .expect("Bismark mapping");
    assert_eq!(mapping.auxiliary_mode(), AlignmentAuxiliaryMode::Bismark);
    let record = AlignmentRecord::new(
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
    .expect("record");
    assert_eq!(
        sam_record_bytes(&record, limits).expect("SAM record"),
        b"read1\t0\tchr1\t1\t42\t4M\t*\t0\t0\tACGT\tIIII\tNM:i:0\tXG:Z:CT\tMD:Z:4\tXM:Z:....\tXR:Z:CT\n"
    );
}

#[test]
fn compact_borrowed_record_is_shared_by_sam_and_bam_without_owned_conversion() {
    let (header, _) = fixture();
    let cigar = [AlignmentCigarRun::new(AlignmentCigarOp::Match, 4).expect("CIGAR")];
    let record = BorrowedAlignmentRecord::new(
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
    .expect("compact record");
    assert_eq!(
        sam_borrowed_record_bytes(&record, &header, AlignmentRecordLimits::default())
            .expect("SAM compact record"),
        b"read1\t0\tchr1\t1\t42\t4M\t*\t0\t0\tACGT\tIIII\tNM:i:0\tXG:Z:CT\tMD:Z:4\tXM:Z:....\tXR:Z:CT\n"
    );

    let staging = unique_path("borrowed-sam-stage");
    let target = staging.with_extension("sam");
    let mut writer =
        SamFileWriter::create_new(&target, &staging, &header, AlignmentRecordLimits::default())
            .expect("SAM staging");
    writer
        .write_borrowed_record(&record, AlignmentRecordLimits::default())
        .expect("SAM compact record writes");
    match writer.finish() {
        Ok(publication) => {
            assert_eq!(publication.records_written(), 1);
            let bytes = fs::read(&target).expect("SAM target");
            assert!(bytes.ends_with(b"\tXR:Z:CT\n"));
            fs::remove_file(target).expect("target cleanup");
        }
        Err(error) if error.kind() == Some(std::io::ErrorKind::Unsupported) => {}
        Err(error) => panic!("publish SAM: {error}"),
    }
}

#[test]
fn compact_sam_rnext_uses_reference_ordinal_not_a_duplicate_name() {
    let limits = AlignmentRecordLimits::default();
    let header = SamHeader::new(
        vec![
            SamHeaderReference::new(0, b"dup", 8).expect("first dictionary entry"),
            SamHeaderReference::new(1, b"dup", 8).expect("second dictionary entry"),
        ],
        limits,
    )
    .expect("header");
    let cigar = [AlignmentCigarRun::new(AlignmentCigarOp::Match, 4).expect("CIGAR")];
    let record = BorrowedAlignmentRecord::new(
        b"read1",
        1,
        Some(0),
        1,
        42,
        &cigar,
        Some(1),
        2,
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
    .expect("compact record");
    let sam = sam_borrowed_record_bytes(&record, &header, limits).expect("SAM compact record");
    assert!(
        sam.windows(b"\tdup\t2\t".len())
            .any(|window| window == b"\tdup\t2\t")
    );
    assert!(
        !sam.windows(b"\t=\t2\t".len())
            .any(|window| window == b"\t=\t2\t")
    );
}

#[test]
fn sam_writer_delegates_file_lifecycle_to_generic_io() {
    let (header, record) = fixture();
    let target = unique_path("record.sam");
    let staging = target.with_extension("sam.tmp");
    let mut writer =
        SamFileWriter::create_new(&target, &staging, &header, AlignmentRecordLimits::default())
            .expect("SAM staging");
    writer
        .write_record(&record, AlignmentRecordLimits::default())
        .expect("record");
    match writer.finish() {
        Ok(publication) => {
            assert_eq!(publication.records_written(), 1);
            assert!(
                fs::read(&target)
                    .expect("SAM target")
                    .ends_with(b"XG:Z:CT\n")
            );
            fs::remove_file(target).expect("target cleanup");
        }
        Err(error) if error.kind() == Some(std::io::ErrorKind::Unsupported) => {}
        Err(error) => panic!("publish SAM: {error}"),
    }
}

#[test]
fn reference_names_reject_sam_reserved_first_bytes_only_at_the_first_byte() {
    for (name, expected_byte) in [(b"*chr".as_slice(), b'*'), (b"=chr".as_slice(), b'=')] {
        let error = SamHeaderReference::new(0, name, 8).expect_err("reserved first byte");
        assert!(matches!(
            error,
            AlignmentRecordError::InvalidReferenceNameByte {
                ordinal: 0,
                offset: 0,
                byte: Some(byte),
            } if byte == expected_byte
        ));
    }
    assert!(SamHeaderReference::new(0, b"chr*", 8).is_ok());
    assert!(SamHeaderReference::new(0, b"chr=", 8).is_ok());
}

#[test]
fn sam_encoding_rechecks_cigar_limits_at_the_writer_boundary() {
    let (_, record) = fixture();
    let limits = AlignmentRecordLimits::default();
    let strict_cigar_text = AlignmentRecordLimits::new(
        limits.max_query_name_bytes(),
        limits.max_read_bases(),
        limits.max_quality_bytes(),
        limits.max_cigar_runs(),
        1,
        limits.max_md_bytes(),
        limits.max_optional_field_bytes(),
        limits.max_sam_line_bytes(),
        limits.max_header_references(),
        limits.max_header_name_bytes(),
        limits.max_header_bytes(),
    );
    let error = sam_record_bytes(&record, strict_cigar_text)
        .expect_err("the writer enforces its own CIGAR text cap");
    assert!(matches!(
        error,
        bsbit_hts::AlignmentRecordError::LimitExceeded {
            resource: AlignmentRecordResource::CigarTextBytes,
            ..
        }
    ));
}
