//! Content-detected plain/gzip transport contracts for neutral text records.

use std::fs;
use std::io::{Cursor, Read};
use std::path::PathBuf;
use std::time::{SystemTime, UNIX_EPOCH};

use bsbit_hts::{
    Compression, DecodedFastaReader, DecodedFastqReader, DecodedPairedFastqReader, DecodedReader,
    FastaReader, FastqReader, HtsError, HtsErrorKind, PairSourceSide, RecordField,
    TextRecordErrorKind, TextRecordFormat, TextRecordLimits,
};

const FASTQ_PAYLOAD: &[u8] = b"@r1\nACGT\n+\nIIII\n";
const FASTA_PAYLOAD: &[u8] = b">chr description\r\nAC\r\nGT\n";
const GZIP_FASTQ_PAYLOAD: &[u8] = &[
    0x1f, 0x8b, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x03, 0x73, 0x28, 0x32, 0xe4, 0x72, 0x74,
    0x76, 0x0f, 0xe1, 0xd2, 0xe6, 0xf2, 0x04, 0x02, 0x2e, 0x00, 0xfe, 0x49, 0x16, 0x27, 0x10, 0x00,
    0x00, 0x00,
];

fn unique_directory(label: &str) -> PathBuf {
    let nonce = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("clock")
        .as_nanos();
    std::env::temp_dir().join(format!(
        "bsbit-hts-decoded-{label}-{}-{nonce}",
        std::process::id()
    ))
}

fn limits() -> TextRecordLimits {
    TextRecordLimits::new(1_024, 10, 128, 256, 1_024, 10_000, 1_024)
}

#[test]
fn plain_and_gzip_decode_to_identical_bytes() {
    let directory = unique_directory("decode");
    fs::create_dir(&directory).expect("directory");
    let plain = directory.join("misleading.gz");
    let gzip = directory.join("neutral.data");
    let concatenated = directory.join("concatenated.fastq.gz");
    fs::write(&plain, FASTQ_PAYLOAD).expect("plain fixture");
    fs::write(&gzip, GZIP_FASTQ_PAYLOAD).expect("gzip fixture");
    fs::write(
        &concatenated,
        [GZIP_FASTQ_PAYLOAD, GZIP_FASTQ_PAYLOAD].concat(),
    )
    .expect("concatenated gzip fixture");

    let mut plain_reader = DecodedReader::open(&plain).expect("plain opens");
    let mut gzip_reader = DecodedReader::open(&gzip).expect("gzip opens");
    let mut concatenated_reader = DecodedReader::open(&concatenated).expect("gzip opens");
    assert_eq!(plain_reader.compression(), Compression::Plain);
    assert_eq!(gzip_reader.compression(), Compression::Gzip);
    assert_eq!(concatenated_reader.compression(), Compression::Gzip);

    let mut plain_bytes = Vec::new();
    let mut gzip_bytes = Vec::new();
    let mut concatenated_bytes = Vec::new();
    plain_reader
        .read_to_end(&mut plain_bytes)
        .expect("plain decodes");
    gzip_reader
        .read_to_end(&mut gzip_bytes)
        .expect("gzip decodes");
    concatenated_reader
        .read_to_end(&mut concatenated_bytes)
        .expect("all gzip members decode");
    assert_eq!(plain_bytes, FASTQ_PAYLOAD);
    assert_eq!(gzip_bytes, FASTQ_PAYLOAD);
    assert_eq!(concatenated_bytes, [FASTQ_PAYLOAD, FASTQ_PAYLOAD].concat());
    plain_reader.close().expect("plain closes");
    gzip_reader.close().expect("gzip closes");
    concatenated_reader
        .close()
        .expect("concatenated gzip closes");

    fs::remove_dir_all(directory).expect("fixture cleanup");
}

#[test]
fn decoded_fasta_and_fastq_match_the_generic_rust_parsers() {
    let directory = unique_directory("parsed-single");
    fs::create_dir(&directory).expect("directory");
    let fasta_path = directory.join("reference.gz");
    let fastq_path = directory.join("reads.data");
    fs::write(&fasta_path, FASTA_PAYLOAD).expect("plain FASTA fixture");
    fs::write(&fastq_path, GZIP_FASTQ_PAYLOAD).expect("gzip FASTQ fixture");

    let mut reference_source =
        DecodedFastaReader::open(&fasta_path, limits()).expect("FASTA opens");
    assert_eq!(reference_source.path(), fasta_path);
    assert_eq!(reference_source.compression(), Compression::Plain);
    let decoded_reference = reference_source
        .next_record()
        .expect("FASTA parses")
        .expect("one FASTA record");
    let memory_reference = FastaReader::new(Cursor::new(FASTA_PAYLOAD), limits())
        .next_record()
        .expect("memory FASTA parses")
        .expect("one memory FASTA record");
    assert_eq!(decoded_reference, memory_reference);
    assert!(reference_source.next_record().expect("FASTA EOF").is_none());
    reference_source.close().expect("FASTA closes");

    let mut read_source = DecodedFastqReader::open(&fastq_path, limits()).expect("FASTQ opens");
    assert_eq!(read_source.path(), fastq_path);
    assert_eq!(read_source.compression(), Compression::Gzip);
    let decoded_read = read_source
        .next_record()
        .expect("FASTQ parses")
        .expect("one FASTQ record");
    let memory_read = FastqReader::new(Cursor::new(FASTQ_PAYLOAD), limits())
        .next_record()
        .expect("memory FASTQ parses")
        .expect("one memory FASTQ record");
    assert_eq!(decoded_read, memory_read);
    assert!(read_source.next_record().expect("FASTQ EOF").is_none());
    read_source.close().expect("FASTQ closes");

    fs::remove_dir_all(directory).expect("fixture cleanup");
}

#[test]
fn decoded_paired_fastq_preserves_independent_compression_and_pair_semantics() {
    let directory = unique_directory("parsed-pair");
    fs::create_dir(&directory).expect("directory");
    let first_path = directory.join("first.fastq");
    let second_path = directory.join("second.fastq.gz");
    fs::write(&first_path, GZIP_FASTQ_PAYLOAD).expect("gzip first fixture");
    fs::write(&second_path, b"@r1\nTGCA\n+\nJJJJ\n").expect("plain second fixture");

    let mut reader = DecodedPairedFastqReader::open(&first_path, &second_path, limits())
        .expect("paired sources open");
    assert_eq!(
        reader.paths(),
        (first_path.as_path(), second_path.as_path())
    );
    assert_eq!(
        reader.compressions(),
        (Compression::Gzip, Compression::Plain)
    );
    let pair = reader.next_pair().expect("pair parses").expect("one pair");
    assert_eq!(pair.first().record_name().name(), b"r1");
    assert_eq!(pair.first().sequence().to_ascii(), b"ACGT");
    assert_eq!(pair.second().record_name().name(), b"r1");
    assert_eq!(pair.second().sequence().to_ascii(), b"TGCA");
    assert!(reader.next_pair().expect("paired EOF").is_none());
    reader.close().expect("both sources close");

    let missing = directory.join("missing.fastq");
    let open_error = DecodedPairedFastqReader::open(&first_path, &missing, limits())
        .err()
        .expect("second-source open failure is returned");
    assert_eq!(open_error.path(), missing);
    assert_eq!(
        open_error.kind(),
        HtsErrorKind::Native(bsbit_hts::NativeStatus::OpenFailed)
    );

    fs::remove_dir_all(directory).expect("fixture cleanup");
}

#[test]
fn corrupt_gzip_parser_error_retains_native_cause_and_terminal_context() {
    let directory = unique_directory("parsed-corrupt");
    fs::create_dir(&directory).expect("directory");
    let path = directory.join("truncated.fastq");
    fs::write(&path, &GZIP_FASTQ_PAYLOAD[..GZIP_FASTQ_PAYLOAD.len() - 5]).expect("corrupt fixture");
    let mut reader = DecodedFastqReader::open(&path, limits()).expect("gzip header opens");

    let error = reader.next_record().expect_err("decode fails in parser");
    assert_eq!(error.format(), TextRecordFormat::Fastq);
    assert_eq!(error.side(), PairSourceSide::Single);
    assert_eq!(error.ordinal().get(), 0);
    assert_eq!(error.line(), 1);
    assert_eq!(error.field(), RecordField::Header);
    let TextRecordErrorKind::Io { source } = error.kind() else {
        panic!("expected I/O wrapper, got {:?}", error.kind());
    };
    let native = source
        .get_ref()
        .and_then(|source| source.downcast_ref::<HtsError>())
        .expect("I/O source retains HtsError");
    assert_eq!(native.path(), path);
    assert_eq!(
        native.kind(),
        HtsErrorKind::Native(bsbit_hts::NativeStatus::ReadFailed)
    );
    assert!(matches!(
        reader
            .next_record()
            .expect_err("parser remains terminal")
            .kind(),
        TextRecordErrorKind::TerminalState
    ));
    let close = reader
        .close()
        .expect_err("corrupt decoder reports its close failure separately");
    assert_eq!(close.path(), path);
    assert_eq!(
        close.kind(),
        HtsErrorKind::Native(bsbit_hts::NativeStatus::CloseFailed)
    );

    let first_path = directory.join("first.fastq");
    fs::write(&first_path, FASTQ_PAYLOAD).expect("valid first mate");
    let mut paired = DecodedPairedFastqReader::open(&first_path, &path, limits())
        .expect("paired sources open before decode");
    let paired_error = paired.next_pair().expect_err("second decoder fails");
    assert_eq!(paired_error.side(), PairSourceSide::Second);
    assert_eq!(paired_error.ordinal().get(), 0);
    assert_eq!(paired_error.line(), 1);
    assert_eq!(paired_error.field(), RecordField::Header);
    assert!(matches!(
        paired_error.kind(),
        TextRecordErrorKind::Io { .. }
    ));
    let paired_close = paired
        .close()
        .expect_err("second corrupt source reports close failure");
    assert_eq!(paired_close.path(), path);
    assert_eq!(
        paired_close.kind(),
        HtsErrorKind::Native(bsbit_hts::NativeStatus::CloseFailed)
    );

    fs::remove_dir_all(directory).expect("fixture cleanup");
}

#[test]
fn corrupt_gzip_failure_is_terminal_and_publishes_no_failing_bytes() {
    let directory = unique_directory("corrupt-gzip");
    fs::create_dir(&directory).expect("directory");
    let corrupt = directory.join("truncated.data");
    fs::write(
        &corrupt,
        &GZIP_FASTQ_PAYLOAD[..GZIP_FASTQ_PAYLOAD.len() - 5],
    )
    .expect("corrupt fixture");
    let mut reader = DecodedReader::open(&corrupt).expect("gzip header opens");
    assert_eq!(reader.compression(), Compression::Gzip);

    let mut buffer = [0_u8; 7];
    let first_error = loop {
        match reader.read_decoded(&mut buffer) {
            Ok(0) => panic!("truncated gzip reached clean EOF"),
            Ok(_) => {}
            Err(error) => break error,
        }
    };
    assert_eq!(
        first_error.kind(),
        HtsErrorKind::Native(bsbit_hts::NativeStatus::ReadFailed)
    );
    let replay = reader
        .read_decoded(&mut buffer)
        .expect_err("terminal read replays an error");
    assert_eq!(replay.kind(), first_error.kind());
    drop(reader);

    fs::remove_dir_all(directory).expect("fixture cleanup");
}

#[test]
fn unsupported_and_missing_paths_have_typed_boundaries() {
    assert_eq!(
        DecodedReader::open("-")
            .err()
            .expect("stdio rejected")
            .kind(),
        HtsErrorKind::UnsupportedPath
    );
    assert_eq!(
        DecodedReader::open("https://example.invalid/x")
            .err()
            .expect("URL rejected")
            .kind(),
        HtsErrorKind::UnsupportedPath
    );
    assert_eq!(
        DecodedReader::open("embedded\0nul")
            .err()
            .expect("NUL rejected")
            .kind(),
        HtsErrorKind::PathContainsNul
    );
    let missing = unique_directory("missing").join("none.fastq");
    let error = DecodedReader::open(&missing)
        .err()
        .expect("missing path fails natively");
    assert_eq!(
        error.kind(),
        HtsErrorKind::Native(bsbit_hts::NativeStatus::OpenFailed)
    );
    assert_eq!(error.path(), missing);
}
