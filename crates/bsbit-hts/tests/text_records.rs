//! Independent ground-truth and boundary tests for bounded text records.

use std::cell::Cell;
use std::io::{self, BufRead, BufReader, Cursor, Read};
use std::rc::Rc;
use std::thread;

use bsbit_core::sequence::{NormalizationError, normalize_dna};
use bsbit_hts::{
    FastaReader, FastqReader, PairSourceSide, PairedFastqReader, RecordField, TextRecordErrorKind,
    TextRecordLimits, TextRecordResource,
};

fn limits() -> TextRecordLimits {
    TextRecordLimits::new(1_024, 100, 128, 256, 1_024, 10_000, 1_024)
}

fn fasta(input: &[u8]) -> FastaReader<Cursor<&[u8]>> {
    FastaReader::new(Cursor::new(input), limits())
}

fn fastq(input: &[u8]) -> FastqReader<Cursor<&[u8]>> {
    FastqReader::new(Cursor::new(input), limits())
}

#[test]
fn empty_inputs_yield_no_records() {
    assert!(fasta(b"").next_record().expect("empty FASTA").is_none());
    assert!(fastq(b"").next_record().expect("empty FASTQ").is_none());
}

#[test]
fn fastq_batch_is_contiguous_and_fully_validated() {
    let mut reader = FastqReader::new(
        Cursor::new(b"@read description bytes\nACGT\n+\n!~!!\n"),
        limits(),
    );
    let batch = reader.next_batch(8).expect("validated batch parses");
    let record = batch.get(0).expect("record is present");
    assert_eq!(record.name(), b"read");
    assert_eq!(
        record
            .sequence()
            .iter()
            .map(|base| base.as_ascii())
            .collect::<Vec<_>>(),
        b"ACGT"
    );
    assert_eq!(record.quality(), b"!~!!");

    let header_error = FastqReader::new(
        Cursor::new(b"@read bad\x1fdescription\nA\n+\n!\n"),
        limits(),
    )
    .next_batch(1)
    .expect_err("control byte in description is rejected");
    assert_eq!(header_error.field(), RecordField::Description);
    assert!(matches!(
        header_error.kind(),
        TextRecordErrorKind::InvalidHeaderByte { byte: 0x1f, .. }
    ));

    let name_error = FastqReader::new(Cursor::new(b"@re\x1fad\nA\n+\n!\n"), limits())
        .next_batch(1)
        .expect_err("control byte in name is rejected");
    assert_eq!(name_error.field(), RecordField::Name);
    assert!(matches!(
        name_error.kind(),
        TextRecordErrorKind::InvalidHeaderByte { byte: 0x1f, .. }
    ));

    let quality_error = FastqReader::new(Cursor::new(b"@read\nA\n+\n \n"), limits())
        .next_batch(1)
        .expect_err("quality byte below Phred+33 is rejected");
    assert_eq!(quality_error.field(), RecordField::Quality);
    assert!(matches!(
        quality_error.kind(),
        TextRecordErrorKind::InvalidQualityByte { byte: b' ', .. }
    ));
}

#[test]
fn fastq_batch_and_owned_record_paths_have_identical_validation_outcomes() {
    let cases: &[&[u8]] = &[
        b"@read description\nAcgTN\n+\n!#$%&\n",
        b"read\nA\n+\n!\n",
        b"@\nA\n+\n!\n",
        b"@re\x1fad\nA\n+\n!\n",
        b"@read bad\x1fdescription\nA\n+\n!\n",
        b"@read\n\n+\n\n",
        b"@read\nR\n+\n!\n",
        b"@read description\nA\n+other\n!\n",
        b"@read\nAA\n+\n!\n",
        b"@read\nA\n+\n \n",
    ];

    for &input in cases {
        let owned = FastqReader::new(Cursor::new(input), limits()).next_record();
        let batch = FastqReader::new(Cursor::new(input), limits()).next_batch(1);
        match (owned, batch) {
            (Ok(Some(owned)), Ok(batch)) => {
                let borrowed = batch.get(0).expect("batch record");
                assert_eq!(borrowed.ordinal(), owned.ordinal(), "input {input:?}");
                assert_eq!(
                    borrowed.name(),
                    owned.record_name().name(),
                    "input {input:?}"
                );
                assert_eq!(
                    borrowed
                        .sequence()
                        .iter()
                        .map(|base| base.as_ascii())
                        .collect::<Vec<_>>(),
                    owned.sequence().to_ascii(),
                    "input {input:?}"
                );
                assert_eq!(borrowed.quality(), owned.quality(), "input {input:?}");
            }
            (Err(owned), Err(batch)) => {
                assert_eq!(batch.format(), owned.format(), "input {input:?}");
                assert_eq!(batch.side(), owned.side(), "input {input:?}");
                assert_eq!(batch.ordinal(), owned.ordinal(), "input {input:?}");
                assert_eq!(batch.line(), owned.line(), "input {input:?}");
                assert_eq!(batch.field(), owned.field(), "input {input:?}");
                assert_eq!(
                    format!("{:?}", batch.kind()),
                    format!("{:?}", owned.kind()),
                    "input {input:?}"
                );
            }
            (owned, batch) => {
                panic!("record/batch outcome diverged for {input:?}: {owned:?} {batch:?}")
            }
        }
    }
}

#[test]
fn fastq_batch_capacity_failure_is_typed_and_terminal() {
    let mut reader = FastqReader::new(Cursor::new(b""), limits());
    let error = reader
        .next_batch(usize::MAX)
        .expect_err("impossible batch capacity is rejected");
    assert!(matches!(
        error.kind(),
        TextRecordErrorKind::AllocationFailed {
            allocation: bsbit_hts::TextRecordAllocation::Record,
            ..
        }
    ));
    assert!(matches!(
        reader
            .next_batch(1)
            .expect_err("allocation failure is terminal")
            .kind(),
        TextRecordErrorKind::TerminalState
    ));
}

#[test]
fn fasta_retains_identity_and_normalizes_multiline_crlf() {
    let mut reader = fasta(b">chr1\tfirst description\r\nacN\r\ngt\r\n>chr2\nT");
    let first = reader
        .next_record()
        .expect("valid first record")
        .expect("first record");
    assert_eq!(first.ordinal().get(), 0);
    assert_eq!(first.record_name().name(), b"chr1");
    assert_eq!(first.record_name().description(), b"first description");
    assert_eq!(first.sequence().to_ascii(), b"ACNGT");

    let second = reader
        .next_record()
        .expect("valid second record")
        .expect("second record");
    assert_eq!(second.ordinal().get(), 1);
    assert_eq!(second.sequence().to_ascii(), b"T");
    assert!(reader.next_record().expect("end of FASTA").is_none());
}

#[test]
fn fasta_retains_name_and_sequence_without_index_coupling() {
    let record = fasta(b">contig description\nAC\n")
        .next_record()
        .expect("valid record")
        .expect("one record");
    assert_eq!(record.record_name().name(), b"contig");
    assert_eq!(record.record_name().description(), b"description");
    assert_eq!(record.sequence().to_ascii(), b"AC");
}

#[test]
fn fasta_reports_concatenated_and_physical_sequence_offsets() {
    let error = fasta(b">x\nAC\nTz\n")
        .next_record()
        .expect_err("invalid sequence byte");
    assert_eq!(error.line(), 3);
    assert_eq!(error.field(), RecordField::Sequence);
    match error.kind() {
        TextRecordErrorKind::InvalidSequence {
            column,
            record_offset,
            source,
        } => {
            assert_eq!(*column, 2);
            assert_eq!(*record_offset, 3);
            assert_eq!(
                *source,
                NormalizationError::InvalidBaseByte {
                    byte: b'z',
                    offset: 1
                }
            );
        }
        other => panic!("unexpected error: {other:?}"),
    }
}

#[test]
fn blank_fasta_sequence_lines_are_not_skipped() {
    let error = fasta(b">x\n\nA\n")
        .next_record()
        .expect_err("blank sequence line");
    assert_eq!(error.line(), 2);
    assert!(matches!(error.kind(), TextRecordErrorKind::EmptySequence));
}

#[test]
fn fasta_header_without_sequence_fails_closed() {
    let error = fasta(b">x\n>y\nA\n")
        .next_record()
        .expect_err("first record has no bases");
    assert_eq!(error.line(), 1);
    assert!(matches!(error.kind(), TextRecordErrorKind::EmptySequence));
}

#[test]
fn crlf_does_not_consume_line_content_budget_even_when_split() {
    let limits = TextRecordLimits::new(2, 1, 1, 0, 2, 2, 2);
    let input = b">x\r\nAC\r\n";
    let buffered = BufReader::with_capacity(1, Cursor::new(input));
    let record = FastaReader::new(buffered, limits)
        .next_record()
        .expect("CRLF at exact cap")
        .expect("one record");
    assert_eq!(record.sequence().to_ascii(), b"AC");
}

#[test]
fn bare_carriage_return_is_field_data_and_is_rejected() {
    let error = fasta(b">x\nA\r")
        .next_record()
        .expect_err("bare CR is invalid sequence data");
    assert_eq!(error.line(), 2);
    assert!(matches!(
        error.kind(),
        TextRecordErrorKind::InvalidSequence {
            source: NormalizationError::InvalidBaseByte { byte: b'\r', .. },
            ..
        }
    ));
}

#[test]
fn line_limit_failure_is_terminal() {
    let limits = TextRecordLimits::new(2, 10, 10, 10, 10, 10, 10);
    let mut reader = FastaReader::new(Cursor::new(b">abc\nA\n"), limits);
    let error = reader.next_record().expect_err("header exceeds line cap");
    assert!(matches!(
        error.kind(),
        TextRecordErrorKind::LimitExceeded {
            resource: TextRecordResource::LineBytes,
            observed: 4,
            limit: 2,
        }
    ));
    let terminal = reader.next_record().expect_err("terminal state");
    assert!(matches!(
        terminal.kind(),
        TextRecordErrorKind::TerminalState
    ));
    assert_eq!(terminal.line(), error.line());
}

#[test]
fn header_controls_non_ascii_and_empty_names_are_rejected() {
    for input in [
        b">\nA\n".as_slice(),
        b"> name\nA\n".as_slice(),
        b">x\0y\nA\n".as_slice(),
        b">x\x7fy\nA\n".as_slice(),
        b">x\xffy\nA\n".as_slice(),
        b">x d\0e\nA\n".as_slice(),
    ] {
        let error = fasta(input).next_record().expect_err("invalid header");
        assert!(matches!(
            error.kind(),
            TextRecordErrorKind::EmptyName | TextRecordErrorKind::InvalidHeaderByte { .. }
        ));
    }
}

#[test]
fn fastq_accepts_empty_or_exact_plus_suffix_and_final_line_without_lf() {
    let input = b"@r1 description\nacn\n+\n!~!\n@r2\tsecond\r\nT\r\n+r2\tsecond\r\nI";
    let mut reader = fastq(input);
    let first = reader
        .next_record()
        .expect("valid first FASTQ")
        .expect("first record");
    assert_eq!(first.record_name().name(), b"r1");
    assert_eq!(first.record_name().description(), b"description");
    assert_eq!(first.sequence().to_ascii(), b"ACN");
    assert_eq!(first.quality(), b"!~!");

    let second = reader
        .next_record()
        .expect("valid second FASTQ")
        .expect("second record");
    assert_eq!(second.record_name().description(), b"second");
    assert_eq!(second.quality(), b"I");
    assert!(reader.next_record().expect("end of FASTQ").is_none());
}

#[test]
fn canonical_fasta_round_trip_preserves_record_semantics() {
    let fasta_record = fasta(b">ref\tretained description\r\nacN\r\nGT")
        .next_record()
        .expect("source FASTA")
        .expect("one FASTA record");
    let mut canonical_fasta = Vec::new();
    fasta_record
        .write_canonical(&mut canonical_fasta)
        .expect("write canonical FASTA");
    assert_eq!(canonical_fasta, b">ref retained description\nACNGT\n");
    let reparsed_fasta = fasta(&canonical_fasta)
        .next_record()
        .expect("canonical FASTA")
        .expect("one canonical FASTA record");
    assert_eq!(reparsed_fasta, fasta_record);
}

#[test]
fn canonical_fastq_round_trip_preserves_record_semantics() {
    let fastq_record =
        fastq(b"@read\tretained description\r\nacn\r\n+read\tretained description\r\n!~!\r\n")
            .next_record()
            .expect("source FASTQ")
            .expect("one FASTQ record");
    let mut canonical_fastq = Vec::new();
    fastq_record
        .write_canonical(&mut canonical_fastq)
        .expect("write canonical FASTQ");
    assert_eq!(
        canonical_fastq,
        b"@read retained description\nACN\n+\n!~!\n"
    );
    let reparsed_fastq = fastq(&canonical_fastq)
        .next_record()
        .expect("canonical FASTQ")
        .expect("one canonical FASTQ record");
    assert_eq!(reparsed_fastq, fastq_record);
}

#[test]
fn fastq_plus_suffix_compares_the_entire_raw_header_tail() {
    let mut valid = fastq(b"@r  description  \nA\n+r  description  \n!\n");
    assert!(valid.next_record().expect("valid exact suffix").is_some());

    let error = fastq(b"@r  description\nA\n+r description\n!\n")
        .next_record()
        .expect_err("spacing is part of exact suffix");
    assert!(matches!(
        error.kind(),
        TextRecordErrorKind::PlusHeaderMismatch
    ));
}

#[test]
fn every_fastq_truncation_reports_the_missing_field() {
    let fixtures: [(&[u8], RecordField, u64); 3] = [
        (b"@r\n", RecordField::Sequence, 2),
        (b"@r\nA\n", RecordField::Plus, 3),
        (b"@r\nA\n+\n", RecordField::Quality, 4),
    ];
    for (input, field, line) in fixtures {
        let error = fastq(input).next_record().expect_err("truncated FASTQ");
        assert_eq!(error.field(), field);
        assert_eq!(error.line(), line);
        assert!(matches!(error.kind(), TextRecordErrorKind::UnexpectedEof));
    }
}

#[test]
fn fastq_rejects_wrong_markers_empty_sequence_and_wrapping() {
    let cases = [
        b">r\nA\n+\n!\n".as_slice(),
        b"@r\n\n+\n\n".as_slice(),
        b"@r\nA\nA\n+\n!\n".as_slice(),
    ];
    for input in cases {
        assert!(fastq(input).next_record().is_err());
    }
}

#[test]
fn fastq_quality_checks_length_before_byte_range() {
    let short = fastq(b"@r\nAC\n+\n \n")
        .next_record()
        .expect_err("short quality");
    assert!(matches!(
        short.kind(),
        TextRecordErrorKind::QualityLengthMismatch {
            sequence: 2,
            quality: 1
        }
    ));

    let invalid = fastq(b"@r\nA\n+\n \n")
        .next_record()
        .expect_err("quality byte 32");
    assert!(matches!(
        invalid.kind(),
        TextRecordErrorKind::InvalidQualityByte {
            byte: 32,
            column: 1
        }
    ));
    let valid = fastq(b"@r\nAA\n+\n!~\n")
        .next_record()
        .expect("boundary quality bytes")
        .expect("one record");
    assert_eq!(valid.quality(), b"!~");
}

#[test]
fn configured_limits_accept_boundaries_and_reject_the_next_value() {
    let exact = TextRecordLimits::new(4, 1, 1, 1, 2, 2, 2);
    let mut reader = FastqReader::new(Cursor::new(b"@r d\nAC\n+\n!!\n"), exact);
    assert!(reader.next_record().expect("all exact limits").is_some());
    assert!(
        reader
            .next_record()
            .expect("EOF at exact record cap")
            .is_none()
    );

    let name_error = FastaReader::new(
        Cursor::new(b">xy\nA\n"),
        TextRecordLimits::new(3, 1, 1, 0, 1, 1, 1),
    )
    .next_record()
    .expect_err("name cap");
    assert!(matches!(
        name_error.kind(),
        TextRecordErrorKind::LimitExceeded {
            resource: TextRecordResource::NameBytes,
            observed: 2,
            limit: 1
        }
    ));

    let total_error = fastq_with_limits(
        b"@r\nAC\nnot-plus\n!!\n",
        TextRecordLimits::new(16, 1, 1, 0, 2, 1, 2),
    )
    .next_record()
    .expect_err("total cap has priority once sequence is known");
    assert!(matches!(
        total_error.kind(),
        TextRecordErrorKind::LimitExceeded {
            resource: TextRecordResource::TotalBases,
            observed: 2,
            limit: 1
        }
    ));
}

#[test]
fn every_logical_limit_is_reported_at_its_own_boundary() {
    let cases = [
        (
            b">r dd\nA\n".as_slice(),
            TextRecordLimits::new(8, 1, 1, 1, 1, 1, 1),
            TextRecordResource::DescriptionBytes,
            2,
            1,
        ),
        (
            b">r\nAC\n".as_slice(),
            TextRecordLimits::new(8, 1, 1, 0, 1, 2, 2),
            TextRecordResource::BasesPerRecord,
            2,
            1,
        ),
    ];
    for (input, limits, resource, observed, limit) in cases {
        let error = FastaReader::new(Cursor::new(input), limits)
            .next_record()
            .expect_err("configured cap must fail");
        assert!(matches!(
            error.kind(),
            TextRecordErrorKind::LimitExceeded {
                resource: actual,
                observed: actual_observed,
                limit: actual_limit,
            } if *actual == resource && *actual_observed == observed && *actual_limit == limit
        ));
    }

    let quality = FastqReader::new(
        Cursor::new(b"@r\nAC\n+\n!!\n"),
        TextRecordLimits::new(8, 1, 1, 0, 2, 2, 1),
    )
    .next_record()
    .expect_err("quality cap");
    assert!(matches!(
        quality.kind(),
        TextRecordErrorKind::LimitExceeded {
            resource: TextRecordResource::QualityBytes,
            observed: 2,
            limit: 1,
        }
    ));

    let records = FastaReader::new(
        Cursor::new(b">r\nA\n"),
        TextRecordLimits::new(8, 0, 1, 0, 1, 1, 1),
    )
    .next_record()
    .expect_err("zero record cap");
    assert!(matches!(
        records.kind(),
        TextRecordErrorKind::LimitExceeded {
            resource: TextRecordResource::Records,
            observed: 1,
            limit: 0,
        }
    ));
}

fn fastq_with_limits(input: &[u8], limits: TextRecordLimits) -> FastqReader<Cursor<&[u8]>> {
    FastqReader::new(Cursor::new(input), limits)
}

#[test]
fn record_limit_rejects_a_real_following_record() {
    let limits = TextRecordLimits::new(8, 1, 2, 0, 1, 2, 1);
    let mut reader = FastqReader::new(Cursor::new(b"@a\nA\n+\n!\n@b\nC\n+\n!\n"), limits);
    assert!(reader.next_record().expect("first record").is_some());
    let error = reader.next_record().expect_err("second record exceeds cap");
    assert!(matches!(
        error.kind(),
        TextRecordErrorKind::LimitExceeded {
            resource: TextRecordResource::Records,
            observed: 2,
            limit: 1
        }
    ));
}

#[test]
fn chunked_and_contiguous_inputs_have_identical_results() {
    let input = b"@r description\r\nacn\r\n+r description\r\n!~!\r\n";
    let contiguous = FastqReader::new(Cursor::new(input), limits())
        .next_record()
        .expect("contiguous parse");
    for capacity in 1..=7 {
        let chunked = FastqReader::new(
            BufReader::with_capacity(capacity, Cursor::new(input)),
            limits(),
        )
        .next_record()
        .expect("chunked parse");
        assert_eq!(chunked, contiguous, "buffer capacity {capacity}");
    }
}

#[test]
fn line_endings_and_fasta_rewrapping_are_semantic_invariants() {
    let fasta_variants: [&[u8]; 4] = [
        b">r description\nACNGT\n",
        b">r description\r\nACNGT\r\n",
        b">r description\nAC\nNG\nT",
        b">r description\r\nA\r\nCNGT",
    ];
    let expected = fasta(fasta_variants[0])
        .next_record()
        .expect("reference FASTA")
        .expect("one reference record");
    for input in fasta_variants {
        let actual = fasta(input)
            .next_record()
            .expect("variant FASTA")
            .expect("one variant record");
        assert_eq!(actual, expected);
    }
}

#[test]
fn fastq_line_endings_and_final_newline_are_semantic_invariants() {
    let fastq_variants: [&[u8]; 3] = [
        b"@r description\nACN\n+\n!~!\n",
        b"@r description\r\nACN\r\n+\r\n!~!\r\n",
        b"@r description\nACN\n+\n!~!",
    ];
    let expected = fastq(fastq_variants[0])
        .next_record()
        .expect("reference FASTQ")
        .expect("one reference record");
    for input in fastq_variants {
        let actual = fastq(input)
            .next_record()
            .expect("variant FASTQ")
            .expect("one variant record");
        assert_eq!(actual, expected);
    }
}

#[test]
fn paired_fastq_accepts_exact_and_slash_names_without_clipping() {
    let first = b"@same\nA\n+\n!\n@read/1 one\nACG\n+\n!!!\n";
    let second = b"@same two\nTT\n+\n!!\n@read/2 two\nT\n+\n!\n";
    let mut reader = PairedFastqReader::new(Cursor::new(first), Cursor::new(second), limits());
    let exact = reader
        .next_pair()
        .expect("exact-name pair")
        .expect("first pair");
    assert_eq!(exact.first().sequence().len(), 1);
    assert_eq!(exact.second().sequence().len(), 2);
    let suffixed = reader
        .next_pair()
        .expect("slash pair")
        .expect("second pair");
    assert_eq!(suffixed.first().record_name().name(), b"read/1");
    assert_eq!(suffixed.second().record_name().name(), b"read/2");
    assert_eq!(suffixed.shared_name(), b"read");
    assert_eq!(suffixed.first().sequence().len(), 3);
    assert_eq!(suffixed.second().sequence().len(), 1);
}

#[test]
fn paired_fastq_reports_name_and_count_mismatches_then_stops() {
    let mut names = PairedFastqReader::new(
        Cursor::new(b"@a\nA\n+\n!\n"),
        Cursor::new(b"@b\nA\n+\n!\n"),
        limits(),
    );
    let error = names.next_pair().expect_err("name mismatch");
    match error.kind() {
        TextRecordErrorKind::PairNameMismatch { first, second } => {
            assert_eq!(&**first, b"a");
            assert_eq!(&**second, b"b");
        }
        other => panic!("unexpected error: {other:?}"),
    }
    assert!(matches!(
        names.next_pair().expect_err("paired terminal").kind(),
        TextRecordErrorKind::TerminalState
    ));

    let mut count =
        PairedFastqReader::new(Cursor::new(b"@a\nA\n+\n!\n"), Cursor::new(b""), limits());
    let error = count.next_pair().expect_err("second source ended early");
    assert_eq!(error.side(), PairSourceSide::Second);
    assert!(matches!(
        error.kind(),
        TextRecordErrorKind::PairCountMismatch {
            missing: PairSourceSide::Second
        }
    ));
}

#[test]
fn paired_fastq_reports_first_source_ending_early() {
    let mut count =
        PairedFastqReader::new(Cursor::new(b""), Cursor::new(b"@a\nA\n+\n!\n"), limits());
    let error = count.next_pair().expect_err("first source ended early");
    assert_eq!(error.side(), PairSourceSide::First);
    assert!(matches!(
        error.kind(),
        TextRecordErrorKind::PairCountMismatch {
            missing: PairSourceSide::First
        }
    ));
}

#[test]
fn paired_fastq_rejects_reversed_or_one_sided_mate_suffixes() {
    for (first_name, second_name) in [
        (b"read/2".as_slice(), b"read/1".as_slice()),
        (b"read/1".as_slice(), b"read".as_slice()),
        (b"read".as_slice(), b"read/2".as_slice()),
    ] {
        let mut first = Vec::new();
        first.extend_from_slice(b"@");
        first.extend_from_slice(first_name);
        first.extend_from_slice(b"\nA\n+\n!\n");
        let mut second = Vec::new();
        second.extend_from_slice(b"@");
        second.extend_from_slice(second_name);
        second.extend_from_slice(b"\nA\n+\n!\n");
        let error = PairedFastqReader::new(Cursor::new(first), Cursor::new(second), limits())
            .next_pair()
            .expect_err("incompatible suffixes");
        assert!(matches!(
            error.kind(),
            TextRecordErrorKind::PairNameMismatch { .. }
        ));
    }
}

#[derive(Debug)]
struct FailingBufRead<'a> {
    bytes: &'a [u8],
    position: usize,
    fail_at: usize,
    fill_calls: Rc<Cell<u64>>,
}

impl Read for FailingBufRead<'_> {
    fn read(&mut self, output: &mut [u8]) -> io::Result<usize> {
        let available = self.fill_buf()?;
        let length = available.len().min(output.len());
        output[..length].copy_from_slice(&available[..length]);
        self.consume(length);
        Ok(length)
    }
}

impl BufRead for FailingBufRead<'_> {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        self.fill_calls.set(self.fill_calls.get().saturating_add(1));
        if self.position >= self.fail_at {
            return Err(io::Error::other("injected read failure"));
        }
        let end = self.bytes.len().min(self.fail_at);
        Ok(&self.bytes[self.position..end])
    }

    fn consume(&mut self, amount: usize) {
        self.position = self.position.saturating_add(amount).min(self.bytes.len());
    }
}

#[test]
fn injected_io_errors_preserve_field_context_and_terminal_reads_stop() {
    let input = b"@r\nAC\n+\n!!\n";
    let cases = [
        (0, RecordField::Header, 1),
        (1, RecordField::Header, 1),
        (3, RecordField::Sequence, 2),
        (4, RecordField::Sequence, 2),
        (6, RecordField::Plus, 3),
        (7, RecordField::Plus, 3),
        (8, RecordField::Quality, 4),
        (9, RecordField::Quality, 4),
    ];
    for (fail_at, field, line) in cases {
        let fill_calls = Rc::new(Cell::new(0));
        let source = FailingBufRead {
            bytes: input,
            position: 0,
            fail_at,
            fill_calls: Rc::clone(&fill_calls),
        };
        let mut reader = FastqReader::new(source, limits());
        let error = reader.next_record().expect_err("injected I/O error");
        assert_eq!(error.field(), field, "fail offset {fail_at}");
        assert_eq!(error.line(), line, "fail offset {fail_at}");
        assert!(matches!(error.kind(), TextRecordErrorKind::Io { .. }));
        let calls_after_error = fill_calls.get();
        assert!(matches!(
            reader.next_record().expect_err("terminal parser").kind(),
            TextRecordErrorKind::TerminalState
        ));
        assert_eq!(fill_calls.get(), calls_after_error);
    }
}

#[test]
fn local_fallible_normalization_matches_core_for_every_embeddable_byte() {
    for byte in u8::MIN..=u8::MAX {
        if matches!(byte, b'\n' | b'\r') {
            continue;
        }
        let input = [b'@', b'r', b'\n', byte, b'\n', b'+', b'\n', b'!', b'\n'];
        let parsed = fastq(&input).next_record();
        match normalize_dna(&[byte]) {
            Ok(expected) => {
                let record = parsed.expect("same accepted alphabet").expect("one record");
                assert_eq!(record.sequence(), &expected, "byte 0x{byte:02X}");
            }
            Err(expected) => {
                let error = parsed.expect_err("same rejected alphabet");
                match error.kind() {
                    TextRecordErrorKind::InvalidSequence { source, .. } => {
                        assert_eq!(*source, expected, "byte 0x{byte:02X}");
                    }
                    other => panic!("byte 0x{byte:02X}: unexpected error {other:?}"),
                }
            }
        }
    }
}

#[test]
fn independent_threads_parse_identically() {
    let input = b"@r/1 d\nACNT\n+\n!!!!\n".to_vec();
    let expected = fastq(&input)
        .next_record()
        .expect("reference parse")
        .expect("one record");
    let handles: Vec<_> = (0..8)
        .map(|_| {
            let input = input.clone();
            thread::spawn(move || {
                FastqReader::new(Cursor::new(input), limits())
                    .next_record()
                    .expect("thread parse")
                    .expect("one record")
            })
        })
        .collect();
    for handle in handles {
        assert_eq!(handle.join().expect("thread completed"), expected);
    }
}
