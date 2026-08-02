//! Independent byte-slice oracle for the accepted FASTA/FASTQ subset.

use std::io::Cursor;

use bsbit_core::sequence::NormalizationError;
use bsbit_hts::{
    FastaReader, FastqReader, RecordField, TextRecordError, TextRecordErrorKind, TextRecordLimits,
};

#[derive(Clone, Debug, Eq, PartialEq)]
struct RecordSnapshot {
    ordinal: u64,
    name: Vec<u8>,
    description: Vec<u8>,
    sequence: Vec<u8>,
    quality: Option<Vec<u8>>,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum SequenceErrorCategory {
    UnsupportedIupac,
    InvalidByte,
}

#[derive(Clone, Debug, Eq, PartialEq)]
enum ErrorCategory {
    UnexpectedEof,
    InvalidMarker {
        expected: u8,
        found: Option<u8>,
    },
    EmptyName,
    InvalidHeaderByte {
        byte: u8,
        column: u64,
    },
    EmptySequence,
    InvalidSequence {
        category: SequenceErrorCategory,
        byte: u8,
        column: u64,
        record_offset: u64,
    },
    PlusHeaderMismatch,
    QualityLengthMismatch {
        sequence: u64,
        quality: u64,
    },
    InvalidQualityByte {
        byte: u8,
        column: u64,
    },
}

#[derive(Clone, Debug, Eq, PartialEq)]
struct ErrorSnapshot {
    ordinal: u64,
    line: u64,
    field: RecordField,
    category: ErrorCategory,
}

#[derive(Clone, Debug, Eq, PartialEq)]
enum Outcome {
    Records(Vec<RecordSnapshot>),
    Error(ErrorSnapshot),
}

#[derive(Clone, Copy)]
struct SliceLine<'a> {
    number: u64,
    bytes: &'a [u8],
}

fn limits() -> TextRecordLimits {
    TextRecordLimits::new(1_024, 100, 128, 256, 1_024, 10_000, 1_024)
}

fn slice_lines(input: &[u8]) -> Vec<SliceLine<'_>> {
    let mut result = Vec::new();
    let mut start = 0;
    let mut number = 1_u64;
    while start < input.len() {
        let relative_newline = input[start..].iter().position(|&byte| byte == b'\n');
        let (end, next, terminated) = relative_newline
            .map_or((input.len(), input.len(), false), |relative| {
                (start + relative, start + relative + 1, true)
            });
        let mut content = &input[start..end];
        if terminated && content.last() == Some(&b'\r') {
            content = &content[..content.len() - 1];
        }
        result.push(SliceLine {
            number,
            bytes: content,
        });
        number = number.checked_add(1).expect("tiny oracle line count");
        start = next;
    }
    result
}

fn error(ordinal: u64, line: u64, field: RecordField, category: ErrorCategory) -> ErrorSnapshot {
    ErrorSnapshot {
        ordinal,
        line,
        field,
        category,
    }
}

fn oracle_header(
    line: SliceLine<'_>,
    marker: u8,
    ordinal: u64,
) -> Result<(Vec<u8>, Vec<u8>), ErrorSnapshot> {
    if line.bytes.first().copied() != Some(marker) {
        return Err(error(
            ordinal,
            line.number,
            RecordField::Header,
            ErrorCategory::InvalidMarker {
                expected: marker,
                found: line.bytes.first().copied(),
            },
        ));
    }
    let tail = &line.bytes[1..];
    let split = tail
        .iter()
        .position(|&byte| matches!(byte, b' ' | b'\t'))
        .unwrap_or(tail.len());
    let name = &tail[..split];
    if name.is_empty() {
        return Err(error(
            ordinal,
            line.number,
            RecordField::Name,
            ErrorCategory::EmptyName,
        ));
    }
    for (offset, &byte) in name.iter().enumerate() {
        if !byte.is_ascii_graphic() || byte.is_ascii_whitespace() {
            return Err(error(
                ordinal,
                line.number,
                RecordField::Name,
                ErrorCategory::InvalidHeaderByte {
                    byte,
                    column: u64::try_from(offset).expect("tiny offset") + 2,
                },
            ));
        }
    }
    let description_start = tail[split..]
        .iter()
        .position(|&byte| !matches!(byte, b' ' | b'\t'))
        .map_or(tail.len(), |offset| split + offset);
    let description = &tail[description_start..];
    for (offset, &byte) in description.iter().enumerate() {
        if !(matches!(byte, b' ' | b'\t') || byte.is_ascii_graphic()) {
            return Err(error(
                ordinal,
                line.number,
                RecordField::Description,
                ErrorCategory::InvalidHeaderByte {
                    byte,
                    column: u64::try_from(description_start + offset).expect("tiny offset") + 2,
                },
            ));
        }
    }
    Ok((name.to_vec(), description.to_vec()))
}

fn oracle_base(byte: u8) -> Result<u8, (SequenceErrorCategory, u8)> {
    match byte {
        b'A' | b'a' => Ok(b'A'),
        b'C' | b'c' => Ok(b'C'),
        b'G' | b'g' => Ok(b'G'),
        b'T' | b't' => Ok(b'T'),
        b'N' | b'n' => Ok(b'N'),
        b'R' | b'r' | b'Y' | b'y' | b'S' | b's' | b'W' | b'w' | b'K' | b'k' | b'M' | b'm'
        | b'B' | b'b' | b'D' | b'd' | b'H' | b'h' | b'V' | b'v' => {
            Err((SequenceErrorCategory::UnsupportedIupac, byte))
        }
        _ => Err((SequenceErrorCategory::InvalidByte, byte)),
    }
}

fn oracle_sequence_line(
    line: SliceLine<'_>,
    ordinal: u64,
    record_start: u64,
) -> Result<Vec<u8>, ErrorSnapshot> {
    let mut sequence = Vec::new();
    for (offset, &byte) in line.bytes.iter().enumerate() {
        let offset = u64::try_from(offset).expect("tiny offset");
        match oracle_base(byte) {
            Ok(base) => sequence.push(base),
            Err((category, byte)) => {
                return Err(error(
                    ordinal,
                    line.number,
                    RecordField::Sequence,
                    ErrorCategory::InvalidSequence {
                        category,
                        byte,
                        column: offset + 1,
                        record_offset: record_start + offset,
                    },
                ));
            }
        }
    }
    Ok(sequence)
}

fn oracle_fasta(input: &[u8]) -> Outcome {
    let lines = slice_lines(input);
    let mut cursor = 0;
    let mut records = Vec::new();
    while cursor < lines.len() {
        let ordinal = u64::try_from(records.len()).expect("tiny record count");
        let header = lines[cursor];
        cursor += 1;
        let (name, description) = match oracle_header(header, b'>', ordinal) {
            Ok(value) => value,
            Err(error) => return Outcome::Error(error),
        };
        let mut sequence = Vec::new();
        while cursor < lines.len() && lines[cursor].bytes.first() != Some(&b'>') {
            let line = lines[cursor];
            cursor += 1;
            if line.bytes.is_empty() {
                return Outcome::Error(error(
                    ordinal,
                    line.number,
                    RecordField::Sequence,
                    ErrorCategory::EmptySequence,
                ));
            }
            let record_start = u64::try_from(sequence.len()).expect("tiny sequence");
            match oracle_sequence_line(line, ordinal, record_start) {
                Ok(bases) => sequence.extend_from_slice(&bases),
                Err(error) => return Outcome::Error(error),
            }
        }
        if sequence.is_empty() {
            return Outcome::Error(error(
                ordinal,
                header.number,
                RecordField::Sequence,
                ErrorCategory::EmptySequence,
            ));
        }
        records.push(RecordSnapshot {
            ordinal,
            name,
            description,
            sequence,
            quality: None,
        });
    }
    Outcome::Records(records)
}

fn required_line<'a>(
    lines: &'a [SliceLine<'a>],
    cursor: &mut usize,
    ordinal: u64,
    field: RecordField,
) -> Result<SliceLine<'a>, ErrorSnapshot> {
    let Some(&line) = lines.get(*cursor) else {
        return Err(error(
            ordinal,
            u64::try_from(lines.len()).expect("tiny line count") + 1,
            field,
            ErrorCategory::UnexpectedEof,
        ));
    };
    *cursor += 1;
    Ok(line)
}

fn oracle_fastq(input: &[u8]) -> Outcome {
    let lines = slice_lines(input);
    let mut cursor = 0;
    let mut records = Vec::new();
    while cursor < lines.len() {
        let ordinal = u64::try_from(records.len()).expect("tiny record count");
        let header = lines[cursor];
        cursor += 1;
        let (name, description) = match oracle_header(header, b'@', ordinal) {
            Ok(value) => value,
            Err(error) => return Outcome::Error(error),
        };
        let sequence_line = match required_line(&lines, &mut cursor, ordinal, RecordField::Sequence)
        {
            Ok(line) => line,
            Err(error) => return Outcome::Error(error),
        };
        if sequence_line.bytes.is_empty() {
            return Outcome::Error(error(
                ordinal,
                sequence_line.number,
                RecordField::Sequence,
                ErrorCategory::EmptySequence,
            ));
        }
        let sequence = match oracle_sequence_line(sequence_line, ordinal, 0) {
            Ok(value) => value,
            Err(error) => return Outcome::Error(error),
        };
        let plus = match required_line(&lines, &mut cursor, ordinal, RecordField::Plus) {
            Ok(line) => line,
            Err(error) => return Outcome::Error(error),
        };
        if plus.bytes.first().copied() != Some(b'+') {
            return Outcome::Error(error(
                ordinal,
                plus.number,
                RecordField::Plus,
                ErrorCategory::InvalidMarker {
                    expected: b'+',
                    found: plus.bytes.first().copied(),
                },
            ));
        }
        if plus.bytes.len() > 1 && plus.bytes[1..] != header.bytes[1..] {
            return Outcome::Error(error(
                ordinal,
                plus.number,
                RecordField::Plus,
                ErrorCategory::PlusHeaderMismatch,
            ));
        }
        let quality = match required_line(&lines, &mut cursor, ordinal, RecordField::Quality) {
            Ok(line) => line,
            Err(error) => return Outcome::Error(error),
        };
        if quality.bytes.len() != sequence.len() {
            return Outcome::Error(error(
                ordinal,
                quality.number,
                RecordField::Quality,
                ErrorCategory::QualityLengthMismatch {
                    sequence: u64::try_from(sequence.len()).expect("tiny sequence"),
                    quality: u64::try_from(quality.bytes.len()).expect("tiny quality"),
                },
            ));
        }
        if let Some((offset, &byte)) = quality
            .bytes
            .iter()
            .enumerate()
            .find(|(_, byte)| !(33..=126).contains(*byte))
        {
            return Outcome::Error(error(
                ordinal,
                quality.number,
                RecordField::Quality,
                ErrorCategory::InvalidQualityByte {
                    byte,
                    column: u64::try_from(offset).expect("tiny offset") + 1,
                },
            ));
        }
        records.push(RecordSnapshot {
            ordinal,
            name,
            description,
            sequence,
            quality: Some(quality.bytes.to_vec()),
        });
    }
    Outcome::Records(records)
}

fn implementation_error(error: &TextRecordError) -> ErrorSnapshot {
    let category = match error.kind() {
        TextRecordErrorKind::UnexpectedEof => ErrorCategory::UnexpectedEof,
        TextRecordErrorKind::InvalidMarker { expected, found } => ErrorCategory::InvalidMarker {
            expected: *expected,
            found: *found,
        },
        TextRecordErrorKind::EmptyName => ErrorCategory::EmptyName,
        TextRecordErrorKind::InvalidHeaderByte { byte, column } => {
            ErrorCategory::InvalidHeaderByte {
                byte: *byte,
                column: *column,
            }
        }
        TextRecordErrorKind::EmptySequence => ErrorCategory::EmptySequence,
        TextRecordErrorKind::InvalidSequence {
            column,
            record_offset,
            source,
        } => {
            let (category, byte) = match source {
                NormalizationError::UnsupportedIupac { byte, .. } => {
                    (SequenceErrorCategory::UnsupportedIupac, *byte)
                }
                NormalizationError::InvalidBaseByte { byte, .. } => {
                    (SequenceErrorCategory::InvalidByte, *byte)
                }
            };
            ErrorCategory::InvalidSequence {
                category,
                byte,
                column: *column,
                record_offset: *record_offset,
            }
        }
        TextRecordErrorKind::PlusHeaderMismatch => ErrorCategory::PlusHeaderMismatch,
        TextRecordErrorKind::QualityLengthMismatch { sequence, quality } => {
            ErrorCategory::QualityLengthMismatch {
                sequence: *sequence,
                quality: *quality,
            }
        }
        TextRecordErrorKind::InvalidQualityByte { byte, column } => {
            ErrorCategory::InvalidQualityByte {
                byte: *byte,
                column: *column,
            }
        }
        other => panic!("unexpected implementation error under generous limits: {other:?}"),
    };
    ErrorSnapshot {
        ordinal: error.ordinal().get(),
        line: error.line(),
        field: error.field(),
        category,
    }
}

fn implementation_fasta(input: &[u8]) -> Outcome {
    let mut reader = FastaReader::new(Cursor::new(input), limits());
    let mut records = Vec::new();
    loop {
        match reader.next_record() {
            Ok(Some(record)) => records.push(RecordSnapshot {
                ordinal: record.ordinal().get(),
                name: record.record_name().name().to_vec(),
                description: record.record_name().description().to_vec(),
                sequence: record.sequence().to_ascii(),
                quality: None,
            }),
            Ok(None) => return Outcome::Records(records),
            Err(error) => return Outcome::Error(implementation_error(&error)),
        }
    }
}

fn implementation_fastq(input: &[u8]) -> Outcome {
    let mut reader = FastqReader::new(Cursor::new(input), limits());
    let mut records = Vec::new();
    loop {
        match reader.next_record() {
            Ok(Some(record)) => records.push(RecordSnapshot {
                ordinal: record.ordinal().get(),
                name: record.record_name().name().to_vec(),
                description: record.record_name().description().to_vec(),
                sequence: record.sequence().to_ascii(),
                quality: Some(record.quality().to_vec()),
            }),
            Ok(None) => return Outcome::Records(records),
            Err(error) => return Outcome::Error(implementation_error(&error)),
        }
    }
}

fn joined(lines: &[&[u8]], crlf: bool, final_newline: bool) -> Vec<u8> {
    let mut result = Vec::new();
    for (index, line) in lines.iter().enumerate() {
        result.extend_from_slice(line);
        if index + 1 < lines.len() || final_newline {
            result.extend_from_slice(if crlf { b"\r\n" } else { b"\n" });
        }
    }
    result
}

#[test]
fn independent_oracle_matches_exhaustive_tiny_fastq_cross_product() {
    let headers: &[&[u8]] = &[b"@r", b"@r d", b"@r\td", b"@", b">r", b"@r\0x"];
    let sequences: &[&[u8]] = &[b"A", b"a", b"N", b"R", b" ", b"", b"\xff"];
    let plus_lines: &[&[u8]] = &[b"+", b"+r", b"+r d", b"-", b""];
    let qualities: &[&[u8]] = &[b"!", b"~", b" ", b"", b"!!", b"\x7f", b"\xff"];
    let mut cases = 0_u64;
    for header in headers {
        for sequence in sequences {
            for plus in plus_lines {
                for quality in qualities {
                    for crlf in [false, true] {
                        for final_newline in [false, true] {
                            let input =
                                joined(&[*header, *sequence, *plus, *quality], crlf, final_newline);
                            assert_eq!(
                                implementation_fastq(&input),
                                oracle_fastq(&input),
                                "input bytes {input:?}"
                            );
                            cases += 1;
                        }
                    }
                }
            }
        }
    }
    assert_eq!(cases, 5_880);
}

#[test]
fn independent_oracle_matches_exhaustive_tiny_fasta_cross_product() {
    let headers: &[&[u8]] = &[b">r", b">r d", b">r\td", b">", b"@r", b">r\x7f"];
    let sequences: &[&[u8]] = &[b"A", b"aN", b"R", b" ", b"", b"\xff"];
    let continuations: &[&[u8]] = &[b"C", b"nT", b"Y", b"", b">next", b"\r"];
    let mut cases = 0_u64;
    for header in headers {
        for sequence in sequences {
            for continuation in continuations {
                for crlf in [false, true] {
                    for final_newline in [false, true] {
                        let input =
                            joined(&[*header, *sequence, *continuation], crlf, final_newline);
                        assert_eq!(
                            implementation_fasta(&input),
                            oracle_fasta(&input),
                            "input bytes {input:?}"
                        );
                        cases += 1;
                    }
                }
            }
        }
    }
    assert_eq!(cases, 864);
}

#[test]
fn independent_fasta_oracle_matches_multi_record_ordinals_and_offsets() {
    let fasta = b">a d\nAC\nN\n>b\nT\n>c\nG";
    assert_eq!(implementation_fasta(fasta), oracle_fasta(fasta));
}

#[test]
fn independent_fastq_oracle_matches_multi_record_ordinals() {
    let fastq = b"@a d\nAC\n+a d\n!~\n@b\nN\n+\n!";
    assert_eq!(implementation_fastq(fastq), oracle_fastq(fastq));
}
