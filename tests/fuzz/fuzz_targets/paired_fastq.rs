#![no_main]

use std::io::{BufReader, Cursor};

use bsbit_hts::{PairedFastqReader, PairedFastqRecord, TextRecordErrorKind};
use libfuzzer_sys::fuzz_target;

mod support;

fn parse(
    first: &[u8],
    second: &[u8],
    first_capacity: usize,
    second_capacity: usize,
) -> Result<Vec<PairedFastqRecord>, String> {
    let limits = support::text_limits();
    let mut reader = PairedFastqReader::new(
        BufReader::with_capacity(first_capacity, Cursor::new(first)),
        BufReader::with_capacity(second_capacity, Cursor::new(second)),
        limits,
    );
    let mut records = Vec::new();
    for _ in 0..=limits.max_records() {
        match reader.next_pair() {
            Ok(Some(record)) => records.push(record),
            Ok(None) => return Ok(records),
            Err(error) => {
                let diagnostic = format!("{error:?}");
                assert!(matches!(
                    reader
                        .next_pair()
                        .expect_err("paired FASTQ failure must make the reader terminal")
                        .kind(),
                    TextRecordErrorKind::TerminalState
                ));
                return Err(diagnostic);
            }
        }
    }
    panic!("bounded paired FASTQ parser did not reach EOF or a terminal error")
}

fuzz_target!(|data: &[u8]| {
    let first_control = data.first().copied().unwrap_or(0);
    let second_control = data.get(1).copied().unwrap_or(0);
    let split_control = u16::from_le_bytes([
        data.get(2).copied().unwrap_or(0),
        data.get(3).copied().unwrap_or(0),
    ]);
    let payload = data.get(4..).unwrap_or(&[]);
    let split = usize::from(split_control) % (payload.len() + 1);
    let (first_input, second_input) = payload.split_at(split);
    let first_capacity = support::buffer_capacity(first_control);
    let second_capacity = support::buffer_capacity(second_control);

    let first = parse(first_input, second_input, first_capacity, second_capacity);
    let repeated = parse(first_input, second_input, first_capacity, second_capacity);
    assert_eq!(
        first, repeated,
        "identical paired FASTQ bytes must be deterministic"
    );

    if let Ok(records) = first {
        let mut canonical_first = Vec::new();
        let mut canonical_second = Vec::new();
        for pair in &records {
            pair.first()
                .write_canonical(&mut canonical_first)
                .expect("Vec writes cannot fail");
            pair.second()
                .write_canonical(&mut canonical_second)
                .expect("Vec writes cannot fail");
        }
        let reparsed = parse(
            &canonical_first,
            &canonical_second,
            first_capacity,
            second_capacity,
        )
        .expect("canonical paired FASTQ must parse successfully");
        assert_eq!(
            records, reparsed,
            "canonical paired FASTQ must preserve semantics"
        );
    }
});
