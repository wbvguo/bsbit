#![no_main]

use std::io::{BufReader, Cursor};

use bsbit_hts::{FastqReader, FastqRecord, TextRecordErrorKind};
use libfuzzer_sys::fuzz_target;

mod support;

fn parse(input: &[u8], capacity: usize) -> Result<Vec<FastqRecord>, String> {
    let limits = support::text_limits();
    let mut reader = FastqReader::new(
        BufReader::with_capacity(capacity, Cursor::new(input)),
        limits,
    );
    let mut records = Vec::new();
    for _ in 0..=limits.max_records() {
        match reader.next_record() {
            Ok(Some(record)) => records.push(record),
            Ok(None) => return Ok(records),
            Err(error) => {
                let diagnostic = format!("{error:?}");
                assert!(matches!(
                    reader
                        .next_record()
                        .expect_err("FASTQ failure must make the reader terminal")
                        .kind(),
                    TextRecordErrorKind::TerminalState
                ));
                return Err(diagnostic);
            }
        }
    }
    panic!("bounded FASTQ parser did not reach EOF or a terminal error")
}

fuzz_target!(|data: &[u8]| {
    let (control, input) = data.split_first().unwrap_or((&0, &[]));
    let capacity = support::buffer_capacity(*control);
    let first = parse(input, capacity);
    let repeated = parse(input, capacity);
    assert_eq!(
        first, repeated,
        "identical FASTQ bytes must be deterministic"
    );

    if let Ok(records) = first {
        let mut canonical = Vec::new();
        for record in &records {
            record
                .write_canonical(&mut canonical)
                .expect("Vec writes cannot fail");
        }
        let reparsed =
            parse(&canonical, capacity).expect("canonical FASTQ must parse successfully");
        assert_eq!(records, reparsed, "canonical FASTQ must preserve semantics");
    }
});
