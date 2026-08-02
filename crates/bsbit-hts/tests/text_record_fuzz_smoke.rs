//! Deterministic mutation smoke tests for parser state and panic safety.

use std::io::{BufReader, Cursor};

use bsbit_hts::{
    FastaReader, FastqReader, PairedFastqReader, TextRecordErrorKind, TextRecordLimits,
};

const DEFAULT_CASES: u64 = 4_096;
const MAX_CASES: u64 = 1_000_000;

fn case_count() -> u64 {
    let Some(raw) = std::env::var_os("BSBIT_TEXT_FUZZ_CASES") else {
        return DEFAULT_CASES;
    };
    let cases = raw
        .to_str()
        .expect("BSBIT_TEXT_FUZZ_CASES must be UTF-8")
        .parse::<u64>()
        .expect("BSBIT_TEXT_FUZZ_CASES must be an integer");
    assert!(cases > 0, "BSBIT_TEXT_FUZZ_CASES must be positive");
    assert!(
        cases <= MAX_CASES,
        "BSBIT_TEXT_FUZZ_CASES exceeds the bounded soak cap"
    );
    cases
}

fn limits() -> TextRecordLimits {
    TextRecordLimits::new(96, 8, 32, 48, 64, 256, 64)
}

fn bytes(seed: u64) -> Vec<u8> {
    let mut state = seed ^ 0x4253_4249_545f_494f;
    let length = usize::try_from(state % 128).expect("bounded fixture length");
    let mut result = Vec::with_capacity(length);
    for _ in 0..length {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        result.push(state.to_le_bytes()[0]);
    }
    result
}

fn exhaust_fasta(input: &[u8], capacity: usize) {
    let source = BufReader::with_capacity(capacity, Cursor::new(input));
    let mut reader = FastaReader::new(source, limits());
    for _ in 0..=limits().max_records() {
        match reader.next_record() {
            Ok(Some(_)) => {}
            Ok(None) => return,
            Err(_) => {
                assert!(matches!(
                    reader
                        .next_record()
                        .expect_err("stable terminal state")
                        .kind(),
                    TextRecordErrorKind::TerminalState
                ));
                return;
            }
        }
    }
    panic!("bounded FASTA parser did not reach a terminal outcome");
}

fn exhaust_fastq(input: &[u8], capacity: usize) {
    let source = BufReader::with_capacity(capacity, Cursor::new(input));
    let mut reader = FastqReader::new(source, limits());
    for _ in 0..=limits().max_records() {
        match reader.next_record() {
            Ok(Some(_)) => {}
            Ok(None) => return,
            Err(_) => {
                assert!(matches!(
                    reader
                        .next_record()
                        .expect_err("stable terminal state")
                        .kind(),
                    TextRecordErrorKind::TerminalState
                ));
                return;
            }
        }
    }
    panic!("bounded FASTQ parser did not reach a terminal outcome");
}

#[test]
fn deterministic_arbitrary_bytes_reach_stable_single_source_outcomes() {
    for seed in 0..case_count() {
        let input = bytes(seed);
        let capacity = usize::try_from(seed % 17 + 1).expect("bounded buffer capacity");
        exhaust_fasta(&input, capacity);
        exhaust_fastq(&input, capacity);
    }
}

#[test]
fn deterministic_arbitrary_pairs_reach_stable_synchronized_outcomes() {
    for seed in 0..case_count() {
        let first = bytes(seed);
        let second = bytes(seed ^ 0x9e37_79b9_7f4a_7c15);
        let mut reader = PairedFastqReader::new(
            BufReader::with_capacity(1, Cursor::new(first)),
            BufReader::with_capacity(7, Cursor::new(second)),
            limits(),
        );
        for _ in 0..=limits().max_records() {
            match reader.next_pair() {
                Ok(Some(_)) => {}
                Ok(None) => break,
                Err(_) => {
                    assert!(matches!(
                        reader
                            .next_pair()
                            .expect_err("stable paired terminal")
                            .kind(),
                        TextRecordErrorKind::TerminalState
                    ));
                    break;
                }
            }
        }
    }
}
