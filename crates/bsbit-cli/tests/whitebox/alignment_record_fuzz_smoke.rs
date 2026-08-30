//! Deterministic mutation smoke tests for external record metadata.

use std::mem::discriminant;

use super::{
    RecordBuildError as AlignmentRecordError, build_sam_header, build_single_alignment_record,
};
use bsbit_align::extension::VerifiedAlignment;
use bsbit_align::materialize::traceback_read_placement;
use bsbit_core::bisulfite::BisulfiteStrand;
use bsbit_core::coordinate::{ReferenceInterval, ReferenceLength};
use bsbit_core::sequence::{NormalizedSequence, normalize_dna};
use bsbit_hts::{
    AlignmentRead, AlignmentRecordError as HtsAlignmentRecordError, AlignmentRecordLimits,
    RecordMappingQuality, sam_record_bytes,
};
use bsbit_index::reference::{ContigInput, ReferenceBuildLimits, ReferenceIndex};

fn normalized(raw: &[u8]) -> NormalizedSequence {
    normalize_dna(raw).expect("fixed test sequence is valid")
}

fn reference(name: &[u8]) -> ReferenceIndex {
    ReferenceIndex::build(
        vec![ContigInput::new(name.to_vec(), normalized(b"GGACCTAA"))],
        ReferenceBuildLimits::MAX,
    )
    .expect("one nonempty contig builds")
}

fn exact_alignment(reference: &ReferenceIndex, read: &NormalizedSequence) -> VerifiedAlignment {
    let contig = reference
        .contig_by_ordinal(0)
        .expect("fixture contig exists");
    let interval = ReferenceInterval::new(
        2,
        2 + read.len(),
        ReferenceLength::new(contig.sequence().len()),
    )
    .expect("fixture interval is bounded");
    let contig_id = reference.contig_id(0).expect("fixture contig id exists");
    traceback_read_placement(
        reference,
        read,
        &contig_id,
        interval,
        BisulfiteStrand::OT,
        0,
    )
    .expect("fixture alignment materializes")
}

fn next(state: &mut u64) -> u8 {
    *state = state
        .wrapping_mul(6_364_136_223_846_793_005)
        .wrapping_add(1_442_695_040_888_963_407);
    state.to_le_bytes()[3]
}

fn valid_query_name(name: &[u8]) -> bool {
    !name.is_empty()
        && name.len() <= 254
        && name
            .iter()
            .all(|byte| matches!(byte, b'!'..=b'?' | b'A'..=b'~'))
}

fn valid_quality(quality: &[u8]) -> bool {
    quality.len() == 4 && quality.iter().all(|byte| (b'!'..=b'~').contains(byte))
}

fn valid_reference_byte(byte: u8) -> bool {
    byte.is_ascii_graphic()
        && !matches!(
            byte,
            b'\\' | b',' | b'"' | b'\'' | b'(' | b')' | b'[' | b']' | b'{' | b'}' | b'<' | b'>'
        )
}

#[test]
fn every_byte_has_the_pinned_qname_and_reference_name_classification() {
    let canonical_reference = reference(b"chr");
    let read = normalized(b"ACCT");
    let alignment = exact_alignment(&canonical_reference, &read);
    for byte in u8::MIN..=u8::MAX {
        let name = [byte];
        let built = build_single_alignment_record(
            &canonical_reference,
            &name,
            AlignmentRead::new(&read, None),
            Some(&alignment),
            RecordMappingQuality::Calibrated(40),
            AlignmentRecordLimits::default(),
        );
        assert_eq!(
            built.is_ok(),
            valid_query_name(&name),
            "QNAME byte {byte:#04x}"
        );

        let tail_name = [b'a', byte];
        let tail_reference = reference(&tail_name);
        assert_eq!(
            build_sam_header(&tail_reference, AlignmentRecordLimits::default()).is_ok(),
            valid_reference_byte(byte),
            "reference tail byte {byte:#04x}"
        );

        let first_name = [byte, b'a'];
        let first_reference = reference(&first_name);
        assert_eq!(
            build_sam_header(&first_reference, AlignmentRecordLimits::default()).is_ok(),
            valid_reference_byte(byte) && !matches!(byte, b'*' | b'='),
            "reference first byte {byte:#04x}"
        );
    }
}

#[test]
fn four_thousand_ninety_six_metadata_mutations_are_stable_and_panic_free() {
    let reference = reference(b"chr");
    let read = normalized(b"ACCT");
    let alignment = exact_alignment(&reference, &read);
    let mut state = 0x4253_4249_545f_3643_u64;
    for seed in 0..4_096_u64 {
        let name_length = usize::from(next(&mut state) % 12);
        let quality_length = usize::from(next(&mut state) % 8);
        let name = (0..name_length)
            .map(|_| next(&mut state))
            .collect::<Vec<_>>();
        let quality = (0..quality_length)
            .map(|_| next(&mut state))
            .collect::<Vec<_>>();
        let first = build_single_alignment_record(
            &reference,
            &name,
            AlignmentRead::new(&read, Some(&quality)),
            Some(&alignment),
            RecordMappingQuality::Calibrated(40),
            AlignmentRecordLimits::default(),
        );
        let second = build_single_alignment_record(
            &reference,
            &name,
            AlignmentRead::new(&read, Some(&quality)),
            Some(&alignment),
            RecordMappingQuality::Calibrated(40),
            AlignmentRecordLimits::default(),
        );
        assert_eq!(
            first.is_ok(),
            valid_query_name(&name) && valid_quality(&quality),
            "seed {seed} acceptance"
        );
        match (first, second) {
            (Ok(first), Ok(second)) => {
                assert_eq!(first, second, "seed {seed} record determinism");
                assert_eq!(
                    sam_record_bytes(&first, AlignmentRecordLimits::default())
                        .expect("accepted record serializes"),
                    sam_record_bytes(&second, AlignmentRecordLimits::default())
                        .expect("repeated accepted record serializes"),
                    "seed {seed} byte determinism"
                );
            }
            (Err(first), Err(second)) => {
                assert_eq!(discriminant(&first), discriminant(&second));
                assert_eq!(first.to_string(), second.to_string());
            }
            _ => panic!("seed {seed} repeated construction changed outcome"),
        }
    }
}

#[test]
fn qname_and_header_size_boundaries_fail_at_exact_next_byte() {
    let reference = reference(b"chr");
    let read = normalized(b"ACCT");
    let alignment = exact_alignment(&reference, &read);
    let accepted = vec![b'q'; 254];
    let rejected = vec![b'q'; 255];
    assert!(
        build_single_alignment_record(
            &reference,
            &accepted,
            AlignmentRead::new(&read, None),
            Some(&alignment),
            RecordMappingQuality::Calibrated(40),
            AlignmentRecordLimits::default(),
        )
        .is_ok()
    );
    assert!(matches!(
        build_single_alignment_record(
            &reference,
            &rejected,
            AlignmentRead::new(&read, None),
            Some(&alignment),
            RecordMappingQuality::Calibrated(40),
            AlignmentRecordLimits::default(),
        ),
        Err(AlignmentRecordError::LimitExceeded { .. })
    ));

    let header =
        build_sam_header(&reference, AlignmentRecordLimits::default()).expect("header builds");
    let exact_length = u64::try_from(
        bsbit_hts::sam_header_bytes(&header, AlignmentRecordLimits::default())
            .expect("header encodes")
            .len(),
    )
    .expect("header length is portable");
    let exact =
        AlignmentRecordLimits::new(254, 10, 10, 10, 100, 100, 100, 1_000, 10, 100, exact_length);
    let short = AlignmentRecordLimits::new(
        254,
        10,
        10,
        10,
        100,
        100,
        100,
        1_000,
        10,
        100,
        exact_length - 1,
    );
    assert!(bsbit_hts::sam_header_bytes(&header, exact).is_ok());
    assert!(matches!(
        bsbit_hts::sam_header_bytes(&header, short),
        Err(HtsAlignmentRecordError::LimitExceeded { .. })
    ));
}
