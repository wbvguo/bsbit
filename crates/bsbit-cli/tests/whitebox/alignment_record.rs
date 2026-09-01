//! Ground-truth tests for immutable records and canonical SAM text.

use std::io::{self, Write};
use std::sync::Arc;
use std::thread;

use super::record_fixture::{SingleFixture, single_fixture};
use super::{
    PairedRecordComposer, RecordBuildError as AlignmentRecordError, append_u64,
    bismark_methylation_call, build_indexed_single_alignment_record, build_sam_header,
    build_single_alignment_record, build_single_alignment_record_with_auxiliary_mode,
    checked_add_resource, decimal_digits, storage_len,
    try_build_indexed_single_ungapped_alignment_record,
};
use bsbit_align::extension::VerifiedAlignment;
use bsbit_align::materialize::traceback_read_placement;
use bsbit_core::alphabet::Base;
use bsbit_core::bisulfite::{AlignmentOrientation, BisulfiteStrand, CytosineStrand};
use bsbit_core::coordinate::{ReferenceInterval, ReferenceLength};
use bsbit_core::sequence::{NormalizedSequence, normalize_dna};
use bsbit_hts::{
    AlignmentAuxiliaryMode, AlignmentCigarOp, AlignmentPlacement, AlignmentRead, AlignmentRecord,
    AlignmentRecordBatch, AlignmentRecordError as HtsAlignmentRecordError, AlignmentRecordLimits,
    AlignmentRecordResource, BorrowedAlignmentRead, RecordMappingQuality, SamWriteError,
    SamWritePhase, sam_flag, sam_header_bytes, sam_record_bytes, write_sam_header,
    write_sam_record,
};
use bsbit_index::reference::{ContigInput, ReferenceBuildLimits, ReferenceIndex};

fn normalized(raw: &[u8]) -> NormalizedSequence {
    normalize_dna(raw).expect("test input is normalized DNA")
}

fn reference(catalog: &[(&[u8], &[u8])]) -> ReferenceIndex {
    ReferenceIndex::build(
        catalog
            .iter()
            .map(|(name, sequence)| ContigInput::new(name.to_vec(), normalized(sequence)))
            .collect(),
        ReferenceBuildLimits::MAX,
    )
    .expect("bounded test reference builds")
}

fn single(reference: &ReferenceIndex, raw: &[u8], budget: u64) -> SingleFixture {
    single_fixture(reference, raw, budget)
}

fn exact_alignment(
    reference: &ReferenceIndex,
    read: &[u8],
    start: u64,
    strand: BisulfiteStrand,
) -> VerifiedAlignment {
    let query = normalized(read);
    let contig = reference
        .contig_by_ordinal(0)
        .expect("fixture contig exists");
    let end = start
        .checked_add(u64::try_from(query.bases().len()).expect("fixture read length fits u64"))
        .expect("fixture endpoint fits u64");
    let interval =
        ReferenceInterval::new(start, end, ReferenceLength::new(contig.sequence().len()))
            .expect("fixture placement is bounded");
    let contig_id = reference.contig_id(0).expect("fixture contig id exists");
    traceback_read_placement(reference, &query, &contig_id, interval, strand, 0)
        .expect("exact fixture placement materializes")
}

fn assert_single_ungapped_fast_path_matches_traceback(
    reference_raw: &[u8],
    read_raw: &[u8],
    strand: BisulfiteStrand,
    distance: u8,
) {
    let reference = reference(&[(b"chr1", reference_raw)]);
    let read = normalized(read_raw);
    let contig = reference
        .contig_by_ordinal(0)
        .expect("fixture contig exists");
    let interval =
        ReferenceInterval::new(0, read.len(), ReferenceLength::new(contig.sequence().len()))
            .expect("fixture interval is bounded");
    let placement = AlignmentPlacement::new(0, interval, strand, distance);
    let mapping_quality = RecordMappingQuality::Calibrated(37);
    let limits = AlignmentRecordLimits::default();
    let fast = try_build_indexed_single_ungapped_alignment_record(
        &reference,
        b"read",
        AlignmentRead::new(
            &read,
            Some(b"IIII".get(..read_raw.len()).expect("short quality")),
        ),
        placement,
        mapping_quality,
        limits,
    )
    .expect("fast-path evaluation succeeds")
    .expect("fixture is certified ungapped");
    let contig_id = reference.contig_id(0).expect("fixture contig id");
    let alignment =
        traceback_read_placement(&reference, &read, &contig_id, interval, strand, distance)
            .expect("full traceback succeeds");
    let expected = build_indexed_single_alignment_record(
        &reference,
        b"read",
        AlignmentRead::new(
            &read,
            Some(b"IIII".get(..read_raw.len()).expect("short quality")),
        ),
        Some(&alignment),
        mapping_quality,
        limits,
    )
    .expect("traceback record builds");
    assert_eq!(fast, expected);
}

#[test]
fn indexed_single_ungapped_fast_path_is_record_identical() {
    assert_single_ungapped_fast_path_matches_traceback(b"ACGT", b"ACGT", BisulfiteStrand::OT, 0);
    assert_single_ungapped_fast_path_matches_traceback(b"AAAA", b"AACA", BisulfiteStrand::OT, 1);
    assert_single_ungapped_fast_path_matches_traceback(b"ACGT", b"AGGA", BisulfiteStrand::OT, 2);
    assert_single_ungapped_fast_path_matches_traceback(b"AACG", b"CGTT", BisulfiteStrand::OB, 0);
}

#[test]
fn indexed_single_ungapped_fast_path_declines_shifted_gap_ties() {
    let reference = reference(&[(b"chr1", b"AC")]);
    let read = normalized(b"CA");
    let interval = ReferenceInterval::new(0, 2, ReferenceLength::new(2)).expect("interval");
    let observed = try_build_indexed_single_ungapped_alignment_record(
        &reference,
        b"tie",
        AlignmentRead::new(&read, Some(b"II")),
        AlignmentPlacement::new(0, interval, BisulfiteStrand::OT, 2),
        RecordMappingQuality::Calibrated(37),
        AlignmentRecordLimits::default(),
    )
    .expect("certificate evaluation succeeds");
    assert!(observed.is_none());
}

fn single_record(
    reference: &ReferenceIndex,
    name: &[u8],
    raw: &[u8],
    quality: Option<&[u8]>,
    budget: u64,
) -> AlignmentRecord {
    single_record_with_mode(
        reference,
        name,
        raw,
        quality,
        budget,
        AlignmentAuxiliaryMode::Minimal,
    )
}

fn single_record_with_mode(
    reference: &ReferenceIndex,
    name: &[u8],
    raw: &[u8],
    quality: Option<&[u8]>,
    budget: u64,
    auxiliary_mode: AlignmentAuxiliaryMode,
) -> AlignmentRecord {
    let fixture = single(reference, raw, budget);
    build_single_alignment_record_with_auxiliary_mode(
        reference,
        name,
        AlignmentRead::new(&fixture.query, quality),
        fixture.alignment.as_ref(),
        fixture.mapping_quality,
        AlignmentRecordLimits::default(),
        auxiliary_mode,
    )
    .expect("valid record builds")
}

fn fields(line: &[u8]) -> Vec<&[u8]> {
    assert_eq!(line.last(), Some(&b'\n'));
    line[..line.len() - 1]
        .split(|byte| *byte == b'\t')
        .collect()
}

fn number(raw: &[u8]) -> u64 {
    std::str::from_utf8(raw)
        .expect("numeric SAM field is UTF-8")
        .parse()
        .expect("numeric SAM field parses")
}

fn tag<'a>(decoded: &'a [&'a [u8]], prefix: &[u8]) -> &'a [u8] {
    decoded
        .iter()
        .copied()
        .find_map(|field| field.strip_prefix(prefix))
        .expect("required tag exists")
}

fn independent_literal_oracle(
    forward_reference: &[u8],
    start: usize,
    cigar: &[u8],
    stored_sequence: &[u8],
) -> (u64, Vec<u8>) {
    let mut reference_index = start;
    let mut query_index = 0_usize;
    let mut length = 0_usize;
    let mut saw_length_digit = false;
    let mut matches = 0_u64;
    let mut nm = 0_u64;
    let mut md = Vec::new();
    for &byte in cigar {
        if byte.is_ascii_digit() {
            length = length
                .checked_mul(10)
                .and_then(|value| value.checked_add(usize::from(byte - b'0')))
                .expect("tiny CIGAR length is representable");
            saw_length_digit = true;
            continue;
        }
        assert!(saw_length_digit && length > 0, "canonical CIGAR run");
        match byte {
            b'M' => {
                for _ in 0..length {
                    let reference_base = forward_reference[reference_index].to_ascii_uppercase();
                    let query_base = stored_sequence[query_index].to_ascii_uppercase();
                    if reference_base == query_base && reference_base != b'N' {
                        matches = matches.checked_add(1).expect("tiny match count");
                    } else {
                        md.extend_from_slice(matches.to_string().as_bytes());
                        md.push(reference_base);
                        matches = 0;
                        nm += 1;
                    }
                    reference_index += 1;
                    query_index += 1;
                }
            }
            b'I' => {
                query_index += length;
                nm += u64::try_from(length).expect("tiny insertion length");
            }
            b'D' => {
                md.extend_from_slice(matches.to_string().as_bytes());
                md.push(b'^');
                md.extend_from_slice(&forward_reference[reference_index..reference_index + length]);
                matches = 0;
                reference_index += length;
                nm += u64::try_from(length).expect("tiny deletion length");
            }
            _ => panic!("unsupported test CIGAR operation"),
        }
        length = 0;
        saw_length_digit = false;
    }
    assert_eq!(length, 0);
    assert_eq!(query_index, stored_sequence.len());
    md.extend_from_slice(matches.to_string().as_bytes());
    (nm, md)
}

#[test]
fn named_forward_reverse_conversion_and_unmapped_records_are_exact() {
    let exact_reference = reference(&[(b"chr", b"GGACCTAA")]);
    let exact = single_record(&exact_reference, b"exact", b"ACCT", Some(b"ABCD"), 0);
    assert_eq!(
        sam_record_bytes(&exact, AlignmentRecordLimits::default()).expect("SAM encodes"),
        b"exact\t0\tchr\t3\t255\t4M\t*\t0\t0\tACCT\tABCD\tNM:i:0\tXG:Z:CT\n"
    );
    assert_eq!(exact.mapping_quality(), RecordMappingQuality::Unavailable);

    let reverse_reference = reference(&[(b"chr", b"TTAACGAA")]);
    let input_quality = *b"ABCD";
    let reverse = single_record(
        &reverse_reference,
        b"reverse",
        b"CGTT",
        Some(&input_quality),
        0,
    );
    assert_eq!(reverse.sequence(), b"AACG");
    assert_eq!(reverse.quality(), Some(b"DCBA".as_slice()));
    assert_eq!(sam_flag(&reverse), 16);
    assert_eq!(input_quality, *b"ABCD", "input quality is not mutated");
    assert_eq!(
        sam_record_bytes(&reverse, AlignmentRecordLimits::default()).expect("SAM encodes"),
        b"reverse\t16\tchr\t3\t255\t4M\t*\t0\t0\tAACG\tDCBA\tNM:i:0\tXG:Z:GA\n"
    );

    let conversion_reference = reference(&[(b"chr", b"GGACCGAA")]);
    let conversion = single_record(&conversion_reference, b"ct", b"ATTG", None, 0);
    let mapping = conversion.mapping().expect("conversion maps");
    assert_eq!(mapping.literal_nm(), 2);
    assert_eq!(mapping.md(), None);
    assert_eq!(
        sam_record_bytes(&conversion, AlignmentRecordLimits::default()).expect("SAM encodes"),
        b"ct\t0\tchr\t3\t255\t4M\t*\t0\t0\tATTG\t*\tNM:i:2\tXG:Z:CT\n"
    );

    let bottom_reference = reference(&[(b"chr", b"AGGTCC")]);
    let bottom = single_record(&bottom_reference, b"ga", b"ATTT", None, 0);
    let mapping = bottom.mapping().expect("bottom conversion maps");
    assert_eq!(mapping.orientation(), AlignmentOrientation::Reverse);
    assert_eq!(mapping.literal_nm(), 2);
    assert_eq!(mapping.md(), None);
    assert_eq!(bottom.sequence(), b"AAAT");
    assert_eq!(
        sam_record_bytes(&bottom, AlignmentRecordLimits::default()).expect("SAM encodes"),
        b"ga\t16\tchr\t1\t255\t4M\t*\t0\t0\tAAAT\t*\tNM:i:2\tXG:Z:GA\n"
    );

    let unmapped_reference = reference(&[(b"chr", b"AAAA")]);
    let unmapped = single_record(&unmapped_reference, b"none", b"GGGG", Some(b"!!!!"), 0);
    assert!(!unmapped.is_mapped());
    assert_eq!(sam_flag(&unmapped), 4);
    assert_eq!(
        sam_record_bytes(&unmapped, AlignmentRecordLimits::default()).expect("SAM encodes"),
        b"none\t4\t*\t0\t0\t*\t*\t0\t0\tGGGG\t!!!!\n"
    );
}

#[test]
fn literal_nm_md_and_cigar_equal_an_independent_forward_reference_oracle() {
    struct LiteralCase {
        reference: &'static [u8],
        read: &'static [u8],
        budget: u64,
        cigar: &'static [u8],
        nm: u64,
        md: &'static [u8],
    }
    let cases = [
        LiteralCase {
            reference: b"GGACCGAA",
            read: b"ATTG",
            budget: 0,
            cigar: b"4M",
            nm: 2,
            md: b"1C0C1",
        },
        LiteralCase {
            reference: b"AGGTCC",
            read: b"ATTT",
            budget: 0,
            cigar: b"4M",
            nm: 2,
            md: b"1G0G1",
        },
        LiteralCase {
            reference: b"AAGC",
            read: b"AATC",
            budget: 1,
            cigar: b"4M",
            nm: 1,
            md: b"2G1",
        },
        LiteralCase {
            reference: b"AGTC",
            read: b"ACGTC",
            budget: 1,
            cigar: b"1M1I3M",
            nm: 1,
            md: b"4",
        },
        LiteralCase {
            reference: b"ACGAC",
            read: b"ACAC",
            budget: 1,
            cigar: b"2M1D2M",
            nm: 1,
            md: b"2^G2",
        },
        LiteralCase {
            reference: b"AANCG",
            read: b"AATCG",
            budget: 1,
            cigar: b"5M",
            nm: 1,
            md: b"2N2",
        },
    ];
    for case in cases {
        let reference = reference(&[(b"chr", case.reference)]);
        let record = single_record_with_mode(
            &reference,
            b"oracle",
            case.read,
            None,
            case.budget,
            AlignmentAuxiliaryMode::Bismark,
        );
        let line =
            sam_record_bytes(&record, AlignmentRecordLimits::default()).expect("SAM encodes");
        let decoded = fields(&line);
        assert_eq!(decoded.len(), 16);
        assert_eq!(decoded[5], case.cigar);
        let start = usize::try_from(number(decoded[3]) - 1).expect("tiny start");
        let (oracle_nm, oracle_md) =
            independent_literal_oracle(case.reference, start, decoded[5], decoded[9]);
        assert_eq!(oracle_nm, case.nm);
        assert_eq!(oracle_md, case.md);
        assert_eq!(number(tag(&decoded, b"NM:i:")), oracle_nm);
        assert_eq!(
            record.mapping().and_then(|mapping| mapping.md()),
            Some(oracle_md.as_slice())
        );
        assert_eq!(tag(&decoded, b"MD:Z:"), oracle_md);
        assert_eq!(
            tag(&decoded, b"XM:Z:"),
            record.mapping().unwrap().bismark_xm().unwrap()
        );
        assert_eq!(
            tag(&decoded, b"XR:Z:"),
            record.mapping().unwrap().bismark_xr()
        );
        assert_eq!(
            tag(&decoded, b"XG:Z:"),
            record.mapping().unwrap().bismark_xg()
        );
    }
}

#[test]
fn path_ambiguity_and_reference_boundaries_preserve_primary_coordinates() {
    let path_reference = reference(&[(b"chr", b"AA")]);
    let path = single_record(&path_reference, b"path", b"AAA", None, 1);
    let path_line = sam_record_bytes(&path, AlignmentRecordLimits::default()).expect("SAM");
    let decoded = fields(&path_line);
    assert_eq!(decoded[3], b"1");
    assert_eq!(decoded[4], b"255");
    assert_eq!(decoded[5], b"1I2M");
    assert_eq!(tag(&decoded, b"NM:i:"), b"1");
    assert_eq!(tag(&decoded, b"XG:Z:"), b"CT");

    for (reference_raw, expected_position) in
        [(b"AGTCGGGG".as_slice(), 1_u64), (b"GGGGAGTC".as_slice(), 5)]
    {
        let reference = reference(&[(b"chr", reference_raw)]);
        let record = single_record(&reference, b"edge", b"AGTC", None, 0);
        let line = sam_record_bytes(&record, AlignmentRecordLimits::default()).expect("SAM");
        assert_eq!(number(fields(&line)[3]), expected_position);
    }
}

#[test]
fn ambiguous_primary_is_canonical_and_truthfully_uses_zero_mapq() {
    let reference = reference(&[(b"chr", b"GGACCTAAACCTCC")]);
    let result = single(&reference, b"ACCT", 0);
    assert_eq!(result.mapping_quality, RecordMappingQuality::Tied);
    let first = build_single_alignment_record(
        &reference,
        b"tie",
        AlignmentRead::new(&result.query, None),
        result.alignment.as_ref(),
        result.mapping_quality,
        AlignmentRecordLimits::default(),
    )
    .expect("ambiguous record builds");
    let second = build_single_alignment_record(
        &reference,
        b"tie",
        AlignmentRead::new(&result.query, None),
        result.alignment.as_ref(),
        result.mapping_quality,
        AlignmentRecordLimits::default(),
    )
    .expect("ambiguous record rebuilds");
    assert_eq!(first, second);
    assert_eq!(first.mapping_quality(), RecordMappingQuality::Tied);
    assert_eq!(
        fields(&sam_record_bytes(&first, AlignmentRecordLimits::default()).expect("SAM"))[4],
        b"0"
    );
}

#[test]
fn direct_pair_preserves_full_reads_and_orients_3prime_soft_clips() {
    let reference = reference(&[(b"chr", b"AACCGTGATCTAGGCTTACGGAAT")]);
    let first_retained = normalized(b"CCGTGA");
    let second_retained = normalized(b"TCCGTA");
    let first_alignment = exact_alignment(&reference, b"CCGTGA", 2, BisulfiteStrand::OT);
    let second_alignment = exact_alignment(&reference, b"TCCGTA", 16, BisulfiteStrand::CTOT);
    let first_full = normalized(b"CCGTGAAAAA");
    let second_full = normalized(b"TCCGTACCCC");
    let mut batch = AlignmentRecordBatch::new();
    let mut composer = PairedRecordComposer::new();
    composer
        .push_soft_clipped_retained_unique_pair(
            &reference,
            b"clipped",
            BorrowedAlignmentRead::new(first_full.bases(), b"ABCDEFGHIJ"),
            BorrowedAlignmentRead::new(second_full.bases(), b"123456789:"),
            0..first_retained.bases().len(),
            0..second_retained.bases().len(),
            &first_retained,
            &second_retained,
            &first_alignment,
            &second_alignment,
            AlignmentRecordLimits::default(),
            AlignmentAuxiliaryMode::Bismark,
            20,
        )
        .expect("soft-clipped direct pair builds");
    composer
        .flush_into(&mut batch, AlignmentRecordLimits::default())
        .expect("soft-clipped pair flushes");
    let records = batch.records().collect::<Vec<_>>();
    assert_eq!(records.len(), 2);
    assert_eq!(records[0].mapping_quality(), 20);
    assert_eq!(records[1].mapping_quality(), 20);
    assert_eq!(records[0].sequence(), b"CCGTGAAAAA");
    assert_eq!(records[1].sequence(), b"GGGGTACGGA");
    assert_eq!(records[1].quality(), Some(b":987654321".as_slice()));
    assert_eq!(
        records[0]
            .cigar()
            .iter()
            .map(|run| (run.operation(), run.length()))
            .collect::<Vec<_>>(),
        vec![
            (AlignmentCigarOp::Match, 6),
            (AlignmentCigarOp::SoftClip, 4)
        ]
    );
    assert_eq!(
        records[1]
            .cigar()
            .iter()
            .map(|run| (run.operation(), run.length()))
            .collect::<Vec<_>>(),
        vec![
            (AlignmentCigarOp::SoftClip, 4),
            (AlignmentCigarOp::Match, 6)
        ]
    );
    assert_eq!(records[0].md(), Some(b"6".as_slice()));
    assert_eq!(records[1].md(), Some(b"6".as_slice()));
    assert_eq!(records[0].bismark_xm(), Some(b"XZ........".as_slice()));
    assert_eq!(records[1].bismark_xm(), Some(b"......Z...".as_slice()));
    assert_eq!(records[0].bismark_xr(), b"CT");
    assert_eq!(records[0].bismark_xg(), b"CT");
    assert_eq!(records[1].bismark_xr(), b"GA");
    assert_eq!(records[1].bismark_xg(), b"CT");
}

#[test]
fn unmapped_pair_preserves_primary_records_and_input_orientation() {
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

    let records = batch.records().collect::<Vec<_>>();
    assert_eq!(records.len(), 2);
    assert_eq!(records[0].query_name(), b"unmapped-pair");
    assert_eq!(records[1].query_name(), b"unmapped-pair");
    assert_eq!(records[0].flag(), 77);
    assert_eq!(records[1].flag(), 141);
    for record in &records {
        assert_eq!(record.reference_ordinal(), None);
        assert_eq!(record.position(), 0);
        assert_eq!(record.mapping_quality(), 0);
        assert!(record.cigar().is_empty());
        assert_eq!(record.mate_reference_ordinal(), None);
        assert_eq!(record.mate_position(), 0);
        assert_eq!(record.template_length(), 0);
        assert_eq!(record.md(), None);
    }
    assert_eq!(records[0].sequence(), b"ACGTN");
    assert_eq!(records[0].quality(), Some(b"ABCDE".as_slice()));
    assert_eq!(records[1].sequence(), b"TGCAN");
    assert_eq!(records[1].quality(), Some(b"12345".as_slice()));
}

#[test]
fn direct_pair_orients_bidirectional_soft_clips() {
    let reference = reference(&[(b"chr", b"AACCGTGATCTAGGCTTACGGAAT")]);
    let first_retained = normalized(b"CCGTGA");
    let second_retained = normalized(b"TCCGTA");
    let first_alignment = exact_alignment(&reference, b"CCGTGA", 2, BisulfiteStrand::OT);
    let second_alignment = exact_alignment(&reference, b"TCCGTA", 16, BisulfiteStrand::CTOT);
    let first_full = normalized(b"TCCGTGAAAA");
    let second_full = normalized(b"GGTCCGTACCCC");
    let mut batch = AlignmentRecordBatch::new();
    let mut composer = PairedRecordComposer::new();
    composer
        .push_soft_clipped_retained_unique_pair(
            &reference,
            b"two-ended",
            BorrowedAlignmentRead::new(first_full.bases(), b"ABCDEFGHIJ"),
            BorrowedAlignmentRead::new(second_full.bases(), b"123456789:;<"),
            1..7,
            2..8,
            &first_retained,
            &second_retained,
            &first_alignment,
            &second_alignment,
            AlignmentRecordLimits::default(),
            AlignmentAuxiliaryMode::Minimal,
            20,
        )
        .expect("bidirectionally clipped direct pair builds");
    composer
        .flush_into(&mut batch, AlignmentRecordLimits::default())
        .expect("bidirectionally clipped pair flushes");
    let records = batch.records().collect::<Vec<_>>();
    assert_eq!(records.len(), 2);
    assert_eq!(records[0].sequence(), b"TCCGTGAAAA");
    assert_eq!(records[1].sequence(), b"GGGGTACGGACC");
    assert_eq!(
        records[0]
            .cigar()
            .iter()
            .map(|run| (run.operation(), run.length()))
            .collect::<Vec<_>>(),
        vec![
            (AlignmentCigarOp::SoftClip, 1),
            (AlignmentCigarOp::Match, 6),
            (AlignmentCigarOp::SoftClip, 3),
        ]
    );
    assert_eq!(
        records[1]
            .cigar()
            .iter()
            .map(|run| (run.operation(), run.length()))
            .collect::<Vec<_>>(),
        vec![
            (AlignmentCigarOp::SoftClip, 4),
            (AlignmentCigarOp::Match, 6),
            (AlignmentCigarOp::SoftClip, 2),
        ]
    );
}

#[test]
fn direct_fast_path_orients_bidirectional_soft_clips() {
    let reference = reference(&[(b"chr", b"AACCGTGATCTAGGCTTACGGAAT")]);
    let first_alignment = exact_alignment(&reference, b"CCGTGA", 2, BisulfiteStrand::OT);
    let second_alignment = exact_alignment(&reference, b"TCCGTA", 16, BisulfiteStrand::CTOT);
    let first_full = normalized(b"TCCGTGAAAA");
    let second_full = normalized(b"GGTCCGTACCCC");
    let mut batch = AlignmentRecordBatch::new();
    let mut composer = PairedRecordComposer::new();
    assert!(
        composer
            .try_push_soft_clipped_ungapped_pair(
            &reference,
            b"two-ended-fast",
            BorrowedAlignmentRead::new(first_full.bases(), b"ABCDEFGHIJ"),
            BorrowedAlignmentRead::new(second_full.bases(), b"123456789:;<"),
            1..7,
            2..8,
            AlignmentPlacement::new(0, first_alignment.interval(), first_alignment.strand(), 0,),
            AlignmentPlacement::new(0, second_alignment.interval(), second_alignment.strand(), 0,),
            AlignmentRecordLimits::default(),
            20,
        )
        .expect("slab soft-clipped fast path succeeds")
    );
    composer
        .flush_into(&mut batch, AlignmentRecordLimits::default())
        .expect("fast-path pair flushes");
    let records = batch.records().collect::<Vec<_>>();
    assert_eq!(records.len(), 2);
    assert_eq!(records[0].sequence(), b"TCCGTGAAAA");
    assert_eq!(records[1].sequence(), b"GGGGTACGGACC");
    assert_eq!(records[0].mapping_quality(), 20);
    assert_eq!(records[1].mapping_quality(), 20);
    assert_eq!(records[0].literal_nm(), 0);
    assert_eq!(records[1].literal_nm(), 0);
    assert_eq!(records[0].md(), None);
    assert_eq!(records[1].md(), None);
    assert_eq!(
        records[0]
            .cigar()
            .iter()
            .map(|run| (run.operation(), run.length()))
            .collect::<Vec<_>>(),
        vec![
            (AlignmentCigarOp::SoftClip, 1),
            (AlignmentCigarOp::Match, 6),
            (AlignmentCigarOp::SoftClip, 3),
        ]
    );
    assert_eq!(
        records[1]
            .cigar()
            .iter()
            .map(|run| (run.operation(), run.length()))
            .collect::<Vec<_>>(),
        vec![
            (AlignmentCigarOp::SoftClip, 4),
            (AlignmentCigarOp::Match, 6),
            (AlignmentCigarOp::SoftClip, 2),
        ]
    );
}

#[test]
fn header_bytes_are_exact_ordered_and_bounded() {
    let reference = reference(&[(b"alpha", b"ACGT"), (b"beta", b"NN")]);
    let header =
        build_sam_header(&reference, AlignmentRecordLimits::default()).expect("header builds");
    assert_eq!(header.references()[0].name(), b"alpha");
    assert_eq!(header.references()[1].length(), 2);
    assert_eq!(
        sam_header_bytes(&header, AlignmentRecordLimits::default()).expect("header encodes"),
        concat!(
            "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:alpha\tLN:4\n@SQ\tSN:beta\tLN:2\n@PG\tID:bsbit\tPN:bsbit\tVN:",
            env!("CARGO_PKG_VERSION"),
            "\n"
        )
        .as_bytes()
    );

    let limits = AlignmentRecordLimits::new(254, 10, 10, 10, 100, 100, 100, 1_000, 1, 100, 1_000);
    assert!(matches!(
        build_sam_header(&reference, limits),
        Err(AlignmentRecordError::LimitExceeded {
            resource: AlignmentRecordResource::HeaderReferences,
            observed: 2,
            limit: 1,
        })
    ));
}

#[test]
#[allow(clippy::too_many_lines)]
fn lexical_owner_and_resource_failures_are_typed_and_fail_closed() {
    let primary_reference = reference(&[(b"chr", b"GGACCTAA")]);
    let read = normalized(b"ACCT");
    let result = single(&primary_reference, b"ACCT", 0);
    let default = AlignmentRecordLimits::default();
    assert!(matches!(
        build_single_alignment_record(
            &primary_reference,
            b"",
            AlignmentRead::new(&read, None),
            result.alignment.as_ref(),
            result.mapping_quality,
            default,
        ),
        Err(AlignmentRecordError::Format {
            source: HtsAlignmentRecordError::EmptyQueryName,
        })
    ));
    assert!(matches!(
        build_single_alignment_record(
            &primary_reference,
            b"bad name",
            AlignmentRead::new(&read, None),
            result.alignment.as_ref(),
            result.mapping_quality,
            default,
        ),
        Err(AlignmentRecordError::Format {
            source: HtsAlignmentRecordError::InvalidQueryNameByte {
                offset: 3,
                byte: b' '
            }
        })
    ));
    assert!(matches!(
        build_single_alignment_record(
            &primary_reference,
            b"q",
            AlignmentRead::new(&read, Some(b"!!!")),
            result.alignment.as_ref(),
            result.mapping_quality,
            default,
        ),
        Err(AlignmentRecordError::Format {
            source: HtsAlignmentRecordError::QualityLengthMismatch {
                sequence: 4,
                quality: 3,
            }
        })
    ));
    assert!(matches!(
        build_single_alignment_record(
            &primary_reference,
            b"q",
            AlignmentRead::new(&read, Some(b"!!\x7f!")),
            result.alignment.as_ref(),
            result.mapping_quality,
            default,
        ),
        Err(AlignmentRecordError::Format {
            source: HtsAlignmentRecordError::InvalidQualityByte {
                offset: 2,
                byte: 0x7f
            }
        })
    ));

    let tiny = AlignmentRecordLimits::new(254, 3, 4, 10, 100, 100, 100, 1_000, 10, 100, 1_000);
    assert!(matches!(
        build_single_alignment_record(
            &primary_reference,
            b"q",
            AlignmentRead::new(&read, None),
            result.alignment.as_ref(),
            result.mapping_quality,
            tiny,
        ),
        Err(AlignmentRecordError::LimitExceeded {
            resource: AlignmentRecordResource::ReadBases,
            observed: 4,
            limit: 3,
        })
    ));

    for (limits, expected_resource) in [
        (
            AlignmentRecordLimits::new(1, 4, 4, 10, 100, 100, 100, 1_000, 10, 100, 1_000),
            AlignmentRecordResource::QueryNameBytes,
        ),
        (
            AlignmentRecordLimits::new(254, 4, 4, 0, 100, 100, 100, 1_000, 10, 100, 1_000),
            AlignmentRecordResource::CigarRuns,
        ),
    ] {
        let error = build_single_alignment_record(
            &primary_reference,
            b"qq",
            AlignmentRead::new(&read, Some(b"!!!!")),
            result.alignment.as_ref(),
            result.mapping_quality,
            limits,
        )
        .expect_err("each exact resource cap fails closed");
        assert!(matches!(
            error,
            AlignmentRecordError::LimitExceeded { resource, .. }
                if resource == expected_resource
        ));
    }

    let encoded_resource_record = build_single_alignment_record(
        &primary_reference,
        b"qq",
        AlignmentRead::new(&read, Some(b"!!!!")),
        result.alignment.as_ref(),
        result.mapping_quality,
        default,
    )
    .expect("shared alignment record builds independently of text encoding caps");
    let short_cigar_text =
        AlignmentRecordLimits::new(254, 4, 4, 10, 1, 100, 100, 1_000, 10, 100, 1_000);
    assert!(matches!(
        sam_record_bytes(&encoded_resource_record, short_cigar_text),
        Err(HtsAlignmentRecordError::LimitExceeded {
            resource: AlignmentRecordResource::CigarTextBytes,
            ..
        })
    ));
    let no_optional_fields =
        AlignmentRecordLimits::new(254, 4, 4, 10, 100, 100, 0, 1_000, 10, 100, 1_000);
    assert!(matches!(
        build_single_alignment_record_with_auxiliary_mode(
            &primary_reference,
            b"qq",
            AlignmentRead::new(&read, Some(b"!!!!")),
            result.alignment.as_ref(),
            result.mapping_quality,
            no_optional_fields,
            AlignmentAuxiliaryMode::Bismark,
        ),
        Err(AlignmentRecordError::LimitExceeded {
            resource: AlignmentRecordResource::OptionalFieldBytes,
            ..
        })
    ));

    let no_md_capacity =
        AlignmentRecordLimits::new(254, 4, 4, 10, 100, 0, 100, 1_000, 10, 100, 1_000);
    let minimal = build_single_alignment_record(
        &primary_reference,
        b"qq",
        AlignmentRead::new(&read, Some(b"!!!!")),
        result.alignment.as_ref(),
        result.mapping_quality,
        no_md_capacity,
    )
    .expect("minimal mode allocates no MD bytes");
    assert_eq!(minimal.mapping().and_then(|mapping| mapping.md()), None);
    assert!(matches!(
        build_single_alignment_record_with_auxiliary_mode(
            &primary_reference,
            b"qq",
            AlignmentRead::new(&read, Some(b"!!!!")),
            result.alignment.as_ref(),
            result.mapping_quality,
            no_md_capacity,
            AlignmentAuxiliaryMode::Bismark,
        ),
        Err(AlignmentRecordError::LimitExceeded {
            resource: AlignmentRecordResource::MdBytes,
            ..
        })
    ));

    let other = reference(&[(b"chr", b"GGACCTAA")]);
    assert!(matches!(
        build_single_alignment_record(
            &other,
            b"q",
            AlignmentRead::new(&read, None),
            result.alignment.as_ref(),
            result.mapping_quality,
            default,
        ),
        Err(AlignmentRecordError::ReferenceAccess { .. })
    ));

    let invalid_name = reference(&[(b"*invalid", b"A")]);
    assert!(matches!(
        build_sam_header(&invalid_name, default),
        Err(AlignmentRecordError::Format {
            source: HtsAlignmentRecordError::InvalidReferenceNameByte {
                ordinal: 0,
                offset: 0,
                byte: Some(b'*'),
            }
        })
    ));

    let valid_punctuation = reference(&[(b"a:b|c", b"A")]);
    assert!(build_sam_header(&valid_punctuation, default).is_ok());
    for invalid in [
        b"bad\\name".as_slice(),
        b"bad,name".as_slice(),
        b"bad(name".as_slice(),
        b"bad[name".as_slice(),
        b"bad{name".as_slice(),
        b"bad<name".as_slice(),
    ] {
        let invalid_reference = reference(&[(invalid, b"A")]);
        assert!(matches!(
            build_sam_header(&invalid_reference, default),
            Err(AlignmentRecordError::Format {
                source: HtsAlignmentRecordError::InvalidReferenceNameByte { .. },
            })
        ));
    }

    let record = single_record(&primary_reference, b"q", b"ACCT", None, 0);
    let line_limit = AlignmentRecordLimits::new(254, 10, 10, 10, 100, 100, 100, 1, 10, 100, 100);
    assert!(matches!(
        sam_record_bytes(&record, line_limit),
        Err(HtsAlignmentRecordError::LimitExceeded {
            resource: AlignmentRecordResource::SamLineBytes,
            limit: 1,
            ..
        })
    ));
}

#[derive(Default)]
struct ChunkedWriter {
    bytes: Vec<u8>,
    chunk: usize,
}

impl Write for ChunkedWriter {
    fn write(&mut self, buffer: &[u8]) -> io::Result<usize> {
        let count = buffer.len().min(self.chunk.max(1));
        self.bytes.extend_from_slice(&buffer[..count]);
        Ok(count)
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

struct FailingWriter;

impl Write for FailingWriter {
    fn write(&mut self, _buffer: &[u8]) -> io::Result<usize> {
        Err(io::Error::other("injected writer failure"))
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

#[test]
fn generic_writers_preserve_exact_bytes_and_report_phase() {
    let reference = reference(&[(b"chr", b"GGACCTAA")]);
    let record = single_record(&reference, b"q", b"ACCT", None, 0);
    let header =
        build_sam_header(&reference, AlignmentRecordLimits::default()).expect("header builds");
    let expected_record =
        sam_record_bytes(&record, AlignmentRecordLimits::default()).expect("record encodes");
    let expected_header =
        sam_header_bytes(&header, AlignmentRecordLimits::default()).expect("header encodes");
    for chunk in 1..=7 {
        let mut writer = ChunkedWriter {
            bytes: Vec::new(),
            chunk,
        };
        write_sam_header(&mut writer, &header, AlignmentRecordLimits::default())
            .expect("chunked header writes");
        assert_eq!(writer.bytes, expected_header);
        writer.bytes.clear();
        write_sam_record(&mut writer, &record, AlignmentRecordLimits::default())
            .expect("chunked record writes");
        assert_eq!(writer.bytes, expected_record);
    }

    assert!(matches!(
        write_sam_header(
            &mut FailingWriter,
            &header,
            AlignmentRecordLimits::default(),
        ),
        Err(SamWriteError::Io {
            phase: SamWritePhase::Header,
            ..
        })
    ));
    assert!(matches!(
        write_sam_record(
            &mut FailingWriter,
            &record,
            AlignmentRecordLimits::default(),
        ),
        Err(SamWriteError::Io {
            phase: SamWritePhase::Record,
            ..
        })
    ));
}

#[test]
fn repeated_and_parallel_serialization_is_byte_identical() {
    let reference = reference(&[(b"chr", b"TTAACGAA")]);
    let record = Arc::new(single_record(
        &reference,
        b"parallel",
        b"CGTT",
        Some(b"ABCD"),
        0,
    ));
    let expected = sam_record_bytes(&record, AlignmentRecordLimits::default()).expect("SAM");
    for _ in 0..8 {
        assert_eq!(
            sam_record_bytes(&record, AlignmentRecordLimits::default()).expect("SAM"),
            expected
        );
    }
    let workers = (0..8)
        .map(|_| {
            let record = Arc::clone(&record);
            thread::spawn(move || {
                sam_record_bytes(&record, AlignmentRecordLimits::default()).expect("worker SAM")
            })
        })
        .collect::<Vec<_>>();
    for worker in workers {
        assert_eq!(worker.join().expect("worker does not panic"), expected);
    }
}

#[test]
fn reverse_orientation_is_present_in_the_named_fixture() {
    let reference = reference(&[(b"chr", b"TTAACGAA")]);
    let record = single_record(&reference, b"reverse", b"CGTT", None, 0);
    assert_eq!(
        record.mapping().expect("mapped").orientation(),
        AlignmentOrientation::Reverse
    );
}

#[test]
fn composition_arithmetic_and_decimal_helpers_are_exact() {
    assert!(matches!(
        checked_add_resource(u64::MAX, 1, AlignmentRecordResource::MdBytes),
        Err(AlignmentRecordError::ArithmeticOverflow {
            resource: AlignmentRecordResource::MdBytes,
            current: u64::MAX,
            increment: 1,
        })
    ));
    for (value, expected) in [
        (0, b"0".as_slice()),
        (9, b"9".as_slice()),
        (10, b"10".as_slice()),
        (u64::MAX, b"18446744073709551615".as_slice()),
    ] {
        let mut output = Vec::new();
        append_u64(&mut output, value);
        assert_eq!(output, expected);
        assert_eq!(storage_len(output.len()), decimal_digits(value));
    }
}

#[test]
fn bismark_methylation_calls_cover_context_case_and_boundaries() {
    for (reference, index, cytosine_strand, methylated, converted, upper, lower) in [
        (
            [Base::C, Base::G, Base::A],
            0,
            CytosineStrand::Top,
            Base::C,
            Base::T,
            b'Z',
            b'z',
        ),
        (
            [Base::C, Base::A, Base::G],
            0,
            CytosineStrand::Top,
            Base::C,
            Base::T,
            b'X',
            b'x',
        ),
        (
            [Base::C, Base::A, Base::A],
            0,
            CytosineStrand::Top,
            Base::C,
            Base::T,
            b'H',
            b'h',
        ),
        (
            [Base::A, Base::C, Base::G],
            2,
            CytosineStrand::Bottom,
            Base::G,
            Base::A,
            b'Z',
            b'z',
        ),
        (
            [Base::C, Base::A, Base::G],
            2,
            CytosineStrand::Bottom,
            Base::G,
            Base::A,
            b'X',
            b'x',
        ),
        (
            [Base::A, Base::A, Base::G],
            2,
            CytosineStrand::Bottom,
            Base::G,
            Base::A,
            b'H',
            b'h',
        ),
    ] {
        assert_eq!(
            bismark_methylation_call(
                &reference,
                index,
                reference[index],
                methylated,
                cytosine_strand,
            ),
            upper,
        );
        assert_eq!(
            bismark_methylation_call(
                &reference,
                index,
                reference[index],
                converted,
                cytosine_strand,
            ),
            lower,
        );
    }

    assert_eq!(
        bismark_methylation_call(&[Base::C], 0, Base::C, Base::C, CytosineStrand::Top),
        b'U',
    );
    assert_eq!(
        bismark_methylation_call(&[Base::G], 0, Base::G, Base::A, CytosineStrand::Bottom),
        b'u',
    );
    assert_eq!(
        bismark_methylation_call(
            &[Base::C, Base::N, Base::G],
            0,
            Base::C,
            Base::A,
            CytosineStrand::Top,
        ),
        b'.',
    );
}
