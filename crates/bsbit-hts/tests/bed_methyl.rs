//! Public extended-bedMethyl parse and encode contracts.

use bsbit_hts::{BedMethylContext, BedMethylError, BedMethylRecord, BedMethylStrand};

const ROW: &[u8] = b"chr1\t9\t10\tm,CHG,0\t10\t-\t9\t10\t255,0,0\t10\t70.00\t7\t3\t1\t2\t3\t4\t5";

#[test]
fn parsed_record_encodes_one_canonical_eighteen_column_row() {
    let parsed = BedMethylRecord::parse(ROW).expect("strict row");
    assert_eq!(parsed.contig(), b"chr1");
    assert_eq!(parsed.start(), 9);
    assert_eq!(parsed.end(), 10);
    assert_eq!(parsed.context(), BedMethylContext::Chg);
    assert_eq!(parsed.strand(), BedMethylStrand::Reverse);
    assert_eq!(parsed.coverage(), 10);
    assert_eq!(parsed.methylated(), 7);
    assert_eq!(parsed.unmethylated(), 3);
    assert_eq!(parsed.other_modification(), 1);
    assert_eq!(parsed.deleted(), 2);
    assert_eq!(parsed.failed(), 3);
    assert_eq!(parsed.different(), 4);
    assert_eq!(parsed.no_call(), 5);

    let mut encoded = Vec::new();
    parsed.encode(&mut encoded).expect("canonical row");
    assert_eq!(encoded, [ROW, b"\n"].concat());
    assert_eq!(
        BedMethylRecord::parse(encoded.strip_suffix(b"\n").expect("newline")).expect("reparse"),
        parsed
    );
}

#[test]
fn semantic_constructor_computes_span_coverage_and_rounded_percent() {
    let record = BedMethylRecord::new(
        b"chr2",
        3,
        BedMethylContext::Cg,
        BedMethylStrand::Forward,
        b"0,0,255",
        1,
        2,
        0,
        4,
        0,
        5,
        0,
    )
    .expect("record");
    let mut encoded = Vec::new();
    record.encode(&mut encoded).expect("row");
    assert_eq!(
        encoded,
        b"chr2\t3\t4\tm,CG,0\t3\t+\t3\t4\t0,0,255\t3\t33.33\t1\t2\t0\t4\t0\t5\t0\n"
    );
}

#[test]
fn parser_rejects_column_vocabulary_and_cross_column_inconsistency() {
    for invalid in [
        b"chr1\t9\t10".as_slice(),
        b"chr1\t9\t10\tm,CNN,0\t10\t-\t9\t10\t255,0,0\t10\t70.00\t7\t3\t1\t2\t3\t4\t5".as_slice(),
        b"chr1\t9\t11\tm,CHG,0\t10\t-\t9\t11\t255,0,0\t10\t70.00\t7\t3\t1\t2\t3\t4\t5".as_slice(),
        b"chr1\t9\t10\tm,CHG,0\t10\t.\t9\t10\t255,0,0\t10\t70.00\t7\t3\t1\t2\t3\t4\t5".as_slice(),
        b"chr1\t9\t10\tm,CHG,0\t11\t-\t9\t10\t255,0,0\t10\t70.00\t7\t3\t1\t2\t3\t4\t5".as_slice(),
        b"chr1\t9\t10\tm,CHG,0\t10\t-\t9\t10\t255,0,0\t10\t100.01\t7\t3\t1\t2\t3\t4\t5".as_slice(),
    ] {
        assert!(BedMethylRecord::parse(invalid).is_err(), "{invalid:?}");
    }

    assert!(matches!(
        BedMethylRecord::new(
            b"chr",
            0,
            BedMethylContext::Chh,
            BedMethylStrand::Forward,
            b"255,0,0",
            u64::MAX,
            1,
            0,
            0,
            0,
            0,
            0,
        ),
        Err(BedMethylError::CoverageOverflow)
    ));
}
