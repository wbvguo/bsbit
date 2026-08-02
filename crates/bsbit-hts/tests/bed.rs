//! Public BED3+ interval-line syntax contracts.

use bsbit_hts::{BedError, BedInterval};

#[test]
fn bed3_and_bed3_plus_lines_borrow_the_first_three_fields() {
    let line = b"chr2\t17\t23\tfeature\t99";
    let interval = BedInterval::parse_line(line)
        .expect("BED3+ syntax")
        .expect("data line");
    assert_eq!(interval.contig(), b"chr2");
    assert_eq!(interval.start(), 17);
    assert_eq!(interval.end(), 23);
    assert_eq!(interval.contig().as_ptr(), line.as_ptr());

    let zero_length = BedInterval::parse_line(b"chr2\t23\t23")
        .expect("span policy is not codec syntax")
        .expect("data line");
    assert_eq!((zero_length.start(), zero_length.end()), (23, 23));

    let non_utf8 = BedInterval::parse_line(b"chr\xff\t1\t2")
        .expect("text policy belongs to the consumer")
        .expect("data line");
    assert_eq!(non_utf8.contig(), b"chr\xff");
}

#[test]
fn blank_comment_track_and_browser_lines_are_ignored() {
    for line in [
        b"".as_slice(),
        b"# comment".as_slice(),
        b"track name=targets".as_slice(),
        b"track\tname=targets".as_slice(),
        b"browser position chr1:1-10".as_slice(),
        b"browser\tposition chr1:1-10".as_slice(),
    ] {
        assert_eq!(BedInterval::parse_line(line), Ok(None), "{line:?}");
    }
}

#[test]
fn missing_invalid_and_overflowing_coordinates_are_typed() {
    assert_eq!(
        BedInterval::parse_line(b"chr1"),
        Err(BedError::ColumnCount { observed: 1 })
    );
    assert_eq!(
        BedInterval::parse_line(b"chr1\t1"),
        Err(BedError::ColumnCount { observed: 2 })
    );
    for (line, column) in [
        (b"chr1\t\t2".as_slice(), 2),
        (b"chr1\t-1\t2".as_slice(), 2),
        (b"chr1\t1\t2x".as_slice(), 3),
    ] {
        assert_eq!(
            BedInterval::parse_line(line),
            Err(BedError::InvalidInteger { column })
        );
    }
    assert_eq!(
        BedInterval::parse_line(b"chr1\t18446744073709551616\t20"),
        Err(BedError::IntegerOverflow { column: 2 })
    );
    assert_eq!(
        BedInterval::parse_line(b"chr1\t1\t18446744073709551616"),
        Err(BedError::IntegerOverflow { column: 3 })
    );
}
