#![allow(missing_docs)]

use bsbit_align::score::EditDistance;

#[test]
fn edit_distance_is_checked() {
    assert_eq!(
        EditDistance::new(41).checked_add(1),
        Ok(EditDistance::new(42))
    );
    let overflow = EditDistance::new(u64::MAX)
        .checked_add(1)
        .expect_err("the next edit is not representable");
    assert_eq!(overflow.accumulated(), u64::MAX);
    assert_eq!(overflow.increment(), 1);
    assert_eq!(
        overflow.to_string(),
        "edit distance addition 18446744073709551615 + 1 overflowed"
    );
}
