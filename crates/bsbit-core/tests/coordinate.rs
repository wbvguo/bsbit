//! Exhaustive and boundary tests for checked coordinate-domain primitives.

use std::thread;

use bsbit_core::coordinate::{
    CoordinateDomain, CoordinateError, CoordinateOperation, CoordinateShift, OneBasedPosition,
    PositionConvention, QueryInterval, QueryLength, QueryPosition, ReferenceInterval,
    ReferenceLength, ReferencePosition,
};

fn assert_copy<T: Copy>() {}
fn assert_error<T: std::error::Error>() {}

#[test]
fn coordinate_values_have_copy_semantics_and_errors_are_standard_errors() {
    assert_copy::<ReferenceLength>();
    assert_copy::<QueryLength>();
    assert_copy::<ReferencePosition>();
    assert_copy::<QueryPosition>();
    assert_copy::<ReferenceInterval>();
    assert_copy::<QueryInterval>();
    assert_copy::<CoordinateShift>();
    assert_copy::<OneBasedPosition>();
    assert_error::<CoordinateError>();
}

#[test]
fn logical_display_forms_are_stable_and_show_half_open_bounds() {
    let reference_length = ReferenceLength::new(10);
    let query_length = QueryLength::new(7);
    let reference_position = ReferencePosition::new(2, reference_length).unwrap();
    let query_position = QueryPosition::new(4, query_length).unwrap();
    let one_based = reference_position.to_one_based(reference_length).unwrap();
    let reference_interval = ReferenceInterval::new(2, 5, reference_length).unwrap();
    let query_interval = QueryInterval::new(0, 7, query_length).unwrap();

    assert_eq!(reference_length.to_string(), "reference-length:10");
    assert_eq!(query_length.to_string(), "query-length:7");
    assert_eq!(reference_position.to_string(), "reference:2");
    assert_eq!(query_position.to_string(), "query:4");
    assert_eq!(one_based.to_string(), "reference-1based:3");
    assert_eq!(reference_interval.to_string(), "reference:[2,5)");
    assert_eq!(query_interval.to_string(), "query:[0,7)");
    assert_eq!(CoordinateShift::Zero.to_string(), "0");
    assert_eq!(CoordinateShift::forward(3).to_string(), "+3");
    assert_eq!(CoordinateShift::backward(3).to_string(), "-3");
}

#[test]
fn position_construction_distinguishes_existing_bases_from_boundaries() {
    for length in 0_u64..=12 {
        let reference_length = ReferenceLength::new(length);
        let query_length = QueryLength::new(length);

        for value in 0_u64..=13 {
            let reference = ReferencePosition::new(value, reference_length);
            let query = QueryPosition::new(value, query_length);
            assert_eq!(reference.is_ok(), value < length);
            assert_eq!(query.is_ok(), value < length);

            if value < length {
                assert_eq!(reference.unwrap().get(), value);
                assert_eq!(query.unwrap().get(), value);
            } else {
                assert_eq!(
                    reference.unwrap_err(),
                    CoordinateError::PositionOutOfBounds {
                        domain: CoordinateDomain::Reference,
                        convention: PositionConvention::ZeroBased,
                        value,
                        length,
                    }
                );
                assert_eq!(
                    query.unwrap_err(),
                    CoordinateError::PositionOutOfBounds {
                        domain: CoordinateDomain::Query,
                        convention: PositionConvention::ZeroBased,
                        value,
                        length,
                    }
                );
            }
        }

        // An interval boundary may equal the length even though a base
        // position with that value may not.
        assert_eq!(
            ReferenceInterval::new(length, length, reference_length)
                .unwrap()
                .len(),
            0
        );
        assert_eq!(
            QueryInterval::new(length, length, query_length)
                .unwrap()
                .len(),
            0
        );
    }
}

#[test]
fn reference_and_query_domains_remain_distinct_in_values_and_errors() {
    let reference = ReferencePosition::new(2, ReferenceLength::new(4)).unwrap();
    let query = QueryPosition::new(2, QueryLength::new(4)).unwrap();
    assert_eq!(reference.get(), query.get());
    assert_ne!(
        std::any::type_name::<ReferencePosition>(),
        std::any::type_name::<QueryPosition>()
    );
    assert_ne!(
        std::any::type_name::<ReferenceInterval>(),
        std::any::type_name::<QueryInterval>()
    );

    let reference_error = ReferenceInterval::new(0, 5, ReferenceLength::new(4)).unwrap_err();
    let query_error = QueryInterval::new(0, 5, QueryLength::new(4)).unwrap_err();
    assert!(matches!(
        reference_error,
        CoordinateError::OutOfBounds {
            domain: CoordinateDomain::Reference,
            ..
        }
    ));
    assert!(matches!(
        query_error,
        CoordinateError::OutOfBounds {
            domain: CoordinateDomain::Query,
            ..
        }
    ));
}

#[test]
fn interval_construction_is_exhaustive_for_small_lengths() {
    for length in 0_u64..=10 {
        for start in 0_u64..=12 {
            for end in 0_u64..=12 {
                let reference = ReferenceInterval::new(start, end, ReferenceLength::new(length));
                let query = QueryInterval::new(start, end, QueryLength::new(length));

                if start > end {
                    assert_eq!(
                        reference.unwrap_err(),
                        CoordinateError::InvertedInterval {
                            domain: CoordinateDomain::Reference,
                            start,
                            end,
                        }
                    );
                    assert_eq!(
                        query.unwrap_err(),
                        CoordinateError::InvertedInterval {
                            domain: CoordinateDomain::Query,
                            start,
                            end,
                        }
                    );
                } else if end > length {
                    assert_eq!(
                        reference.unwrap_err(),
                        CoordinateError::OutOfBounds {
                            domain: CoordinateDomain::Reference,
                            operation: CoordinateOperation::IntervalConstruction,
                            start,
                            end,
                            length,
                        }
                    );
                    assert_eq!(
                        query.unwrap_err(),
                        CoordinateError::OutOfBounds {
                            domain: CoordinateDomain::Query,
                            operation: CoordinateOperation::IntervalConstruction,
                            start,
                            end,
                            length,
                        }
                    );
                } else {
                    let reference = reference.unwrap();
                    let query = query.unwrap();
                    assert_eq!((reference.start(), reference.end()), (start, end));
                    assert_eq!((query.start(), query.end()), (start, end));
                    assert_eq!(reference.len(), end - start);
                    assert_eq!(query.len(), end - start);
                    assert_eq!(reference.is_empty(), start == end);
                    assert_eq!(query.is_empty(), start == end);
                }
            }
        }
    }
}

#[test]
fn ground_truth_reverse_transform_is_correct() {
    let interval = ReferenceInterval::new(2, 5, ReferenceLength::new(10)).unwrap();
    let reversed = interval.reverse(ReferenceLength::new(10)).unwrap();
    assert_eq!((reversed.start(), reversed.end()), (5, 8));
}

#[test]
fn reverse_transform_is_an_involution_for_all_small_intervals() {
    for length in 0_u64..=16 {
        for start in 0_u64..=length {
            for end in start..=length {
                let reference =
                    ReferenceInterval::new(start, end, ReferenceLength::new(length)).unwrap();
                let query = QueryInterval::new(start, end, QueryLength::new(length)).unwrap();

                let reference_reversed = reference.reverse(ReferenceLength::new(length)).unwrap();
                let query_reversed = query.reverse(QueryLength::new(length)).unwrap();
                assert_eq!(
                    (reference_reversed.start(), reference_reversed.end()),
                    (length - end, length - start)
                );
                assert_eq!(
                    (query_reversed.start(), query_reversed.end()),
                    (length - end, length - start)
                );
                assert_eq!(
                    reference_reversed
                        .reverse(ReferenceLength::new(length))
                        .unwrap(),
                    reference
                );
                assert_eq!(
                    query_reversed.reverse(QueryLength::new(length)).unwrap(),
                    query
                );
            }
        }
    }
}

#[test]
fn reverse_revalidates_the_interval_against_the_supplied_length() {
    let reference = ReferenceInterval::new(4, 8, ReferenceLength::new(8)).unwrap();
    assert_eq!(
        reference.reverse(ReferenceLength::new(7)).unwrap_err(),
        CoordinateError::OutOfBounds {
            domain: CoordinateDomain::Reference,
            operation: CoordinateOperation::ReverseTransform,
            start: 4,
            end: 8,
            length: 7,
        }
    );

    let query = QueryInterval::new(4, 8, QueryLength::new(8)).unwrap();
    assert_eq!(
        query.reverse(QueryLength::new(7)).unwrap_err(),
        CoordinateError::OutOfBounds {
            domain: CoordinateDomain::Query,
            operation: CoordinateOperation::ReverseTransform,
            start: 4,
            end: 8,
            length: 7,
        }
    );
}

#[test]
fn shift_constructors_have_only_one_zero_representation() {
    assert_eq!(CoordinateShift::forward(0), CoordinateShift::Zero);
    assert_eq!(CoordinateShift::backward(0), CoordinateShift::Zero);
    assert_eq!(CoordinateShift::Zero.magnitude(), 0);
    assert_eq!(CoordinateShift::forward(7).magnitude(), 7);
    assert_eq!(CoordinateShift::backward(9).magnitude(), 9);
}

fn assert_backward_translations(
    reference: ReferenceInterval,
    query: QueryInterval,
    start: u64,
    end: u64,
    length: u64,
    amount: u64,
) {
    let reference_result = reference.translate(
        CoordinateShift::backward(amount),
        ReferenceLength::new(length),
    );
    let query_result = query.translate(CoordinateShift::backward(amount), QueryLength::new(length));
    if amount <= start {
        let expected = (start - amount, end - amount);
        let reference_shifted = reference_result.unwrap();
        let query_shifted = query_result.unwrap();
        assert_eq!(
            (reference_shifted.start(), reference_shifted.end()),
            expected
        );
        assert_eq!((query_shifted.start(), query_shifted.end()), expected);
    } else {
        for (actual, domain) in [
            (reference_result.unwrap_err(), CoordinateDomain::Reference),
            (query_result.unwrap_err(), CoordinateDomain::Query),
        ] {
            assert_eq!(
                actual,
                CoordinateError::CoordinateUnderflow {
                    domain,
                    operation: CoordinateOperation::BackwardTranslation,
                    lhs: start,
                    rhs: amount,
                }
            );
        }
    }
}

fn assert_forward_translations(
    reference: ReferenceInterval,
    query: QueryInterval,
    start: u64,
    end: u64,
    length: u64,
    amount: u64,
) {
    let reference_result = reference.translate(
        CoordinateShift::forward(amount),
        ReferenceLength::new(length),
    );
    let query_result = query.translate(CoordinateShift::forward(amount), QueryLength::new(length));
    if amount <= length - end {
        let expected = (start + amount, end + amount);
        let reference_shifted = reference_result.unwrap();
        let query_shifted = query_result.unwrap();
        assert_eq!(
            (reference_shifted.start(), reference_shifted.end()),
            expected
        );
        assert_eq!((query_shifted.start(), query_shifted.end()), expected);
    } else {
        for (actual, domain) in [
            (reference_result.unwrap_err(), CoordinateDomain::Reference),
            (query_result.unwrap_err(), CoordinateDomain::Query),
        ] {
            assert_eq!(
                actual,
                CoordinateError::OutOfBounds {
                    domain,
                    operation: CoordinateOperation::ForwardTranslation,
                    start: start + amount,
                    end: end + amount,
                    length,
                }
            );
        }
    }
}

#[test]
fn translations_match_checked_arithmetic_exhaustively_on_small_inputs() {
    for length in 0_u64..=12 {
        for start in 0_u64..=length {
            for end in start..=length {
                let reference =
                    ReferenceInterval::new(start, end, ReferenceLength::new(length)).unwrap();
                let query = QueryInterval::new(start, end, QueryLength::new(length)).unwrap();

                assert_eq!(
                    reference
                        .translate(CoordinateShift::Zero, ReferenceLength::new(length))
                        .unwrap(),
                    reference
                );
                assert_eq!(
                    query
                        .translate(CoordinateShift::Zero, QueryLength::new(length))
                        .unwrap(),
                    query
                );

                for amount in 1_u64..=14 {
                    assert_backward_translations(reference, query, start, end, length, amount);
                    assert_forward_translations(reference, query, start, end, length, amount);
                }
            }
        }
    }
}

#[test]
fn translation_uses_documented_error_priority() {
    let interval = ReferenceInterval::new(8, 10, ReferenceLength::new(10)).unwrap();

    // Input revalidation wins over both a later underflow and a later overflow.
    for shift in [
        CoordinateShift::backward(u64::MAX),
        CoordinateShift::forward(u64::MAX),
    ] {
        assert_eq!(
            interval
                .translate(shift, ReferenceLength::new(9))
                .unwrap_err(),
            CoordinateError::OutOfBounds {
                domain: CoordinateDomain::Reference,
                operation: CoordinateOperation::TranslationInput,
                start: 8,
                end: 10,
                length: 9,
            }
        );
    }

    // For a valid input, representational overflow wins over enclosing bounds.
    let near_max =
        ReferenceInterval::new(u64::MAX - 2, u64::MAX, ReferenceLength::new(u64::MAX)).unwrap();
    assert_eq!(
        near_max
            .translate(CoordinateShift::forward(1), ReferenceLength::new(u64::MAX))
            .unwrap_err(),
        CoordinateError::CoordinateOverflow {
            domain: CoordinateDomain::Reference,
            operation: CoordinateOperation::ForwardTranslation,
            lhs: u64::MAX,
            rhs: 1,
        }
    );

    // If addition is representable, an enclosing-bound violation is distinct.
    let bounded = ReferenceInterval::new(3, 5, ReferenceLength::new(10)).unwrap();
    assert_eq!(
        bounded
            .translate(CoordinateShift::forward(6), ReferenceLength::new(10))
            .unwrap_err(),
        CoordinateError::OutOfBounds {
            domain: CoordinateDomain::Reference,
            operation: CoordinateOperation::ForwardTranslation,
            start: 9,
            end: 11,
            length: 10,
        }
    );
}

#[test]
fn one_based_reference_conversion_round_trips_all_small_positions() {
    for length in 0_u64..=32 {
        for one_based in 0_u64..=length.saturating_add(1) {
            let result = ReferencePosition::from_one_based(one_based, ReferenceLength::new(length));
            if one_based == 0 || one_based > length {
                assert_eq!(
                    result.unwrap_err(),
                    CoordinateError::PositionOutOfBounds {
                        domain: CoordinateDomain::Reference,
                        convention: PositionConvention::OneBased,
                        value: one_based,
                        length,
                    }
                );
            } else {
                let zero_based = result.unwrap();
                assert_eq!(zero_based.get(), one_based - 1);
                let external = zero_based
                    .to_one_based(ReferenceLength::new(length))
                    .unwrap();
                assert_eq!(external.get(), one_based);
                assert_eq!(
                    external
                        .to_zero_based(ReferenceLength::new(length))
                        .unwrap(),
                    zero_based
                );
            }
        }
    }
}

#[test]
fn one_based_conversion_handles_the_maximum_representable_length() {
    let length = ReferenceLength::new(u64::MAX);
    let last = ReferencePosition::new(u64::MAX - 1, length).unwrap();
    let one_based = last.to_one_based(length).unwrap();
    assert_eq!(one_based.get(), u64::MAX);
    assert_eq!(one_based.to_zero_based(length).unwrap(), last);
}

#[test]
fn intervals_handle_the_maximum_representable_length_without_wrapping() {
    let reference_length = ReferenceLength::new(u64::MAX);
    let query_length = QueryLength::new(u64::MAX);
    let reference = ReferenceInterval::new(u64::MAX - 2, u64::MAX, reference_length).unwrap();
    let query = QueryInterval::new(u64::MAX - 2, u64::MAX, query_length).unwrap();

    assert_eq!(
        reference.translate(CoordinateShift::Zero, reference_length),
        Ok(reference)
    );
    assert_eq!(
        query.translate(CoordinateShift::Zero, query_length),
        Ok(query)
    );
    assert_eq!(
        reference.reverse(reference_length).unwrap(),
        ReferenceInterval::new(0, 2, reference_length).unwrap()
    );
    assert_eq!(
        query.reverse(query_length).unwrap(),
        QueryInterval::new(0, 2, query_length).unwrap()
    );
}

#[test]
fn coordinate_operations_are_concurrently_deterministic() {
    let mut workers = Vec::new();
    for _ in 0..8 {
        workers.push(thread::spawn(|| {
            let length = ReferenceLength::new(10);
            let interval = ReferenceInterval::new(2, 5, length).unwrap();
            for _ in 0..1_000 {
                assert_eq!(
                    interval.reverse(length).unwrap(),
                    ReferenceInterval::new(5, 8, length).unwrap()
                );
                assert_eq!(
                    interval
                        .translate(CoordinateShift::forward(3), length)
                        .unwrap(),
                    ReferenceInterval::new(5, 8, length).unwrap()
                );
            }
        }));
    }
    for worker in workers {
        worker.join().expect("coordinate worker did not panic");
    }
}

#[test]
fn to_one_based_revalidates_against_the_supplied_contig_length() {
    let position = ReferencePosition::new(4, ReferenceLength::new(5)).unwrap();
    assert_eq!(
        position.to_one_based(ReferenceLength::new(4)).unwrap_err(),
        CoordinateError::PositionOutOfBounds {
            domain: CoordinateDomain::Reference,
            convention: PositionConvention::ZeroBased,
            value: 4,
            length: 4,
        }
    );
}
