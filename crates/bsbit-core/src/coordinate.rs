//! Checked coordinate-domain primitives.
//!
//! Reference and query coordinates intentionally use different types.  Base
//! positions name existing bases, whereas interval endpoints are half-open
//! boundaries and may equal the enclosing sequence length.  Every arithmetic
//! operation is checked and preserves the coordinate domain.
//!
//! Reference and query positions cannot be interchanged accidentally:
//!
//! ```compile_fail
//! use bsbit_core::coordinate::{QueryLength, QueryPosition, ReferencePosition};
//!
//! fn needs_reference(_: ReferencePosition) {}
//!
//! let query = QueryPosition::new(0, QueryLength::new(1)).unwrap();
//! needs_reference(query);
//! ```
//!
//! The same separation applies to intervals:
//!
//! ```compile_fail
//! use bsbit_core::coordinate::{QueryInterval, QueryLength, ReferenceInterval};
//!
//! fn needs_reference(_: ReferenceInterval) {}
//!
//! let query = QueryInterval::new(0, 1, QueryLength::new(1)).unwrap();
//! needs_reference(query);
//! ```

use core::fmt;
use core::num::NonZeroU64;

/// The coordinate space in which an error occurred.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum CoordinateDomain {
    /// A position or interval on one forward-reference contig.
    Reference,
    /// A position or interval on one oriented query.
    Query,
}

/// The external convention used by an existing-base position.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum PositionConvention {
    /// An internal, zero-based position.
    ZeroBased,
    /// An external, one-based reference position.
    OneBased,
}

/// The operation being performed when checked coordinate arithmetic failed.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum CoordinateOperation {
    /// Construction of a validated half-open interval.
    IntervalConstruction,
    /// Validation of the input to a translation.
    TranslationInput,
    /// Addition for a forward translation.
    ForwardTranslation,
    /// Subtraction for a backward translation.
    BackwardTranslation,
    /// Validation or subtraction for a reverse-coordinate transform.
    ReverseTransform,
    /// Conversion of a zero-based reference position to one-based form.
    ZeroToOneBased,
    /// Conversion of a one-based reference position to zero-based form.
    OneToZeroBased,
}

/// A structured failure from coordinate validation or arithmetic.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum CoordinateError {
    /// A value does not name an existing base under the stated convention.
    PositionOutOfBounds {
        /// Coordinate domain of the supplied value.
        domain: CoordinateDomain,
        /// Whether the supplied value is zero- or one-based.
        convention: PositionConvention,
        /// Offending position value.
        value: u64,
        /// Length of the enclosing sequence or contig.
        length: u64,
    },
    /// A half-open interval has `start > end`.
    InvertedInterval {
        /// Coordinate domain of the interval.
        domain: CoordinateDomain,
        /// Supplied start boundary.
        start: u64,
        /// Supplied end boundary.
        end: u64,
    },
    /// A half-open interval extends beyond its enclosing sequence.
    OutOfBounds {
        /// Coordinate domain of the interval.
        domain: CoordinateDomain,
        /// Operation that exposed the invalid bounds.
        operation: CoordinateOperation,
        /// Supplied or attempted start boundary.
        start: u64,
        /// Supplied or attempted end boundary.
        end: u64,
        /// Length of the enclosing sequence.
        length: u64,
    },
    /// Unsigned coordinate addition could not be represented.
    CoordinateOverflow {
        /// Coordinate domain of the operands.
        domain: CoordinateDomain,
        /// Operation whose addition overflowed.
        operation: CoordinateOperation,
        /// Left operand.
        lhs: u64,
        /// Right operand.
        rhs: u64,
    },
    /// Unsigned coordinate subtraction would produce a negative value.
    CoordinateUnderflow {
        /// Coordinate domain of the operands.
        domain: CoordinateDomain,
        /// Operation whose subtraction underflowed.
        operation: CoordinateOperation,
        /// Left operand.
        lhs: u64,
        /// Right operand.
        rhs: u64,
    },
}

impl fmt::Display for CoordinateError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match *self {
            Self::PositionOutOfBounds {
                domain,
                convention,
                value,
                length,
            } => write!(
                formatter,
                "{convention:?} {domain:?} position {value} is outside length {length}"
            ),
            Self::InvertedInterval { domain, start, end } => write!(
                formatter,
                "{domain:?} half-open interval [{start}, {end}) is inverted"
            ),
            Self::OutOfBounds {
                domain,
                operation,
                start,
                end,
                length,
            } => write!(
                formatter,
                "{domain:?} interval [{start}, {end}) is outside length {length} during {operation:?}"
            ),
            Self::CoordinateOverflow {
                domain,
                operation,
                lhs,
                rhs,
            } => write!(
                formatter,
                "{domain:?} coordinate addition {lhs} + {rhs} overflowed during {operation:?}"
            ),
            Self::CoordinateUnderflow {
                domain,
                operation,
                lhs,
                rhs,
            } => write!(
                formatter,
                "{domain:?} coordinate subtraction {lhs} - {rhs} underflowed during {operation:?}"
            ),
        }
    }
}

impl std::error::Error for CoordinateError {}

/// Length of one forward-reference contig.
#[derive(Clone, Copy, Debug, Default, Eq, Ord, PartialEq, PartialOrd)]
pub struct ReferenceLength(u64);

impl ReferenceLength {
    /// Constructs a reference length.  Zero is a valid mathematical length.
    #[must_use]
    pub const fn new(value: u64) -> Self {
        Self(value)
    }

    /// Returns the length as an unsigned base count.
    #[must_use]
    pub const fn get(self) -> u64 {
        self.0
    }
}

impl fmt::Display for ReferenceLength {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "reference-length:{}", self.0)
    }
}

/// Length of one normalized query.
#[derive(Clone, Copy, Debug, Default, Eq, Ord, PartialEq, PartialOrd)]
pub struct QueryLength(u64);

impl QueryLength {
    /// Constructs a query length.  Zero is a valid mathematical length.
    #[must_use]
    pub const fn new(value: u64) -> Self {
        Self(value)
    }

    /// Returns the length as an unsigned base count.
    #[must_use]
    pub const fn get(self) -> u64 {
        self.0
    }
}

impl fmt::Display for QueryLength {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "query-length:{}", self.0)
    }
}

/// A zero-based position of an existing base on one reference contig.
#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub struct ReferencePosition(u64);

impl ReferencePosition {
    /// Validates a zero-based reference-base position against `length`.
    ///
    /// No position is valid when `length` is zero.  This constructor is
    /// allocation-free and runs in `O(1)` time.
    ///
    /// # Errors
    ///
    /// Returns [`CoordinateError::PositionOutOfBounds`] when `value` is not
    /// smaller than `length`.
    pub const fn new(value: u64, length: ReferenceLength) -> Result<Self, CoordinateError> {
        if value < length.0 {
            Ok(Self(value))
        } else {
            Err(position_out_of_bounds(
                CoordinateDomain::Reference,
                PositionConvention::ZeroBased,
                value,
                length.0,
            ))
        }
    }

    /// Returns the zero-based position value.
    #[must_use]
    pub const fn get(self) -> u64 {
        self.0
    }

    /// Converts this reference position to one-based form after revalidating
    /// it against `length`.
    ///
    /// Requiring the contig length prevents a position validated for one
    /// contig from being silently reused at a shorter boundary.  This is an
    /// `O(1)`, allocation-free operation.
    ///
    /// # Errors
    ///
    /// Returns [`CoordinateError::PositionOutOfBounds`] if this position is
    /// not valid for `length`, or [`CoordinateError::CoordinateOverflow`] if
    /// the checked addition cannot be represented.
    pub const fn to_one_based(
        self,
        length: ReferenceLength,
    ) -> Result<OneBasedPosition, CoordinateError> {
        if self.0 >= length.0 {
            return Err(position_out_of_bounds(
                CoordinateDomain::Reference,
                PositionConvention::ZeroBased,
                self.0,
                length.0,
            ));
        }

        match self.0.checked_add(1) {
            Some(value) => Ok(OneBasedPosition(value)),
            None => Err(CoordinateError::CoordinateOverflow {
                domain: CoordinateDomain::Reference,
                operation: CoordinateOperation::ZeroToOneBased,
                lhs: self.0,
                rhs: 1,
            }),
        }
    }

    /// Converts a raw one-based reference position to zero-based form.
    ///
    /// Both zero and values greater than `length` return
    /// [`CoordinateError::PositionOutOfBounds`].  The operation is `O(1)` and
    /// allocation-free.
    ///
    /// # Errors
    ///
    /// Returns [`CoordinateError::PositionOutOfBounds`] when `value` is zero
    /// or greater than `length`.
    pub const fn from_one_based(
        value: u64,
        length: ReferenceLength,
    ) -> Result<Self, CoordinateError> {
        if value == 0 || value > length.0 {
            return Err(position_out_of_bounds(
                CoordinateDomain::Reference,
                PositionConvention::OneBased,
                value,
                length.0,
            ));
        }

        match value.checked_sub(1) {
            Some(position) => Ok(Self(position)),
            None => Err(CoordinateError::CoordinateUnderflow {
                domain: CoordinateDomain::Reference,
                operation: CoordinateOperation::OneToZeroBased,
                lhs: value,
                rhs: 1,
            }),
        }
    }
}

impl fmt::Display for ReferencePosition {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "reference:{}", self.0)
    }
}

/// A positive one-based reference-base position at an external boundary.
///
/// Values of this type are created only through checked reference-position
/// conversion.  Query positions deliberately have no corresponding API.
#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub struct OneBasedPosition(u64);

impl OneBasedPosition {
    /// Returns the positive one-based value.
    #[must_use]
    pub const fn get(self) -> u64 {
        self.0
    }

    /// Converts this external position to a zero-based reference position,
    /// explicitly validating it against `length`.
    ///
    /// # Errors
    ///
    /// Returns [`CoordinateError::PositionOutOfBounds`] if this one-based
    /// value is not within `1..=length`.
    pub const fn to_zero_based(
        self,
        length: ReferenceLength,
    ) -> Result<ReferencePosition, CoordinateError> {
        ReferencePosition::from_one_based(self.0, length)
    }
}

impl fmt::Display for OneBasedPosition {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "reference-1based:{}", self.0)
    }
}

/// A zero-based position of an existing base in one query.
#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub struct QueryPosition(u64);

impl QueryPosition {
    /// Validates a zero-based query-base position against `length`.
    ///
    /// No position is valid when `length` is zero.  This constructor is
    /// allocation-free and runs in `O(1)` time.
    ///
    /// # Errors
    ///
    /// Returns [`CoordinateError::PositionOutOfBounds`] when `value` is not
    /// smaller than `length`.
    pub const fn new(value: u64, length: QueryLength) -> Result<Self, CoordinateError> {
        if value < length.0 {
            Ok(Self(value))
        } else {
            Err(position_out_of_bounds(
                CoordinateDomain::Query,
                PositionConvention::ZeroBased,
                value,
                length.0,
            ))
        }
    }

    /// Returns the zero-based position value.
    #[must_use]
    pub const fn get(self) -> u64 {
        self.0
    }
}

impl fmt::Display for QueryPosition {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "query:{}", self.0)
    }
}

/// A normalized displacement that has exactly one representation of zero.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum CoordinateShift {
    /// No displacement.
    Zero,
    /// Move toward larger coordinates by a positive magnitude.
    Forward(NonZeroU64),
    /// Move toward smaller coordinates by a positive magnitude.
    Backward(NonZeroU64),
}

impl CoordinateShift {
    /// Constructs a forward displacement, normalizing magnitude zero to
    /// [`CoordinateShift::Zero`].
    #[must_use]
    pub const fn forward(magnitude: u64) -> Self {
        match NonZeroU64::new(magnitude) {
            Some(value) => Self::Forward(value),
            None => Self::Zero,
        }
    }

    /// Constructs a backward displacement, normalizing magnitude zero to
    /// [`CoordinateShift::Zero`] so negative zero cannot be represented.
    #[must_use]
    pub const fn backward(magnitude: u64) -> Self {
        match NonZeroU64::new(magnitude) {
            Some(value) => Self::Backward(value),
            None => Self::Zero,
        }
    }

    /// Returns the unsigned magnitude of this displacement.
    #[must_use]
    pub const fn magnitude(self) -> u64 {
        match self {
            Self::Zero => 0,
            Self::Forward(value) | Self::Backward(value) => value.get(),
        }
    }
}

impl fmt::Display for CoordinateShift {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Zero => formatter.write_str("0"),
            Self::Forward(value) => write!(formatter, "+{value}"),
            Self::Backward(value) => write!(formatter, "-{value}"),
        }
    }
}

/// A validated zero-based half-open interval on one reference contig.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ReferenceInterval {
    start: u64,
    end: u64,
}

impl ReferenceInterval {
    /// Constructs `[start, end)` after validation against `length`.
    ///
    /// Empty intervals and the end boundary `length` are valid.  An inverted
    /// interval is reported before any out-of-bounds condition.  Construction
    /// is `O(1)` and allocation-free.
    ///
    /// # Errors
    ///
    /// Returns [`CoordinateError::InvertedInterval`] when `start > end`, or
    /// [`CoordinateError::OutOfBounds`] when `end > length`, in that order.
    pub const fn new(
        start: u64,
        end: u64,
        length: ReferenceLength,
    ) -> Result<Self, CoordinateError> {
        match validate_interval(
            CoordinateDomain::Reference,
            CoordinateOperation::IntervalConstruction,
            start,
            end,
            length.0,
        ) {
            Ok(()) => Ok(Self { start, end }),
            Err(error) => Err(error),
        }
    }

    /// Returns the inclusive start boundary.
    #[must_use]
    pub const fn start(self) -> u64 {
        self.start
    }

    /// Returns the exclusive end boundary.
    #[must_use]
    pub const fn end(self) -> u64 {
        self.end
    }

    /// Returns `end - start`, which is representable by construction.
    #[must_use]
    pub const fn len(self) -> u64 {
        self.end - self.start
    }

    /// Reports whether this interval contains no bases.
    #[must_use]
    pub const fn is_empty(self) -> bool {
        self.start == self.end
    }

    /// Translates this interval while keeping it within `length`.
    ///
    /// The input is first revalidated against the supplied length.  Forward
    /// translation checks representational overflow before the enclosing
    /// bound; backward translation reports underflow when `start` is smaller
    /// than the magnitude.  This is an `O(1)`, allocation-free operation.
    ///
    /// # Errors
    ///
    /// Returns [`CoordinateError::OutOfBounds`] if the input or translated
    /// interval exceeds `length`, [`CoordinateError::CoordinateUnderflow`] for
    /// a backward shift past zero, or [`CoordinateError::CoordinateOverflow`]
    /// if forward arithmetic is not representable.
    pub const fn translate(
        self,
        shift: CoordinateShift,
        length: ReferenceLength,
    ) -> Result<Self, CoordinateError> {
        match translate_bounds(
            CoordinateDomain::Reference,
            self.start,
            self.end,
            shift,
            length.0,
        ) {
            Ok((start, end)) => Ok(Self { start, end }),
            Err(error) => Err(error),
        }
    }

    /// Maps `[start, end)` to `[length - end, length - start)`.
    ///
    /// The input is revalidated against `length`.  The transform is
    /// allocation-free, runs in `O(1)`, and applying it twice with the same
    /// length returns the original interval.
    ///
    /// # Errors
    ///
    /// Returns [`CoordinateError::OutOfBounds`] if this interval is not valid
    /// for `length`.  Checked subtraction failures are represented by
    /// [`CoordinateError::CoordinateUnderflow`].
    pub const fn reverse(self, length: ReferenceLength) -> Result<Self, CoordinateError> {
        match reverse_bounds(CoordinateDomain::Reference, self.start, self.end, length.0) {
            Ok((start, end)) => Ok(Self { start, end }),
            Err(error) => Err(error),
        }
    }
}

impl fmt::Display for ReferenceInterval {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "reference:[{},{})", self.start, self.end)
    }
}

/// A validated zero-based half-open interval in one query.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct QueryInterval {
    start: u64,
    end: u64,
}

impl QueryInterval {
    /// Constructs `[start, end)` after validation against `length`.
    ///
    /// Empty intervals and the end boundary `length` are valid.  An inverted
    /// interval is reported before any out-of-bounds condition.  Construction
    /// is `O(1)` and allocation-free.
    ///
    /// # Errors
    ///
    /// Returns [`CoordinateError::InvertedInterval`] when `start > end`, or
    /// [`CoordinateError::OutOfBounds`] when `end > length`, in that order.
    pub const fn new(start: u64, end: u64, length: QueryLength) -> Result<Self, CoordinateError> {
        match validate_interval(
            CoordinateDomain::Query,
            CoordinateOperation::IntervalConstruction,
            start,
            end,
            length.0,
        ) {
            Ok(()) => Ok(Self { start, end }),
            Err(error) => Err(error),
        }
    }

    /// Returns the inclusive start boundary.
    #[must_use]
    pub const fn start(self) -> u64 {
        self.start
    }

    /// Returns the exclusive end boundary.
    #[must_use]
    pub const fn end(self) -> u64 {
        self.end
    }

    /// Returns `end - start`, which is representable by construction.
    #[must_use]
    pub const fn len(self) -> u64 {
        self.end - self.start
    }

    /// Reports whether this interval contains no bases.
    #[must_use]
    pub const fn is_empty(self) -> bool {
        self.start == self.end
    }

    /// Translates this interval while keeping it within `length`.
    ///
    /// The input is first revalidated against the supplied length.  Forward
    /// translation checks representational overflow before the enclosing
    /// bound; backward translation reports underflow when `start` is smaller
    /// than the magnitude.  This is an `O(1)`, allocation-free operation.
    ///
    /// # Errors
    ///
    /// Returns [`CoordinateError::OutOfBounds`] if the input or translated
    /// interval exceeds `length`, [`CoordinateError::CoordinateUnderflow`] for
    /// a backward shift past zero, or [`CoordinateError::CoordinateOverflow`]
    /// if forward arithmetic is not representable.
    pub const fn translate(
        self,
        shift: CoordinateShift,
        length: QueryLength,
    ) -> Result<Self, CoordinateError> {
        match translate_bounds(
            CoordinateDomain::Query,
            self.start,
            self.end,
            shift,
            length.0,
        ) {
            Ok((start, end)) => Ok(Self { start, end }),
            Err(error) => Err(error),
        }
    }

    /// Maps `[start, end)` to `[length - end, length - start)`.
    ///
    /// The input is revalidated against `length`.  The transform is
    /// allocation-free, runs in `O(1)`, and applying it twice with the same
    /// length returns the original interval.
    ///
    /// # Errors
    ///
    /// Returns [`CoordinateError::OutOfBounds`] if this interval is not valid
    /// for `length`.  Checked subtraction failures are represented by
    /// [`CoordinateError::CoordinateUnderflow`].
    pub const fn reverse(self, length: QueryLength) -> Result<Self, CoordinateError> {
        match reverse_bounds(CoordinateDomain::Query, self.start, self.end, length.0) {
            Ok((start, end)) => Ok(Self { start, end }),
            Err(error) => Err(error),
        }
    }
}

impl fmt::Display for QueryInterval {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "query:[{},{})", self.start, self.end)
    }
}

const fn position_out_of_bounds(
    domain: CoordinateDomain,
    convention: PositionConvention,
    value: u64,
    length: u64,
) -> CoordinateError {
    CoordinateError::PositionOutOfBounds {
        domain,
        convention,
        value,
        length,
    }
}

const fn validate_interval(
    domain: CoordinateDomain,
    operation: CoordinateOperation,
    start: u64,
    end: u64,
    length: u64,
) -> Result<(), CoordinateError> {
    if start > end {
        return Err(CoordinateError::InvertedInterval { domain, start, end });
    }
    if end > length {
        return Err(CoordinateError::OutOfBounds {
            domain,
            operation,
            start,
            end,
            length,
        });
    }
    Ok(())
}

const fn translate_bounds(
    domain: CoordinateDomain,
    start: u64,
    end: u64,
    shift: CoordinateShift,
    length: u64,
) -> Result<(u64, u64), CoordinateError> {
    match validate_interval(
        domain,
        CoordinateOperation::TranslationInput,
        start,
        end,
        length,
    ) {
        Ok(()) => {}
        Err(error) => return Err(error),
    }

    match shift {
        CoordinateShift::Zero => Ok((start, end)),
        CoordinateShift::Backward(magnitude) => {
            let amount = magnitude.get();
            let Some(new_start) = start.checked_sub(amount) else {
                return Err(CoordinateError::CoordinateUnderflow {
                    domain,
                    operation: CoordinateOperation::BackwardTranslation,
                    lhs: start,
                    rhs: amount,
                });
            };
            let Some(new_end) = end.checked_sub(amount) else {
                return Err(CoordinateError::CoordinateUnderflow {
                    domain,
                    operation: CoordinateOperation::BackwardTranslation,
                    lhs: end,
                    rhs: amount,
                });
            };
            Ok((new_start, new_end))
        }
        CoordinateShift::Forward(magnitude) => {
            let amount = magnitude.get();
            let Some(new_end) = end.checked_add(amount) else {
                return Err(CoordinateError::CoordinateOverflow {
                    domain,
                    operation: CoordinateOperation::ForwardTranslation,
                    lhs: end,
                    rhs: amount,
                });
            };
            let Some(new_start) = start.checked_add(amount) else {
                return Err(CoordinateError::CoordinateOverflow {
                    domain,
                    operation: CoordinateOperation::ForwardTranslation,
                    lhs: start,
                    rhs: amount,
                });
            };
            if new_end > length {
                return Err(CoordinateError::OutOfBounds {
                    domain,
                    operation: CoordinateOperation::ForwardTranslation,
                    start: new_start,
                    end: new_end,
                    length,
                });
            }
            Ok((new_start, new_end))
        }
    }
}

const fn reverse_bounds(
    domain: CoordinateDomain,
    start: u64,
    end: u64,
    length: u64,
) -> Result<(u64, u64), CoordinateError> {
    match validate_interval(
        domain,
        CoordinateOperation::ReverseTransform,
        start,
        end,
        length,
    ) {
        Ok(()) => {}
        Err(error) => return Err(error),
    }

    let Some(new_start) = length.checked_sub(end) else {
        return Err(CoordinateError::CoordinateUnderflow {
            domain,
            operation: CoordinateOperation::ReverseTransform,
            lhs: length,
            rhs: end,
        });
    };
    let Some(new_end) = length.checked_sub(start) else {
        return Err(CoordinateError::CoordinateUnderflow {
            domain,
            operation: CoordinateOperation::ReverseTransform,
            lhs: length,
            rhs: start,
        });
    };
    Ok((new_start, new_end))
}
