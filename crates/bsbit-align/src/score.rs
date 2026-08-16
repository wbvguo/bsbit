//! Checked score values produced by sequence-alignment algorithms.

use core::fmt;

/// A nonnegative bisulfite-aware unit edit cost.
///
/// This is the unit-cost score used by alignment kernels, distinct from MAPQ,
/// affine scores, and quality-derived penalties.
#[derive(Clone, Copy, Debug, Default, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct EditDistance(u64);

impl EditDistance {
    /// Constructs an edit-distance value.
    #[must_use]
    pub const fn new(value: u64) -> Self {
        Self(value)
    }

    /// Returns the unit edit cost.
    #[must_use]
    pub const fn get(self) -> u64 {
        self.0
    }

    /// Adds a nonnegative cost with overflow checking.
    ///
    /// # Errors
    ///
    /// Returns the exact operands when the sum is not representable.
    pub const fn checked_add(self, increment: u64) -> Result<Self, EditDistanceOverflow> {
        match self.0.checked_add(increment) {
            Some(value) => Ok(Self(value)),
            None => Err(EditDistanceOverflow {
                accumulated: self.0,
                increment,
            }),
        }
    }
}

impl fmt::Display for EditDistance {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.fmt(formatter)
    }
}

/// Exact operands from an edit-distance addition that overflowed.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub struct EditDistanceOverflow {
    accumulated: u64,
    increment: u64,
}

impl EditDistanceOverflow {
    /// Returns the value before addition.
    #[must_use]
    pub const fn accumulated(self) -> u64 {
        self.accumulated
    }

    /// Returns the requested nonnegative increment.
    #[must_use]
    pub const fn increment(self) -> u64 {
        self.increment
    }
}

impl fmt::Display for EditDistanceOverflow {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "edit distance addition {} + {} overflowed",
            self.accumulated, self.increment
        )
    }
}

impl std::error::Error for EditDistanceOverflow {}
