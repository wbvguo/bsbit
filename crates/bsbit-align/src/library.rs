//! Paired-library and template-span values used by read alignment.

use core::fmt;

/// Paired-library profile.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub enum PairedLibraryProfile {
    /// Conventional directional R1 OT/OB and R2 CTOT/CTOB pairing.
    Directional,
    /// Non-directional pairing over OT/OB/CTOT/CTOB for both mates.
    NonDirectional,
}

/// A reference-consuming outer template span in bases.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct TemplateSpan(u64);

impl TemplateSpan {
    /// Constructs a typed template-span base count.
    #[must_use]
    pub const fn new(value: u64) -> Self {
        Self(value)
    }

    /// Returns the base count.
    #[must_use]
    pub const fn get(self) -> u64 {
        self.0
    }
}

/// Inclusive accepted template-span bounds.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct TemplateSpanBounds {
    minimum: TemplateSpan,
    maximum: TemplateSpan,
}

impl TemplateSpanBounds {
    /// Validates inclusive minimum/maximum bounds.
    ///
    /// # Errors
    ///
    /// Returns [`PairConstraintError`] when `minimum > maximum`.
    pub const fn new(
        minimum: TemplateSpan,
        maximum: TemplateSpan,
    ) -> Result<Self, PairConstraintError> {
        if minimum.0 > maximum.0 {
            Err(PairConstraintError::InvertedSpanBounds { minimum, maximum })
        } else {
            Ok(Self { minimum, maximum })
        }
    }

    /// Returns the inclusive minimum.
    #[must_use]
    pub const fn minimum(self) -> TemplateSpan {
        self.minimum
    }

    /// Returns the inclusive maximum.
    #[must_use]
    pub const fn maximum(self) -> TemplateSpan {
        self.maximum
    }

    /// Returns whether a span lies inside both inclusive bounds.
    #[must_use]
    pub const fn contains(self, span: TemplateSpan) -> bool {
        self.minimum.0 <= span.0 && span.0 <= self.maximum.0
    }
}

/// Invalid typed pair constraints.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum PairConstraintError {
    /// The inclusive minimum exceeds the maximum.
    InvertedSpanBounds {
        /// Supplied minimum.
        minimum: TemplateSpan,
        /// Supplied maximum.
        maximum: TemplateSpan,
    },
}

impl fmt::Display for PairConstraintError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvertedSpanBounds { minimum, maximum } => write!(
                formatter,
                "minimum template span {} exceeds maximum {}",
                minimum.get(),
                maximum.get()
            ),
        }
    }
}

impl std::error::Error for PairConstraintError {}

/// Pure profile and span constraints for one pair.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct PairConstraints {
    profile: PairedLibraryProfile,
    span_bounds: TemplateSpanBounds,
}

impl PairConstraints {
    /// Constructs already validated pair constraints.
    #[must_use]
    pub const fn new(profile: PairedLibraryProfile, span_bounds: TemplateSpanBounds) -> Self {
        Self {
            profile,
            span_bounds,
        }
    }

    /// Returns the selected profile.
    #[must_use]
    pub const fn profile(self) -> PairedLibraryProfile {
        self.profile
    }

    /// Returns the inclusive span bounds.
    #[must_use]
    pub const fn span_bounds(self) -> TemplateSpanBounds {
        self.span_bounds
    }
}
