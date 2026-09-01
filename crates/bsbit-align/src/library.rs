//! Library-profile, conversion-pass, and template-span alignment policy.

use core::fmt;

use bsbit_core::bisulfite::BisulfiteStrand;

/// Bisulfite library profile shared by single-end and paired-end alignment.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub enum LibraryProfile {
    /// Conventional directional search over one original conversion pass.
    Directional,
    /// Non-directional search over original and complementary conversion passes.
    NonDirectional,
}

/// Backward-compatible name for callers that constructed paired options.
///
/// New code should use [`LibraryProfile`], which also describes single-end
/// alignment.
pub type PairedLibraryProfile = LibraryProfile;

/// One directional two-strand pass inside a library-profile search.
///
/// A single read uses the pass to select its query projection and hit labels.
/// A read pair uses the original pass in input mate order and the complementary
/// pass in swapped mate order before restoring result order. Layout-specific
/// candidate representation and evidence reduction remain outside this type.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum ConversionPass {
    Original,
    Complementary,
}

const DIRECTIONAL_PASSES: [ConversionPass; 1] = [ConversionPass::Original];
const NON_DIRECTIONAL_PASSES: [ConversionPass; 2] =
    [ConversionPass::Original, ConversionPass::Complementary];

impl LibraryProfile {
    /// Returns the number of directional two-strand passes required by this
    /// profile.
    #[must_use]
    pub const fn conversion_pass_count(self) -> u8 {
        match self {
            Self::Directional => 1,
            Self::NonDirectional => 2,
        }
    }

    pub(crate) const fn conversion_passes(self) -> &'static [ConversionPass] {
        match self {
            Self::Directional => &DIRECTIONAL_PASSES,
            Self::NonDirectional => &NON_DIRECTIONAL_PASSES,
        }
    }
}

impl ConversionPass {
    /// Reports whether a single query uses the complementary projection.
    pub(crate) const fn reverse_complement_query(self) -> bool {
        matches!(self, Self::Complementary)
    }

    /// Reports whether a paired pass swaps input mates and restores them after
    /// mapping through the canonical directional pair executor.
    pub(crate) const fn swaps_mates(self) -> bool {
        matches!(self, Self::Complementary)
    }

    /// Converts one combined-index hit into the molecular strand represented
    /// by this pass.
    ///
    /// The combined index yields OT/OB labels for a projected pass. The
    /// complementary projection reinterprets them as CTOT/CTOB. Unexpected
    /// already-complementary labels are rejected in that pass rather than
    /// being relabelled twice.
    pub(crate) const fn relabel_combined_hit(
        self,
        strand: BisulfiteStrand,
    ) -> Option<BisulfiteStrand> {
        match self {
            Self::Original => Some(strand),
            Self::Complementary => match strand {
                BisulfiteStrand::OT => Some(BisulfiteStrand::CTOT),
                BisulfiteStrand::OB => Some(BisulfiteStrand::CTOB),
                BisulfiteStrand::CTOT | BisulfiteStrand::CTOB => None,
            },
        }
    }
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
    profile: LibraryProfile,
    span_bounds: TemplateSpanBounds,
}

impl PairConstraints {
    /// Constructs already validated pair constraints.
    #[must_use]
    pub const fn new(profile: LibraryProfile, span_bounds: TemplateSpanBounds) -> Self {
        Self {
            profile,
            span_bounds,
        }
    }

    /// Returns the selected profile.
    #[must_use]
    pub const fn profile(self) -> LibraryProfile {
        self.profile
    }

    /// Returns the inclusive span bounds.
    #[must_use]
    pub const fn span_bounds(self) -> TemplateSpanBounds {
        self.span_bounds
    }
}

#[cfg(test)]
mod tests {
    use super::{ConversionPass, LibraryProfile};
    use bsbit_core::bisulfite::BisulfiteStrand;

    #[test]
    fn profiles_expand_to_static_conversion_passes() {
        assert_eq!(
            LibraryProfile::Directional.conversion_passes(),
            [ConversionPass::Original]
        );
        assert_eq!(
            LibraryProfile::NonDirectional.conversion_passes(),
            [ConversionPass::Original, ConversionPass::Complementary]
        );
        assert_eq!(LibraryProfile::Directional.conversion_pass_count(), 1);
        assert_eq!(LibraryProfile::NonDirectional.conversion_pass_count(), 2);
    }

    #[test]
    fn conversion_pass_owns_projection_and_hit_relabelling() {
        assert!(!ConversionPass::Original.reverse_complement_query());
        assert!(ConversionPass::Complementary.reverse_complement_query());
        assert!(!ConversionPass::Original.swaps_mates());
        assert!(ConversionPass::Complementary.swaps_mates());
        assert_eq!(
            ConversionPass::Original.relabel_combined_hit(BisulfiteStrand::OT),
            Some(BisulfiteStrand::OT)
        );
        assert_eq!(
            ConversionPass::Original.relabel_combined_hit(BisulfiteStrand::OB),
            Some(BisulfiteStrand::OB)
        );
        assert_eq!(
            ConversionPass::Complementary.relabel_combined_hit(BisulfiteStrand::OT),
            Some(BisulfiteStrand::CTOT)
        );
        assert_eq!(
            ConversionPass::Complementary.relabel_combined_hit(BisulfiteStrand::OB),
            Some(BisulfiteStrand::CTOB)
        );
        assert_eq!(
            ConversionPass::Complementary.relabel_combined_hit(BisulfiteStrand::CTOT),
            None
        );
        assert_eq!(
            ConversionPass::Complementary.relabel_combined_hit(BisulfiteStrand::CTOB),
            None
        );
    }
}
