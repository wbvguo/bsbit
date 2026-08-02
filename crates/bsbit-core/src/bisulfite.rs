//! Bisulfite strand, conversion, and matching semantics.
//!
//! Alignment orientation, cytosine evidence strand, molecular bisulfite
//! strand, and three-letter conversion are intentionally separate axes. The
//! total four-strand table in this module is the only supported derivation
//! between those stable chemistry concepts. Library and mate admissibility are
//! alignment policy and live in `bsbit-align`.

use core::fmt;

use crate::alphabet::Base;

/// Query orientation relative to forward FASTA coordinates.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub enum AlignmentOrientation {
    /// The oriented query runs with the forward reference.
    Forward,
    /// The oriented query is the reverse complement of the input read.
    Reverse,
}

impl fmt::Display for AlignmentOrientation {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::Forward => "Forward",
            Self::Reverse => "Reverse",
        })
    }
}

/// The genomic cytosine strand that supplies bisulfite evidence.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub enum CytosineStrand {
    /// Top/reference C sites: reference C may match query T.
    Top,
    /// Bottom/reference G sites: reference G may match query A.
    Bottom,
}

impl fmt::Display for CytosineStrand {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::Top => "Top",
            Self::Bottom => "Bottom",
        })
    }
}

/// A molecular/derived strand in the four-strand bisulfite vocabulary.
#[allow(clippy::upper_case_acronyms)]
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub enum BisulfiteStrand {
    /// Original top strand.
    OT,
    /// Original bottom strand.
    OB,
    /// Complementary to original top strand.
    CTOT,
    /// Complementary to original bottom strand.
    CTOB,
}

impl BisulfiteStrand {
    /// Every strand in canonical diagnostic order.
    pub const ALL: [Self; 4] = [Self::OT, Self::OB, Self::CTOT, Self::CTOB];
}

/// A three-letter sequence transformation used only for candidate search.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub enum ThreeLetterConversion {
    /// Convert C to T, preserving every other base.
    CToT,
    /// Convert G to A, preserving every other base.
    GToA,
}

impl fmt::Display for ThreeLetterConversion {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::CToT => "CToT",
            Self::GToA => "GToA",
        })
    }
}

/// The validated semantic axes associated with one bisulfite strand.
///
/// Its search conversion applies only after the raw query has been oriented
/// left-to-right in forward-reference coordinates. A representation that keeps
/// a reverse-lane query in sequencing order must derive the dual conversion;
/// it must not reuse this value unchanged.
///
/// Fields are private so an inconsistent combination cannot be constructed.
/// Use [`strand_semantics`] or [`validate_strand_assignment`].
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub struct StrandSemantics {
    strand: BisulfiteStrand,
    orientation: AlignmentOrientation,
    cytosine_strand: CytosineStrand,
    search_conversion: ThreeLetterConversion,
}

impl StrandSemantics {
    /// Returns the molecular/derived strand.
    #[must_use]
    pub const fn strand(self) -> BisulfiteStrand {
        self.strand
    }

    /// Returns the alignment orientation relative to the forward reference.
    #[must_use]
    pub const fn orientation(self) -> AlignmentOrientation {
        self.orientation
    }

    /// Returns the genomic cytosine evidence strand.
    #[must_use]
    pub const fn cytosine_strand(self) -> CytosineStrand {
        self.cytosine_strand
    }

    /// Returns the oriented-query transformation used for candidate search.
    ///
    /// The query must already be oriented left-to-right in forward-reference
    /// coordinates. A raw sequencing-order reverse-lane view uses the dual.
    #[must_use]
    pub const fn search_conversion(self) -> ThreeLetterConversion {
        self.search_conversion
    }
}

/// Which axis failed [`validate_strand_assignment`].
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum StrandAssignmentAxis {
    /// The supplied alignment orientation was inconsistent.
    Orientation,
    /// The supplied cytosine evidence strand was inconsistent.
    CytosineStrand,
}

/// A structured inconsistent-strand-assignment error.
///
/// Orientation is validated before cytosine strand, so [`Self::axis`] is
/// deterministic when both supplied axes are wrong.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct InconsistentStrandAssignment {
    /// The strand being validated.
    pub strand: BisulfiteStrand,
    /// The first inconsistent axis in the specified validation order.
    pub axis: StrandAssignmentAxis,
    /// The supplied orientation.
    pub supplied_orientation: AlignmentOrientation,
    /// The expected orientation.
    pub expected_orientation: AlignmentOrientation,
    /// The supplied cytosine evidence strand.
    pub supplied_cytosine_strand: CytosineStrand,
    /// The expected cytosine evidence strand.
    pub expected_cytosine_strand: CytosineStrand,
}

impl fmt::Display for InconsistentStrandAssignment {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "inconsistent {:?} assignment for {}: supplied ({:?}, {:?}), expected ({:?}, {:?})",
            self.axis,
            self.strand,
            self.supplied_orientation,
            self.supplied_cytosine_strand,
            self.expected_orientation,
            self.expected_cytosine_strand,
        )
    }
}

impl std::error::Error for InconsistentStrandAssignment {}

/// Returns the exact orientation, cytosine strand, and oriented-query search
/// conversion for a bisulfite strand.
///
/// The returned conversion applies after orienting the query to
/// forward-reference coordinates. It is not a raw reverse-lane projection.
///
#[must_use]
pub const fn strand_semantics(strand: BisulfiteStrand) -> StrandSemantics {
    match strand {
        BisulfiteStrand::OT => StrandSemantics {
            strand,
            orientation: AlignmentOrientation::Forward,
            cytosine_strand: CytosineStrand::Top,
            search_conversion: ThreeLetterConversion::CToT,
        },
        BisulfiteStrand::OB => StrandSemantics {
            strand,
            orientation: AlignmentOrientation::Reverse,
            cytosine_strand: CytosineStrand::Bottom,
            search_conversion: ThreeLetterConversion::GToA,
        },
        BisulfiteStrand::CTOT => StrandSemantics {
            strand,
            orientation: AlignmentOrientation::Reverse,
            cytosine_strand: CytosineStrand::Top,
            search_conversion: ThreeLetterConversion::CToT,
        },
        BisulfiteStrand::CTOB => StrandSemantics {
            strand,
            orientation: AlignmentOrientation::Forward,
            cytosine_strand: CytosineStrand::Bottom,
            search_conversion: ThreeLetterConversion::GToA,
        },
    }
}

/// Validates explicitly supplied axes against the canonical four-strand table.
///
/// Orientation is checked before cytosine strand. The returned value cannot
/// contain an inconsistent combination.
///
/// # Errors
///
/// Returns [`InconsistentStrandAssignment`] with both supplied and expected
/// axes when either input axis disagrees with the strand table.
pub const fn validate_strand_assignment(
    strand: BisulfiteStrand,
    orientation: AlignmentOrientation,
    cytosine_strand: CytosineStrand,
) -> Result<StrandSemantics, InconsistentStrandAssignment> {
    let expected = strand_semantics(strand);
    let axis = if !orientation_eq(orientation, expected.orientation) {
        Some(StrandAssignmentAxis::Orientation)
    } else if !cytosine_strand_eq(cytosine_strand, expected.cytosine_strand) {
        Some(StrandAssignmentAxis::CytosineStrand)
    } else {
        None
    };

    match axis {
        None => Ok(expected),
        Some(axis) => Err(InconsistentStrandAssignment {
            strand,
            axis,
            supplied_orientation: orientation,
            expected_orientation: expected.orientation,
            supplied_cytosine_strand: cytosine_strand,
            expected_cytosine_strand: expected.cytosine_strand,
        }),
    }
}

const fn orientation_eq(left: AlignmentOrientation, right: AlignmentOrientation) -> bool {
    matches!(
        (left, right),
        (AlignmentOrientation::Forward, AlignmentOrientation::Forward)
            | (AlignmentOrientation::Reverse, AlignmentOrientation::Reverse)
    )
}

const fn cytosine_strand_eq(left: CytosineStrand, right: CytosineStrand) -> bool {
    matches!(
        (left, right),
        (CytosineStrand::Top, CytosineStrand::Top)
            | (CytosineStrand::Bottom, CytosineStrand::Bottom)
    )
}

/// The conversion-aware relation between a forward-reference base and an
/// oriented-query base.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum BaseRelation {
    /// Canonical bases are literally equal; cost zero.
    LiteralMatch,
    /// An allowed asymmetric C-to-T or G-to-A observation; cost zero.
    BisulfiteCompatible,
    /// A canonical substitution not allowed by bisulfite chemistry; cost one.
    Mismatch,
    /// At least one base is `N`; cost one, including `N`/`N`.
    Unknown,
}

impl BaseRelation {
    /// Returns the reference unit substitution cost.
    #[must_use]
    pub const fn cost(self) -> u64 {
        match self {
            Self::LiteralMatch | Self::BisulfiteCompatible => 0,
            Self::Mismatch | Self::Unknown => 1,
        }
    }

    /// Returns whether this relation has zero substitution cost.
    #[must_use]
    pub const fn is_zero_cost(self) -> bool {
        matches!(self, Self::LiteralMatch | Self::BisulfiteCompatible)
    }
}

/// Classifies `(reference_base, oriented_query_base)` for a cytosine strand.
///
/// Argument order is semantically significant. The allowed conversion pairs
/// are reference C/query T on [`CytosineStrand::Top`] and reference G/query A
/// on [`CytosineStrand::Bottom`]. The reverse pairs remain mismatches.
#[must_use]
pub const fn classify_bases(
    reference_base: Base,
    oriented_query_base: Base,
    cytosine_strand: CytosineStrand,
) -> BaseRelation {
    if reference_base.is_unknown() || oriented_query_base.is_unknown() {
        BaseRelation::Unknown
    } else if base_eq(reference_base, oriented_query_base) {
        BaseRelation::LiteralMatch
    } else {
        match (cytosine_strand, reference_base, oriented_query_base) {
            (CytosineStrand::Top, Base::C, Base::T)
            | (CytosineStrand::Bottom, Base::G, Base::A) => BaseRelation::BisulfiteCompatible,
            _ => BaseRelation::Mismatch,
        }
    }
}

const fn base_eq(left: Base, right: Base) -> bool {
    matches!(
        (left, right),
        (Base::A, Base::A)
            | (Base::C, Base::C)
            | (Base::G, Base::G)
            | (Base::T, Base::T)
            | (Base::N, Base::N)
    )
}

impl fmt::Display for BisulfiteStrand {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::OT => "OT",
            Self::OB => "OB",
            Self::CTOT => "CTOT",
            Self::CTOB => "CTOB",
        })
    }
}
