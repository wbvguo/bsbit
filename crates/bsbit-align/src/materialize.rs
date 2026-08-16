//! Verified alignment facts for output-format composition.

use core::fmt;

use crate::extension::{ExtensionError, VerifiedAlignment, traceback_retained_placement_banded};
use crate::score::EditDistance;
use crate::verification::cigar::{CigarEvaluationError, evaluate_cigar};
use bsbit_core::alphabet::Base;
use bsbit_core::bisulfite::{AlignmentOrientation, BisulfiteStrand, strand_semantics};
use bsbit_core::cigar::CoreCigar;
use bsbit_core::coordinate::ReferenceInterval;
use bsbit_core::sequence::NormalizedSequence;
use bsbit_index::reference::{ContigId, ReferenceAccessError, ReferenceIndex};

/// Conversion-aware and literal edit facts for one alignment.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct AlignmentEvaluation {
    distance: EditDistance,
    literal_nm: u64,
}

impl AlignmentEvaluation {
    /// Returns the verified conversion-aware unit edit distance.
    #[must_use]
    pub const fn distance(self) -> EditDistance {
        self.distance
    }

    /// Returns literal SAM NM, including conversion-compatible substitutions.
    #[must_use]
    pub const fn literal_nm(self) -> u64 {
        self.literal_nm
    }
}

/// Failure while deriving output facts from a retained alignment.
#[derive(Debug)]
pub enum AlignmentEvaluationError {
    /// The owner-bound reference contig could not be resolved.
    Reference(ReferenceAccessError),
    /// A checked reference interval could not address local storage.
    IntervalNotRepresentable,
    /// Structural or arithmetic CIGAR replay failed.
    Cigar(CigarEvaluationError),
    /// Replay disagreed with the mapper's verified distance.
    DistanceMismatch {
        /// Distance retained by mapping.
        expected: EditDistance,
        /// Distance obtained by authoritative replay.
        observed: EditDistance,
    },
    /// Ungapped evaluation requires equal nonzero sequence lengths.
    InvalidUngappedLengths,
}

impl fmt::Display for AlignmentEvaluationError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Reference(error) => {
                write!(formatter, "alignment reference access failed: {error}")
            }
            Self::IntervalNotRepresentable => {
                formatter.write_str("alignment interval is not locally addressable")
            }
            Self::Cigar(error) => write!(formatter, "alignment CIGAR evaluation failed: {error}"),
            Self::DistanceMismatch { expected, observed } => write!(
                formatter,
                "alignment distance mismatch: expected {}, observed {}",
                expected.get(),
                observed.get()
            ),
            Self::InvalidUngappedLengths => {
                formatter.write_str("ungapped alignment sequences must have equal nonzero lengths")
            }
        }
    }
}

impl std::error::Error for AlignmentEvaluationError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Reference(error) => Some(error),
            Self::Cigar(error) => Some(error),
            Self::IntervalNotRepresentable
            | Self::DistanceMismatch { .. }
            | Self::InvalidUngappedLengths => None,
        }
    }
}

/// Recovers the canonical traceback for a retained read placement.
///
/// Search stores a conservative distance. Canonical materialization checks
/// every admissible exact band from that bound down to zero and returns the
/// first verified path.
///
/// # Errors
///
/// Returns the final [`ExtensionError`] when no exact band can materialize the
/// retained placement.
///
/// # Panics
///
/// The inclusive `u8` distance range always performs at least one attempt; a
/// panic would therefore indicate an internal control-flow invariant failure.
pub fn traceback_read_placement(
    reference: &ReferenceIndex,
    query: &NormalizedSequence,
    contig: &ContigId,
    interval: ReferenceInterval,
    strand: BisulfiteStrand,
    conservative_distance: u8,
) -> Result<VerifiedAlignment, ExtensionError> {
    let mut last_error = None;
    for distance in (0..=conservative_distance).rev() {
        match traceback_retained_placement_banded(
            reference,
            query,
            contig,
            interval,
            strand,
            EditDistance::new(u64::from(distance)),
        ) {
            Ok(alignment) => return Ok(alignment),
            Err(error) => last_error = Some(error),
        }
    }
    Err(last_error.expect("the inclusive placement-distance range is non-empty"))
}

/// Evaluates one canonical verified alignment against its owner-bound contig.
///
/// # Errors
///
/// Returns [`AlignmentEvaluationError`] for reference access, interval,
/// structural CIGAR, counter, or verified-distance disagreement.
pub fn evaluate_verified_alignment(
    reference: &ReferenceIndex,
    raw_query: &NormalizedSequence,
    alignment: &VerifiedAlignment,
) -> Result<AlignmentEvaluation, AlignmentEvaluationError> {
    let contig = reference
        .resolve_contig(alignment.contig())
        .map_err(AlignmentEvaluationError::Reference)?;
    let start = usize::try_from(alignment.interval().start())
        .map_err(|_| AlignmentEvaluationError::IntervalNotRepresentable)?;
    let end = usize::try_from(alignment.interval().end())
        .map_err(|_| AlignmentEvaluationError::IntervalNotRepresentable)?;
    let reference_bases = contig
        .sequence()
        .bases()
        .get(start..end)
        .ok_or(AlignmentEvaluationError::IntervalNotRepresentable)?;
    let reference_interval = NormalizedSequence::from_bases(reference_bases.iter().copied());
    let reversed;
    let oriented_query = match alignment.orientation() {
        AlignmentOrientation::Forward => raw_query,
        AlignmentOrientation::Reverse => {
            reversed = raw_query.reverse_complement();
            &reversed
        }
    };
    let evaluation = evaluate_cigar(
        alignment.cigar(),
        &reference_interval,
        oriented_query,
        alignment.cytosine_strand(),
    )
    .map_err(AlignmentEvaluationError::Cigar)?;
    if evaluation.distance() != alignment.distance() {
        return Err(AlignmentEvaluationError::DistanceMismatch {
            expected: alignment.distance(),
            observed: evaluation.distance(),
        });
    }
    Ok(AlignmentEvaluation {
        distance: evaluation.distance(),
        literal_nm: evaluation.literal_nm(),
    })
}

/// Evaluates an ungapped paired-end placement without output-format policy.
///
/// # Errors
///
/// Returns [`AlignmentEvaluationError`] when the two sequences do not have the
/// same nonzero length or the authoritative all-match replay fails.
pub fn evaluate_ungapped_alignment(
    reference_bases: &[Base],
    raw_query: &[Base],
    strand: BisulfiteStrand,
) -> Result<AlignmentEvaluation, AlignmentEvaluationError> {
    if reference_bases.is_empty() || reference_bases.len() != raw_query.len() {
        return Err(AlignmentEvaluationError::InvalidUngappedLengths);
    }
    let reference = NormalizedSequence::from_bases(reference_bases.iter().copied());
    let raw_query = NormalizedSequence::from_bases(raw_query.iter().copied());
    let reversed;
    let oriented_query = match strand_semantics(strand).orientation() {
        AlignmentOrientation::Forward => &raw_query,
        AlignmentOrientation::Reverse => {
            reversed = raw_query.reverse_complement();
            &reversed
        }
    };
    let cigar = CoreCigar::all_matches(raw_query.len());
    let evaluation = evaluate_cigar(
        &cigar,
        &reference,
        oriented_query,
        strand_semantics(strand).cytosine_strand(),
    )
    .map_err(AlignmentEvaluationError::Cigar)?;
    Ok(AlignmentEvaluation {
        distance: evaluation.distance(),
        literal_nm: evaluation.literal_nm(),
    })
}
