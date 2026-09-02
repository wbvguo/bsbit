//! Verified alignment facts for output-format composition.

use core::fmt;

use crate::extension::{ExtensionError, VerifiedAlignment, traceback_retained_placement_banded};
use crate::score::EditDistance;
use crate::verification::cigar::{CigarEvaluationError, evaluate_cigar};
use crate::verification::distance::DistanceError;
use crate::verification::prefix_filter::ungapped_traceback_at_most_two_certified_cached_nm;
use crate::verification::ungapped::MAX_UNGAPPED_QUERY_BASES;
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
    /// Certified ungapped traceback construction failed.
    Distance(DistanceError),
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
            Self::Distance(error) => {
                write!(formatter, "certified ungapped alignment failed: {error}")
            }
        }
    }
}

impl std::error::Error for AlignmentEvaluationError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Reference(error) => Some(error),
            Self::Cigar(error) => Some(error),
            Self::Distance(error) => Some(error),
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

/// Evaluates a canonical all-match alignment when the exact ungapped path can
/// be certified without an equally scoring shifted-gap path.
///
/// `Ok(None)` asks the caller to use the full canonical traceback. This occurs
/// above distance two, for an equal-distance shifted-gap tie, or above the
/// fixed read length supported by the aligner hot path.
///
/// # Errors
///
/// Returns [`AlignmentEvaluationError`] for unequal or empty spans and for a
/// failed certified traceback result.
pub fn evaluate_certified_ungapped_alignment(
    reference_bases: &[Base],
    raw_query: &[Base],
    strand: BisulfiteStrand,
) -> Result<Option<AlignmentEvaluation>, AlignmentEvaluationError> {
    if reference_bases.is_empty() || reference_bases.len() != raw_query.len() {
        return Err(AlignmentEvaluationError::InvalidUngappedLengths);
    }
    if raw_query.len() > MAX_UNGAPPED_QUERY_BASES {
        return Ok(None);
    }
    let semantics = strand_semantics(strand);
    let mut reverse_storage = [Base::N; MAX_UNGAPPED_QUERY_BASES];
    let oriented_query = match semantics.orientation() {
        AlignmentOrientation::Forward => raw_query,
        AlignmentOrientation::Reverse => {
            for (output, base) in reverse_storage.iter_mut().zip(raw_query.iter().rev()) {
                *output = base.complement();
            }
            &reverse_storage[..raw_query.len()]
        }
    };
    ungapped_traceback_at_most_two_certified_cached_nm(
        reference_bases,
        oriented_query,
        semantics.cytosine_strand(),
    )
    .map_err(AlignmentEvaluationError::Distance)
    .map(|certified| {
        certified.map(|(traceback, literal_nm)| AlignmentEvaluation {
            distance: traceback.distance(),
            literal_nm,
        })
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn raw_query_certification_matches_the_oriented_certificate_exhaustively() {
        for strand in [
            BisulfiteStrand::OT,
            BisulfiteStrand::OB,
            BisulfiteStrand::CTOT,
            BisulfiteStrand::CTOB,
        ] {
            let semantics = strand_semantics(strand);
            for reference_first in Base::ALL {
                for reference_second in Base::ALL {
                    let reference = [reference_first, reference_second];
                    for query_first in Base::ALL {
                        for query_second in Base::ALL {
                            let raw_query = [query_first, query_second];
                            let oriented = match semantics.orientation() {
                                AlignmentOrientation::Forward => raw_query.to_vec(),
                                AlignmentOrientation::Reverse => raw_query
                                    .iter()
                                    .rev()
                                    .map(|base| base.complement())
                                    .collect(),
                            };
                            let expected = ungapped_traceback_at_most_two_certified_cached_nm(
                                &reference,
                                &oriented,
                                semantics.cytosine_strand(),
                            )
                            .expect("short certificate is representable")
                            .map(|(traceback, literal_nm)| AlignmentEvaluation {
                                distance: traceback.distance(),
                                literal_nm,
                            });
                            assert_eq!(
                                evaluate_certified_ungapped_alignment(
                                    &reference, &raw_query, strand,
                                )
                                .expect("valid short span"),
                                expected,
                                "strand={strand:?}, reference={reference:?}, query={raw_query:?}",
                            );
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn raw_query_certification_declines_unsupported_lengths() {
        assert!(matches!(
            evaluate_certified_ungapped_alignment(&[], &[], BisulfiteStrand::OT),
            Err(AlignmentEvaluationError::InvalidUngappedLengths)
        ));
        assert!(matches!(
            evaluate_certified_ungapped_alignment(
                &[Base::A],
                &[Base::A, Base::A],
                BisulfiteStrand::OT,
            ),
            Err(AlignmentEvaluationError::InvalidUngappedLengths)
        ));
        assert_eq!(
            evaluate_certified_ungapped_alignment(
                &[Base::A; MAX_UNGAPPED_QUERY_BASES + 1],
                &[Base::A; MAX_UNGAPPED_QUERY_BASES + 1],
                BisulfiteStrand::OT,
            )
            .expect("long equal span is a fallback"),
            None
        );
    }
}
