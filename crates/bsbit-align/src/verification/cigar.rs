//! Bisulfite-aware replay of structural core CIGAR values.

use core::fmt;

use crate::score::EditDistance;
use bsbit_core::bisulfite::{BaseRelation, CytosineStrand, classify_bases};
use bsbit_core::cigar::{CigarDomain, CigarError, CoreCigar, CoreCigarOp, validate_cigar};
use bsbit_core::sequence::NormalizedSequence;

/// A counted field that overflowed while replaying a CIGAR.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum CigarEvaluationField {
    /// Literal-match columns.
    LiteralMatches,
    /// Bisulfite-compatible columns.
    BisulfiteCompatible,
    /// Literal mismatch columns.
    Mismatches,
    /// Columns containing at least one unknown base.
    UnknownColumns,
    /// Inserted query bases.
    InsertedBases,
    /// Deleted reference bases.
    DeletedBases,
    /// Coalesced insertion/deletion runs.
    GapRuns,
    /// Unit edit cost.
    Distance,
    /// Literal SAM NM count.
    LiteralNm,
}

/// A structural or counter failure while replaying a CIGAR against sequences.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum CigarEvaluationError {
    /// The structural CIGAR did not consume the supplied sequences exactly.
    Structural {
        /// Underlying structural validation failure.
        error: CigarError,
    },
    /// A checked relation or gap counter overflowed.
    CounterOverflow {
        /// Field whose increment or sum overflowed.
        field: CigarEvaluationField,
    },
}

impl fmt::Display for CigarEvaluationError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Structural { error } => write!(formatter, "CIGAR replay rejected: {error}"),
            Self::CounterOverflow { field } => {
                write!(formatter, "CIGAR evaluation counter {field:?} overflowed")
            }
        }
    }
}

impl std::error::Error for CigarEvaluationError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Structural { error } => Some(error),
            Self::CounterOverflow { .. } => None,
        }
    }
}

impl From<CigarError> for CigarEvaluationError {
    fn from(error: CigarError) -> Self {
        Self::Structural { error }
    }
}

/// Counts obtained by replaying a validated CIGAR against two sequences.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct CigarEvaluation {
    distance: EditDistance,
    literal_nm: u64,
    literal_matches: u64,
    bisulfite_compatible: u64,
    mismatches: u64,
    unknown_columns: u64,
    inserted_bases: u64,
    deleted_bases: u64,
    gap_runs: u64,
}

impl CigarEvaluation {
    /// Returns the replayed bisulfite-aware unit edit distance.
    #[must_use]
    pub const fn distance(self) -> EditDistance {
        self.distance
    }

    /// Returns literal SAM NM, including bisulfite-compatible substitutions.
    #[must_use]
    pub const fn literal_nm(self) -> u64 {
        self.literal_nm
    }

    /// Returns literal zero-cost `M` columns.
    #[must_use]
    pub const fn literal_matches(self) -> u64 {
        self.literal_matches
    }

    /// Returns asymmetric bisulfite-compatible zero-cost `M` columns.
    #[must_use]
    pub const fn bisulfite_compatible(self) -> u64 {
        self.bisulfite_compatible
    }

    /// Returns canonical substitution columns with unit cost.
    #[must_use]
    pub const fn mismatches(self) -> u64 {
        self.mismatches
    }

    /// Returns `M` columns containing one or two `N` bases.
    #[must_use]
    pub const fn unknown_columns(self) -> u64 {
        self.unknown_columns
    }

    /// Returns query bases consumed by insertion runs.
    #[must_use]
    pub const fn inserted_bases(self) -> u64 {
        self.inserted_bases
    }

    /// Returns reference bases consumed by deletion runs.
    #[must_use]
    pub const fn deleted_bases(self) -> u64 {
        self.deleted_bases
    }

    /// Returns the number of coalesced insertion/deletion runs.
    #[must_use]
    pub const fn gap_runs(self) -> u64 {
        self.gap_runs
    }
}

/// Replays a CIGAR against `(reference, oriented_query, cytosine_strand)`.
///
/// Consumption is validated before replay. Every `M` column is classified by
/// the bisulfite model, preserving the distinction between literal equality
/// and a conversion-compatible zero-cost column.
///
/// # Errors
///
/// Returns structured consumption, length, or evaluation-counter errors.
pub fn evaluate_cigar(
    cigar: &CoreCigar,
    reference: &NormalizedSequence,
    oriented_query: &NormalizedSequence,
    cytosine_strand: CytosineStrand,
) -> Result<CigarEvaluation, CigarEvaluationError> {
    validate_cigar(cigar, reference.len(), oriented_query.len())?;

    let mut evaluation = EvaluationCounts::default();
    let mut reference_index = 0_usize;
    let mut query_index = 0_usize;

    for (run_index, run) in cigar.runs().iter().copied().enumerate() {
        let operation = run.operation();
        let run_length = run.length();
        let length = usize::try_from(run_length).map_err(|_| CigarEvaluationError::Structural {
            error: CigarError::CigarConsumptionOverflow {
                run_index,
                operation,
                domain: coalescing_domain(operation),
                accumulated: 0,
                run_length,
            },
        })?;
        match operation {
            CoreCigarOp::M => {
                for _ in 0..length {
                    let relation = classify_bases(
                        reference.bases()[reference_index],
                        oriented_query.bases()[query_index],
                        cytosine_strand,
                    );
                    evaluation.increment_relation(relation)?;
                    reference_index += 1;
                    query_index += 1;
                }
            }
            CoreCigarOp::I => {
                evaluation.inserted_bases = checked_count_add(
                    evaluation.inserted_bases,
                    run_length,
                    CigarEvaluationField::InsertedBases,
                )?;
                evaluation.gap_runs =
                    checked_count_add(evaluation.gap_runs, 1, CigarEvaluationField::GapRuns)?;
                query_index += length;
            }
            CoreCigarOp::D => {
                evaluation.deleted_bases = checked_count_add(
                    evaluation.deleted_bases,
                    run_length,
                    CigarEvaluationField::DeletedBases,
                )?;
                evaluation.gap_runs =
                    checked_count_add(evaluation.gap_runs, 1, CigarEvaluationField::GapRuns)?;
                reference_index += length;
            }
        }
    }

    evaluation.finish()
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
struct EvaluationCounts {
    literal_matches: u64,
    bisulfite_compatible: u64,
    mismatches: u64,
    unknown_columns: u64,
    inserted_bases: u64,
    deleted_bases: u64,
    gap_runs: u64,
}

impl EvaluationCounts {
    fn increment_relation(&mut self, relation: BaseRelation) -> Result<(), CigarEvaluationError> {
        let (field, counter) = match relation {
            BaseRelation::LiteralMatch => (
                CigarEvaluationField::LiteralMatches,
                &mut self.literal_matches,
            ),
            BaseRelation::BisulfiteCompatible => (
                CigarEvaluationField::BisulfiteCompatible,
                &mut self.bisulfite_compatible,
            ),
            BaseRelation::Mismatch => (CigarEvaluationField::Mismatches, &mut self.mismatches),
            BaseRelation::Unknown => (
                CigarEvaluationField::UnknownColumns,
                &mut self.unknown_columns,
            ),
        };
        *counter = checked_count_add(*counter, 1, field)?;
        Ok(())
    }

    fn finish(self) -> Result<CigarEvaluation, CigarEvaluationError> {
        let distance = self
            .mismatches
            .checked_add(self.unknown_columns)
            .and_then(|value| value.checked_add(self.inserted_bases))
            .and_then(|value| value.checked_add(self.deleted_bases))
            .ok_or(CigarEvaluationError::CounterOverflow {
                field: CigarEvaluationField::Distance,
            })?;
        let literal_nm = distance.checked_add(self.bisulfite_compatible).ok_or(
            CigarEvaluationError::CounterOverflow {
                field: CigarEvaluationField::LiteralNm,
            },
        )?;
        Ok(CigarEvaluation {
            distance: EditDistance::new(distance),
            literal_nm,
            literal_matches: self.literal_matches,
            bisulfite_compatible: self.bisulfite_compatible,
            mismatches: self.mismatches,
            unknown_columns: self.unknown_columns,
            inserted_bases: self.inserted_bases,
            deleted_bases: self.deleted_bases,
            gap_runs: self.gap_runs,
        })
    }
}

fn checked_count_add(
    accumulated: u64,
    increment: u64,
    field: CigarEvaluationField,
) -> Result<u64, CigarEvaluationError> {
    accumulated
        .checked_add(increment)
        .ok_or(CigarEvaluationError::CounterOverflow { field })
}

const fn coalescing_domain(operation: CoreCigarOp) -> CigarDomain {
    match operation {
        CoreCigarOp::D | CoreCigarOp::M => CigarDomain::Reference,
        CoreCigarOp::I => CigarDomain::Query,
    }
}
