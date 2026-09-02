//! Application-layer mapping-result to SAM/BAM record composition.
//!
//! The aligner and index stay independent of output formats, while bsbit-hts
//! owns the shared SAM/BAM format model. This module is the explicit composition
//! boundary between those domains.

use core::fmt;
use core::mem::size_of;

use bsbit_align::extension::VerifiedAlignment;
use bsbit_align::materialize::{
    AlignmentEvaluationError, evaluate_certified_ungapped_alignment, evaluate_ungapped_alignment,
    evaluate_verified_alignment,
};
use bsbit_core::alphabet::Base;
use bsbit_core::bisulfite::{
    AlignmentOrientation, BisulfiteStrand, CytosineStrand, strand_semantics,
};
use bsbit_core::cigar::{CoreCigarOp, CoreCigarRun};
use bsbit_core::coordinate::ReferenceInterval;
use bsbit_core::sequence::NormalizedSequence;
use bsbit_hts::{
    AlignmentAuxiliaryMode, AlignmentCigarOp, AlignmentCigarRun, AlignmentPlacement, AlignmentRead,
    AlignmentRecordAllocation, AlignmentRecordBatch,
    AlignmentRecordError as HtsAlignmentRecordError, AlignmentRecordField, AlignmentRecordLimits,
    AlignmentRecordResource, BorrowedAlignmentRead, BorrowedAlignmentRecord, RecordSegment,
    SAM_MAX_REFERENCE_LENGTH, SamHeader, SamHeaderReference,
};
#[cfg(test)]
use bsbit_hts::{
    AlignmentRecord, MappedAlignmentRecord, RecordMappingQuality, RecordMateLocation,
    RecordReference,
};
use bsbit_index::reference::{ReferenceAccessError, ReferenceIndex};

#[derive(Debug)]
pub(crate) enum RecordBuildError {
    /// A configured cap was exceeded.
    LimitExceeded {
        /// Controlled resource.
        resource: AlignmentRecordResource,
        /// First known complete value outside the cap.
        observed: u64,
        /// Configured cap.
        limit: u64,
    },
    /// A checked logical sum overflowed.
    ArithmeticOverflow {
        /// Resource being counted.
        resource: AlignmentRecordResource,
        /// Value before addition.
        current: u64,
        /// Attempted increment.
        increment: u64,
    },
    /// A bounded allocation could not be reserved.
    AllocationFailed {
        /// Allocation site.
        allocation: AlignmentRecordAllocation,
        /// Requested bytes or elements.
        requested: u64,
    },
    /// A retained soft-clipped query was not an exact prefix of its full read.
    SoftClippedSequenceMismatch {
        /// One-based mate number.
        mate: u8,
    },
    /// An owner-bound alignment did not resolve through the supplied reference.
    ReferenceAccess {
        /// Underlying exact-owner access failure.
        source: ReferenceAccessError,
    },
    /// A record integer cannot be represented in its SAM field.
    FieldOutOfRange {
        /// Failed field.
        field: AlignmentRecordField,
        /// Observed unsigned magnitude.
        value: u64,
    },
    /// Authoritative alignment-fact evaluation failed.
    AlignmentEvaluation {
        /// Alignment-layer failure.
        source: AlignmentEvaluationError,
    },
    /// Literal SAM replay disagreed with the alignment-layer NM fact.
    AlignmentLiteralNmMismatch {
        /// Alignment-layer literal NM.
        expected: u64,
        /// Format replay literal NM.
        observed: u64,
    },
    /// A concordant pair did not satisfy the same-reference record invariant.
    ConcordantReferenceMismatch,
    /// The shared HTS alignment-record model rejected composed format fields.
    Format { source: HtsAlignmentRecordError },
}

impl fmt::Display for RecordBuildError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::LimitExceeded {
                resource,
                observed,
                limit,
            } => write!(
                formatter,
                "alignment-record resource {resource:?} observed {observed}, exceeding {limit}"
            ),
            Self::ArithmeticOverflow {
                resource,
                current,
                increment,
            } => write!(
                formatter,
                "alignment-record resource {resource:?} overflowed: {current} + {increment}"
            ),
            Self::AllocationFailed {
                allocation,
                requested,
            } => write!(
                formatter,
                "failed to reserve {requested} bytes/elements for {allocation:?}"
            ),
            Self::SoftClippedSequenceMismatch { mate } => write!(
                formatter,
                "mate {mate} soft-clipped query is not an exact nonempty prefix of the full read"
            ),
            Self::ReferenceAccess { source } => {
                write!(formatter, "alignment reference access failed: {source}")
            }
            Self::FieldOutOfRange { field, value } => {
                write!(
                    formatter,
                    "alignment-record field {field:?} cannot represent {value}"
                )
            }
            Self::AlignmentEvaluation { source } => {
                write!(formatter, "alignment evaluation failed: {source}")
            }
            Self::AlignmentLiteralNmMismatch { expected, observed } => write!(
                formatter,
                "alignment literal NM {expected} differs from format replay NM {observed}"
            ),
            Self::ConcordantReferenceMismatch => {
                formatter.write_str("concordant pair resolved to different reference sequences")
            }
            Self::Format { source } => {
                write!(formatter, "alignment record validation failed: {source}")
            }
        }
    }
}

impl std::error::Error for RecordBuildError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::ReferenceAccess { source } => Some(source),
            Self::AlignmentEvaluation { source } => Some(source),
            Self::Format { source } => Some(source),
            _ => None,
        }
    }
}

impl From<HtsAlignmentRecordError> for RecordBuildError {
    fn from(source: HtsAlignmentRecordError) -> Self {
        match source {
            HtsAlignmentRecordError::LimitExceeded {
                resource,
                observed,
                limit,
            } => Self::LimitExceeded {
                resource,
                observed,
                limit,
            },
            HtsAlignmentRecordError::ArithmeticOverflow {
                resource,
                current,
                increment,
            } => Self::ArithmeticOverflow {
                resource,
                current,
                increment,
            },
            HtsAlignmentRecordError::AllocationFailed {
                allocation,
                requested,
            } => Self::AllocationFailed {
                allocation,
                requested,
            },
            HtsAlignmentRecordError::FieldOutOfRange { field, value } => {
                Self::FieldOutOfRange { field, value }
            }
            source => Self::Format { source },
        }
    }
}

mod auxiliary;
mod direct;
#[cfg(test)]
mod materialized;

#[cfg(test)]
pub(crate) use auxiliary::bismark_methylation_call;
pub(crate) use direct::{PairedRecordComposer, SingleRecordComposer};
#[cfg(test)]
pub(crate) use materialized::{
    build_single_alignment_record, build_single_alignment_record_with_auxiliary_mode,
};

pub(crate) fn build_sam_header(
    reference: &ReferenceIndex,
    limits: AlignmentRecordLimits,
) -> Result<SamHeader, RecordBuildError> {
    let capacity = usize::try_from(reference.contig_count()).map_err(|_| {
        RecordBuildError::AllocationFailed {
            allocation: AlignmentRecordAllocation::HeaderReferences,
            requested: reference.contig_count(),
        }
    })?;
    let mut entries = Vec::new();
    entries
        .try_reserve_exact(capacity)
        .map_err(|_| RecordBuildError::AllocationFailed {
            allocation: AlignmentRecordAllocation::HeaderReferences,
            requested: reference.contig_count(),
        })?;
    for ordinal in 0..reference.contig_count() {
        let id = reference
            .contig_id(ordinal)
            .map_err(|source| RecordBuildError::ReferenceAccess { source })?;
        let contig = reference
            .resolve_contig(&id)
            .map_err(|source| RecordBuildError::ReferenceAccess { source })?;
        entries.push(SamHeaderReference::new(
            ordinal,
            contig.name(),
            contig.sequence().len(),
        )?);
    }
    SamHeader::new(entries, limits).map_err(Into::into)
}

#[cfg(test)]
pub(crate) fn check_limit(
    observed: u64,
    limit: u64,
    resource: AlignmentRecordResource,
) -> Result<(), RecordBuildError> {
    if observed > limit {
        Err(RecordBuildError::LimitExceeded {
            resource,
            observed,
            limit,
        })
    } else {
        Ok(())
    }
}

pub(crate) fn checked_add_resource(
    current: u64,
    increment: u64,
    resource: AlignmentRecordResource,
) -> Result<u64, RecordBuildError> {
    current
        .checked_add(increment)
        .ok_or(RecordBuildError::ArithmeticOverflow {
            resource,
            current,
            increment,
        })
}

pub(crate) const fn decimal_digits(mut value: u64) -> u64 {
    let mut digits = 1;
    while value >= 10 {
        value /= 10;
        digits += 1;
    }
    digits
}

pub(crate) fn append_u64(output: &mut Vec<u8>, mut value: u64) {
    let mut digits = [0_u8; 20];
    let mut start = digits.len();
    loop {
        start -= 1;
        digits[start] = b'0' + u8::try_from(value % 10).expect("single decimal digit");
        value /= 10;
        if value == 0 {
            break;
        }
    }
    output.extend_from_slice(&digits[start..]);
}

pub(crate) const fn storage_len(length: usize) -> u64 {
    length as u64
}

fn storage_count(
    value: u64,
    allocation: AlignmentRecordAllocation,
) -> Result<usize, RecordBuildError> {
    usize::try_from(value).map_err(|_| RecordBuildError::AllocationFailed {
        allocation,
        requested: value.saturating_mul(size_of::<u8>() as u64),
    })
}

#[cfg(test)]
#[path = "../../tests/whitebox/alignment_record.rs"]
mod tests;

#[cfg(test)]
#[path = "../../tests/whitebox/alignment_record_fuzz_smoke.rs"]
mod fuzz_tests;

#[cfg(test)]
#[path = "../../tests/whitebox/bam_fields.rs"]
mod bam_tests;

#[cfg(test)]
#[path = "../../tests/whitebox/record_fixture.rs"]
mod record_fixture;
