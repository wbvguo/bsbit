//! Shared record-composition errors and checked storage arithmetic.

use core::fmt;
use core::mem::size_of;

use bsbit_align::materialize::AlignmentEvaluationError;
use bsbit_hts::{
    AlignmentRecordAllocation, AlignmentRecordError as HtsAlignmentRecordError,
    AlignmentRecordField, AlignmentRecordResource,
};
use bsbit_index::reference::ReferenceAccessError;

#[derive(Debug)]
pub(crate) enum RecordBuildError {
    LimitExceeded {
        resource: AlignmentRecordResource,
        observed: u64,
        limit: u64,
    },
    ArithmeticOverflow {
        resource: AlignmentRecordResource,
        current: u64,
        increment: u64,
    },
    AllocationFailed {
        allocation: AlignmentRecordAllocation,
        requested: u64,
    },
    SoftClippedSequenceMismatch {
        mate: u8,
    },
    ReferenceAccess {
        source: ReferenceAccessError,
    },
    FieldOutOfRange {
        field: AlignmentRecordField,
        value: u64,
    },
    AlignmentEvaluation {
        source: AlignmentEvaluationError,
    },
    AlignmentLiteralNmMismatch {
        expected: u64,
        observed: u64,
    },
    ConcordantReferenceMismatch,
    Format {
        source: HtsAlignmentRecordError,
    },
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
            Self::FieldOutOfRange { field, value } => write!(
                formatter,
                "alignment-record field {field:?} cannot represent {value}"
            ),
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

pub(crate) fn storage_count(
    value: u64,
    allocation: AlignmentRecordAllocation,
) -> Result<usize, RecordBuildError> {
    usize::try_from(value).map_err(|_| RecordBuildError::AllocationFailed {
        allocation,
        requested: value.saturating_mul(size_of::<u8>() as u64),
    })
}
