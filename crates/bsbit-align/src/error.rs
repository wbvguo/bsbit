use core::fmt;

use crate::verification::NarrowBandedError;

use crate::read_mapping_limits::{MAX_READ_BASES, VERIFICATION_BATCH};

/// Failure while searching or verifying a read-alignment batch.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum AlignmentError {
    /// A read length falls outside the supported mapping range.
    UnsupportedReadLength {
        /// Observed read length.
        length: usize,
    },
    /// Narrow-band alignment verification failed.
    Verification(NarrowBandedError),
    /// A located contig ordinal was not present in the reference.
    InvalidContigOrdinal {
        /// Missing zero-based contig ordinal.
        ordinal: u64,
    },
    /// Candidate-start arithmetic overflowed.
    CandidateCoordinateOverflow {
        /// Candidate start involved in the overflow.
        start: u64,
    },
    /// Candidate endpoint arithmetic overflowed.
    CandidateEndpointOverflow {
        /// Candidate start.
        start: u64,
        /// Candidate reference length.
        length: usize,
    },
    /// The minimum template span exceeds the maximum.
    InvertedTemplateSpan {
        /// Requested minimum span.
        minimum: u64,
        /// Requested maximum span.
        maximum: u64,
    },
    /// The requested edit budget exceeds the supported maximum.
    UnsupportedEditDistance {
        /// Requested edit-distance budget.
        requested: u8,
        /// Largest supported budget.
        maximum: u8,
    },
    /// A verification batch exceeded the fixed vectorized batch size.
    VerificationBatchSize {
        /// Observed batch size.
        observed: usize,
    },
    /// Candidate and verification output lengths did not agree.
    VerificationOutputSize {
        /// Candidate count.
        candidates: usize,
        /// Output slot count.
        output: usize,
    },
    /// A verification batch mixed bisulfite strands.
    MixedVerificationStrands,
    /// The located-coordinate counter overflowed.
    LocatedCountOverflow,
    /// Combined-index search failed.
    CombinedIndex,
}

impl From<NarrowBandedError> for AlignmentError {
    fn from(error: NarrowBandedError) -> Self {
        Self::Verification(error)
    }
}

impl fmt::Display for AlignmentError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnsupportedReadLength { length } => write!(
                formatter,
                "alignment supports read lengths 3 through {MAX_READ_BASES}, got {length}"
            ),
            Self::Verification(error) => error.fmt(formatter),
            Self::InvalidContigOrdinal { ordinal } => {
                write!(
                    formatter,
                    "alignment candidate contig ordinal {ordinal} is invalid"
                )
            }
            Self::CandidateCoordinateOverflow { start } => {
                write!(formatter, "alignment candidate start {start} exceeds usize")
            }
            Self::CandidateEndpointOverflow { start, length } => write!(
                formatter,
                "alignment candidate endpoint {start}+{length} overflowed"
            ),
            Self::InvertedTemplateSpan { minimum, maximum } => write!(
                formatter,
                "paired-end minimum template span {minimum} exceeds maximum {maximum}"
            ),
            Self::UnsupportedEditDistance { requested, maximum } => write!(
                formatter,
                "alignment maximum edit distance {requested} exceeds supported maximum {maximum}"
            ),
            Self::VerificationBatchSize { observed } => write!(
                formatter,
                "alignment verification batch size {observed} is outside 1..={VERIFICATION_BATCH}"
            ),
            Self::VerificationOutputSize { candidates, output } => write!(
                formatter,
                "alignment verification candidates/output differ: {candidates}/{output}"
            ),
            Self::MixedVerificationStrands => {
                formatter.write_str("alignment verification batch mixes bisulfite strands")
            }
            Self::LocatedCountOverflow => {
                formatter.write_str("alignment located-row count overflowed")
            }
            Self::CombinedIndex => formatter.write_str("alignment combined-index query failed"),
        }
    }
}

impl std::error::Error for AlignmentError {}
