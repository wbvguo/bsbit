//! Stable semantic identity for an ordered normalized reference catalog.
//!
//! This module defines only the value and hashing contract shared by index
//! construction, alignment provenance, and callers. It owns no index or file
//! format and performs no I/O.

use core::fmt;

use sha2::{Digest, Sha256};

use crate::alphabet::Base;

/// Domain separator for the frozen reference semantic digest contract.
pub const REFERENCE_SEMANTIC_DIGEST_DOMAIN: &[u8] = b"BSBIT-REFERENCE-SEMANTIC-SHA256-V1\0";

/// Frozen descriptor hashed into every reference semantic digest.
pub const REFERENCE_SEMANTICS_DESCRIPTOR: [u8; 64] = [
    b'B', b'S', b'B', b'S', b'E', b'M', b'0', b'1', // magic
    1, 0, // descriptor major
    0, 0, // descriptor minor
    1, 0, // normalization
    1, 0, // alphabet
    1, 0, // ambiguity
    1, 0, // barriers
    1, 0, // lanes
    1, 0, // FM/search
    1, 0, // coordinates
    1, 0, // catalog/name identity
    64, 0, 0, 0, // descriptor length
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // reserved
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // reserved
];

/// A verified semantic reference digest.
#[derive(Clone, Copy, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct ReferenceSemanticDigest([u8; 32]);

impl ReferenceSemanticDigest {
    /// Constructs a digest from its exact 32 bytes.
    #[must_use]
    pub const fn from_bytes(bytes: [u8; 32]) -> Self {
        Self(bytes)
    }

    /// Returns the exact digest bytes.
    #[must_use]
    pub const fn as_bytes(&self) -> &[u8; 32] {
        &self.0
    }

    /// Consumes the value and returns its exact bytes.
    #[must_use]
    pub const fn into_bytes(self) -> [u8; 32] {
        self.0
    }
}

impl fmt::Debug for ReferenceSemanticDigest {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "ReferenceSemanticDigest({self})")
    }
}

impl fmt::Display for ReferenceSemanticDigest {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        for byte in self.0 {
            write!(formatter, "{byte:02x}")?;
        }
        Ok(())
    }
}

/// A streaming failure while deriving the semantic digest from normalized
/// reference records.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ReferenceSemanticDigestBuildError {
    /// More records were supplied than the declared catalog size.
    TooManyContigs {
        /// Declared number of contigs.
        expected: u64,
    },
    /// The stream ended before the declared catalog size was reached.
    ContigCountMismatch {
        /// Declared number of contigs.
        expected: u64,
        /// Number of records actually supplied.
        observed: u64,
    },
    /// A name or sequence length could not be represented by the contract.
    LengthOverflow,
    /// A contig stream operation was called in the wrong order.
    StreamState {
        /// Stable explanation of the invalid transition.
        reason: &'static str,
    },
    /// Supplied chunks did not equal the declared sequence length.
    SequenceLengthMismatch {
        /// Zero-based contig ordinal.
        contig: u64,
        /// Declared normalized sequence length.
        expected: u64,
        /// Sequence bytes supplied so far.
        observed: u64,
    },
    /// A sequence byte was outside the normalized `A/C/G/T/N` alphabet.
    InvalidBase {
        /// Zero-based contig ordinal.
        contig: u64,
        /// Zero-based byte offset within the contig.
        offset: u64,
        /// Rejected original byte.
        byte: u8,
    },
}

impl fmt::Display for ReferenceSemanticDigestBuildError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match *self {
            Self::TooManyContigs { expected } => {
                write!(
                    formatter,
                    "semantic digest received more than {expected} contigs"
                )
            }
            Self::ContigCountMismatch { expected, observed } => write!(
                formatter,
                "semantic digest expected {expected} contigs, observed {observed}"
            ),
            Self::LengthOverflow => formatter.write_str("semantic digest length exceeds u64"),
            Self::StreamState { reason } => {
                write!(
                    formatter,
                    "semantic digest stream state is invalid: {reason}"
                )
            }
            Self::SequenceLengthMismatch {
                contig,
                expected,
                observed,
            } => write!(
                formatter,
                "semantic digest contig {contig} expected {expected} bases, observed {observed}"
            ),
            Self::InvalidBase {
                contig,
                offset,
                byte,
            } => write!(
                formatter,
                "semantic digest contig {contig} has invalid base 0x{byte:02x} at offset {offset}"
            ),
        }
    }
}

impl std::error::Error for ReferenceSemanticDigestBuildError {}

/// Streaming builder for the exact semantic digest shared by bsbit products.
///
/// Records must be supplied in reference/BAM dictionary order. ASCII input
/// accepts lowercase bases but hashes its uppercase normalized form.
#[derive(Debug)]
pub struct ReferenceSemanticDigestBuilder {
    hasher: Sha256,
    expected_contigs: u64,
    observed_contigs: u64,
    active_contig: Option<ActiveSemanticContig>,
}

#[derive(Clone, Copy, Debug)]
struct ActiveSemanticContig {
    ordinal: u64,
    expected_bases: u64,
    observed_bases: u64,
}

impl ReferenceSemanticDigestBuilder {
    /// Starts a semantic digest for an ordered catalog of `contig_count`
    /// records.
    #[must_use]
    pub fn new(contig_count: u64) -> Self {
        let mut hasher = Sha256::new();
        hasher.update(REFERENCE_SEMANTIC_DIGEST_DOMAIN);
        hasher.update(REFERENCE_SEMANTICS_DESCRIPTOR);
        hasher.update(contig_count.to_le_bytes());
        Self {
            hasher,
            expected_contigs: contig_count,
            observed_contigs: 0,
            active_contig: None,
        }
    }

    /// Starts one contig whose sequence will be supplied in bounded chunks.
    ///
    /// # Errors
    ///
    /// Rejects a nested/excess contig or an unrepresentable name length.
    pub fn begin_ascii_contig(
        &mut self,
        name: &[u8],
        sequence_length: u64,
    ) -> Result<(), ReferenceSemanticDigestBuildError> {
        if self.active_contig.is_some() {
            return Err(ReferenceSemanticDigestBuildError::StreamState {
                reason: "a contig is already active",
            });
        }
        if self.observed_contigs >= self.expected_contigs {
            return Err(ReferenceSemanticDigestBuildError::TooManyContigs {
                expected: self.expected_contigs,
            });
        }
        let name_length = u64::try_from(name.len())
            .map_err(|_| ReferenceSemanticDigestBuildError::LengthOverflow)?;
        self.hasher.update(name_length.to_le_bytes());
        self.hasher.update(sequence_length.to_le_bytes());
        self.hasher.update(name);
        self.active_contig = Some(ActiveSemanticContig {
            ordinal: self.observed_contigs,
            expected_bases: sequence_length,
            observed_bases: 0,
        });
        Ok(())
    }

    /// Hashes the next case-insensitive `A/C/G/T/N` chunk of the active contig.
    ///
    /// # Errors
    ///
    /// Rejects missing setup, invalid bases, overflow, and excess sequence.
    pub fn push_ascii_bases(
        &mut self,
        sequence: &[u8],
    ) -> Result<(), ReferenceSemanticDigestBuildError> {
        let active = self
            .active_contig
            .ok_or(ReferenceSemanticDigestBuildError::StreamState {
                reason: "no contig is active",
            })?;
        let chunk_length = u64::try_from(sequence.len())
            .map_err(|_| ReferenceSemanticDigestBuildError::LengthOverflow)?;
        let final_length = active
            .observed_bases
            .checked_add(chunk_length)
            .ok_or(ReferenceSemanticDigestBuildError::LengthOverflow)?;
        if final_length > active.expected_bases {
            return Err(ReferenceSemanticDigestBuildError::SequenceLengthMismatch {
                contig: active.ordinal,
                expected: active.expected_bases,
                observed: final_length,
            });
        }
        let mut normalized = [0_u8; 8192];
        for (chunk_ordinal, chunk) in sequence.chunks(normalized.len()).enumerate() {
            for (within, (&source, target)) in chunk.iter().zip(normalized.iter_mut()).enumerate() {
                *target = match source {
                    b'A' | b'a' => b'A',
                    b'C' | b'c' => b'C',
                    b'G' | b'g' => b'G',
                    b'T' | b't' => b'T',
                    b'N' | b'n' => b'N',
                    _ => {
                        let within_chunk = chunk_ordinal
                            .checked_mul(normalized.len())
                            .and_then(|value| value.checked_add(within))
                            .ok_or(ReferenceSemanticDigestBuildError::LengthOverflow)?;
                        let offset =
                            active
                                .observed_bases
                                .checked_add(u64::try_from(within_chunk).map_err(|_| {
                                    ReferenceSemanticDigestBuildError::LengthOverflow
                                })?)
                                .ok_or(ReferenceSemanticDigestBuildError::LengthOverflow)?;
                        return Err(ReferenceSemanticDigestBuildError::InvalidBase {
                            contig: active.ordinal,
                            offset,
                            byte: source,
                        });
                    }
                };
            }
            self.hasher.update(&normalized[..chunk.len()]);
        }
        self.active_contig = Some(ActiveSemanticContig {
            observed_bases: final_length,
            ..active
        });
        Ok(())
    }

    /// Finishes the active contig after requiring its declared length.
    ///
    /// # Errors
    ///
    /// Rejects a missing active contig or an incomplete sequence stream.
    pub fn end_ascii_contig(&mut self) -> Result<(), ReferenceSemanticDigestBuildError> {
        let active = self
            .active_contig
            .ok_or(ReferenceSemanticDigestBuildError::StreamState {
                reason: "no contig is active",
            })?;
        if active.observed_bases != active.expected_bases {
            return Err(ReferenceSemanticDigestBuildError::SequenceLengthMismatch {
                contig: active.ordinal,
                expected: active.expected_bases,
                observed: active.observed_bases,
            });
        }
        self.active_contig = None;
        self.observed_contigs += 1;
        Ok(())
    }

    /// Adds one named ASCII sequence in catalog order.
    ///
    /// # Errors
    ///
    /// Rejects excess records, unrepresentable lengths, and bytes outside
    /// case-insensitive `A/C/G/T/N`.
    pub fn push_ascii_contig(
        &mut self,
        name: &[u8],
        sequence: &[u8],
    ) -> Result<(), ReferenceSemanticDigestBuildError> {
        let sequence_length = u64::try_from(sequence.len())
            .map_err(|_| ReferenceSemanticDigestBuildError::LengthOverflow)?;
        self.begin_ascii_contig(name, sequence_length)?;
        self.push_ascii_bases(sequence)?;
        self.end_ascii_contig()
    }

    /// Adds one named sequence already represented by normalized core bases.
    ///
    /// # Errors
    ///
    /// Rejects excess records or unrepresentable lengths.
    pub fn push_normalized_contig(
        &mut self,
        name: &[u8],
        sequence: &[Base],
    ) -> Result<(), ReferenceSemanticDigestBuildError> {
        let sequence_length = u64::try_from(sequence.len())
            .map_err(|_| ReferenceSemanticDigestBuildError::LengthOverflow)?;
        self.begin_ascii_contig(name, sequence_length)?;
        let mut ascii = [0_u8; 8192];
        for chunk in sequence.chunks(ascii.len()) {
            for (&base, target) in chunk.iter().zip(ascii.iter_mut()) {
                *target = base.as_ascii();
            }
            self.push_ascii_bases(&ascii[..chunk.len()])?;
        }
        self.end_ascii_contig()
    }

    /// Finishes the digest after requiring the declared record count.
    ///
    /// # Errors
    ///
    /// Rejects an incomplete record stream.
    pub fn finish(self) -> Result<ReferenceSemanticDigest, ReferenceSemanticDigestBuildError> {
        if self.active_contig.is_some() {
            return Err(ReferenceSemanticDigestBuildError::StreamState {
                reason: "the final contig was not ended",
            });
        }
        if self.observed_contigs != self.expected_contigs {
            return Err(ReferenceSemanticDigestBuildError::ContigCountMismatch {
                expected: self.expected_contigs,
                observed: self.observed_contigs,
            });
        }
        Ok(ReferenceSemanticDigest::from_bytes(
            self.hasher.finalize().into(),
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ascii_stream_is_case_insensitive_and_chunk_stable() {
        let mut whole = ReferenceSemanticDigestBuilder::new(1);
        whole.push_ascii_contig(b"chr1", b"ACGTN").unwrap();

        let mut chunks = ReferenceSemanticDigestBuilder::new(1);
        chunks.begin_ascii_contig(b"chr1", 5).unwrap();
        chunks.push_ascii_bases(b"ac").unwrap();
        chunks.push_ascii_bases(b"gtn").unwrap();
        chunks.end_ascii_contig().unwrap();

        assert_eq!(whole.finish().unwrap(), chunks.finish().unwrap());
    }
}
