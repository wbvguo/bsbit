//! Shared exact 3' adapter evidence and clipping policy for read-layout output.

use core::fmt;

use bsbit_core::alphabet::Base;

use crate::alignment_policy::{
    DEFAULT_ADAPTER_MAX_CLIP_BASES, DEFAULT_MAX_SOFT_CLIP_BASES, ILLUMINA_UNIVERSAL_ADAPTER,
    MIN_ADAPTER_SUPPORT_BASES,
};
use crate::read_mapping_limits::MAX_READ_BASES;

/// Which endpoint-repair paths an alignment run may use.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum SoftClipMode {
    /// Preserve the mode-qualified behavior: exact adapter repair for both
    /// layouts and candidate-local semi-global completion in sensitive PE.
    #[default]
    Auto,
    /// Emit only whole-read alignments.
    None,
    /// Permit exact adapter-supported clipping but no generic semi-global
    /// endpoint completion.
    Adapter,
}

/// Invalid user-configurable adapter or clipping policy.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum AlignmentOutputPolicyError {
    /// An enabled adapter sequence was empty.
    EmptyAdapter,
    /// An adapter sequence exceeded the mapper's fixed read domain.
    AdapterTooLong {
        /// Observed adapter length.
        length: usize,
        /// Largest adapter supported by the fixed read domain.
        maximum: usize,
    },
    /// An adapter sequence contained a non-ACGTN byte.
    InvalidAdapterBase {
        /// Zero-based byte offset.
        position: usize,
        /// Invalid byte value.
        byte: u8,
    },
    /// The minimum exact overlap was zero or longer than the adapter.
    InvalidAdapterMinimumOverlap {
        /// Requested overlap.
        requested: usize,
        /// Adapter length and therefore largest valid overlap.
        maximum: usize,
    },
    /// An adapter scan bound exceeded the fixed read domain.
    AdapterClipTooLarge {
        /// Requested scan bound.
        requested: usize,
        /// Largest supported scan bound.
        maximum: usize,
    },
    /// A soft-clip bound exceeded the fixed read domain.
    SoftClipTooLarge {
        /// Requested clipping bound.
        requested: usize,
        /// Largest supported clipping bound.
        maximum: usize,
    },
    /// Adapter clipping was enabled but its effective clip domain could not
    /// contain the requested exact overlap.
    AdapterOverlapOutsideClipDomain {
        /// Requested exact adapter overlap.
        overlap: usize,
        /// Effective clipping bound after both limits are applied.
        maximum_clip: usize,
    },
}

impl fmt::Display for AlignmentOutputPolicyError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match *self {
            Self::EmptyAdapter => formatter.write_str("adapter sequence must not be empty"),
            Self::AdapterTooLong { length, maximum } => {
                write!(formatter, "adapter length {length} exceeds {maximum}")
            }
            Self::InvalidAdapterBase { position, byte } => write!(
                formatter,
                "adapter byte 0x{byte:02x} at offset {position} is outside A/C/G/T/N"
            ),
            Self::InvalidAdapterMinimumOverlap { requested, maximum } => write!(
                formatter,
                "adapter minimum overlap {requested} is outside 1..={maximum}"
            ),
            Self::AdapterClipTooLarge { requested, maximum } => write!(
                formatter,
                "adapter maximum clip {requested} exceeds {maximum}"
            ),
            Self::SoftClipTooLarge { requested, maximum } => {
                write!(formatter, "maximum soft clip {requested} exceeds {maximum}")
            }
            Self::AdapterOverlapOutsideClipDomain {
                overlap,
                maximum_clip,
            } => write!(
                formatter,
                "adapter minimum overlap {overlap} exceeds effective clip bound {maximum_clip}"
            ),
        }
    }
}

impl std::error::Error for AlignmentOutputPolicyError {}

/// Validated run-level policy shared by SE and PE output qualification.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct AlignmentOutputPolicy {
    adapter: Option<Box<[u8]>>,
    adapter_minimum_overlap: usize,
    adapter_maximum_clip_bases: usize,
    soft_clip_mode: SoftClipMode,
    maximum_soft_clip_bases: usize,
}

impl Default for AlignmentOutputPolicy {
    fn default() -> Self {
        Self::new(
            Some(ILLUMINA_UNIVERSAL_ADAPTER),
            MIN_ADAPTER_SUPPORT_BASES,
            DEFAULT_ADAPTER_MAX_CLIP_BASES,
            SoftClipMode::Auto,
            DEFAULT_MAX_SOFT_CLIP_BASES,
        )
        .expect("the built-in alignment output policy is valid")
    }
}

impl AlignmentOutputPolicy {
    /// Validates and owns one run-level adapter/clipping policy.
    ///
    /// `None` disables adapter recognition. Adapter bytes are normalized to
    /// uppercase and must belong to the core A/C/G/T/N alphabet.
    ///
    /// # Errors
    ///
    /// Returns [`AlignmentOutputPolicyError`] when a sequence or bound lies
    /// outside the fixed read domain, or when enabled adapter clipping cannot
    /// satisfy its requested exact overlap.
    pub fn new(
        adapter: Option<&[u8]>,
        adapter_minimum_overlap: usize,
        adapter_maximum_clip_bases: usize,
        soft_clip_mode: SoftClipMode,
        maximum_soft_clip_bases: usize,
    ) -> Result<Self, AlignmentOutputPolicyError> {
        if !(1..=MAX_READ_BASES).contains(&adapter_minimum_overlap) {
            return Err(AlignmentOutputPolicyError::InvalidAdapterMinimumOverlap {
                requested: adapter_minimum_overlap,
                maximum: MAX_READ_BASES,
            });
        }
        if adapter_maximum_clip_bases > MAX_READ_BASES {
            return Err(AlignmentOutputPolicyError::AdapterClipTooLarge {
                requested: adapter_maximum_clip_bases,
                maximum: MAX_READ_BASES,
            });
        }
        if maximum_soft_clip_bases > MAX_READ_BASES {
            return Err(AlignmentOutputPolicyError::SoftClipTooLarge {
                requested: maximum_soft_clip_bases,
                maximum: MAX_READ_BASES,
            });
        }
        let adapter = adapter
            .map(|sequence| {
                if sequence.is_empty() {
                    return Err(AlignmentOutputPolicyError::EmptyAdapter);
                }
                if sequence.len() > MAX_READ_BASES {
                    return Err(AlignmentOutputPolicyError::AdapterTooLong {
                        length: sequence.len(),
                        maximum: MAX_READ_BASES,
                    });
                }
                let mut normalized = Vec::with_capacity(sequence.len());
                for (position, &byte) in sequence.iter().enumerate() {
                    let byte = byte.to_ascii_uppercase();
                    if !matches!(byte, b'A' | b'C' | b'G' | b'T' | b'N') {
                        return Err(AlignmentOutputPolicyError::InvalidAdapterBase {
                            position,
                            byte,
                        });
                    }
                    normalized.push(byte);
                }
                if !(1..=normalized.len()).contains(&adapter_minimum_overlap) {
                    return Err(AlignmentOutputPolicyError::InvalidAdapterMinimumOverlap {
                        requested: adapter_minimum_overlap,
                        maximum: normalized.len(),
                    });
                }
                Ok(normalized.into_boxed_slice())
            })
            .transpose()?;
        let effective_adapter_clip = adapter_maximum_clip_bases.min(maximum_soft_clip_bases);
        if adapter.is_some()
            && !matches!(soft_clip_mode, SoftClipMode::None)
            && effective_adapter_clip != 0
            && adapter_minimum_overlap > effective_adapter_clip
        {
            return Err(
                AlignmentOutputPolicyError::AdapterOverlapOutsideClipDomain {
                    overlap: adapter_minimum_overlap,
                    maximum_clip: effective_adapter_clip,
                },
            );
        }
        Ok(Self {
            adapter,
            adapter_minimum_overlap,
            adapter_maximum_clip_bases,
            soft_clip_mode,
            maximum_soft_clip_bases,
        })
    }

    /// Returns the normalized configured adapter, or `None` when disabled.
    #[must_use]
    pub fn adapter(&self) -> Option<&[u8]> {
        self.adapter.as_deref()
    }

    /// Returns the minimum exact adapter overlap.
    #[must_use]
    pub const fn adapter_minimum_overlap(&self) -> usize {
        self.adapter_minimum_overlap
    }

    /// Returns the configured adapter scan bound.
    #[must_use]
    pub const fn adapter_maximum_clip_bases(&self) -> usize {
        self.adapter_maximum_clip_bases
    }

    /// Returns the configured endpoint-repair mode.
    #[must_use]
    pub const fn soft_clip_mode(&self) -> SoftClipMode {
        self.soft_clip_mode
    }

    /// Returns the maximum total soft-clipped query bases admitted per read.
    #[must_use]
    pub const fn maximum_soft_clip_bases(&self) -> usize {
        self.maximum_soft_clip_bases
    }

    #[must_use]
    pub(crate) const fn adapter_clipping_enabled(&self) -> bool {
        self.adapter.is_some()
            && !matches!(self.soft_clip_mode, SoftClipMode::None)
            && self.adapter_maximum_clip_bases != 0
            && self.maximum_soft_clip_bases != 0
    }

    #[must_use]
    pub(crate) const fn semi_global_clipping_enabled(&self) -> bool {
        matches!(self.soft_clip_mode, SoftClipMode::Auto) && self.maximum_soft_clip_bases != 0
    }

    #[must_use]
    const fn effective_adapter_maximum_clip_bases(&self) -> usize {
        if self.adapter_maximum_clip_bases < self.maximum_soft_clip_bases {
            self.adapter_maximum_clip_bases
        } else {
            self.maximum_soft_clip_bases
        }
    }
}

pub(crate) fn sequencing_three_prime_adapter_supported(
    read: &[Base],
    retained_end: usize,
    policy: &AlignmentOutputPolicy,
) -> bool {
    let Some(adapter) = policy.adapter() else {
        return false;
    };
    let clipped = read.get(retained_end..).unwrap_or_default();
    let supported = clipped.len().min(adapter.len());
    supported >= policy.adapter_minimum_overlap()
        && clipped
            .iter()
            .take(supported)
            .zip(adapter.iter().take(supported))
            .all(|(observed, expected)| observed.as_ascii() == *expected)
}

#[must_use]
pub(crate) fn supported_three_prime_adapter_start(
    read: &[Base],
    policy: &AlignmentOutputPolicy,
) -> Option<usize> {
    if !policy.adapter_clipping_enabled() {
        return None;
    }
    let earliest = read
        .len()
        .saturating_sub(policy.effective_adapter_maximum_clip_bases());
    let latest = read.len().checked_sub(policy.adapter_minimum_overlap())?;
    (earliest..=latest).find(|&start| sequencing_three_prime_adapter_supported(read, start, policy))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn policy_normalizes_custom_adapter_and_exposes_independent_bounds() {
        let policy =
            AlignmentOutputPolicy::new(Some(b"acgtnacg"), 5, 17, SoftClipMode::Adapter, 19)
                .expect("valid custom policy");
        assert_eq!(policy.adapter(), Some(b"ACGTNACG".as_slice()));
        assert_eq!(policy.adapter_minimum_overlap(), 5);
        assert_eq!(policy.adapter_maximum_clip_bases(), 17);
        assert_eq!(policy.soft_clip_mode(), SoftClipMode::Adapter);
        assert_eq!(policy.maximum_soft_clip_bases(), 19);
        assert!(policy.adapter_clipping_enabled());
        assert!(!policy.semi_global_clipping_enabled());
    }

    #[test]
    fn zero_clip_bound_cleanly_disables_clipping() {
        let policy = AlignmentOutputPolicy::new(
            Some(ILLUMINA_UNIVERSAL_ADAPTER),
            MIN_ADAPTER_SUPPORT_BASES,
            DEFAULT_ADAPTER_MAX_CLIP_BASES,
            SoftClipMode::Auto,
            0,
        )
        .expect("zero is a valid disabling bound");
        assert!(!policy.adapter_clipping_enabled());
        assert!(!policy.semi_global_clipping_enabled());
    }

    #[test]
    fn policy_rejects_invalid_sequence_overlap_and_bounds() {
        assert!(matches!(
            AlignmentOutputPolicy::new(Some(b"AC-X"), 2, 3, SoftClipMode::Auto, 3),
            Err(AlignmentOutputPolicyError::InvalidAdapterBase { position: 2, .. })
        ));
        assert!(matches!(
            AlignmentOutputPolicy::new(Some(b"ACGT"), 5, 5, SoftClipMode::Auto, 5),
            Err(AlignmentOutputPolicyError::InvalidAdapterMinimumOverlap { .. })
        ));
        assert!(matches!(
            AlignmentOutputPolicy::new(Some(b"ACGTACGT"), 8, 7, SoftClipMode::Adapter, 7,),
            Err(AlignmentOutputPolicyError::AdapterOverlapOutsideClipDomain { .. })
        ));
        assert!(matches!(
            AlignmentOutputPolicy::new(None, 8, MAX_READ_BASES + 1, SoftClipMode::None, 0),
            Err(AlignmentOutputPolicyError::AdapterClipTooLarge { .. })
        ));
    }
}
