//! Immutable normalized DNA sequences and transformations.

#[cfg(not(any(
    target_pointer_width = "16",
    target_pointer_width = "32",
    target_pointer_width = "64"
)))]
compile_error!("bsbit-core requires a target whose pointer width does not exceed 64 bits");

use core::fmt;
use core::hash::{Hash, Hasher};
use core::ops::Deref;

use crate::alphabet::Base;
use crate::bisulfite::ThreeLetterConversion;

/// A failure to normalize an external DNA byte string.
///
/// Both variants retain the original byte and its zero-based byte offset. The
/// first invalid byte is always reported and no partial sequence is returned.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum NormalizationError {
    /// A recognized IUPAC ambiguity code other than `N` was supplied.
    UnsupportedIupac {
        /// The original byte, preserving its case.
        byte: u8,
        /// The zero-based byte offset in the input.
        offset: u64,
    },
    /// A byte outside the accepted DNA alphabet was supplied.
    InvalidBaseByte {
        /// The original byte.
        byte: u8,
        /// The zero-based byte offset in the input.
        offset: u64,
    },
}

impl NormalizationError {
    /// Returns the original offending byte.
    #[must_use]
    pub const fn byte(self) -> u8 {
        match self {
            Self::UnsupportedIupac { byte, .. } | Self::InvalidBaseByte { byte, .. } => byte,
        }
    }

    /// Returns the zero-based byte offset of the error.
    #[must_use]
    pub const fn offset(self) -> u64 {
        match self {
            Self::UnsupportedIupac { offset, .. } | Self::InvalidBaseByte { offset, .. } => offset,
        }
    }
}

impl fmt::Display for NormalizationError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        let (kind, byte, offset) = match *self {
            Self::UnsupportedIupac { byte, offset } => {
                ("unsupported IUPAC ambiguity code", byte, offset)
            }
            Self::InvalidBaseByte { byte, offset } => ("invalid DNA base byte", byte, offset),
        };

        if byte.is_ascii_graphic() {
            write!(
                formatter,
                "{kind} '{}' at byte offset {offset}",
                char::from(byte)
            )
        } else {
            write!(formatter, "{kind} 0x{byte:02X} at byte offset {offset}")
        }
    }
}

impl std::error::Error for NormalizationError {}

/// An owned, immutable sequence over [`Base`].
///
/// Construction from external text goes through [`normalize_dna`]. Direct
/// construction from `Base` values is infallible because every `Base` is a
/// valid normalized value. Empty sequences are valid at this mathematical
/// layer.
/// An owned immutable sequence over [`Base`].
#[derive(Clone, Debug)]
pub struct NormalizedSequence {
    bases: Box<[Base]>,
    length: u64,
}

impl NormalizedSequence {
    /// Constructs an immutable sequence from already normalized bases.
    #[must_use]
    pub fn from_bases(bases: impl IntoIterator<Item = Base>) -> Self {
        let bases: Box<[Base]> = bases.into_iter().collect();
        let length = storage_len_to_u64(bases.len());
        Self { bases, length }
    }

    /// Returns the normalized bases as an immutable slice.
    #[must_use]
    pub fn bases(&self) -> &[Base] {
        &self.bases
    }

    /// Returns the logical number of bases as an architecture-independent value.
    #[must_use]
    pub const fn len(&self) -> u64 {
        self.length
    }

    /// Returns whether this sequence is empty.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.length == 0
    }

    /// Returns the base at `position`, or `None` when it is out of bounds.
    #[must_use]
    pub fn get(&self, position: u64) -> Option<Base> {
        let Ok(storage_position) = usize::try_from(position) else {
            return None;
        };
        self.bases().get(storage_position).copied()
    }

    /// Iterates over normalized bases by value.
    #[must_use]
    pub fn iter(&self) -> impl ExactSizeIterator<Item = Base> + DoubleEndedIterator + '_ {
        self.bases().iter().copied()
    }

    /// Returns a new reverse-complemented sequence without modifying `self`.
    #[must_use]
    pub fn reverse_complement(&self) -> Self {
        reverse_complement(self)
    }

    /// Returns a new three-letter search projection without modifying `self`.
    #[must_use]
    pub fn three_letter_convert(&self, conversion: ThreeLetterConversion) -> Self {
        three_letter_convert(self, conversion)
    }

    /// Returns a newly allocated uppercase ASCII representation.
    #[must_use]
    pub fn to_ascii(&self) -> Vec<u8> {
        self.iter().map(Base::as_ascii).collect()
    }
}

impl Default for NormalizedSequence {
    fn default() -> Self {
        Self::from_bases([])
    }
}

impl PartialEq for NormalizedSequence {
    fn eq(&self, other: &Self) -> bool {
        self.bases() == other.bases()
    }
}

impl Eq for NormalizedSequence {}

impl PartialOrd for NormalizedSequence {
    fn partial_cmp(&self, other: &Self) -> Option<core::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for NormalizedSequence {
    fn cmp(&self, other: &Self) -> core::cmp::Ordering {
        self.bases().cmp(other.bases())
    }
}

impl Hash for NormalizedSequence {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.bases().hash(state);
    }
}

impl Deref for NormalizedSequence {
    type Target = [Base];

    fn deref(&self) -> &Self::Target {
        self.bases()
    }
}

impl AsRef<[Base]> for NormalizedSequence {
    fn as_ref(&self) -> &[Base] {
        self.bases()
    }
}

impl FromIterator<Base> for NormalizedSequence {
    fn from_iter<T: IntoIterator<Item = Base>>(iter: T) -> Self {
        Self::from_bases(iter)
    }
}

impl From<Vec<Base>> for NormalizedSequence {
    fn from(bases: Vec<Base>) -> Self {
        let length = storage_len_to_u64(bases.len());
        Self {
            bases: bases.into_boxed_slice(),
            length,
        }
    }
}

impl<'a> IntoIterator for &'a NormalizedSequence {
    type Item = Base;
    type IntoIter = core::iter::Copied<core::slice::Iter<'a, Base>>;

    fn into_iter(self) -> Self::IntoIter {
        self.bases().iter().copied()
    }
}

impl fmt::Display for NormalizedSequence {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        for base in self {
            base.fmt(formatter)?;
        }
        Ok(())
    }
}

/// Normalizes an ASCII DNA sequence into an owned immutable value.
///
/// Accepted bytes are uppercase or lowercase `A`, `C`, `G`, `T`, and `N`.
/// Recognized IUPAC ambiguity codes are distinguished from all other invalid
/// bytes. Whitespace is not trimmed because record parsing belongs to the I/O
/// layer.
///
/// # Errors
///
/// Returns the first [`NormalizationError`] from left to right. No partial
/// sequence is returned.
pub fn normalize_dna(input: &[u8]) -> Result<NormalizedSequence, NormalizationError> {
    let mut bases = Vec::with_capacity(input.len());
    for (storage_offset, &byte) in input.iter().enumerate() {
        let offset = storage_len_to_u64(storage_offset);
        let base = match byte {
            b'A' | b'a' => Base::A,
            b'C' | b'c' => Base::C,
            b'G' | b'g' => Base::G,
            b'T' | b't' => Base::T,
            b'N' | b'n' => Base::N,
            b'R' | b'r' | b'Y' | b'y' | b'S' | b's' | b'W' | b'w' | b'K' | b'k' | b'M' | b'm'
            | b'B' | b'b' | b'D' | b'd' | b'H' | b'h' | b'V' | b'v' => {
                return Err(NormalizationError::UnsupportedIupac { byte, offset });
            }
            _ => return Err(NormalizationError::InvalidBaseByte { byte, offset }),
        };
        bases.push(base);
    }
    Ok(bases.into())
}

/// Returns a new sequence in reverse order with every base complemented.
#[must_use]
pub fn reverse_complement(sequence: &NormalizedSequence) -> NormalizedSequence {
    sequence.iter().rev().map(Base::complement).collect()
}

/// Returns a new three-letter search projection of `sequence`.
///
/// `CToT` maps only C to T and `GToA` maps only G to A. This transform is for
/// candidate search; final verification must use the original four-letter
/// sequences and [`crate::bisulfite::classify_bases`].
#[must_use]
pub fn three_letter_convert(
    sequence: &NormalizedSequence,
    conversion: ThreeLetterConversion,
) -> NormalizedSequence {
    sequence
        .iter()
        .map(|base| match (conversion, base) {
            (ThreeLetterConversion::CToT, Base::C) => Base::T,
            (ThreeLetterConversion::GToA, Base::G) => Base::A,
            _ => base,
        })
        .collect()
}

/// Converts a storage index into the public logical width.
///
/// The module-level target guard makes this cast non-narrowing. Keeping the
/// conversion here prevents storage-sized integers from leaking into the
/// scientific API.
const fn storage_len_to_u64(length: usize) -> u64 {
    length as u64
}
