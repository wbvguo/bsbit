//! DNA alphabet primitives.
//!
//! [`Base`] deliberately has only five values. IUPAC ambiguity codes other
//! than `N` are rejected while normalizing external bytes; they are not hidden
//! inside this type as sets of possible bases.

use core::fmt;

/// A normalized DNA base used by the scientific core.
///
/// `N` denotes an unknown observation. It is neither a wildcard nor a fifth
/// nucleotide that compares at zero cost; bisulfite-aware comparison gives any
/// column containing `N` unit cost.
#[repr(transparent)]
#[derive(Clone, Copy, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct Base(u8);

impl Base {
    /// Adenine.
    pub const A: Self = Self(0);
    /// Cytosine.
    pub const C: Self = Self(1);
    /// Guanine.
    pub const G: Self = Self(2);
    /// Thymine.
    pub const T: Self = Self(3);
    /// An unknown base.
    pub const N: Self = Self(4);

    /// Every base in the canonical diagnostic order `A`, `C`, `G`, `T`, `N`.
    pub const ALL: [Self; 5] = [Self::A, Self::C, Self::G, Self::T, Self::N];

    /// The four concrete DNA bases in canonical diagnostic order.
    pub const CANONICAL: [Self; 4] = [Self::A, Self::C, Self::G, Self::T];

    /// Returns this base's uppercase ASCII representation.
    #[must_use]
    pub const fn as_ascii(self) -> u8 {
        match self {
            Self::A => b'A',
            Self::C => b'C',
            Self::G => b'G',
            Self::T => b'T',
            _ => b'N',
        }
    }

    /// Returns the Watson-Crick complement, preserving unknown `N`.
    #[must_use]
    pub const fn complement(self) -> Self {
        match self {
            Self::A => Self::T,
            Self::C => Self::G,
            Self::G => Self::C,
            Self::T => Self::A,
            _ => Self::N,
        }
    }

    /// Returns whether this is the unknown base `N`.
    #[must_use]
    pub const fn is_unknown(self) -> bool {
        self.0 >= Self::N.0
    }

    /// Returns whether this value is one of the five normalized DNA codes.
    #[doc(hidden)]
    #[must_use]
    pub const fn is_normalized(self) -> bool {
        self.0 <= Self::N.0
    }

    /// Returns the stable storage code, mapping invalid backing bytes to `N`.
    #[doc(hidden)]
    #[must_use]
    pub const fn storage_code(self) -> u8 {
        if self.is_normalized() {
            self.0
        } else {
            Self::N.0
        }
    }
}

impl fmt::Debug for Base {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Display::fmt(self, formatter)
    }
}

impl fmt::Display for Base {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match *self {
            Self::A => "A",
            Self::C => "C",
            Self::G => "G",
            Self::T => "T",
            _ => "N",
        })
    }
}
