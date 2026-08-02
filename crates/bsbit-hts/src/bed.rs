//! Borrowed BED3+ interval-line syntax.

use core::fmt;

/// BED3+ line syntax failure.
#[non_exhaustive]
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum BedError {
    /// A data line contained fewer than the three required tab-separated fields.
    ColumnCount {
        /// Observed field count.
        observed: usize,
    },
    /// A coordinate field was empty or contained a non-decimal byte.
    InvalidInteger {
        /// One-based BED column number.
        column: u8,
    },
    /// A coordinate field exceeded `u64`.
    IntegerOverflow {
        /// One-based BED column number.
        column: u8,
    },
}

impl fmt::Display for BedError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ColumnCount { observed } => write!(
                formatter,
                "expected at least 3 tab-separated BED columns, observed {observed}"
            ),
            Self::InvalidInteger { column } => {
                write!(
                    formatter,
                    "BED column {column} must be a nonnegative integer"
                )
            }
            Self::IntegerOverflow { column } => {
                write!(formatter, "BED column {column} overflows u64")
            }
        }
    }
}

impl std::error::Error for BedError {}

/// Borrowed projection of the first three fields in one BED3+ data line.
///
/// This type validates only BED line and decimal syntax. Consumers retain
/// ownership of text encoding, nonempty-span, dictionary, and coordinate-bound
/// policy.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct BedInterval<'a> {
    contig: &'a [u8],
    start: u64,
    end: u64,
}

impl<'a> BedInterval<'a> {
    /// Parses one logical BED3+ line without its line terminator.
    ///
    /// Empty lines, comments, and UCSC `track`/`browser` directives return
    /// `Ok(None)`. Data lines may contain additional columns, which are left
    /// uninterpreted.
    ///
    /// # Errors
    ///
    /// Returns a typed syntax error when a data line lacks a required field or
    /// either coordinate is not a `u64` decimal integer.
    pub fn parse_line(line: &'a [u8]) -> Result<Option<Self>, BedError> {
        if line.is_empty()
            || line.starts_with(b"#")
            || line.starts_with(b"track\t")
            || line.starts_with(b"track ")
            || line.starts_with(b"browser\t")
            || line.starts_with(b"browser ")
        {
            return Ok(None);
        }

        let mut fields = line.split(|byte| *byte == b'\t');
        let contig = fields.next().unwrap_or_default();
        let start = fields.next().ok_or(BedError::ColumnCount { observed: 1 })?;
        let end = fields.next().ok_or(BedError::ColumnCount { observed: 2 })?;
        Ok(Some(Self {
            contig,
            start: parse_u64(start, 2)?,
            end: parse_u64(end, 3)?,
        }))
    }

    /// Returns the borrowed chromosome or reference-sequence field.
    #[must_use]
    pub const fn contig(self) -> &'a [u8] {
        self.contig
    }

    /// Returns the zero-based BED start coordinate.
    #[must_use]
    pub const fn start(self) -> u64 {
        self.start
    }

    /// Returns the zero-based BED end coordinate.
    #[must_use]
    pub const fn end(self) -> u64 {
        self.end
    }
}

fn parse_u64(value: &[u8], column: u8) -> Result<u64, BedError> {
    if value.is_empty() {
        return Err(BedError::InvalidInteger { column });
    }
    let mut parsed = 0_u64;
    for byte in value {
        if !byte.is_ascii_digit() {
            return Err(BedError::InvalidInteger { column });
        }
        parsed = parsed
            .checked_mul(10)
            .and_then(|current| current.checked_add(u64::from(*byte - b'0')))
            .ok_or(BedError::IntegerOverflow { column })?;
    }
    Ok(parsed)
}
