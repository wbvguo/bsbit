//! Strict eighteen-column extended bedMethyl records.

use core::fmt;
use std::io::{self, Write};

/// Cytosine context encoded in bedMethyl column four.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub enum BedMethylContext {
    /// `CpG` context (`m,CG,0`).
    Cg,
    /// CHG context (`m,CHG,0`).
    Chg,
    /// CHH context (`m,CHH,0`).
    Chh,
}

impl BedMethylContext {
    /// Returns the complete canonical column-four value.
    #[must_use]
    pub const fn modification(self) -> &'static [u8] {
        match self {
            Self::Cg => b"m,CG,0",
            Self::Chg => b"m,CHG,0",
            Self::Chh => b"m,CHH,0",
        }
    }

    fn parse(value: &[u8]) -> Result<Self, BedMethylError> {
        match value {
            b"m,CG,0" => Ok(Self::Cg),
            b"m,CHG,0" => Ok(Self::Chg),
            b"m,CHH,0" => Ok(Self::Chh),
            _ => Err(BedMethylError::InvalidModification),
        }
    }
}

/// Genomic strand encoded in bedMethyl column six.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub enum BedMethylStrand {
    /// Forward (`+`) strand.
    Forward,
    /// Reverse (`-`) strand.
    Reverse,
}

impl BedMethylStrand {
    /// Returns the canonical single-byte field.
    #[must_use]
    pub const fn as_bytes(self) -> &'static [u8] {
        match self {
            Self::Forward => b"+",
            Self::Reverse => b"-",
        }
    }

    fn parse(value: &[u8]) -> Result<Self, BedMethylError> {
        match value {
            b"+" => Ok(Self::Forward),
            b"-" => Ok(Self::Reverse),
            _ => Err(BedMethylError::InvalidStrand),
        }
    }
}

/// Strict extended-bedMethyl syntax or consistency failure.
#[non_exhaustive]
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum BedMethylError {
    /// The row did not contain exactly eighteen tab-separated fields.
    ColumnCount {
        /// Observed field count.
        observed: usize,
    },
    /// Column one was empty.
    EmptyContig,
    /// A numeric column was empty or contained a non-decimal byte.
    InvalidInteger {
        /// One-based column number.
        column: u8,
    },
    /// A numeric column exceeded `u64`.
    IntegerOverflow {
        /// One-based column number.
        column: u8,
    },
    /// Columns two and three did not describe one half-open base.
    InvalidSpan,
    /// Column four was not a supported methyl-cytosine context.
    InvalidModification,
    /// Column six was neither `+` nor `-`.
    InvalidStrand,
    /// Thick coordinates did not equal the site coordinates.
    ThickCoordinatesMismatch,
    /// Column nine was empty.
    EmptyDisplayColor,
    /// Column eleven was not a decimal percentage in `0..=100`.
    InvalidPercent,
    /// Methylated plus unmethylated coverage overflowed `u64`.
    CoverageOverflow,
    /// The two coverage columns or count-derived coverage disagreed.
    CoverageMismatch,
}

impl fmt::Display for BedMethylError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ColumnCount { observed } => write!(
                formatter,
                "expected exactly 18 tab-separated extended bedMethyl columns, observed {observed}"
            ),
            Self::EmptyContig => formatter.write_str("contig in column 1 must not be empty"),
            Self::InvalidInteger { column } => {
                write!(formatter, "column {column} must be a nonnegative integer")
            }
            Self::IntegerOverflow { column } => write!(formatter, "column {column} overflows u64"),
            Self::InvalidSpan => {
                formatter.write_str("columns 2 and 3 must describe one 0-based half-open base")
            }
            Self::InvalidModification => {
                formatter.write_str("column 4 must be `m,CG,0`, `m,CHG,0`, or `m,CHH,0`")
            }
            Self::InvalidStrand => formatter.write_str("column 6 must be `+` or `-`"),
            Self::ThickCoordinatesMismatch => formatter
                .write_str("thick coordinates in columns 7 and 8 must match columns 2 and 3"),
            Self::EmptyDisplayColor => {
                formatter.write_str("display color in column 9 must not be empty")
            }
            Self::InvalidPercent => formatter
                .write_str("percent methylated in column 11 must be a decimal within 0..=100"),
            Self::CoverageOverflow => {
                formatter.write_str("methylated plus unmethylated coverage overflows u64")
            }
            Self::CoverageMismatch => formatter.write_str(
                "coverage columns 5 and 10 must equal methylated plus unmethylated counts",
            ),
        }
    }
}

impl std::error::Error for BedMethylError {}

/// One borrowed, validated eighteen-column extended bedMethyl record.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct BedMethylRecord<'a> {
    contig: &'a [u8],
    start: u64,
    context: BedMethylContext,
    strand: BedMethylStrand,
    display_color: &'a [u8],
    methylated: u64,
    unmethylated: u64,
    other_modification: u64,
    deleted: u64,
    failed: u64,
    different: u64,
    no_call: u64,
}

impl<'a> BedMethylRecord<'a> {
    /// Creates a validated borrowed record from semantic values.
    ///
    /// # Errors
    ///
    /// Returns an empty-field, coordinate-overflow, or coverage-overflow error.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        contig: &'a [u8],
        start: u64,
        context: BedMethylContext,
        strand: BedMethylStrand,
        display_color: &'a [u8],
        methylated: u64,
        unmethylated: u64,
        other_modification: u64,
        deleted: u64,
        failed: u64,
        different: u64,
        no_call: u64,
    ) -> Result<Self, BedMethylError> {
        if contig.is_empty() {
            return Err(BedMethylError::EmptyContig);
        }
        if display_color.is_empty() {
            return Err(BedMethylError::EmptyDisplayColor);
        }
        start.checked_add(1).ok_or(BedMethylError::InvalidSpan)?;
        methylated
            .checked_add(unmethylated)
            .ok_or(BedMethylError::CoverageOverflow)?;
        Ok(Self {
            contig,
            start,
            context,
            strand,
            display_color,
            methylated,
            unmethylated,
            other_modification,
            deleted,
            failed,
            different,
            no_call,
        })
    }

    /// Parses exactly one row without allocating on success.
    ///
    /// # Errors
    ///
    /// Returns the first structural, numeric, vocabulary, or consistency error.
    pub fn parse(line: &'a [u8]) -> Result<Self, BedMethylError> {
        let mut columns = [&[][..]; 18];
        let mut fields = line.split(|byte| *byte == b'\t');
        for (index, column) in columns.iter_mut().enumerate() {
            let Some(value) = fields.next() else {
                return Err(BedMethylError::ColumnCount { observed: index });
            };
            *column = value;
        }
        if fields.next().is_some() {
            return Err(BedMethylError::ColumnCount {
                observed: 19 + fields.count(),
            });
        }
        if columns[0].is_empty() {
            return Err(BedMethylError::EmptyContig);
        }
        let start = parse_u64(columns[1], 2)?;
        let end = parse_u64(columns[2], 3)?;
        if start.checked_add(1) != Some(end) {
            return Err(BedMethylError::InvalidSpan);
        }
        let context = BedMethylContext::parse(columns[3])?;
        let first_coverage = parse_u64(columns[4], 5)?;
        let strand = BedMethylStrand::parse(columns[5])?;
        if parse_u64(columns[6], 7)? != start || parse_u64(columns[7], 8)? != end {
            return Err(BedMethylError::ThickCoordinatesMismatch);
        }
        if columns[8].is_empty() {
            return Err(BedMethylError::EmptyDisplayColor);
        }
        let second_coverage = parse_u64(columns[9], 10)?;
        validate_percent(columns[10])?;
        let methylated = parse_u64(columns[11], 12)?;
        let unmethylated = parse_u64(columns[12], 13)?;
        let other_modification = parse_u64(columns[13], 14)?;
        let deleted = parse_u64(columns[14], 15)?;
        let failed = parse_u64(columns[15], 16)?;
        let different = parse_u64(columns[16], 17)?;
        let no_call = parse_u64(columns[17], 18)?;
        let coverage = methylated
            .checked_add(unmethylated)
            .ok_or(BedMethylError::CoverageOverflow)?;
        if first_coverage != coverage || second_coverage != coverage {
            return Err(BedMethylError::CoverageMismatch);
        }
        Ok(Self {
            contig: columns[0],
            start,
            context,
            strand,
            display_color: columns[8],
            methylated,
            unmethylated,
            other_modification,
            deleted,
            failed,
            different,
            no_call,
        })
    }

    /// Encodes the canonical eighteen-column row with two percent decimals.
    ///
    /// # Errors
    ///
    /// Returns the first error from `writer`.
    pub fn encode<W: Write + ?Sized>(&self, writer: &mut W) -> io::Result<()> {
        let end = self.start + 1;
        let coverage = self.methylated + self.unmethylated;
        let percent = rounded_scaled_ratio(self.methylated, coverage, 10_000);
        writer.write_all(self.contig)?;
        write!(writer, "\t{}\t{end}\t", self.start)?;
        writer.write_all(self.context.modification())?;
        write!(writer, "\t{coverage}\t")?;
        writer.write_all(self.strand.as_bytes())?;
        write!(writer, "\t{}\t{end}\t", self.start)?;
        writer.write_all(self.display_color)?;
        writeln!(
            writer,
            "\t{coverage}\t{}.{:02}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            percent / 100,
            percent % 100,
            self.methylated,
            self.unmethylated,
            self.other_modification,
            self.deleted,
            self.failed,
            self.different,
            self.no_call,
        )
    }

    /// Returns the contig name bytes.
    #[must_use]
    pub const fn contig(self) -> &'a [u8] {
        self.contig
    }
    /// Returns the zero-based site start.
    #[must_use]
    pub const fn start(self) -> u64 {
        self.start
    }
    /// Returns the exclusive site end.
    #[must_use]
    pub const fn end(self) -> u64 {
        self.start + 1
    }
    /// Returns the cytosine context.
    #[must_use]
    pub const fn context(self) -> BedMethylContext {
        self.context
    }
    /// Returns the genomic strand.
    #[must_use]
    pub const fn strand(self) -> BedMethylStrand {
        self.strand
    }
    /// Returns the display-color field.
    #[must_use]
    pub const fn display_color(self) -> &'a [u8] {
        self.display_color
    }
    /// Returns methylated observations.
    #[must_use]
    pub const fn methylated(self) -> u64 {
        self.methylated
    }
    /// Returns unmethylated observations.
    #[must_use]
    pub const fn unmethylated(self) -> u64 {
        self.unmethylated
    }
    /// Returns valid methylated-plus-unmethylated coverage.
    #[must_use]
    pub const fn coverage(self) -> u64 {
        self.methylated + self.unmethylated
    }
    /// Returns other-modification observations.
    #[must_use]
    pub const fn other_modification(self) -> u64 {
        self.other_modification
    }
    /// Returns deletion observations.
    #[must_use]
    pub const fn deleted(self) -> u64 {
        self.deleted
    }
    /// Returns failed observations.
    #[must_use]
    pub const fn failed(self) -> u64 {
        self.failed
    }
    /// Returns different-base observations.
    #[must_use]
    pub const fn different(self) -> u64 {
        self.different
    }
    /// Returns no-call observations.
    #[must_use]
    pub const fn no_call(self) -> u64 {
        self.no_call
    }
}

fn parse_u64(value: &[u8], column: u8) -> Result<u64, BedMethylError> {
    if value.is_empty() {
        return Err(BedMethylError::InvalidInteger { column });
    }
    let mut parsed = 0_u64;
    for &byte in value {
        if !byte.is_ascii_digit() {
            return Err(BedMethylError::InvalidInteger { column });
        }
        parsed = parsed
            .checked_mul(10)
            .and_then(|current| current.checked_add(u64::from(byte - b'0')))
            .ok_or(BedMethylError::IntegerOverflow { column })?;
    }
    Ok(parsed)
}

fn validate_percent(value: &[u8]) -> Result<(), BedMethylError> {
    let mut parts = value.split(|byte| *byte == b'.');
    let whole = parts.next().ok_or(BedMethylError::InvalidPercent)?;
    let whole = parse_u64(whole, 11).map_err(|_| BedMethylError::InvalidPercent)?;
    let fraction = parts.next();
    if parts.next().is_some()
        || fraction
            .is_some_and(|digits| digits.is_empty() || !digits.iter().all(u8::is_ascii_digit))
        || whole > 100
        || (whole == 100
            && fraction.is_some_and(|digits| digits.iter().any(|digit| *digit != b'0')))
    {
        return Err(BedMethylError::InvalidPercent);
    }
    Ok(())
}

fn rounded_scaled_ratio(numerator: u64, denominator: u64, scale: u64) -> u64 {
    if denominator == 0 {
        return 0;
    }
    let denominator = u128::from(denominator);
    let scaled = (u128::from(numerator) * u128::from(scale) + denominator / 2) / denominator;
    u64::try_from(scaled).expect("bounded methylation percentage fits u64")
}
