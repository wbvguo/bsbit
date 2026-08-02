//! Typed core CIGAR construction, validation, parsing, and evaluation.
//!
//! Level 1 deliberately supports only `M`, `I`, and `D`. A validated
//! [`CoreCigar`] has positive run lengths and no adjacent equal operations.
//! Construction and replay are deterministic, use checked arithmetic, and do
//! not mutate either input sequence.

use core::fmt;
use core::str::FromStr;

/// A Level 1 core CIGAR operation.
///
/// The declaration order is the canonical lexicographic order `D < I < M`.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub enum CoreCigarOp {
    /// Consume one reference base and no query base.
    D,
    /// Consume one query base and no reference base.
    I,
    /// Consume one reference base and one query base.
    M,
}

impl CoreCigarOp {
    /// All operations in the canonical lexicographic order `D < I < M`.
    pub const LEXICOGRAPHIC: [Self; 3] = [Self::D, Self::I, Self::M];

    /// Returns whether this operation consumes reference bases.
    #[must_use]
    pub const fn consumes_reference(self) -> bool {
        matches!(self, Self::D | Self::M)
    }

    /// Returns whether this operation consumes query bases.
    #[must_use]
    pub const fn consumes_query(self) -> bool {
        matches!(self, Self::I | Self::M)
    }

    /// Returns whether this operation is an insertion or deletion.
    #[must_use]
    pub const fn is_gap(self) -> bool {
        matches!(self, Self::D | Self::I)
    }

    const fn as_char(self) -> char {
        match self {
            Self::D => 'D',
            Self::I => 'I',
            Self::M => 'M',
        }
    }

    const fn from_byte(byte: u8) -> Option<Self> {
        match byte {
            b'D' => Some(Self::D),
            b'I' => Some(Self::I),
            b'M' => Some(Self::M),
            _ => None,
        }
    }
}

impl fmt::Display for CoreCigarOp {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_fmt(format_args!("{}", self.as_char()))
    }
}

/// One untrusted run at a CIGAR construction or decoding boundary.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub struct RawCigarRun {
    operation: CoreCigarOp,
    length: u64,
}

impl RawCigarRun {
    /// Creates a raw run. Zero length is retained for subsequent validation.
    #[must_use]
    pub const fn new(operation: CoreCigarOp, length: u64) -> Self {
        Self { operation, length }
    }

    /// Returns this run's operation.
    #[must_use]
    pub const fn operation(self) -> CoreCigarOp {
        self.operation
    }

    /// Returns this untrusted run's length.
    #[must_use]
    pub const fn length(self) -> u64 {
        self.length
    }
}

/// An ordered, potentially noncanonical core CIGAR.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct RawCoreCigar {
    runs: Vec<RawCigarRun>,
}

impl RawCoreCigar {
    /// Collects raw runs without validating or normalizing them.
    #[must_use]
    pub fn new(runs: impl IntoIterator<Item = RawCigarRun>) -> Self {
        Self {
            runs: runs.into_iter().collect(),
        }
    }

    /// Returns the raw runs in forward-reference order.
    #[must_use]
    pub fn runs(&self) -> &[RawCigarRun] {
        &self.runs
    }
}

impl From<Vec<RawCigarRun>> for RawCoreCigar {
    fn from(runs: Vec<RawCigarRun>) -> Self {
        Self { runs }
    }
}

/// One validated, positive-length core CIGAR run.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub struct CoreCigarRun {
    operation: CoreCigarOp,
    length: u64,
}

impl CoreCigarRun {
    fn positive(operation: CoreCigarOp, length: u64) -> Self {
        debug_assert!(length > 0);
        Self { operation, length }
    }

    /// Returns this run's operation.
    #[must_use]
    pub const fn operation(self) -> CoreCigarOp {
        self.operation
    }

    /// Returns this run's positive length.
    #[must_use]
    pub const fn length(self) -> u64 {
        self.length
    }
}

/// A canonical run-length encoded Level 1 CIGAR.
///
/// Every run is positive and adjacent runs always have different operations.
/// Expected sequence lengths are checked separately by [`validate_cigar`].
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct CoreCigar {
    runs: Vec<CoreCigarRun>,
}

impl CoreCigar {
    /// Constructs the unique ungapped all-match path for a nonempty interval.
    #[doc(hidden)]
    #[must_use]
    pub fn all_matches(length: u64) -> Self {
        debug_assert!(length > 0);
        Self {
            runs: vec![CoreCigarRun::positive(CoreCigarOp::M, length)],
        }
    }

    /// Returns the canonical runs in forward-reference order.
    #[must_use]
    pub fn runs(&self) -> &[CoreCigarRun] {
        &self.runs
    }

    /// Returns the number of coalesced runs.
    #[must_use]
    pub fn run_count(&self) -> usize {
        self.runs.len()
    }

    /// Returns whether this CIGAR has no operations.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.runs.is_empty()
    }
}

impl fmt::Display for CoreCigar {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        for run in &self.runs {
            write!(formatter, "{}{}", run.length, run.operation)?;
        }
        Ok(())
    }
}

impl FromStr for CoreCigar {
    type Err = CigarError;

    fn from_str(input: &str) -> Result<Self, Self::Err> {
        parse_core_cigar(input)
    }
}

/// The sequence domain involved in CIGAR consumption arithmetic.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum CigarDomain {
    /// Reference-sequence consumption.
    Reference,
    /// Oriented-query-sequence consumption.
    Query,
}

/// A structured CIGAR construction, parsing, or validation error.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum CigarError {
    /// A raw run had length zero.
    ZeroLengthCigarRun {
        /// Zero-based run index.
        run_index: usize,
        /// Operation carried by the run.
        operation: CoreCigarOp,
    },
    /// Two adjacent raw runs had the same operation in a strict input.
    NonCanonicalCigarRuns {
        /// Index of the preceding run.
        previous_run_index: usize,
        /// Index of the repeated run.
        run_index: usize,
        /// Repeated operation.
        operation: CoreCigarOp,
    },
    /// Consumption arithmetic exceeded `u64`.
    CigarConsumptionOverflow {
        /// Zero-based run index.
        run_index: usize,
        /// Run operation.
        operation: CoreCigarOp,
        /// Sequence domain whose count overflowed.
        domain: CigarDomain,
        /// Count before adding this run.
        accumulated: u64,
        /// Run length being added.
        run_length: u64,
    },
    /// Final reference or query consumption did not equal the expected length.
    CigarLengthMismatch {
        /// Expected reference length.
        expected_reference: u64,
        /// Observed reference consumption.
        observed_reference: u64,
        /// Expected query length.
        expected_query: u64,
        /// Observed query consumption.
        observed_query: u64,
    },
    /// A decimal run length exceeded `u64` while parsing.
    CigarRunLengthOverflow {
        /// Zero-based run index.
        run_index: usize,
        /// Zero-based byte offset of the digit that overflowed.
        byte_offset: usize,
    },
    /// A decimal run length used a noncanonical leading zero.
    NonCanonicalCigarRunLength {
        /// Zero-based run index.
        run_index: usize,
        /// Zero-based byte offset of the leading zero.
        byte_offset: usize,
    },
    /// A parsed run did not begin with an ASCII decimal digit.
    ExpectedCigarRunLength {
        /// Zero-based byte offset.
        byte_offset: usize,
        /// Byte found at that offset.
        found: u8,
    },
    /// A decimal run length reached end-of-input without an operation.
    MissingCigarOperation {
        /// Zero-based run index.
        run_index: usize,
        /// End-of-input byte offset.
        byte_offset: usize,
    },
    /// A run used an operation outside `M`, `I`, and `D`.
    UnknownCigarOperation {
        /// Zero-based run index.
        run_index: usize,
        /// Zero-based byte offset of the operation.
        byte_offset: usize,
        /// Unsupported byte.
        found: u8,
    },
}

impl fmt::Display for CigarError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ZeroLengthCigarRun {
                run_index,
                operation,
            } => write!(
                formatter,
                "CIGAR run {run_index} has zero length for operation {operation}"
            ),
            Self::NonCanonicalCigarRuns {
                previous_run_index,
                run_index,
                operation,
            } => write!(
                formatter,
                "CIGAR runs {previous_run_index} and {run_index} repeat operation {operation}"
            ),
            Self::CigarConsumptionOverflow {
                run_index,
                operation,
                domain,
                accumulated,
                run_length,
            } => write!(
                formatter,
                "{domain:?} consumption overflow at CIGAR run {run_index} ({operation}): {accumulated} + {run_length}"
            ),
            Self::CigarLengthMismatch {
                expected_reference,
                observed_reference,
                expected_query,
                observed_query,
            } => write!(
                formatter,
                "CIGAR consumes reference/query {observed_reference}/{observed_query}, expected {expected_reference}/{expected_query}"
            ),
            Self::CigarRunLengthOverflow {
                run_index,
                byte_offset,
            } => write!(
                formatter,
                "CIGAR run {run_index} length overflows at byte offset {byte_offset}"
            ),
            Self::NonCanonicalCigarRunLength {
                run_index,
                byte_offset,
            } => write!(
                formatter,
                "CIGAR run {run_index} has a noncanonical leading zero at byte offset {byte_offset}"
            ),
            Self::ExpectedCigarRunLength { byte_offset, found } => write!(
                formatter,
                "expected CIGAR run length at byte offset {byte_offset}, found 0x{found:02X}"
            ),
            Self::MissingCigarOperation {
                run_index,
                byte_offset,
            } => write!(
                formatter,
                "CIGAR run {run_index} is missing an operation at byte offset {byte_offset}"
            ),
            Self::UnknownCigarOperation {
                run_index,
                byte_offset,
                found,
            } => write!(
                formatter,
                "unknown CIGAR operation 0x{found:02X} for run {run_index} at byte offset {byte_offset}"
            ),
        }
    }
}

impl std::error::Error for CigarError {}

/// Checked reference and query consumption derived from a CIGAR.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct CigarConsumption {
    reference: u64,
    query: u64,
}

impl CigarConsumption {
    /// Returns total consumed reference bases.
    #[must_use]
    pub const fn reference(self) -> u64 {
        self.reference
    }

    /// Returns total consumed oriented-query bases.
    #[must_use]
    pub const fn query(self) -> u64 {
        self.query
    }
}

/// Coalesces adjacent equal raw runs while retaining every consumed base.
///
/// Empty input produces an empty [`CoreCigar`]. Zero raw runs and checked
/// coalescing overflow are errors. This function does not validate sequence
/// consumption and does not choose among alternative alignment paths.
///
/// # Errors
///
/// Returns [`CigarError::ZeroLengthCigarRun`] or
/// [`CigarError::CigarConsumptionOverflow`] at the first offending run.
pub fn canonicalize_cigar(raw: &RawCoreCigar) -> Result<CoreCigar, CigarError> {
    let mut runs: Vec<CoreCigarRun> = Vec::with_capacity(raw.runs.len());
    for (run_index, raw_run) in raw.runs.iter().copied().enumerate() {
        if raw_run.length == 0 {
            return Err(CigarError::ZeroLengthCigarRun {
                run_index,
                operation: raw_run.operation,
            });
        }

        if let Some(last) = runs.last_mut()
            && last.operation == raw_run.operation
        {
            last.length = last.length.checked_add(raw_run.length).ok_or(
                CigarError::CigarConsumptionOverflow {
                    run_index,
                    operation: raw_run.operation,
                    domain: coalescing_domain(raw_run.operation),
                    accumulated: last.length,
                    run_length: raw_run.length,
                },
            )?;
            continue;
        }

        runs.push(CoreCigarRun::positive(raw_run.operation, raw_run.length));
    }
    Ok(CoreCigar { runs })
}

/// Canonicalizes an expanded forward operation stream.
///
/// Every item denotes one consumed column/base. The output is equivalent to
/// converting items to length-one raw runs and calling [`canonicalize_cigar`].
///
/// # Errors
///
/// Returns a checked overflow if one coalesced run would exceed `u64`.
pub fn canonicalize_operations(
    operations: impl IntoIterator<Item = CoreCigarOp>,
) -> Result<CoreCigar, CigarError> {
    let iterator = operations.into_iter();
    // `IntoIterator` is an untrusted public boundary. Its `size_hint` may be
    // arbitrarily large or dishonest, so it must not drive an eager allocation.
    let mut runs: Vec<CoreCigarRun> = Vec::new();
    for (operation_index, operation) in iterator.enumerate() {
        if let Some(last) = runs.last_mut()
            && last.operation == operation
        {
            last.length =
                last.length
                    .checked_add(1)
                    .ok_or(CigarError::CigarConsumptionOverflow {
                        run_index: operation_index,
                        operation,
                        domain: coalescing_domain(operation),
                        accumulated: last.length,
                        run_length: 1,
                    })?;
        } else {
            runs.push(CoreCigarRun::positive(operation, 1));
        }
    }
    Ok(CoreCigar { runs })
}

/// Strictly constructs a validated CIGAR for expected sequence lengths.
///
/// Each raw run is checked in this order: zero length, adjacent equality, then
/// checked reference/query consumption. Exact final consumption is checked
/// after the last run. No invalid raw state is returned as a [`CoreCigar`].
///
/// # Errors
///
/// Returns the first structured [`CigarError`] in the specified order.
pub fn try_core_cigar(
    raw: &RawCoreCigar,
    expected_reference: u64,
    expected_query: u64,
) -> Result<CoreCigar, CigarError> {
    let mut runs = Vec::with_capacity(raw.runs.len());
    let mut reference = 0_u64;
    let mut query = 0_u64;
    let mut previous = None;

    for (run_index, raw_run) in raw.runs.iter().copied().enumerate() {
        if raw_run.length == 0 {
            return Err(CigarError::ZeroLengthCigarRun {
                run_index,
                operation: raw_run.operation,
            });
        }
        if previous == Some(raw_run.operation) {
            return Err(CigarError::NonCanonicalCigarRuns {
                previous_run_index: run_index - 1,
                run_index,
                operation: raw_run.operation,
            });
        }

        add_consumption(
            &mut reference,
            &mut query,
            run_index,
            raw_run.operation,
            raw_run.length,
        )?;
        runs.push(CoreCigarRun::positive(raw_run.operation, raw_run.length));
        previous = Some(raw_run.operation);
    }

    exact_consumption(reference, query, expected_reference, expected_query)?;
    Ok(CoreCigar { runs })
}

/// Checks a structurally valid CIGAR against expected sequence lengths.
///
/// # Errors
///
/// Returns checked consumption overflow or exact-length mismatch.
pub fn validate_cigar(
    cigar: &CoreCigar,
    expected_reference: u64,
    expected_query: u64,
) -> Result<CigarConsumption, CigarError> {
    let consumption = cigar_consumption(cigar)?;
    exact_consumption(
        consumption.reference,
        consumption.query,
        expected_reference,
        expected_query,
    )?;
    Ok(consumption)
}

/// Parses the stable logical form of a canonical Level 1 CIGAR.
///
/// The empty string represents the empty CIGAR. Nonempty input must consist
/// entirely of positive decimal lengths followed by `M`, `I`, or `D`.
/// Adjacent equal operations, length overflow, unknown operations, and trailing
/// incomplete/garbage data are rejected. Sequence consumption is not embedded
/// in this logical representation and is checked separately.
///
/// # Errors
///
/// Returns the first structured syntax or canonicality error.
pub fn parse_core_cigar(input: &str) -> Result<CoreCigar, CigarError> {
    let bytes = input.as_bytes();
    let mut offset = 0_usize;
    let mut run_index = 0_usize;
    let mut previous = None;
    let mut runs = Vec::new();

    while offset < bytes.len() {
        if !bytes[offset].is_ascii_digit() {
            return Err(CigarError::ExpectedCigarRunLength {
                byte_offset: offset,
                found: bytes[offset],
            });
        }

        if bytes[offset] == b'0' && bytes.get(offset + 1).is_some_and(u8::is_ascii_digit) {
            return Err(CigarError::NonCanonicalCigarRunLength {
                run_index,
                byte_offset: offset,
            });
        }

        let mut length = 0_u64;
        while offset < bytes.len() && bytes[offset].is_ascii_digit() {
            length = length
                .checked_mul(10)
                .and_then(|value| value.checked_add(u64::from(bytes[offset] - b'0')))
                .ok_or(CigarError::CigarRunLengthOverflow {
                    run_index,
                    byte_offset: offset,
                })?;
            offset += 1;
        }

        if offset == bytes.len() {
            return Err(CigarError::MissingCigarOperation {
                run_index,
                byte_offset: offset,
            });
        }
        let operation =
            CoreCigarOp::from_byte(bytes[offset]).ok_or(CigarError::UnknownCigarOperation {
                run_index,
                byte_offset: offset,
                found: bytes[offset],
            })?;
        offset += 1;

        if length == 0 {
            return Err(CigarError::ZeroLengthCigarRun {
                run_index,
                operation,
            });
        }
        if previous == Some(operation) {
            return Err(CigarError::NonCanonicalCigarRuns {
                previous_run_index: run_index - 1,
                run_index,
                operation,
            });
        }

        runs.push(CoreCigarRun::positive(operation, length));
        previous = Some(operation);
        run_index += 1;
    }

    Ok(CoreCigar { runs })
}

fn coalescing_domain(operation: CoreCigarOp) -> CigarDomain {
    if operation == CoreCigarOp::I {
        CigarDomain::Query
    } else {
        CigarDomain::Reference
    }
}

fn cigar_consumption(cigar: &CoreCigar) -> Result<CigarConsumption, CigarError> {
    let mut reference = 0_u64;
    let mut query = 0_u64;
    for (run_index, run) in cigar.runs.iter().copied().enumerate() {
        add_consumption(
            &mut reference,
            &mut query,
            run_index,
            run.operation,
            run.length,
        )?;
    }
    Ok(CigarConsumption { reference, query })
}

fn add_consumption(
    reference: &mut u64,
    query: &mut u64,
    run_index: usize,
    operation: CoreCigarOp,
    run_length: u64,
) -> Result<(), CigarError> {
    if operation.consumes_reference() {
        *reference =
            reference
                .checked_add(run_length)
                .ok_or(CigarError::CigarConsumptionOverflow {
                    run_index,
                    operation,
                    domain: CigarDomain::Reference,
                    accumulated: *reference,
                    run_length,
                })?;
    }
    if operation.consumes_query() {
        *query = query
            .checked_add(run_length)
            .ok_or(CigarError::CigarConsumptionOverflow {
                run_index,
                operation,
                domain: CigarDomain::Query,
                accumulated: *query,
                run_length,
            })?;
    }
    Ok(())
}

fn exact_consumption(
    observed_reference: u64,
    observed_query: u64,
    expected_reference: u64,
    expected_query: u64,
) -> Result<(), CigarError> {
    if observed_reference == expected_reference && observed_query == expected_query {
        Ok(())
    } else {
        Err(CigarError::CigarLengthMismatch {
            expected_reference,
            observed_reference,
            expected_query,
            observed_query,
        })
    }
}
