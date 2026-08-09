//! Exact in-memory FM-index reference backend.
//!
//! This module owns a full suffix array, an explicit BWT, and a prefix rank
//! table. It is the safe scalar oracle for later compressed or sampled
//! backends. Input is already projected and contains only canonical bases.

use core::fmt;
use core::mem::size_of;

use bsbit_core::alphabet::Base;

const SENTINEL_CODE: u8 = 0;

/// A canonical base accepted by exact FM search.
///
/// Unknown bases, separators, and the terminal sentinel are deliberately not
/// representable.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
#[repr(u8)]
pub enum SearchBase {
    /// Adenine.
    A = b'A',
    /// Cytosine.
    C = b'C',
    /// Guanine.
    G = b'G',
    /// Thymine.
    T = b'T',
}

/// A projected query symbol in the persisted three-letter order.
///
/// This is intentionally separate from [`SearchBase`]: its representation is
/// the combined-index digit, so mapping searches can
/// consume projected reads without translating every symbol in their hot
/// loops.
#[doc(hidden)]
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
#[repr(u8)]
pub enum ProjectedBase {
    /// Guanine, digit zero in the frozen index.
    G = 0,
    /// Thymine, digit one in the frozen index.
    T = 1,
    /// Adenine, digit two in the frozen index.
    A = 2,
}

impl ProjectedBase {
    /// Returns the digit consumed directly by the frozen combined index.
    #[must_use]
    pub const fn digit(self) -> u8 {
        self as u8
    }
}

impl SearchBase {
    /// Every searchable base in canonical order.
    pub const ALL: [Self; 4] = [Self::A, Self::C, Self::G, Self::T];

    /// Converts a normalized biological base when it is canonical.
    #[must_use]
    pub const fn from_base(base: Base) -> Option<Self> {
        match base {
            Base::A => Some(Self::A),
            Base::C => Some(Self::C),
            Base::G => Some(Self::G),
            Base::T => Some(Self::T),
            _ => None,
        }
    }

    /// Converts this exact-search symbol to the general normalized alphabet.
    #[must_use]
    pub const fn as_base(self) -> Base {
        match self {
            Self::A => Base::A,
            Self::C => Base::C,
            Self::G => Base::G,
            Self::T => Base::T,
        }
    }

    /// Returns the uppercase ASCII byte for diagnostics.
    #[must_use]
    pub const fn as_ascii(self) -> u8 {
        self.as_base().as_ascii()
    }

    const fn rank_index(self) -> usize {
        match self {
            Self::A => 0,
            Self::C => 1,
            Self::G => 2,
            Self::T => 3,
        }
    }

    pub(crate) const fn bwt_code(self) -> u8 {
        match self {
            Self::A => 1,
            Self::C => 2,
            Self::G => 3,
            Self::T => 4,
        }
    }
}

impl fmt::Display for SearchBase {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        char::from(self.as_ascii()).fmt(formatter)
    }
}

/// Maximum permitted number of suffix rows, including the terminal suffix.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct FmBuildLimit(u64);

impl FmBuildLimit {
    /// A limit that admits every representable suffix count.
    pub const MAX: Self = Self(u64::MAX);

    /// Constructs an explicit suffix-row limit.
    ///
    /// Zero rejects every index because even an empty text has one terminal
    /// suffix.
    #[must_use]
    pub const fn new(max_suffix_count: u64) -> Self {
        Self(max_suffix_count)
    }

    /// Returns the maximum suffix count.
    #[must_use]
    pub const fn get(self) -> u64 {
        self.0
    }
}

/// A diagnostic allocation site in the current reference backend.
///
/// These variants explain failures; they are not a stable storage-layout
/// contract for later backends.
#[non_exhaustive]
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum FmAllocation {
    /// Full suffix-array rows.
    SuffixArray,
    /// Current prefix-doubling equivalence classes.
    RankClasses,
    /// Next prefix-doubling equivalence classes.
    NextRankClasses,
    /// Explicit BWT symbols.
    Bwt,
    /// Prefix occurrence-count rows.
    RankPrefixes,
    /// Owned locate output.
    LocateResults,
}

/// A checked FM-index construction or query error.
#[non_exhaustive]
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum FmError {
    /// A physical input length cannot be represented by the logical width.
    TextLengthNotRepresentable {
        /// Physical input length.
        text_len: usize,
    },
    /// Adding the unique terminal suffix overflowed.
    SuffixCountOverflow {
        /// Logical text length before adding the terminal suffix.
        text_len: u64,
    },
    /// Adding the final rank boundary row overflowed.
    RankRowCountOverflow {
        /// Logical suffix count before adding the boundary row.
        suffix_count: u64,
    },
    /// The explicit construction limit rejects the requested suffix count.
    BuildLimitExceeded {
        /// Requested suffix rows, including the terminal suffix.
        suffix_count: u64,
        /// Configured maximum suffix rows.
        max_suffix_count: u64,
    },
    /// A component cannot fit the architecture or allocation byte limit.
    AllocationSizeOverflow {
        /// Component being sized.
        component: FmAllocation,
        /// Requested element count.
        elements: u64,
        /// Size of one element in bytes.
        element_size: u64,
    },
    /// Fallible reservation failed after size validation.
    AllocationFailed {
        /// Component being reserved.
        component: FmAllocation,
        /// Requested element count.
        elements: u64,
    },
    /// A rank boundary lies outside the BWT boundary domain.
    RankBoundaryOutOfBounds {
        /// Requested prefix boundary.
        boundary: u64,
        /// Number of BWT rows; the same value is the largest legal boundary.
        suffix_count: u64,
    },
    /// A half-open interval has its upper bound before its lower bound.
    InvertedInterval {
        /// Requested inclusive lower row.
        lower: u64,
        /// Requested exclusive upper row.
        upper: u64,
    },
    /// An interval was created for a different suffix-count domain.
    IntervalDomainMismatch {
        /// Suffix count recorded by the interval.
        interval_suffix_count: u64,
        /// Suffix count of the receiving index.
        index_suffix_count: u64,
    },
    /// A half-open interval lies outside the receiving index.
    IntervalOutOfBounds {
        /// Requested inclusive lower row.
        lower: u64,
        /// Requested exclusive upper row.
        upper: u64,
        /// Number of suffix rows in the receiving index.
        suffix_count: u64,
    },
}

impl fmt::Display for FmError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::TextLengthNotRepresentable { text_len } => {
                write!(
                    formatter,
                    "text length {text_len} is not representable as u64"
                )
            }
            Self::SuffixCountOverflow { text_len } => {
                write!(
                    formatter,
                    "text length {text_len} cannot include a terminal suffix"
                )
            }
            Self::RankRowCountOverflow { suffix_count } => write!(
                formatter,
                "suffix count {suffix_count} cannot include the final rank boundary"
            ),
            Self::BuildLimitExceeded {
                suffix_count,
                max_suffix_count,
            } => write!(
                formatter,
                "requested {suffix_count} suffix rows exceeds build limit {max_suffix_count}"
            ),
            Self::AllocationSizeOverflow {
                component,
                elements,
                element_size,
            } => write!(
                formatter,
                "cannot size {component:?}: {elements} elements of {element_size} bytes"
            ),
            Self::AllocationFailed {
                component,
                elements,
            } => write!(
                formatter,
                "failed to reserve {elements} elements for {component:?}"
            ),
            Self::RankBoundaryOutOfBounds {
                boundary,
                suffix_count,
            } => write!(
                formatter,
                "rank boundary {boundary} exceeds BWT length {suffix_count}"
            ),
            Self::InvertedInterval { lower, upper } => {
                write!(formatter, "FM interval [{lower}, {upper}) is inverted")
            }
            Self::IntervalDomainMismatch {
                interval_suffix_count,
                index_suffix_count,
            } => write!(
                formatter,
                "FM interval domain has {interval_suffix_count} suffix rows; index has {index_suffix_count}"
            ),
            Self::IntervalOutOfBounds {
                lower,
                upper,
                suffix_count,
            } => write!(
                formatter,
                "FM interval [{lower}, {upper}) exceeds suffix count {suffix_count}"
            ),
        }
    }
}

impl std::error::Error for FmError {}

/// A checked half-open suffix-array row interval.
///
/// Empty intervals preserve their exact insertion-point coordinate. The private
/// suffix-count tag rejects use with a differently sized index; full reference
/// identity is a later projected-reference responsibility.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub struct FmInterval {
    lower: u64,
    upper: u64,
    suffix_count: u64,
}

impl FmInterval {
    const fn new_unchecked(lower: u64, upper: u64, suffix_count: u64) -> Self {
        Self {
            lower,
            upper,
            suffix_count,
        }
    }

    /// Constructs an interval for a checked private FM backend.
    ///
    /// This workspace-internal boundary lets an independently implemented
    /// packed backend reuse the owner-bound reference search contract without
    /// exposing unchecked interval construction to ordinary callers.
    ///
    /// # Errors
    ///
    /// Rejects inverted bounds or rows outside `suffix_count`.
    #[doc(hidden)]
    pub fn private_checked(lower: u64, upper: u64, suffix_count: u64) -> Result<Self, FmError> {
        if upper < lower {
            return Err(FmError::InvertedInterval { lower, upper });
        }
        if upper > suffix_count {
            return Err(FmError::IntervalOutOfBounds {
                lower,
                upper,
                suffix_count,
            });
        }
        Ok(Self::new_unchecked(lower, upper, suffix_count))
    }

    /// Returns the inclusive lower suffix-array row.
    #[must_use]
    pub const fn lower(self) -> u64 {
        self.lower
    }

    /// Returns the exclusive upper suffix-array row.
    #[must_use]
    pub const fn upper(self) -> u64 {
        self.upper
    }

    /// Returns the number of rows.
    #[must_use]
    pub const fn len(self) -> u64 {
        self.upper - self.lower
    }

    /// Returns whether this interval contains no rows.
    #[must_use]
    pub const fn is_empty(self) -> bool {
        self.lower == self.upper
    }

    /// Returns the private suffix-count domain tag used by packed backends.
    #[doc(hidden)]
    #[must_use]
    pub const fn private_suffix_count(self) -> u64 {
        self.suffix_count
    }
}

impl fmt::Display for FmInterval {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "[{}, {})", self.lower, self.upper)
    }
}

/// An offset into the indexed text.
///
/// The boundary value equal to the text length denotes the terminal suffix and
/// is present when locating the empty pattern.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct TextOffset(u64);

impl TextOffset {
    /// Returns the zero-based text offset, which may equal the text length.
    #[must_use]
    pub const fn get(self) -> u64 {
        self.0
    }
}

impl fmt::Display for TextOffset {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.fmt(formatter)
    }
}

/// Immutable exact FM-index reference representation.
///
/// The implementation stores a complete suffix array and full occurrence
/// prefixes. It is intended as a correctness oracle, not the persisted
/// compressed layout.
pub struct FmIndex {
    text_len: u64,
    suffix_count: u64,
    suffix_array: Vec<usize>,
    bwt: Vec<u8>,
    rank_prefixes: Vec<[u64; 4]>,
    first_occurrence: [u64; 4],
}

impl FmIndex {
    /// Builds an exact reference index over an already-projected canonical text.
    ///
    /// The input is not retained or modified. The index owns every derived
    /// structure and is immutable after successful construction.
    ///
    /// # Errors
    ///
    /// Returns a structured logical-size, explicit-limit, physical-size, or
    /// allocation error before publishing a partial index.
    pub fn build_reference(text: &[SearchBase], limit: FmBuildLimit) -> Result<Self, FmError> {
        let text_len =
            u64::try_from(text.len()).map_err(|_| FmError::TextLengthNotRepresentable {
                text_len: text.len(),
            })?;
        let dimensions = validate_logical_dimensions(text_len, limit)?;
        let capacities = preflight_build(dimensions)?;

        let suffix_array = build_suffix_array(text, dimensions, capacities)?;
        let bwt = build_bwt(text, &suffix_array, dimensions, capacities)?;
        let (rank_prefixes, first_occurrence) = build_rank_prefixes(&bwt, dimensions, capacities)?;

        debug_assert_eq!(suffix_array.len(), capacities.suffix_storage);
        debug_assert_eq!(bwt.len(), capacities.suffix_storage);
        debug_assert_eq!(rank_prefixes.len(), capacities.rank_row_storage);

        Ok(Self {
            text_len,
            suffix_count: dimensions.suffix_count,
            suffix_array,
            bwt,
            rank_prefixes,
            first_occurrence,
        })
    }

    /// Returns the indexed text length, excluding the terminal suffix.
    #[must_use]
    pub const fn text_len(&self) -> u64 {
        self.text_len
    }

    /// Returns the suffix-array and BWT row count.
    #[must_use]
    pub const fn suffix_count(&self) -> u64 {
        self.suffix_count
    }

    /// Returns the interval containing every suffix, including the terminal
    /// suffix.
    #[must_use]
    pub const fn full_interval(&self) -> FmInterval {
        FmInterval::new_unchecked(0, self.suffix_count, self.suffix_count)
    }

    /// Constructs a checked row interval for this index.
    ///
    /// Inversion is diagnosed before bounds.
    ///
    /// # Errors
    ///
    /// Returns a structured interval error when the bounds are inverted or
    /// exceed this index's suffix count.
    pub const fn interval(&self, lower: u64, upper: u64) -> Result<FmInterval, FmError> {
        match validate_interval(lower, upper, self.suffix_count) {
            Ok(()) => Ok(FmInterval::new_unchecked(lower, upper, self.suffix_count)),
            Err(error) => Err(error),
        }
    }

    /// Counts one base in a BWT prefix.
    ///
    /// Position zero counts an empty prefix. Position equal to
    /// `suffix_count()` counts the complete BWT.
    ///
    /// # Errors
    ///
    /// Returns a structured bounds error when the boundary exceeds the BWT
    /// length.
    pub fn rank(&self, base: SearchBase, boundary: u64) -> Result<u64, FmError> {
        debug_assert_eq!(self.bwt.len() + 1, self.rank_prefixes.len());
        if boundary > self.suffix_count {
            return Err(FmError::RankBoundaryOutOfBounds {
                boundary,
                suffix_count: self.suffix_count,
            });
        }
        let storage = usize::try_from(boundary).map_err(|_| FmError::RankBoundaryOutOfBounds {
            boundary,
            suffix_count: self.suffix_count,
        })?;
        Ok(self.rank_prefixes[storage][base.rank_index()])
    }

    /// Searches for an exact canonical pattern.
    ///
    /// The empty pattern returns the full interval. Remaining characters are
    /// processed even after an interval becomes empty so its exact insertion
    /// point remains correct.
    #[must_use]
    pub fn exact_search(&self, pattern: &[SearchBase]) -> FmInterval {
        let mut interval = self.full_interval();
        for &base in pattern.iter().rev() {
            let base_index = base.rank_index();
            let lower_rank = self.rank_unchecked(base_index, interval.lower);
            let upper_rank = self.rank_unchecked(base_index, interval.upper);
            interval = FmInterval::new_unchecked(
                self.first_occurrence[base_index] + lower_rank,
                self.first_occurrence[base_index] + upper_rank,
                self.suffix_count,
            );
            debug_assert!(interval.upper <= self.suffix_count);
        }
        interval
    }

    /// Locates every suffix row in an interval.
    ///
    /// This reference backend emits suffix-array row order. Callers must treat
    /// the offset set, not this backend-specific order, as semantic identity.
    /// The full interval includes the terminal offset equal to `text_len()`.
    ///
    /// # Errors
    ///
    /// Returns a structured interval, size, or output-allocation error.
    pub fn locate(&self, interval: FmInterval) -> Result<Vec<TextOffset>, FmError> {
        if interval.suffix_count != self.suffix_count {
            return Err(FmError::IntervalDomainMismatch {
                interval_suffix_count: interval.suffix_count,
                index_suffix_count: self.suffix_count,
            });
        }
        validate_interval(interval.lower, interval.upper, self.suffix_count)?;
        let logical_elements = interval.len();
        let storage_elements =
            preflight_allocation::<TextOffset>(FmAllocation::LocateResults, logical_elements)?;
        let lower = usize::try_from(interval.lower).map_err(|_| FmError::IntervalOutOfBounds {
            lower: interval.lower,
            upper: interval.upper,
            suffix_count: self.suffix_count,
        })?;
        let upper = usize::try_from(interval.upper).map_err(|_| FmError::IntervalOutOfBounds {
            lower: interval.lower,
            upper: interval.upper,
            suffix_count: self.suffix_count,
        })?;

        let mut offsets = try_vector(
            FmAllocation::LocateResults,
            logical_elements,
            storage_elements,
        )?;
        for &position in &self.suffix_array[lower..upper] {
            let logical = u64::try_from(position)
                .map_err(|_| FmError::TextLengthNotRepresentable { text_len: position })?;
            offsets.push(TextOffset(logical));
        }
        Ok(offsets)
    }

    /// Streams suffix offsets in FM-row order without allocating a result
    /// vector. Returning `false` stops after the current offset.
    pub(crate) fn visit_locate(
        &self,
        interval: FmInterval,
        visitor: &mut dyn FnMut(u64) -> bool,
    ) -> Result<u64, FmError> {
        if interval.suffix_count != self.suffix_count {
            return Err(FmError::IntervalDomainMismatch {
                interval_suffix_count: interval.suffix_count,
                index_suffix_count: self.suffix_count,
            });
        }
        validate_interval(interval.lower, interval.upper, self.suffix_count)?;
        let lower = usize::try_from(interval.lower).map_err(|_| FmError::IntervalOutOfBounds {
            lower: interval.lower,
            upper: interval.upper,
            suffix_count: self.suffix_count,
        })?;
        let upper = usize::try_from(interval.upper).map_err(|_| FmError::IntervalOutOfBounds {
            lower: interval.lower,
            upper: interval.upper,
            suffix_count: self.suffix_count,
        })?;
        let mut located = 0_u64;
        for &position in &self.suffix_array[lower..upper] {
            let logical = u64::try_from(position)
                .map_err(|_| FmError::TextLengthNotRepresentable { text_len: position })?;
            located = located
                .checked_add(1)
                .ok_or(FmError::TextLengthNotRepresentable { text_len: position })?;
            if !visitor(logical) {
                break;
            }
        }
        Ok(located)
    }

    fn rank_unchecked(&self, base_index: usize, boundary: u64) -> u64 {
        let Ok(storage) = usize::try_from(boundary) else {
            unreachable!("validated FM boundary must fit physical storage");
        };
        self.rank_prefixes[storage][base_index]
    }

    #[cfg(test)]
    fn bwt(&self) -> &[u8] {
        &self.bwt
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct BuildDimensions {
    suffix_count: u64,
    rank_rows: u64,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct BuildCapacities {
    suffix_storage: usize,
    rank_row_storage: usize,
}

const fn validate_logical_dimensions(
    text_len: u64,
    limit: FmBuildLimit,
) -> Result<BuildDimensions, FmError> {
    let Some(suffix_count) = text_len.checked_add(1) else {
        return Err(FmError::SuffixCountOverflow { text_len });
    };
    if suffix_count > limit.get() {
        return Err(FmError::BuildLimitExceeded {
            suffix_count,
            max_suffix_count: limit.get(),
        });
    }
    let Some(rank_rows) = suffix_count.checked_add(1) else {
        return Err(FmError::RankRowCountOverflow { suffix_count });
    };
    Ok(BuildDimensions {
        suffix_count,
        rank_rows,
    })
}

fn preflight_build(dimensions: BuildDimensions) -> Result<BuildCapacities, FmError> {
    let suffix_storage =
        preflight_allocation::<usize>(FmAllocation::SuffixArray, dimensions.suffix_count)?;
    preflight_allocation::<u64>(FmAllocation::RankClasses, dimensions.suffix_count)?;
    preflight_allocation::<u64>(FmAllocation::NextRankClasses, dimensions.suffix_count)?;
    preflight_allocation::<u8>(FmAllocation::Bwt, dimensions.suffix_count)?;
    let rank_row_storage =
        preflight_allocation::<[u64; 4]>(FmAllocation::RankPrefixes, dimensions.rank_rows)?;
    Ok(BuildCapacities {
        suffix_storage,
        rank_row_storage,
    })
}

fn preflight_allocation<T>(component: FmAllocation, elements: u64) -> Result<usize, FmError> {
    let element_size = size_of::<T>();
    let logical_element_size =
        u64::try_from(element_size).map_err(|_| FmError::AllocationSizeOverflow {
            component,
            elements,
            element_size: u64::MAX,
        })?;
    elements
        .checked_mul(logical_element_size)
        .ok_or(FmError::AllocationSizeOverflow {
            component,
            elements,
            element_size: logical_element_size,
        })?;
    let storage = usize::try_from(elements).map_err(|_| FmError::AllocationSizeOverflow {
        component,
        elements,
        element_size: logical_element_size,
    })?;
    if element_size != 0 && storage > isize::MAX.unsigned_abs() / element_size {
        return Err(FmError::AllocationSizeOverflow {
            component,
            elements,
            element_size: logical_element_size,
        });
    }
    Ok(storage)
}

fn try_vector<T>(
    component: FmAllocation,
    logical_elements: u64,
    storage_elements: usize,
) -> Result<Vec<T>, FmError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(storage_elements)
        .map_err(|_| FmError::AllocationFailed {
            component,
            elements: logical_elements,
        })?;
    Ok(values)
}

fn build_suffix_array(
    text: &[SearchBase],
    dimensions: BuildDimensions,
    capacities: BuildCapacities,
) -> Result<Vec<usize>, FmError> {
    let mut suffix_array = try_vector(
        FmAllocation::SuffixArray,
        dimensions.suffix_count,
        capacities.suffix_storage,
    )?;
    suffix_array.extend(0..capacities.suffix_storage);

    let mut ranks = try_vector(
        FmAllocation::RankClasses,
        dimensions.suffix_count,
        capacities.suffix_storage,
    )?;
    ranks.resize(capacities.suffix_storage, 0_u64);
    let mut next_ranks = try_vector(
        FmAllocation::NextRankClasses,
        dimensions.suffix_count,
        capacities.suffix_storage,
    )?;
    next_ranks.resize(capacities.suffix_storage, 0_u64);

    suffix_array.sort_unstable_by_key(|&position| (initial_symbol_key(text, position), position));

    let first_position = suffix_array[0];
    ranks[first_position] = 0;
    let mut classes = 1_u64;
    for row in 1..suffix_array.len() {
        let previous = suffix_array[row - 1];
        let current = suffix_array[row];
        if initial_symbol_key(text, previous) != initial_symbol_key(text, current) {
            classes += 1;
        }
        ranks[current] = classes - 1;
    }

    let mut width = 1_usize;
    while classes < dimensions.suffix_count {
        suffix_array.sort_unstable_by(|left, right| {
            rank_pair(*left, width, &ranks)
                .cmp(&rank_pair(*right, width, &ranks))
                .then_with(|| left.cmp(right))
        });

        next_ranks[suffix_array[0]] = 0;
        classes = 1;
        for row in 1..suffix_array.len() {
            let previous = suffix_array[row - 1];
            let current = suffix_array[row];
            if rank_pair(previous, width, &ranks) != rank_pair(current, width, &ranks) {
                classes += 1;
            }
            next_ranks[current] = classes - 1;
        }
        core::mem::swap(&mut ranks, &mut next_ranks);
        if classes < dimensions.suffix_count {
            width = width.saturating_mul(2).min(capacities.suffix_storage);
        }
    }

    Ok(suffix_array)
}

fn initial_symbol_key(text: &[SearchBase], position: usize) -> u8 {
    text.get(position)
        .map_or(SENTINEL_CODE, |base| base.bwt_code())
}

fn rank_pair(position: usize, width: usize, ranks: &[u64]) -> (u64, Option<u64>) {
    let second = position
        .checked_add(width)
        .and_then(|next| ranks.get(next))
        .copied();
    (ranks[position], second)
}

fn build_bwt(
    text: &[SearchBase],
    suffix_array: &[usize],
    dimensions: BuildDimensions,
    capacities: BuildCapacities,
) -> Result<Vec<u8>, FmError> {
    let mut bwt = try_vector(
        FmAllocation::Bwt,
        dimensions.suffix_count,
        capacities.suffix_storage,
    )?;
    for &start in suffix_array {
        let code = if start == 0 {
            SENTINEL_CODE
        } else {
            text[start - 1].bwt_code()
        };
        bwt.push(code);
    }
    Ok(bwt)
}

fn build_rank_prefixes(
    bwt: &[u8],
    dimensions: BuildDimensions,
    capacities: BuildCapacities,
) -> Result<(Vec<[u64; 4]>, [u64; 4]), FmError> {
    let mut prefixes = try_vector(
        FmAllocation::RankPrefixes,
        dimensions.rank_rows,
        capacities.rank_row_storage,
    )?;
    let mut counts = [0_u64; 4];
    prefixes.push(counts);
    for &symbol in bwt {
        if symbol != SENTINEL_CODE {
            counts[usize::from(symbol - 1)] += 1;
        }
        prefixes.push(counts);
    }

    let first_occurrence = [
        1,
        1 + counts[0],
        1 + counts[0] + counts[1],
        1 + counts[0] + counts[1] + counts[2],
    ];
    Ok((prefixes, first_occurrence))
}

const fn validate_interval(lower: u64, upper: u64, suffix_count: u64) -> Result<(), FmError> {
    if lower > upper {
        return Err(FmError::InvertedInterval { lower, upper });
    }
    if upper > suffix_count {
        return Err(FmError::IntervalOutOfBounds {
            lower,
            upper,
            suffix_count,
        });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn bases(text: &str) -> Vec<SearchBase> {
        text.bytes()
            .map(|byte| match byte {
                b'A' => SearchBase::A,
                b'C' => SearchBase::C,
                b'G' => SearchBase::G,
                b'T' => SearchBase::T,
                _ => panic!("invalid test base"),
            })
            .collect()
    }

    fn built(text: &str) -> FmIndex {
        FmIndex::build_reference(&bases(text), FmBuildLimit::MAX)
            .expect("small test index should build")
    }

    fn offsets(index: &FmIndex, pattern: &str) -> Vec<u64> {
        index
            .locate(index.exact_search(&bases(pattern)))
            .expect("search interval must locate")
            .iter()
            .map(|offset| offset.get())
            .collect()
    }

    #[test]
    fn logical_dimension_errors_are_ordered() {
        assert_eq!(
            validate_logical_dimensions(u64::MAX, FmBuildLimit::new(0)),
            Err(FmError::SuffixCountOverflow { text_len: u64::MAX })
        );
        assert_eq!(
            validate_logical_dimensions(4, FmBuildLimit::new(4)),
            Err(FmError::BuildLimitExceeded {
                suffix_count: 5,
                max_suffix_count: 4,
            })
        );
        assert_eq!(
            validate_logical_dimensions(u64::MAX - 1, FmBuildLimit::new(u64::MAX - 1)),
            Err(FmError::BuildLimitExceeded {
                suffix_count: u64::MAX,
                max_suffix_count: u64::MAX - 1,
            })
        );
        assert_eq!(
            validate_logical_dimensions(u64::MAX - 1, FmBuildLimit::MAX),
            Err(FmError::RankRowCountOverflow {
                suffix_count: u64::MAX,
            })
        );
        assert_eq!(
            validate_logical_dimensions(4, FmBuildLimit::new(5)),
            Ok(BuildDimensions {
                suffix_count: 5,
                rank_rows: 6,
            })
        );
    }

    #[test]
    fn physical_preflight_reports_the_first_component() {
        let too_many_u64 = u64::try_from(isize::MAX.unsigned_abs() / size_of::<usize>())
            .expect("usize capacity boundary fits u64")
            + 1;
        assert_eq!(
            preflight_build(BuildDimensions {
                suffix_count: too_many_u64,
                rank_rows: too_many_u64,
            }),
            Err(FmError::AllocationSizeOverflow {
                component: FmAllocation::SuffixArray,
                elements: too_many_u64,
                element_size: u64::try_from(size_of::<usize>()).expect("element size fits u64"),
            })
        );
        assert_eq!(
            preflight_allocation::<[u64; 4]>(FmAllocation::RankPrefixes, u64::MAX),
            Err(FmError::AllocationSizeOverflow {
                component: FmAllocation::RankPrefixes,
                elements: u64::MAX,
                element_size: 32,
            })
        );
    }

    #[test]
    fn private_named_suffix_arrays_and_bwts_are_exact() {
        let fixtures: [(&str, &[usize], &[u8]); 6] = [
            ("", &[0], &[0]),
            ("A", &[1, 0], &[1, 0]),
            ("T", &[1, 0], &[4, 0]),
            ("AAAA", &[4, 3, 2, 1, 0], &[1, 1, 1, 1, 0]),
            ("ACAC", &[4, 2, 0, 3, 1], &[2, 2, 0, 1, 1]),
            ("ACGT", &[4, 0, 1, 2, 3], &[4, 0, 1, 2, 3]),
        ];
        for (text, expected_sa, expected_bwt) in fixtures {
            let index = built(text);
            assert_eq!(index.suffix_array, expected_sa, "{text}");
            assert_eq!(index.bwt(), expected_bwt, "{text}");
        }
    }

    #[test]
    fn empty_interval_keeps_moving_to_its_exact_insertion_point() {
        let index = built("C");
        let single_base_interval = index.exact_search(&bases("A"));
        assert_eq!(
            (single_base_interval.lower(), single_base_interval.upper()),
            (1, 1)
        );

        let two_base_interval = index.exact_search(&bases("CA"));
        assert_eq!(
            (two_base_interval.lower(), two_base_interval.upper()),
            (2, 2)
        );
        assert_eq!(offsets(&index, "CA"), Vec::<u64>::new());
    }

    #[test]
    fn empty_pattern_and_terminal_suffix_are_exact() {
        let empty = built("");
        assert_eq!(empty.full_interval(), FmInterval::new_unchecked(0, 1, 1));
        assert_eq!(offsets(&empty, ""), vec![0]);
        assert_eq!(
            empty.exact_search(&bases("A")),
            FmInterval::new_unchecked(1, 1, 1)
        );

        let singleton = built("A");
        assert_eq!(offsets(&singleton, ""), vec![1, 0]);
        assert_eq!(offsets(&singleton, "A"), vec![0]);
        assert_eq!(singleton.suffix_array[0], 1);
        assert_eq!(
            singleton
                .bwt()
                .iter()
                .position(|&symbol| symbol == SENTINEL_CODE),
            Some(1)
        );
    }
}
