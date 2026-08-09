//! Checked index-construction boundaries.

#[cfg(feature = "index-construction")]
#[allow(unsafe_code)]
pub mod combined;
#[cfg(feature = "index-construction")]
#[allow(unsafe_code)]
pub(crate) mod combined_blocks;
#[cfg(feature = "index-construction")]
#[allow(unsafe_code)]
pub(crate) mod libsais;

#[cfg(test)]
use core::fmt;

#[cfg(test)]
use crate::storage::fm::SearchBase;

#[cfg(test)]
const MIN_PARALLEL_VALIDATION_ROWS: usize = 262_144;

/// Failure to construct or validate one complete suffix array.
#[cfg(test)]
#[derive(Clone, Debug, Eq, PartialEq)]
#[doc(hidden)]
#[non_exhaustive]
pub(crate) enum SuffixArrayBuildError {
    /// The text plus its terminal suffix exceeds the packed `u32` domain.
    #[cfg(test)]
    TextExceedsU32,
    /// A checked physical dimension overflowed.
    #[cfg(test)]
    SizeOverflow,
    /// A fallible allocation failed.
    AllocationFailed {
        /// Logical allocation site.
        component: &'static str,
        /// Requested element count.
        elements: usize,
    },
    /// An isolated implementation reported a backend-specific failure.
    #[cfg(feature = "index-construction")]
    Backend {
        /// Stable backend identifier.
        backend: &'static str,
        /// Owned diagnostic without borrowed native or process storage.
        message: String,
    },
    /// The builder returned the wrong number of suffix rows.
    #[cfg(test)]
    InvalidLength {
        /// Expected rows, including the terminal suffix.
        expected: usize,
        /// Returned rows.
        observed: usize,
    },
    /// A suffix offset lies outside the terminal-inclusive domain.
    #[cfg(test)]
    OffsetOutOfBounds {
        /// Row containing the invalid offset.
        row: usize,
        /// Returned suffix offset.
        offset: u64,
        /// Largest permitted offset.
        text_len: usize,
    },
    /// Two rows contain the same suffix offset.
    #[cfg(test)]
    DuplicateOffset {
        /// Later row containing the duplicate.
        row: usize,
        /// Duplicated suffix offset.
        offset: u32,
    },
    /// The terminal suffix was not the first lexicographic row.
    #[cfg(test)]
    TerminalNotFirst {
        /// Observed first suffix offset.
        observed: u32,
        /// Required terminal offset.
        expected: u32,
    },
    /// Adjacent rows violate suffix lexicographic order.
    #[cfg(test)]
    InvalidOrder {
        /// Later row in the invalid adjacent pair.
        row: usize,
        /// Earlier suffix offset.
        previous: u32,
        /// Later suffix offset.
        current: u32,
    },
}

#[cfg(test)]
impl fmt::Display for SuffixArrayBuildError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            #[cfg(test)]
            Self::TextExceedsU32 => formatter.write_str("suffix-array text exceeds u32 offsets"),
            #[cfg(test)]
            Self::SizeOverflow => formatter.write_str("suffix-array dimensions overflow"),
            Self::AllocationFailed {
                component,
                elements,
            } => write!(
                formatter,
                "allocate {elements} suffix-array elements for {component}: failed"
            ),
            #[cfg(feature = "index-construction")]
            Self::Backend { backend, message } => {
                write!(formatter, "suffix-array backend {backend}: {message}")
            }
            #[cfg(test)]
            Self::InvalidLength { expected, observed } => write!(
                formatter,
                "suffix-array row count {observed} differs from expected {expected}"
            ),
            #[cfg(test)]
            Self::OffsetOutOfBounds {
                row,
                offset,
                text_len,
            } => write!(
                formatter,
                "suffix-array row {row} has offset {offset} outside 0..={text_len}"
            ),
            #[cfg(test)]
            Self::DuplicateOffset { row, offset } => {
                write!(formatter, "suffix-array row {row} repeats offset {offset}")
            }
            #[cfg(test)]
            Self::TerminalNotFirst { observed, expected } => write!(
                formatter,
                "suffix-array first row is {observed}, expected terminal offset {expected}"
            ),
            #[cfg(test)]
            Self::InvalidOrder {
                row,
                previous,
                current,
            } => write!(
                formatter,
                "suffix-array rows {} and {row} are not ordered: {previous} before {current}",
                row - 1
            ),
        }
    }
}

#[cfg(test)]
impl std::error::Error for SuffixArrayBuildError {}

/// One interchangeable suffix-array constructor.
///
/// Implementations return every suffix of `text` plus the empty terminal
/// suffix, ordered lexicographically. Packed image construction validates an
/// alternate builder's complete result before consuming it.
#[doc(hidden)]
#[cfg(test)]
pub(crate) trait SuffixArrayBuilder: Send + Sync {
    /// Stable identifier recorded by build output.
    #[cfg(feature = "index-construction")]
    fn backend_name(&self) -> &'static str;

    /// Builds terminal-inclusive `u32` suffix positions.
    ///
    /// # Errors
    ///
    /// Returns a checked resource, backend, or dimension failure without a
    /// partial suffix array.
    fn build_suffix_array(&self, text: &[SearchBase]) -> Result<Vec<u32>, SuffixArrayBuildError>;
}

/// Safe Rust prefix-doubling constructor retained as the validation baseline.
#[derive(Clone, Copy, Debug, Default)]
#[doc(hidden)]
#[cfg(test)]
pub(crate) struct PrefixDoublingSuffixArrayBuilder;

#[cfg(test)]
impl SuffixArrayBuilder for PrefixDoublingSuffixArrayBuilder {
    #[cfg(feature = "index-construction")]
    fn backend_name(&self) -> &'static str {
        "rust-prefix-doubling"
    }

    fn build_suffix_array(&self, text: &[SearchBase]) -> Result<Vec<u32>, SuffixArrayBuildError> {
        build_prefix_doubling(text)
    }
}

/// Validates a terminal-inclusive suffix array in linear time and memory.
///
/// This is the trust boundary between the pinned libsais FFI and the Rust
/// correctness oracle. Inverse ranks prove adjacent equal-prefix ordering without
/// comparing whole suffix slices.
///
/// # Errors
///
/// Returns the first dimension, permutation, terminal, or ordering defect.
#[doc(hidden)]
#[cfg(test)]
pub(crate) fn validate_suffix_array(
    text: &[SearchBase],
    suffix_array: &[u32],
) -> Result<(), SuffixArrayBuildError> {
    validate_suffix_array_with_threads(text, suffix_array, 1)
}

/// Validates a terminal-inclusive suffix array with parallel adjacent-order
/// checks after the deterministic permutation scan.
///
/// The inverse-rank construction remains serial because it establishes the
/// unique-offset trust boundary. Adjacent order checks are independent and
/// return the earliest invalid row regardless of worker scheduling.
///
/// # Errors
///
/// Returns the first dimension, permutation, terminal, or ordering defect.
#[doc(hidden)]
#[cfg(test)]
pub(crate) fn validate_suffix_array_with_threads(
    text: &[SearchBase],
    suffix_array: &[u32],
    threads: u32,
) -> Result<(), SuffixArrayBuildError> {
    let suffix_storage = checked_suffix_storage(text.len())?;
    if suffix_array.len() != suffix_storage {
        return Err(SuffixArrayBuildError::InvalidLength {
            expected: suffix_storage,
            observed: suffix_array.len(),
        });
    }
    let terminal = u32::try_from(text.len()).map_err(|_| SuffixArrayBuildError::TextExceedsU32)?;
    let mut inverse = reserved_vec(suffix_storage, "inverse suffix ranks")?;
    inverse.resize(suffix_storage, u32::MAX);
    for (row, &offset) in suffix_array.iter().enumerate() {
        let storage =
            usize::try_from(offset).map_err(|_| SuffixArrayBuildError::OffsetOutOfBounds {
                row,
                offset: u64::from(offset),
                text_len: text.len(),
            })?;
        let Some(slot) = inverse.get_mut(storage) else {
            return Err(SuffixArrayBuildError::OffsetOutOfBounds {
                row,
                offset: u64::from(offset),
                text_len: text.len(),
            });
        };
        if *slot != u32::MAX {
            return Err(SuffixArrayBuildError::DuplicateOffset { row, offset });
        }
        *slot = u32::try_from(row).map_err(|_| SuffixArrayBuildError::TextExceedsU32)?;
    }
    if suffix_array.first().copied() != Some(terminal) {
        return Err(SuffixArrayBuildError::TerminalNotFirst {
            observed: suffix_array[0],
            expected: terminal,
        });
    }
    let pair_count = suffix_array.len().saturating_sub(1);
    if pair_count == 0 {
        return Ok(());
    }
    let requested_workers = if suffix_array.len() < MIN_PARALLEL_VALIDATION_ROWS {
        1
    } else {
        threads.max(1)
    };
    let worker_count = usize::try_from(requested_workers)
        .unwrap_or(usize::MAX)
        .min(pair_count);
    let pairs_per_worker = pair_count.div_ceil(worker_count);
    let failure = std::thread::scope(|scope| {
        let mut workers = Vec::with_capacity(worker_count);
        let inverse = &inverse;
        for first_pair in (0..pair_count).step_by(pairs_per_worker) {
            let end_pair = first_pair.saturating_add(pairs_per_worker).min(pair_count);
            workers.push(scope.spawn(move || {
                for pair in first_pair..end_pair {
                    let row = pair + 1;
                    let previous = suffix_array[pair];
                    let current = suffix_array[row];
                    if !ordered_pair(text, inverse, previous, current) {
                        return Some((row, previous, current));
                    }
                }
                None
            }));
        }
        workers.into_iter().find_map(|worker| {
            worker
                .join()
                .expect("suffix validation worker cannot panic")
        })
    });
    if let Some((row, previous, current)) = failure {
        return Err(SuffixArrayBuildError::InvalidOrder {
            row,
            previous,
            current,
        });
    }
    Ok(())
}

#[cfg(test)]
fn ordered_pair(text: &[SearchBase], inverse: &[u32], previous: u32, current: u32) -> bool {
    let previous = usize::try_from(previous).expect("validated u32 offset fits usize");
    let current = usize::try_from(current).expect("validated u32 offset fits usize");
    let previous_key = initial_key(text, previous);
    let current_key = initial_key(text, current);
    if previous_key != current_key {
        return previous_key < current_key;
    }
    if previous == text.len() || current == text.len() {
        return false;
    }
    inverse[previous + 1] < inverse[current + 1]
}

#[cfg(test)]
fn build_prefix_doubling(text: &[SearchBase]) -> Result<Vec<u32>, SuffixArrayBuildError> {
    let suffix_storage = checked_suffix_storage(text.len())?;
    let suffix_count =
        u64::try_from(suffix_storage).map_err(|_| SuffixArrayBuildError::SizeOverflow)?;
    let mut suffix_array = reserved_vec(suffix_storage, "suffix array")?;
    for position in 0..suffix_storage {
        suffix_array
            .push(u32::try_from(position).map_err(|_| SuffixArrayBuildError::TextExceedsU32)?);
    }
    let mut ranks = reserved_vec(suffix_storage, "rank classes")?;
    ranks.resize(suffix_storage, 0_u32);
    let mut next_ranks = reserved_vec(suffix_storage, "next rank classes")?;
    next_ranks.resize(suffix_storage, 0_u32);

    suffix_array.sort_unstable_by_key(|&position| {
        (
            initial_key(
                text,
                usize::try_from(position).expect("u32 position fits usize"),
            ),
            position,
        )
    });
    let first_position = usize::try_from(suffix_array[0]).expect("u32 position fits usize");
    ranks[first_position] = 0;
    let mut classes = 1_u64;
    for row in 1..suffix_array.len() {
        let previous = suffix_array[row - 1];
        let current = suffix_array[row];
        if initial_key(
            text,
            usize::try_from(previous).expect("u32 position fits usize"),
        ) != initial_key(
            text,
            usize::try_from(current).expect("u32 position fits usize"),
        ) {
            classes += 1;
        }
        ranks[usize::try_from(current).expect("u32 position fits usize")] =
            u32::try_from(classes - 1).expect("rank class fits u32");
    }

    let mut width = 1_usize;
    while classes < suffix_count {
        suffix_array.sort_unstable_by(|left, right| {
            rank_pair(*left, width, &ranks)
                .cmp(&rank_pair(*right, width, &ranks))
                .then_with(|| left.cmp(right))
        });
        next_ranks[usize::try_from(suffix_array[0]).expect("u32 position fits usize")] = 0;
        classes = 1;
        for row in 1..suffix_array.len() {
            let previous = suffix_array[row - 1];
            let current = suffix_array[row];
            if rank_pair(previous, width, &ranks) != rank_pair(current, width, &ranks) {
                classes += 1;
            }
            next_ranks[usize::try_from(current).expect("u32 position fits usize")] =
                u32::try_from(classes - 1).expect("rank class fits u32");
        }
        core::mem::swap(&mut ranks, &mut next_ranks);
        if classes < suffix_count {
            width = width.saturating_mul(2).min(suffix_storage);
        }
    }
    Ok(suffix_array)
}

#[cfg(test)]
fn checked_suffix_storage(text_len: usize) -> Result<usize, SuffixArrayBuildError> {
    let logical = u64::try_from(text_len).map_err(|_| SuffixArrayBuildError::TextExceedsU32)?;
    if logical > u64::from(u32::MAX) {
        return Err(SuffixArrayBuildError::TextExceedsU32);
    }
    text_len
        .checked_add(1)
        .ok_or(SuffixArrayBuildError::SizeOverflow)
}

#[cfg(test)]
fn initial_key(text: &[SearchBase], position: usize) -> u8 {
    text.get(position).map_or(0, |base| base.bwt_code())
}

#[cfg(test)]
fn rank_pair(position: u32, width: usize, ranks: &[u32]) -> (u32, Option<u32>) {
    let position = usize::try_from(position).expect("u32 position fits usize");
    let second = position
        .checked_add(width)
        .and_then(|next| ranks.get(next))
        .copied();
    (ranks[position], second)
}

#[cfg(test)]
fn reserved_vec<T>(
    elements: usize,
    component: &'static str,
) -> Result<Vec<T>, SuffixArrayBuildError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(elements)
        .map_err(|_| SuffixArrayBuildError::AllocationFailed {
            component,
            elements,
        })?;
    Ok(values)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn text(bytes: &[u8]) -> Vec<SearchBase> {
        bytes
            .iter()
            .map(|byte| match byte {
                b'A' => SearchBase::A,
                b'C' => SearchBase::C,
                b'G' => SearchBase::G,
                b'T' => SearchBase::T,
                _ => panic!("canonical fixture"),
            })
            .collect()
    }

    #[test]
    fn prefix_doubling_and_linear_validator_are_exact() {
        let builder = PrefixDoublingSuffixArrayBuilder;
        let input = text(b"ACGACGTA");
        let suffixes = builder.build_suffix_array(&input).expect("builds");
        validate_suffix_array(&input, &suffixes).expect("validates");
        let mut expected: Vec<u32> = (0..=u32::try_from(input.len()).unwrap()).collect();
        expected.sort_unstable_by(|left, right| {
            let left = usize::try_from(*left).unwrap();
            let right = usize::try_from(*right).unwrap();
            input[left..].cmp(&input[right..])
        });
        assert_eq!(suffixes, expected);
    }

    #[test]
    fn validator_rejects_length_bounds_duplicates_terminal_and_order() {
        let input = text(b"ACG");
        assert!(matches!(
            validate_suffix_array(&input, &[3, 0, 1]),
            Err(SuffixArrayBuildError::InvalidLength { .. })
        ));
        assert!(matches!(
            validate_suffix_array(&input, &[3, 0, 1, 4]),
            Err(SuffixArrayBuildError::OffsetOutOfBounds { row: 3, .. })
        ));
        assert!(matches!(
            validate_suffix_array(&input, &[3, 0, 1, 1]),
            Err(SuffixArrayBuildError::DuplicateOffset { row: 3, .. })
        ));
        assert!(matches!(
            validate_suffix_array(&input, &[0, 3, 1, 2]),
            Err(SuffixArrayBuildError::TerminalNotFirst { .. })
        ));
        assert!(matches!(
            validate_suffix_array(&input, &[3, 1, 0, 2]),
            Err(SuffixArrayBuildError::InvalidOrder { row: 2, .. })
        ));
    }
}
