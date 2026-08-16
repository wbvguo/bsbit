//! Private, semantics-preserving extension filters.
//!
//! The public Level 1 scalar APIs remain the exact oracle. This module supplies
//! an ungapped exact/one-mismatch shortcut and one-/two-word Myers distance
//! sweeps for the candidate-extension layer.

use super::MAX_QUERY_BASES;

#[cfg(test)]
use super::{KernelFlavor, myers_distance_batch};

use crate::score::EditDistance;
use crate::verification::distance::{DistanceError, MatrixAllocation, TracebackResult};
use bsbit_core::alphabet::Base;
use bsbit_core::bisulfite::{CytosineStrand, classify_bases};

/// Smallest query for which the bit-vector prefix filter is worthwhile.
pub const MIN_FILTER_QUERY_BASES: u64 = 4;
const MAX_PREFIX_MYERS_QUERY_BASES: usize = 128;

/// Prepared bisulfite equality masks for one- or two-word Myers filtering.
pub struct WordMyersQuery {
    equality_masks: [u128; 5],
    query_length: usize,
}

impl WordMyersQuery {
    /// Prepares a query of at most 128 bases for prefix-distance filtering.
    #[must_use]
    pub fn new(query: &[Base], strand: CytosineStrand) -> Option<Self> {
        if !(1..=MAX_PREFIX_MYERS_QUERY_BASES).contains(&query.len()) {
            return None;
        }
        let mut equality_masks = [0_u128; 5];
        for (position, &query_base) in query.iter().enumerate() {
            let bit = 1_u128 << position;
            for reference_base in Base::ALL {
                if classify_bases(reference_base, query_base, strand).is_zero_cost() {
                    equality_masks[base_code(reference_base)] |= bit;
                }
            }
        }
        Some(Self {
            equality_masks,
            query_length: query.len(),
        })
    }

    #[cfg(test)]
    // The test-only path is entered only for one-word queries, so all upper
    // 64 bits are known to be zero before the conversion.
    #[allow(clippy::cast_possible_truncation)]
    pub(crate) fn distances(
        &self,
        reference_codes: &[u8],
        starts: &[usize],
        lengths: &[usize],
        output: &mut [u64],
    ) -> Option<KernelFlavor> {
        if self.query_length > MAX_QUERY_BASES {
            return None;
        }
        let equality_masks = self.equality_masks.map(|mask| mask as u64);
        myers_distance_batch(
            &equality_masks,
            self.query_length,
            reference_codes,
            starts,
            lengths,
            output,
        )
        .ok()
    }

    /// Computes the edit distance to every nonempty reference prefix.
    pub fn prefix_distances(&self, reference: &[Base], output: &mut [u64]) -> bool {
        self.prefix_distances_from_codes(
            reference.len(),
            reference.iter().copied().map(base_code_u8),
            output,
        )
    }

    #[cfg(test)]
    fn prefix_distances_encoded(&self, reference_codes: &[u8], output: &mut [u64]) -> bool {
        self.prefix_distances_from_codes(
            reference_codes.len(),
            reference_codes.iter().copied(),
            output,
        )
    }

    fn prefix_distances_from_codes(
        &self,
        reference_length: usize,
        reference_codes: impl Iterator<Item = u8>,
        output: &mut [u64],
    ) -> bool {
        if output.len() != reference_length {
            return false;
        }
        if self.query_length <= MAX_QUERY_BASES {
            self.prefix_distances_u64(reference_codes, output);
        } else {
            self.prefix_distances_u128(reference_codes, output);
        }
        true
    }

    // Dispatch reaches this recurrence only for queries of at most 64 bases,
    // so every equality mask is exactly representable as u64.
    #[allow(clippy::cast_possible_truncation)]
    fn prefix_distances_u64(&self, reference_codes: impl Iterator<Item = u8>, output: &mut [u64]) {
        let mut positive = !0_u64;
        let mut negative = 0_u64;
        let mut score = u64::try_from(self.query_length).expect("query length fits u64");
        let high_bit = 1_u64 << (self.query_length - 1);
        for (code, result) in reference_codes.zip(output) {
            let equal = equality_mask_u128(&self.equality_masks, code) as u64;
            let horizontal_input = equal | negative;
            let horizontal = (((equal & positive).wrapping_add(positive)) ^ positive) | equal;
            let positive_horizontal = negative | !(horizontal | positive);
            let negative_horizontal = positive & horizontal;
            if positive_horizontal & high_bit != 0 {
                score += 1;
            } else if negative_horizontal & high_bit != 0 {
                score -= 1;
            }
            let shifted_positive = (positive_horizontal << 1) | 1;
            let shifted_negative = negative_horizontal << 1;
            positive = shifted_negative | !(horizontal_input | shifted_positive);
            negative = shifted_positive & horizontal_input;
            *result = score;
        }
    }

    fn prefix_distances_u128(&self, reference_codes: impl Iterator<Item = u8>, output: &mut [u64]) {
        let mut positive = !0_u128;
        let mut negative = 0_u128;
        let mut score = u64::try_from(self.query_length).expect("query length fits u64");
        let high_bit = 1_u128 << (self.query_length - 1);
        for (code, result) in reference_codes.zip(output) {
            let equal = equality_mask_u128(&self.equality_masks, code);
            let horizontal_input = equal | negative;
            let horizontal = (((equal & positive).wrapping_add(positive)) ^ positive) | equal;
            let positive_horizontal = negative | !(horizontal | positive);
            let negative_horizontal = positive & horizontal;
            if positive_horizontal & high_bit != 0 {
                score += 1;
            } else if negative_horizontal & high_bit != 0 {
                score -= 1;
            }
            let shifted_positive = (positive_horizontal << 1) | 1;
            let shifted_negative = negative_horizontal << 1;
            positive = shifted_negative | !(horizontal_input | shifted_positive);
            negative = shifted_positive & horizontal_input;
            *result = score;
        }
    }
}

/// Reusable scalar workspace for bounded prefix-distance sweeps.
pub struct BandedPrefixDistanceWorkspace {
    query_length: u64,
    max_edit_distance: u64,
    capped_distance: u64,
    previous: Vec<u64>,
    current: Vec<u64>,
}

impl BandedPrefixDistanceWorkspace {
    /// Allocates a checked workspace for one query length and edit budget.
    pub fn new(query_length: u64, max_edit_distance: EditDistance) -> Result<Self, DistanceError> {
        let query_storage = usize::try_from(query_length).map_err(|_| {
            matrix_size_error(
                0,
                query_length,
                MatrixAllocation::LogicalQueryExtent,
                Some(query_length),
            )
        })?;
        let columns = query_storage.checked_add(1).ok_or_else(|| {
            matrix_size_error(0, query_length, MatrixAllocation::LogicalQueryExtent, None)
        })?;
        let row_elements = columns.checked_mul(2).ok_or_else(|| {
            matrix_size_error(0, query_length, MatrixAllocation::DistanceRows, None)
        })?;
        let row_elements_u64 = u64::try_from(row_elements).map_err(|_| {
            matrix_size_error(0, query_length, MatrixAllocation::DistanceRows, None)
        })?;
        let row_bytes = row_elements_u64.checked_mul(8).ok_or_else(|| {
            matrix_size_error(
                0,
                query_length,
                MatrixAllocation::DistanceRows,
                Some(row_elements_u64),
            )
        })?;
        if usize::try_from(row_bytes).is_err()
            || row_bytes > u64::try_from(usize::MAX >> 1).unwrap()
        {
            return Err(matrix_size_error(
                0,
                query_length,
                MatrixAllocation::DistanceRows,
                Some(row_elements_u64),
            ));
        }

        let mut previous = Vec::new();
        previous.try_reserve_exact(columns).map_err(|_| {
            matrix_size_error(
                0,
                query_length,
                MatrixAllocation::DistanceRows,
                Some(row_elements_u64),
            )
        })?;
        previous.extend(0..=query_length);
        let mut current = Vec::new();
        current.try_reserve_exact(columns).map_err(|_| {
            matrix_size_error(
                0,
                query_length,
                MatrixAllocation::DistanceRows,
                Some(row_elements_u64),
            )
        })?;
        current.resize(columns, 0);
        Ok(Self {
            query_length,
            max_edit_distance: max_edit_distance.get(),
            capped_distance: max_edit_distance.get().saturating_add(1),
            previous,
            current,
        })
    }

    /// Computes capped distances to every supplied reference prefix.
    pub fn prefix_distances(
        &mut self,
        reference: &[Base],
        query: &[Base],
        strand: CytosineStrand,
        output: &mut [u64],
    ) -> Result<u64, DistanceError> {
        let (reference_length, query_length) =
            self.validate_prefix_request(reference, query, output)?;

        self.previous.fill(self.capped_distance);
        let initial_maximum = self.max_edit_distance.min(query_length);
        for column in 0..=initial_maximum {
            let column_storage = usize::try_from(column).map_err(|_| {
                matrix_size_error(
                    reference_length,
                    query_length,
                    MatrixAllocation::LogicalQueryExtent,
                    Some(column),
                )
            })?;
            self.previous[column_storage] = column;
        }
        let mut previous_minimum = 0_u64;
        let mut previous_maximum = initial_maximum;
        let mut distance_updates = 0_u64;
        for (reference_index, (&reference_base, result)) in reference.iter().zip(output).enumerate()
        {
            let row = u64::try_from(reference_index + 1).map_err(|_| {
                matrix_size_error(
                    reference_length,
                    query_length,
                    MatrixAllocation::LogicalReferenceExtent,
                    None,
                )
            })?;
            let current_minimum = row.saturating_sub(self.max_edit_distance);
            let current_maximum = row.saturating_add(self.max_edit_distance).min(query_length);
            if current_minimum > current_maximum {
                previous_minimum = 1;
                previous_maximum = 0;
                *result = self.capped_distance;
                continue;
            }
            if current_minimum == 0 {
                self.current[0] = row.min(self.capped_distance);
            }
            let recurrence_minimum = current_minimum.max(1);
            for column in recurrence_minimum..=current_maximum {
                let column_storage = usize::try_from(column).map_err(|_| {
                    matrix_size_error(
                        reference_length,
                        query_length,
                        MatrixAllocation::LogicalQueryExtent,
                        Some(column),
                    )
                })?;
                let query_index = column_storage - 1;
                let deletion = if in_band(column, previous_minimum, previous_maximum) {
                    capped_distance_add(self.previous[column_storage], 1, self.capped_distance)
                } else {
                    self.capped_distance
                };
                let insertion =
                    if in_band(column.saturating_sub(1), current_minimum, current_maximum) {
                        capped_distance_add(
                            self.current[column_storage - 1],
                            1,
                            self.capped_distance,
                        )
                    } else {
                        self.capped_distance
                    };
                let diagonal =
                    if in_band(column.saturating_sub(1), previous_minimum, previous_maximum) {
                        capped_distance_add(
                            self.previous[column_storage - 1],
                            classify_bases(reference_base, query[query_index], strand).cost(),
                            self.capped_distance,
                        )
                    } else {
                        self.capped_distance
                    };
                self.current[column_storage] = deletion.min(insertion).min(diagonal);
                distance_updates = checked_distance_add(distance_updates, 1)?;
            }
            *result = if in_band(query_length, current_minimum, current_maximum) {
                self.current[query.len()]
            } else {
                self.capped_distance
            };
            core::mem::swap(&mut self.previous, &mut self.current);
            previous_minimum = current_minimum;
            previous_maximum = current_maximum;
        }
        Ok(distance_updates)
    }

    fn validate_prefix_request(
        &self,
        reference: &[Base],
        query: &[Base],
        output: &[u64],
    ) -> Result<(u64, u64), DistanceError> {
        let reference_length = u64::try_from(reference.len()).map_err(|_| {
            matrix_size_error(
                u64::MAX,
                self.query_length,
                MatrixAllocation::LogicalReferenceExtent,
                None,
            )
        })?;
        let query_length = u64::try_from(query.len()).map_err(|_| {
            matrix_size_error(
                reference_length,
                u64::MAX,
                MatrixAllocation::LogicalQueryExtent,
                None,
            )
        })?;
        if query_length != self.query_length {
            return Err(matrix_size_error(
                reference_length,
                query_length,
                MatrixAllocation::LogicalQueryExtent,
                Some(self.query_length),
            ));
        }
        if output.len() != reference.len() {
            return Err(matrix_size_error(
                reference_length,
                query_length,
                MatrixAllocation::LogicalReferenceExtent,
                Some(u64::try_from(output.len()).unwrap_or(u64::MAX)),
            ));
        }
        Ok((reference_length, query_length))
    }
}

const fn in_band(index: u64, minimum: u64, maximum: u64) -> bool {
    minimum <= index && index <= maximum
}

fn capped_distance_add(accumulated: u64, increment: u64, cap: u64) -> u64 {
    if accumulated >= cap {
        cap
    } else {
        accumulated.saturating_add(increment).min(cap)
    }
}

#[cfg(test)]
pub(crate) fn encode_bases(bases: &[Base]) -> Vec<u8> {
    bases.iter().copied().map(base_code_u8).collect()
}

/// Returns an exact ungapped traceback when there is at most one mismatch.
pub fn ungapped_traceback_at_most_one(
    reference: &[Base],
    query: &[Base],
    strand: CytosineStrand,
) -> Result<Option<TracebackResult>, DistanceError> {
    ungapped_traceback_with_literal_nm::<false, false>(reference, query, strand)
        .map(|result| result.map(|(traceback, _)| traceback))
}

/// Returns an exact ungapped traceback for a certified distance up to two.
pub fn ungapped_traceback_at_most_two_certified_cached_nm(
    reference: &[Base],
    query: &[Base],
    strand: CytosineStrand,
) -> Result<Option<(TracebackResult, u64)>, DistanceError> {
    ungapped_traceback_with_literal_nm::<true, true>(reference, query, strand)
}

fn ungapped_traceback_with_literal_nm<
    const TRACK_LITERAL_NM: bool,
    const CERTIFY_TWO_MISMATCHES: bool,
>(
    reference: &[Base],
    query: &[Base],
    strand: CytosineStrand,
) -> Result<Option<(TracebackResult, u64)>, DistanceError> {
    if reference.len() != query.len() {
        return Ok(None);
    }
    let mut mismatch_count = 0_u32;
    let mut mismatch_positions = [0_usize; 2];
    let mut recorded_mismatches = 0_usize;
    let mut literal_nm = 0_u64;
    for (chunk_start, (reference_chunk, query_chunk)) in reference
        .chunks(64)
        .zip(query.chunks(64))
        .enumerate()
        .map(|(chunk, pair)| (chunk * 64, pair))
    {
        let mut mismatch_bits = 0_u64;
        for (position, (&reference_base, &query_base)) in
            reference_chunk.iter().zip(query_chunk).enumerate()
        {
            if TRACK_LITERAL_NM
                && (reference_base != query_base
                    || !matches!(reference_base, Base::A | Base::C | Base::G | Base::T))
            {
                literal_nm = literal_nm.saturating_add(1);
            }
            if !classify_bases(reference_base, query_base, strand).is_zero_cost() {
                mismatch_bits |= 1_u64 << position;
            }
        }
        let mut remaining_mismatches = mismatch_bits;
        while remaining_mismatches != 0 && recorded_mismatches < mismatch_positions.len() {
            mismatch_positions[recorded_mismatches] = chunk_start
                + usize::try_from(remaining_mismatches.trailing_zeros())
                    .expect("trailing-zero count fits usize");
            recorded_mismatches += 1;
            remaining_mismatches &= remaining_mismatches - 1;
        }
        mismatch_count += mismatch_bits.count_ones();
        let maximum_mismatches = if CERTIFY_TWO_MISMATCHES { 2 } else { 1 };
        if mismatch_count > maximum_mismatches {
            return Ok(None);
        }
    }
    if mismatch_count == 2 {
        let [first, second] = mismatch_positions;
        debug_assert!(first < second && second < reference.len());
        let forward_shift_tie = reference[first..second]
            .iter()
            .zip(&query[first + 1..=second])
            .all(|(&reference_base, &query_base)| {
                classify_bases(reference_base, query_base, strand).is_zero_cost()
            });
        let reverse_shift_tie = reference[first + 1..=second]
            .iter()
            .zip(&query[first..second])
            .all(|(&reference_base, &query_base)| {
                classify_bases(reference_base, query_base, strand).is_zero_cost()
            });
        if forward_shift_tie || reverse_shift_tie {
            return Ok(None);
        }
    }
    let length = u64::try_from(reference.len()).expect("supported pointer width fits u64");
    TracebackResult::ungapped(EditDistance::new(u64::from(mismatch_count)), length)
        .map(|traceback| Some((traceback, literal_nm)))
}

fn base_code(base: Base) -> usize {
    usize::from(base_code_u8(base))
}

const fn equality_mask_u128(masks: &[u128; 5], code: u8) -> u128 {
    match code {
        0 => masks[0],
        1 => masks[1],
        2 => masks[2],
        3 => masks[3],
        _ => 0,
    }
}

fn checked_distance_add(accumulated: u64, increment: u64) -> Result<u64, DistanceError> {
    EditDistance::new(accumulated)
        .checked_add(increment)
        .map_err(DistanceError::from)
        .map(EditDistance::get)
}

const fn matrix_size_error(
    reference_length: u64,
    query_length: u64,
    allocation: MatrixAllocation,
    requested_elements: Option<u64>,
) -> DistanceError {
    DistanceError::MatrixSizeOverflow {
        reference_length,
        query_length,
        allocation,
        requested_elements,
        requested_bytes: None,
    }
}

const fn base_code_u8(base: Base) -> u8 {
    match base {
        Base::A => 0,
        Base::C => 1,
        Base::G => 2,
        Base::T => 3,
        _ => 4,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::verification::distance::{DpCellLimit, global_bs_alignment, global_bs_distance};
    use bsbit_core::sequence::NormalizedSequence;

    #[test]
    fn word_myers_equals_scalar_global_distance_exhaustively() {
        let references = sequences_through(3);
        let queries = sequences_through(3)
            .into_iter()
            .filter(|sequence| !sequence.is_empty())
            .collect::<Vec<_>>();
        for strand in [CytosineStrand::Top, CytosineStrand::Bottom] {
            for query in &queries {
                let kernel = WordMyersQuery::new(query.bases(), strand).expect("short query");
                for reference in &references {
                    let encoded = encode_bases(reference.bases());
                    let mut observed = [u64::MAX];
                    kernel
                        .distances(&encoded, &[0], &[encoded.len()], &mut observed)
                        .expect("kernel dispatch");
                    let expected = global_bs_distance(reference, query, strand, DpCellLimit::MAX)
                        .expect("scalar distance");
                    assert_eq!(
                        observed[0],
                        expected.get(),
                        "{strand:?} {reference:?} {query:?}"
                    );
                }
            }
        }
    }

    #[test]
    fn word_myers_prefix_sweep_equals_every_scalar_prefix() {
        let references = sequences_through(4);
        let queries = sequences_through(3)
            .into_iter()
            .filter(|sequence| !sequence.is_empty())
            .collect::<Vec<_>>();
        for strand in [CytosineStrand::Top, CytosineStrand::Bottom] {
            for query in &queries {
                let kernel = WordMyersQuery::new(query.bases(), strand).expect("short query");
                for reference in &references {
                    let encoded = encode_bases(reference.bases());
                    let mut observed = vec![u64::MAX; encoded.len()];
                    assert!(kernel.prefix_distances_encoded(&encoded, &mut observed));
                    for (prefix_length, &distance) in observed.iter().enumerate() {
                        let prefix = NormalizedSequence::from_bases(
                            reference.bases()[..=prefix_length].iter().copied(),
                        );
                        let expected = global_bs_distance(&prefix, query, strand, DpCellLimit::MAX)
                            .expect("scalar distance");
                        assert_eq!(
                            distance,
                            expected.get(),
                            "{strand:?} {reference:?} {query:?} prefix={}",
                            prefix_length + 1
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn two_word_myers_prefix_sweep_equals_every_scalar_prefix() {
        for query_length in [65_usize, 100, 127, 128] {
            let query = patterned_sequence(query_length, 3);
            let reference = patterned_sequence(query_length + 2, 11);
            for strand in [CytosineStrand::Top, CytosineStrand::Bottom] {
                let kernel = WordMyersQuery::new(query.bases(), strand).expect("two-word query");
                let mut observed = vec![u64::MAX; reference.bases().len()];
                assert!(kernel.prefix_distances(reference.bases(), &mut observed));
                for (prefix_length, &distance) in observed.iter().enumerate() {
                    let prefix = NormalizedSequence::from_bases(
                        reference.bases()[..=prefix_length].iter().copied(),
                    );
                    let expected = global_bs_distance(&prefix, &query, strand, DpCellLimit::MAX)
                        .expect("scalar distance");
                    assert_eq!(
                        distance,
                        expected.get(),
                        "length={query_length} {strand:?} prefix={}",
                        prefix_length + 1
                    );
                }
            }
        }
    }

    #[test]
    fn banded_prefix_sweep_is_exact_inside_budget_and_rejects_outside() {
        let references = sequences_through(4);
        let queries = sequences_through(3)
            .into_iter()
            .filter(|sequence| !sequence.is_empty())
            .collect::<Vec<_>>();
        for strand in [CytosineStrand::Top, CytosineStrand::Bottom] {
            for query in &queries {
                for budget in 0..=3 {
                    let mut workspace =
                        BandedPrefixDistanceWorkspace::new(query.len(), EditDistance::new(budget))
                            .unwrap();
                    for reference in &references {
                        let mut observed = vec![u64::MAX; reference.bases().len()];
                        let updates = workspace
                            .prefix_distances(
                                reference.bases(),
                                query.bases(),
                                strand,
                                &mut observed,
                            )
                            .expect("banded sweep");
                        assert_eq!(
                            updates,
                            expected_banded_updates(reference.len(), query.len(), budget)
                        );
                        for (prefix_length, &distance) in observed.iter().enumerate() {
                            let prefix = NormalizedSequence::from_bases(
                                reference.bases()[..=prefix_length].iter().copied(),
                            );
                            let expected =
                                global_bs_distance(&prefix, query, strand, DpCellLimit::MAX)
                                    .expect("scalar distance");
                            if expected.get() <= budget {
                                assert_eq!(
                                    distance,
                                    expected.get(),
                                    "{strand:?} {reference:?} {query:?} budget={budget} prefix={}",
                                    prefix_length + 1
                                );
                            } else {
                                assert!(
                                    distance > budget,
                                    "{strand:?} {reference:?} {query:?} budget={budget} prefix={} observed={distance} expected={}",
                                    prefix_length + 1,
                                    expected.get()
                                );
                            }
                        }
                    }
                }
            }
        }
    }

    fn expected_banded_updates(reference_length: u64, query_length: u64, budget: u64) -> u64 {
        (1..=reference_length)
            .map(|row| {
                let first = row.saturating_sub(budget).max(1);
                let last = row.saturating_add(budget).min(query_length);
                last.saturating_sub(first)
                    .saturating_add(u64::from(first <= last))
            })
            .sum()
    }

    fn patterned_sequence(length: usize, offset: usize) -> NormalizedSequence {
        let alphabet = [Base::A, Base::C, Base::G, Base::T, Base::N];
        NormalizedSequence::from_bases(
            (0..length).map(|position| alphabet[(position * 7 + offset) % alphabet.len()]),
        )
    }

    #[test]
    fn ungapped_shortcut_exactly_matches_scalar_zero_and_one_paths() {
        let sequences = sequences_through(3);
        for strand in [CytosineStrand::Top, CytosineStrand::Bottom] {
            for reference in &sequences {
                for query in &sequences {
                    let observed =
                        ungapped_traceback_at_most_one(reference.bases(), query.bases(), strand)
                            .expect("shortcut");
                    let expected = global_bs_alignment(reference, query, strand, DpCellLimit::MAX)
                        .expect("scalar alignment");
                    if reference.len() == query.len() && expected.distance().get() <= 1 {
                        assert_eq!(observed.as_ref(), Some(&expected));
                    } else {
                        assert!(observed.is_none());
                    }
                }
            }
        }
    }

    #[test]
    fn cached_literal_nm_counts_non_acgt_self_pairs() {
        let reference = [Base::A, Base::N, Base::C, Base::T];
        let observed = ungapped_traceback_at_most_two_certified_cached_nm(
            &reference,
            &reference,
            CytosineStrand::Top,
        )
        .expect("cached shortcut")
        .expect("zero-distance ungapped alignment");
        assert_eq!(observed.0.distance(), EditDistance::new(1));
        assert_eq!(observed.1, 1);
    }

    #[test]
    fn certified_two_mismatch_shortcut_matches_scalar_traceback() {
        let sequences = sequences_through(3);
        let mut admitted_two_mismatch_paths = 0_usize;
        for strand in [CytosineStrand::Top, CytosineStrand::Bottom] {
            for reference in &sequences {
                for query in &sequences {
                    let observed = ungapped_traceback_at_most_two_certified_cached_nm(
                        reference.bases(),
                        query.bases(),
                        strand,
                    )
                    .expect("certified shortcut");
                    let expected = global_bs_alignment(reference, query, strand, DpCellLimit::MAX)
                        .expect("scalar alignment");
                    if let Some((traceback, literal_nm)) = observed {
                        assert_eq!(traceback, expected);
                        let expected_literal_nm = reference
                            .bases()
                            .iter()
                            .zip(query.bases())
                            .filter(|(reference_base, query_base)| {
                                reference_base != query_base
                                    || !matches!(
                                        **reference_base,
                                        Base::A | Base::C | Base::G | Base::T
                                    )
                            })
                            .count();
                        assert_eq!(literal_nm, u64::try_from(expected_literal_nm).unwrap());
                        if traceback.distance().get() == 2 {
                            admitted_two_mismatch_paths += 1;
                        }
                    }
                }
            }
        }
        assert!(admitted_two_mismatch_paths > 0);
    }

    #[test]
    fn ungapped_shortcut_crosses_multiple_words_without_losing_one_mismatch() {
        let reference = NormalizedSequence::from_bases(
            (0..130).map(|position| Base::CANONICAL[position % Base::CANONICAL.len()]),
        );
        let exact = reference.clone();
        let mut one_mismatch = reference.bases().to_vec();
        one_mismatch[96] = one_mismatch[96].complement();
        let one_mismatch = NormalizedSequence::from_bases(one_mismatch);
        for query in [&exact, &one_mismatch] {
            let observed = ungapped_traceback_at_most_one(
                reference.bases(),
                query.bases(),
                CytosineStrand::Top,
            )
            .unwrap()
            .expect("zero/one mismatch shortcut");
            let expected =
                global_bs_alignment(&reference, query, CytosineStrand::Top, DpCellLimit::MAX)
                    .unwrap();
            assert_eq!(observed, expected);
        }
    }

    fn sequences_through(maximum_length: usize) -> Vec<NormalizedSequence> {
        let mut sequences = Vec::new();
        for length in 0..=maximum_length {
            let count = 5_usize.pow(u32::try_from(length).unwrap());
            for ordinal in 0..count {
                let mut remaining = ordinal;
                let mut bases = Vec::with_capacity(length);
                for _ in 0..length {
                    bases.push(Base::ALL[remaining % Base::ALL.len()]);
                    remaining /= Base::ALL.len();
                }
                sequences.push(NormalizedSequence::from_bases(bases));
            }
        }
        sequences
    }
}
