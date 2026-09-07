//! Banded affine alignment scoring for an already selected reference slice.

use core::ops::Range;

use bsbit_core::alphabet::Base;
use bsbit_core::bisulfite::{
    AlignmentOrientation, BisulfiteStrand, CytosineStrand, classify_bases, strand_semantics,
};
use bsbit_index::reference::ReferenceIndex;

use crate::AlignmentError;
use crate::placement::ReadPlacement;
use crate::read_mapping_limits::{MAX_EDIT_DISTANCE, MAX_READ_BASES};

const MATCH_SCORE: i16 = 1;
const MISMATCH_PENALTY: i16 = 4;
const GAP_OPEN_PENALTY: i16 = 6;
const GAP_EXTENSION_PENALTY: i16 = 1;
const AFFINE_BAND: usize = 2 * MAX_EDIT_DISTANCE as usize;
const NEGATIVE_INFINITY: i16 = i16::MIN / 4;
const ROW_CELLS: usize = MAX_READ_BASES + 2 * MAX_EDIT_DISTANCE as usize + 1;
const BASE_COUNT: usize = Base::ALL.len();
const SUBSTITUTION_CELLS: usize = BASE_COUNT * BASE_COUNT;
const TOP_SUBSTITUTION_SCORES: [i16; SUBSTITUTION_CELLS] =
    substitution_score_table(CytosineStrand::Top);
const BOTTOM_SUBSTITUTION_SCORES: [i16; SUBSTITUTION_CELLS] =
    substitution_score_table(CytosineStrand::Bottom);

const fn substitution_score_table(strand: CytosineStrand) -> [i16; SUBSTITUTION_CELLS] {
    let mut scores = [-MISMATCH_PENALTY; SUBSTITUTION_CELLS];
    let mut query = 0;
    while query < BASE_COUNT {
        let mut reference = 0;
        while reference < BASE_COUNT {
            if classify_bases(Base::ALL[reference], Base::ALL[query], strand).is_zero_cost() {
                scores[query * BASE_COUNT + reference] = MATCH_SCORE;
            }
            reference += 1;
        }
        query += 1;
    }
    scores
}

/// Reusable fixed-capacity rows for bounded affine scoring.
#[derive(Clone)]
pub struct AffineScoreWorkspace {
    reference_codes: [u8; ROW_CELLS],
    match_previous: [i16; ROW_CELLS],
    insertion_previous: [i16; ROW_CELLS],
    deletion_previous: [i16; ROW_CELLS],
    match_current: [i16; ROW_CELLS],
    insertion_current: [i16; ROW_CELLS],
    deletion_current: [i16; ROW_CELLS],
}

impl Default for AffineScoreWorkspace {
    fn default() -> Self {
        Self {
            reference_codes: [Base::N.storage_code(); ROW_CELLS],
            match_previous: [NEGATIVE_INFINITY; ROW_CELLS],
            insertion_previous: [NEGATIVE_INFINITY; ROW_CELLS],
            deletion_previous: [NEGATIVE_INFINITY; ROW_CELLS],
            match_current: [NEGATIVE_INFINITY; ROW_CELLS],
            insertion_current: [NEGATIVE_INFINITY; ROW_CELLS],
            deletion_current: [NEGATIVE_INFINITY; ROW_CELLS],
        }
    }
}

/// Scores one already selected reference/query placement with a fixed band.
///
/// `retained_query` is expressed in sequencing orientation. The bisulfite
/// strand determines both query orientation and zero-cost conversion policy.
/// The returned score includes the supplied linear penalty for clipped query
/// bases.
///
/// Returns `None` when the query, reference slice, retained interval, or band
/// exceeds this bounded kernel's contract.
#[must_use]
// The three affine states and band boundaries advance as one coupled dynamic
// program; keeping the recurrence together makes its invariants auditable.
#[allow(clippy::too_many_lines)]
pub fn banded_affine_score(
    reference: &[Base],
    read: &[Base],
    retained_query: Range<usize>,
    strand: BisulfiteStrand,
    clip_penalty: u8,
    workspace: &mut AffineScoreWorkspace,
) -> Option<i16> {
    if retained_query.start > retained_query.end || retained_query.end > read.len() {
        return None;
    }
    let query_len = retained_query.end - retained_query.start;
    let reference_len = reference.len();
    if reference_len >= ROW_CELLS
        || query_len > MAX_READ_BASES
        || query_len.abs_diff(reference_len) > AFFINE_BAND
    {
        return None;
    }
    for (code, base) in workspace.reference_codes[..reference_len]
        .iter_mut()
        .zip(reference)
    {
        *code = base.storage_code();
    }

    let initial_upper = reference_len.min(AFFINE_BAND);
    workspace.match_previous[..=initial_upper].fill(NEGATIVE_INFINITY);
    workspace.insertion_previous[..=initial_upper].fill(NEGATIVE_INFINITY);
    workspace.deletion_previous[..=initial_upper].fill(NEGATIVE_INFINITY);
    workspace.match_previous[0] = 0;
    for reference_position in 1..=initial_upper {
        workspace.deletion_previous[reference_position] =
            -GAP_OPEN_PENALTY - GAP_EXTENSION_PENALTY * i16::try_from(reference_position).ok()?;
    }
    if initial_upper < reference_len {
        let boundary = initial_upper + 1;
        workspace.match_previous[boundary] = NEGATIVE_INFINITY;
        workspace.insertion_previous[boundary] = NEGATIVE_INFINITY;
        workspace.deletion_previous[boundary] = NEGATIVE_INFINITY;
    }

    let semantics = strand_semantics(strand);
    let substitution_scores = match semantics.cytosine_strand() {
        CytosineStrand::Top => &TOP_SUBSTITUTION_SCORES,
        CytosineStrand::Bottom => &BOTTOM_SUBSTITUTION_SCORES,
    };
    for query_position in 1..=query_len {
        let query_base = match semantics.orientation() {
            AlignmentOrientation::Forward => read[retained_query.start + query_position - 1],
            AlignmentOrientation::Reverse => read[retained_query.end - query_position].complement(),
        };
        let substitution_row = usize::from(query_base.storage_code()) * BASE_COUNT;
        let lower = query_position.saturating_sub(AFFINE_BAND).max(1);
        let upper = reference_len.min(query_position.saturating_add(AFFINE_BAND));
        let left_boundary = lower - 1;
        workspace.match_current[left_boundary] = NEGATIVE_INFINITY;
        workspace.insertion_current[left_boundary] =
            if left_boundary == 0 && query_position <= AFFINE_BAND {
                -GAP_OPEN_PENALTY - GAP_EXTENSION_PENALTY * i16::try_from(query_position).ok()?
            } else {
                NEGATIVE_INFINITY
            };
        workspace.deletion_current[left_boundary] = NEGATIVE_INFINITY;
        if upper < reference_len {
            let right_boundary = upper + 1;
            workspace.match_current[right_boundary] = NEGATIVE_INFINITY;
            workspace.insertion_current[right_boundary] = NEGATIVE_INFINITY;
            workspace.deletion_current[right_boundary] = NEGATIVE_INFINITY;
        }
        for reference_position in lower..=upper {
            let substitution = substitution_scores
                [substitution_row + usize::from(workspace.reference_codes[reference_position - 1])];
            workspace.match_current[reference_position] = workspace.match_previous
                [reference_position - 1]
                .max(workspace.insertion_previous[reference_position - 1])
                .max(workspace.deletion_previous[reference_position - 1])
                .saturating_add(substitution);
            workspace.insertion_current[reference_position] = workspace.match_previous
                [reference_position]
                .max(workspace.deletion_previous[reference_position])
                .saturating_sub(GAP_OPEN_PENALTY + GAP_EXTENSION_PENALTY)
                .max(
                    workspace.insertion_previous[reference_position]
                        .saturating_sub(GAP_EXTENSION_PENALTY),
                );
            workspace.deletion_current[reference_position] = workspace.match_current
                [reference_position - 1]
                .max(workspace.insertion_current[reference_position - 1])
                .saturating_sub(GAP_OPEN_PENALTY + GAP_EXTENSION_PENALTY)
                .max(
                    workspace.deletion_current[reference_position - 1]
                        .saturating_sub(GAP_EXTENSION_PENALTY),
                );
        }
        core::mem::swap(&mut workspace.match_previous, &mut workspace.match_current);
        core::mem::swap(
            &mut workspace.insertion_previous,
            &mut workspace.insertion_current,
        );
        core::mem::swap(
            &mut workspace.deletion_previous,
            &mut workspace.deletion_current,
        );
    }

    let alignment_score = workspace.match_previous[reference_len]
        .max(workspace.insertion_previous[reference_len])
        .max(workspace.deletion_previous[reference_len]);
    let clipped = read.len().saturating_sub(query_len);
    Some(
        alignment_score
            .saturating_sub(i16::from(clip_penalty).saturating_mul(i16::try_from(clipped).ok()?)),
    )
}

pub(crate) fn affine_placement_score(
    reference: &ReferenceIndex,
    read: &[Base],
    placement: ReadPlacement,
    clip_penalty: u8,
    workspace: &mut AffineScoreWorkspace,
) -> Result<i16, AlignmentError> {
    let contig = reference
        .contig_by_ordinal(placement.contig_ordinal())
        .ok_or(AlignmentError::InvalidContigOrdinal {
            ordinal: placement.contig_ordinal(),
        })?;
    let start = usize::try_from(placement.start()).map_err(|_| {
        AlignmentError::CandidateCoordinateOverflow {
            start: placement.start(),
        }
    })?;
    let end = usize::try_from(placement.end()).map_err(|_| {
        AlignmentError::CandidateCoordinateOverflow {
            start: placement.end(),
        }
    })?;
    let reference_bases = contig.sequence().bases().get(start..end).ok_or(
        AlignmentError::CandidateCoordinateOverflow {
            start: placement.end(),
        },
    )?;
    let retained = placement.retained_query_interval(read.len());
    let retained_length = retained.end.saturating_sub(retained.start);
    banded_affine_score(
        reference_bases,
        read,
        retained,
        placement.strand(),
        clip_penalty,
        workspace,
    )
    .ok_or(AlignmentError::UnsupportedReadLength {
        length: retained_length,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn substitution_tables_match_the_canonical_bisulfite_relation() {
        for (strand, scores) in [
            (CytosineStrand::Top, &TOP_SUBSTITUTION_SCORES),
            (CytosineStrand::Bottom, &BOTTOM_SUBSTITUTION_SCORES),
        ] {
            for query in Base::ALL {
                for reference in Base::ALL {
                    let expected = if classify_bases(reference, query, strand).is_zero_cost() {
                        MATCH_SCORE
                    } else {
                        -MISMATCH_PENALTY
                    };
                    let actual = scores[usize::from(query.storage_code()) * BASE_COUNT
                        + usize::from(reference.storage_code())];
                    assert_eq!(actual, expected, "{strand:?}: {reference}/{query}");
                }
            }
        }
    }
}
