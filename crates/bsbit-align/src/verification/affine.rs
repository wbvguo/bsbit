//! Banded affine alignment scoring for an already selected reference slice.

use core::ops::Range;

use bsbit_core::alphabet::Base;
use bsbit_core::bisulfite::{
    AlignmentOrientation, BisulfiteStrand, classify_bases, strand_semantics,
};

const MATCH_SCORE: i16 = 1;
const MISMATCH_PENALTY: i16 = 4;
const GAP_OPEN_PENALTY: i16 = 6;
const GAP_EXTENSION_PENALTY: i16 = 1;
const MAX_EDIT_DISTANCE: usize = 5;
const AFFINE_BAND: usize = 2 * MAX_EDIT_DISTANCE;
const MAX_QUERY_BASES: usize = 192;
const NEGATIVE_INFINITY: i16 = i16::MIN / 4;
const ROW_CELLS: usize = MAX_QUERY_BASES + 2 * MAX_EDIT_DISTANCE + 1;

/// Reusable fixed-capacity rows for bounded affine scoring.
#[derive(Clone)]
pub struct AffineScoreWorkspace {
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
        || query_len > MAX_QUERY_BASES
        || query_len.abs_diff(reference_len) > AFFINE_BAND
    {
        return None;
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
    for query_position in 1..=query_len {
        let query_base = match semantics.orientation() {
            AlignmentOrientation::Forward => read[retained_query.start + query_position - 1],
            AlignmentOrientation::Reverse => read[retained_query.end - query_position].complement(),
        };
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
            let substitution = if classify_bases(
                reference[reference_position - 1],
                query_base,
                semantics.cytosine_strand(),
            )
            .is_zero_cost()
            {
                MATCH_SCORE
            } else {
                -MISMATCH_PENALTY
            };
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
