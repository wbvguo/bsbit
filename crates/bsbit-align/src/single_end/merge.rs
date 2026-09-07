//! Origin grouping, fair representative selection, and library-pass merging.

use bsbit_core::alphabet::Base;
use bsbit_core::bisulfite::BisulfiteStrand;
use bsbit_index::reference::ReferenceIndex;

use super::result::{SingleAlignmentResult, SingleMappingStatus};
use crate::AlignmentError;
use crate::placement::{ReadPlacement, placement_net_gap_bases};
use crate::read_mapping::ReadWorkspace;
use crate::reporting_tie_break::{
    ReportingTieBreak, compare_placements_without_reference_order, placement_origin_hash,
};
use crate::single_end::mapq::{
    affine_rerank_origin_count_supported, cross_pass_mapping_quality_cap,
    merged_incomplete_repeat_mapping_quality,
};
use crate::verification::affine::{AffineScoreWorkspace, affine_placement_score};

pub(super) fn origins_share_local_locus(
    left: (u64, BisulfiteStrand, i128),
    right: (u64, BisulfiteStrand, i128),
    edit_radius: u8,
) -> bool {
    left.0 == right.0 && left.1 == right.1 && left.2.abs_diff(right.2) <= u128::from(edit_radius)
}

pub(super) fn local_origin_locus_count(
    sorted_origins: &[(u64, BisulfiteStrand, i128)],
    edit_radius: u8,
) -> usize {
    let mut clusters = 0usize;
    let mut cluster_start: Option<(u64, BisulfiteStrand, i128)> = None;
    for &origin in sorted_origins {
        if cluster_start.is_none_or(|start| {
            start.0 != origin.0
                || start.1 != origin.1
                || origin.2.saturating_sub(start.2) > i128::from(edit_radius)
        }) {
            clusters = clusters.saturating_add(1);
            cluster_start = Some(origin);
        }
    }
    clusters
}

#[cold]
#[inline(never)]
pub(super) fn prefer_minimum_net_gap_representative(
    placements: &[ReadPlacement],
    best_distance: u8,
    read_length: usize,
) -> Option<ReadPlacement> {
    placements
        .iter()
        .copied()
        .filter(|placement| placement.distance() == best_distance)
        .min_by_key(|placement| (placement_net_gap_bases(*placement, read_length), *placement))
}

pub(super) fn prefer_fair_ambiguous_representative(
    reference: &ReferenceIndex,
    read: &[Base],
    workspace: &ReadWorkspace,
    mut result: SingleAlignmentResult,
    tie_break: ReportingTieBreak,
) -> Result<SingleAlignmentResult, AlignmentError> {
    if !matches!(result.status, SingleMappingStatus::Ambiguous) || result.mapping_quality != 0 {
        return Ok(result);
    }
    let Some(current) = result.placement else {
        return Ok(result);
    };
    let bounded_affine_tie = affine_rerank_origin_count_supported(result.best_origin_count);
    if bounded_affine_tie
        && let Some((_, retained_affine_score)) = workspace
            .affine_scores
            .iter()
            .find(|(placement, _)| *placement == current)
    {
        let mut best = current;
        let mut best_hash = placement_origin_hash(reference, tie_break, current, read.len())?;
        for &(placement, score) in &workspace.affine_scores {
            if score != *retained_affine_score {
                continue;
            }
            let hash = placement_origin_hash(reference, tie_break, placement, read.len())?;
            if hash < best_hash
                || (hash == best_hash
                    && compare_placements_without_reference_order(
                        reference,
                        placement,
                        best,
                        read.len(),
                    )?
                    .is_lt())
            {
                best = placement;
                best_hash = hash;
            }
        }
        result.placement = Some(best);
        return Ok(result);
    }
    let mut affine_workspace = AffineScoreWorkspace::default();
    let retained_affine_score = if bounded_affine_tie {
        Some(affine_placement_score(
            reference,
            read,
            current,
            0,
            &mut affine_workspace,
        )?)
    } else {
        None
    };
    let mut best = current;
    let mut best_hash = placement_origin_hash(reference, tie_break, current, read.len())?;
    // Reporting remains a seeded lottery only after all bounded alignment
    // objectives tie. This preserves reference-order independence without
    // discarding a real affine distinction between equal-edit origins.
    for placement in workspace
        .placements
        .iter()
        .copied()
        .filter(|placement| placement.distance() == current.distance())
    {
        if let Some(retained_score) = retained_affine_score
            && affine_placement_score(reference, read, placement, 0, &mut affine_workspace)?
                != retained_score
        {
            continue;
        }
        let hash = placement_origin_hash(reference, tie_break, placement, read.len())?;
        if hash < best_hash
            || (hash == best_hash
                && compare_placements_without_reference_order(
                    reference,
                    placement,
                    best,
                    read.len(),
                )?
                .is_lt())
        {
            best = placement;
            best_hash = hash;
        }
    }
    result.placement = Some(best);
    Ok(result)
}

#[cfg(test)]
pub(super) fn merge_non_directional_results(
    original: SingleAlignmentResult,
    complementary: SingleAlignmentResult,
) -> SingleAlignmentResult {
    merge_non_directional_results_with_tie_break(None, original, complementary, false)
        .expect("a disabled tie-break does not access a reference")
}

#[cfg(test)]
pub(super) fn merge_non_directional_completed_frontiers(
    original: SingleAlignmentResult,
    complementary: SingleAlignmentResult,
) -> SingleAlignmentResult {
    merge_non_directional_results_with_tie_break(None, original, complementary, true)
        .expect("a disabled tie-break does not access a reference")
}

pub(super) fn merge_non_directional_results_with_tie_break(
    tie_break: Option<(&ReferenceIndex, &[Base], ReportingTieBreak)>,
    original: SingleAlignmentResult,
    complementary: SingleAlignmentResult,
    frontiers_complete: bool,
) -> Result<SingleAlignmentResult, AlignmentError> {
    let original_distance = original.placement.map(ReadPlacement::distance);
    let complementary_distance = complementary.placement.map(ReadPlacement::distance);
    let (mut selected, other, tied) = match (original_distance, complementary_distance) {
        (Some(left), Some(right)) if left < right => (original, complementary, false),
        (Some(left), Some(right)) if right < left => (complementary, original, false),
        (Some(_), Some(_)) => {
            let select_complementary = if let (
                Some((reference, read, tie_break)),
                Some(original_placement),
                Some(complementary_placement),
            ) =
                (tie_break, original.placement, complementary.placement)
            {
                let original_hash =
                    placement_origin_hash(reference, tie_break, original_placement, read.len())?;
                let complementary_hash = placement_origin_hash(
                    reference,
                    tie_break,
                    complementary_placement,
                    read.len(),
                )?;
                complementary_hash < original_hash
                    || (complementary_hash == original_hash
                        && compare_placements_without_reference_order(
                            reference,
                            complementary_placement,
                            original_placement,
                            read.len(),
                        )?
                        .is_lt())
            } else {
                false
            };
            if select_complementary {
                (complementary, original, true)
            } else {
                (original, complementary, true)
            }
        }
        (None, Some(_)) => (complementary, original, false),
        (_, None) => (original, complementary, false),
    };
    selected.located_rows = original
        .located_rows
        .saturating_add(complementary.located_rows);
    selected.distinct_candidate_starts = original
        .distinct_candidate_starts
        .saturating_add(complementary.distinct_candidate_starts);
    selected.verified_placements = original
        .verified_placements
        .saturating_add(complementary.verified_placements);

    if tied {
        selected.best_origin_count = original
            .best_origin_count
            .saturating_add(complementary.best_origin_count);
        selected.status = SingleMappingStatus::Ambiguous;
        selected.mapping_quality = 0;
    } else if matches!(selected.status, SingleMappingStatus::Unique) {
        if let (Some(best), Some(runner_up)) = (
            selected.placement.map(ReadPlacement::distance),
            other.placement.map(ReadPlacement::distance),
        ) {
            selected.mapping_quality = selected
                .mapping_quality
                .min(cross_pass_mapping_quality_cap(best, runner_up));
        }
        selected.mapping_quality = merged_incomplete_repeat_mapping_quality(
            selected.mapping_quality,
            frontiers_complete,
            selected.located_rows,
            selected.distinct_candidate_starts,
            selected.verified_placements,
        );
    } else {
        selected.mapping_quality = 0;
    }
    Ok(selected)
}
