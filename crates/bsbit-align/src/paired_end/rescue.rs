//! Bounded mate-rescue window construction and exact-anchor completion.

use super::adapter::best_ungapped_semi_global_placement;
use super::frontier::{append_local_flexible_proof_candidates, balanced_rescue_blocks};
use super::selection::{counterpart_strand, expected_mate2_strand, is_inward};
use super::{
    AlignmentError, Base, CombinedSearchReferenceExt, CombinedSeedHit, ConversionPass,
    FLEXIBLE_NOMINAL_PROOF, INITIAL_EDIT_DISTANCE, MAX_READ_BASES, MateRescueWindow,
    PairedPlacement, ProjectedBase, RESCUE_BLOCKS, ReadAlignmentMetrics, ReadCandidate,
    ReadPlacement, ReadWorkspace, ReferenceIndex, SearchBase,
};

fn rescue_window_contains_candidate(
    windows: &[MateRescueWindow],
    candidate: ReadCandidate,
    maximum_edit_distance: u8,
) -> bool {
    let edit_budget = u64::from(maximum_edit_distance);
    windows.iter().any(|window| {
        window.contig_ordinal == candidate.contig_ordinal()
            && window.strand == candidate.strand()
            && (window.start.saturating_sub(edit_budget)..=window.end.saturating_add(edit_budget))
                .contains(&candidate.start())
    })
}

/// Completes the missing-mate frontier inside every window induced by a
/// fully enumerated anchor frontier. Unlike the initial rescue path, the
/// block count follows the requested edit budget.
#[allow(clippy::too_many_arguments)]
pub(super) fn rescue_from_ranked_anchor_windows(
    workspace: &mut ReadWorkspace,
    rescue_windows: &mut Vec<MateRescueWindow>,
    reference: &ReferenceIndex,
    read: &[Base],
    anchors: &[ReadPlacement],
    rescuing_mate1: bool,
    maximum_template_span: u64,
    maximum_edit_distance: u8,
) -> Result<ReadAlignmentMetrics, AlignmentError> {
    workspace.candidates.clear();
    workspace.candidate_nominals.clear();
    workspace.placements.clear();
    prepare_rescue_windows(
        rescue_windows,
        reference,
        anchors,
        rescuing_mate1,
        maximum_template_span,
    )?;
    for &window in rescue_windows.iter() {
        let contig = reference.contig_by_ordinal(window.contig_ordinal).ok_or(
            AlignmentError::InvalidContigOrdinal {
                ordinal: window.contig_ordinal,
            },
        )?;
        append_local_flexible_proof_candidates(
            read,
            contig.sequence().bases(),
            window,
            maximum_edit_distance,
            &mut workspace.candidate_nominals,
        );
    }
    let (_, metrics) = workspace.verify_candidates_with_budget(
        reference,
        read,
        ReadAlignmentMetrics::default(),
        maximum_edit_distance,
    )?;
    Ok(metrics)
}

// Exact-block enumeration and bounded rescue-window completion share one
// proof budget and one metrics transaction.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub(super) fn rescue_from_combined_exact_blocks(
    workspace: &mut ReadWorkspace,
    rescue_windows: &mut Vec<MateRescueWindow>,
    reference: &ReferenceIndex,
    read: &[Base],
    reversed_projected: &[ProjectedBase],
    anchors: &[ReadPlacement],
    rescuing_mate1: bool,
    maximum_template_span: u64,
    maximum_edit_distance: u8,
    incremental_fallback_requested: bool,
    maximum_located_hits: u64,
) -> Result<ReadAlignmentMetrics, AlignmentError> {
    workspace.candidates.clear();
    workspace.candidate_nominals.clear();
    workspace.placements.clear();
    prepare_rescue_windows(
        rescue_windows,
        reference,
        anchors,
        rescuing_mate1,
        maximum_template_span,
    )?;
    if rescue_windows.is_empty() {
        return Ok(ReadAlignmentMetrics::default());
    }

    let blocks = balanced_rescue_blocks(read.len());
    let mut projected_search = [SearchBase::A; MAX_READ_BASES];
    for (destination, &source) in projected_search.iter_mut().zip(reversed_projected) {
        *destination = match source {
            ProjectedBase::A => SearchBase::A,
            ProjectedBase::G => SearchBase::G,
            ProjectedBase::T => SearchBase::T,
        };
    }
    let mut exact = [None; RESCUE_BLOCKS];
    let mut projected_hits = 0_u64;
    for (ordinal, block) in blocks.into_iter().enumerate() {
        let query_start = usize::from(block.query_start());
        let query_end = usize::from(block.query_end());
        let source = if rescuing_mate1 {
            &read[query_start..query_end]
        } else {
            &read[read.len() - query_end..read.len() - query_start]
        };
        if source.contains(&Base::N) {
            continue;
        }
        let reversed_start = read.len() - query_end;
        let reversed_end = read.len() - query_start;
        let pattern = &projected_search[reversed_start..reversed_end];
        let Some(matches) = reference
            .combined_exact_seed(pattern)
            .map_err(|_| AlignmentError::CombinedIndex)?
        else {
            continue;
        };
        if matches.matched_bases()
            != u64::try_from(query_end - query_start).expect("bounded block length fits u64")
        {
            return Err(AlignmentError::CombinedIndex);
        }
        projected_hits = projected_hits.saturating_add(matches.exact_hit_count());
        if projected_hits > maximum_located_hits {
            if maximum_edit_distance <= INITIAL_EDIT_DISTANCE && !incremental_fallback_requested {
                return Ok(ReadAlignmentMetrics::default());
            }
            // The exact blocks are useful proofs inside the already
            // bounded mate windows even when their whole-genome FM
            // intervals are too repetitive to locate.  Falling back to
            // a local rolling scan preserves the hit cap while avoiding
            // a false-negative return for repeat-rich mates.
            // Expand only the anchor's minimum verified edit tier. The
            // higher-distance alternatives remain available to the
            // ordinary capped-FM path, but expanding all of them locally
            // creates broad low-confidence repeat windows.
            prepare_best_distance_rescue_windows(
                rescue_windows,
                reference,
                anchors,
                rescuing_mate1,
                maximum_template_span,
            )?;
            for &window in rescue_windows.iter() {
                let contig = reference.contig_by_ordinal(window.contig_ordinal).ok_or(
                    AlignmentError::InvalidContigOrdinal {
                        ordinal: window.contig_ordinal,
                    },
                )?;
                append_local_flexible_proof_candidates(
                    read,
                    contig.sequence().bases(),
                    window,
                    maximum_edit_distance,
                    &mut workspace.candidate_nominals,
                );
            }
            let (_, metrics) = workspace.verify_candidates_with_budget(
                reference,
                read,
                ReadAlignmentMetrics::default(),
                maximum_edit_distance,
            )?;
            return Ok(metrics);
        }
        exact[ordinal] = Some((
            matches,
            u64::try_from(query_start).expect("bounded query start fits u64"),
            1_u8 << ordinal,
        ));
    }

    let conversion_pass = if rescuing_mate1 {
        ConversionPass::Original
    } else {
        ConversionPass::Complementary
    };
    let query_len = u64::try_from(read.len()).expect("bounded read length fits u64");
    let mut located_rows = 0_u64;
    for (matches, query_offset, proof_mask) in exact.into_iter().flatten() {
        let metrics = reference
            .visit_combined_seed(matches, query_offset, query_len, &mut |hit| {
                let Some(strand) = conversion_pass.relabel_combined_hit(hit.strand()) else {
                    return true;
                };
                let candidate = ReadCandidate {
                    contig_ordinal: hit.contig_ordinal(),
                    start: hit.start(),
                    strand,
                    // The exact block establishes the nominal start;
                    // the flexible d3 verifier covers every start and
                    // endpoint displacement around it directly. This
                    // avoids constructing the whole-reference local
                    // filter planes for a sparse mate-rescue frontier.
                    proof_mask: FLEXIBLE_NOMINAL_PROOF | proof_mask,
                };
                if rescue_window_contains_candidate(
                    rescue_windows,
                    candidate,
                    maximum_edit_distance,
                ) {
                    workspace.candidate_nominals.push(candidate);
                }
                true
            })
            .map_err(|_| AlignmentError::CombinedIndex)?;
        located_rows = located_rows.saturating_add(metrics.located_coordinates());
    }
    let (_, metrics) = workspace.verify_candidates_with_budget(
        reference,
        read,
        ReadAlignmentMetrics {
            located_rows,
            ..ReadAlignmentMetrics::default()
        },
        maximum_edit_distance,
    )?;
    Ok(metrics)
}

pub(super) fn prepare_rescue_windows(
    rescue_windows: &mut Vec<MateRescueWindow>,
    reference: &ReferenceIndex,
    anchors: &[ReadPlacement],
    rescuing_mate1: bool,
    maximum_template_span: u64,
) -> Result<(), AlignmentError> {
    rescue_windows.clear();
    for &anchor in anchors {
        let Some(strand) = counterpart_strand(anchor.strand(), rescuing_mate1) else {
            continue;
        };
        let Some(contig) = reference.contig_by_ordinal(anchor.contig_ordinal()) else {
            return Err(AlignmentError::InvalidContigOrdinal {
                ordinal: anchor.contig_ordinal(),
            });
        };
        let lower = anchor.end().saturating_sub(maximum_template_span);
        let upper = anchor
            .start()
            .saturating_add(maximum_template_span)
            .min(contig.sequence().len().saturating_sub(1));
        rescue_windows.push(MateRescueWindow {
            contig_ordinal: anchor.contig_ordinal(),
            strand,
            start: lower,
            end: upper,
        });
    }
    merge_overlapping_rescue_windows(rescue_windows);
    Ok(())
}

pub(super) fn prepare_best_distance_rescue_windows(
    rescue_windows: &mut Vec<MateRescueWindow>,
    reference: &ReferenceIndex,
    anchors: &[ReadPlacement],
    rescuing_mate1: bool,
    maximum_template_span: u64,
) -> Result<(), AlignmentError> {
    rescue_windows.clear();
    let Some(best_distance) = anchors.iter().map(|anchor| anchor.distance()).min() else {
        return Ok(());
    };
    if best_distance > 1
        || anchors
            .iter()
            .filter(|anchor| anchor.distance() == best_distance)
            .take(2)
            .count()
            != 1
    {
        return Ok(());
    }
    for &anchor in anchors
        .iter()
        .filter(|anchor| anchor.distance() == best_distance)
    {
        let Some(strand) = counterpart_strand(anchor.strand(), rescuing_mate1) else {
            continue;
        };
        let Some(contig) = reference.contig_by_ordinal(anchor.contig_ordinal()) else {
            return Err(AlignmentError::InvalidContigOrdinal {
                ordinal: anchor.contig_ordinal(),
            });
        };
        let lower = anchor.end().saturating_sub(maximum_template_span);
        let upper = anchor
            .start()
            .saturating_add(maximum_template_span)
            .min(contig.sequence().len().saturating_sub(1));
        rescue_windows.push(MateRescueWindow {
            contig_ordinal: anchor.contig_ordinal(),
            strand,
            start: lower,
            end: upper,
        });
    }
    merge_overlapping_rescue_windows(rescue_windows);
    Ok(())
}

pub(super) fn merge_overlapping_rescue_windows(rescue_windows: &mut Vec<MateRescueWindow>) {
    rescue_windows.sort_unstable();
    let mut retained = 0_usize;
    for index in 0..rescue_windows.len() {
        let incoming = rescue_windows[index];
        if retained != 0 {
            let previous = &mut rescue_windows[retained - 1];
            if previous.contig_ordinal == incoming.contig_ordinal
                && previous.strand == incoming.strand
                && incoming.start <= previous.end.saturating_add(1)
            {
                previous.end = previous.end.max(incoming.end);
                continue;
            }
        }
        rescue_windows[retained] = incoming;
        retained += 1;
    }
    rescue_windows.truncate(retained);
}

/// Adds candidate-local ungapped semi-global placements without another
/// FM-index traversal. The existing seed frontier supplies full-read
/// nominal origins; this pass only chooses bounded terminal endpoints.
pub(super) fn append_ungapped_semi_global_placements(
    workspace: &mut ReadWorkspace,
    reference: &ReferenceIndex,
    read: &[Base],
    maximum_edit_distance: u8,
    clip_penalty: u8,
) {
    for &candidate in &workspace.candidate_nominals {
        if let Some(placement) = best_ungapped_semi_global_placement(
            reference,
            read,
            candidate,
            maximum_edit_distance,
            clip_penalty,
        ) {
            workspace.placements.push(placement);
        }
    }
    workspace.placements.sort_unstable_by_key(|placement| {
        (
            placement.contig_ordinal,
            placement.strand,
            placement.start,
            placement.end,
            placement.distance,
            placement.query_start,
            placement.query_end,
            placement.fallback_score,
        )
    });
    workspace.placements.dedup();
}

pub(super) fn relabel_exact_retained_hit(
    hit: CombinedSeedHit,
    lane: usize,
) -> Option<ReadCandidate> {
    let conversion_pass = if lane == 1 {
        ConversionPass::Complementary
    } else {
        ConversionPass::Original
    };
    let strand = conversion_pass.relabel_combined_hit(hit.strand())?;
    Some(ReadCandidate {
        contig_ordinal: hit.contig_ordinal(),
        start: hit.start(),
        strand,
        proof_mask: 0,
    })
}

pub(super) fn exact_retained_placement(
    candidate: ReadCandidate,
    selected: ReadPlacement,
    retained_length: usize,
) -> Option<ReadPlacement> {
    Some(ReadPlacement {
        contig_ordinal: candidate.contig_ordinal(),
        start: candidate.start(),
        end: candidate
            .start()
            .checked_add(u64::try_from(retained_length).ok()?)?,
        strand: candidate.strand(),
        distance: 0,
        query_start: selected.query_start,
        query_end: selected.query_end,
        fallback_score: selected.fallback_score,
    })
}

pub(super) fn exact_compatible_pair(
    first: ReadPlacement,
    second: ReadPlacement,
    minimum_template_span: u64,
    maximum_template_span: u64,
) -> Option<PairedPlacement> {
    if expected_mate2_strand(first.strand()) != Some(second.strand())
        || first.contig_ordinal() != second.contig_ordinal()
    {
        return None;
    }
    let template_start = first.start().min(second.start());
    let template_end = first.end().max(second.end());
    let span = template_end.checked_sub(template_start)?;
    if !(minimum_template_span..=maximum_template_span).contains(&span) || !is_inward(first, second)
    {
        return None;
    }
    Some(PairedPlacement {
        mate1: first,
        mate2: second,
        template_start,
        template_end,
        distance: 0,
        score: first.fallback_score.saturating_add(second.fallback_score),
    })
}
