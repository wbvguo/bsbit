//! Paired mate-rescue policy and genomic-window completion.
//!
//! This module is a mechanical responsibility split from `paired_end`; it
//! shares the parent's private bounded-policy vocabulary and does not alter the
//! crate's public API or algorithmic ordering.

use crate::AlignmentError;
use crate::alignment_policy::RESCUE_BLOCKS;
use crate::placement::ReadPlacement;
use crate::read_mapping::{ReadAlignmentMetrics, ReadCandidate, ReadWorkspace};
use crate::read_mapping_limits::{INITIAL_EDIT_DISTANCE, MAX_EDIT_DISTANCE, MAX_READ_BASES};
use crate::search::combined_adaptive::{DIRECT_SINGLETON_PROOF, FLEXIBLE_NOMINAL_PROOF};
use crate::search::combined_query::{CombinedSearchReferenceExt, CombinedSeedMatches};
use bsbit_core::alphabet::Base;
use bsbit_core::bisulfite::{
    AlignmentOrientation, BisulfiteStrand, CytosineStrand, strand_semantics,
};
use bsbit_index::reference::ReferenceIndex;
use bsbit_index::storage::fm::{ProjectedBase, SearchBase};

use super::selection::counterpart_strand;

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub(super) struct MateRescueWindow {
    pub(super) contig_ordinal: u64,
    pub(super) strand: BisulfiteStrand,
    pub(super) start: u64,
    pub(super) end: u64,
}

/// One disjoint exact block used to prove a bounded mate-rescue frontier.
#[derive(Clone, Copy, Debug)]
struct RescueProofBlock {
    query_start: u16,
    query_end: u16,
}

impl RescueProofBlock {
    const fn query_start(self) -> u16 {
        self.query_start
    }

    const fn query_end(self) -> u16 {
        self.query_end
    }
}

/// Emits a complete flexible candidate proof frontier inside one mate window.
///
/// `maximum_edit_distance + 1` disjoint query blocks guarantee that every
/// in-budget gapped placement leaves at least one exact projected block. The
/// exact block can be displaced from the true alignment origin by at most the
/// edit budget, which is precisely the start domain covered by the flexible
/// verifier. This is used only after a whole-genome rescue interval exceeds
/// its locate cap, so work is proportional to the paired genomic window rather
/// than to the global repeat count.
pub(super) fn append_local_flexible_proof_candidates(
    read: &[Base],
    contig: &[Base],
    window: MateRescueWindow,
    maximum_edit_distance: u8,
    candidates: &mut Vec<ReadCandidate>,
) {
    let block_count = usize::from(maximum_edit_distance) + 1;
    debug_assert!(block_count <= usize::from(MAX_EDIT_DISTANCE) + 1);
    let short = read.len() / block_count;
    let long_count = read.len() % block_count;
    let semantics = strand_semantics(window.strand);
    let edit_budget = u64::from(maximum_edit_distance);
    let mut query_start = 0_usize;

    for block_ordinal in 0..block_count {
        let block_len = short + usize::from(block_ordinal < long_count);
        let query_end = query_start + block_len;
        let mut query_code = 0_u128;
        let mut canonical = true;
        for oriented_position in query_start..query_end {
            let query_base = match semantics.orientation() {
                AlignmentOrientation::Forward => read[oriented_position],
                AlignmentOrientation::Reverse => {
                    read[read.len() - oriented_position - 1].complement()
                }
            };
            let Some(code) = rescue_projected_code(query_base, semantics.cytosine_strand()) else {
                canonical = false;
                break;
            };
            query_code = (query_code << 2) | u128::from(code);
        }
        if !canonical {
            query_start = query_end;
            continue;
        }

        let query_offset = u64::try_from(query_start).expect("bounded query offset fits u64");
        let scan_start = window
            .start
            .saturating_add(query_offset)
            .saturating_sub(edit_budget);
        let scan_end = window
            .end
            .saturating_add(query_offset)
            .saturating_add(edit_budget);
        let Some(mut position) = usize::try_from(scan_start).ok() else {
            query_start = query_end;
            continue;
        };
        let Some(last) = usize::try_from(scan_end).ok() else {
            query_start = query_end;
            continue;
        };
        let last = last.min(contig.len().saturating_sub(block_len));
        if position > last || position.saturating_add(block_len) > contig.len() {
            query_start = query_end;
            continue;
        }
        let bits = block_len * 2;
        let mask = if bits == u128::BITS as usize {
            u128::MAX
        } else {
            (1_u128 << bits) - 1
        };
        let mut reference_code = pack_rescue_projection(
            &contig[position..position + block_len],
            semantics.cytosine_strand(),
        );
        loop {
            if reference_code == query_code {
                let observed = u64::try_from(position).expect("reference position fits u64");
                if let Some(nominal) = observed.checked_sub(query_offset) {
                    let extended_start = window.start.saturating_sub(edit_budget);
                    let extended_end = window.end.saturating_add(edit_budget);
                    if (extended_start..=extended_end).contains(&nominal) {
                        candidates.push(ReadCandidate {
                            contig_ordinal: window.contig_ordinal,
                            start: nominal,
                            strand: window.strand,
                            proof_mask: FLEXIBLE_NOMINAL_PROOF | (1_u16 << block_ordinal),
                        });
                    }
                }
            }
            if position == last {
                break;
            }
            position += 1;
            let incoming = rescue_projected_code(
                contig[position + block_len - 1],
                semantics.cytosine_strand(),
            )
            .unwrap_or(3);
            reference_code = ((reference_code << 2) & mask) | u128::from(incoming);
        }
        query_start = query_end;
    }
}

fn balanced_rescue_blocks(read_len: usize) -> [RescueProofBlock; RESCUE_BLOCKS] {
    let short = read_len / RESCUE_BLOCKS;
    let long_count = read_len % RESCUE_BLOCKS;
    let mut cursor = 0_usize;
    core::array::from_fn(|ordinal| {
        let length = short + usize::from(ordinal < long_count);
        let start = cursor;
        cursor += length;
        RescueProofBlock {
            query_start: u16::try_from(start).expect("bounded rescue block start fits u16"),
            query_end: u16::try_from(cursor).expect("bounded rescue block end fits u16"),
        }
    })
}

fn pack_rescue_projection(bases: &[Base], strand: CytosineStrand) -> u128 {
    bases.iter().fold(0_u128, |packed, &base| {
        (packed << 2) | u128::from(rescue_projected_code(base, strand).unwrap_or(3))
    })
}

const fn rescue_projected_code(base: Base, strand: CytosineStrand) -> Option<u8> {
    use CytosineStrand::{Bottom, Top};
    match (strand, base) {
        (_, Base::A) | (Bottom, Base::G) => Some(0),
        (Top, Base::C | Base::T) | (Bottom, Base::C) => Some(1),
        (Top, Base::G) | (Bottom, Base::T) => Some(2),
        _ => None,
    }
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
            1_u16 << ordinal,
        ));
    }

    let query_len = u64::try_from(read.len()).expect("bounded read length fits u64");
    let mut located_rows = 0_u64;
    for (matches, query_offset, proof_mask) in exact.into_iter().flatten() {
        let metrics = reference
            .visit_combined_seed(matches, query_offset, query_len, &mut |hit| {
                let strand = if rescuing_mate1 {
                    hit.strand()
                } else {
                    match hit.strand() {
                        BisulfiteStrand::OT => BisulfiteStrand::CTOT,
                        BisulfiteStrand::OB => BisulfiteStrand::CTOB,
                        BisulfiteStrand::CTOT | BisulfiteStrand::CTOB => return true,
                    }
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

fn prepare_rescue_windows(
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

fn prepare_best_distance_rescue_windows(
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

fn merge_overlapping_rescue_windows(rescue_windows: &mut Vec<MateRescueWindow>) {
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

pub(super) fn select_combined_window_rescue_anchor(
    first_seeds: [Option<CombinedSeedMatches>; 2],
    mate1: &[ReadCandidate],
    mate2: &[ReadCandidate],
) -> Option<usize> {
    let pools = [mate1, mate2];
    let evidence = core::array::from_fn::<_, 2, _>(|mate| {
        let seed = first_seeds[mate]?;
        if seed.exact_hit_count() != 1 || pools[mate].is_empty() {
            return None;
        }
        let direct = pools[mate]
            .iter()
            .any(|candidate| candidate.proof_mask & DIRECT_SINGLETON_PROOF != 0);
        Some((direct, seed.matched_bases(), pools[mate].len()))
    });
    match (evidence[0], evidence[1]) {
        (None, None) => None,
        (Some(_), None) => Some(0),
        (None, Some(_)) => Some(1),
        (Some(first), Some(second)) => {
            if first.0 != second.0 {
                Some(usize::from(!first.0))
            } else if first.1 != second.1 {
                Some(usize::from(first.1 < second.1))
            } else {
                Some(usize::from(first.2 > second.2))
            }
        }
    }
}

pub(super) fn nominal_pair_geometry_exists(
    mate1: &[ReadCandidate],
    mate2: &[ReadCandidate],
    read1_len: usize,
    read2_len: usize,
    maximum_span: u64,
    maximum_edit_distance: u8,
) -> bool {
    mate1.iter().any(|&left| {
        nominal_partner_exists(
            left,
            mate2,
            true,
            read1_len,
            read2_len,
            maximum_span,
            maximum_edit_distance,
        )
    })
}

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

pub(super) fn retain_nominal_pair_geometry(
    mate1: &mut Vec<ReadCandidate>,
    mate2: &mut Vec<ReadCandidate>,
    read1_len: usize,
    read2_len: usize,
    maximum_span: u64,
    maximum_edit_distance: u8,
) {
    mate1.retain(|left| {
        nominal_partner_exists(
            *left,
            mate2,
            true,
            read1_len,
            read2_len,
            maximum_span,
            maximum_edit_distance,
        )
    });
    mate2.retain(|right| {
        nominal_partner_exists(
            *right,
            mate1,
            false,
            read1_len,
            read2_len,
            maximum_span,
            maximum_edit_distance,
        )
    });
}

fn nominal_partner_exists(
    candidate: ReadCandidate,
    pool: &[ReadCandidate],
    candidate_is_mate1: bool,
    read1_len: usize,
    read2_len: usize,
    maximum_span: u64,
    maximum_edit_distance: u8,
) -> bool {
    let edit_budget = u64::from(maximum_edit_distance);
    let target = match (candidate_is_mate1, candidate.strand()) {
        (true, BisulfiteStrand::OT) => BisulfiteStrand::CTOT,
        (true, BisulfiteStrand::OB) => BisulfiteStrand::CTOB,
        (false, BisulfiteStrand::CTOT) => BisulfiteStrand::OT,
        (false, BisulfiteStrand::CTOB) => BisulfiteStrand::OB,
        _ => return false,
    };
    let lower_start = candidate.start().saturating_sub(maximum_span);
    let upper_start = candidate.start().saturating_add(maximum_span);
    let lower = pool.partition_point(|partner| {
        (partner.strand(), partner.contig_ordinal(), partner.start())
            < (target, candidate.contig_ordinal(), lower_start)
    });
    let upper = pool.partition_point(|partner| {
        (partner.strand(), partner.contig_ordinal(), partner.start())
            <= (target, candidate.contig_ordinal(), upper_start)
    });
    pool[lower..upper].iter().any(|partner| {
        let (left, right) = if candidate_is_mate1 {
            (candidate, *partner)
        } else {
            (*partner, candidate)
        };
        match (left.strand(), right.strand()) {
            (BisulfiteStrand::OT, BisulfiteStrand::CTOT) => {
                left.start()
                    < right
                        .start()
                        .saturating_add(u64::try_from(read2_len).unwrap_or(u64::MAX) + edit_budget)
            }
            (BisulfiteStrand::OB, BisulfiteStrand::CTOB) => {
                right.start()
                    < left
                        .start()
                        .saturating_add(u64::try_from(read1_len).unwrap_or(u64::MAX) + edit_budget)
            }
            _ => false,
        }
    })
}
