//! Canonical paired-end read-to-reference alignment.
//!
//! Information-first maximal-suffix seeds feed an integer-only candidate path
//! and a worker-owned edit-distance-three verifier.

use crate::AlignmentError;
use crate::placement::{
    ReadPlacement, SEMI_GLOBAL_EDIT_PENALTY, placement_net_gap_bases, placement_origin_key,
};
#[cfg(test)]
use crate::read_mapping::{
    LOCAL_FILTER_BLOCKS, LocalCandidateFilter, PlacementVerifier, VerificationCacheEntry,
};
use crate::read_mapping::{
    ReadAlignmentMetrics, ReadCandidate, ReadWorkspace, sort_nominal_candidates, strand_index,
};
use crate::read_mapping_limits::{
    INITIAL_EDIT_DISTANCE, MAX_EDIT_DISTANCE, MAX_READ_BASES, MIN_SUFFIX_BASES,
};
use crate::verification::affine::{AffineScoreWorkspace, banded_affine_score};
use crate::verification::ungapped::UngappedEndpoint;
use crate::verification::ungapped::{BoundedSemiglobalConfig, UngappedProfile};

#[cfg(test)]
use self::mapq::{
    PARSIMONY_MAX_LOCATED_ROWS as SENSITIVE_PARSIMONY_MAX_LOCATED_ROWS,
    PARSIMONY_REQUIRED_SCORE_GAP as SENSITIVE_PARSIMONY_REQUIRED_SCORE_GAP,
};
pub use self::mapq::{SENSITIVE_MAPQ_REPEAT_RISK_ROWS, bwa_pair_mapping_quality_from_evidence};
use self::mapq::{
    ambiguity_q10_certified as sensitive_ambiguity_q10_certified,
    effective_mapping_quality as sensitive_effective_mapping_quality,
    incomplete_sparse_completion_required as sensitive_incomplete_sparse_completion_required,
    paired_mapping_quality, stable_rescue_q20_certified as sensitive_stable_rescue_q20_certified,
    two_way_parsimony_q20_certified as sensitive_two_way_parsimony_q20_certified,
};
use crate::library::PairedLibraryProfile;
use crate::search::combined_adaptive::{
    CombinedSearchLimits, CombinedTwoLaneSearchState, DIRECT_SINGLETON_PROOF,
    FLEXIBLE_NOMINAL_PROOF, INITIAL_SEARCH_LIMITS, continue_combined_two_lane_search,
    prepare_combined_projection, prepare_combined_search_projection,
    start_combined_two_lane_search,
};
use crate::search::combined_query::{
    CombinedSearchReferenceExt, CombinedSeedHit, CombinedSeedMatches,
};
use bsbit_core::alphabet::Base;
use bsbit_core::bisulfite::{AlignmentOrientation, BisulfiteStrand, strand_semantics};
use bsbit_index::reference::ReferenceIndex;
use bsbit_index::storage::fm::{ProjectedBase, SearchBase};

mod mapq;
mod options;
mod result;

pub use options::{PairedAlignmentOptions, PairedSearchMode};
pub use result::{PairMappingStatus, PairedAlignmentResult, PairedPlacement};

/// Maximum number of paired reads mapped in one paired-end worker batch.
pub const PAIRED_ALIGNMENT_BATCH_SIZE: usize = 32;
/// Fixed edit-distance budget used by every supported paired-end search mode.
pub const PAIRED_MAX_EDIT_DISTANCE: u8 = MAX_EDIT_DISTANCE;
/// Largest per-mate edit-distance budget supported by the paired-end mapper.
/// Distances four and five use the generic narrow-band AVX2 kernel only in an
/// incremental fallback; the distance-three specialization is retained for
/// the common first pass inside both supported modes.
// Four disjoint blocks make the local mate-rescue candidate frontier complete
// for the paired-end edit-distance-three budget: at most three edits can
// disturb at most three blocks, leaving at least one exact proof block.
const RESCUE_BLOCKS: usize = INITIAL_EDIT_DISTANCE as usize + 1;
const SENSITIVE_RANKED_BLOCK_HITS: u64 = 512;
const SENSITIVE_UNMAPPED_RANKED_BLOCK_HITS: u64 = SENSITIVE_RANKED_BLOCK_HITS.saturating_mul(2);
const SENSITIVE_SELECTIVE_UNMAPPED_RANKED_BLOCK_HITS: u64 =
    SENSITIVE_RANKED_BLOCK_HITS.saturating_mul(8);
const SENSITIVE_SELECTIVE_UNMAPPED_MIN_RETAINED_HITS: u64 = 32;
const SENSITIVE_POSITIVE_MAPQ_REPORTING_MIN_RETAINED_HITS: u64 = 128;
const SENSITIVE_SELECTIVE_UNMAPPED_MAX_RETAINED_HITS: u64 = 2_049;
const SENSITIVE_POSITIVE_MAPQ_REPORTING_MAX_RETAINED_HITS: u64 = 512;
/// Candidate-row pressure that triggers a bounded second-best completion.
/// This starts below the MAPQ risk threshold so search can prove away false
/// uniqueness before the reporting layer needs to lower confidence.
const SENSITIVE_REPEAT_RECHECK_ROWS: u64 = 256;
const SENSITIVE_PROOF_BLOCKS: usize = MAX_EDIT_DISTANCE as usize + 1;
const SENSITIVE_ADAPTIVE_MIN_BLOCK_BASES: usize = 19;
const SENSITIVE_BALANCED_BOUNDARY_SHIFTS: [i8; SENSITIVE_PROOF_BLOCKS - 1] =
    [0; SENSITIVE_PROOF_BLOCKS - 1];
const SENSITIVE_ADAPTIVE_BOUNDARY_SHIFTS: [i8; 3] = [-3, 0, 3];
// The qualified endpoint search is bounded to a 30-base terminal domain so it
// remains a candidate-local operation rather than an unrestricted local
// aligner.
const SEMI_GLOBAL_MAX_CLIP_BASES: usize = 30;
const SEMI_GLOBAL_MIN_ALIGNED_BASES: usize = 50;
const ADAPTER_STABILITY_DELTA: usize = 8;
const SEMI_GLOBAL_ADMISSION_EDIT_PENALTY: u8 = 2;
const SEMI_GLOBAL_CLIP_PENALTY: u8 = 1;
// Endpoint representation is selected independently from genomic-locus
// ranking in the origin-grouped policy. Unsupported clipping must not win
// merely because several sequencing errors happen to be terminal: its affine
// extension equals the mismatch penalty. Explicit adapter evidence receives a
// separate, favorable clipping prior below.
const ORIGIN_ENDPOINT_CLIP_OPEN_PENALTY: u16 = 8;
const ORIGIN_ENDPOINT_CLIP_EXTENSION_PENALTY: u16 = 7;
const ORIGIN_ENDPOINT_ADAPTER_CLIP_OPEN_PENALTY: u16 = 2;
const ORIGIN_ENDPOINT_ADAPTER_CLIP_EXTENSION_PENALTY: u16 = 0;
const ORIGIN_ENDPOINT_MIN_ADAPTER_SUPPORT: usize = 8;
const ILLUMINA_UNIVERSAL_ADAPTER: &[u8] = b"AGATCGGAAGAGC";
const SENSITIVE_CLIP_PENALTY: u8 = 4;
// A complete terminal mismatch costs seven. Eight makes endpoint selection
// prefer the full-read placement over removing that mismatch;
// seven would still prefer the clipped placement through the retained-edit
// tie-break in pair selection.
const SENSITIVE_MIN_EVENT_PENALTY: u8 = if SENSITIVE_CLIP_PENALTY < SEMI_GLOBAL_EDIT_PENALTY {
    SENSITIVE_CLIP_PENALTY
} else {
    SEMI_GLOBAL_EDIT_PENALTY
};
const SEMI_GLOBAL_MAX_EXACT_ANCHOR_HITS: u64 = 256;

// BWA-MEM-compatible score units used only by the residual sensitive
// selector.  Candidate discovery and the qualified d3/d5 verifier continue to
// use conversion-aware edit distance, so the common path does not pay dynamic
// programming cost.
const BWA_MATCH_SCORE: i16 = 1;
const BWA_MISMATCH_PENALTY: i16 = 4;
const BWA_NEAR_SUBOPTIMAL_DELTA: i16 = 7;

#[derive(Clone, Copy)]
struct RankedBlockSeed {
    matches: CombinedSeedMatches,
    query_offset: u64,
    proof_mask: u8,
}

#[derive(Clone, Copy)]
struct RankedBlockSelection {
    retained_hits: u64,
    complete: bool,
}

type RankedBlockPartition = Option<(u64, [Option<RankedBlockSeed>; SENSITIVE_PROOF_BLOCKS])>;
type EndpointKey = (u16, usize, u8, usize, usize, u8, usize, usize);

fn selective_unmapped_frontier_deepening_required(
    selections: [Option<RankedBlockSelection>; 2],
    ordinary_frontier_metrics: Option<PairAlignmentMetrics>,
) -> bool {
    let [Some(first), Some(second)] = selections else {
        return false;
    };
    let retained_hits = first.retained_hits.saturating_add(second.retained_hits);
    let both_incomplete = !first.complete && !second.complete;
    let inside_bounded_window = both_incomplete
        && (SENSITIVE_SELECTIVE_UNMAPPED_MIN_RETAINED_HITS
            ..SENSITIVE_SELECTIVE_UNMAPPED_MAX_RETAINED_HITS)
            .contains(&retained_hits);
    {
        let _ = ordinary_frontier_metrics;
        inside_bounded_window
    }
}

/// One disjoint exact block and its already-counted FM interval.
#[derive(Clone, Copy, Debug)]
struct ProofBlock {
    query_start: u16,
    query_end: u16,
}

impl ProofBlock {
    #[must_use]
    const fn query_start(self) -> u16 {
        self.query_start
    }

    #[must_use]
    const fn query_end(self) -> u16 {
        self.query_end
    }
}

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
struct MateRescueWindow {
    contig_ordinal: u64,
    strand: BisulfiteStrand,
    start: u64,
    end: u64,
}

/// Completes the missing-mate frontier inside every window induced by a
/// fully enumerated anchor frontier. Unlike the initial rescue path, the
/// block count follows the requested edit budget.
#[allow(clippy::too_many_arguments)]
fn rescue_from_ranked_anchor_windows(
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
fn rescue_from_combined_exact_blocks(
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

/// Adds candidate-local ungapped semi-global placements without another
/// FM-index traversal. The existing seed frontier supplies full-read
/// nominal origins; this pass only chooses bounded terminal endpoints.
fn append_ungapped_semi_global_placements(
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

fn relabel_exact_retained_hit(hit: CombinedSeedHit, lane: usize) -> Option<ReadCandidate> {
    let strand = if lane == 1 {
        match hit.strand() {
            BisulfiteStrand::OT => BisulfiteStrand::CTOT,
            BisulfiteStrand::OB => BisulfiteStrand::CTOB,
            BisulfiteStrand::CTOT | BisulfiteStrand::CTOB => return None,
        }
    } else {
        hit.strand()
    };
    Some(ReadCandidate {
        contig_ordinal: hit.contig_ordinal(),
        start: hit.start(),
        strand,
        proof_mask: 0,
    })
}

fn exact_retained_placement(
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

fn exact_compatible_pair(
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

type OriginPairStorageKey = ((u64, u8, i128), (u64, u8, i128));

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct OriginPairEvidence {
    // Larger is better. Each biological origin contributes only its best
    // endpoint score to MAPQ, while raw endpoint statistics remain available
    // to control the search pipeline.
    mapq_score: i16,
}

fn candidate_for_origin_endpoint(
    placement: ReadPlacement,
    read_length: usize,
) -> Option<ReadCandidate> {
    let (_, _, five_prime) = placement_origin_key(placement, read_length);
    let nominal_start = match strand_semantics(placement.strand()).orientation() {
        AlignmentOrientation::Forward => five_prime,
        AlignmentOrientation::Reverse => five_prime.checked_sub(
            i128::try_from(read_length.saturating_sub(1)).expect("bounded read length fits i128"),
        )?,
    };
    Some(ReadCandidate {
        contig_ordinal: placement.contig_ordinal(),
        start: u64::try_from(nominal_start).ok()?,
        strand: placement.strand(),
        proof_mask: 0,
    })
}

fn origin_endpoint_variant(
    reference: &ReferenceIndex,
    read: &[Base],
    placement: ReadPlacement,
    maximum_edit_distance: u8,
) -> ReadPlacement {
    if !placement.is_soft_clipped(read.len()) && !read_has_supported_three_prime_adapter(read) {
        return placement;
    }
    let Some(candidate) = candidate_for_origin_endpoint(placement, read.len()) else {
        return placement;
    };
    let endpoint_edit_limit =
        maximum_edit_distance.max(u8::try_from(read.len() / 5).unwrap_or(u8::MAX));
    let Some(endpoint) = best_ungapped_origin_endpoint_placement(
        reference,
        read,
        candidate,
        endpoint_edit_limit,
        SENSITIVE_CLIP_PENALTY,
    ) else {
        return placement;
    };
    if placement_origin_key(endpoint, read.len()) != placement_origin_key(placement, read.len()) {
        return placement;
    }
    if placement_endpoint_cost(read, endpoint) < placement_endpoint_cost(read, placement) {
        endpoint
    } else {
        placement
    }
}

/// Chooses the reported endpoint/CIGAR inside an already selected biological
/// locus. Mapping rank, ambiguity, and MAPQ must be frozen before this runs.
/// The endpoint objective uses conversion-aware edit evidence, affine terminal
/// clipping, and explicit 3-prime adapter support. It cannot move either mate
/// to a different five-prime origin.
#[must_use]
fn select_reported_origin_endpoint(
    reference: &ReferenceIndex,
    reads: [&[Base]; 2],
    selected: PairedPlacement,
    maximum_edit_distance: u8,
    minimum_template_span: u64,
    maximum_template_span: u64,
) -> PairedPlacement {
    let mate1_variant =
        origin_endpoint_variant(reference, reads[0], selected.mate1(), maximum_edit_distance);
    let mate2_variant =
        origin_endpoint_variant(reference, reads[1], selected.mate2(), maximum_edit_distance);
    // Almost every selected pair is already a whole-read endpoint and has no
    // adapter evidence. Avoid constructing and rescoring four identical pair
    // combinations on that common path.
    if mate1_variant == selected.mate1() && mate2_variant == selected.mate2() {
        return selected;
    }
    let alternatives = [
        [selected.mate1(), mate1_variant],
        [selected.mate2(), mate2_variant],
    ];
    let selected_origin = pair_origin_key(selected, reads[0].len(), reads[1].len());
    let mut best = selected;
    let mut best_key = pair_endpoint_key(reads, selected);
    for mate1 in alternatives[0] {
        for mate2 in alternatives[1] {
            let template_start = mate1.start().min(mate2.start());
            let template_end = mate1.end().max(mate2.end());
            let span = template_end.saturating_sub(template_start);
            if mate1.contig_ordinal() != mate2.contig_ordinal()
                || expected_mate2_strand(mate1.strand()) != Some(mate2.strand())
                || !(minimum_template_span..=maximum_template_span).contains(&span)
                || !is_inward(mate1, mate2)
            {
                continue;
            }
            let candidate = PairedPlacement {
                mate1,
                mate2,
                template_start,
                template_end,
                distance: mate1.distance().saturating_add(mate2.distance()),
                // Preserve the score that selected the locus. Endpoint choice
                // is downstream of confidence and cannot rerank candidates.
                score: selected.score(),
            };
            if pair_origin_key(candidate, reads[0].len(), reads[1].len()) != selected_origin {
                continue;
            }
            let key = pair_endpoint_key(reads, candidate);
            if key < best_key {
                best = candidate;
                best_key = key;
            }
        }
    }
    best
}

// This ordering is applied only after pair selection has retained an ambiguous
// best-score tie. It chooses the BAM representative without removing a tied
// placement or changing the combined-index pair class.
fn prefer_minimum_net_gap_representative(
    pairs: &mut [PairedPlacement],
    read1_len: usize,
    read2_len: usize,
) {
    pairs.sort_unstable_by_key(|pair| {
        (
            placement_net_gap_bases(pair.mate1(), read1_len)
                .saturating_add(placement_net_gap_bases(pair.mate2(), read2_len)),
            *pair,
        )
    });
}

fn pair_net_gap_profile(
    pairs: &[PairedPlacement],
    read1_len: usize,
    read2_len: usize,
) -> (Option<u64>, Option<u64>, u64) {
    let mut minimum = None;
    let mut second = None;
    let mut minimum_count = 0_u64;
    for pair in pairs {
        let gap = placement_net_gap_bases(pair.mate1(), read1_len)
            .saturating_add(placement_net_gap_bases(pair.mate2(), read2_len));
        match minimum {
            None => {
                minimum = Some(gap);
                minimum_count = 1;
            }
            Some(current) if gap < current => {
                second = minimum;
                minimum = Some(gap);
                minimum_count = 1;
            }
            Some(current) if gap == current => {
                minimum_count = minimum_count.saturating_add(1);
            }
            Some(_) if second.is_none_or(|current| gap < current) => {
                second = Some(gap);
            }
            Some(_) => {}
        }
    }
    (minimum, second, minimum_count)
}

/// Result counters for one paired-end alignment.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct PairAlignmentMetrics {
    pub(crate) mate1: ReadAlignmentMetrics,
    pub(crate) mate2: ReadAlignmentMetrics,
    pub(crate) compatible_pairs: u64,
    pub(crate) best_pair_placements: u64,
    pub(crate) window_rescue_attempted: bool,
    pub(crate) semi_global_attempted: bool,
    /// Best compatible pair score in BWA score units (larger is better).
    pub(crate) best_pair_score: Option<i16>,
    /// Best strictly lower compatible pair score, when one was observed.
    pub(crate) second_best_pair_score: Option<i16>,
    /// Number of alternative pairings within the BWA near-suboptimal window.
    pub(crate) near_best_pairings: u64,
    /// Confidence evidence collapsed to distinct biological pair origins.
    /// These fields are consumed only by MAPQ; raw fields above continue to
    /// control affine rescoring, rescue, and candidate-search decisions.
    pub(crate) mapq_compatible_pairs: u64,
    pub(crate) mapq_best_pair_score: Option<i16>,
    pub(crate) mapq_second_best_pair_score: Option<i16>,
    pub(crate) mapq_near_best_pairings: u64,
    /// Whether all candidate work required by the active bounded search ended.
    pub(crate) frontier_complete: bool,
}

/// Reusable state for one paired-end mapping worker.
struct PairWorkspace {
    mate1: ReadWorkspace,
    mate2: ReadWorkspace,
    rescue_windows: Vec<MateRescueWindow>,
    best_pairs: Vec<PairedPlacement>,
    exact_anchor_candidates: Vec<ReadCandidate>,
    ranked_anchor_placements: Vec<ReadPlacement>,
    mate1_affine_scores: Vec<i16>,
    mate2_affine_scores: Vec<i16>,
    affine: AffineScoreWorkspace,
    semi_global_clip_penalty: u8,
    prefer_minimum_net_gap: bool,
    origin_pair_evidence: std::collections::HashMap<OriginPairStorageKey, OriginPairEvidence>,
    combined_search_state: CombinedTwoLaneSearchState,
    fallback_mate1_nominals: Vec<ReadCandidate>,
    fallback_mate2_nominals: Vec<ReadCandidate>,
    ranked_extension_selections: [Option<RankedBlockSelection>; 2],
}

/// One copied result from a cross-pair combined first-seed wavefront.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
// These booleans are independent output certificates, not interchangeable
// state flags; a bitfield would make their meanings less explicit.
#[allow(clippy::struct_excessive_bools)]
pub(crate) struct PairedBatchResult {
    class: PairMappingStatus,
    metrics: PairAlignmentMetrics,
    best_pair: Option<PairedPlacement>,
    second_best_distance: Option<u8>,
    repeat_risk_q20_certified: bool,
    parsimony_q20_certified: bool,
    ambiguity_q10_certified: bool,
    requires_positive_mapq_for_reporting: bool,
}

impl PairedBatchResult {
    #[must_use]
    pub(crate) const fn class(self) -> PairMappingStatus {
        self.class
    }

    #[must_use]
    pub(crate) const fn metrics(self) -> PairAlignmentMetrics {
        self.metrics
    }

    #[must_use]
    pub(crate) const fn best_pair(self) -> Option<PairedPlacement> {
        self.best_pair
    }

    /// Returns the selected pair score in BWA score units (larger is better).
    #[must_use]
    #[cfg(test)]
    pub(crate) const fn best_pair_score(self) -> Option<i16> {
        self.metrics.best_pair_score
    }

    /// Returns the best strictly lower pair score, when observed.
    #[must_use]
    #[cfg(test)]
    pub(crate) const fn second_best_pair_score(self) -> Option<i16> {
        self.metrics.second_best_pair_score
    }

    /// Returns the number of alternative pairings close enough to penalize
    /// BWA-style mapping quality.
    #[must_use]
    #[cfg(test)]
    pub(crate) const fn near_best_pairings(self) -> u64 {
        self.metrics.near_best_pairings
    }

    /// Returns the pair-level BWA-style score-gap mapping quality before any
    /// reporting-layer cap for clipping or repeat-risk provenance.
    #[must_use]
    pub(crate) fn evidence_mapping_quality(self) -> u8 {
        bwa_pair_mapping_quality_from_evidence(
            self.class,
            self.metrics.frontier_complete,
            self.metrics.best_pair_score,
            self.metrics.second_best_pair_score,
            self.metrics.near_best_pairings,
        )
    }

    /// Reports whether an independent, bounded endpoint pass supplied enough
    /// stable one-mate-rescue evidence to clear only the MAPQ-20 boundary.
    #[must_use]
    pub(crate) const fn repeat_risk_q20_certified(self) -> bool {
        self.repeat_risk_q20_certified
    }

    /// Reports whether a complete two-way semi-global tie had one uniquely
    /// parsimonious representative inside the qualified Q20 envelope.
    #[must_use]
    pub(crate) const fn parsimony_q20_certified(self) -> bool {
        self.parsimony_q20_certified
    }
}

#[derive(Clone, Copy)]
struct AdapterFallbackResult {
    result: PairedBatchResult,
    stability_result: Option<PairedBatchResult>,
    final_class: PairMappingStatus,
    retained_bases: [usize; 2],
}

/// Worker-owned storage for one combined cross-read-pair seed wavefront.
pub struct PairedBatchAligner {
    pair: PairWorkspace,
    projections: Vec<[[ProjectedBase; MAX_READ_BASES]; 2]>,
    first_seeds: Vec<Option<CombinedSeedMatches>>,
    results: Vec<PairedBatchResult>,
}

impl PairedBatchAligner {
    /// Allocates reusable mapping storage for at least `pair_capacity` pairs.
    #[must_use]
    pub fn with_capacity(pair_capacity: usize) -> Self {
        Self {
            pair: PairWorkspace::with_capacity(4096, 1024, 32),
            projections: Vec::with_capacity(pair_capacity),
            first_seeds: Vec::with_capacity(pair_capacity.saturating_mul(2)),
            results: Vec::with_capacity(pair_capacity),
        }
    }

    /// Maps a directional or non-directional paired-end batch with one of the
    /// qualified paired-end strategies.
    ///
    /// # Errors
    ///
    /// Returns [`AlignmentError`] for an unsupported library profile,
    /// invalid template spans, or any mapping failure.
    fn map_pairs_combined<'a>(
        &'a mut self,
        reference: &ReferenceIndex,
        reads: &[[&[Base]; 2]],
        options: PairedAlignmentOptions,
    ) -> Result<&'a [PairedBatchResult], AlignmentError> {
        let (maximum_edit_distance, window_rescue, semi_global) = options.derived_policy();
        match options.library_profile {
            PairedLibraryProfile::Directional => self.map_directional_pairs_combined_inner(
                reference,
                reads,
                maximum_edit_distance,
                options.minimum_template_span,
                options.maximum_template_span,
                window_rescue,
                semi_global,
                options.search_mode,
            ),
            PairedLibraryProfile::NonDirectional => self
                .map_non_directional_pairs_combined_with_search_mode(
                    reference,
                    reads,
                    maximum_edit_distance,
                    options.minimum_template_span,
                    options.maximum_template_span,
                    window_rescue,
                    semi_global,
                    options.search_mode,
                ),
        }
    }

    /// Maps a paired-read batch through the complete qualified output policy.
    ///
    /// Adapter-supported trimming, stability remapping, MAPQ certificates,
    /// positive-MAPQ admission, and endpoint representation are resolved here
    /// so serialization callers receive facts rather than policy controls.
    ///
    /// # Errors
    ///
    /// Returns [`AlignmentError`] when the template span is invalid
    /// or any qualified mapping phase fails.
    ///
    /// # Panics
    ///
    /// Panics only if internally generated adapter-stability metadata loses
    /// its matching adapter result, which would violate this method's local
    /// construction invariant.
    // Adapter repair, stability proof, MAPQ certification, and reporting
    // admission form one ordered output-policy transaction.
    #[allow(clippy::too_many_lines)]
    pub fn map_pairs_for_output(
        &mut self,
        reference: &ReferenceIndex,
        reads: &[[&[Base]; 2]],
        options: PairedAlignmentOptions,
    ) -> Result<Vec<PairedAlignmentResult>, AlignmentError> {
        let primary = self.map_pairs_combined(reference, reads, options)?.to_vec();
        let mut adapter_results = vec![None; reads.len()];
        let mut adapter_classes = vec![None; reads.len()];
        let mut adapter_attempted = vec![false; reads.len()];
        let mut adapter_clipped_mates = vec![0_u8; reads.len()];
        let mut adapter_clipped_bases = vec![0_usize; reads.len()];
        let mut clipped_reads = Vec::with_capacity(reads.len());
        let mut clipped_metadata = Vec::with_capacity(reads.len());

        for (offset, (pair, result)) in reads.iter().zip(&primary).enumerate() {
            let should_attempt = matches!(result.class(), PairMappingStatus::Unmapped)
                || (options.search_mode.is_sensitive()
                    && result.metrics().window_rescue_attempted
                    && matches!(result.class(), PairMappingStatus::Ambiguous));
            if !should_attempt {
                continue;
            }
            let retained = [
                supported_three_prime_adapter_start(pair[0])
                    .filter(|&start| start >= SEMI_GLOBAL_MIN_ALIGNED_BASES)
                    .unwrap_or(pair[0].len()),
                supported_three_prime_adapter_start(pair[1])
                    .filter(|&start| start >= SEMI_GLOBAL_MIN_ALIGNED_BASES)
                    .unwrap_or(pair[1].len()),
            ];
            if retained == [pair[0].len(), pair[1].len()] {
                continue;
            }
            adapter_attempted[offset] = true;
            adapter_clipped_mates[offset] = u8::from(retained[0] != pair[0].len())
                .saturating_add(u8::from(retained[1] != pair[1].len()));
            adapter_clipped_bases[offset] = pair[0]
                .len()
                .saturating_sub(retained[0])
                .saturating_add(pair[1].len().saturating_sub(retained[1]));
            clipped_reads.push([&pair[0][..retained[0]], &pair[1][..retained[1]]]);
            clipped_metadata.push((offset, retained));
        }

        if !clipped_reads.is_empty() {
            let adapter_options = PairedAlignmentOptions::adapter_trimmed(
                options.library_profile,
                options.search_mode,
                options.minimum_template_span,
                options.maximum_template_span,
            );
            let remapped = self
                .map_pairs_combined(reference, &clipped_reads, adapter_options)?
                .to_vec();
            for ((offset, retained_bases), result) in clipped_metadata.iter().copied().zip(remapped)
            {
                adapter_results[offset] = Some(AdapterFallbackResult {
                    result,
                    stability_result: None,
                    final_class: result.class(),
                    retained_bases,
                });
            }

            let mut stability_reads = Vec::with_capacity(clipped_reads.len());
            let mut stability_metadata = Vec::with_capacity(clipped_reads.len());
            for (offset, fallback) in adapter_results.iter_mut().enumerate() {
                let Some(fallback) = fallback else {
                    continue;
                };
                if !matches!(fallback.final_class, PairMappingStatus::Unique) {
                    continue;
                }
                let full_lengths = [reads[offset][0].len(), reads[offset][1].len()];
                let mut retained = fallback.retained_bases;
                let mut stable_domain = true;
                for mate in 0..2 {
                    if retained[mate] == full_lengths[mate] {
                        continue;
                    }
                    if retained[mate]
                        < SEMI_GLOBAL_MIN_ALIGNED_BASES.saturating_add(ADAPTER_STABILITY_DELTA)
                    {
                        stable_domain = false;
                        break;
                    }
                    retained[mate] -= ADAPTER_STABILITY_DELTA;
                }
                if !stable_domain {
                    fallback.final_class = PairMappingStatus::Ambiguous;
                    continue;
                }
                stability_reads.push([
                    &reads[offset][0][..retained[0]],
                    &reads[offset][1][..retained[1]],
                ]);
                stability_metadata.push(offset);
            }

            if !stability_reads.is_empty() {
                let stability = self
                    .map_pairs_combined(reference, &stability_reads, adapter_options)?
                    .to_vec();
                for (offset, stability_result) in stability_metadata.into_iter().zip(stability) {
                    let fallback = adapter_results[offset]
                        .as_mut()
                        .expect("stability metadata refers to an adapter result");
                    fallback.stability_result = Some(stability_result);
                    let same_origin = fallback
                        .result
                        .best_pair()
                        .zip(stability_result.best_pair())
                        .is_some_and(|(primary, stability)| {
                            pair_origin_key(
                                primary,
                                fallback.retained_bases[0],
                                fallback.retained_bases[1],
                            ) == pair_origin_key(
                                stability,
                                fallback.retained_bases[0],
                                fallback.retained_bases[1],
                            )
                        });
                    if !matches!(stability_result.class(), PairMappingStatus::Unique)
                        || !same_origin
                    {
                        fallback.final_class = PairMappingStatus::Ambiguous;
                    }
                }
            }

            for (offset, fallback) in adapter_results.iter_mut().enumerate() {
                let Some(candidate) = *fallback else {
                    continue;
                };
                if matches!(primary[offset].class(), PairMappingStatus::Ambiguous)
                    && primary[offset].metrics().window_rescue_attempted
                    && !matches!(candidate.final_class, PairMappingStatus::Unique)
                {
                    adapter_classes[offset] = Some(PairMappingStatus::Ambiguous);
                    *fallback = None;
                } else {
                    adapter_classes[offset] = Some(candidate.final_class);
                }
            }
        }

        let mut outputs = Vec::with_capacity(reads.len());
        for (offset, (pair, strict_result)) in reads.iter().zip(primary).enumerate() {
            let adapter = adapter_results[offset];
            let result = adapter.map_or(strict_result, |fallback| fallback.result);
            let mut class = adapter.map_or(result.class(), |fallback| fallback.final_class);
            let mate_rescue_attempted = result.metrics().window_rescue_attempted
                || adapter.is_some_and(|fallback| {
                    fallback
                        .stability_result
                        .is_some_and(|stability| stability.metrics().window_rescue_attempted)
                });
            let semi_global_attempted = result.metrics().semi_global_attempted;
            let mut semi_global_clipped_mates = 0_u8;
            let mut semi_global_clipped_bases = 0_usize;
            if semi_global_attempted && let Some(selected) = result.best_pair() {
                for (placement, read) in [(selected.mate1(), pair[0]), (selected.mate2(), pair[1])]
                {
                    let retained = placement.retained_query_interval(read.len());
                    let clipped = read
                        .len()
                        .saturating_sub(retained.end.saturating_sub(retained.start));
                    semi_global_clipped_mates =
                        semi_global_clipped_mates.saturating_add(u8::from(clipped != 0));
                    semi_global_clipped_bases = semi_global_clipped_bases.saturating_add(clipped);
                }
            }

            let report_ambiguous =
                matches!(class, PairMappingStatus::Ambiguous) && result.best_pair().is_some();
            let Some(selected) = result
                .best_pair()
                .filter(|_| matches!(class, PairMappingStatus::Unique) || report_ambiguous)
            else {
                outputs.push(PairedAlignmentResult {
                    class,
                    placement: None,
                    retained_query_intervals: [0..pair[0].len(), 0..pair[1].len()],
                    mapping_quality: 0,
                    adapter_attempted: adapter_attempted[offset],
                    adapter_class: adapter_classes[offset],
                    adapter_clipped_mates: adapter_clipped_mates[offset],
                    adapter_clipped_bases: adapter_clipped_bases[offset],
                    semi_global_attempted,
                    semi_global_clipped_mates,
                    semi_global_clipped_bases,
                    mate_rescue_attempted,
                });
                continue;
            };
            let mut retained_query_intervals = adapter.map_or_else(
                || {
                    [
                        selected.mate1().retained_query_interval(pair[0].len()),
                        selected.mate2().retained_query_interval(pair[1].len()),
                    ]
                },
                |fallback| [0..fallback.retained_bases[0], 0..fallback.retained_bases[1]],
            );
            let mapping_quality = paired_mapping_quality(
                result,
                adapter.and_then(|fallback| fallback.stability_result),
                class,
                options.search_mode,
                [pair[0].len(), pair[1].len()],
                [&retained_query_intervals[0], &retained_query_intervals[1]],
            );
            let requires_positive_mapq = strict_result.requires_positive_mapq_for_reporting
                || result.requires_positive_mapq_for_reporting;
            if requires_positive_mapq && mapping_quality == 0 {
                class = PairMappingStatus::Unmapped;
                outputs.push(PairedAlignmentResult {
                    class,
                    placement: None,
                    retained_query_intervals,
                    mapping_quality,
                    adapter_attempted: adapter_attempted[offset],
                    adapter_class: adapter_classes[offset],
                    adapter_clipped_mates: adapter_clipped_mates[offset],
                    adapter_clipped_bases: adapter_clipped_bases[offset],
                    semi_global_attempted,
                    semi_global_clipped_mates,
                    semi_global_clipped_bases,
                    mate_rescue_attempted,
                });
                continue;
            }
            let selected = if adapter.is_none() && options.search_mode.is_sensitive() {
                select_reported_origin_endpoint(
                    reference,
                    *pair,
                    selected,
                    PAIRED_MAX_EDIT_DISTANCE,
                    options.minimum_template_span,
                    options.maximum_template_span,
                )
            } else {
                selected
            };
            retained_query_intervals = adapter.map_or_else(
                || {
                    [
                        selected.mate1().retained_query_interval(pair[0].len()),
                        selected.mate2().retained_query_interval(pair[1].len()),
                    ]
                },
                |fallback| [0..fallback.retained_bases[0], 0..fallback.retained_bases[1]],
            );
            outputs.push(PairedAlignmentResult {
                class,
                placement: Some(selected),
                retained_query_intervals,
                mapping_quality,
                adapter_attempted: adapter_attempted[offset],
                adapter_class: adapter_classes[offset],
                adapter_clipped_mates: adapter_clipped_mates[offset],
                adapter_clipped_bases: adapter_clipped_bases[offset],
                semi_global_attempted,
                semi_global_clipped_mates,
                semi_global_clipped_bases,
                mate_rescue_attempted,
            });
        }
        Ok(outputs)
    }

    #[allow(clippy::too_many_arguments)]
    fn map_non_directional_pairs_combined_with_search_mode<'a>(
        &'a mut self,
        reference: &ReferenceIndex,
        reads: &[[&[Base]; 2]],
        maximum_edit_distance: u8,
        minimum_template_span: u64,
        maximum_template_span: u64,
        window_rescue: bool,
        semi_global: bool,
        search_mode: PairedSearchMode,
    ) -> Result<&'a [PairedBatchResult], AlignmentError> {
        let original = self
            .map_directional_pairs_combined_inner(
                reference,
                reads,
                maximum_edit_distance,
                minimum_template_span,
                maximum_template_span,
                window_rescue,
                semi_global,
                search_mode,
            )?
            .to_vec();
        let swapped_reads = reads
            .iter()
            .map(|pair| [pair[1], pair[0]])
            .collect::<Vec<_>>();
        let complementary = self
            .map_directional_pairs_combined_inner(
                reference,
                &swapped_reads,
                maximum_edit_distance,
                minimum_template_span,
                maximum_template_span,
                window_rescue,
                semi_global,
                search_mode,
            )?
            .iter()
            .copied()
            .map(swap_batch_result_mates)
            .collect::<Vec<_>>();
        self.results.clear();
        self.results.extend(original.iter().zip(&complementary).map(
            |(original, complementary)| {
                merge_non_directional_batch_results(original, complementary)
            },
        ));
        Ok(&self.results)
    }

    #[allow(clippy::too_many_arguments)]
    // Cross-pair wavefront seeding and per-pair completion intentionally share
    // one batch workspace so projected reads and seed states are not copied.
    #[allow(clippy::too_many_lines)]
    fn map_directional_pairs_combined_inner<'a>(
        &'a mut self,
        reference: &ReferenceIndex,
        reads: &[[&[Base]; 2]],
        maximum_edit_distance: u8,
        minimum_template_span: u64,
        maximum_template_span: u64,
        window_rescue: bool,
        semi_global: bool,
        search_mode: PairedSearchMode,
    ) -> Result<&'a [PairedBatchResult], AlignmentError> {
        if maximum_edit_distance > MAX_EDIT_DISTANCE {
            return Err(AlignmentError::UnsupportedEditDistance {
                requested: maximum_edit_distance,
                maximum: MAX_EDIT_DISTANCE,
            });
        }
        if minimum_template_span > maximum_template_span {
            return Err(AlignmentError::InvertedTemplateSpan {
                minimum: minimum_template_span,
                maximum: maximum_template_span,
            });
        }
        if reads.len().saturating_mul(2) > 64 {
            return Err(AlignmentError::LocatedCountOverflow);
        }
        self.projections.clear();
        self.projections
            .resize(reads.len(), [[ProjectedBase::A; MAX_READ_BASES]; 2]);
        for (projection, pair) in self.projections.iter_mut().zip(reads) {
            prepare_combined_projection(pair[0], false, &mut projection[0])?;
            prepare_combined_projection(pair[1], true, &mut projection[1])?;
        }
        let patterns = self
            .projections
            .iter()
            .zip(reads)
            .flat_map(|(projection, pair)| {
                [
                    &projection[0][..pair[0].len()],
                    &projection[1][..pair[1].len()],
                ]
            })
            .collect::<Vec<_>>();
        self.first_seeds = reference
            .combined_maximal_suffix_projected_wavefront(&patterns, MIN_SUFFIX_BASES)
            .map_err(|_| AlignmentError::CombinedIndex)?;
        self.pair.semi_global_clip_penalty = search_mode.semi_global_clip_penalty();
        self.pair.prefer_minimum_net_gap = search_mode.is_sensitive();
        // Sensitive mode uses semi-global alignment as a confidence repair,
        // not as an eager replacement objective. Run the proof-oriented
        // strict search first, preserve every already-high-confidence result,
        // and revisit only the small residual low-confidence frontier below.
        let eager_semi_global = semi_global && !search_mode.is_sensitive();
        self.results.clear();
        for (ordinal, pair) in reads.iter().enumerate() {
            let mut repeat_risk_q20_certified = false;
            let mut parsimony_q20_certified = false;
            let mut ambiguity_q10_certified = false;
            // The final paired-end frontier may append bounded completion work.
            let mut requires_positive_mapq_for_reporting = false;
            let first_seeds = [
                self.first_seeds[ordinal * 2],
                self.first_seeds[ordinal * 2 + 1],
            ];
            let projection = &self.projections[ordinal];
            let initial_edit_distance = maximum_edit_distance.min(INITIAL_EDIT_DISTANCE);
            let (mut class, mut metrics, mut second_best_distance) =
                self.pair.map_directional_pair_combined_prepared(
                    reference,
                    pair[0],
                    pair[1],
                    [
                        &projection[0][..pair[0].len()],
                        &projection[1][..pair[1].len()],
                    ],
                    first_seeds,
                    initial_edit_distance,
                    minimum_template_span,
                    maximum_template_span,
                    window_rescue,
                    eager_semi_global,
                    INITIAL_SEARCH_LIMITS,
                    true,
                )?;
            if matches!(class, PairMappingStatus::Unmapped) {
                if maximum_edit_distance > initial_edit_distance {
                    let reverified = self.pair.reverify_directional_pair_combined_candidates(
                        reference,
                        pair[0],
                        pair[1],
                        maximum_edit_distance,
                        minimum_template_span,
                        maximum_template_span,
                        eager_semi_global,
                        metrics,
                    )?;
                    if !matches!(reverified.0, PairMappingStatus::Unmapped) {
                        (class, metrics, second_best_distance) = reverified;
                    }
                }
                if matches!(class, PairMappingStatus::Unmapped) {
                    (class, metrics, second_best_distance) =
                        self.pair.continue_directional_pair_combined_incremental(
                            reference,
                            pair[0],
                            pair[1],
                            [
                                &projection[0][..pair[0].len()],
                                &projection[1][..pair[1].len()],
                            ],
                            first_seeds,
                            maximum_edit_distance,
                            minimum_template_span,
                            maximum_template_span,
                            window_rescue,
                            eager_semi_global,
                            metrics,
                        )?;
                }
            }
            let complete_suspicious_unique = sensitive_repeat_recheck_required(class, metrics);
            if search_mode.is_sensitive()
                && window_rescue
                && (matches!(class, PairMappingStatus::Unmapped) || complete_suspicious_unique)
            {
                let original = (class, metrics, second_best_distance);
                let original_best = complete_suspicious_unique
                    .then(|| self.pair.best_pairs().first().copied())
                    .flatten();
                let mut completed = self.pair.extend_directional_pair_from_ranked_blocks(
                    reference,
                    pair[0],
                    pair[1],
                    [
                        &projection[0][..pair[0].len()],
                        &projection[1][..pair[1].len()],
                    ],
                    maximum_edit_distance,
                    minimum_template_span,
                    maximum_template_span,
                    eager_semi_global,
                    SENSITIVE_RANKED_BLOCK_HITS,
                )?;
                if matches!(class, PairMappingStatus::Unmapped)
                    && matches!(
                        completed.as_ref().map(|candidate| candidate.0),
                        None | Some(PairMappingStatus::Unmapped)
                    )
                {
                    completed = self.pair.extend_directional_pair_from_ranked_blocks(
                        reference,
                        pair[0],
                        pair[1],
                        [
                            &projection[0][..pair[0].len()],
                            &projection[1][..pair[1].len()],
                        ],
                        maximum_edit_distance,
                        minimum_template_span,
                        maximum_template_span,
                        eager_semi_global,
                        SENSITIVE_UNMAPPED_RANKED_BLOCK_HITS,
                    )?;
                    if matches!(
                        completed.as_ref().map(|candidate| candidate.0),
                        None | Some(PairMappingStatus::Unmapped)
                    ) && self.pair.selective_unmapped_frontier_deepening_required(
                        completed.as_ref().map(|candidate| candidate.1),
                    ) {
                        let extended_frontier_requires_positive_mapq = self
                            .pair
                            .selective_unmapped_frontier_requires_positive_mapq();
                        completed = self.pair.extend_directional_pair_from_ranked_blocks(
                            reference,
                            pair[0],
                            pair[1],
                            [
                                &projection[0][..pair[0].len()],
                                &projection[1][..pair[1].len()],
                            ],
                            maximum_edit_distance,
                            minimum_template_span,
                            maximum_template_span,
                            eager_semi_global,
                            SENSITIVE_SELECTIVE_UNMAPPED_RANKED_BLOCK_HITS,
                        )?;
                        {
                            requires_positive_mapq_for_reporting =
                                completed.as_ref().is_some_and(|candidate| {
                                    !matches!(candidate.0, PairMappingStatus::Unmapped)
                                }) && extended_frontier_requires_positive_mapq;
                        }
                    }
                }
                match (original_best, completed) {
                    (_, Some(completed)) if !matches!(completed.0, PairMappingStatus::Unmapped) => {
                        let completed_best = self.pair.best_pairs().first().copied();
                        if original_best
                            .zip(completed_best)
                            .is_some_and(|(original, completed)| {
                                completed.score() > original.score()
                            })
                        {
                            self.pair.best_pairs.clear();
                            self.pair
                                .best_pairs
                                .push(original_best.expect("rescued unique has a best pair"));
                            (class, metrics, second_best_distance) = original;
                        } else {
                            (class, metrics, second_best_distance) = completed;
                            if let Some((original, completed)) = original_best.zip(completed_best)
                                && original.score() == completed.score()
                                && pair_origin_key(original, pair[0].len(), pair[1].len())
                                    != pair_origin_key(completed, pair[0].len(), pair[1].len())
                            {
                                self.pair.best_pairs.push(original);
                                collapse_equivalent_pair_origins(
                                    &mut self.pair.best_pairs,
                                    pair[0].len(),
                                    pair[1].len(),
                                );
                                prefer_minimum_net_gap_representative(
                                    &mut self.pair.best_pairs,
                                    pair[0].len(),
                                    pair[1].len(),
                                );
                                class = PairMappingStatus::Ambiguous;
                                metrics.best_pair_placements = metrics
                                    .best_pair_placements
                                    .max(
                                        u64::try_from(self.pair.best_pairs.len())
                                            .unwrap_or(u64::MAX),
                                    )
                                    .max(2);
                                metrics.near_best_pairings = metrics.near_best_pairings.max(1);
                            }
                        }
                    }
                    (Some(original_best), _) => {
                        self.pair.best_pairs.clear();
                        self.pair.best_pairs.push(original_best);
                        (class, metrics, second_best_distance) = original;
                    }
                    (None, _) => {}
                }
            }
            if search_mode.is_sensitive()
                && self
                    .pair
                    .should_affine_rescore(class, metrics, pair[0].len(), pair[1].len())
            {
                (class, metrics, second_best_distance) =
                    self.pair.affine_rescore_directional_pair(
                        reference,
                        pair[0],
                        pair[1],
                        class,
                        maximum_edit_distance,
                        minimum_template_span,
                        maximum_template_span,
                        metrics,
                    )?;
            }
            if search_mode.is_sensitive()
                && semi_global
                && sensitive_targeted_semi_global_required(
                    class,
                    metrics,
                    self.pair
                        .best_pairs()
                        .first()
                        .copied()
                        .map(PairedPlacement::distance),
                )
            {
                let original = (class, metrics, second_best_distance);
                let original_best = self.pair.best_pairs().first().copied();
                let original_confidence = sensitive_effective_mapping_quality(class, metrics);
                let mut candidate = self.pair.finish_directional_pair_combined(
                    reference,
                    pair[0],
                    pair[1],
                    maximum_edit_distance,
                    minimum_template_span,
                    maximum_template_span,
                    true,
                    metrics.mate1,
                    metrics.mate2,
                    metrics.window_rescue_attempted,
                )?;
                let candidate_best = self.pair.best_pairs().first().copied();
                let candidate_net_gap_profile =
                    pair_net_gap_profile(self.pair.best_pairs(), pair[0].len(), pair[1].len());
                let same_origin =
                    original_best
                        .zip(candidate_best)
                        .is_some_and(|(original, candidate)| {
                            pair_origin_key(original, pair[0].len(), pair[1].len())
                                == pair_origin_key(candidate, pair[0].len(), pair[1].len())
                        });
                ambiguity_q10_certified = sensitive_ambiguity_q10_certified(
                    original.0,
                    original.1,
                    original_best.map(PairedPlacement::score),
                    original_best.map(PairedPlacement::distance),
                    candidate.0,
                    candidate.1,
                    candidate_best.map(PairedPlacement::distance),
                    candidate_net_gap_profile,
                    same_origin,
                );
                let bounded_parsimony = sensitive_two_way_parsimony_q20_certified(
                    original.0,
                    original.1,
                    original_best,
                    candidate.0,
                    candidate.1,
                    candidate_best,
                    self.pair.best_pairs().len(),
                    candidate_net_gap_profile,
                    pair[0].len(),
                    pair[1].len(),
                    same_origin,
                );
                if bounded_parsimony {
                    // `prefer_minimum_net_gap_representative` has already put
                    // the uniquely parsimonious member first.  Keep the
                    // primary-score tie in the score metrics, but make the
                    // secondary decision explicit in class/cardinality.
                    self.pair.best_pairs.truncate(1);
                    candidate.0 = PairMappingStatus::Unique;
                    candidate.1.best_pair_placements = 1;
                }
                let candidate_confidence = if bounded_parsimony {
                    20
                } else {
                    sensitive_effective_mapping_quality(candidate.0, candidate.1)
                };
                repeat_risk_q20_certified = sensitive_stable_rescue_q20_certified(
                    original.0,
                    original.1,
                    original_confidence,
                    candidate.0,
                    candidate.1,
                    candidate_confidence,
                    same_origin,
                );
                // A first coordinate discovered only after endpoint clipping
                // has no independent full-read origin to stabilize it. Keep
                // it as unresolved; targeted semi-global is allowed to repair
                // confidence only at an origin already supported by the
                // strict/affine search.
                let completed_incomplete_q10 = !original.1.frontier_complete
                    && ambiguity_q10_certified
                    && candidate.1.frontier_complete
                    && same_origin;
                if completed_incomplete_q10 {
                    // Completion proved a stable representative, but this
                    // The qualified cell is only a Q10 ambiguity
                    // certificate.  Do not let the candidate's temporary
                    // Unique class enter the unique-result Q30/Q40 LUT.
                    candidate.0 = PairMappingStatus::Ambiguous;
                    candidate.1.best_pair_placements = candidate.1.best_pair_placements.max(2);
                }
                let accepted = same_origin
                    && ((original.1.frontier_complete
                        && (bounded_parsimony || candidate_confidence > original_confidence))
                        || completed_incomplete_q10);
                parsimony_q20_certified = bounded_parsimony && accepted;
                if accepted {
                    if bounded_parsimony {
                        // The targeted pass is secondary confidence evidence
                        // at the already supported biological origin.  Keep
                        // the strict representative for traceback/CIGAR: its
                        // full-read placement remains the primary alignment,
                        // while the uniquely parsimonious endpoint tie-break
                        // only certifies class and the Q20 boundary.
                        self.pair.best_pairs.clear();
                        self.pair
                            .best_pairs
                            .push(original_best.expect("certified parsimony has a strict origin"));
                    }
                    (class, metrics, second_best_distance) = candidate;
                } else {
                    self.pair.best_pairs.clear();
                    if let Some(original_best) = original_best {
                        self.pair.best_pairs.push(original_best);
                    }
                    (class, metrics, second_best_distance) = original;
                }
            }
            self.results.push(PairedBatchResult {
                class,
                metrics,
                best_pair: self.pair.best_pairs().first().copied(),
                second_best_distance,
                repeat_risk_q20_certified,
                parsimony_q20_certified,
                ambiguity_q10_certified,
                requires_positive_mapq_for_reporting,
            });
        }
        Ok(&self.results)
    }
}

fn swap_batch_result_mates(mut result: PairedBatchResult) -> PairedBatchResult {
    result.metrics = PairAlignmentMetrics {
        mate1: result.metrics.mate2,
        mate2: result.metrics.mate1,
        ..result.metrics
    };
    result.best_pair = result.best_pair.map(|pair| PairedPlacement {
        mate1: pair.mate2,
        mate2: pair.mate1,
        ..pair
    });
    result
}

// Classification, score confidence, and MAPQ evidence must be merged under
// the same directional tie decision.
#[allow(clippy::too_many_lines)]
fn merge_non_directional_batch_results(
    original: &PairedBatchResult,
    complementary: &PairedBatchResult,
) -> PairedBatchResult {
    let original_score = batch_best_score(original);
    let complementary_score = batch_best_score(complementary);
    let (mut selected, other, tied) = match (original_score, complementary_score) {
        (Some(left), Some(right)) if left > right => (*original, complementary, false),
        (Some(left), Some(right)) if right > left => (*complementary, original, false),
        (Some(_), Some(_)) => (*original, complementary, true),
        (None, Some(_)) => (*complementary, original, false),
        (_, None) => (*original, complementary, false),
    };
    let selected_score = batch_best_score(&selected);
    let other_score = batch_best_score(other);
    let selected_metrics = selected.metrics;
    let other_metrics = other.metrics;

    selected.metrics = PairAlignmentMetrics {
        mate1: merge_read_metrics(selected_metrics.mate1, other_metrics.mate1),
        mate2: merge_read_metrics(selected_metrics.mate2, other_metrics.mate2),
        compatible_pairs: selected_metrics
            .compatible_pairs
            .saturating_add(other_metrics.compatible_pairs),
        best_pair_placements: if tied {
            selected_metrics
                .best_pair_placements
                .saturating_add(other_metrics.best_pair_placements)
                .max(2)
        } else {
            selected_metrics.best_pair_placements
        },
        window_rescue_attempted: selected_metrics.window_rescue_attempted
            || other_metrics.window_rescue_attempted,
        semi_global_attempted: selected_metrics.semi_global_attempted
            || other_metrics.semi_global_attempted,
        best_pair_score: selected_score,
        second_best_pair_score: if tied {
            selected_score
        } else {
            maximum_optional_score(selected_metrics.second_best_pair_score, other_score)
        },
        near_best_pairings: selected_metrics
            .near_best_pairings
            .saturating_add(other_score.map_or(0, |other_score| {
                selected_score.map_or(0, |selected_score| {
                    if selected_score.saturating_sub(other_score) <= BWA_NEAR_SUBOPTIMAL_DELTA {
                        other_metrics.near_best_pairings.saturating_add(1)
                    } else {
                        0
                    }
                })
            })),
        mapq_compatible_pairs: selected_metrics
            .mapq_compatible_pairs
            .saturating_add(other_metrics.mapq_compatible_pairs),
        mapq_best_pair_score: maximum_optional_score(
            selected_metrics.mapq_best_pair_score,
            other_metrics.mapq_best_pair_score,
        ),
        mapq_second_best_pair_score: if tied {
            maximum_optional_score(
                selected_metrics.mapq_best_pair_score,
                other_metrics.mapq_best_pair_score,
            )
        } else {
            maximum_optional_score(
                selected_metrics.mapq_second_best_pair_score,
                other_metrics.mapq_best_pair_score,
            )
        },
        mapq_near_best_pairings: selected_metrics.mapq_near_best_pairings.saturating_add(
            other_metrics.mapq_best_pair_score.map_or(0, |other_score| {
                selected_metrics
                    .mapq_best_pair_score
                    .map_or(0, |selected_score| {
                        if selected_score.saturating_sub(other_score) <= BWA_NEAR_SUBOPTIMAL_DELTA {
                            other_metrics.mapq_near_best_pairings.saturating_add(1)
                        } else {
                            0
                        }
                    })
            }),
        ),
        frontier_complete: selected_metrics.frontier_complete && other_metrics.frontier_complete,
    };
    selected.second_best_distance = if tied {
        selected.best_pair.map(PairedPlacement::score)
    } else {
        minimum_optional_distance(
            selected.second_best_distance,
            other.best_pair.map(PairedPlacement::score),
        )
    };
    if tied {
        selected.class = PairMappingStatus::Ambiguous;
        selected.repeat_risk_q20_certified = false;
        selected.parsimony_q20_certified = false;
    } else if !selected.metrics.frontier_complete
        && matches!(selected.class, PairMappingStatus::Unique)
    {
        selected.class = PairMappingStatus::Ambiguous;
        selected.metrics.best_pair_placements = selected.metrics.best_pair_placements.max(2);
        selected.repeat_risk_q20_certified = false;
        selected.parsimony_q20_certified = false;
    } else if other_score.is_some() {
        // Cross-conversion evidence did not participate in the directional
        // pass's specialized confidence repair, so do not carry that repair
        // certificate across the global four-strand decision.
        selected.repeat_risk_q20_certified = false;
        selected.parsimony_q20_certified = false;
    }
    selected
}

fn batch_best_score(result: &PairedBatchResult) -> Option<i16> {
    result.metrics.best_pair_score.or_else(|| {
        result
            .best_pair
            .map(|pair| -i16::from(pair.score()).saturating_mul(BWA_MISMATCH_PENALTY))
    })
}

const fn maximum_optional_score(left: Option<i16>, right: Option<i16>) -> Option<i16> {
    match (left, right) {
        (Some(left), Some(right)) => Some(if left > right { left } else { right }),
        (Some(value), None) | (None, Some(value)) => Some(value),
        (None, None) => None,
    }
}

const fn minimum_optional_distance(left: Option<u8>, right: Option<u8>) -> Option<u8> {
    match (left, right) {
        (Some(left), Some(right)) => Some(if left < right { left } else { right }),
        (Some(value), None) | (None, Some(value)) => Some(value),
        (None, None) => None,
    }
}

const fn merge_read_metrics(
    left: ReadAlignmentMetrics,
    right: ReadAlignmentMetrics,
) -> ReadAlignmentMetrics {
    ReadAlignmentMetrics {
        located_rows: left.located_rows.saturating_add(right.located_rows),
        emitted_candidate_starts: left
            .emitted_candidate_starts
            .saturating_add(right.emitted_candidate_starts),
        distinct_candidate_starts: left
            .distinct_candidate_starts
            .saturating_add(right.distinct_candidate_starts),
        verified_placements: left
            .verified_placements
            .saturating_add(right.verified_placements),
    }
}

fn sensitive_repeat_recheck_required(
    class: PairMappingStatus,
    metrics: PairAlignmentMetrics,
) -> bool {
    matches!(class, PairMappingStatus::Unique)
        && (metrics.window_rescue_attempted
            || metrics.mate1.located_rows.max(metrics.mate2.located_rows)
                >= SENSITIVE_REPEAT_RECHECK_ROWS)
}

fn sensitive_targeted_semi_global_required(
    class: PairMappingStatus,
    metrics: PairAlignmentMetrics,
    pair_distance: Option<u8>,
) -> bool {
    !matches!(class, PairMappingStatus::Unmapped)
        && ((metrics.frontier_complete && sensitive_effective_mapping_quality(class, metrics) < 20)
            || sensitive_incomplete_sparse_completion_required(class, metrics, pair_distance))
}

fn conservatively_mark_incomplete_frontier(
    result: &mut (PairMappingStatus, PairAlignmentMetrics, Option<u8>),
    complete: bool,
) {
    result.1.frontier_complete = complete;
    if !complete && matches!(result.0, PairMappingStatus::Unique) {
        result.0 = PairMappingStatus::Ambiguous;
        result.1.best_pair_placements = result.1.best_pair_placements.max(2);
        result.2 = None;
    }
}

impl PairWorkspace {
    #[must_use]
    fn with_capacity(
        mate_candidate_capacity: usize,
        mate_placement_capacity: usize,
        pair_capacity: usize,
    ) -> Self {
        Self {
            mate1: ReadWorkspace::with_capacity(mate_candidate_capacity, mate_placement_capacity),
            mate2: ReadWorkspace::with_capacity(mate_candidate_capacity, mate_placement_capacity),
            rescue_windows: Vec::with_capacity(mate_placement_capacity),
            best_pairs: Vec::with_capacity(pair_capacity),
            exact_anchor_candidates: Vec::with_capacity(
                usize::try_from(SEMI_GLOBAL_MAX_EXACT_ANCHOR_HITS)
                    .expect("exact anchor limit fits usize"),
            ),
            ranked_anchor_placements: Vec::with_capacity(mate_placement_capacity),
            mate1_affine_scores: Vec::with_capacity(mate_placement_capacity),
            mate2_affine_scores: Vec::with_capacity(mate_placement_capacity),
            affine: AffineScoreWorkspace::default(),
            semi_global_clip_penalty: SEMI_GLOBAL_CLIP_PENALTY,
            prefer_minimum_net_gap: false,
            origin_pair_evidence: std::collections::HashMap::with_capacity(pair_capacity),
            combined_search_state: CombinedTwoLaneSearchState::new(),
            fallback_mate1_nominals: Vec::with_capacity(mate_candidate_capacity / 5 + 1),
            fallback_mate2_nominals: Vec::with_capacity(mate_candidate_capacity / 5 + 1),
            ranked_extension_selections: [None; 2],
        }
    }

    fn selective_unmapped_frontier_deepening_required(
        &self,
        ordinary_frontier_metrics: Option<PairAlignmentMetrics>,
    ) -> bool {
        selective_unmapped_frontier_deepening_required(
            self.ranked_extension_selections,
            ordinary_frontier_metrics,
        )
    }

    fn selective_unmapped_frontier_requires_positive_mapq(&self) -> bool {
        let [Some(first), Some(second)] = self.ranked_extension_selections else {
            return false;
        };
        !(SENSITIVE_POSITIVE_MAPQ_REPORTING_MIN_RETAINED_HITS
            ..SENSITIVE_POSITIVE_MAPQ_REPORTING_MAX_RETAINED_HITS)
            .contains(&first.retained_hits.saturating_add(second.retained_hits))
    }

    #[allow(clippy::too_many_arguments)]
    fn collect_ranked_block_seeds_for_lane(
        reference: &ReferenceIndex,
        read: &[Base],
        reversed_projected: &[ProjectedBase],
        lane: usize,
        maximum_edit_distance: u8,
        maximum_ranked_block_hits: u64,
        output: &mut [Option<RankedBlockSeed>; SENSITIVE_PROOF_BLOCKS],
    ) -> Result<Option<RankedBlockSelection>, AlignmentError> {
        let budget = usize::from(maximum_edit_distance);
        debug_assert!(lane < 2);
        debug_assert!(budget < SENSITIVE_PROOF_BLOCKS);
        let selection = collect_ranked_block_seeds(
            reference,
            read,
            reversed_projected,
            lane == 0,
            maximum_edit_distance,
            maximum_ranked_block_hits,
            output,
        )?;
        Ok(selection)
    }

    fn append_ranked_block_candidates_for_lane(
        &mut self,
        reference: &ReferenceIndex,
        read_len: usize,
        lane: usize,
        maximum_edit_distance: u8,
        seeds: &[Option<RankedBlockSeed>; SENSITIVE_PROOF_BLOCKS],
    ) -> Result<u64, AlignmentError> {
        let budget = usize::from(maximum_edit_distance);
        debug_assert!(lane < 2);
        debug_assert!(budget < SENSITIVE_PROOF_BLOCKS);
        let candidates = if lane == 0 {
            &mut self.mate1.candidate_nominals
        } else {
            &mut self.mate2.candidate_nominals
        };
        let located_rows =
            append_ranked_block_candidates(reference, read_len, lane == 0, seeds, candidates)?;
        Ok(located_rows)
    }

    // This is the internal handoff of one fully prepared pair; a parameter
    // object would obscure which search proof each value belongs to.
    #[allow(clippy::too_many_arguments)]
    fn map_directional_pair_combined_prepared(
        &mut self,
        reference: &ReferenceIndex,
        read1: &[Base],
        read2: &[Base],
        projected: [&[ProjectedBase]; 2],
        first_seeds: [Option<CombinedSeedMatches>; 2],
        maximum_edit_distance: u8,
        minimum_template_span: u64,
        maximum_template_span: u64,
        window_rescue: bool,
        semi_global: bool,
        search_limits: CombinedSearchLimits,
        preserve_fallback_frontier: bool,
    ) -> Result<(PairMappingStatus, PairAlignmentMetrics, Option<u8>), AlignmentError> {
        {
            self.mate1.begin_verification_cache_read();
            self.mate2.begin_verification_cache_read();
        }
        self.best_pairs.clear();
        self.mate1.candidates.clear();
        self.mate1.candidate_nominals.clear();
        self.mate1.placements.clear();
        self.mate2.candidates.clear();
        self.mate2.candidate_nominals.clear();
        self.mate2.placements.clear();
        self.combined_search_state = CombinedTwoLaneSearchState::new();
        self.fallback_mate1_nominals.clear();
        self.fallback_mate2_nominals.clear();
        if !semi_global
            && (read1.iter().filter(|base| base.is_unknown()).count()
                > usize::from(maximum_edit_distance)
                || read2.iter().filter(|base| base.is_unknown()).count()
                    > usize::from(maximum_edit_distance))
        {
            return Ok((PairMappingStatus::Unmapped, empty_pair_metrics(), None));
        }
        self.combined_search_state = start_combined_two_lane_search(
            reference,
            [read1, read2],
            projected,
            first_seeds,
            [false, true],
            search_limits,
            &mut self.mate1.candidate_nominals,
            &mut self.mate2.candidate_nominals,
        )?;
        let located_rows = self.combined_search_state.located;
        let mate1_metrics = ReadAlignmentMetrics {
            located_rows: located_rows[0],
            ..ReadAlignmentMetrics::default()
        };
        let mate2_metrics = ReadAlignmentMetrics {
            located_rows: located_rows[1],
            ..ReadAlignmentMetrics::default()
        };
        self.verify_directional_pair_combined_frontier(
            reference,
            read1,
            read2,
            projected,
            first_seeds,
            maximum_edit_distance,
            minimum_template_span,
            maximum_template_span,
            window_rescue,
            semi_global,
            search_limits,
            preserve_fallback_frontier,
            mate1_metrics,
            mate2_metrics,
        )
    }

    // Geometry filtering, rescue selection, and verification share one
    // candidate frontier and must update its metrics atomically.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn verify_directional_pair_combined_frontier(
        &mut self,
        reference: &ReferenceIndex,
        read1: &[Base],
        read2: &[Base],
        projected: [&[ProjectedBase]; 2],
        first_seeds: [Option<CombinedSeedMatches>; 2],
        maximum_edit_distance: u8,
        minimum_template_span: u64,
        maximum_template_span: u64,
        window_rescue: bool,
        semi_global: bool,
        search_limits: CombinedSearchLimits,
        preserve_fallback_frontier: bool,
        mate1_metrics: ReadAlignmentMetrics,
        mate2_metrics: ReadAlignmentMetrics,
    ) -> Result<(PairMappingStatus, PairAlignmentMetrics, Option<u8>), AlignmentError> {
        sort_nominal_candidates(&mut self.mate1.candidate_nominals);
        sort_nominal_candidates(&mut self.mate2.candidate_nominals);
        let nominal_geometry = nominal_pair_geometry_exists(
            &self.mate1.candidate_nominals,
            &self.mate2.candidate_nominals,
            read1.len(),
            read2.len(),
            maximum_template_span,
            maximum_edit_distance,
        );
        if preserve_fallback_frontier && !nominal_geometry {
            self.fallback_mate1_nominals
                .clone_from(&self.mate1.candidate_nominals);
            self.fallback_mate2_nominals
                .clone_from(&self.mate2.candidate_nominals);
        }
        let rescue_anchor = if window_rescue && !nominal_geometry {
            select_combined_window_rescue_anchor(
                first_seeds,
                &self.mate1.candidate_nominals,
                &self.mate2.candidate_nominals,
            )
        } else {
            None
        };
        let (mate1_metrics, mate2_metrics, window_rescue_attempted) = match rescue_anchor {
            Some(0) => {
                let (anchors, mate1_metrics) = self.mate1.verify_sorted_candidates_with_budget(
                    reference,
                    read1,
                    mate1_metrics,
                    maximum_edit_distance,
                )?;
                if anchors.is_empty() {
                    self.mate2.candidate_nominals.clear();
                    self.mate2.candidates.clear();
                    self.mate2.placements.clear();
                    (mate1_metrics, mate2_metrics, false)
                } else {
                    let mut rescued_metrics = rescue_from_combined_exact_blocks(
                        &mut self.mate2,
                        &mut self.rescue_windows,
                        reference,
                        read2,
                        projected[1],
                        anchors,
                        false,
                        maximum_template_span,
                        maximum_edit_distance,
                        preserve_fallback_frontier,
                        search_limits.maximum_combined_rescue_hits,
                    )?;
                    rescued_metrics.located_rows = mate2_metrics.located_rows;
                    (mate1_metrics, rescued_metrics, true)
                }
            }
            Some(1) => {
                let (anchors, mate2_metrics) = self.mate2.verify_sorted_candidates_with_budget(
                    reference,
                    read2,
                    mate2_metrics,
                    maximum_edit_distance,
                )?;
                if anchors.is_empty() {
                    self.mate1.candidate_nominals.clear();
                    self.mate1.candidates.clear();
                    self.mate1.placements.clear();
                    (mate1_metrics, mate2_metrics, false)
                } else {
                    let mut rescued_metrics = rescue_from_combined_exact_blocks(
                        &mut self.mate1,
                        &mut self.rescue_windows,
                        reference,
                        read1,
                        projected[0],
                        anchors,
                        true,
                        maximum_template_span,
                        maximum_edit_distance,
                        preserve_fallback_frontier,
                        search_limits.maximum_combined_rescue_hits,
                    )?;
                    rescued_metrics.located_rows = mate1_metrics.located_rows;
                    (rescued_metrics, mate2_metrics, true)
                }
            }
            Some(_) => unreachable!("a pair has exactly two mates"),
            None => {
                retain_nominal_pair_geometry(
                    &mut self.mate1.candidate_nominals,
                    &mut self.mate2.candidate_nominals,
                    read1.len(),
                    read2.len(),
                    maximum_template_span,
                    maximum_edit_distance,
                );
                let (_, mate1_metrics) = self.mate1.verify_sorted_candidates_with_budget(
                    reference,
                    read1,
                    mate1_metrics,
                    maximum_edit_distance,
                )?;
                let (_, mate2_metrics) = self.mate2.verify_sorted_candidates_with_budget(
                    reference,
                    read2,
                    mate2_metrics,
                    maximum_edit_distance,
                )?;
                (mate1_metrics, mate2_metrics, false)
            }
        };
        self.finish_directional_pair_combined(
            reference,
            read1,
            read2,
            maximum_edit_distance,
            minimum_template_span,
            maximum_template_span,
            semi_global,
            mate1_metrics,
            mate2_metrics,
            window_rescue_attempted,
        )
    }

    /// Reuses the candidate frontier left by the initial d3 pass and only reruns
    /// bounded verification at the requested edit budget. A pair that remains
    /// unmapped can still fall through to the deeper incremental seed search.
    #[allow(clippy::too_many_arguments)]
    fn reverify_directional_pair_combined_candidates(
        &mut self,
        reference: &ReferenceIndex,
        read1: &[Base],
        read2: &[Base],
        maximum_edit_distance: u8,
        minimum_template_span: u64,
        maximum_template_span: u64,
        semi_global: bool,
        previous_metrics: PairAlignmentMetrics,
    ) -> Result<(PairMappingStatus, PairAlignmentMetrics, Option<u8>), AlignmentError> {
        self.best_pairs.clear();
        self.mate1.candidates.clear();
        self.mate1.placements.clear();
        self.mate2.candidates.clear();
        self.mate2.placements.clear();
        self.mate1
            .candidate_nominals
            .extend_from_slice(&self.fallback_mate1_nominals);
        self.mate2
            .candidate_nominals
            .extend_from_slice(&self.fallback_mate2_nominals);
        sort_nominal_candidates(&mut self.mate1.candidate_nominals);
        sort_nominal_candidates(&mut self.mate2.candidate_nominals);
        retain_nominal_pair_geometry(
            &mut self.mate1.candidate_nominals,
            &mut self.mate2.candidate_nominals,
            read1.len(),
            read2.len(),
            maximum_template_span,
            maximum_edit_distance,
        );
        let mate1_metrics = ReadAlignmentMetrics {
            located_rows: previous_metrics.mate1.located_rows,
            ..ReadAlignmentMetrics::default()
        };
        let mate2_metrics = ReadAlignmentMetrics {
            located_rows: previous_metrics.mate2.located_rows,
            ..ReadAlignmentMetrics::default()
        };
        let (_, mate1_metrics) = self.mate1.verify_sorted_candidates_with_budget(
            reference,
            read1,
            mate1_metrics,
            maximum_edit_distance,
        )?;
        let (_, mate2_metrics) = self.mate2.verify_sorted_candidates_with_budget(
            reference,
            read2,
            mate2_metrics,
            maximum_edit_distance,
        )?;
        self.finish_directional_pair_combined(
            reference,
            read1,
            read2,
            maximum_edit_distance,
            minimum_template_span,
            maximum_template_span,
            semi_global,
            mate1_metrics,
            mate2_metrics,
            previous_metrics.window_rescue_attempted,
        )
    }

    /// Replays only initial-pass seeds that become admissible in the incremental
    /// fallback,
    /// then continues from the saved seed offset into the additional round.
    /// This preserves the deeper frontier without repeating rounds 0 through 4.
    #[allow(clippy::too_many_arguments)]
    fn continue_directional_pair_combined_incremental(
        &mut self,
        reference: &ReferenceIndex,
        read1: &[Base],
        read2: &[Base],
        projected: [&[ProjectedBase]; 2],
        first_seeds: [Option<CombinedSeedMatches>; 2],
        maximum_edit_distance: u8,
        minimum_template_span: u64,
        maximum_template_span: u64,
        window_rescue: bool,
        semi_global: bool,
        previous_metrics: PairAlignmentMetrics,
    ) -> Result<(PairMappingStatus, PairAlignmentMetrics, Option<u8>), AlignmentError> {
        if !self.combined_search_state.initialized {
            return self.map_directional_pair_combined_prepared(
                reference,
                read1,
                read2,
                projected,
                first_seeds,
                maximum_edit_distance,
                minimum_template_span,
                maximum_template_span,
                window_rescue,
                semi_global,
                PairedSearchMode::Default.limits(),
                false,
            );
        }
        self.best_pairs.clear();
        self.mate1.candidates.clear();
        self.mate1.placements.clear();
        self.mate2.candidates.clear();
        self.mate2.placements.clear();
        self.mate1
            .candidate_nominals
            .append(&mut self.fallback_mate1_nominals);
        self.mate2
            .candidate_nominals
            .append(&mut self.fallback_mate2_nominals);
        let additional_located = continue_combined_two_lane_search(
            reference,
            [read1, read2],
            projected,
            [false, true],
            &mut self.combined_search_state,
            &mut self.mate1.candidate_nominals,
            &mut self.mate2.candidate_nominals,
        )?;
        let mate1_metrics = ReadAlignmentMetrics {
            located_rows: previous_metrics
                .mate1
                .located_rows
                .saturating_add(additional_located[0]),
            ..ReadAlignmentMetrics::default()
        };
        let mate2_metrics = ReadAlignmentMetrics {
            located_rows: previous_metrics
                .mate2
                .located_rows
                .saturating_add(additional_located[1]),
            ..ReadAlignmentMetrics::default()
        };
        self.verify_directional_pair_combined_frontier(
            reference,
            read1,
            read2,
            projected,
            first_seeds,
            maximum_edit_distance,
            minimum_template_span,
            maximum_template_span,
            window_rescue,
            semi_global,
            PairedSearchMode::Default.limits(),
            false,
            mate1_metrics,
            mate2_metrics,
        )
    }

    /// Sensitive failed-pair completion.
    ///
    /// The `d + 1` disjoint exact blocks are ranked by occurrence count and the
    /// rarest intervals fitting the global enumeration budget form the anchor
    /// frontier. Every verified anchor then induces a bounded window in which
    /// the partner is completed with the full disjoint-block proof. A result
    /// recovered from an incomplete global block frontier is reported only as
    /// ambiguous; no truth coordinate participates in paired-end
    /// classification.
    // The ranked proof pipeline intentionally retains its seed, verification,
    // and fallback frontiers in one worker-owned transaction.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn extend_directional_pair_from_ranked_blocks(
        &mut self,
        reference: &ReferenceIndex,
        read1: &[Base],
        read2: &[Base],
        projected: [&[ProjectedBase]; 2],
        maximum_edit_distance: u8,
        minimum_template_span: u64,
        maximum_template_span: u64,
        semi_global: bool,
        maximum_ranked_block_hits: u64,
    ) -> Result<Option<(PairMappingStatus, PairAlignmentMetrics, Option<u8>)>, AlignmentError> {
        let mut seed_sets = [[None; SENSITIVE_PROOF_BLOCKS]; 2];
        let first_selection = Self::collect_ranked_block_seeds_for_lane(
            reference,
            read1,
            projected[0],
            0,
            maximum_edit_distance,
            maximum_ranked_block_hits,
            &mut seed_sets[0],
        )?;
        let second_selection = Self::collect_ranked_block_seeds_for_lane(
            reference,
            read2,
            projected[1],
            1,
            maximum_edit_distance,
            maximum_ranked_block_hits,
            &mut seed_sets[1],
        )?;
        {
            self.ranked_extension_selections = [first_selection, second_selection];
        }
        // The paired geometry is much more selective than either short
        // bisulfite block by itself.  Intersect both bounded frontiers before
        // invoking the d5 verifier; the anchor/window path below remains the
        // fallback when only one mate has an informative retained block.
        if let (Some(first_ranked), Some(second_ranked)) = (first_selection, second_selection) {
            self.best_pairs.clear();
            self.mate1.candidates.clear();
            self.mate1.candidate_nominals.clear();
            self.mate1.placements.clear();
            self.mate2.candidates.clear();
            self.mate2.candidate_nominals.clear();
            self.mate2.placements.clear();
            let mate1_rows = self.append_ranked_block_candidates_for_lane(
                reference,
                read1.len(),
                0,
                maximum_edit_distance,
                &seed_sets[0],
            )?;
            let mate2_rows = self.append_ranked_block_candidates_for_lane(
                reference,
                read2.len(),
                1,
                maximum_edit_distance,
                &seed_sets[1],
            )?;
            sort_nominal_candidates(&mut self.mate1.candidate_nominals);
            sort_nominal_candidates(&mut self.mate2.candidate_nominals);
            retain_nominal_pair_geometry(
                &mut self.mate1.candidate_nominals,
                &mut self.mate2.candidate_nominals,
                read1.len(),
                read2.len(),
                maximum_template_span,
                maximum_edit_distance,
            );
            if !self.mate1.candidate_nominals.is_empty()
                && !self.mate2.candidate_nominals.is_empty()
            {
                let mate1_seed_metrics = ReadAlignmentMetrics {
                    located_rows: mate1_rows,
                    ..ReadAlignmentMetrics::default()
                };
                let mate2_seed_metrics = ReadAlignmentMetrics {
                    located_rows: mate2_rows,
                    ..ReadAlignmentMetrics::default()
                };
                let (_, mate1_metrics) = self.mate1.verify_sorted_candidates_with_budget(
                    reference,
                    read1,
                    mate1_seed_metrics,
                    maximum_edit_distance,
                )?;
                let (_, mate2_metrics) = self.mate2.verify_sorted_candidates_with_budget(
                    reference,
                    read2,
                    mate2_seed_metrics,
                    maximum_edit_distance,
                )?;
                let mut completed = self.finish_directional_pair_combined(
                    reference,
                    read1,
                    read2,
                    maximum_edit_distance,
                    minimum_template_span,
                    maximum_template_span,
                    semi_global,
                    mate1_metrics,
                    mate2_metrics,
                    false,
                )?;
                if !matches!(completed.0, PairMappingStatus::Unmapped) {
                    let complete_frontier = first_ranked.complete && second_ranked.complete;
                    if !complete_frontier
                        && let Some(certified) = self.certify_ranked_pair_frontier(
                            reference,
                            read1,
                            read2,
                            projected,
                            maximum_edit_distance,
                            minimum_template_span,
                            maximum_template_span,
                            semi_global,
                            SENSITIVE_RANKED_BLOCK_HITS.saturating_mul(2),
                            completed,
                        )?
                    {
                        return Ok(Some(certified));
                    }
                    conservatively_mark_incomplete_frontier(&mut completed, complete_frontier);
                    return Ok(Some(completed));
                }
            }
        }

        let anchor_lane = match (first_selection, second_selection) {
            (None, None) => return Ok(None),
            (Some(_), None) => 0,
            (None, Some(_)) => 1,
            (Some(first), Some(second)) => match (first.complete, second.complete) {
                (true, false) => 0,
                (false, true) => 1,
                _ => usize::from(second.retained_hits < first.retained_hits),
            },
        };
        let anchor_frontier_complete = [first_selection, second_selection][anchor_lane]
            .expect("selected anchor has retained hits")
            .complete;

        self.best_pairs.clear();
        self.mate1.candidates.clear();
        self.mate1.candidate_nominals.clear();
        self.mate1.placements.clear();
        self.mate2.candidates.clear();
        self.mate2.candidate_nominals.clear();
        self.mate2.placements.clear();
        self.ranked_anchor_placements.clear();

        let located_rows = if anchor_lane == 0 {
            self.append_ranked_block_candidates_for_lane(
                reference,
                read1.len(),
                0,
                maximum_edit_distance,
                &seed_sets[0],
            )?
        } else {
            self.append_ranked_block_candidates_for_lane(
                reference,
                read2.len(),
                1,
                maximum_edit_distance,
                &seed_sets[1],
            )?
        };
        let anchor_metrics = ReadAlignmentMetrics {
            located_rows,
            ..ReadAlignmentMetrics::default()
        };
        let anchor_metrics = if anchor_lane == 0 {
            let (_, metrics) = self.mate1.verify_candidates_with_budget(
                reference,
                read1,
                anchor_metrics,
                maximum_edit_distance,
            )?;
            self.ranked_anchor_placements
                .extend_from_slice(&self.mate1.placements);
            metrics
        } else {
            let (_, metrics) = self.mate2.verify_candidates_with_budget(
                reference,
                read2,
                anchor_metrics,
                maximum_edit_distance,
            )?;
            self.ranked_anchor_placements
                .extend_from_slice(&self.mate2.placements);
            metrics
        };

        let (mate1_metrics, mate2_metrics) = if self.ranked_anchor_placements.is_empty() {
            if anchor_lane == 0 {
                (anchor_metrics, ReadAlignmentMetrics::default())
            } else {
                (ReadAlignmentMetrics::default(), anchor_metrics)
            }
        } else if anchor_lane == 0 {
            let partner_metrics = rescue_from_ranked_anchor_windows(
                &mut self.mate2,
                &mut self.rescue_windows,
                reference,
                read2,
                &self.ranked_anchor_placements,
                false,
                maximum_template_span,
                maximum_edit_distance,
            )?;
            (anchor_metrics, partner_metrics)
        } else {
            let partner_metrics = rescue_from_ranked_anchor_windows(
                &mut self.mate1,
                &mut self.rescue_windows,
                reference,
                read1,
                &self.ranked_anchor_placements,
                true,
                maximum_template_span,
                maximum_edit_distance,
            )?;
            (partner_metrics, anchor_metrics)
        };
        let mut completed = self.finish_directional_pair_combined(
            reference,
            read1,
            read2,
            maximum_edit_distance,
            minimum_template_span,
            maximum_template_span,
            semi_global,
            mate1_metrics,
            mate2_metrics,
            !self.ranked_anchor_placements.is_empty(),
        )?;
        if !anchor_frontier_complete
            && let Some(certified) = self.certify_ranked_pair_frontier(
                reference,
                read1,
                read2,
                projected,
                maximum_edit_distance,
                minimum_template_span,
                maximum_template_span,
                semi_global,
                SENSITIVE_RANKED_BLOCK_HITS.saturating_mul(2),
                completed,
            )?
        {
            return Ok(Some(certified));
        }
        conservatively_mark_incomplete_frontier(&mut completed, anchor_frontier_complete);
        Ok(Some(completed))
    }

    /// Attempts a complete uniqueness proof at the score already discovered by
    /// an incomplete maximum-distance frontier.
    ///
    /// Under strict scoring, a pair with score `k` can disrupt at most `k`
    /// query blocks. Under semi-global scoring every retained mismatch costs
    /// seven and every clipped base has the configured sensitive
    /// penalty, so dividing by the smaller event penalty is a conservative
    /// bound. A `k + 1` partition therefore covers every pair
    /// that could tie or beat the selected score while using longer, much less
    /// repetitive exact blocks than the original distance-five search.
    // Certification replays the complete ranked frontier while preserving the
    // originally selected pair for conservative fallback.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn certify_ranked_pair_frontier(
        &mut self,
        reference: &ReferenceIndex,
        read1: &[Base],
        read2: &[Base],
        projected: [&[ProjectedBase]; 2],
        maximum_edit_distance: u8,
        minimum_template_span: u64,
        maximum_template_span: u64,
        semi_global: bool,
        maximum_combined_proof_hits: u64,
        original: (PairMappingStatus, PairAlignmentMetrics, Option<u8>),
    ) -> Result<Option<(PairMappingStatus, PairAlignmentMetrics, Option<u8>)>, AlignmentError> {
        if !matches!(original.0, PairMappingStatus::Unique) {
            return Ok(None);
        }
        let Some(original_best) = self.best_pairs.first().copied() else {
            return Ok(None);
        };
        let proof_budget = if original.1.semi_global_attempted {
            original_best.score() / SENSITIVE_MIN_EVENT_PENALTY
        } else {
            original_best.score()
        };
        if proof_budget >= maximum_edit_distance {
            return Ok(None);
        }

        let mut seed_sets = [[None; SENSITIVE_PROOF_BLOCKS]; 2];
        let first_selection = Self::collect_ranked_block_seeds_for_lane(
            reference,
            read1,
            projected[0],
            0,
            proof_budget,
            SENSITIVE_RANKED_BLOCK_HITS,
            &mut seed_sets[0],
        )?;
        let second_selection = Self::collect_ranked_block_seeds_for_lane(
            reference,
            read2,
            projected[1],
            1,
            proof_budget,
            SENSITIVE_RANKED_BLOCK_HITS,
            &mut seed_sets[1],
        )?;
        let Some((first_selection, second_selection)) = first_selection.zip(second_selection)
        else {
            return Ok(None);
        };
        if !first_selection.complete || !second_selection.complete {
            return Ok(None);
        }
        if first_selection
            .retained_hits
            .saturating_add(second_selection.retained_hits)
            > maximum_combined_proof_hits
        {
            return Ok(None);
        }

        // The proof partition is complete for every per-mate placement whose
        // distance can contribute to a pair tying or beating `original_best`.
        // Verification still uses the paired-end maximum-distance budget so
        // candidates outside the proof partition retain ordinary semantics.
        let certification_budget = maximum_edit_distance;

        self.best_pairs.clear();
        self.mate1.candidates.clear();
        self.mate1.candidate_nominals.clear();
        self.mate1.placements.clear();
        self.mate2.candidates.clear();
        self.mate2.candidate_nominals.clear();
        self.mate2.placements.clear();
        let mate1_rows = self.append_ranked_block_candidates_for_lane(
            reference,
            read1.len(),
            0,
            proof_budget,
            &seed_sets[0],
        )?;
        let mate2_rows = self.append_ranked_block_candidates_for_lane(
            reference,
            read2.len(),
            1,
            proof_budget,
            &seed_sets[1],
        )?;
        sort_nominal_candidates(&mut self.mate1.candidate_nominals);
        sort_nominal_candidates(&mut self.mate2.candidate_nominals);
        retain_nominal_pair_geometry(
            &mut self.mate1.candidate_nominals,
            &mut self.mate2.candidate_nominals,
            read1.len(),
            read2.len(),
            maximum_template_span,
            certification_budget,
        );
        if self.mate1.candidate_nominals.is_empty() || self.mate2.candidate_nominals.is_empty() {
            self.best_pairs.push(original_best);
            return Ok(None);
        }

        let mate1_seed_metrics = ReadAlignmentMetrics {
            located_rows: mate1_rows,
            ..ReadAlignmentMetrics::default()
        };
        let mate2_seed_metrics = ReadAlignmentMetrics {
            located_rows: mate2_rows,
            ..ReadAlignmentMetrics::default()
        };
        let (_, mate1_metrics) = self.mate1.verify_sorted_candidates_with_budget(
            reference,
            read1,
            mate1_seed_metrics,
            certification_budget,
        )?;
        let (_, mate2_metrics) = self.mate2.verify_sorted_candidates_with_budget(
            reference,
            read2,
            mate2_seed_metrics,
            certification_budget,
        )?;
        let certified = self.finish_directional_pair_combined(
            reference,
            read1,
            read2,
            certification_budget,
            minimum_template_span,
            maximum_template_span,
            semi_global,
            mate1_metrics,
            mate2_metrics,
            false,
        )?;
        if matches!(certified.0, PairMappingStatus::Unmapped)
            || self
                .best_pairs
                .first()
                .is_none_or(|best| best.score() > original_best.score())
        {
            self.best_pairs.clear();
            self.best_pairs.push(original_best);
            return Ok(None);
        }
        Ok(Some(certified))
    }

    // Pair selection, optional rescoring, and confidence aggregation consume
    // the same workspace state and therefore remain one finishing pass.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn finish_directional_pair_combined(
        &mut self,
        reference: &ReferenceIndex,
        read1: &[Base],
        read2: &[Base],
        maximum_edit_distance: u8,
        minimum_template_span: u64,
        maximum_template_span: u64,
        semi_global: bool,
        mut mate1_metrics: ReadAlignmentMetrics,
        mut mate2_metrics: ReadAlignmentMetrics,
        window_rescue_attempted: bool,
    ) -> Result<(PairMappingStatus, PairAlignmentMetrics, Option<u8>), AlignmentError> {
        let mut selection = if self.prefer_minimum_net_gap {
            select_best_pair_origins_with_endpoint_policy(
                &self.mate1.placements,
                &self.mate2.placements,
                [read1, read2],
                maximum_edit_distance,
                minimum_template_span,
                maximum_template_span,
                false,
                &mut self.origin_pair_evidence,
                &mut self.best_pairs,
            )
        } else {
            select_best_pairs(
                &self.mate1.placements,
                &self.mate2.placements,
                maximum_edit_distance,
                minimum_template_span,
                maximum_template_span,
                &mut self.best_pairs,
            )
        };
        collapse_equivalent_pair_origins(&mut self.best_pairs, read1.len(), read2.len());
        if self.prefer_minimum_net_gap && self.best_pairs.len() > 1 {
            prefer_minimum_net_gap_representative(&mut self.best_pairs, read1.len(), read2.len());
        }
        let semi_global_attempted = semi_global
            && (self.best_pairs.len() != 1 || self.best_pairs[0].distance() != 0)
            && !self.mate1.candidate_nominals.is_empty()
            && !self.mate2.candidate_nominals.is_empty();
        if semi_global_attempted {
            append_ungapped_semi_global_placements(
                &mut self.mate1,
                reference,
                read1,
                maximum_edit_distance,
                self.semi_global_clip_penalty,
            );
            mate1_metrics.verified_placements =
                u64::try_from(self.mate1.placements.len()).unwrap_or(u64::MAX);
            append_ungapped_semi_global_placements(
                &mut self.mate2,
                reference,
                read2,
                maximum_edit_distance,
                self.semi_global_clip_penalty,
            );
            mate2_metrics.verified_placements =
                u64::try_from(self.mate2.placements.len()).unwrap_or(u64::MAX);
            // The pair join uses partition points over spatially sorted mate-2
            // placements. Appending a second, individually ordered frontier
            // does not preserve that global order.
            self.mate1
                .placements
                .sort_unstable_by_key(|placement| spatial_key(*placement));
            self.mate2
                .placements
                .sort_unstable_by_key(|placement| spatial_key(*placement));
            selection = if self.prefer_minimum_net_gap {
                select_best_pair_origins_with_endpoint_policy(
                    &self.mate1.placements,
                    &self.mate2.placements,
                    [read1, read2],
                    maximum_edit_distance,
                    minimum_template_span,
                    maximum_template_span,
                    true,
                    &mut self.origin_pair_evidence,
                    &mut self.best_pairs,
                )
            } else {
                select_best_pairs_with_fallback_score(
                    &self.mate1.placements,
                    &self.mate2.placements,
                    maximum_edit_distance,
                    minimum_template_span,
                    maximum_template_span,
                    &mut self.best_pairs,
                )
            };
            collapse_equivalent_pair_origins(&mut self.best_pairs, read1.len(), read2.len());
            if self.prefer_minimum_net_gap && self.best_pairs.len() > 1 {
                prefer_minimum_net_gap_representative(
                    &mut self.best_pairs,
                    read1.len(),
                    read2.len(),
                );
            }
        }
        let exact_retained_alternative = if semi_global_attempted && self.best_pairs.len() == 1 {
            self.exact_retained_pair_has_alternative(
                reference,
                read1,
                read2,
                self.best_pairs[0],
                minimum_template_span,
                maximum_template_span,
            )?
        } else {
            false
        };
        let class = match (self.best_pairs.len(), exact_retained_alternative) {
            (0, _) => PairMappingStatus::Unmapped,
            (1, false) => PairMappingStatus::Unique,
            _ => PairMappingStatus::Ambiguous,
        };
        Ok((
            class,
            PairAlignmentMetrics {
                mate1: mate1_metrics,
                mate2: mate2_metrics,
                compatible_pairs: selection.compatible_pairs,
                best_pair_placements: if exact_retained_alternative {
                    2
                } else {
                    u64::try_from(self.best_pairs.len()).unwrap_or(u64::MAX)
                },
                window_rescue_attempted,
                semi_global_attempted,
                best_pair_score: selection.best_pair_score,
                second_best_pair_score: selection.second_best_pair_score,
                near_best_pairings: selection
                    .near_best_pairings
                    .max(u64::from(exact_retained_alternative)),
                mapq_compatible_pairs: selection.mapq_compatible_pairs,
                mapq_best_pair_score: selection.mapq_best_pair_score,
                mapq_second_best_pair_score: selection.mapq_second_best_pair_score,
                mapq_near_best_pairings: selection
                    .mapq_near_best_pairings
                    .max(u64::from(exact_retained_alternative)),
                frontier_complete: true,
            },
            selection.second_best_distance,
        ))
    }

    // Both retained mates must be re-enumerated under one exact-origin proof;
    // splitting the scan would duplicate its shared anchor state.
    #[allow(clippy::too_many_lines)]
    fn exact_retained_pair_has_alternative(
        &mut self,
        reference: &ReferenceIndex,
        read1: &[Base],
        read2: &[Base],
        selected: PairedPlacement,
        minimum_template_span: u64,
        maximum_template_span: u64,
    ) -> Result<bool, AlignmentError> {
        let placements = [selected.mate1(), selected.mate2()];
        if placements.iter().any(|placement| placement.distance() != 0)
            || !placements[0].is_soft_clipped(read1.len())
                && !placements[1].is_soft_clipped(read2.len())
        {
            return Ok(false);
        }
        let ranges = [
            placements[0].retained_query_interval(read1.len()),
            placements[1].retained_query_interval(read2.len()),
        ];
        let retained = [&read1[ranges[0].clone()], &read2[ranges[1].clone()]];
        let mut first_projected = [SearchBase::A; MAX_READ_BASES];
        let mut second_projected = [SearchBase::A; MAX_READ_BASES];
        prepare_combined_search_projection(retained[0], false, &mut first_projected)?;
        prepare_combined_search_projection(retained[1], true, &mut second_projected)?;
        let projected = [
            &first_projected[..retained[0].len()],
            &second_projected[..retained[1].len()],
        ];
        let Some(first_seed) = reference
            .combined_exact_seed(projected[0])
            .map_err(|_| AlignmentError::CombinedIndex)?
        else {
            return Ok(true);
        };
        let Some(second_seed) = reference
            .combined_exact_seed(projected[1])
            .map_err(|_| AlignmentError::CombinedIndex)?
        else {
            return Ok(true);
        };
        let seeds = [first_seed, second_seed];
        for lane in 0..2 {
            if seeds[lane].matched_bases()
                != u64::try_from(retained[lane].len()).expect("bounded retained length fits u64")
            {
                return Err(AlignmentError::CombinedIndex);
            }
        }
        let hits = seeds.map(CombinedSeedMatches::exact_hit_count);
        if hits == [1, 1] {
            return Ok(false);
        }
        let anchor_lane = usize::from(hits[1] < hits[0]);
        if hits[anchor_lane] > SEMI_GLOBAL_MAX_EXACT_ANCHOR_HITS {
            return Ok(true);
        }
        let other_lane = 1 - anchor_lane;
        self.exact_anchor_candidates.clear();
        reference
            .visit_combined_seed(
                seeds[anchor_lane],
                0,
                u64::try_from(retained[anchor_lane].len())
                    .expect("bounded retained length fits u64"),
                &mut |hit| {
                    if let Some(candidate) = relabel_exact_retained_hit(hit, anchor_lane) {
                        self.exact_anchor_candidates.push(candidate);
                    }
                    true
                },
            )
            .map_err(|_| AlignmentError::CombinedIndex)?;
        if self.exact_anchor_candidates.is_empty() {
            return Ok(true);
        }

        let selected_origin = pair_origin_key(selected, read1.len(), read2.len());
        let anchors = &self.exact_anchor_candidates;
        let mut alternative = false;
        reference
            .visit_combined_seed(
                seeds[other_lane],
                0,
                u64::try_from(retained[other_lane].len())
                    .expect("bounded retained length fits u64"),
                &mut |hit| {
                    let Some(other_candidate) = relabel_exact_retained_hit(hit, other_lane) else {
                        return true;
                    };
                    let Some(other) = exact_retained_placement(
                        other_candidate,
                        placements[other_lane],
                        retained[other_lane].len(),
                    ) else {
                        return true;
                    };
                    for &anchor_candidate in anchors {
                        let Some(anchor) = exact_retained_placement(
                            anchor_candidate,
                            placements[anchor_lane],
                            retained[anchor_lane].len(),
                        ) else {
                            continue;
                        };
                        let pair = if anchor_lane == 0 {
                            exact_compatible_pair(
                                anchor,
                                other,
                                minimum_template_span,
                                maximum_template_span,
                            )
                        } else {
                            exact_compatible_pair(
                                other,
                                anchor,
                                minimum_template_span,
                                maximum_template_span,
                            )
                        };
                        if pair.is_some_and(|pair| {
                            pair_origin_key(pair, read1.len(), read2.len()) != selected_origin
                        }) {
                            alternative = true;
                            return false;
                        }
                    }
                    true
                },
            )
            .map_err(|_| AlignmentError::CombinedIndex)?;
        Ok(alternative)
    }

    fn should_affine_rescore(
        &self,
        class: PairMappingStatus,
        metrics: PairAlignmentMetrics,
        read1_len: usize,
        read2_len: usize,
    ) -> bool {
        if !matches!(class, PairMappingStatus::Ambiguous) || self.best_pairs.is_empty() {
            return false;
        }
        metrics.compatible_pairs <= 64
            && self.best_pairs.iter().any(|pair| {
                pair.mate1().is_soft_clipped(read1_len)
                    || pair.mate2().is_soft_clipped(read2_len)
                    || placement_net_gap_bases(pair.mate1(), read1_len) != 0
                    || placement_net_gap_bases(pair.mate2(), read2_len) != 0
            })
    }

    #[allow(clippy::too_many_arguments)]
    fn affine_rescore_directional_pair(
        &mut self,
        reference: &ReferenceIndex,
        read1: &[Base],
        read2: &[Base],
        original_class: PairMappingStatus,
        maximum_edit_distance: u8,
        minimum_template_span: u64,
        maximum_template_span: u64,
        mut metrics: PairAlignmentMetrics,
    ) -> Result<(PairMappingStatus, PairAlignmentMetrics, Option<u8>), AlignmentError> {
        let prior_best_count = metrics.best_pair_placements;
        let prior_retained_count = u64::try_from(self.best_pairs.len()).unwrap_or(u64::MAX);
        let concealed_or_incomplete =
            !metrics.frontier_complete || prior_best_count > prior_retained_count;

        self.mate1_affine_scores.clear();
        for &placement in &self.mate1.placements {
            self.mate1_affine_scores.push(affine_placement_score(
                reference,
                read1,
                placement,
                self.semi_global_clip_penalty,
                &mut self.affine,
            )?);
        }
        self.mate2_affine_scores.clear();
        for &placement in &self.mate2.placements {
            self.mate2_affine_scores.push(affine_placement_score(
                reference,
                read2,
                placement,
                self.semi_global_clip_penalty,
                &mut self.affine,
            )?);
        }

        let selection = select_best_pair_origins_with_affine_score(
            &self.mate1.placements,
            &self.mate1_affine_scores,
            &self.mate2.placements,
            &self.mate2_affine_scores,
            [read1, read2],
            maximum_edit_distance,
            minimum_template_span,
            maximum_template_span,
            &mut self.origin_pair_evidence,
            &mut self.best_pairs,
        );
        collapse_equivalent_pair_origins(&mut self.best_pairs, read1.len(), read2.len());
        if self.best_pairs.len() > 1 {
            prefer_minimum_net_gap_representative(&mut self.best_pairs, read1.len(), read2.len());
        }
        let insufficient_affine_separation = matches!(original_class, PairMappingStatus::Ambiguous)
            && selection
                .best_pair_score
                .zip(selection.second_best_pair_score)
                .is_some_and(|(best, second)| {
                    best.saturating_sub(second)
                        < BWA_MATCH_SCORE.saturating_add(BWA_MISMATCH_PENALTY)
                });
        let must_remain_ambiguous = concealed_or_incomplete || insufficient_affine_separation;
        let class = match (self.best_pairs.len(), must_remain_ambiguous) {
            (0, _) => PairMappingStatus::Unmapped,
            (1, false) => PairMappingStatus::Unique,
            _ => PairMappingStatus::Ambiguous,
        };
        metrics.compatible_pairs = selection.compatible_pairs;
        metrics.best_pair_placements = if must_remain_ambiguous {
            u64::try_from(self.best_pairs.len())
                .unwrap_or(u64::MAX)
                .max(2)
        } else {
            u64::try_from(self.best_pairs.len()).unwrap_or(u64::MAX)
        };
        metrics.best_pair_score = selection.best_pair_score;
        metrics.second_best_pair_score = selection.second_best_pair_score;
        metrics.near_best_pairings = selection
            .near_best_pairings
            .max(u64::from(must_remain_ambiguous));
        metrics.mapq_compatible_pairs = selection.mapq_compatible_pairs;
        metrics.mapq_best_pair_score = selection.mapq_best_pair_score;
        metrics.mapq_second_best_pair_score = selection.mapq_second_best_pair_score;
        metrics.mapq_near_best_pairings = selection
            .mapq_near_best_pairings
            .max(u64::from(must_remain_ambiguous));
        Ok((class, metrics, None))
    }

    #[must_use]
    pub fn best_pairs(&self) -> &[PairedPlacement] {
        &self.best_pairs
    }
}

const fn empty_pair_metrics() -> PairAlignmentMetrics {
    PairAlignmentMetrics {
        mate1: ReadAlignmentMetrics {
            located_rows: 0,
            emitted_candidate_starts: 0,
            distinct_candidate_starts: 0,
            verified_placements: 0,
        },
        mate2: ReadAlignmentMetrics {
            located_rows: 0,
            emitted_candidate_starts: 0,
            distinct_candidate_starts: 0,
            verified_placements: 0,
        },
        compatible_pairs: 0,
        best_pair_placements: 0,
        window_rescue_attempted: false,
        semi_global_attempted: false,
        best_pair_score: None,
        second_best_pair_score: None,
        near_best_pairings: 0,
        mapq_compatible_pairs: 0,
        mapq_best_pair_score: None,
        mapq_second_best_pair_score: None,
        mapq_near_best_pairings: 0,
        frontier_complete: false,
    }
}

fn select_combined_window_rescue_anchor(
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

fn nominal_pair_geometry_exists(
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

fn retain_nominal_pair_geometry(
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

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
struct PairSelection {
    compatible_pairs: u64,
    second_best_distance: Option<u8>,
    best_pair_score: Option<i16>,
    second_best_pair_score: Option<i16>,
    near_best_pairings: u64,
    mapq_compatible_pairs: u64,
    mapq_best_pair_score: Option<i16>,
    mapq_second_best_pair_score: Option<i16>,
    mapq_near_best_pairings: u64,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
struct PairScoreConfidence {
    best: Option<i16>,
    second: Option<i16>,
    counts_by_delta: [u64; BWA_NEAR_SUBOPTIMAL_DELTA as usize + 1],
}

impl PairScoreConfidence {
    fn observe(&mut self, score: i16) {
        let Some(best) = self.best else {
            self.best = Some(score);
            self.counts_by_delta[0] = 1;
            return;
        };
        if score > best {
            self.second = Some(self.second.map_or(best, |second| second.max(best)));
            let shift = usize::try_from(score - best).unwrap_or(usize::MAX);
            if shift >= self.counts_by_delta.len() {
                self.counts_by_delta.fill(0);
            } else {
                let retained = self.counts_by_delta.len() - shift;
                self.counts_by_delta.copy_within(..retained, shift);
                self.counts_by_delta[..shift].fill(0);
            }
            self.counts_by_delta[0] = 1;
            self.best = Some(score);
        } else {
            let delta = usize::try_from(best - score).unwrap_or(usize::MAX);
            if delta < self.counts_by_delta.len() {
                self.counts_by_delta[delta] = self.counts_by_delta[delta].saturating_add(1);
            }
            if score < best {
                self.second = Some(self.second.map_or(score, |second| second.max(score)));
            }
        }
    }

    fn near_best_alternatives(self) -> u64 {
        self.counts_by_delta
            .into_iter()
            .fold(0_u64, u64::saturating_add)
            .saturating_sub(u64::from(self.best.is_some()))
    }
}

fn affine_placement_score(
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

fn select_best_pairs(
    mate1: &[ReadPlacement],
    mate2: &[ReadPlacement],
    maximum_edit_distance: u8,
    minimum_span: u64,
    maximum_span: u64,
    best_pairs: &mut Vec<PairedPlacement>,
) -> PairSelection {
    select_best_pairs_with_objective(
        mate1,
        mate2,
        maximum_edit_distance,
        minimum_span,
        maximum_span,
        false,
        best_pairs,
    )
}

/// Returns false only when every placement is proven to have a distinct
/// biological origin. The common one-placement frontier exits immediately;
/// large repeat frontiers conservatively use grouped scoring instead of
/// paying for a separate allocation or quadratic duplicate proof.
fn placements_may_share_origin(placements: &[ReadPlacement], read_length: usize) -> bool {
    const EXACT_SCAN_LIMIT: usize = 64;
    if placements.len() < 2 {
        return false;
    }
    if placements.len() > EXACT_SCAN_LIMIT {
        return true;
    }
    placements.iter().enumerate().any(|(index, placement)| {
        let origin = placement_origin_key(*placement, read_length);
        placements[..index]
            .iter()
            .any(|previous| placement_origin_key(*previous, read_length) == origin)
    })
}

/// Preserves the established mapping selection while collapsing MAPQ evidence
/// to distinct biological origins. Raw endpoint counts and scores continue to
/// drive search control; only the dedicated `mapq_*` result fields are grouped.
// Endpoint selection and origin-collapsed MAPQ evidence must observe the exact
// same compatible-pair stream.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn select_best_pair_origins_with_endpoint_policy(
    mate1: &[ReadPlacement],
    mate2: &[ReadPlacement],
    reads: [&[Base]; 2],
    maximum_edit_distance: u8,
    minimum_span: u64,
    maximum_span: u64,
    fallback_scoring: bool,
    origin_evidence: &mut std::collections::HashMap<OriginPairStorageKey, OriginPairEvidence>,
    best_pairs: &mut Vec<PairedPlacement>,
) -> PairSelection {
    origin_evidence.clear();
    if !placements_may_share_origin(mate1, reads[0].len())
        && !placements_may_share_origin(mate2, reads[1].len())
    {
        return if fallback_scoring {
            select_best_pairs_with_fallback_score(
                mate1,
                mate2,
                maximum_edit_distance,
                minimum_span,
                maximum_span,
                best_pairs,
            )
        } else {
            select_best_pairs(
                mate1,
                mate2,
                maximum_edit_distance,
                minimum_span,
                maximum_span,
                best_pairs,
            )
        };
    }
    best_pairs.clear();
    let mut best_score = u8::MAX;
    let mut best_distance = u8::MAX;
    let mut second_best_score = u8::MAX;
    let mut compatible_count = 0_u64;
    let mut raw_confidence = PairScoreConfidence::default();
    for &first in mate1 {
        if first.distance() > maximum_edit_distance {
            continue;
        }
        let Some(expected) = expected_mate2_strand(first.strand()) else {
            continue;
        };
        let lower_start = first.end().saturating_sub(maximum_span);
        let upper_start = first.start().saturating_add(maximum_span);
        let lower = mate2.partition_point(|second| {
            spatial_key(*second) < (first.contig_ordinal(), expected, lower_start, 0, 0)
        });
        let upper = mate2.partition_point(|second| {
            spatial_key(*second)
                <= (
                    first.contig_ordinal(),
                    expected,
                    upper_start,
                    u64::MAX,
                    u8::MAX,
                )
        });
        for &second in &mate2[lower..upper] {
            if second.distance() > maximum_edit_distance {
                continue;
            }
            let template_start = first.start().min(second.start());
            let template_end = first.end().max(second.end());
            let span = template_end.saturating_sub(template_start);
            if !(minimum_span..=maximum_span).contains(&span) || !is_inward(first, second) {
                continue;
            }
            compatible_count = compatible_count.saturating_add(1);
            let distance = first.distance().saturating_add(second.distance());
            let score = if fallback_scoring {
                first.fallback_score.saturating_add(second.fallback_score)
            } else {
                distance
            };
            let confidence_score = -i16::from(score) * BWA_MISMATCH_PENALTY;
            raw_confidence.observe(confidence_score);
            let pair = PairedPlacement {
                mate1: first,
                mate2: second,
                template_start,
                template_end,
                distance,
                score,
            };
            let origin = pair_origin_storage_key(pair, reads[0].len(), reads[1].len());
            match origin_evidence.entry(origin) {
                std::collections::hash_map::Entry::Vacant(entry) => {
                    entry.insert(OriginPairEvidence {
                        mapq_score: confidence_score,
                    });
                }
                std::collections::hash_map::Entry::Occupied(mut entry) => {
                    let evidence = entry.get_mut();
                    if confidence_score > evidence.mapq_score {
                        evidence.mapq_score = confidence_score;
                    }
                }
            }
            let objective_distance = if fallback_scoring { distance } else { 0 };
            if (score, objective_distance) < (best_score, best_distance) {
                second_best_score = best_score;
                best_score = score;
                best_distance = objective_distance;
                best_pairs.clear();
                best_pairs.push(pair);
            } else if (score, objective_distance) == (best_score, best_distance) {
                best_pairs.push(pair);
            } else if score < second_best_score {
                second_best_score = score;
            }
        }
    }

    best_pairs.sort_unstable();
    best_pairs.dedup();
    let mut mapq_confidence = PairScoreConfidence::default();
    for evidence in origin_evidence.values() {
        mapq_confidence.observe(evidence.mapq_score);
    }

    PairSelection {
        compatible_pairs: compatible_count,
        second_best_distance: (second_best_score != u8::MAX).then_some(second_best_score),
        best_pair_score: raw_confidence.best,
        second_best_pair_score: raw_confidence.second,
        near_best_pairings: raw_confidence.near_best_alternatives(),
        mapq_compatible_pairs: u64::try_from(origin_evidence.len()).unwrap_or(u64::MAX),
        mapq_best_pair_score: mapq_confidence.best,
        mapq_second_best_pair_score: mapq_confidence.second,
        mapq_near_best_pairings: mapq_confidence.near_best_alternatives(),
    }
}

fn select_best_pairs_with_fallback_score(
    mate1: &[ReadPlacement],
    mate2: &[ReadPlacement],
    maximum_edit_distance: u8,
    minimum_span: u64,
    maximum_span: u64,
    best_pairs: &mut Vec<PairedPlacement>,
) -> PairSelection {
    select_best_pairs_with_objective(
        mate1,
        mate2,
        maximum_edit_distance,
        minimum_span,
        maximum_span,
        true,
        best_pairs,
    )
}

#[allow(clippy::too_many_arguments)]
fn select_best_pairs_with_objective(
    mate1: &[ReadPlacement],
    mate2: &[ReadPlacement],
    maximum_edit_distance: u8,
    minimum_span: u64,
    maximum_span: u64,
    fallback_scoring: bool,
    best_pairs: &mut Vec<PairedPlacement>,
) -> PairSelection {
    best_pairs.clear();
    let mut best_score = u8::MAX;
    let mut best_distance = u8::MAX;
    let mut second_best_score = u8::MAX;
    let mut compatible_count = 0_u64;
    let mut confidence = PairScoreConfidence::default();
    for &first in mate1 {
        if first.distance() > maximum_edit_distance {
            continue;
        }
        let Some(expected) = expected_mate2_strand(first.strand()) else {
            continue;
        };
        let lower_start = first.end().saturating_sub(maximum_span);
        let upper_start = first.start().saturating_add(maximum_span);
        let lower = mate2.partition_point(|second| {
            spatial_key(*second) < (first.contig_ordinal(), expected, lower_start, 0, 0)
        });
        let upper = mate2.partition_point(|second| {
            spatial_key(*second)
                <= (
                    first.contig_ordinal(),
                    expected,
                    upper_start,
                    u64::MAX,
                    u8::MAX,
                )
        });
        for &second in &mate2[lower..upper] {
            if second.distance() > maximum_edit_distance {
                continue;
            }
            let template_start = first.start().min(second.start());
            let template_end = first.end().max(second.end());
            let span = template_end - template_start;
            if !(minimum_span..=maximum_span).contains(&span) || !is_inward(first, second) {
                continue;
            }
            compatible_count = compatible_count.saturating_add(1);
            let distance = first.distance().saturating_add(second.distance());
            let score = if fallback_scoring {
                first.fallback_score.saturating_add(second.fallback_score)
            } else {
                distance
            };
            confidence.observe(-i16::from(score) * BWA_MISMATCH_PENALTY);
            let pair = PairedPlacement {
                mate1: first,
                mate2: second,
                template_start,
                template_end,
                distance,
                score,
            };
            let objective_distance = if fallback_scoring { distance } else { 0 };
            if (score, objective_distance) < (best_score, best_distance) {
                second_best_score = best_score;
                best_score = score;
                best_distance = objective_distance;
                best_pairs.clear();
                best_pairs.push(pair);
            } else if (score, objective_distance) == (best_score, best_distance) {
                best_pairs.push(pair);
            } else if score < second_best_score {
                second_best_score = score;
            }
        }
    }
    best_pairs.sort_unstable();
    best_pairs.dedup();
    PairSelection {
        compatible_pairs: compatible_count,
        second_best_distance: (second_best_score != u8::MAX).then_some(second_best_score),
        best_pair_score: confidence.best,
        second_best_pair_score: confidence.second,
        near_best_pairings: confidence.near_best_alternatives(),
        mapq_compatible_pairs: compatible_count,
        mapq_best_pair_score: confidence.best,
        mapq_second_best_pair_score: confidence.second,
        mapq_near_best_pairings: confidence.near_best_alternatives(),
    }
}

#[allow(clippy::too_many_arguments)]
fn select_best_pairs_with_affine_score(
    mate1: &[ReadPlacement],
    mate1_scores: &[i16],
    mate2: &[ReadPlacement],
    mate2_scores: &[i16],
    maximum_edit_distance: u8,
    minimum_span: u64,
    maximum_span: u64,
    best_pairs: &mut Vec<PairedPlacement>,
) -> PairSelection {
    debug_assert_eq!(mate1.len(), mate1_scores.len());
    debug_assert_eq!(mate2.len(), mate2_scores.len());
    best_pairs.clear();
    let mut best_score = i16::MIN;
    let mut compatible_count = 0_u64;
    let mut confidence = PairScoreConfidence::default();
    for (first_index, &first) in mate1.iter().enumerate() {
        if first.distance() > maximum_edit_distance {
            continue;
        }
        let Some(expected) = expected_mate2_strand(first.strand()) else {
            continue;
        };
        let lower_start = first.end().saturating_sub(maximum_span);
        let upper_start = first.start().saturating_add(maximum_span);
        let lower = mate2.partition_point(|second| {
            spatial_key(*second) < (first.contig_ordinal(), expected, lower_start, 0, 0)
        });
        let upper = mate2.partition_point(|second| {
            spatial_key(*second)
                <= (
                    first.contig_ordinal(),
                    expected,
                    upper_start,
                    u64::MAX,
                    u8::MAX,
                )
        });
        for second_index in lower..upper {
            let second = mate2[second_index];
            if second.distance() > maximum_edit_distance {
                continue;
            }
            let template_start = first.start().min(second.start());
            let template_end = first.end().max(second.end());
            let span = template_end - template_start;
            if !(minimum_span..=maximum_span).contains(&span) || !is_inward(first, second) {
                continue;
            }
            compatible_count = compatible_count.saturating_add(1);
            let pair_score = mate1_scores[first_index].saturating_add(mate2_scores[second_index]);
            confidence.observe(pair_score);
            let distance = first.distance().saturating_add(second.distance());
            let pair = PairedPlacement {
                mate1: first,
                mate2: second,
                template_start,
                template_end,
                distance,
                score: first.fallback_score.saturating_add(second.fallback_score),
            };
            if pair_score > best_score {
                best_score = pair_score;
                best_pairs.clear();
                best_pairs.push(pair);
            } else if pair_score == best_score {
                best_pairs.push(pair);
            }
        }
    }
    best_pairs.sort_unstable();
    best_pairs.dedup();
    PairSelection {
        compatible_pairs: compatible_count,
        second_best_distance: None,
        best_pair_score: confidence.best,
        second_best_pair_score: confidence.second,
        near_best_pairings: confidence.near_best_alternatives(),
        mapq_compatible_pairs: compatible_count,
        mapq_best_pair_score: confidence.best,
        mapq_second_best_pair_score: confidence.second,
        mapq_near_best_pairings: confidence.near_best_alternatives(),
    }
}

// Affine endpoint selection and origin-collapsed MAPQ evidence must observe
// the exact same compatible-pair stream.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn select_best_pair_origins_with_affine_score(
    mate1: &[ReadPlacement],
    mate1_scores: &[i16],
    mate2: &[ReadPlacement],
    mate2_scores: &[i16],
    reads: [&[Base]; 2],
    maximum_edit_distance: u8,
    minimum_span: u64,
    maximum_span: u64,
    origin_evidence: &mut std::collections::HashMap<OriginPairStorageKey, OriginPairEvidence>,
    best_pairs: &mut Vec<PairedPlacement>,
) -> PairSelection {
    debug_assert_eq!(mate1.len(), mate1_scores.len());
    debug_assert_eq!(mate2.len(), mate2_scores.len());
    origin_evidence.clear();
    if !placements_may_share_origin(mate1, reads[0].len())
        && !placements_may_share_origin(mate2, reads[1].len())
    {
        return select_best_pairs_with_affine_score(
            mate1,
            mate1_scores,
            mate2,
            mate2_scores,
            maximum_edit_distance,
            minimum_span,
            maximum_span,
            best_pairs,
        );
    }
    best_pairs.clear();
    let mut best_score = i16::MIN;
    let mut compatible_count = 0_u64;
    let mut raw_confidence = PairScoreConfidence::default();
    for (first_index, &first) in mate1.iter().enumerate() {
        if first.distance() > maximum_edit_distance {
            continue;
        }
        let Some(expected) = expected_mate2_strand(first.strand()) else {
            continue;
        };
        let lower_start = first.end().saturating_sub(maximum_span);
        let upper_start = first.start().saturating_add(maximum_span);
        let lower = mate2.partition_point(|second| {
            spatial_key(*second) < (first.contig_ordinal(), expected, lower_start, 0, 0)
        });
        let upper = mate2.partition_point(|second| {
            spatial_key(*second)
                <= (
                    first.contig_ordinal(),
                    expected,
                    upper_start,
                    u64::MAX,
                    u8::MAX,
                )
        });
        for second_index in lower..upper {
            let second = mate2[second_index];
            if second.distance() > maximum_edit_distance {
                continue;
            }
            let template_start = first.start().min(second.start());
            let template_end = first.end().max(second.end());
            let span = template_end.saturating_sub(template_start);
            if !(minimum_span..=maximum_span).contains(&span) || !is_inward(first, second) {
                continue;
            }
            compatible_count = compatible_count.saturating_add(1);
            let pair_score = mate1_scores[first_index].saturating_add(mate2_scores[second_index]);
            raw_confidence.observe(pair_score);
            let distance = first.distance().saturating_add(second.distance());
            let pair = PairedPlacement {
                mate1: first,
                mate2: second,
                template_start,
                template_end,
                distance,
                score: first.fallback_score.saturating_add(second.fallback_score),
            };
            let origin = pair_origin_storage_key(pair, reads[0].len(), reads[1].len());
            match origin_evidence.entry(origin) {
                std::collections::hash_map::Entry::Vacant(entry) => {
                    entry.insert(OriginPairEvidence {
                        mapq_score: pair_score,
                    });
                }
                std::collections::hash_map::Entry::Occupied(mut entry) => {
                    let evidence = entry.get_mut();
                    if pair_score > evidence.mapq_score {
                        evidence.mapq_score = pair_score;
                    }
                }
            }
            if pair_score > best_score {
                best_score = pair_score;
                best_pairs.clear();
                best_pairs.push(pair);
            } else if pair_score == best_score {
                best_pairs.push(pair);
            }
        }
    }

    best_pairs.sort_unstable();
    best_pairs.dedup();
    let mut mapq_confidence = PairScoreConfidence::default();
    for evidence in origin_evidence.values() {
        mapq_confidence.observe(evidence.mapq_score);
    }
    PairSelection {
        compatible_pairs: compatible_count,
        second_best_distance: None,
        best_pair_score: raw_confidence.best,
        second_best_pair_score: raw_confidence.second,
        near_best_pairings: raw_confidence.near_best_alternatives(),
        mapq_compatible_pairs: u64::try_from(origin_evidence.len()).unwrap_or(u64::MAX),
        mapq_best_pair_score: mapq_confidence.best,
        mapq_second_best_pair_score: mapq_confidence.second,
        mapq_near_best_pairings: mapq_confidence.near_best_alternatives(),
    }
}

fn collapse_equivalent_pair_origins(
    best_pairs: &mut Vec<PairedPlacement>,
    mate1_read_length: usize,
    mate2_read_length: usize,
) {
    if best_pairs.len() < 2 {
        return;
    }
    best_pairs.sort_unstable_by_key(|pair| {
        (
            pair_origin_key(*pair, mate1_read_length, mate2_read_length),
            *pair,
        )
    });
    best_pairs.dedup_by_key(|pair| pair_origin_key(*pair, mate1_read_length, mate2_read_length));
    best_pairs.sort_unstable();
}

fn pair_origin_key(
    pair: PairedPlacement,
    mate1_read_length: usize,
    mate2_read_length: usize,
) -> ((u64, BisulfiteStrand, i128), (u64, BisulfiteStrand, i128)) {
    (
        placement_origin_key(pair.mate1(), mate1_read_length),
        placement_origin_key(pair.mate2(), mate2_read_length),
    )
}

fn pair_origin_storage_key(
    pair: PairedPlacement,
    mate1_read_length: usize,
    mate2_read_length: usize,
) -> OriginPairStorageKey {
    let encode = |(contig, strand, five_prime): (u64, BisulfiteStrand, i128)| {
        (
            contig,
            u8::try_from(strand_index(strand)).expect("four strands fit u8"),
            five_prime,
        )
    };
    let (mate1, mate2) = pair_origin_key(pair, mate1_read_length, mate2_read_length);
    (encode(mate1), encode(mate2))
}

fn sequencing_three_prime_adapter_supported(read: &[Base], retained_end: usize) -> bool {
    let clipped = read.get(retained_end..).unwrap_or_default();
    let supported = clipped.len().min(ILLUMINA_UNIVERSAL_ADAPTER.len());
    supported >= ORIGIN_ENDPOINT_MIN_ADAPTER_SUPPORT
        && clipped
            .iter()
            .take(supported)
            .zip(ILLUMINA_UNIVERSAL_ADAPTER.iter().take(supported))
            .all(|(observed, expected)| observed.as_ascii() == *expected)
}

#[must_use]
fn supported_three_prime_adapter_start(read: &[Base]) -> Option<usize> {
    let earliest = read.len().saturating_sub(SEMI_GLOBAL_MAX_CLIP_BASES);
    let latest = read
        .len()
        .checked_sub(ORIGIN_ENDPOINT_MIN_ADAPTER_SUPPORT)?;
    (earliest..=latest).find(|&start| sequencing_three_prime_adapter_supported(read, start))
}

fn read_has_supported_three_prime_adapter(read: &[Base]) -> bool {
    supported_three_prime_adapter_start(read).is_some()
}

fn affine_terminal_clip_cost(length: usize, adapter_supported: bool) -> u16 {
    if length == 0 {
        return 0;
    }
    let (open, extension) = if adapter_supported {
        (
            ORIGIN_ENDPOINT_ADAPTER_CLIP_OPEN_PENALTY,
            ORIGIN_ENDPOINT_ADAPTER_CLIP_EXTENSION_PENALTY,
        )
    } else {
        (
            ORIGIN_ENDPOINT_CLIP_OPEN_PENALTY,
            ORIGIN_ENDPOINT_CLIP_EXTENSION_PENALTY,
        )
    };
    open.saturating_add(extension.saturating_mul(u16::try_from(length - 1).unwrap_or(u16::MAX)))
}

fn placement_endpoint_cost(read: &[Base], placement: ReadPlacement) -> u16 {
    let retained = placement.retained_query_interval(read.len());
    let five_prime_clip = retained.start;
    let three_prime_clip = read.len().saturating_sub(retained.end);
    u16::from(placement.distance())
        .saturating_mul(u16::from(SEMI_GLOBAL_EDIT_PENALTY))
        .saturating_add(affine_terminal_clip_cost(five_prime_clip, false))
        .saturating_add(affine_terminal_clip_cost(
            three_prime_clip,
            sequencing_three_prime_adapter_supported(read, retained.end),
        ))
}

fn pair_endpoint_key(
    reads: [&[Base]; 2],
    pair: PairedPlacement,
) -> (u16, u64, u8, u64, PairedPlacement) {
    let retained = [
        pair.mate1().retained_query_interval(reads[0].len()),
        pair.mate2().retained_query_interval(reads[1].len()),
    ];
    let clipped = reads[0]
        .len()
        .saturating_sub(retained[0].end.saturating_sub(retained[0].start))
        .saturating_add(
            reads[1]
                .len()
                .saturating_sub(retained[1].end.saturating_sub(retained[1].start)),
        );
    (
        placement_endpoint_cost(reads[0], pair.mate1())
            .saturating_add(placement_endpoint_cost(reads[1], pair.mate2())),
        u64::try_from(clipped).unwrap_or(u64::MAX),
        pair.distance(),
        placement_net_gap_bases(pair.mate1(), reads[0].len())
            .saturating_add(placement_net_gap_bases(pair.mate2(), reads[1].len())),
        pair,
    )
}

const fn spatial_key(placement: ReadPlacement) -> (u64, BisulfiteStrand, u64, u64, u8) {
    (
        placement.contig_ordinal,
        placement.strand,
        placement.start,
        placement.end,
        placement.distance,
    )
}

const fn expected_mate2_strand(strand: BisulfiteStrand) -> Option<BisulfiteStrand> {
    match strand {
        BisulfiteStrand::OT => Some(BisulfiteStrand::CTOT),
        BisulfiteStrand::OB => Some(BisulfiteStrand::CTOB),
        BisulfiteStrand::CTOT | BisulfiteStrand::CTOB => None,
    }
}

const fn counterpart_strand(
    anchor: BisulfiteStrand,
    rescuing_mate1: bool,
) -> Option<BisulfiteStrand> {
    match (anchor, rescuing_mate1) {
        (BisulfiteStrand::CTOT, true) => Some(BisulfiteStrand::OT),
        (BisulfiteStrand::CTOB, true) => Some(BisulfiteStrand::OB),
        (BisulfiteStrand::OT, false) => Some(BisulfiteStrand::CTOT),
        (BisulfiteStrand::OB, false) => Some(BisulfiteStrand::CTOB),
        _ => None,
    }
}

const fn is_inward(mate1: ReadPlacement, mate2: ReadPlacement) -> bool {
    match (mate1.strand, mate2.strand) {
        (BisulfiteStrand::OT, BisulfiteStrand::CTOT) => mate1.start < mate2.end,
        (BisulfiteStrand::OB, BisulfiteStrand::CTOB) => mate2.start < mate1.end,
        _ => false,
    }
}

/// Emits a complete flexible candidate proof frontier inside one mate window.
///
/// `maximum_edit_distance + 1` disjoint query blocks guarantee that every
/// in-budget gapped placement leaves at least one exact projected block.  The
/// exact block can be displaced from the true alignment origin by at most the
/// edit budget, which is precisely the start domain covered by the flexible
/// verifier.  This is used only after a whole-genome rescue interval exceeds
/// its locate cap, so work is proportional to the paired genomic window rather
/// than to the global repeat count.
fn append_local_flexible_proof_candidates(
    read: &[Base],
    contig: &[Base],
    window: MateRescueWindow,
    maximum_edit_distance: u8,
    candidates: &mut Vec<ReadCandidate>,
) {
    append_scalar_local_flexible_proof_candidates(
        read,
        contig,
        window,
        maximum_edit_distance,
        candidates,
    );
}

fn append_scalar_local_flexible_proof_candidates(
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
            let Some(code) = projected_code(query_base, semantics.cytosine_strand()) else {
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
        let mut reference_code = pack_projected(
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
                            proof_mask: FLEXIBLE_NOMINAL_PROOF | (1_u8 << block_ordinal),
                        });
                    }
                }
            }
            if position == last {
                break;
            }
            position += 1;
            let incoming = projected_code(
                contig[position + block_len - 1],
                semantics.cytosine_strand(),
            )
            .unwrap_or(3);
            reference_code = ((reference_code << 2) & mask) | u128::from(incoming);
        }
        query_start = query_end;
    }
}

/// Builds a rarity-first disjoint-block seed set for one mate.
///
/// The balanced `d + 1` partition is inspected before any suffix-array rows are
/// located. If it is empty or over limit, a small boundary-state dynamic
/// program chooses the complete shifted partition with the smallest
/// FM-interval sum. Every path remains a disjoint cover and independently
/// retains the pigeonhole completeness guarantee.
///
/// A repetitive block must not make a different, informative block unusable.
/// Retaining only the rarest prefix also keeps the global locate and verifier
/// work bounded by `SENSITIVE_RANKED_BLOCK_HITS` per attempted anchor.
fn shifted_ranked_boundary(nominal: usize, shift: i8) -> Option<usize> {
    if shift.is_negative() {
        nominal.checked_sub(usize::from(shift.unsigned_abs()))
    } else {
        nominal.checked_add(usize::from(shift.unsigned_abs()))
    }
}

fn ranked_block_boundaries(
    read_len: usize,
    block_count: usize,
    boundary_shifts: [i8; SENSITIVE_PROOF_BLOCKS - 1],
    minimum_block_bases: usize,
) -> Option<[usize; SENSITIVE_PROOF_BLOCKS + 1]> {
    if block_count == 0 || block_count > SENSITIVE_PROOF_BLOCKS {
        return None;
    }
    let short = read_len / block_count;
    let long_count = read_len % block_count;
    let mut boundaries = [0_usize; SENSITIVE_PROOF_BLOCKS + 1];
    for ordinal in 1..block_count {
        let nominal = short
            .checked_mul(ordinal)?
            .checked_add(long_count.min(ordinal))?;
        boundaries[ordinal] = shifted_ranked_boundary(nominal, boundary_shifts[ordinal - 1])?;
    }
    boundaries[block_count] = read_len;
    for ordinal in 0..block_count {
        if boundaries[ordinal + 1].checked_sub(boundaries[ordinal])? < minimum_block_bases {
            return None;
        }
    }
    Some(boundaries)
}

#[allow(clippy::too_many_arguments)]
fn ranked_block_seed_for_range(
    reference: &ReferenceIndex,
    read: &[Base],
    projected_search: &[SearchBase; MAX_READ_BASES],
    mate1: bool,
    block_ordinal: usize,
    query_start: usize,
    query_end: usize,
) -> Result<Option<RankedBlockSeed>, AlignmentError> {
    let block_len = query_end - query_start;
    let source = if mate1 {
        &read[query_start..query_end]
    } else {
        &read[read.len() - query_end..read.len() - query_start]
    };
    if source.contains(&Base::N) {
        return Ok(None);
    }
    let reversed_start = read.len() - query_end;
    let reversed_end = read.len() - query_start;
    let pattern = &projected_search[reversed_start..reversed_end];
    let Some(matches) = reference
        .combined_exact_seed(pattern)
        .map_err(|_| AlignmentError::CombinedIndex)?
    else {
        return Ok(None);
    };
    if matches.matched_bases() != u64::try_from(block_len).expect("bounded block length fits u64") {
        return Err(AlignmentError::CombinedIndex);
    }
    Ok(Some(RankedBlockSeed {
        matches,
        query_offset: u64::try_from(query_start).expect("bounded query offset fits u64"),
        proof_mask: 1_u8 << block_ordinal,
    }))
}

#[allow(clippy::too_many_arguments)]
fn fill_ranked_block_seed_partition(
    reference: &ReferenceIndex,
    read: &[Base],
    projected_search: &[SearchBase; MAX_READ_BASES],
    mate1: bool,
    block_count: usize,
    boundaries: &[usize; SENSITIVE_PROOF_BLOCKS + 1],
    output: &mut [Option<RankedBlockSeed>; SENSITIVE_PROOF_BLOCKS],
) -> Result<u64, AlignmentError> {
    output.fill(None);
    for (block_ordinal, slot) in output[..block_count].iter_mut().enumerate() {
        let query_start = boundaries[block_ordinal];
        let query_end = boundaries[block_ordinal + 1];
        *slot = ranked_block_seed_for_range(
            reference,
            read,
            projected_search,
            mate1,
            block_ordinal,
            query_start,
            query_end,
        )?;
    }
    Ok(output[..block_count]
        .iter()
        .flatten()
        .fold(0_u64, |total, seed| {
            total.saturating_add(seed.matches.exact_hit_count())
        }))
}

// This small dynamic program deliberately keeps edge evaluation and
// predecessor reconstruction together so completeness is easy to audit.
#[allow(clippy::too_many_lines)]
fn optimal_adaptive_ranked_block_partition(
    reference: &ReferenceIndex,
    read: &[Base],
    projected_search: &[SearchBase; MAX_READ_BASES],
    mate1: bool,
    block_count: usize,
    maximum_ranked_block_hits: u64,
) -> Result<RankedBlockPartition, AlignmentError> {
    const STATES: usize = SENSITIVE_ADAPTIVE_BOUNDARY_SHIFTS.len();
    if !(2..=SENSITIVE_PROOF_BLOCKS).contains(&block_count)
        || read.len() < SENSITIVE_ADAPTIVE_MIN_BLOCK_BASES * block_count
    {
        return Ok(None);
    }

    let short = read.len() / block_count;
    let long_count = read.len() % block_count;
    let mut positions = [[0_usize; STATES]; SENSITIVE_PROOF_BLOCKS + 1];
    positions[block_count][0] = read.len();
    for (ordinal, layer) in positions.iter_mut().enumerate().take(block_count).skip(1) {
        let nominal = short * ordinal + long_count.min(ordinal);
        for (position, &shift) in layer.iter_mut().zip(&SENSITIVE_ADAPTIVE_BOUNDARY_SHIFTS) {
            *position = shifted_ranked_boundary(nominal, shift)
                .expect("qualified adaptive boundary remains inside the read");
        }
    }

    let mut edges = [[[None::<Option<RankedBlockSeed>>; STATES]; STATES]; SENSITIVE_PROOF_BLOCKS];
    let mut costs = [[u64::MAX; 2]; STATES];
    costs[0][0] = 0;
    let mut predecessors = [[[None::<(u8, u8)>; 2]; STATES]; SENSITIVE_PROOF_BLOCKS + 1];

    for block_ordinal in 0..block_count {
        let from_states = if block_ordinal == 0 { 1 } else { STATES };
        let to_states = if block_ordinal + 1 == block_count {
            1
        } else {
            STATES
        };
        let mut next_costs = [[u64::MAX; 2]; STATES];
        for from_state in 0..from_states {
            // Edge weights are non-negative and the caller only accepts a
            // complete partition within the locate cap. Once both hit states
            // exceed that cap, no continuation through this boundary can
            // become admissible, so avoid its FM interval queries entirely.
            if costs[from_state]
                .iter()
                .all(|&cost| cost > maximum_ranked_block_hits)
            {
                continue;
            }
            for to_state in 0..to_states {
                let query_start = positions[block_ordinal][from_state];
                let query_end = positions[block_ordinal + 1][to_state];
                if query_end <= query_start
                    || query_end - query_start < SENSITIVE_ADAPTIVE_MIN_BLOCK_BASES
                {
                    continue;
                }
                let seed = ranked_block_seed_for_range(
                    reference,
                    read,
                    projected_search,
                    mate1,
                    block_ordinal,
                    query_start,
                    query_end,
                )?;
                let hits = seed.map_or(0, |seed| seed.matches.exact_hit_count());
                edges[block_ordinal][from_state][to_state] = Some(seed);
                for (had_hits, &cost) in costs[from_state].iter().enumerate() {
                    if cost == u64::MAX {
                        continue;
                    }
                    let has_hits = usize::from(had_hits != 0 || hits != 0);
                    let candidate_cost = cost.saturating_add(hits);
                    if candidate_cost <= maximum_ranked_block_hits
                        && candidate_cost < next_costs[to_state][has_hits]
                    {
                        next_costs[to_state][has_hits] = candidate_cost;
                        predecessors[block_ordinal + 1][to_state][has_hits] = Some((
                            u8::try_from(from_state).expect("adaptive state fits u8"),
                            u8::try_from(had_hits).expect("hit state fits u8"),
                        ));
                    }
                }
            }
        }
        if next_costs
            .iter()
            .flatten()
            .all(|&cost| cost > maximum_ranked_block_hits)
        {
            return Ok(None);
        }
        costs = next_costs;
    }

    let all_hits = costs[0][1];
    if all_hits == u64::MAX {
        return Ok(None);
    }
    let mut seeds = [None; SENSITIVE_PROOF_BLOCKS];
    let mut state = 0_usize;
    let mut had_hits = 1_usize;
    for block_ordinal in (0..block_count).rev() {
        let (from_state, previous_had_hits) = predecessors[block_ordinal + 1][state][had_hits]
            .expect("reachable adaptive state has a predecessor");
        let from_state = usize::from(from_state);
        seeds[block_ordinal] =
            edges[block_ordinal][from_state][state].expect("selected adaptive edge was evaluated");
        state = from_state;
        had_hits = usize::from(previous_had_hits);
    }
    debug_assert_eq!(state, 0);
    debug_assert_eq!(had_hits, 0);
    Ok(Some((all_hits, seeds)))
}

fn collect_ranked_block_seeds(
    reference: &ReferenceIndex,
    read: &[Base],
    reversed_projected: &[ProjectedBase],
    mate1: bool,
    maximum_edit_distance: u8,
    maximum_ranked_block_hits: u64,
    output: &mut [Option<RankedBlockSeed>; SENSITIVE_PROOF_BLOCKS],
) -> Result<Option<RankedBlockSelection>, AlignmentError> {
    output.fill(None);
    let block_count = usize::from(maximum_edit_distance) + 1;
    let mut projected_search = [SearchBase::A; MAX_READ_BASES];
    for (destination, &source) in projected_search.iter_mut().zip(reversed_projected) {
        *destination = match source {
            ProjectedBase::A => SearchBase::A,
            ProjectedBase::G => SearchBase::G,
            ProjectedBase::T => SearchBase::T,
        };
    }
    let balanced_boundaries = ranked_block_boundaries(
        read.len(),
        block_count,
        SENSITIVE_BALANCED_BOUNDARY_SHIFTS,
        0,
    )
    .expect("a bounded read has a balanced d+1 partition");
    let mut best_seeds = [None; SENSITIVE_PROOF_BLOCKS];
    let mut all_hits = fill_ranked_block_seed_partition(
        reference,
        read,
        &projected_search,
        mate1,
        block_count,
        &balanced_boundaries,
        &mut best_seeds,
    )?;
    if (all_hits == 0 || all_hits > maximum_ranked_block_hits)
        && let Some((candidate_hits, candidate_seeds)) = optimal_adaptive_ranked_block_partition(
            reference,
            read,
            &projected_search,
            mate1,
            block_count,
            maximum_ranked_block_hits,
        )?
        && candidate_hits <= maximum_ranked_block_hits
    {
        best_seeds = candidate_seeds;
        all_hits = candidate_hits;
    }
    *output = best_seeds;
    output[..block_count]
        .sort_by_key(|seed| seed.map_or(u64::MAX, |seed| seed.matches.exact_hit_count()));
    let mut retained_hits = 0_u64;
    for slot in &mut output[..block_count] {
        let Some(seed) = *slot else {
            break;
        };
        let next = retained_hits.saturating_add(seed.matches.exact_hit_count());
        if next > maximum_ranked_block_hits {
            *slot = None;
            continue;
        }
        retained_hits = next;
    }
    Ok((retained_hits != 0).then_some(RankedBlockSelection {
        retained_hits,
        complete: retained_hits == all_hits,
    }))
}

fn append_ranked_block_candidates(
    reference: &ReferenceIndex,
    read_len: usize,
    mate1: bool,
    seeds: &[Option<RankedBlockSeed>; SENSITIVE_PROOF_BLOCKS],
    candidates: &mut Vec<ReadCandidate>,
) -> Result<u64, AlignmentError> {
    let query_len = u64::try_from(read_len).expect("bounded read length fits u64");
    let mut located_rows = 0_u64;
    for seed in seeds.iter().flatten() {
        let metrics = reference
            .visit_combined_seed(seed.matches, seed.query_offset, query_len, &mut |hit| {
                let strand = if mate1 {
                    hit.strand()
                } else {
                    match hit.strand() {
                        BisulfiteStrand::OT => BisulfiteStrand::CTOT,
                        BisulfiteStrand::OB => BisulfiteStrand::CTOB,
                        BisulfiteStrand::CTOT | BisulfiteStrand::CTOB => return true,
                    }
                };
                candidates.push(ReadCandidate {
                    contig_ordinal: hit.contig_ordinal(),
                    start: hit.start(),
                    strand,
                    proof_mask: FLEXIBLE_NOMINAL_PROOF | seed.proof_mask,
                });
                true
            })
            .map_err(|_| AlignmentError::CombinedIndex)?;
        located_rows = located_rows.saturating_add(metrics.located_coordinates());
    }
    Ok(located_rows)
}

fn balanced_rescue_blocks(read_len: usize) -> [ProofBlock; RESCUE_BLOCKS] {
    let short = read_len / RESCUE_BLOCKS;
    let long_count = read_len % RESCUE_BLOCKS;
    let mut cursor = 0_usize;
    core::array::from_fn(|ordinal| {
        let length = short + usize::from(ordinal < long_count);
        let start = cursor;
        cursor += length;
        ProofBlock {
            query_start: u16::try_from(start).expect("bounded rescue block start fits u16"),
            query_end: u16::try_from(cursor).expect("bounded rescue block end fits u16"),
        }
    })
}

fn pack_projected(bases: &[Base], strand: bsbit_core::bisulfite::CytosineStrand) -> u128 {
    bases.iter().fold(0_u128, |packed, &base| {
        (packed << 2) | u128::from(projected_code(base, strand).unwrap_or(3))
    })
}

const fn projected_code(base: Base, strand: bsbit_core::bisulfite::CytosineStrand) -> Option<u8> {
    use bsbit_core::bisulfite::CytosineStrand::{Bottom, Top};
    match (strand, base) {
        (_, Base::A) | (Bottom, Base::G) => Some(0),
        (Top, Base::C | Base::T) | (Bottom, Base::C) => Some(1),
        (Top, Base::G) | (Bottom, Base::T) => Some(2),
        _ => None,
    }
}

fn best_ungapped_semi_global_placement(
    reference: &ReferenceIndex,
    read: &[Base],
    candidate: ReadCandidate,
    maximum_edit_distance: u8,
    clip_penalty: u8,
) -> Option<ReadPlacement> {
    let contig = reference.contig_by_ordinal(candidate.contig_ordinal())?;
    let nominal_start = usize::try_from(candidate.start()).ok()?;
    let alignment = UngappedProfile::new(
        contig.sequence().bases(),
        nominal_start,
        read,
        candidate.strand(),
    )?
    .best_bounded_semiglobal(BoundedSemiglobalConfig::new(
        maximum_edit_distance,
        SEMI_GLOBAL_MAX_CLIP_BASES,
        SEMI_GLOBAL_MIN_ALIGNED_BASES,
        SEMI_GLOBAL_EDIT_PENALTY,
        clip_penalty,
        SEMI_GLOBAL_ADMISSION_EDIT_PENALTY,
        SEMI_GLOBAL_CLIP_PENALTY,
        u8::try_from(read.len() / 5).unwrap_or(u8::MAX),
    ))?;
    let endpoint = alignment.endpoint();
    Some(ReadPlacement {
        contig_ordinal: candidate.contig_ordinal(),
        start: u64::try_from(endpoint.reference_start()).ok()?,
        end: u64::try_from(endpoint.reference_end()).ok()?,
        strand: candidate.strand(),
        distance: endpoint.distance(),
        query_start: u16::try_from(endpoint.query_start()).ok()?,
        query_end: u16::try_from(endpoint.query_end()).ok()?,
        fallback_score: alignment.score(),
    })
}
fn best_ungapped_origin_endpoint_placement(
    reference: &ReferenceIndex,
    read: &[Base],
    candidate: ReadCandidate,
    maximum_edit_distance: u8,
    clip_penalty: u8,
) -> Option<ReadPlacement> {
    if read.len() < SEMI_GLOBAL_MIN_ALIGNED_BASES {
        return None;
    }
    let contig = reference.contig_by_ordinal(candidate.contig_ordinal())?;
    let nominal_start = usize::try_from(candidate.start()).ok()?;
    let profile = UngappedProfile::new(
        contig.sequence().bases(),
        nominal_start,
        read,
        candidate.strand(),
    )?;
    let maximum_clip =
        SEMI_GLOBAL_MAX_CLIP_BASES.min(read.len().saturating_sub(SEMI_GLOBAL_MIN_ALIGNED_BASES));
    let mut best: Option<(EndpointKey, UngappedEndpoint)> = None;
    for oriented_left_clip in 0..=maximum_clip {
        for oriented_right_clip in 0..=maximum_clip {
            let clipped = oriented_left_clip.saturating_add(oriented_right_clip);
            if read.len().saturating_sub(clipped) < SEMI_GLOBAL_MIN_ALIGNED_BASES {
                continue;
            }
            let Some(endpoint) = profile.endpoint(oriented_left_clip, oriented_right_clip) else {
                continue;
            };
            if endpoint.distance() > maximum_edit_distance {
                continue;
            }
            let admission_score = endpoint
                .distance()
                .saturating_mul(SEMI_GLOBAL_ADMISSION_EDIT_PENALTY)
                .saturating_add(
                    u8::try_from(clipped)
                        .unwrap_or(u8::MAX)
                        .saturating_mul(SEMI_GLOBAL_CLIP_PENALTY),
                );
            if admission_score > u8::try_from(read.len() / 5).unwrap_or(u8::MAX) {
                continue;
            }
            let endpoint_cost = u16::from(endpoint.distance())
                .saturating_mul(u16::from(SEMI_GLOBAL_EDIT_PENALTY))
                .saturating_add(affine_terminal_clip_cost(endpoint.query_start(), false))
                .saturating_add(affine_terminal_clip_cost(
                    read.len().saturating_sub(endpoint.query_end()),
                    sequencing_three_prime_adapter_supported(read, endpoint.query_end()),
                ));
            let fallback_score = endpoint
                .distance()
                .saturating_mul(SEMI_GLOBAL_EDIT_PENALTY)
                .saturating_add(
                    u8::try_from(clipped)
                        .unwrap_or(u8::MAX)
                        .saturating_mul(clip_penalty),
                );
            let key = (
                endpoint_cost,
                clipped,
                endpoint.distance(),
                endpoint.oriented_left_clip(),
                endpoint.oriented_right_clip(),
                fallback_score,
                endpoint.query_start(),
                endpoint.query_end(),
            );
            if best.as_ref().is_none_or(|(current, _)| key < *current) {
                best = Some((key, endpoint));
            }
        }
    }
    let ((_, _, _, _, _, fallback_score, _, _), endpoint) = best?;
    Some(ReadPlacement {
        contig_ordinal: candidate.contig_ordinal(),
        start: u64::try_from(endpoint.reference_start()).ok()?,
        end: u64::try_from(endpoint.reference_end()).ok()?,
        strand: candidate.strand(),
        distance: endpoint.distance(),
        query_start: u16::try_from(endpoint.query_start()).ok()?,
        query_end: u16::try_from(endpoint.query_end()).ok()?,
        fallback_score,
    })
}
#[cfg(test)]
#[path = "../../tests/whitebox/paired_end.rs"]
mod tests;
