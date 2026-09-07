//! Pair selection, origin grouping, endpoint policy, and affine rescoring.
//!
//! This module is a mechanical responsibility split from `paired_end`; it
//! shares the parent's private bounded-policy vocabulary and does not alter the
//! crate's public API or algorithmic ordering.

use super::mapq::{PAIR_NEAR_SUBOPTIMAL_SCORE_DELTA, confidence_score_from_edit_distance};
use crate::adapter::{
    AlignmentOutputPolicy, SoftClipMode, sequencing_three_prime_adapter_supported,
    supported_three_prime_adapter_start,
};
use crate::alignment_policy::{
    ORIGIN_ENDPOINT_ADAPTER_CLIP_EXTENSION_PENALTY, ORIGIN_ENDPOINT_ADAPTER_CLIP_OPEN_PENALTY,
    ORIGIN_ENDPOINT_CLIP_EXTENSION_PENALTY, ORIGIN_ENDPOINT_CLIP_OPEN_PENALTY,
    PAIR_ORIGIN_EXACT_SCAN_LIMIT, SEMI_GLOBAL_ADMISSION_EDIT_PENALTY, SEMI_GLOBAL_CLIP_PENALTY,
    SEMI_GLOBAL_EDIT_PENALTY, SEMI_GLOBAL_MIN_ALIGNED_BASES, SENSITIVE_CLIP_PENALTY,
};
use crate::placement::{ReadPlacement, placement_net_gap_bases, placement_origin_key};
use crate::read_mapping::{ReadCandidate, ReadWorkspace, strand_index};
use crate::search::combined_query::CombinedSeedHit;
use crate::verification::ungapped::{BoundedSemiglobalConfig, UngappedEndpoint, UngappedProfile};
use bsbit_core::alphabet::Base;
use bsbit_core::bisulfite::{AlignmentOrientation, BisulfiteStrand, strand_semantics};
use bsbit_index::reference::ReferenceIndex;

use super::result::PairedPlacement;

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
struct EndpointKey {
    cost: u16,
    clipped: usize,
    distance: u8,
    oriented_left_clip: usize,
    oriented_right_clip: usize,
    fallback_score: u8,
    query_start: usize,
    query_end: usize,
}

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub(super) struct OriginPairStorageKey {
    mate1: (u64, u8, i128),
    mate2: (u64, u8, i128),
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) struct OriginPairEvidence {
    // Larger is better. Each biological origin contributes only its best
    // endpoint score to MAPQ, while raw endpoint statistics remain available
    // to control the search pipeline.
    mapq_score: i16,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub(super) struct PairSelection {
    pub(super) compatible_pairs: u64,
    pub(super) second_best_distance: Option<u8>,
    pub(super) best_pair_score: Option<i16>,
    pub(super) second_best_pair_score: Option<i16>,
    pub(super) near_best_pairings: u64,
    pub(super) mapq_compatible_pairs: u64,
    pub(super) mapq_best_pair_score: Option<i16>,
    pub(super) mapq_second_best_pair_score: Option<i16>,
    pub(super) mapq_near_best_pairings: u64,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
struct PairScoreConfidence {
    best: Option<i16>,
    second: Option<i16>,
    counts_by_delta: [u64; PAIR_NEAR_SUBOPTIMAL_SCORE_DELTA as usize + 1],
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

/// Adds candidate-local ungapped semi-global placements without another
/// FM-index traversal. The existing seed frontier supplies full-read
/// nominal origins; this pass only chooses bounded terminal endpoints.
pub(super) fn append_ungapped_semi_global_placements(
    workspace: &mut ReadWorkspace,
    reference: &ReferenceIndex,
    read: &[Base],
    maximum_edit_distance: u8,
    clip_penalty: u8,
    output_policy: &AlignmentOutputPolicy,
) {
    for &candidate in &workspace.candidate_nominals {
        if let Some(placement) = best_ungapped_semi_global_placement_with_policy(
            reference,
            read,
            candidate,
            maximum_edit_distance,
            clip_penalty,
            output_policy,
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
    output_policy: &AlignmentOutputPolicy,
) -> ReadPlacement {
    if !placement.is_soft_clipped(read.len())
        && !read_has_supported_three_prime_adapter(read, output_policy)
    {
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
        output_policy,
    ) else {
        return placement;
    };
    if placement_origin_key(endpoint, read.len()) != placement_origin_key(placement, read.len()) {
        return placement;
    }
    if placement_endpoint_cost_with_policy(read, endpoint, output_policy)
        < placement_endpoint_cost_with_policy(read, placement, output_policy)
    {
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
pub(super) fn select_reported_origin_endpoint(
    reference: &ReferenceIndex,
    reads: [&[Base]; 2],
    selected: PairedPlacement,
    maximum_edit_distance: u8,
    minimum_template_span: u64,
    maximum_template_span: u64,
    output_policy: &AlignmentOutputPolicy,
) -> PairedPlacement {
    let mate1_variant = origin_endpoint_variant(
        reference,
        reads[0],
        selected.mate1(),
        maximum_edit_distance,
        output_policy,
    );
    let mate2_variant = origin_endpoint_variant(
        reference,
        reads[1],
        selected.mate2(),
        maximum_edit_distance,
        output_policy,
    );
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
    let mut best_key = pair_endpoint_key(reads, selected, output_policy);
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
            let key = pair_endpoint_key(reads, candidate, output_policy);
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
pub(super) fn prefer_minimum_net_gap_representative(
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

pub(super) fn select_best_pairs(
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
pub(super) fn placements_may_share_origin(
    placements: &[ReadPlacement],
    read_length: usize,
) -> bool {
    if placements.len() < 2 {
        return false;
    }
    if placements.len() > PAIR_ORIGIN_EXACT_SCAN_LIMIT {
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
pub(super) fn select_best_pair_origins_with_endpoint_policy(
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
            let confidence_score = confidence_score_from_edit_distance(score);
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

pub(super) fn select_best_pairs_with_fallback_score(
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
            confidence.observe(confidence_score_from_edit_distance(score));
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
pub(super) fn select_best_pair_origins_with_affine_score(
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

pub(super) fn collapse_equivalent_pair_origins(
    best_pairs: &mut Vec<PairedPlacement>,
    mate1_read_length: usize,
    mate2_read_length: usize,
    prefer_minimum_net_gap: bool,
) {
    if best_pairs.len() < 2 {
        return;
    }
    best_pairs.sort_unstable_by_key(|pair| {
        (
            pair_origin_key(*pair, mate1_read_length, mate2_read_length),
            if prefer_minimum_net_gap {
                placement_net_gap_bases(pair.mate1(), mate1_read_length)
                    .saturating_add(placement_net_gap_bases(pair.mate2(), mate2_read_length))
            } else {
                0
            },
            *pair,
        )
    });
    best_pairs.dedup_by_key(|pair| pair_origin_key(*pair, mate1_read_length, mate2_read_length));
    if prefer_minimum_net_gap {
        let representative = best_pairs
            .iter()
            .enumerate()
            .min_by_key(|(_, pair)| {
                (
                    placement_net_gap_bases(pair.mate1(), mate1_read_length)
                        .saturating_add(placement_net_gap_bases(pair.mate2(), mate2_read_length)),
                    **pair,
                )
            })
            .map_or(0, |(index, _)| index);
        best_pairs.swap(0, representative);
    } else {
        best_pairs.sort_unstable();
    }
}

pub(super) fn pair_origin_key(
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
    OriginPairStorageKey {
        mate1: encode(mate1),
        mate2: encode(mate2),
    }
}

fn read_has_supported_three_prime_adapter(
    read: &[Base],
    output_policy: &AlignmentOutputPolicy,
) -> bool {
    supported_three_prime_adapter_start(read, output_policy).is_some()
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

#[cfg(test)]
pub(super) fn placement_endpoint_cost(read: &[Base], placement: ReadPlacement) -> u16 {
    placement_endpoint_cost_with_policy(read, placement, &AlignmentOutputPolicy::default())
}

fn placement_endpoint_cost_with_policy(
    read: &[Base],
    placement: ReadPlacement,
    output_policy: &AlignmentOutputPolicy,
) -> u16 {
    let retained = placement.retained_query_interval(read.len());
    let five_prime_clip = retained.start;
    let three_prime_clip = read.len().saturating_sub(retained.end);
    u16::from(placement.distance())
        .saturating_mul(u16::from(SEMI_GLOBAL_EDIT_PENALTY))
        .saturating_add(affine_terminal_clip_cost(five_prime_clip, false))
        .saturating_add(affine_terminal_clip_cost(
            three_prime_clip,
            sequencing_three_prime_adapter_supported(read, retained.end, output_policy),
        ))
}

fn pair_endpoint_key(
    reads: [&[Base]; 2],
    pair: PairedPlacement,
    output_policy: &AlignmentOutputPolicy,
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
        placement_endpoint_cost_with_policy(reads[0], pair.mate1(), output_policy).saturating_add(
            placement_endpoint_cost_with_policy(reads[1], pair.mate2(), output_policy),
        ),
        u64::try_from(clipped).unwrap_or(u64::MAX),
        pair.distance(),
        placement_net_gap_bases(pair.mate1(), reads[0].len())
            .saturating_add(placement_net_gap_bases(pair.mate2(), reads[1].len())),
        pair,
    )
}

pub(super) const fn spatial_key(placement: ReadPlacement) -> (u64, BisulfiteStrand, u64, u64, u8) {
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

pub(super) const fn counterpart_strand(
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

#[cfg(test)]
pub(super) fn best_ungapped_semi_global_placement(
    reference: &ReferenceIndex,
    read: &[Base],
    candidate: ReadCandidate,
    maximum_edit_distance: u8,
    clip_penalty: u8,
) -> Option<ReadPlacement> {
    best_ungapped_semi_global_placement_with_policy(
        reference,
        read,
        candidate,
        maximum_edit_distance,
        clip_penalty,
        &AlignmentOutputPolicy::default(),
    )
}

fn best_ungapped_semi_global_placement_with_policy(
    reference: &ReferenceIndex,
    read: &[Base],
    candidate: ReadCandidate,
    maximum_edit_distance: u8,
    clip_penalty: u8,
    output_policy: &AlignmentOutputPolicy,
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
        output_policy.maximum_soft_clip_bases(),
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
    output_policy: &AlignmentOutputPolicy,
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
    let maximum_clip = output_policy
        .maximum_soft_clip_bases()
        .min(read.len().saturating_sub(SEMI_GLOBAL_MIN_ALIGNED_BASES));
    let mut best: Option<(EndpointKey, UngappedEndpoint)> = None;
    for oriented_left_clip in 0..=maximum_clip {
        for oriented_right_clip in 0..=maximum_clip {
            let clipped = oriented_left_clip.saturating_add(oriented_right_clip);
            if clipped > output_policy.maximum_soft_clip_bases() {
                continue;
            }
            if read.len().saturating_sub(clipped) < SEMI_GLOBAL_MIN_ALIGNED_BASES {
                continue;
            }
            let Some(endpoint) = profile.endpoint(oriented_left_clip, oriented_right_clip) else {
                continue;
            };
            if matches!(output_policy.soft_clip_mode(), SoftClipMode::Adapter)
                && (endpoint.query_start() != 0
                    || !sequencing_three_prime_adapter_supported(
                        read,
                        endpoint.query_end(),
                        output_policy,
                    ))
            {
                continue;
            }
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
                    sequencing_three_prime_adapter_supported(
                        read,
                        endpoint.query_end(),
                        output_policy,
                    ),
                ));
            let fallback_score = endpoint
                .distance()
                .saturating_mul(SEMI_GLOBAL_EDIT_PENALTY)
                .saturating_add(
                    u8::try_from(clipped)
                        .unwrap_or(u8::MAX)
                        .saturating_mul(clip_penalty),
                );
            let key = EndpointKey {
                cost: endpoint_cost,
                clipped,
                distance: endpoint.distance(),
                oriented_left_clip: endpoint.oriented_left_clip(),
                oriented_right_clip: endpoint.oriented_right_clip(),
                fallback_score,
                query_start: endpoint.query_start(),
                query_end: endpoint.query_end(),
            };
            if best.as_ref().is_none_or(|(current, _)| key < *current) {
                best = Some((key, endpoint));
            }
        }
    }
    let (key, endpoint) = best?;
    let fallback_score = key.fallback_score;
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
