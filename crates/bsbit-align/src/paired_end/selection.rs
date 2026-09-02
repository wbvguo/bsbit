//! Compatible-pair selection, biological-origin collapse, and score evidence.

use super::endpoint::{
    best_ungapped_origin_endpoint_placement, pair_endpoint_key, placement_endpoint_cost,
    read_has_supported_three_prime_adapter,
};
use super::{
    AffineScoreWorkspace, AlignmentError, AlignmentOrientation, BWA_MISMATCH_PENALTY, Base,
    BisulfiteStrand, MateRole, OriginPairEvidence, OriginPairStorageKey, PairScoreConfidence,
    PairSelection, PairedPlacement, ReadCandidate, ReadPlacement, ReferenceIndex,
    SENSITIVE_CLIP_PENALTY, banded_affine_score, placement_net_gap_bases, placement_origin_key,
    strand_index, strand_semantics,
};

pub(super) fn candidate_for_origin_endpoint(
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

pub(super) fn origin_endpoint_variant(
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
pub(super) fn select_reported_origin_endpoint(
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

pub(super) fn pair_net_gap_profile(
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

pub(super) fn affine_placement_score(
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
pub(super) fn select_best_pairs_with_objective(
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
pub(super) fn select_best_pairs_with_affine_score(
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

pub(super) fn pair_origin_storage_key(
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

pub(super) const fn spatial_key(placement: ReadPlacement) -> (u64, BisulfiteStrand, u64, u64, u8) {
    (
        placement.contig_ordinal,
        placement.strand,
        placement.start,
        placement.end,
        placement.distance,
    )
}

pub(super) const fn expected_mate2_strand(strand: BisulfiteStrand) -> Option<BisulfiteStrand> {
    match strand {
        BisulfiteStrand::OT => Some(BisulfiteStrand::CTOT),
        BisulfiteStrand::OB => Some(BisulfiteStrand::CTOB),
        BisulfiteStrand::CTOT | BisulfiteStrand::CTOB => None,
    }
}

pub(super) const fn counterpart_strand(
    anchor: BisulfiteStrand,
    rescuing_mate: MateRole,
) -> Option<BisulfiteStrand> {
    match (anchor, rescuing_mate) {
        (BisulfiteStrand::CTOT, MateRole::First) => Some(BisulfiteStrand::OT),
        (BisulfiteStrand::CTOB, MateRole::First) => Some(BisulfiteStrand::OB),
        (BisulfiteStrand::OT, MateRole::Second) => Some(BisulfiteStrand::CTOT),
        (BisulfiteStrand::OB, MateRole::Second) => Some(BisulfiteStrand::CTOB),
        _ => None,
    }
}

pub(super) const fn is_inward(mate1: ReadPlacement, mate2: ReadPlacement) -> bool {
    match (mate1.strand, mate2.strand) {
        (BisulfiteStrand::OT, BisulfiteStrand::CTOT) => mate1.start < mate2.end,
        (BisulfiteStrand::OB, BisulfiteStrand::CTOB) => mate2.start < mate1.end,
        _ => false,
    }
}
