//! Bounded endpoint-placement policy with shared adapter evidence.

use crate::adapter::sequencing_three_prime_adapter_supported;

use super::{
    Base, BoundedSemiglobalConfig, EndpointKey, ORIGIN_ENDPOINT_ADAPTER_CLIP_EXTENSION_PENALTY,
    ORIGIN_ENDPOINT_ADAPTER_CLIP_OPEN_PENALTY, ORIGIN_ENDPOINT_CLIP_EXTENSION_PENALTY,
    ORIGIN_ENDPOINT_CLIP_OPEN_PENALTY, PairedPlacement, ReadCandidate, ReadPlacement,
    ReferenceIndex, SEMI_GLOBAL_ADMISSION_EDIT_PENALTY, SEMI_GLOBAL_CLIP_PENALTY,
    SEMI_GLOBAL_EDIT_PENALTY, SEMI_GLOBAL_MAX_CLIP_BASES, SEMI_GLOBAL_MIN_ALIGNED_BASES,
    UngappedEndpoint, UngappedProfile, placement_net_gap_bases,
};

pub(super) fn affine_terminal_clip_cost(length: usize, adapter_supported: bool) -> u16 {
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

pub(super) fn placement_endpoint_cost(read: &[Base], placement: ReadPlacement) -> u16 {
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

pub(super) fn pair_endpoint_key(
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

pub(super) fn best_ungapped_semi_global_placement(
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
pub(super) fn best_ungapped_origin_endpoint_placement(
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
