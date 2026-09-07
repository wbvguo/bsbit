//! Stable, reference-order-independent lottery for unresolved reporting ties.

use core::cmp::Ordering;

use bsbit_index::reference::ReferenceIndex;

use crate::AlignmentError;
use crate::placement::{ReadPlacement, placement_origin_key};

const DOMAIN: u64 = 0x4253_4249_545f_5442;

/// One caller-supplied identity for a reproducible reporting lottery.
#[derive(Clone, Copy)]
pub(crate) struct ReportingTieBreak {
    pub(crate) seed: u64,
    pub(crate) read_key: u64,
}

#[inline]
const fn mix64(mut value: u64) -> u64 {
    value = value.wrapping_add(0x9e37_79b9_7f4a_7c15);
    value = (value ^ (value >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    value = (value ^ (value >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    value ^ (value >> 31)
}

#[inline]
const fn absorb(state: u64, value: u64) -> u64 {
    mix64(state ^ mix64(value))
}

fn absorb_bytes(mut state: u64, bytes: &[u8]) -> u64 {
    state = absorb(state, u64::try_from(bytes.len()).unwrap_or(u64::MAX));
    for chunk in bytes.chunks(8) {
        let mut word = [0_u8; 8];
        word[..chunk.len()].copy_from_slice(chunk);
        state = absorb(state, u64::from_le_bytes(word));
    }
    state
}

/// Hashes a biological single-read origin without consulting its contig ordinal.
pub(crate) fn placement_origin_hash(
    reference: &ReferenceIndex,
    tie_break: ReportingTieBreak,
    placement: ReadPlacement,
    read_length: usize,
) -> Result<u64, AlignmentError> {
    let contig = reference
        .contig_by_ordinal(placement.contig_ordinal())
        .ok_or(AlignmentError::InvalidContigOrdinal {
            ordinal: placement.contig_ordinal(),
        })?;
    let (_, strand, five_prime) = placement_origin_key(placement, read_length);
    let five_prime = five_prime.to_le_bytes();
    let mut state = absorb(absorb(DOMAIN, tie_break.seed), tie_break.read_key);
    state = absorb_bytes(state, contig.name());
    state = absorb(state, strand_code(strand));
    state = absorb_bytes(state, &five_prime);
    Ok(mix64(state))
}

/// Combines two ordered mate origins into one pair lottery value.
pub(crate) fn pair_origin_hash(
    reference: &ReferenceIndex,
    tie_break: ReportingTieBreak,
    mate1: ReadPlacement,
    mate1_read_length: usize,
    mate2: ReadPlacement,
    mate2_read_length: usize,
) -> Result<u64, AlignmentError> {
    let first = placement_origin_hash(reference, tie_break, mate1, mate1_read_length)?;
    let second = placement_origin_hash(reference, tie_break, mate2, mate2_read_length)?;
    Ok(mix64(absorb(absorb(DOMAIN ^ 0x5041_4952, first), second)))
}

/// Stable collision fallback that never compares a reference ordinal.
pub(crate) fn compare_placements_without_reference_order(
    reference: &ReferenceIndex,
    left: ReadPlacement,
    right: ReadPlacement,
    read_length: usize,
) -> Result<Ordering, AlignmentError> {
    let left_contig = reference.contig_by_ordinal(left.contig_ordinal()).ok_or(
        AlignmentError::InvalidContigOrdinal {
            ordinal: left.contig_ordinal(),
        },
    )?;
    let right_contig = reference.contig_by_ordinal(right.contig_ordinal()).ok_or(
        AlignmentError::InvalidContigOrdinal {
            ordinal: right.contig_ordinal(),
        },
    )?;
    let (_, left_strand, left_five_prime) = placement_origin_key(left, read_length);
    let (_, right_strand, right_five_prime) = placement_origin_key(right, read_length);
    Ok((
        left_contig.name(),
        strand_code(left_strand),
        left_five_prime,
        left.start,
        left.end,
        left.distance,
        left.query_start,
        left.query_end,
        left.fallback_score,
    )
        .cmp(&(
            right_contig.name(),
            strand_code(right_strand),
            right_five_prime,
            right.start,
            right.end,
            right.distance,
            right.query_start,
            right.query_end,
            right.fallback_score,
        )))
}

const fn strand_code(strand: bsbit_core::bisulfite::BisulfiteStrand) -> u64 {
    use bsbit_core::bisulfite::BisulfiteStrand;
    match strand {
        BisulfiteStrand::OT => 0,
        BisulfiteStrand::OB => 1,
        BisulfiteStrand::CTOT => 2,
        BisulfiteStrand::CTOB => 3,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bsbit_core::alphabet::Base;
    use bsbit_core::bisulfite::BisulfiteStrand;
    use bsbit_core::sequence::NormalizedSequence;
    use bsbit_index::reference::{ContigInput, ReferenceBuildLimits};

    fn reference(names: [&[u8]; 2]) -> ReferenceIndex {
        ReferenceIndex::build(
            names
                .into_iter()
                .map(|name| {
                    ContigInput::new(
                        name.to_vec(),
                        NormalizedSequence::from_bases([Base::A; 128]),
                    )
                })
                .collect(),
            ReferenceBuildLimits::MAX,
        )
        .expect("two-contig tie-break reference builds")
    }

    #[test]
    fn stable_mixer_is_seed_and_read_dependent() {
        let base = absorb(absorb(DOMAIN, 0), 17);
        assert_eq!(base, absorb(absorb(DOMAIN, 0), 17));
        assert_ne!(base, absorb(absorb(DOMAIN, 1), 17));
        assert_ne!(base, absorb(absorb(DOMAIN, 0), 18));
    }

    #[test]
    fn byte_hash_has_explicit_chunk_boundaries() {
        let state = absorb(DOMAIN, 9);
        assert_ne!(absorb_bytes(state, b"a"), absorb_bytes(state, b"a\0"));
        assert_ne!(
            absorb_bytes(state, b"abcdefgh"),
            absorb_bytes(state, b"abcdefgh\0")
        );
    }

    #[test]
    fn origin_hash_uses_contig_name_instead_of_reference_order() {
        let forward = reference([b"alpha", b"beta"]);
        let reversed = reference([b"beta", b"alpha"]);
        let tie_break = ReportingTieBreak {
            seed: 42,
            read_key: 7,
        };
        let forward_alpha = ReadPlacement::strict(0, 11, 31, BisulfiteStrand::OT, 0);
        let reversed_alpha = ReadPlacement::strict(1, 11, 31, BisulfiteStrand::OT, 0);
        assert_eq!(
            placement_origin_hash(&forward, tie_break, forward_alpha, 20),
            placement_origin_hash(&reversed, tie_break, reversed_alpha, 20),
        );
        assert_ne!(
            placement_origin_hash(&forward, tie_break, forward_alpha, 20),
            placement_origin_hash(
                &forward,
                ReportingTieBreak {
                    seed: 43,
                    ..tie_break
                },
                forward_alpha,
                20,
            ),
        );
    }
}
