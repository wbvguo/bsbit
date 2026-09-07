//! Reference-order-independent selection of one equal-best pair representative.

use bsbit_core::alphabet::Base;
use bsbit_index::reference::ReferenceIndex;

use super::result::PairedPlacement;
use crate::AlignmentError;
use crate::placement::placement_net_gap_bases;
use crate::reporting_tie_break::{
    ReportingTieBreak, compare_placements_without_reference_order, pair_origin_hash,
};

// Classification, score confidence, and MAPQ evidence must be merged under
// the same directional tie decision.
pub(super) fn prefer_fair_pair_representative(
    reference: &ReferenceIndex,
    reads: [&[Base]; 2],
    pairs: &mut [PairedPlacement],
    tie_break: ReportingTieBreak,
    reporting_order_swapped: bool,
) -> Result<(), AlignmentError> {
    let Some(first) = pairs.first().copied() else {
        return Ok(());
    };
    let minimum_net_gap = placement_net_gap_bases(first.mate1(), reads[0].len())
        .saturating_add(placement_net_gap_bases(first.mate2(), reads[1].len()));
    let mut best_index = 0_usize;
    let mut best_hash = pair_origin_hash_in_reporting_order(
        reference,
        tie_break,
        first,
        reads,
        reporting_order_swapped,
    )?;
    for (index, pair) in pairs.iter().copied().enumerate().skip(1) {
        let net_gap = placement_net_gap_bases(pair.mate1(), reads[0].len())
            .saturating_add(placement_net_gap_bases(pair.mate2(), reads[1].len()));
        if net_gap != minimum_net_gap {
            continue;
        }
        let hash = pair_origin_hash_in_reporting_order(
            reference,
            tie_break,
            pair,
            reads,
            reporting_order_swapped,
        )?;
        if hash < best_hash
            || (hash == best_hash
                && pair_less_in_reporting_order(
                    reference,
                    reads,
                    pair,
                    pairs[best_index],
                    reporting_order_swapped,
                )?)
        {
            best_index = index;
            best_hash = hash;
        }
    }
    pairs.swap(0, best_index);
    Ok(())
}

pub(super) fn pair_origin_hash_in_reporting_order(
    reference: &ReferenceIndex,
    tie_break: ReportingTieBreak,
    pair: PairedPlacement,
    reads: [&[Base]; 2],
    reporting_order_swapped: bool,
) -> Result<u64, AlignmentError> {
    if reporting_order_swapped {
        pair_origin_hash(
            reference,
            tie_break,
            pair.mate2(),
            reads[1].len(),
            pair.mate1(),
            reads[0].len(),
        )
    } else {
        pair_origin_hash(
            reference,
            tie_break,
            pair.mate1(),
            reads[0].len(),
            pair.mate2(),
            reads[1].len(),
        )
    }
}

fn pair_less_in_reporting_order(
    reference: &ReferenceIndex,
    reads: [&[Base]; 2],
    left: PairedPlacement,
    right: PairedPlacement,
    reporting_order_swapped: bool,
) -> Result<bool, AlignmentError> {
    if reporting_order_swapped {
        let first = compare_placements_without_reference_order(
            reference,
            left.mate2(),
            right.mate2(),
            reads[1].len(),
        )?;
        if !first.is_eq() {
            return Ok(first.is_lt());
        }
        Ok(compare_placements_without_reference_order(
            reference,
            left.mate1(),
            right.mate1(),
            reads[0].len(),
        )?
        .is_lt())
    } else {
        pair_less_without_reference_order(reference, reads, left, right)
    }
}

pub(super) fn pair_less_without_reference_order(
    reference: &ReferenceIndex,
    reads: [&[Base]; 2],
    left: PairedPlacement,
    right: PairedPlacement,
) -> Result<bool, AlignmentError> {
    let first = compare_placements_without_reference_order(
        reference,
        left.mate1(),
        right.mate1(),
        reads[0].len(),
    )?;
    if !first.is_eq() {
        return Ok(first.is_lt());
    }
    Ok(compare_placements_without_reference_order(
        reference,
        left.mate2(),
        right.mate2(),
        reads[1].len(),
    )?
    .is_lt())
}
