//! Verified per-read placement facts shared by single-end and paired-end mapping.

#[cfg(test)]
use bsbit_core::alphabet::Base;
#[cfg(test)]
use bsbit_core::bisulfite::CytosineStrand;
use bsbit_core::bisulfite::{AlignmentOrientation, BisulfiteStrand, strand_semantics};
#[cfg(test)]
use bsbit_index::reference::ReferenceIndex;

use crate::alignment_policy::SEMI_GLOBAL_EDIT_PENALTY;

pub(crate) const FULL_QUERY_END: u16 = u16::MAX;

/// One verified in-budget placement represented in reference coordinates.
#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub struct ReadPlacement {
    pub(crate) contig_ordinal: u64,
    pub(crate) start: u64,
    pub(crate) end: u64,
    pub(crate) strand: BisulfiteStrand,
    pub(crate) distance: u8,
    pub(crate) query_start: u16,
    pub(crate) query_end: u16,
    pub(crate) fallback_score: u8,
}

impl ReadPlacement {
    pub(crate) const fn strict(
        contig_ordinal: u64,
        start: u64,
        end: u64,
        strand: BisulfiteStrand,
        distance: u8,
    ) -> Self {
        Self {
            contig_ordinal,
            start,
            end,
            strand,
            distance,
            query_start: 0,
            query_end: FULL_QUERY_END,
            fallback_score: distance.saturating_mul(SEMI_GLOBAL_EDIT_PENALTY),
        }
    }

    /// Returns the zero-based contig ordinal.
    #[must_use]
    pub const fn contig_ordinal(self) -> u64 {
        self.contig_ordinal
    }

    /// Returns the zero-based reference start.
    #[must_use]
    pub const fn start(self) -> u64 {
        self.start
    }

    /// Returns the exclusive reference end.
    #[must_use]
    pub const fn end(self) -> u64 {
        self.end
    }

    /// Returns the bisulfite alignment strand.
    #[must_use]
    pub const fn strand(self) -> BisulfiteStrand {
        self.strand
    }

    /// Returns the conversion-aware edit distance.
    #[must_use]
    pub const fn distance(self) -> u8 {
        self.distance
    }

    /// Returns the retained sequencing-orientation query interval.
    #[must_use]
    pub fn retained_query_interval(self, read_length: usize) -> core::ops::Range<usize> {
        let end = if self.query_end == FULL_QUERY_END {
            read_length
        } else {
            usize::from(self.query_end)
        };
        usize::from(self.query_start)..end
    }

    /// Reports whether this placement retained less than the complete read.
    #[must_use]
    pub fn is_soft_clipped(self, read_length: usize) -> bool {
        let retained = self.retained_query_interval(read_length);
        retained.start != 0 || retained.end != read_length
    }
}

pub(crate) fn placement_net_gap_bases(placement: ReadPlacement, read_len: usize) -> u64 {
    let reference_bases = placement.end().saturating_sub(placement.start());
    let retained_query = placement.retained_query_interval(read_len);
    let query_bases = u64::try_from(retained_query.end.saturating_sub(retained_query.start))
        .expect("bounded query span fits u64");
    reference_bases.abs_diff(query_bases)
}

pub(crate) fn placement_origin_key(
    placement: ReadPlacement,
    read_length: usize,
) -> (u64, BisulfiteStrand, i128) {
    let sequencing_five_prime_clip =
        i128::try_from(placement.retained_query_interval(read_length).start)
            .expect("bounded read length fits i128");
    let five_prime = match strand_semantics(placement.strand()).orientation() {
        AlignmentOrientation::Forward => i128::from(placement.start()) - sequencing_five_prime_clip,
        AlignmentOrientation::Reverse => {
            i128::from(placement.end()) - 1 + sequencing_five_prime_clip
        }
    };
    (placement.contig_ordinal(), placement.strand(), five_prime)
}

/// Counts converted and unconverted susceptible bases by CG/CHG/CHH context
/// for a full-length ungapped placement.
#[cfg(test)]
pub(crate) fn placement_conversion_counts(
    reference: &ReferenceIndex,
    read: &[Base],
    placement: ReadPlacement,
) -> Option<([u16; 3], [u16; 3])> {
    let retained = placement.retained_query_interval(read.len());
    let reference_span = placement.end().checked_sub(placement.start())?;
    if retained.start != 0
        || retained.end != read.len()
        || reference_span != u64::try_from(read.len()).ok()?
    {
        return None;
    }
    let contig = reference.contig_by_ordinal(placement.contig_ordinal())?;
    let start = usize::try_from(placement.start()).ok()?;
    let end = usize::try_from(placement.end()).ok()?;
    let reference_bases = contig.sequence().bases();
    let aligned = reference_bases.get(start..end)?;
    let semantics = strand_semantics(placement.strand());
    let mut converted = [0_u16; 3];
    let mut unconverted = [0_u16; 3];
    for (offset, &reference_base) in aligned.iter().enumerate() {
        let query_base = match semantics.orientation() {
            AlignmentOrientation::Forward => read[offset],
            AlignmentOrientation::Reverse => read[read.len() - offset - 1].complement(),
        };
        let absolute = start + offset;
        let (susceptible, converted_base) = match semantics.cytosine_strand() {
            CytosineStrand::Top => (Base::C, Base::T),
            CytosineStrand::Bottom => (Base::G, Base::A),
        };
        if reference_base != susceptible {
            continue;
        }
        let context = cytosine_context(reference_bases, absolute, semantics.cytosine_strand());
        if query_base == converted_base {
            converted[context] = converted[context].saturating_add(1);
        } else if query_base == susceptible {
            unconverted[context] = unconverted[context].saturating_add(1);
        }
    }
    Some((converted, unconverted))
}

#[cfg(test)]
fn cytosine_context(reference: &[Base], position: usize, strand: CytosineStrand) -> usize {
    let (first_is_g, second_is_g) = match strand {
        CytosineStrand::Top => (
            position
                .checked_add(1)
                .and_then(|index| reference.get(index))
                == Some(&Base::G),
            position
                .checked_add(2)
                .and_then(|index| reference.get(index))
                == Some(&Base::G),
        ),
        CytosineStrand::Bottom => (
            position
                .checked_sub(1)
                .and_then(|index| reference.get(index))
                == Some(&Base::C),
            position
                .checked_sub(2)
                .and_then(|index| reference.get(index))
                == Some(&Base::C),
        ),
    };
    if first_is_g {
        0
    } else if second_is_g {
        1
    } else {
        2
    }
}
