//! Verified per-read placement facts shared by single-end and paired-end mapping.

use bsbit_core::bisulfite::{AlignmentOrientation, BisulfiteStrand, strand_semantics};

pub(crate) const FULL_QUERY_END: u16 = u16::MAX;
pub(crate) const SEMI_GLOBAL_EDIT_PENALTY: u8 = 7;

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
