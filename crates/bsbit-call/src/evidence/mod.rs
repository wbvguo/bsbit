//! Portable word-parallel base classification for calling hot paths.

pub(crate) mod fragment;

/// The canonical two-bit base code used by bsbit calling kernels.
///
/// The low and high bits follow `A=00`, `C=01`, `G=10`, and `T=11`.
/// Ambiguous or absent bases are represented by a cleared validity bit rather
/// than by aliasing them with `A`.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub(crate) enum BaseCode {
    A,
    C,
    G,
    T,
}

impl BaseCode {
    pub(crate) const ALL: [Self; 4] = [Self::A, Self::C, Self::G, Self::T];

    pub(crate) const fn from_ascii(base: u8) -> Option<Self> {
        match base {
            b'A' | b'a' => Some(Self::A),
            b'C' | b'c' => Some(Self::C),
            b'G' | b'g' => Some(Self::G),
            b'T' | b't' => Some(Self::T),
            _ => None,
        }
    }

    pub(crate) const fn ascii(self) -> u8 {
        match self {
            Self::A => b'A',
            Self::C => b'C',
            Self::G => b'G',
            Self::T => b'T',
        }
    }

    pub(crate) const fn index(self) -> usize {
        match self {
            Self::A => 0,
            Self::C => 1,
            Self::G => 2,
            Self::T => 3,
        }
    }

    const fn low(self) -> bool {
        matches!(self, Self::C | Self::T)
    }

    const fn high(self) -> bool {
        matches!(self, Self::G | Self::T)
    }
}

/// Two bit-planes and a canonical-base validity plane for at most 64 bases.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub(crate) struct BasePlanes {
    low: u64,
    high: u64,
    valid: u64,
}

impl BasePlanes {
    /// Inserts one base at `offset`; `None` leaves the position invalid.
    pub(crate) fn insert(&mut self, offset: usize, base: Option<BaseCode>) {
        debug_assert!(offset < u64::BITS as usize);
        let Some(base) = base else {
            return;
        };
        let bit = 1_u64 << offset;
        self.valid |= bit;
        if base.low() {
            self.low |= bit;
        }
        if base.high() {
            self.high |= bit;
        }
    }

    /// Returns the mask of positions equal to `base`.
    pub(crate) const fn mask(self, base: BaseCode) -> u64 {
        match base {
            BaseCode::A => self.valid & !self.low & !self.high,
            BaseCode::C => self.valid & self.low & !self.high,
            BaseCode::G => self.valid & !self.low & self.high,
            BaseCode::T => self.valid & self.low & self.high,
        }
    }

    #[cfg(test)]
    pub(crate) const fn valid(self) -> u64 {
        self.valid
    }
}

/// Molecular cytosine strand supplying bisulfite evidence.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub(crate) enum EvidenceStrand {
    Top,
    Bottom,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum ContextClass {
    Cg,
    Chg,
    Chh,
}

impl ContextClass {
    pub(crate) const fn name(self) -> &'static [u8] {
        match self {
            Self::Cg => b"CG",
            Self::Chg => b"CHG",
            Self::Chh => b"CHH",
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct CytosineContext {
    pub(crate) class: ContextClass,
    pub(crate) second: u8,
}

/// One caller-neutral aligned observation after reference reconstruction.
/// Methylation, SNP, and joint modules consume this same representation.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct EvidenceObservation {
    pub(crate) reference: u32,
    pub(crate) position: u32,
    pub(crate) reference_base: u8,
    pub(crate) query_base: Option<u8>,
    pub(crate) base_quality: Option<u8>,
    pub(crate) mapping_quality: u8,
    pub(crate) strand: EvidenceStrand,
    pub(crate) context: Option<CytosineContext>,
}

impl EvidenceObservation {
    pub(crate) const fn key(self) -> (u32, u32) {
        (self.reference, self.position)
    }
}

/// Four disjoint masks produced by one directional methylation classification.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub(crate) struct MethylationMasks {
    pub(crate) methylated: u64,
    pub(crate) unmethylated: u64,
    pub(crate) deleted: u64,
    pub(crate) different: u64,
}

impl MethylationMasks {
    pub(crate) const fn callable(self) -> u64 {
        self.methylated | self.unmethylated | self.deleted | self.different
    }
}

/// Eight one-bit planes storing 64 independent saturating `u8` counters.
///
/// `increment_mask` performs carry-propagating addition for every selected
/// lane at once. Lanes that would wrap from 255 to zero are reported so the
/// caller can promote only those rare high-coverage positions to a sparse
/// wide-counter side table.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub(crate) struct BitSlicedU8 {
    planes: [u64; 8],
    wide: u64,
}

/// Lanes requiring sparse wide-counter handling after one masked increment.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub(crate) struct WideIncrement {
    /// Lanes that were already promoted before this increment.
    pub(crate) already_wide: u64,
    /// Lanes promoted by this increment, whose logical value is now 256.
    pub(crate) newly_wide: u64,
}

impl BitSlicedU8 {
    /// Adds one to every selected lane without scalar per-lane arithmetic.
    pub(crate) fn increment_mask(&mut self, mask: u64) -> WideIncrement {
        let already_wide = mask & self.wide;
        let mut carry = mask & !self.wide;
        for plane in &mut self.planes {
            let next_carry = *plane & carry;
            *plane ^= carry;
            carry = next_carry;
        }
        let newly_wide = carry;
        self.wide |= newly_wide;
        WideIncrement {
            already_wide,
            newly_wide,
        }
    }

    /// Returns a narrow lane value, or `None` after sparse promotion.
    pub(crate) fn narrow_value(self, offset: usize) -> Option<u8> {
        let bit = 1_u64.checked_shl(u32::try_from(offset).ok()?)?;
        if self.wide & bit != 0 {
            return None;
        }
        let mut value = 0_u8;
        for (plane_index, plane) in self.planes.iter().enumerate() {
            if plane & bit != 0 {
                value |= 1_u8 << plane_index;
            }
        }
        Some(value)
    }

    /// Returns the mask of lanes promoted beyond the inline `u8` range.
    #[cfg(test)]
    pub(crate) const fn wide_mask(self) -> u64 {
        self.wide
    }
}

/// Classifies at most 64 aligned reference/query columns with word operations.
///
/// `deletion` must be disjoint from `observed.valid()`. Positions outside the
/// logical block must be clear in every input plane.
pub(crate) const fn classify_methylation(
    reference: BasePlanes,
    observed: BasePlanes,
    observed_present: u64,
    deletion: u64,
    strand: EvidenceStrand,
) -> MethylationMasks {
    let (eligible, methylated_base, unmethylated_base) = match strand {
        EvidenceStrand::Top => (
            reference.mask(BaseCode::C),
            observed.mask(BaseCode::C),
            observed.mask(BaseCode::T),
        ),
        EvidenceStrand::Bottom => (
            reference.mask(BaseCode::G),
            observed.mask(BaseCode::G),
            observed.mask(BaseCode::A),
        ),
    };
    let methylated = eligible & methylated_base;
    let unmethylated = eligible & unmethylated_base;
    let deleted = eligible & deletion;
    let different = eligible & observed_present & !(methylated_base | unmethylated_base);
    MethylationMasks {
        methylated,
        unmethylated,
        deleted,
        different,
    }
}

#[cfg(test)]
mod tests {
    use super::{
        BaseCode, BasePlanes, BitSlicedU8, EvidenceStrand, MethylationMasks, WideIncrement,
        classify_methylation,
    };

    #[test]
    fn two_bit_planes_reconstruct_every_position_across_one_word() {
        let mut planes = BasePlanes::default();
        for offset in 0..64 {
            planes.insert(offset, Some(BaseCode::ALL[offset % BaseCode::ALL.len()]));
        }
        for base in BaseCode::ALL {
            let expected = (0..64).fold(0_u64, |mask, offset| {
                mask | (u64::from(BaseCode::ALL[offset % 4] == base) << offset)
            });
            assert_eq!(planes.mask(base), expected);
        }
    }

    #[test]
    fn invalid_bases_do_not_alias_canonical_a() {
        let mut planes = BasePlanes::default();
        planes.insert(0, Some(BaseCode::A));
        planes.insert(1, None);
        assert_eq!(planes.mask(BaseCode::A), 1);
        assert_eq!(planes.valid(), 1);
    }

    #[test]
    fn word_classifier_matches_directional_truth_table() {
        let mut reference = BasePlanes::default();
        let mut observed = BasePlanes::default();
        let reference_bases = [
            BaseCode::C,
            BaseCode::C,
            BaseCode::C,
            BaseCode::C,
            BaseCode::G,
            BaseCode::G,
            BaseCode::G,
            BaseCode::G,
        ];
        let observed_bases = [
            Some(BaseCode::C),
            Some(BaseCode::T),
            Some(BaseCode::A),
            None,
            Some(BaseCode::G),
            Some(BaseCode::A),
            Some(BaseCode::T),
            None,
        ];
        for (offset, base) in reference_bases.into_iter().enumerate() {
            reference.insert(offset, Some(base));
            observed.insert(offset, observed_bases[offset]);
        }
        let deletion = (1 << 3) | (1 << 7);
        assert_eq!(
            classify_methylation(
                reference,
                observed,
                observed.valid(),
                deletion,
                EvidenceStrand::Top,
            ),
            MethylationMasks {
                methylated: 1 << 0,
                unmethylated: 1 << 1,
                deleted: 1 << 3,
                different: 1 << 2,
            }
        );
        assert_eq!(
            classify_methylation(
                reference,
                observed,
                observed.valid(),
                deletion,
                EvidenceStrand::Bottom,
            ),
            MethylationMasks {
                methylated: 1 << 4,
                unmethylated: 1 << 5,
                deleted: 1 << 7,
                different: 1 << 6,
            }
        );
    }

    #[test]
    fn present_ambiguous_query_base_is_different_not_a() {
        let mut reference = BasePlanes::default();
        reference.insert(0, Some(BaseCode::G));
        let observed = BasePlanes::default();
        assert_eq!(
            classify_methylation(reference, observed, 1, 0, EvidenceStrand::Bottom,),
            MethylationMasks {
                different: 1,
                ..MethylationMasks::default()
            }
        );
    }

    #[test]
    fn bit_sliced_counter_matches_scalar_lanes_and_reports_promotion_once() {
        let mut counters = BitSlicedU8::default();
        let mut scalar = [0_u16; 64];
        for round in 0..300_u16 {
            let mask = if round.is_multiple_of(3) {
                u64::MAX
            } else if round.is_multiple_of(2) {
                0x5555_5555_5555_5555
            } else {
                0x8000_0000_0000_0001
            };
            let increment = counters.increment_mask(mask);
            for (offset, value) in scalar.iter_mut().enumerate() {
                if mask & (1_u64 << offset) != 0 {
                    let was_wide = *value > u16::from(u8::MAX);
                    *value += 1;
                    let is_newly_wide = *value == u16::from(u8::MAX) + 1;
                    assert_eq!(increment.already_wide & (1_u64 << offset) != 0, was_wide);
                    assert_eq!(increment.newly_wide & (1_u64 << offset) != 0, is_newly_wide);
                }
            }
            for (offset, value) in scalar.iter().copied().enumerate() {
                if let Ok(narrow) = u8::try_from(value) {
                    assert_eq!(counters.narrow_value(offset), Some(narrow));
                } else {
                    assert_eq!(counters.narrow_value(offset), None);
                }
            }
        }
        assert_ne!(counters.wide_mask(), 0);
    }

    #[test]
    fn bit_sliced_counter_increments_disjoint_masks_in_parallel() {
        let mut counters = BitSlicedU8::default();
        assert_eq!(
            counters.increment_mask(0xaaaa_aaaa_aaaa_aaaa),
            WideIncrement::default()
        );
        assert_eq!(
            counters.increment_mask(0xffff_0000_ffff_0000),
            WideIncrement::default()
        );
        assert_eq!(counters.narrow_value(0), Some(0));
        assert_eq!(counters.narrow_value(1), Some(1));
        assert_eq!(counters.narrow_value(16), Some(1));
        assert_eq!(counters.narrow_value(17), Some(2));
        assert_eq!(counters.narrow_value(64), None);
    }
}
