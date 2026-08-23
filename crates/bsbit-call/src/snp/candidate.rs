//! Dense first-pass SNP candidate detection.

use std::collections::HashMap;

use super::result::{Base, SnpConfig, VariantCall, filtered_observation, validate_config};
use crate::CallError;
use crate::evidence::{BitSlicedU8, EvidenceObservation, EvidenceStrand};

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
enum CandidateCounter {
    Reference,
    Alternate(Base),
}

impl CandidateCounter {
    const ALL: [Self; 5] = [
        Self::Reference,
        Self::Alternate(Base::A),
        Self::Alternate(Base::C),
        Self::Alternate(Base::G),
        Self::Alternate(Base::T),
    ];

    const fn index(self) -> usize {
        match self {
            Self::Reference => 0,
            Self::Alternate(base) => base.index() + 1,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
struct CandidateBlock {
    reference: BitSlicedU8,
    alternate: [BitSlicedU8; 4],
}

pub(crate) const fn candidate_dense_bytes_per_block() -> usize {
    std::mem::size_of::<CandidateBlock>() + u64::BITS as usize
}

/// Conservative per-64-base regional footprint for SNP first-pass state and
/// worst-case retained candidates/calls waiting for ordered output.
pub(crate) const fn snp_region_bytes_per_block() -> usize {
    candidate_dense_bytes_per_block()
        + u64::BITS as usize
            * (std::mem::size_of::<CandidateSite>() + std::mem::size_of::<(u32, VariantCall)>())
}

impl CandidateBlock {
    fn counter_mut(&mut self, counter: CandidateCounter) -> &mut BitSlicedU8 {
        match counter {
            CandidateCounter::Reference => &mut self.reference,
            CandidateCounter::Alternate(base) => &mut self.alternate[base.index()],
        }
    }

    const fn counter(&self, counter: CandidateCounter) -> BitSlicedU8 {
        match counter {
            CandidateCounter::Reference => self.reference,
            CandidateCounter::Alternate(base) => self.alternate[base.index()],
        }
    }
}

/// One first-pass SNP candidate retained for exact quality likelihoods.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct CandidateSite {
    pub(crate) position: u32,
    pub(super) reference: Base,
}

#[cfg(test)]
impl CandidateSite {
    pub(crate) fn for_test(position: u32, reference: u8) -> Self {
        Self {
            position,
            reference: Base::from_ascii(reference).expect("test reference is canonical"),
        }
    }
}

/// Dense, region-bounded first pass using bit-sliced counters.
pub(crate) struct CandidateRegion {
    start: u32,
    end: u32,
    blocks: Vec<CandidateBlock>,
    reference: Vec<u8>,
    wide: HashMap<(u32, CandidateCounter), u64>,
    config: SnpConfig,
}

impl CandidateRegion {
    pub(crate) fn new(start: u32, end: u32, config: SnpConfig) -> Result<Self, CallError> {
        validate_config(config)?;
        if start >= end {
            return Err(CallError::operation(
                "SNP candidate region must be nonempty",
            ));
        }
        let length = usize::try_from(end - start)
            .map_err(|_| CallError::operation("SNP candidate region is not addressable"))?;
        let block_count = length
            .checked_add(63)
            .ok_or_else(|| CallError::operation("SNP candidate block count overflowed"))?
            / 64;
        let mut blocks = Vec::new();
        blocks.try_reserve_exact(block_count).map_err(|error| {
            CallError::with_source(
                crate::CallErrorKind::Calling,
                format!("allocate {block_count} SNP candidate blocks"),
                error,
            )
        })?;
        blocks.resize(block_count, CandidateBlock::default());
        let mut reference = Vec::new();
        reference.try_reserve_exact(length).map_err(|error| {
            CallError::with_source(
                crate::CallErrorKind::Calling,
                format!("allocate {length} SNP reference bytes"),
                error,
            )
        })?;
        reference.resize(length, 0);
        Ok(Self {
            start,
            end,
            blocks,
            reference,
            wide: HashMap::new(),
            config,
        })
    }

    pub(crate) fn observe_fragment(
        &mut self,
        observations: &[EvidenceObservation],
    ) -> Result<(), CallError> {
        let mut current_block = None;
        let mut count_masks = [0_u64; 5];
        for observation in observations {
            let Some((reference, observed, _, _)) = filtered_observation(*observation, self.config)
            else {
                continue;
            };
            let offset = self.offset(observation.position)?;
            let encoded_reference =
                u8::try_from(reference.index() + 1).expect("four canonical bases fit one byte");
            match self.reference[offset] {
                0 => self.reference[offset] = encoded_reference,
                value if value == encoded_reference => {}
                _ => {
                    return Err(CallError::operation(format!(
                        "inconsistent reconstructed reference at SNP position {}",
                        observation.position
                    )));
                }
            }
            let counter = if observed == reference {
                CandidateCounter::Reference
            } else if is_conversion_confounded(reference, observed, observation.strand) {
                continue;
            } else {
                CandidateCounter::Alternate(observed)
            };
            let block_index = offset / 64;
            if current_block.is_some_and(|current| current != block_index) {
                self.increment_block_masks(
                    current_block.expect("checked present block"),
                    count_masks,
                )?;
                count_masks = [0; 5];
            }
            current_block = Some(block_index);
            count_masks[counter.index()] |= 1_u64 << (offset % 64);
        }
        if let Some(block_index) = current_block {
            self.increment_block_masks(block_index, count_masks)?;
        }
        Ok(())
    }

    pub(crate) fn candidates(&self) -> Result<Vec<CandidateSite>, CallError> {
        let mut candidates = Vec::new();
        for (offset, encoded) in self.reference.iter().copied().enumerate() {
            if encoded == 0 {
                continue;
            }
            let reference = Base::ALL
                .get(usize::from(encoded - 1))
                .copied()
                .ok_or_else(|| CallError::operation("invalid encoded SNP reference base"))?;
            let position = self
                .start
                .checked_add(
                    u32::try_from(offset)
                        .map_err(|_| CallError::operation("SNP candidate offset exceeds u32"))?,
                )
                .ok_or_else(|| CallError::operation("SNP candidate position overflowed"))?;
            let reference_count = self.count(position, offset, CandidateCounter::Reference)?;
            let mut alternate_total = 0_u64;
            let mut maximum_alternate = 0_u64;
            for base in Base::ALL {
                if base == reference {
                    continue;
                }
                let count = self.count(position, offset, CandidateCounter::Alternate(base))?;
                alternate_total = alternate_total
                    .checked_add(count)
                    .ok_or_else(|| CallError::operation("SNP alternate count overflowed u64"))?;
                maximum_alternate = maximum_alternate.max(count);
            }
            let depth = reference_count
                .checked_add(alternate_total)
                .ok_or_else(|| CallError::operation("SNP candidate depth overflowed u64"))?;
            if depth >= u64::from(self.config.minimum_depth)
                && maximum_alternate >= u64::from(self.config.minimum_alternate_count)
                && u128::from(maximum_alternate) * 1_000_000_000_u128
                    >= u128::from(depth)
                        * u128::from(self.config.minimum_alternate_fraction_parts_per_billion)
            {
                candidates.try_reserve(1).map_err(|error| {
                    CallError::with_source(
                        crate::CallErrorKind::Calling,
                        "reserve SNP candidate result",
                        error,
                    )
                })?;
                candidates.push(CandidateSite {
                    position,
                    reference,
                });
            }
        }
        Ok(candidates)
    }

    fn offset(&self, position: u32) -> Result<usize, CallError> {
        if position < self.start || position >= self.end {
            return Err(CallError::operation(format!(
                "SNP observation {position} is outside region {}-{}",
                self.start, self.end
            )));
        }
        usize::try_from(position - self.start)
            .map_err(|_| CallError::operation("SNP observation offset is not addressable"))
    }

    fn increment_block_masks(
        &mut self,
        block_index: usize,
        masks: [u64; 5],
    ) -> Result<(), CallError> {
        for counter in CandidateCounter::ALL {
            let mask = masks[counter.index()];
            if mask == 0 {
                continue;
            }
            let increment = self
                .blocks
                .get_mut(block_index)
                .ok_or_else(|| CallError::operation("SNP candidate block is missing"))?
                .counter_mut(counter)
                .increment_mask(mask);
            let mut newly_wide = increment.newly_wide;
            while newly_wide != 0 {
                let lane = newly_wide.trailing_zeros() as usize;
                newly_wide &= newly_wide - 1;
                let position = self.position_for_block_lane(block_index, lane)?;
                if self.wide.insert((position, counter), 256).is_some() {
                    return Err(CallError::operation(
                        "SNP candidate counter was promoted twice",
                    ));
                }
            }
            let mut already_wide = increment.already_wide;
            while already_wide != 0 {
                let lane = already_wide.trailing_zeros() as usize;
                already_wide &= already_wide - 1;
                let position = self.position_for_block_lane(block_index, lane)?;
                let value = self
                    .wide
                    .get_mut(&(position, counter))
                    .ok_or_else(|| CallError::operation("SNP wide candidate counter is missing"))?;
                *value = value
                    .checked_add(1)
                    .ok_or_else(|| CallError::operation("SNP candidate count overflowed u64"))?;
            }
        }
        Ok(())
    }

    fn position_for_block_lane(&self, block_index: usize, lane: usize) -> Result<u32, CallError> {
        let offset = block_index
            .checked_mul(64)
            .and_then(|base| base.checked_add(lane))
            .ok_or_else(|| CallError::operation("SNP candidate block position overflowed"))?;
        self.start
            .checked_add(
                u32::try_from(offset).map_err(|_| {
                    CallError::operation("SNP candidate block position exceeds u32")
                })?,
            )
            .ok_or_else(|| CallError::operation("SNP candidate block position overflowed u32"))
    }

    fn count(
        &self,
        position: u32,
        offset: usize,
        counter: CandidateCounter,
    ) -> Result<u64, CallError> {
        let block = self
            .blocks
            .get(offset / 64)
            .ok_or_else(|| CallError::operation("SNP candidate block is missing"))?;
        if let Some(value) = block.counter(counter).narrow_value(offset % 64) {
            Ok(u64::from(value))
        } else {
            self.wide
                .get(&(position, counter))
                .copied()
                .ok_or_else(|| CallError::operation("SNP wide candidate counter is missing"))
        }
    }
}

const fn is_conversion_confounded(reference: Base, observed: Base, strand: EvidenceStrand) -> bool {
    matches!(
        (reference, observed, strand),
        (Base::C, Base::T, EvidenceStrand::Top) | (Base::G, Base::A, EvidenceStrand::Bottom)
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn observation(
        reference: u8,
        observed: u8,
        strand: EvidenceStrand,
        quality: u8,
    ) -> EvidenceObservation {
        EvidenceObservation {
            reference: 0,
            position: 10,
            reference_base: reference,
            query_base: Some(observed),
            base_quality: Some(quality),
            mapping_quality: 60,
            strand,
            context: None,
        }
    }
    #[test]
    fn candidate_filter_ignores_only_conversion_confounded_strand() {
        let config = SnpConfig {
            minimum_depth: 2,
            minimum_alternate_count: 2,
            ..SnpConfig::default()
        };
        let mut top = CandidateRegion::new(10, 11, config).unwrap();
        for _ in 0..2 {
            top.observe_fragment(&[observation(b'C', b'T', EvidenceStrand::Top, 40)])
                .unwrap();
        }
        assert!(top.candidates().unwrap().is_empty());

        let mut bottom = CandidateRegion::new(10, 11, config).unwrap();
        for _ in 0..2 {
            bottom
                .observe_fragment(&[observation(b'C', b'T', EvidenceStrand::Bottom, 40)])
                .unwrap();
        }
        assert_eq!(bottom.candidates().unwrap().len(), 1);
    }

    #[test]
    fn candidate_fraction_uses_the_strongest_alternate() {
        let strict = SnpConfig {
            minimum_depth: 1,
            minimum_alternate_count: 1,
            minimum_alternate_fraction_parts_per_billion: 100_000_000,
            ..SnpConfig::default()
        };
        let mut region = CandidateRegion::new(10, 11, strict).unwrap();
        for _ in 0..9 {
            region
                .observe_fragment(&[observation(b'A', b'A', EvidenceStrand::Top, 40)])
                .unwrap();
        }
        region
            .observe_fragment(&[observation(b'A', b'C', EvidenceStrand::Top, 40)])
            .unwrap();
        region
            .observe_fragment(&[observation(b'A', b'G', EvidenceStrand::Top, 40)])
            .unwrap();
        assert!(region.candidates().unwrap().is_empty());

        let sensitive = SnpConfig {
            minimum_alternate_fraction_parts_per_billion: 50_000_000,
            ..strict
        };
        let mut region = CandidateRegion::new(10, 11, sensitive).unwrap();
        for _ in 0..9 {
            region
                .observe_fragment(&[observation(b'A', b'A', EvidenceStrand::Top, 40)])
                .unwrap();
        }
        region
            .observe_fragment(&[observation(b'A', b'C', EvidenceStrand::Top, 40)])
            .unwrap();
        region
            .observe_fragment(&[observation(b'A', b'G', EvidenceStrand::Top, 40)])
            .unwrap();
        assert_eq!(region.candidates().unwrap().len(), 1);
    }

    #[test]
    fn candidate_bit_sliced_count_promotes_above_255() {
        let config = SnpConfig {
            minimum_depth: 4,
            minimum_alternate_count: 2,
            ..SnpConfig::default()
        };
        let mut region = CandidateRegion::new(10, 11, config).unwrap();
        for _ in 0..300 {
            region
                .observe_fragment(&[observation(b'A', b'G', EvidenceStrand::Top, 40)])
                .unwrap();
        }
        assert_eq!(region.candidates().unwrap().len(), 1);
    }
}
