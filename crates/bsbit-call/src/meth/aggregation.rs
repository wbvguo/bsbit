//! Region-local bit-sliced methylation state.

use std::collections::HashMap;

use super::Parameters as MethParameters;
use crate::evidence::fragment::merge_contexts;
use crate::evidence::{
    BaseCode, BasePlanes, BitSlicedU8, ContextClass, CytosineContext, EvidenceObservation,
    EvidenceStrand, MethylationMasks, classify_methylation,
};
use crate::{CallError, CallErrorKind};

#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub(crate) struct SiteKey {
    pub(crate) reference: u32,
    pub(crate) position: u32,
    pub(crate) strand: EvidenceStrand,
}

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub(crate) enum CallKind {
    Methylated,
    Unmethylated,
    Deleted,
    Different,
}

impl CallKind {
    const ALL: [Self; 4] = [
        Self::Methylated,
        Self::Unmethylated,
        Self::Deleted,
        Self::Different,
    ];

    const fn index(self) -> usize {
        match self {
            Self::Methylated => 0,
            Self::Unmethylated => 1,
            Self::Deleted => 2,
            Self::Different => 3,
        }
    }
}

#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub(crate) struct SiteCounts {
    pub(crate) context: Option<CytosineContext>,
    pub(crate) methylated: u64,
    pub(crate) unmethylated: u64,
    pub(crate) deleted: u64,
    pub(crate) different: u64,
}

impl SiteCounts {
    pub(crate) fn valid_coverage(&self) -> Result<u64, CallError> {
        self.methylated
            .checked_add(self.unmethylated)
            .ok_or_else(|| CallError::operation("methylation coverage overflowed u64"))
    }
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
struct BitSlicedSiteBlock {
    occupied: u64,
    methylated: BitSlicedU8,
    unmethylated: BitSlicedU8,
    deleted: BitSlicedU8,
    different: BitSlicedU8,
}

impl BitSlicedSiteBlock {
    fn counter_mut(&mut self, call: CallKind) -> &mut BitSlicedU8 {
        match call {
            CallKind::Methylated => &mut self.methylated,
            CallKind::Unmethylated => &mut self.unmethylated,
            CallKind::Deleted => &mut self.deleted,
            CallKind::Different => &mut self.different,
        }
    }

    const fn counter(&self, call: CallKind) -> BitSlicedU8 {
        match call {
            CallKind::Methylated => self.methylated,
            CallKind::Unmethylated => self.unmethylated,
            CallKind::Deleted => self.deleted,
            CallKind::Different => self.different,
        }
    }
}

/// Exact dense-state footprint per 64-base block, including one metadata byte
/// per genomic position.
pub(crate) const fn meth_dense_bytes_per_block() -> usize {
    std::mem::size_of::<BitSlicedSiteBlock>() + u64::BITS as usize
}

pub(crate) fn accumulate_meth_fragment(
    observations: &[EvidenceObservation],
    sites: &mut DenseMethRegion,
    parameters: MethParameters,
) -> Result<(), CallError> {
    for block in observations.chunks(u64::BITS as usize) {
        let mut reference_planes = BasePlanes::default();
        let mut observed_planes = BasePlanes::default();
        let mut observed_present = 0_u64;
        let mut deletion = 0_u64;
        for (offset, observation) in block.iter().enumerate() {
            let mapping_quality_eligible = observation.mapping_quality != u8::MAX
                && observation.mapping_quality >= parameters.minimum_mapping_quality;
            let query_quality_eligible = match observation.query_base {
                Some(_) => observation
                    .base_quality
                    .is_some_and(|quality| quality >= parameters.minimum_base_quality),
                None => true,
            };
            if !mapping_quality_eligible || !query_quality_eligible {
                continue;
            }
            reference_planes.insert(offset, BaseCode::from_ascii(observation.reference_base));
            match observation.query_base {
                Some(base) => {
                    observed_present |= 1_u64 << offset;
                    observed_planes.insert(offset, BaseCode::from_ascii(base));
                }
                None => deletion |= 1_u64 << offset,
            }
        }
        let strand = block
            .first()
            .map(|observation| observation.strand)
            .expect("chunks are nonempty");
        if block.iter().any(|observation| observation.strand != strand) {
            return Err(CallError::input(
                "one fragment contains mixed bisulfite conversion strands",
            ));
        }
        let masks = classify_methylation(
            reference_planes,
            observed_planes,
            observed_present,
            deletion,
            strand,
        );
        sites.add_classified(block, masks)?;
    }
    Ok(())
}

/// Region-local dense state: four bit-sliced `u8` count categories, one
/// metadata byte per genomic position, occupancy/promotion masks, and sparse
/// `u64` storage only for lanes whose count exceeds 255.
#[derive(Debug)]
pub(crate) struct DenseMethRegion {
    reference: u32,
    start: u32,
    end: u32,
    blocks: Vec<BitSlicedSiteBlock>,
    metadata: Vec<u8>,
    wide: HashMap<(u32, CallKind), u64>,
}

impl DenseMethRegion {
    pub(crate) fn new(reference: u32, start: u32, end: u32) -> Result<Self, CallError> {
        if start >= end {
            return Err(CallError::operation("methylation region must be nonempty"));
        }
        let length = usize::try_from(end - start)
            .map_err(|_| CallError::operation("methylation region length is not addressable"))?;
        let block_count = length
            .checked_add(63)
            .ok_or_else(|| CallError::operation("methylation region block count overflowed"))?
            / 64;
        let mut blocks = Vec::new();
        blocks.try_reserve_exact(block_count).map_err(|error| {
            CallError::with_source(
                CallErrorKind::Calling,
                format!("allocate {block_count} methylation count blocks"),
                error,
            )
        })?;
        blocks.resize(block_count, BitSlicedSiteBlock::default());
        let mut metadata = Vec::new();
        metadata.try_reserve_exact(length).map_err(|error| {
            CallError::with_source(
                CallErrorKind::Calling,
                format!("allocate {length} methylation metadata bytes"),
                error,
            )
        })?;
        metadata.resize(length, 0);
        Ok(Self {
            reference,
            start,
            end,
            blocks,
            metadata,
            wide: HashMap::new(),
        })
    }

    fn offset(&self, key: SiteKey) -> Result<usize, CallError> {
        if key.reference != self.reference || key.position < self.start || key.position >= self.end
        {
            return Err(CallError::operation(format!(
                "site reference {}, position {} is outside region {}:{}-{}",
                key.reference, key.position, self.reference, self.start, self.end
            )));
        }
        usize::try_from(key.position - self.start)
            .map_err(|_| CallError::operation("methylation site offset is not addressable"))
    }

    #[cfg(test)]
    fn increment(&mut self, position: u32, offset: usize, call: CallKind) -> Result<(), CallError> {
        let block_index = offset / 64;
        let lane = offset % 64;
        let bit = 1_u64 << lane;
        let block = self
            .blocks
            .get_mut(block_index)
            .ok_or_else(|| CallError::operation("methylation counter block is missing"))?;
        block.occupied |= bit;
        let wide = block.counter_mut(call).increment_mask(bit);
        if wide.newly_wide != 0 {
            if self.wide.insert((position, call), 256).is_some() {
                return Err(CallError::operation(
                    "methylation wide counter was promoted twice",
                ));
            }
        } else if wide.already_wide != 0 {
            checked_increment(
                self.wide
                    .get_mut(&(position, call))
                    .ok_or_else(|| CallError::operation("methylation wide counter is missing"))?,
                "methylation wide count",
            )?;
        }
        Ok(())
    }

    fn increment_block_masks(
        &mut self,
        block_index: usize,
        masks: [u64; 4],
    ) -> Result<(), CallError> {
        for call in CallKind::ALL {
            let mask = masks[call.index()];
            if mask == 0 {
                continue;
            }
            let promoted = {
                let block = self
                    .blocks
                    .get_mut(block_index)
                    .ok_or_else(|| CallError::operation("methylation counter block is missing"))?;
                block.occupied |= mask;
                block.counter_mut(call).increment_mask(mask)
            };
            let mut newly_wide = promoted.newly_wide;
            while newly_wide != 0 {
                let lane = newly_wide.trailing_zeros() as usize;
                newly_wide &= newly_wide - 1;
                let position = self.position_for_block_lane(block_index, lane)?;
                if self.wide.insert((position, call), 256).is_some() {
                    return Err(CallError::operation(
                        "methylation wide counter was promoted twice",
                    ));
                }
            }
            let mut already_wide = promoted.already_wide;
            while already_wide != 0 {
                let lane = already_wide.trailing_zeros() as usize;
                already_wide &= already_wide - 1;
                let position = self.position_for_block_lane(block_index, lane)?;
                checked_increment(
                    self.wide.get_mut(&(position, call)).ok_or_else(|| {
                        CallError::operation("methylation wide counter is missing")
                    })?,
                    "methylation wide count",
                )?;
            }
        }
        Ok(())
    }

    fn position_for_block_lane(&self, block_index: usize, lane: usize) -> Result<u32, CallError> {
        let offset = block_index
            .checked_mul(64)
            .and_then(|base| base.checked_add(lane))
            .ok_or_else(|| CallError::operation("methylation block position overflowed"))?;
        self.start
            .checked_add(
                u32::try_from(offset)
                    .map_err(|_| CallError::operation("methylation block position exceeds u32"))?,
            )
            .ok_or_else(|| CallError::operation("methylation block position overflowed u32"))
    }

    fn count(&self, position: u32, offset: usize, call: CallKind) -> Result<u64, CallError> {
        let block = self
            .blocks
            .get(offset / 64)
            .ok_or_else(|| CallError::operation("methylation counter block is missing"))?;
        if let Some(value) = block.counter(call).narrow_value(offset % 64) {
            Ok(u64::from(value))
        } else {
            self.wide
                .get(&(position, call))
                .copied()
                .ok_or_else(|| CallError::operation("methylation wide counter is missing"))
        }
    }

    pub(super) fn add_classified(
        &mut self,
        observations: &[EvidenceObservation],
        masks: MethylationMasks,
    ) -> Result<(), CallError> {
        let mut current_block = None;
        let mut count_masks = [0_u64; 4];
        let mut callable = masks.callable();
        while callable != 0 {
            let evidence_offset = callable.trailing_zeros() as usize;
            callable &= callable - 1;
            let evidence = observations[evidence_offset];
            let key = SiteKey {
                reference: evidence.reference,
                position: evidence.position,
                strand: evidence.strand,
            };
            let call = call_kind_at(masks, evidence_offset);
            let region_offset = self.offset(key)?;
            let block_index = region_offset / 64;
            if current_block.is_some_and(|current| current != block_index) {
                self.increment_block_masks(
                    current_block.expect("checked present block"),
                    count_masks,
                )?;
                count_masks = [0; 4];
            }
            current_block = Some(block_index);
            self.metadata[region_offset] = merge_site_metadata(
                key,
                self.metadata[region_offset],
                key.strand,
                evidence.context,
            )?;
            count_masks[call.index()] |= 1_u64 << (region_offset % 64);
        }
        if let Some(block_index) = current_block {
            self.increment_block_masks(block_index, count_masks)?;
        }
        Ok(())
    }

    #[cfg(test)]
    pub(crate) fn add_observation(
        &mut self,
        key: SiteKey,
        context: Option<CytosineContext>,
        call: CallKind,
    ) -> Result<(), CallError> {
        let offset = self.offset(key)?;
        self.metadata[offset] =
            merge_site_metadata(key, self.metadata[offset], key.strand, context)?;
        self.increment(key.position, offset, call)
    }

    pub(crate) fn for_each_site(
        &self,
        mut consume: impl FnMut(SiteKey, SiteCounts) -> Result<(), CallError>,
    ) -> Result<(), CallError> {
        for (block_index, block) in self.blocks.iter().enumerate() {
            let mut occupied = block.occupied;
            while occupied != 0 {
                let lane = occupied.trailing_zeros() as usize;
                occupied &= occupied - 1;
                let offset = block_index
                    .checked_mul(64)
                    .and_then(|base| base.checked_add(lane))
                    .ok_or_else(|| CallError::operation("methylation site offset overflowed"))?;
                if offset >= self.metadata.len() {
                    continue;
                }
                let position = self
                    .start
                    .checked_add(
                        u32::try_from(offset).map_err(|_| {
                            CallError::operation("methylation site offset exceeds u32")
                        })?,
                    )
                    .ok_or_else(|| CallError::operation("methylation site position overflowed"))?;
                let (strand, context) = decode_site_metadata(self.metadata[offset])?;
                let key = SiteKey {
                    reference: self.reference,
                    position,
                    strand,
                };
                consume(
                    key,
                    SiteCounts {
                        context,
                        methylated: self.count(position, offset, CallKind::Methylated)?,
                        unmethylated: self.count(position, offset, CallKind::Unmethylated)?,
                        deleted: self.count(position, offset, CallKind::Deleted)?,
                        different: self.count(position, offset, CallKind::Different)?,
                    },
                )?;
            }
        }
        Ok(())
    }

    #[cfg(test)]
    pub(crate) fn into_sites(self) -> Result<Vec<(SiteKey, SiteCounts)>, CallError> {
        let mut sites = Vec::new();
        self.for_each_site(|key, counts| {
            sites.push((key, counts));
            Ok(())
        })?;
        Ok(sites)
    }
}

fn call_kind_at(masks: MethylationMasks, offset: usize) -> CallKind {
    let bit = 1_u64 << offset;
    if masks.methylated & bit != 0 {
        CallKind::Methylated
    } else if masks.unmethylated & bit != 0 {
        CallKind::Unmethylated
    } else if masks.deleted & bit != 0 {
        CallKind::Deleted
    } else {
        debug_assert_ne!(masks.different & bit, 0);
        CallKind::Different
    }
}

fn merge_site_metadata(
    key: SiteKey,
    encoded: u8,
    strand: EvidenceStrand,
    context: Option<CytosineContext>,
) -> Result<u8, CallError> {
    if encoded == 0 {
        return encode_site_metadata(strand, context);
    }
    let (existing_strand, existing_context) = decode_site_metadata(encoded)?;
    if existing_strand != strand {
        return Err(CallError::input(format!(
            "inconsistent methylation strand at reference {}, position {}",
            key.reference, key.position
        )));
    }
    encode_site_metadata(
        strand,
        merge_contexts(key.reference, key.position, existing_context, context)?,
    )
}

fn encode_site_metadata(
    strand: EvidenceStrand,
    context: Option<CytosineContext>,
) -> Result<u8, CallError> {
    const OCCUPIED: u8 = 0x80;
    const BOTTOM: u8 = 0x40;
    const RESOLVED: u8 = 0x20;
    let mut encoded = OCCUPIED;
    if strand == EvidenceStrand::Bottom {
        encoded |= BOTTOM;
    }
    let Some(context) = context else {
        return Ok(encoded);
    };
    let class = match context.class {
        ContextClass::Cg => 0_u8,
        ContextClass::Chg => 1,
        ContextClass::Chh => 2,
    };
    let second = match context.second {
        b'A' => 0_u8,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        value => {
            return Err(CallError::input(format!(
                "invalid methylation context base 0x{value:02x}"
            )));
        }
    };
    Ok(encoded | RESOLVED | (class << 2) | second)
}

fn decode_site_metadata(
    encoded: u8,
) -> Result<(EvidenceStrand, Option<CytosineContext>), CallError> {
    const OCCUPIED: u8 = 0x80;
    const BOTTOM: u8 = 0x40;
    const RESOLVED: u8 = 0x20;
    if encoded & OCCUPIED == 0 {
        return Err(CallError::operation(
            "methylation site metadata is unoccupied",
        ));
    }
    let strand = if encoded & BOTTOM == 0 {
        EvidenceStrand::Top
    } else {
        EvidenceStrand::Bottom
    };
    if encoded & RESOLVED == 0 {
        return Ok((strand, None));
    }
    let class = match (encoded >> 2) & 0x03 {
        0 => ContextClass::Cg,
        1 => ContextClass::Chg,
        2 => ContextClass::Chh,
        value => {
            return Err(CallError::operation(format!(
                "invalid encoded methylation context class {value}"
            )));
        }
    };
    let second = match encoded & 0x03 {
        0 => b'A',
        1 => b'C',
        2 => b'G',
        3 => b'T',
        _ => unreachable!("two-bit methylation context base"),
    };
    Ok((strand, Some(CytosineContext { class, second })))
}

fn checked_increment(value: &mut u64, label: &str) -> Result<(), CallError> {
    *value = value
        .checked_add(1)
        .ok_or_else(|| CallError::operation(format!("{label} overflowed u64")))?;
    Ok(())
}
