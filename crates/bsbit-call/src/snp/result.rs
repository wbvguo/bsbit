//! Shared SNP configuration and variant-domain values.

use super::Parameters;
use crate::CallError;
use crate::evidence::{BaseCode, EvidenceObservation, EvidenceStrand};

pub(super) type Base = BaseCode;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) struct Genotype {
    pub(super) left: Base,
    pub(super) right: Base,
}

impl Genotype {
    pub(super) const fn contains(self, base: Base) -> bool {
        self.left as u8 == base as u8 || self.right as u8 == base as u8
    }

    pub(super) const fn is_conversion_sensitive_on(self, strand: EvidenceStrand) -> bool {
        match strand {
            EvidenceStrand::Top => self.contains(Base::C),
            EvidenceStrand::Bottom => self.contains(Base::G),
        }
    }

    pub(super) const fn is_conversion_sensitive(self) -> bool {
        self.contains(Base::C) || self.contains(Base::G)
    }
}

pub(super) const GENOTYPES: [Genotype; 10] = [
    Genotype {
        left: Base::A,
        right: Base::A,
    },
    Genotype {
        left: Base::A,
        right: Base::C,
    },
    Genotype {
        left: Base::A,
        right: Base::G,
    },
    Genotype {
        left: Base::A,
        right: Base::T,
    },
    Genotype {
        left: Base::C,
        right: Base::C,
    },
    Genotype {
        left: Base::C,
        right: Base::G,
    },
    Genotype {
        left: Base::C,
        right: Base::T,
    },
    Genotype {
        left: Base::G,
        right: Base::G,
    },
    Genotype {
        left: Base::G,
        right: Base::T,
    },
    Genotype {
        left: Base::T,
        right: Base::T,
    },
];

/// Stable filters and bisulfite chemistry parameters for one SNP run.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct SnpConfig {
    pub(crate) minimum_base_quality: u8,
    pub(crate) minimum_mapping_quality: u8,
    pub(crate) ignore_orphans: bool,
    pub(crate) minimum_depth: u32,
    pub(crate) minimum_alternate_count: u32,
    pub(crate) minimum_alternate_fraction_parts_per_billion: u32,
    pub(crate) minimum_genotype_quality: u8,
    pub(crate) minimum_allele_quality: u8,
    /// Prior probability that the sample genotype differs from the reference.
    pub(crate) heterozygosity_rate: f64,
    /// Fraction of unmethylated cytosines that fail to convert.
    pub(crate) underconversion_rate: f64,
    /// Fraction of methylated cytosines that are over-converted.
    pub(crate) overconversion_rate: f64,
}

impl Default for SnpConfig {
    fn default() -> Self {
        Self::from(Parameters::default())
    }
}

impl From<Parameters> for SnpConfig {
    fn from(parameters: Parameters) -> Self {
        Self {
            minimum_base_quality: parameters.minimum_base_quality,
            minimum_mapping_quality: parameters.minimum_mapping_quality,
            ignore_orphans: parameters.ignore_orphans,
            minimum_depth: parameters.minimum_depth,
            minimum_alternate_count: parameters.minimum_alternate_count,
            minimum_alternate_fraction_parts_per_billion: parameters
                .minimum_alternate_fraction_parts_per_billion,
            minimum_genotype_quality: parameters.minimum_genotype_quality,
            minimum_allele_quality: parameters.minimum_allele_quality,
            heterozygosity_rate: f64::from(parameters.heterozygosity_parts_per_billion)
                / 1_000_000_000.0,
            underconversion_rate: f64::from(parameters.underconversion_parts_per_billion)
                / 1_000_000_000.0,
            overconversion_rate: f64::from(parameters.overconversion_parts_per_billion)
                / 1_000_000_000.0,
        }
    }
}

/// One diploid SNV ready for deterministic VCF rendering.
pub(super) const FILTER_LOW_ALTERNATE_DEPTH: u8 = 1 << 0;
pub(super) const FILTER_LOW_GENOTYPE_QUALITY: u8 = 1 << 1;
pub(super) const FILTER_LOW_ALLELE_QUALITY: u8 = 1 << 2;

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct VariantCall {
    pub(crate) position: u32,
    pub(super) reference: Base,
    pub(super) genotype: Genotype,
    pub(crate) depth: u32,
    pub(crate) genotype_quality: u8,
    pub(super) allele_qualities: [u8; 2],
    pub(crate) quality: f64,
    pub(super) conditional_alternate_frequencies: [f64; 2],
    pub(super) phred_likelihoods: [u16; 6],
    pub(super) strand_counts: [[u32; 4]; 2],
    pub(super) filters: u8,
}

impl VariantCall {
    pub(super) fn alternates(&self) -> ([Base; 2], usize) {
        alternate_alleles(self.reference, self.genotype)
    }

    pub(super) fn genotype_indices(&self, alternates: &[Base]) -> (usize, usize) {
        let allele_index = |base: Base| {
            if base == self.reference {
                0
            } else {
                alternates
                    .iter()
                    .position(|alternate| *alternate == base)
                    .map_or(0, |index| index + 1)
            }
        };
        let mut left = allele_index(self.genotype.left);
        let mut right = allele_index(self.genotype.right);
        if left > right {
            std::mem::swap(&mut left, &mut right);
        }
        (left, right)
    }
}

pub(super) fn alternate_alleles(reference: Base, genotype: Genotype) -> ([Base; 2], usize) {
    let mut alternates = [reference; 2];
    let mut count = 0;
    for allele in [genotype.left, genotype.right] {
        if allele != reference && !alternates[..count].contains(&allele) {
            alternates[count] = allele;
            count += 1;
        }
    }
    if count == 2 && alternates[0] > alternates[1] {
        alternates.swap(0, 1);
    }
    (alternates, count)
}

pub(super) fn has_exact_alternate_set(
    reference: Base,
    genotype: Genotype,
    expected: &[Base],
) -> bool {
    let (alternates, count) = alternate_alleles(reference, genotype);
    alternates[..count] == *expected
}

pub(super) fn selected_alleles(reference: Base, alternates: &[Base]) -> ([Base; 3], usize) {
    let mut alleles = [reference; 3];
    for (index, alternate) in alternates.iter().copied().enumerate() {
        alleles[index + 1] = alternate;
    }
    (alleles, alternates.len() + 1)
}

pub(super) fn genotype_index(left: Base, right: Base) -> usize {
    let genotype = if left <= right {
        Genotype { left, right }
    } else {
        Genotype {
            left: right,
            right: left,
        }
    };
    GENOTYPES
        .iter()
        .position(|candidate| *candidate == genotype)
        .expect("canonical base pair has one diploid genotype")
}

pub(super) fn total_allele_depth(strand_counts: [[u32; 4]; 2], allele: Base) -> u32 {
    strand_counts[0][allele.index()].saturating_add(strand_counts[1][allele.index()])
}

pub(super) fn informative_allele_depth(
    strand_counts: [[u32; 4]; 2],
    allele: Base,
    alleles: &[Base],
) -> u32 {
    let top = strand_counts[0][allele.index()];
    let bottom = strand_counts[1][allele.index()];
    if allele == Base::T && alleles.contains(&Base::C) {
        bottom
    } else if allele == Base::A && alleles.contains(&Base::G) {
        top
    } else {
        top.saturating_add(bottom)
    }
}

pub(super) fn has_conversion_confounded_pair(alleles: &[Base]) -> bool {
    (alleles.contains(&Base::C) && alleles.contains(&Base::T))
        || (alleles.contains(&Base::G) && alleles.contains(&Base::A))
}

pub(super) fn filtered_observation(
    observation: EvidenceObservation,
    config: SnpConfig,
) -> Option<(Base, Base, u8, u8)> {
    let reference = Base::from_ascii(observation.reference_base)?;
    let observed = Base::from_ascii(observation.query_base?)?;
    let base_quality = observation.base_quality?;
    let mapping_quality =
        (observation.mapping_quality != u8::MAX).then_some(observation.mapping_quality)?;
    (base_quality >= config.minimum_base_quality
        && mapping_quality >= config.minimum_mapping_quality)
        .then_some((reference, observed, base_quality, mapping_quality))
}

pub(super) fn validate_config(config: SnpConfig) -> Result<(), CallError> {
    if config.minimum_depth == 0 || config.minimum_alternate_count == 0 {
        return Err(CallError::operation(
            "SNP minimum depth and alternate count must be nonzero",
        ));
    }
    if config.minimum_alternate_fraction_parts_per_billion > 1_000_000_000 {
        return Err(CallError::operation(
            "SNP minimum alternate fraction must be within 0..=1",
        ));
    }
    if config.minimum_genotype_quality > 99 || config.minimum_allele_quality > 99 {
        return Err(CallError::operation(
            "SNP minimum genotype and allele qualities must be within 0..=99",
        ));
    }
    if !(0.0..=1.0).contains(&config.underconversion_rate)
        || !(0.0..=1.0).contains(&config.overconversion_rate)
    {
        return Err(CallError::operation(
            "SNP conversion rates must be within 0..=1",
        ));
    }
    if !(0.0..1.0).contains(&config.heterozygosity_rate) || config.heterozygosity_rate == 0.0 {
        return Err(CallError::operation(
            "SNP heterozygosity must be strictly between 0 and 1",
        ));
    }
    Ok(())
}
