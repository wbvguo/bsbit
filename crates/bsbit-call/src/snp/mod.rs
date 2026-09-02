//! Public SNP-calling configuration and entry point.

pub(crate) mod candidate;
pub(crate) mod likelihood;
pub(crate) mod output;
pub(crate) mod result;
mod run;
pub(crate) mod vcf;

use std::path::PathBuf;

use crate::region::RegionSelection;
use crate::{CallError, CallReport, validate_threads};

/// SNP filtering and bisulfite-chemistry parameters.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Parameters {
    /// Minimum observed-base Phred quality in `0..=93`.
    pub minimum_base_quality: u8,
    /// Minimum mapping quality in `0..=254`.
    pub minimum_mapping_quality: u8,
    /// Ignore paired records that do not carry the SAM proper-pair flag.
    pub ignore_orphans: bool,
    /// Minimum candidate and likelihood depth; must be nonzero.
    pub minimum_depth: u32,
    /// Minimum candidate and selected-ALT informative observations; must be nonzero.
    pub minimum_alternate_count: u32,
    /// Minimum strongest-ALT fraction for candidate discovery, in parts per billion.
    pub minimum_alternate_fraction_parts_per_billion: u32,
    /// Genotype-quality threshold for the VCF `LowGQ` filter in `0..=99`.
    pub minimum_genotype_quality: u8,
    /// Per-ALT presence-quality threshold for the VCF `LowAQ` filter in `0..=99`.
    pub minimum_allele_quality: u8,
    /// Prior probability that a site differs from the reference, in parts per billion.
    pub heterozygosity_parts_per_billion: u32,
    /// Underconversion probability represented exactly in parts per billion.
    pub underconversion_parts_per_billion: u32,
    /// Overconversion probability represented exactly in parts per billion.
    pub overconversion_parts_per_billion: u32,
}

impl Default for Parameters {
    fn default() -> Self {
        Self {
            minimum_base_quality: 20,
            minimum_mapping_quality: 20,
            ignore_orphans: false,
            minimum_depth: 10,
            minimum_alternate_count: 2,
            minimum_alternate_fraction_parts_per_billion: 100_000_000,
            minimum_genotype_quality: 0,
            minimum_allele_quality: 30,
            heterozygosity_parts_per_billion: 1_000_000,
            underconversion_parts_per_billion: 2_500_000,
            overconversion_parts_per_billion: 0,
        }
    }
}

impl Parameters {
    pub(crate) fn validate(self, command: &str) -> Result<(), CallError> {
        if self.minimum_base_quality > 93 {
            return Err(CallError::configuration(format!(
                "{command}: minimum base quality must be within 0..=93"
            )));
        }
        if self.minimum_mapping_quality > 254 {
            return Err(CallError::configuration(format!(
                "{command}: minimum mapping quality must be within 0..=254"
            )));
        }
        if self.minimum_depth == 0 || self.minimum_alternate_count == 0 {
            return Err(CallError::configuration(format!(
                "{command}: minimum depth and alternate count must be nonzero"
            )));
        }
        if self.minimum_alternate_fraction_parts_per_billion > 1_000_000_000 {
            return Err(CallError::configuration(format!(
                "{command}: minimum alternate fraction must be within 0..=1"
            )));
        }
        if self.minimum_genotype_quality > 99 || self.minimum_allele_quality > 99 {
            return Err(CallError::configuration(format!(
                "{command}: minimum genotype and allele qualities must be within 0..=99"
            )));
        }
        if self.heterozygosity_parts_per_billion == 0
            || self.heterozygosity_parts_per_billion >= 1_000_000_000
        {
            return Err(CallError::configuration(format!(
                "{command}: heterozygosity must be strictly between 0 and 1"
            )));
        }
        if self.underconversion_parts_per_billion > 1_000_000_000
            || self.overconversion_parts_per_billion > 1_000_000_000
        {
            return Err(CallError::configuration(format!(
                "{command}: conversion rates must be within 0..=1"
            )));
        }
        Ok(())
    }
}

/// SNP-calling configuration.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Options {
    /// Coordinate-sorted, indexed canonical bsbit BAM input.
    pub input: PathBuf,
    /// Authoritative FASTA; an existing FAI is used, otherwise plain FASTA is scanned.
    pub reference: PathBuf,
    /// VCF sample name; defaults to the unique BAM `SM`, then its filename stem.
    pub sample_name: Option<String>,
    /// Optional interval restriction; empty means the whole BAM dictionary.
    pub regions: RegionSelection,
    /// VCF destination, replacing an existing file after completion.
    pub output: PathBuf,
    /// Encode output as BGZF when true, otherwise plain VCF.
    pub compress: bool,
    /// Regional calling workers in `1..=64`.
    pub threads: u64,
    /// Quality, depth, and conversion parameters.
    pub parameters: Parameters,
}

/// Calls bisulfite-aware diploid SNVs and publishes one VCF.
///
/// # Errors
///
/// Returns an operational error for an invalid configuration, input contract,
/// likelihood failure, or output publication failure.
pub fn call(options: &Options) -> Result<CallReport, CallError> {
    validate_threads("call snp", options.threads)?;
    options.parameters.validate("call snp")?;
    run::run(options)
}
