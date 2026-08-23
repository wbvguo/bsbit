//! Joint methylation and SNP calling over shared first-pass evidence.

mod run;

use std::path::PathBuf;

use crate::meth::OutputFormat;
use crate::region::RegionSelection;
use crate::snp::Parameters;
use crate::{CallError, CallReport, validate_threads};

/// Joint-calling configuration.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Options {
    /// Coordinate-sorted, indexed canonical bsbit BAM input.
    pub input: PathBuf,
    /// Indexed FASTA used for authoritative cytosine context.
    pub reference: PathBuf,
    /// VCF sample name; defaults to the unique BAM `SM`, then its basename.
    pub sample_name: Option<String>,
    /// Optional interval restriction; empty means the whole BAM dictionary.
    pub regions: RegionSelection,
    /// Create-only methylation destination.
    pub meth_output: PathBuf,
    /// Methylation output schema.
    pub meth_format: OutputFormat,
    /// Create-only VCF destination.
    pub vcf_output: PathBuf,
    /// Encode both outputs as BGZF when true, otherwise plain text.
    pub compress: bool,
    /// Regional calling workers in `1..=64`.
    pub threads: u64,
    /// SNP filtering and chemistry parameters.
    pub parameters: Parameters,
}

/// Calls methylation and SNPs and publishes both outputs transactionally.
///
/// # Errors
///
/// Returns an operational error for an invalid configuration, input contract,
/// calling failure, or joint publication failure.
pub fn call(options: &Options) -> Result<CallReport, CallError> {
    validate_threads("call joint", options.threads)?;
    options.parameters.validate("call joint")?;
    run::run(options)
}
