//! Joint methylation and SNP calling over shared first-pass evidence.

mod run;

use std::path::PathBuf;

use crate::meth::OutputFormat;
use crate::region::RegionSelection;
use crate::snp::Parameters;
use crate::{CallError, CallReport, validate_compression_threads, validate_threads};

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
    /// Methylation destination, replacing an existing file after completion.
    pub meth_output: PathBuf,
    /// Methylation output schema.
    pub meth_format: OutputFormat,
    /// VCF destination, replacing an existing file after completion.
    pub vcf_output: PathBuf,
    /// Encode both outputs as BGZF when true, otherwise plain text.
    pub compress: bool,
    /// Positive regional calling worker count.
    pub threads: u64,
    /// Private BGZF workers per output; zero performs compression synchronously.
    pub compression_threads: u32,
    /// Emit only `CpG` sites in the methylation output when true.
    pub cg_only: bool,
    /// SNP filtering and chemistry parameters.
    pub parameters: Parameters,
}

/// Calls methylation and SNPs and writes both outputs directly.
///
/// # Errors
///
/// Returns an operational error for an invalid configuration, input contract,
/// calling or output failure.
pub fn call(options: &Options) -> Result<CallReport, CallError> {
    validate_threads("call joint", options.threads)?;
    validate_compression_threads("call joint", options.compress, options.compression_threads)?;
    options.parameters.validate("call joint")?;
    run::run(options)
}
