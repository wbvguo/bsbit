//! Methylation calling over shared fragment evidence.

pub(crate) mod aggregation;
pub(crate) mod output;
mod run;

use std::path::PathBuf;

use crate::region::RegionSelection;
use crate::{CallError, CallReport, validate_compression_threads, validate_threads};

/// Supported methylation output schema.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum OutputFormat {
    /// Standard eight-column `CGmap`.
    Cgmap,
    /// Eighteen-column extended bedMethyl.
    Bed,
}

/// Methylation evidence and site filters.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Parameters {
    /// Minimum observed-base Phred quality in `0..=93`.
    pub minimum_base_quality: u8,
    /// Minimum mapping quality in `0..=254`.
    pub minimum_mapping_quality: u8,
    /// Minimum valid methylated-plus-unmethylated depth; must be nonzero.
    pub minimum_depth: u32,
    /// Emit only `CpG` sites when true.
    pub cg_only: bool,
    /// Ignore paired records that do not carry the SAM proper-pair flag.
    pub ignore_orphans: bool,
}

impl Default for Parameters {
    fn default() -> Self {
        Self {
            minimum_base_quality: 20,
            minimum_mapping_quality: 20,
            minimum_depth: 10,
            cg_only: false,
            ignore_orphans: false,
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
        if self.minimum_depth == 0 {
            return Err(CallError::configuration(format!(
                "{command}: minimum depth must be nonzero"
            )));
        }
        Ok(())
    }
}

/// Methylation-calling configuration.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Options {
    /// Coordinate-sorted, indexed canonical bsbit BAM input.
    pub input: PathBuf,
    /// Indexed FASTA used for authoritative cytosine context.
    pub reference: PathBuf,
    /// Optional interval restriction; empty means the whole BAM dictionary.
    pub regions: RegionSelection,
    /// `CGmap` or BED destination, replacing an existing file after completion.
    pub output: PathBuf,
    /// Output schema.
    pub format: OutputFormat,
    /// Encode output as BGZF when true, otherwise plain text.
    pub compress: bool,
    /// Positive regional calling worker count.
    pub threads: u64,
    /// Private BGZF workers; zero performs compression synchronously.
    pub compression_threads: u32,
    /// Base and mapping quality thresholds.
    pub parameters: Parameters,
}

/// Calls methylation and writes one output directly.
///
/// # Errors
///
/// Returns an operational error for an invalid configuration, input contract,
/// calling or output failure.
pub fn call(options: &Options) -> Result<CallReport, CallError> {
    validate_threads("call meth", options.threads)?;
    validate_compression_threads("call meth", options.compress, options.compression_threads)?;
    options.parameters.validate("call meth")?;
    run::run(options)
}
