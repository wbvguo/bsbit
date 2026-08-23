//! Methylation calling over shared fragment evidence.

pub(crate) mod aggregation;
pub(crate) mod output;
mod run;

use std::path::PathBuf;

use crate::region::RegionSelection;
use crate::{CallError, CallReport, validate_threads};

/// Supported methylation output schema.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum OutputFormat {
    /// Standard eight-column `CGmap`.
    Cgmap,
    /// Eighteen-column extended bedMethyl.
    Bed,
}

/// Base and mapping quality thresholds for methylation evidence.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Parameters {
    /// Minimum observed-base Phred quality in `0..=93`.
    pub minimum_base_quality: u8,
    /// Minimum mapping quality in `0..=254`.
    pub minimum_mapping_quality: u8,
}

impl Default for Parameters {
    fn default() -> Self {
        Self {
            minimum_base_quality: 15,
            minimum_mapping_quality: 20,
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
    /// Create-only `CGmap` or BED destination.
    pub output: PathBuf,
    /// Output schema.
    pub format: OutputFormat,
    /// Encode output as BGZF when true, otherwise plain text.
    pub compress: bool,
    /// Regional calling workers in `1..=64`.
    pub threads: u64,
    /// Base and mapping quality thresholds.
    pub parameters: Parameters,
}

/// Calls methylation and publishes one create-only output.
///
/// # Errors
///
/// Returns an operational error for an invalid configuration, input contract,
/// calling failure, or output publication failure.
pub fn call(options: &Options) -> Result<CallReport, CallError> {
    validate_threads("call meth", options.threads)?;
    options.parameters.validate("call meth")?;
    run::run(options)
}
