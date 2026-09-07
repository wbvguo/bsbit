//! Public input, output, and filtering configuration for matrix assembly.

use std::path::PathBuf;

/// One named methylation sample.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Input {
    /// Unique matrix-column label.
    pub sample: String,
    /// Plain, gzip, or BGZF `CGmap` or extended bedMethyl path.
    pub path: PathBuf,
}

/// Values emitted for each sample and retained site.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum MatrixFormat {
    /// One methylated fraction in `0..=1` per sample.
    Level,
    /// Methylated and total valid coverage columns per sample.
    Count,
    /// Separate level and count matrices produced by one merge.
    Both,
}

/// Per-sample and per-site matrix filters.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Parameters {
    /// Minimum valid methylated-plus-unmethylated coverage for one sample cell.
    pub minimum_count: u64,
    /// Minimum valid-sample proportion, in parts per billion.
    pub minimum_sample_proportion_parts_per_billion: u32,
    /// Retain only `CpG` sites when true.
    pub cg_only: bool,
}

impl Default for Parameters {
    fn default() -> Self {
        Self {
            minimum_count: 1,
            minimum_sample_proportion_parts_per_billion: 0,
            cg_only: false,
        }
    }
}

/// Methylation matrix assembly configuration.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Options {
    /// Ordered sample inputs. This order defines matrix column order.
    pub inputs: Vec<Input>,
    /// Destination, or filename template for `Both`; existing files are replaced.
    ///
    /// `Both` inserts `.level` and `.count` before the recognized BED/gzip
    /// suffix and does not create this path itself.
    pub output: PathBuf,
    /// Per-sample values to emit.
    pub matrix_format: MatrixFormat,
    /// Encode the output as deterministic BGZF when true.
    pub compress: bool,
    /// Positive input merge worker count.
    pub threads: u64,
    /// Private BGZF workers per output; zero performs compression synchronously.
    pub compression_threads: u32,
    /// Coverage and valid-sample filters.
    pub parameters: Parameters,
}
