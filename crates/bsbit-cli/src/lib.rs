//! Thin, deterministic command-line orchestration.

#![deny(unsafe_code)]

mod command;
mod cpu_placement;
mod parallel;
mod record_composition;
mod report;

pub use report::{CliError, CliWarning, RunReport};

use std::ffi::OsString;
use std::io::Write;

use command::Action;

/// Top-level command help.
pub const GENERAL_HELP: &str = concat!(
    "bsbit ",
    env!("CARGO_PKG_VERSION"),
    "\n\nUSAGE:\n    bsbit index -r PATH -o PATH [-t N] [--index-speed balanced|fast]\n    bsbit align -i PATH -1 PATH [-2 PATH] -o PATH [-t N]\n    bsbit call meth|snp|joint [OPTIONS]\n    bsbit combine --input PATH[,PATH...] [--sample-name NAME[,NAME...]]\n                  --output PATH [--matrix level|count|both]\n                  [--min-count N] [--min-prop P]\n                  [--compress true|false] [--threads N]\n\n`bsbit align` is the single-end and paired-end entry point.\n`bsbit index` creates the complete opaque index consumed by alignment; alignment only opens it.\nRun `bsbit COMMAND --help` for command details.\n"
);

/// Index-command help.
pub const INDEX_HELP: &str = "USAGE:\n    bsbit index -r PATH -o PATH [-t N]\n                [--index-speed balanced|fast]\n\nOPTIONS:\n    -r, --reference PATH    plain or BGZF-compressed reference FASTA\n    -o, --output PATH       generated bsbit alignment index\n    -t, --threads N         indexing threads; default: 1; range: 1..=64\n    --index-speed MODE      balanced|fast; default: balanced\n    -h, --help              print help and exit\n\nBuilds the complete reference index consumed by `bsbit align` from local plain or BGZF-compressed FASTA. Ordinary gzip and other compression formats are rejected; compression is detected from content rather than the filename. OUTPUT is the opaque index handle passed to alignment; its physical layout is internal. An existing output is atomically replaced after the new index is complete. Alignment never constructs or modifies the index.\nIndex speed defaults to balanced. Fast halves the sparse suffix-array stride to reduce mapping locate work, at the cost of a larger index and higher mapping RSS.\nReference size has no policy cap; the current accepted index evidence reaches 5,000,000 normalized bases.\n";

/// Align-command help.
pub const ALIGN_HELP: &str = command::align::HELP;

/// Call-module help.
pub const CALL_HELP: &str = "USAGE:\n    bsbit call meth|snp|joint [OPTIONS]\n\nMODULES:\n    meth     Aggregate strand-aware bisulfite methylation calls\n    snp      Call quality-weighted bisulfite-aware diploid SNVs\n    joint    Produce methylation and SNV outputs from shared fragment evidence\n\nRun `bsbit call MODULE --help` for module details.\n";

/// Methylation-calling help.
pub const CALL_METH_HELP: &str = "USAGE:\n    bsbit call meth -i INPUT.bam -r FASTA -o OUTPUT\n                    -f cgmap|bed [OPTIONS]\n\nAggregates primary BS-seq calls. The BAM must be coordinate sorted and indexed and contain the canonical bsbit @PG header. The required FASTA is the authoritative reference and supplies context across read and region edges; its dictionary must match the BAM. An existing adjacent FAI is used; when FAI is absent, plain FASTA is scanned once to build an in-memory position table without creating a sidecar. BGZF-compressed FASTA requires adjacent FAI/GZI; ordinary gzip FASTA is unsupported. MD is ignored by the caller. Every mapped primary record must carry an XG:Z:CT or XG:Z:GA tag selecting conversion-strand identity. A BAM may contain multiple read groups, but all declared @RG SM values must identify one biological sample. Regions are called in parallel and written in deterministic coordinate order without retaining whole-genome results. Overlapping mates contribute once per fragment and site: evidence able to pass the configured BQ/MAPQ thresholds wins first, then canonical/present evidence and lower combined base-plus-mapping error win, with exact ties resolved in favor of R1. Defaults are base quality 15 and MAPQ 20.\n\nWithout a region option the whole BAM dictionary is called. Repeatable --region uses 1-based inclusive CONTIG:START-END coordinates. --regions-file reads plain/gzip/BGZF BED3+ with 0-based half-open coordinates. Both inputs form a merged union, so overlapping targets are counted once. cgmap writes standard 8-column, 1-based CGmap; bed writes 18-column extended bedMethyl. Output is plain text by default; -c true writes deterministic BGZF. Existing output files are atomically replaced.\n";

/// Bisulfite-aware SNP-calling help.
pub const CALL_SNP_HELP: &str = "USAGE:\n    bsbit call snp -i INPUT.bam -r FASTA -o OUTPUT.vcf\n                   [-c true|false] [-t N] [--sample-name NAME]\n                   [--region CONTIG:START-END ...] [--regions-file BED]\n                   [--min-base-quality N] [--min-mapq N] [--min-depth N]\n                   [--min-alt-count N] [--min-alt-fraction P]\n                   [--min-gq N] [--min-aq N] [--heterozygosity P]\n                   [--underconversion-rate P] [--overconversion-rate P]\n\nCalls diploid SNVs with strand-specific bisulfite chemistry, base quality, mapping quality, and adaptive methylation marginalization. The reference-centered heterozygosity prior controls site and ALT identity; after ALT selection, GT dosage is maximum-likelihood and is not biased by the rare-site prior. The coordinate-sorted BAM must have a BAI/CSI index and the canonical bsbit @PG line. The required FASTA is the authoritative reference and its dictionary must match the BAM. An existing adjacent FAI is used; when FAI is absent, plain FASTA is scanned once to build an in-memory position table without creating a sidecar. BGZF-compressed FASTA requires adjacent FAI/GZI; ordinary gzip FASTA is unsupported. MD is ignored by the caller. Every mapped primary record must carry an XG:Z:CT or XG:Z:GA tag selecting conversion-strand identity. Multiple read groups are allowed only when their declared @RG SM values identify one biological sample; multiple distinct SM values are rejected. The VCF sample defaults to the unique SM, then the BAM filename stem; --sample-name only renames that one sample. Candidate discovery uses dense bit-sliced regional counters, exact likelihoods use bounded candidate windows, overlapping mates are collapsed once per fragment/site, and completed regions stream in deterministic order.\n\nWithout a region option the whole BAM dictionary is called. Repeatable --region uses 1-based inclusive CONTIG:START-END coordinates. --regions-file reads plain/gzip/BGZF BED3+ with 0-based half-open coordinates; both form a merged union. Defaults: base quality 15, MAPQ 20, depth 4, alternate observations 2, alternate fraction 0.1, GQ filter disabled (0), AQ filter 30, heterozygosity 0.001, underconversion 0.0025, overconversion 0. -c true writes deterministic, tabix-compatible BGZF. Use --min-gq 20 when FILTER=PASS must also require a high-confidence complete genotype. Existing output files are atomically replaced.\n";

/// Joint methylation/SNP-calling help.
pub const CALL_JOINT_HELP: &str = "USAGE:\n    bsbit call joint -i INPUT.bam -r FASTA\n                     -m METH_OUTPUT -f cgmap|bed -v OUTPUT.vcf\n                     [-c true|false] [-t N] [--sample-name NAME]\n                     [--region CONTIG:START-END ...] [--regions-file BED]\n                     [--min-base-quality N] [--min-mapq N] [--min-depth N]\n                     [--min-alt-count N] [--min-alt-fraction P]\n                     [--min-gq N] [--min-aq N] [--heterozygosity P]\n                     [--underconversion-rate P] [--overconversion-rate P]\n\nProduces both outputs from the same overlap-collapsed first-pass fragment evidence. The required FASTA is authoritative for methylation context and SNP alleles, and its dictionary must match the BAM. An existing adjacent FAI is used; when FAI is absent, plain FASTA is scanned once to build an in-memory position table without creating a sidecar. BGZF-compressed FASTA requires adjacent FAI/GZI; ordinary gzip FASTA is unsupported. MD is ignored by the caller. Every mapped primary record must carry an XG:Z:CT or XG:Z:GA tag selecting conversion-strand identity. Base-quality and MAPQ thresholds apply to both outputs; remaining quality and chemistry options control SNP calling. Multiple read groups are allowed only for one declared SM. The VCF sample defaults to that unique SM, then the BAM filename stem; --sample-name renames it. Without region options the whole dictionary is called. Repeatable 1-based inclusive --region and 0-based half-open BED3+ --regions-file targets are merged as a union. The coordinate-sorted BAM must be indexed and carry the canonical bsbit @PG line. -c true writes deterministic BGZF. The two destinations must differ; existing files are replaced and roll back together if publication cannot complete.\n";

/// Methylation-matrix combine help.
pub const COMBINE_HELP: &str = "USAGE:\n    bsbit combine -i INPUT[.gz][,INPUT[.gz] ...]\n                  [--sample-name NAME[,NAME...]]\n                  -o OUTPUT.bed[.gz] [-m level|count|both]\n                  [--min-count N] [--min-prop P]\n                  [-c true|false] [-t N]\n\nCombines coordinate-sorted 8-column CGmap and 18-column bsbit extended bedMethyl files into one or two BED6-plus-matrix tables. Formats may be mixed across samples. The row schema and plain, gzip, or BGZF transport are detected from content. Both --input and --sample-name accept comma-separated values and may be repeated; values are expanded in declaration order. If --sample-name is omitted, the exact input path is used as the sample label. The number of supplied names must match the number of inputs, and labels must be unique. Commas therefore cannot be embedded in paths or labels.\n\nA sample cell is valid when methylated-plus-unmethylated coverage is at least --min-count (default 1). A site is retained when at least --min-prop of all samples are valid (default 0); at least one valid sample is always required. Low-coverage and absent cells are written as `.` rather than zero. level emits one fraction per sample. count emits SAMPLE_meth_count and SAMPLE_total_count. both performs one merge and writes separate level and count files: OUTPUT cohort.bed.gz becomes cohort.level.bed.gz and cohort.count.bed.gz.\n\nInput workers perform a bounded-memory hierarchical k-way merge. Memory is proportional to the sample count, not the number of genomic sites. A parallel preflight derives and validates a common contig order because methylation tables have no sequence dictionary. Threads default to 1 and accept 1 through 64. Output is plain text unless -c true selects deterministic BGZF. Existing destinations are replaced; both outputs are rolled back together if publication cannot complete.\n";

/// Parses and executes one command from arguments excluding program name.
///
/// Help and version text is written to `output`. Operational commands publish
/// only through staged replacement lifecycles.
///
/// # Errors
///
/// Returns a stable syntax/unsupported error (exit 2) or an operational error
/// (exit 1). Operational failures before publication preserve prior outputs.
pub fn run(
    arguments: impl IntoIterator<Item = OsString>,
    output: &mut impl Write,
) -> Result<RunReport, CliError> {
    match command::parse(arguments)? {
        Action::Help(help) => {
            output
                .write_all(help.as_bytes())
                .map_err(|error| CliError::operation(format!("write help: {error}")))?;
            Ok(RunReport::default())
        }
        Action::Version => {
            writeln!(output, "bsbit {}", env!("CARGO_PKG_VERSION"))
                .map_err(|error| CliError::operation(format!("write version: {error}")))?;
            Ok(RunReport::default())
        }
        Action::Index(options) => command::index::run(&options),
        Action::Align(options) => command::align::run(options),
        Action::CallMeth(options) => adapt_call_report(bsbit_call::meth::call(&options)),
        Action::CallSnp(options) => adapt_call_report(bsbit_call::snp::call(&options)),
        Action::CallJoint(options) => adapt_call_report(bsbit_call::joint::call(&options)),
        Action::Combine(options) => adapt_combine_report(bsbit_combine::combine(&options)),
    }
}

fn adapt_call_report(
    result: Result<bsbit_call::CallReport, bsbit_call::CallError>,
) -> Result<RunReport, CliError> {
    let report = result.map_err(|error| CliError::operation(error.to_string()))?;
    Ok(RunReport {
        warnings: report
            .warnings()
            .iter()
            .map(|warning| CliWarning::new(warning.message()))
            .collect(),
    })
}

fn adapt_combine_report(
    result: Result<bsbit_combine::CombineReport, bsbit_combine::CombineError>,
) -> Result<RunReport, CliError> {
    let report = result.map_err(|error| CliError::operation(error.to_string()))?;
    Ok(RunReport {
        warnings: report
            .warnings()
            .iter()
            .map(|warning| CliWarning::new(warning.message()))
            .collect(),
    })
}
