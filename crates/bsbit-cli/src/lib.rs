//! Thin, deterministic command-line orchestration.

#![deny(unsafe_code)]

mod command;
mod cpu_placement;
mod parallel;
mod progress;
mod record_composition;
mod record_replay;
mod record_support;
mod report;

pub use report::{CliError, CliWarning, RunReport};

use std::ffi::OsString;
use std::io::Write;

use command::Action;

/// Top-level command help.
pub const GENERAL_HELP: &str = concat!(
    "bsbit ",
    env!("CARGO_PKG_VERSION"),
    r#"

USAGE:
    bsbit index -r PATH -o PATH [-t N] [--memory-mib N]
                [--index-speed fast|compact] [--simd-backend BACKEND]
    bsbit align -x PATH -1 PATH [-2 PATH] -o PATH [OPTIONS]
    bsbit cpu [--simd-backend auto|scalar|sse2|sse4.2|avx2|avx512|neon]
    bsbit call meth|snp|joint [OPTIONS]
    bsbit combine -i PATH[,PATH...] -o OUTPUT
                  [-c true|false] [-t N] [--matrix level|count|both]
                  [--min-count N] [--min-prop P] [--cg-only]
                  [--sample-name NAME[,NAME...]]

`bsbit align` is the single-end and paired-end entry point.
`bsbit index` creates the complete reference index used by alignment.
`bsbit cpu` reports process-visible CPU features and the selected SIMD backend.
Run `bsbit COMMAND --help` for command details.
"#
);

/// CPU capability diagnostic help.
pub const CPU_HELP: &str = "USAGE:\n    bsbit cpu [--simd-backend auto|scalar|sse2|sse4.2|avx2|avx512|neon]\n\nReports the process architecture, optional SIMD feature bits, and the final backend. A forced backend is validated before any optimized instruction executes.\n";

/// Index-command help.
pub const INDEX_HELP: &str = r"USAGE:
    bsbit index -r PATH -o PATH [-t N] [--memory-mib N]
                [--index-speed fast|compact] [--simd-backend BACKEND]

OPTIONS:
    -r, --reference PATH    plain, gzip, or BGZF reference FASTA
    -o, --output PATH       bsbit reference index path; existing regular file is truncated
    -t, --threads N         positive indexing workers; default: 1
    --memory-mib N          positive search-index memory budget; default: 9300
    --index-speed MODE      fast|compact; default: fast
    --simd-backend BACKEND  auto|scalar|sse2|sse4.2|avx2|avx512|neon; default: auto
    -h, --help              print help and exit

Fast uses sparse suffix-array stride 8 to reduce locate work during alignment. Compact uses stride 16 to reduce index size and mapping RSS. Both layouts produce the same alignment decisions and are readable by the same aligner. Compression is detected from content rather than the filename. CPU/backend selection, reference dimensions, build phases, elapsed time, and completion are logged to stdout. Alignment only opens the completed index and never rebuilds or modifies it.
";

/// Align-command help.
pub const ALIGN_HELP: &str = command::align::HELP;

/// Call-module help.
pub const CALL_HELP: &str = "USAGE:\n    bsbit call meth|snp|joint [OPTIONS]\n\nMODULES:\n    meth     Aggregate strand-aware bisulfite methylation calls\n    snp      Call quality-weighted bisulfite-aware diploid SNVs\n    joint    Produce methylation and SNV outputs from shared fragment evidence\n\nRun `bsbit call MODULE --help` for module details.\n";

/// Methylation-calling help.
pub const CALL_METH_HELP: &str = r"USAGE:
    bsbit call meth -i INPUT.bam -r FASTA -o OUTPUT -f cgmap|bed [OPTIONS]

OPTIONS:
    -i, --input PATH             coordinate-sorted, indexed bsbit BAM
    -r, --reference FASTA        authoritative reference FASTA
    -o, --output PATH            output path; existing regular file is truncated
    -f, --format FORMAT          cgmap|bed
        --region REGION          one 1-based inclusive CONTIG:START-END target
        --regions-bed BED        BED3+ file for multiple regions
    -c, --compress BOOL          true|false; default: false
    -t, --threads N              positive calling-worker count; default: 1
        --compression-threads N  private BGZF workers; default: 0, or 1 with compressed multi-worker output
        --min-bq N               0..=93; default: 20
        --min-mapq N             0..=254; default: 20
        --min-depth N            positive depth; default: 10
        --cg-only                omit CHG and CHH sites
        --ignore-orphan          skip paired records lacking proper-pair
    -h, --help                   print help and exit

Aggregates primary BS-seq methylation calls. The coordinate-sorted BAM must be indexed and contain a standard M5 checksum on every @SQ record. FASTA is authoritative; its dictionary and per-contig M5 values must match the BAM. An uncompressed FASTA uses FAI when available and otherwise scans once and reports a warning; BGZF requires FAI/GZI, and ordinary gzip is unsupported. MD is ignored and every mapped primary record must carry XG:Z:CT or XG:Z:GA.

Defaults are base quality 20, MAPQ 20, and minimum depth 10. --ignore-orphan skips paired records without the SAM proper-pair flag. --cg-only omits CHG and CHH sites. Use --region for one region and --regions-bed for multiple regions. cgmap writes standard 8-column CGmap and bed writes 18-column extended bedMethyl. Output is plain text by default; -c true writes deterministic BGZF. --compression-threads selects private BGZF workers; 0 compresses synchronously and it must be 0 with -c false. The output path is opened directly; an existing regular file is truncated.
";

/// Bisulfite-aware SNP-calling help.
pub const CALL_SNP_HELP: &str = r"USAGE:
    bsbit call snp -i INPUT.bam -r FASTA -o OUTPUT.vcf [OPTIONS]

OPTIONS:
    -i, --input PATH                coordinate-sorted, indexed bsbit BAM
    -r, --reference FASTA           authoritative reference FASTA
    -o, --output PATH               VCF output path; existing regular file is truncated
        --sample-name NAME          override the BAM-derived sample name
        --region REGION             one 1-based inclusive CONTIG:START-END target
        --regions-bed BED           BED3+ file for multiple regions
    -c, --compress BOOL             true|false; default: false
    -t, --threads N                 positive calling-worker count; default: 1
        --compression-threads N     private BGZF workers; default: 0, or 1 with compressed multi-worker output
        --min-bq N                  0..=93; default: 20
        --min-mapq N                0..=254; default: 20
        --min-depth N               positive depth; default: 10
        --min-alt-count N           positive count; default: 2
        --min-alt-fraction P        probability; default: 0.1
        --min-gq N                  0..=99; default: 0
        --min-aq N                  0..=99; default: 30
        --heterozygosity P          probability strictly between 0 and 1
        --underconversion-rate P    probability; default: 0.0025
        --overconversion-rate P     probability; default: 0
        --ignore-orphan             skip paired records lacking proper-pair
    -h, --help                      print help and exit

Calls diploid SNVs with strand-specific bisulfite chemistry. The coordinate-sorted BAM must have a BAI/CSI index and a standard M5 checksum on every @SQ record. FASTA is authoritative; its dictionary and per-contig M5 values must match the BAM. An uncompressed FASTA uses FAI when available and otherwise scans once and reports a warning; BGZF requires FAI/GZI, and ordinary gzip is unsupported. MD is ignored and every mapped primary record must carry XG:Z:CT or XG:Z:GA.

Defaults: base quality 20, MAPQ 20, depth 10, alternate observations 2, alternate fraction 0.1, GQ filter 0, AQ filter 30, heterozygosity 0.001, underconversion 0.0025, and overconversion 0. --ignore-orphan skips paired records without the SAM proper-pair flag. Use --region for one region and --regions-bed for multiple regions. The VCF sample defaults to the unique BAM SM, then the BAM basename; --sample-name overrides it. Output is plain VCF by default; -c true writes deterministic tabix-compatible BGZF. --compression-threads selects private BGZF workers; 0 compresses synchronously and it must be 0 with -c false. The output path is opened directly; an existing regular file is truncated.
";

/// Joint methylation/SNP-calling help.
pub const CALL_JOINT_HELP: &str = r"USAGE:
    bsbit call joint -i INPUT.bam -r FASTA -p PREFIX -f cgmap|bed [OPTIONS]

OPTIONS:
    -i, --input PATH                coordinate-sorted, indexed bsbit BAM
    -r, --reference FASTA           authoritative reference FASTA
    -p, --prefix PATH               output prefix for .CGmap/.bed and .vcf
    -f, --meth-format FORMAT        cgmap|bed
    -c, --compress BOOL             true|false; default: false
    -t, --threads N                 positive calling-worker count; default: 1
        --sample-name NAME          override the BAM-derived VCF sample name
        --region REGION             one 1-based inclusive CONTIG:START-END target
        --regions-bed BED           BED3+ file for multiple regions
        --min-bq N                  0..=93; default: 20
        --min-mapq N                0..=254; default: 20
        --min-depth N               positive depth; default: 10
        --min-alt-count N           positive count; default: 2
        --min-alt-fraction P        probability; default: 0.1
        --min-gq N                  0..=99; default: 0
        --min-aq N                  0..=99; default: 30
        --heterozygosity P          probability strictly between 0 and 1
        --underconversion-rate P    probability; default: 0.0025
        --overconversion-rate P     probability; default: 0
        --cg-only                   omit CHG and CHH methylation sites
        --ignore-orphan             skip paired records lacking proper-pair
        --compression-threads N     private BGZF workers per output; default: 0, or 1 with compressed multi-worker output
    -h, --help                      print help and exit

Produces methylation and SNP outputs from the same overlap-collapsed fragment evidence. The BAM, FASTA, reference-identity, reference-index, region, quality, and XG contracts are the same as the separate callers. An uncompressed FASTA uses FAI when available and otherwise scans once and reports a warning; BGZF requires FAI/GZI, and ordinary gzip is unsupported. Base quality, MAPQ, minimum depth, and --ignore-orphan apply to both outputs. --cg-only limits only the methylation output. Use --region for one region and --regions-bed for multiple regions.

PREFIX produces PREFIX.CGmap or PREFIX.bed together with PREFIX.vcf. Output is plain text by default; -c true appends .gz and writes deterministic BGZF. --compression-threads is the private worker count for each BGZF output; 0 compresses synchronously and it must be 0 with -c false. Both output paths are opened directly, truncating existing regular files; a failed call may leave empty or partial outputs.
";

/// Methylation-matrix combine help.
pub const COMBINE_HELP: &str = r"USAGE:
    bsbit combine -i INPUT[.gz][,INPUT[.gz] ...] -o OUTPUT.bed[.gz]
                  [-c true|false] [-t N] [--matrix level|count|both]
                  [--min-count N] [--min-prop P] [--cg-only]
                  [--sample-name NAME[,NAME...]] [--compression-threads N]

Combines coordinate-sorted 8-column CGmap and 18-column bsbit extended bedMethyl files into BED6-plus-matrix tables. Formats may be mixed across samples; plain, gzip, and BGZF transport are detected from content. Supply exactly one --input comma-separated list and, when naming samples, exactly one matching --sample-name comma-separated list.

--cg-only excludes CHG and CHH sites. A sample cell is valid at --min-count coverage (default 1); a site is retained at --min-prop valid samples (default 0) with at least one valid sample. level emits fractions, count emits methylated and total counts, and both performs one merge. --output is the exact destination for a single matrix and the filename template for --matrix both. Output is plain text by default; -c true writes deterministic BGZF. --compression-threads is the private worker count for each BGZF output; 0 compresses synchronously and it must be 0 with -c false. Existing destination files are replaced atomically after completion.
";

/// Parses and executes one command from arguments excluding program name.
///
/// Help, version text, and human-readable index/alignment progress are written
/// to `output`. Indexing, alignment, and calling write directly to final paths;
/// matrix outputs use staged replacement.
///
/// # Errors
///
/// Returns a stable syntax/unsupported error (exit 2) or an operational error
/// (exit 1). Operational failures may leave empty or partial direct outputs.
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
        Action::Cpu(options) => {
            command::cpu::run(options, output)?;
            Ok(RunReport::default())
        }
        Action::Index(options) => command::index::run(&options, output),
        Action::Align(options) => command::align::run(options, output)
            .map(|()| RunReport::default())
            .map_err(|error| CliError::operation(error.to_string())),
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
