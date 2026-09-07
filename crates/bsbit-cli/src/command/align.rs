//! Standard FASTQ-to-BAM alignment command.
//!
//! Single-end and paired-end input share the persisted combined index, bounded
//! d3/d5 verification core, canonical traceback, record construction, BAM
//! compression/finalization, and direct output path.

use std::error::Error;
use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};

use bsbit_cpu::initialize;
use bsbit_index::storage::combined::{CombinedIndexSaStride, combined_index_sa_stride_hint};

#[cfg(test)]
use super::align_options::SearchMode;
use super::align_options::{Options, ReadLayout, ReadOutputMode, parse_options_from};
use super::alignment_metrics::MetricsTimer;
#[cfg(test)]
use super::alignment_metrics::{
    sensitive_mapq_zero_strategy_id, sensitive_read_complete_strategy_id, strategy_id,
};
use super::internal_search_file_prefix;
use super::single_end::{SingleEndCommandOptions, run_single_end};
use crate::progress::ProgressLog;

pub(crate) const HELP: &str = r"bsbit align - standard bisulfite read alignment

USAGE:
  bsbit align -x PATH -1 PATH [-2 PATH] -o PATH [OPTIONS]

REQUIRED:
  -x, --index PATH                   complete index created by `bsbit index`
  -1, --read1 PATH                   single-end FASTQ, or R1 FASTQ when paired
  -o, --output PATH                  BAM path; an existing regular file is truncated

INPUT LAYOUT:
  --read1 only                       directional single-end alignment
                                      (add --non-directional for four-strand SE)
  --read1 and --read2                synchronized directional paired-end alignment
                                      (add --non-directional for four-strand PE)

OPTIONAL INPUT:
  -2, --read2 PATH                   R2 FASTQ; requires --read1

OPTIONS FOR BOTH LAYOUTS:
  --sensitive                        enable the layout-specific sensitive search contract
  --non-directional                  search all four bisulfite strands
  --max-edit-distance N              per-read edit budget in 0..=5; default: 5
  --adapter auto|none|illumina|SEQ   exact 3' adapter; default: auto (Illumina universal)
  --adapter-min-overlap N            minimum exact adapter support; default: 8
  --adapter-max-clip N               maximum 3' adapter scan; default: 30
  --soft-clip auto|none|adapter      endpoint repair policy; default: auto
  --max-soft-clip N                  maximum clipped bases per read; default: 30
  --output-contract CONTRACT         minimal|bismark; default: minimal
  --mapped-only                      omit primary records without a placement
  -t, --threads N                    mapping workers; default: 1
  --compression-threads N            BGZF workers; 0 is synchronous; default: 1
  --total-threads N                  split one positive core budget between mapping
                                      and BAM output; conflicts with both thread flags
  --compression-level LEVEL          default|0..9; default: 1
  --batch-size N                     reads or read pairs per batch;
                                      default: 1000 (single), 16384 (paired)
  --tie-break-seed N                 seed for fair MAPQ-0 coordinate ties; default: 0
  --queue-batches N                  bounded batches between pipeline stages;
                                      default: 2
  --simd-backend BACKEND             auto|scalar|sse2|sse4.2|avx2|avx512|neon; default: auto
  --metrics                          suppress human logs and write profiling TSV to stdout

PAIRED-END OPTIONS:
  --min-template-span N              default: 0
  --max-template-span N              default: 1000

Single-end alignment uses the same persisted combined index and bounded d3/d5
verification core as paired-end alignment. Unique single reads receive numeric
MAPQ from their existing score-separation and repeat evidence; unresolved tied
placements use MAPQ 0. The BAM records its single-end layout and directional or
non-directional library profile as informational header metadata. `bsbit call`
accepts it after coordinate sorting, duplicate handling, indexing, and the
same record and reference-identity checks used for paired-end input.

Without --sensitive, default mode runs the low-latency d3 pass plus an
incremental d5 fallback. For single-end input, --sensitive preserves that
result as an incumbent and audits every unique placement against a ten-round
bounded seed frontier. A different-origin replacement must have a positive
policy-declared score; a lower-confidence conflict retains the incumbent at
MAPQ 0. A parsimonious low-edit origin with independent seed support may enter
the completed frontier audit. Final confidence is otherwise derived from the
completed best/second-best score frontier; sensitive mode does not apply
benchmark- or library-profile-specific MAPQ promotion tables.
Single reads with verified equal-best origins retain a deterministic,
read-dependent hash-selected representative at MAPQ 0, including highly repetitive reads;
only reads without a verified placement remain unmapped. For paired-end input, --sensitive
enables the qualified pair-specific recovery policy. A recovered pair that cannot earn
positive MAPQ remains coordinate-bearing, ambiguous, and MAPQ 0; only pairs without a
verified placement remain unmapped. Pair MAPQ uses origin-grouped score separation;
stability, repeat-risk, and clipping checks are hard caps and can never raise confidence.
MAPQ 40 and above require a completed frontier plus an independent origin certificate.
Inputs may remain gzip-compressed; pre-decompression is not required or recommended.
Every FASTQ sequence and quality must contain 3..=192 bytes; both layouts use
the same limit and reject a record before it reaches an alignment backend.
Normal runs log CPU/backend selection, phases, processed reads/pairs, throughput,
and elapsed time to stdout. Progress is rate-limited to one update per five seconds.
";

pub(super) fn parse(arguments: &[std::ffi::OsString]) -> Result<super::Action, crate::CliError> {
    if matches!(arguments, [value] if value == "--help" || value == "-h") {
        return Ok(super::Action::Help(HELP));
    }
    parse_options_from(arguments.iter().cloned())
        .map(super::Action::Align)
        .map_err(|error| crate::CliError::usage(error.to_string()))
}

pub(crate) fn run(mut options: Options, output: &mut impl Write) -> Result<(), Box<dyn Error>> {
    let output_file = open_output(&options)?;
    let configuration = initialize(options.simd_backend)?;
    apply_total_thread_budget(&mut options);
    let process_started = MetricsTimer::start(options.emit_metrics);
    let summary = {
        let mut progress =
            ProgressLog::start(output, "align", configuration, !options.emit_metrics);
        log_alignment_configuration(&options, &mut progress);
        if matches!(options.layout(), ReadLayout::SingleEnd) {
            return run_standard_single_from_options(
                options,
                configuration,
                output_file,
                &mut progress,
            );
        }
        super::paired_end::run_paired_alignment(&options, output_file, &mut progress)?
    };
    super::paired_end::write_metrics(
        output,
        &options,
        configuration,
        &summary,
        process_started.elapsed_ns(),
    )?;
    Ok(())
}

fn open_output(options: &Options) -> Result<File, crate::CliError> {
    let internal_prefix = internal_search_file_prefix(&options.index);
    let internal_components = [
        internal_prefix.clone(),
        append_path_suffix(&internal_prefix, ".bwt"),
        append_path_suffix(&internal_prefix, ".sa"),
        append_path_suffix(&internal_prefix, ".occ"),
    ];
    let mut protected = vec![options.index.as_path(), options.read1.as_path()];
    if let Some(read2) = options.read2.as_deref() {
        protected.push(read2);
    }
    protected.extend(internal_components.iter().map(PathBuf::as_path));
    bsbit_io::open_direct_output_distinct_from(&options.output, &protected).map_err(|error| {
        crate::CliError::operation(format!(
            "align: open output {}: {error}",
            options.output.display()
        ))
    })
}

fn append_path_suffix(path: &Path, suffix: &str) -> PathBuf {
    let mut value = path.as_os_str().to_owned();
    value.push(suffix);
    PathBuf::from(value)
}

fn apply_total_thread_budget(options: &mut Options) {
    let Some(total_threads) = options.total_thread_budget else {
        return;
    };
    let internal_prefix = internal_search_file_prefix(&options.index);
    let fast_index = matches!(
        combined_index_sa_stride_hint(&internal_prefix),
        Some(CombinedIndexSaStride::Eight)
    );
    let (mapping_threads, compression_threads) = throughput_thread_split(total_threads, fast_index);
    options.threads = mapping_threads;
    options.compression_threads = compression_threads;
}

fn throughput_thread_split(total_threads: usize, fast_index: bool) -> (usize, u32) {
    debug_assert!(total_threads > 0);
    if total_threads == 1 {
        return (1, 0);
    }
    let output_share = if fast_index { 4 } else { 5 };
    let output_threads = total_threads.div_ceil(output_share).min(total_threads - 1);
    (
        total_threads - output_threads,
        u32::try_from(output_threads).expect("validated total thread count fits u32"),
    )
}

fn log_alignment_configuration(options: &Options, progress: &mut ProgressLog<'_>) {
    let layout = match options.layout() {
        ReadLayout::SingleEnd => "single-end",
        ReadLayout::PairedEnd => "paired-end",
    };
    progress.phase(
        "configuration",
        format_args!(
            "layout={layout} mapping_threads={} compression_threads={} tie_break_seed={} index={}",
            options.threads,
            options.compression_threads,
            options.tie_break_seed,
            options.index.display(),
        ),
    );
}

fn run_standard_single_from_options(
    options: Options,
    configuration: bsbit_cpu::Configuration,
    output_file: File,
    progress: &mut ProgressLog<'_>,
) -> Result<(), Box<dyn Error>> {
    let align_options = SingleEndCommandOptions {
        index: options.index,
        read1: options.read1,
        output: options.output,
        max_edit_distance: options.maximum_edit_distance,
        output_policy: options.output_policy,
        batch_records: u64::try_from(options.batch_size).expect("validated batch size fits u64"),
        tie_break_seed: options.tie_break_seed,
        queue_batches: options.queue_batches,
        threads: u64::try_from(options.threads).expect("validated thread count fits u64"),
        compression_threads: options.compression_threads,
        compression_level: options.compression_level,
        search_mode: options.search_mode.single(),
        library_profile: options.library_profile,
        output_contract: options.output_contract,
        mapped_only: matches!(options.read_output, ReadOutputMode::MappedOnly),
        configuration,
        emit_metrics: options.emit_metrics,
    };
    run_single_end(&align_options, output_file, progress)
        .map(|_| ())
        .map_err(|error| Box::new(error) as Box<dyn Error>)
}

#[cfg(test)]
#[path = "../../tests/whitebox/align.rs"]
mod tests;
