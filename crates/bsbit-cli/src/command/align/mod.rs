//! Standard single-end and paired-end alignment command.
//!
//! This module owns only common argument parsing and layout dispatch. Each
//! layout keeps its own validated options and FASTQ-to-BAM runtime in a
//! parallel child module.

mod paired;
mod single;

use std::collections::BTreeSet;
use std::io;
use std::path::PathBuf;

use bsbit_align::library::LibraryProfile;
use bsbit_align::paired_end::PairedSearchMode;
use bsbit_align::single_end::{SINGLE_MAX_EDIT_DISTANCE, SingleSearchMode};
use bsbit_hts::AlignmentAuxiliaryMode;

use self::paired::throughput_thread_split;
use super::ReadLayout;
use crate::{CliError, RunReport};

pub(crate) const HELP: &str = r"bsbit align - standard bisulfite read alignment

USAGE:
  bsbit align --index PATH --read1 PATH [--read2 PATH] --output-bam PATH [OPTIONS]

REQUIRED:
  --index PATH                       complete index created by `bsbit index`
  -1, --read1 PATH                   single-end FASTQ, or R1 FASTQ when paired
  --output-bam PATH                  create-only published BAM path

INPUT LAYOUT:
  --read1 only                       directional single-end alignment
                                      (add --non-directional for four-strand SE)
  --read1 and --read2                synchronized directional paired-end alignment
                                      (add --non-directional for four-strand PE)

OPTIONAL INPUT:
  -2, --read2 PATH                   R2 FASTQ; requires --read1

OPTIONS FOR BOTH LAYOUTS:
  --sensitive                        audit a wider bounded candidate frontier
  --non-directional                  search all four bisulfite strands
  --output-contract CONTRACT         minimal|bismark; default: minimal
  --mapped-only                      omit primary records without a placement
  --threads N                        mapping workers; default: 1
  --bam-threads N                    BGZF workers; default: 1
  --bam-compression-level LEVEL      default|0..9; default: 1
  --metrics                          write the layout-specific profiling TSV

PAIRED-END OPTIONS:
  --total-threads N                  split one 1..64 core budget between mapping
                                      and output; conflicts with both thread flags
  --batch-pairs N                    default: 16384
  --alignment-queue-batches N        default: 2
  --min-template-span N              default: 0
  --max-template-span N              default: 1000

Single-end alignment uses the same persisted combined index and bounded d3/d5
verification core as paired-end alignment. Unique single reads receive numeric
MAPQ from their existing score-separation and repeat evidence; tied best
placements use MAPQ 0. Directional and non-directional single-end BAM contracts
are accepted by `bsbit call` after coordinate sorting, duplicate handling, and
indexing.

Without --sensitive, default mode runs the low-latency d3 pass plus an
incremental d5 fallback. For single-end input, --sensitive preserves that
result as an incumbent and audits it against the wider bounded seed frontier.
A different-origin replacement or new rescue must be unique at MAPQ 20 or
above; a lower-confidence conflict retains the incumbent at MAPQ 0. For
paired-end input --sensitive enables the qualified pair-specific recovery policy.
Inputs may remain gzip-compressed; pre-decompression is not required or recommended.
";

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum SearchMode {
    Default,
    Sensitive,
}

impl SearchMode {
    const fn single(self) -> SingleSearchMode {
        match self {
            Self::Default => SingleSearchMode::Default,
            Self::Sensitive => SingleSearchMode::Sensitive,
        }
    }

    const fn paired(self) -> PairedSearchMode {
        match self {
            Self::Default => PairedSearchMode::Default,
            Self::Sensitive => PairedSearchMode::Sensitive,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub(super) enum ReadOutputMode {
    #[default]
    Complete,
    MappedOnly,
}

impl ReadOutputMode {
    pub(super) const fn name(self) -> &'static str {
        match self {
            Self::Complete => "complete",
            Self::MappedOnly => "mapped-only",
        }
    }
}

pub(super) const fn output_contract_name(mode: AlignmentAuxiliaryMode) -> &'static str {
    match mode {
        AlignmentAuxiliaryMode::Minimal => "minimal",
        AlignmentAuxiliaryMode::Bismark => "bismark",
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
struct ParsedOptions {
    index: PathBuf,
    layout: ReadLayout,
    read1: PathBuf,
    read2: Option<PathBuf>,
    output_bam: PathBuf,
    batch_pairs: usize,
    alignment_queue_batches: usize,
    threads: usize,
    bam_threads: u32,
    auxiliary_core_budget: Option<usize>,
    total_thread_budget: Option<usize>,
    bam_compression_level: Option<u8>,
    output_contract: AlignmentAuxiliaryMode,
    library_profile: LibraryProfile,
    search_mode: SearchMode,
    read_output: ReadOutputMode,
    minimum_template_span: u64,
    maximum_template_span: u64,
    emit_metrics: bool,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) enum Options {
    Single(single::Options),
    Paired(paired::Options),
}

impl ParsedOptions {
    fn into_options(self) -> Options {
        match self.layout {
            ReadLayout::SingleEnd => Options::Single(single::Options {
                index: self.index,
                read1: self.read1,
                output_bam: self.output_bam,
                search_mode: self.search_mode.single(),
                library_profile: self.library_profile,
                max_edit_distance: u64::from(SINGLE_MAX_EDIT_DISTANCE),
                batch_records: 1_000,
                threads: u64::try_from(self.threads).expect("validated thread count fits u64"),
                bam_threads: self.bam_threads,
                bam_compression_level: self.bam_compression_level,
                output_contract: self.output_contract,
                read_output: self.read_output,
                emit_metrics: self.emit_metrics,
            }),
            ReadLayout::PairedEnd => Options::Paired(paired::Options {
                index: self.index,
                read1: self.read1,
                read2: self.read2.expect("paired layout was validated with read 2"),
                output_bam: self.output_bam,
                batch_pairs: self.batch_pairs,
                alignment_queue_batches: self.alignment_queue_batches,
                threads: self.threads,
                bam_threads: self.bam_threads,
                auxiliary_core_budget: self.auxiliary_core_budget,
                total_thread_budget: self.total_thread_budget,
                bam_compression_level: self.bam_compression_level,
                output_contract: self.output_contract,
                library_profile: self.library_profile,
                search_mode: self.search_mode.paired(),
                read_output: self.read_output,
                minimum_template_span: self.minimum_template_span,
                maximum_template_span: self.maximum_template_span,
                emit_metrics: self.emit_metrics,
            }),
        }
    }
}

pub(super) fn parse(arguments: &[String]) -> Result<super::Action, CliError> {
    if matches!(arguments, [value] if value == "--help" || value == "-h") {
        return Ok(super::Action::Help(HELP));
    }
    parse_options_from(arguments.iter().map(std::ffi::OsString::from))
        .map(ParsedOptions::into_options)
        .map(super::Action::Align)
        .map_err(|error| CliError::usage(error.to_string()))
}

pub(crate) fn run(options: Options) -> Result<RunReport, CliError> {
    match options {
        Options::Single(options) => single::run(&options),
        Options::Paired(options) => {
            paired::run(options).map_err(|error| CliError::operation(error.to_string()))
        }
    }
}

fn option_takes_value(flag: &str) -> bool {
    matches!(
        flag,
        "--index"
            | "--read1"
            | "--read2"
            | "--output-bam"
            | "--batch-pairs"
            | "--alignment-queue-batches"
            | "--total-threads"
            | "--threads"
            | "--bam-threads"
            | "--bam-compression-level"
            | "--output-contract"
            | "--min-template-span"
            | "--max-template-span"
    )
}

// Keeping option collection and cross-option validation together makes every
// accepted flag combination auditable without introducing a second state model.
#[allow(clippy::too_many_lines)]
fn parse_options_from(
    args: impl IntoIterator<Item = std::ffi::OsString>,
) -> Result<ParsedOptions, io::Error> {
    let mut index = None;
    let mut read1 = None;
    let mut read2 = None;
    let mut output_bam = None;
    let mut batch_pairs = 16_384_usize;
    let mut alignment_queue_batches = 2_usize;
    let mut threads = 1_usize;
    let mut total_threads = None;
    // One BGZF worker lets record compression overlap mapping.  Zero remains
    // available for callers that require a strictly synchronous writer.
    let mut bam_threads = 1_u32;
    let mut bam_compression_level = Some(1_u8);
    let mut output_contract = AlignmentAuxiliaryMode::Minimal;
    let mut library_profile = LibraryProfile::Directional;
    let mut search_mode = SearchMode::Default;
    let mut explicit_search_mode = None;
    let mut read_output = ReadOutputMode::Complete;
    let mut read_output_explicit = false;
    let mut minimum_template_span = 0_u64;
    let mut maximum_template_span = 1_000_u64;
    let mut emit_metrics = false;
    let mut seen_value_flags = BTreeSet::new();
    let mut args = args.into_iter();
    while let Some(flag) = args.next() {
        let flag = flag
            .to_str()
            .ok_or_else(|| invalid("argument name is not UTF-8"))?;
        let flag = match flag {
            "-1" => "--read1",
            "-2" => "--read2",
            flag => flag,
        };
        if flag == "--sensitive" {
            let requested = SearchMode::Sensitive;
            if let Some(previous) = explicit_search_mode {
                if previous == requested {
                    return Err(invalid(format!("{flag} may be specified only once")));
                }
                return Err(invalid(
                    "--sensitive conflicts with the selected search mode",
                ));
            }
            explicit_search_mode = Some(requested);
            search_mode = requested;
            continue;
        }
        if flag == "--non-directional" {
            if matches!(library_profile, LibraryProfile::NonDirectional) {
                return Err(invalid("--non-directional may be specified only once"));
            }
            library_profile = LibraryProfile::NonDirectional;
            continue;
        }
        if flag == "--mapped-only" {
            if read_output_explicit {
                return Err(invalid("--mapped-only may be specified only once"));
            }
            read_output = ReadOutputMode::MappedOnly;
            read_output_explicit = true;
            continue;
        }
        if flag == "--metrics" {
            if emit_metrics {
                return Err(invalid("--metrics may be specified only once"));
            }
            emit_metrics = true;
            continue;
        }
        if !option_takes_value(flag) {
            return Err(invalid(format!("unknown option {flag}")));
        }
        if !seen_value_flags.insert(flag.to_owned()) {
            return Err(invalid(format!("{flag} may be specified only once")));
        }
        let value = args
            .next()
            .ok_or_else(|| invalid(format!("{flag} requires a value")))?;
        match flag {
            "--index" => index = Some(PathBuf::from(value)),
            "--read1" => read1 = Some(PathBuf::from(value)),
            "--read2" => read2 = Some(PathBuf::from(value)),
            "--output-bam" => output_bam = Some(PathBuf::from(value)),
            "--batch-pairs" => batch_pairs = parse_usize(flag, &value)?,
            "--alignment-queue-batches" => {
                alignment_queue_batches = parse_usize(flag, &value)?;
            }
            "--total-threads" => total_threads = Some(parse_usize(flag, &value)?),
            "--threads" => threads = parse_usize(flag, &value)?,
            "--bam-threads" => bam_threads = parse_u32(flag, &value)?,
            "--bam-compression-level" => {
                if value == "default" {
                    bam_compression_level = None;
                } else {
                    let level = parse_u32(flag, &value)?;
                    if level > 9 {
                        return Err(invalid(
                            "--bam-compression-level must be default or in 0..=9",
                        ));
                    }
                    bam_compression_level =
                        Some(u8::try_from(level).expect("level is at most nine"));
                }
            }
            "--output-contract" => {
                output_contract = parse_output_contract(flag, &value)?;
            }
            "--min-template-span" => minimum_template_span = parse_u64(flag, &value)?,
            "--max-template-span" => maximum_template_span = parse_u64(flag, &value)?,
            _ => unreachable!("value-bearing option was validated above"),
        }
    }
    let auxiliary_core_budget = if let Some(total_threads) = total_threads {
        if seen_value_flags.contains("--threads") || seen_value_flags.contains("--bam-threads") {
            return Err(invalid(
                "--total-threads conflicts with --threads and --bam-threads",
            ));
        }
        if total_threads == 0 || total_threads > 64 {
            return Err(invalid("--total-threads must be in 1..=64"));
        }
        let split = throughput_thread_split(total_threads, false);
        threads = split.0;
        bam_threads = split.1;
        Some(usize::try_from(bam_threads).expect("BGZF thread count fits usize"))
    } else {
        None
    };
    if threads == 0 || threads > 64 {
        return Err(invalid("--threads must be in 1..=64"));
    }
    if bam_threads > 64 {
        return Err(invalid("--bam-threads must be in 0..=64"));
    }
    if batch_pairs == 0 {
        return Err(invalid("--batch-pairs must be positive"));
    }
    if alignment_queue_batches == 0 || alignment_queue_batches > 64 {
        return Err(invalid("--alignment-queue-batches must be in 1..=64"));
    }
    if minimum_template_span > maximum_template_span {
        return Err(invalid(
            "--min-template-span must not exceed --max-template-span",
        ));
    }
    let index = required(index, "--index")?;
    let (layout, read1, read2) = match (read1, read2) {
        (Some(read1), Some(read2)) => (ReadLayout::PairedEnd, read1, Some(read2)),
        (Some(read1), None) => (ReadLayout::SingleEnd, read1, None),
        (None, Some(_)) => return Err(invalid("--read2 requires --read1")),
        (None, None) => return Err(invalid("missing --read1")),
    };
    let output_bam = required(output_bam, "--output-bam")?;
    if matches!(layout, ReadLayout::SingleEnd) {
        let unsupported_flag = [
            "--batch-pairs",
            "--alignment-queue-batches",
            "--total-threads",
            "--min-template-span",
            "--max-template-span",
        ]
        .into_iter()
        .find(|flag| seen_value_flags.contains(*flag));
        if let Some(flag) = unsupported_flag {
            return Err(invalid(format!("{flag} requires paired input via --read2")));
        }
    }
    Ok(ParsedOptions {
        index,
        layout,
        read1,
        read2,
        output_bam,
        batch_pairs,
        alignment_queue_batches,
        threads,
        bam_threads,
        auxiliary_core_budget,
        total_thread_budget: total_threads,
        bam_compression_level,
        output_contract,
        library_profile,
        search_mode,
        read_output,
        minimum_template_span,
        maximum_template_span,
        emit_metrics,
    })
}

fn parse_output_contract(
    flag: &str,
    value: &std::ffi::OsStr,
) -> Result<AlignmentAuxiliaryMode, io::Error> {
    match value.to_str() {
        Some("minimal") => Ok(AlignmentAuxiliaryMode::Minimal),
        Some("bismark") => Ok(AlignmentAuxiliaryMode::Bismark),
        _ => Err(invalid(format!(
            "invalid {flag}; expected minimal or bismark"
        ))),
    }
}

fn required(value: Option<PathBuf>, name: &str) -> Result<PathBuf, io::Error> {
    value.ok_or_else(|| invalid(format!("missing {name}")))
}

fn parse_usize(flag: &str, value: &std::ffi::OsStr) -> Result<usize, io::Error> {
    value
        .to_str()
        .ok_or_else(|| invalid(format!("{flag} value is not UTF-8")))?
        .parse()
        .map_err(|_| invalid(format!("invalid {flag}")))
}

fn parse_u64(flag: &str, value: &std::ffi::OsStr) -> Result<u64, io::Error> {
    value
        .to_str()
        .ok_or_else(|| invalid(format!("{flag} value is not UTF-8")))?
        .parse()
        .map_err(|_| invalid(format!("invalid {flag}")))
}

fn parse_u32(flag: &str, value: &std::ffi::OsStr) -> Result<u32, io::Error> {
    value
        .to_str()
        .ok_or_else(|| invalid(format!("{flag} value is not UTF-8")))?
        .parse()
        .map_err(|_| invalid(format!("invalid {flag}")))
}

fn invalid(message: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidInput, message.into())
}

#[cfg(test)]
use self::paired::{
    MetricsTimer, sensitive_mapq_zero_strategy_id, sensitive_read_complete_strategy_id,
    strategy_id_for,
};

#[cfg(test)]
fn strategy_id(options: &ParsedOptions) -> &'static str {
    strategy_id_for(options.search_mode.paired(), options.read_output)
}

#[cfg(test)]
#[path = "../../../tests/whitebox/align.rs"]
mod tests;
