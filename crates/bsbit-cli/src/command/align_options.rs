//! Stable option grammar for the standard alignment command.
//!
//! Parsing remains separate from execution so single-end and paired-end
//! capabilities cannot drift as their pipelines evolve independently.

use std::collections::BTreeSet;
use std::io;
use std::path::PathBuf;

use bsbit_align::library::LibraryProfile;
use bsbit_align::paired_end::PairedSearchMode;
use bsbit_align::single_end::SingleSearchMode;
use bsbit_align::{AlignmentOutputPolicy, MAX_EDIT_DISTANCE, SoftClipMode};
use bsbit_cpu::BackendRequest;
use bsbit_hts::AlignmentAuxiliaryMode;

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub(super) enum SearchMode {
    #[default]
    Default,
    Sensitive,
}

impl SearchMode {
    pub(super) const fn single(self) -> SingleSearchMode {
        match self {
            Self::Default => SingleSearchMode::Default,
            Self::Sensitive => SingleSearchMode::Sensitive,
        }
    }

    pub(super) const fn paired(self) -> PairedSearchMode {
        match self {
            Self::Default => PairedSearchMode::Default,
            Self::Sensitive => PairedSearchMode::Sensitive,
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) enum ReadLayout {
    SingleEnd,
    PairedEnd,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub(super) enum ReadOutputMode {
    #[default]
    Complete,
    MappedOnly,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct Options {
    pub(super) index: PathBuf,
    pub(super) read1: PathBuf,
    pub(super) read2: Option<PathBuf>,
    pub(super) output: PathBuf,
    pub(super) maximum_edit_distance: u8,
    pub(super) output_policy: AlignmentOutputPolicy,
    pub(super) batch_size: usize,
    pub(super) tie_break_seed: u64,
    pub(super) queue_batches: usize,
    pub(super) threads: usize,
    pub(super) total_thread_budget: Option<usize>,
    pub(super) compression_threads: u32,
    pub(super) compression_level: Option<u8>,
    pub(super) output_contract: AlignmentAuxiliaryMode,
    pub(super) library_profile: LibraryProfile,
    pub(super) search_mode: SearchMode,
    pub(super) read_output: ReadOutputMode,
    pub(super) minimum_template_span: u64,
    pub(super) maximum_template_span: u64,
    pub(super) emit_metrics: bool,
    pub(super) simd_backend: BackendRequest,
}

impl Options {
    pub(super) const fn layout(&self) -> ReadLayout {
        if self.read2.is_some() {
            ReadLayout::PairedEnd
        } else {
            ReadLayout::SingleEnd
        }
    }
}

pub(super) fn invalid(message: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidInput, message.into())
}

fn option_takes_value(flag: &str) -> bool {
    matches!(
        flag,
        "--index"
            | "--read1"
            | "--read2"
            | "--output"
            | "--max-edit-distance"
            | "--adapter"
            | "--adapter-min-overlap"
            | "--adapter-max-clip"
            | "--soft-clip"
            | "--max-soft-clip"
            | "--batch-size"
            | "--tie-break-seed"
            | "--queue-batches"
            | "--total-threads"
            | "--threads"
            | "--compression-threads"
            | "--compression-level"
            | "--simd-backend"
            | "--output-contract"
            | "--min-template-span"
            | "--max-template-span"
    )
}

// Keeping option collection and cross-option validation together makes every
// accepted flag combination auditable without introducing a second state model.
#[allow(clippy::too_many_lines)]
pub(super) fn parse_options_from(
    args: impl IntoIterator<Item = std::ffi::OsString>,
) -> Result<Options, io::Error> {
    let mut index = None;
    let mut read1 = None;
    let mut read2 = None;
    let mut output = None;
    let mut maximum_edit_distance = MAX_EDIT_DISTANCE;
    let default_output_policy = AlignmentOutputPolicy::default();
    let mut adapter = default_output_policy.adapter().map(<[u8]>::to_vec);
    let mut adapter_minimum_overlap = default_output_policy.adapter_minimum_overlap();
    let mut adapter_maximum_clip_bases = default_output_policy.adapter_maximum_clip_bases();
    let mut soft_clip_mode = default_output_policy.soft_clip_mode();
    let mut maximum_soft_clip_bases = default_output_policy.maximum_soft_clip_bases();
    let mut batch_size = None;
    let mut tie_break_seed = 0_u64;
    let mut queue_batches = 2_usize;
    let mut threads = 1_usize;
    let mut total_thread_budget = None;
    // One BGZF worker lets record compression overlap mapping.  Zero remains
    // available for callers that require a strictly synchronous writer.
    let mut compression_threads = 1_u32;
    let mut compression_level = Some(1_u8);
    let mut output_contract = AlignmentAuxiliaryMode::Minimal;
    let mut library_profile = LibraryProfile::Directional;
    let mut search_mode = SearchMode::Default;
    let mut sensitive_explicit = false;
    let mut read_output = ReadOutputMode::Complete;
    let mut read_output_explicit = false;
    let mut minimum_template_span = 0_u64;
    let mut maximum_template_span = 1_000_u64;
    let mut emit_metrics = false;
    let mut simd_backend = BackendRequest::Auto;
    let mut seen_value_flags = BTreeSet::new();
    let mut args = args.into_iter();
    while let Some(flag) = args.next() {
        let flag = flag
            .to_str()
            .ok_or_else(|| invalid("argument name is not UTF-8"))?;
        let flag = match flag {
            "-x" => "--index",
            "-1" => "--read1",
            "-2" => "--read2",
            "-o" => "--output",
            "-t" => "--threads",
            flag => flag,
        };
        if flag == "--sensitive" {
            if sensitive_explicit {
                return Err(invalid(format!("{flag} may be specified only once")));
            }
            sensitive_explicit = true;
            search_mode = SearchMode::Sensitive;
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
            "--output" => output = Some(PathBuf::from(value)),
            "--max-edit-distance" => {
                let requested = parse_u32(flag, &value)?;
                if requested > u32::from(MAX_EDIT_DISTANCE) {
                    return Err(invalid(format!(
                        "--max-edit-distance must be in 0..={MAX_EDIT_DISTANCE}"
                    )));
                }
                maximum_edit_distance =
                    u8::try_from(requested).expect("validated edit distance fits u8");
            }
            "--adapter" => {
                let text = value
                    .to_str()
                    .ok_or_else(|| invalid("--adapter value is not UTF-8"))?;
                adapter = match text {
                    "auto" | "illumina" => default_output_policy.adapter().map(<[u8]>::to_vec),
                    "none" => None,
                    sequence => Some(sequence.as_bytes().to_vec()),
                };
            }
            "--adapter-min-overlap" => {
                adapter_minimum_overlap = parse_usize(flag, &value)?;
            }
            "--adapter-max-clip" => {
                adapter_maximum_clip_bases = parse_usize(flag, &value)?;
            }
            "--soft-clip" => {
                soft_clip_mode = match value.to_str() {
                    Some("auto") => SoftClipMode::Auto,
                    Some("none") => SoftClipMode::None,
                    Some("adapter") => SoftClipMode::Adapter,
                    _ => return Err(invalid("--soft-clip must be auto, none, or adapter")),
                };
            }
            "--max-soft-clip" => {
                maximum_soft_clip_bases = parse_usize(flag, &value)?;
            }
            "--batch-size" => batch_size = Some(parse_usize(flag, &value)?),
            "--tie-break-seed" => tie_break_seed = parse_u64(flag, &value)?,
            "--queue-batches" => {
                queue_batches = parse_usize(flag, &value)?;
            }
            "--total-threads" => total_thread_budget = Some(parse_usize(flag, &value)?),
            "--threads" => threads = parse_usize(flag, &value)?,
            "--compression-threads" => compression_threads = parse_u32(flag, &value)?,
            "--compression-level" => {
                if value == "default" {
                    compression_level = None;
                } else {
                    let level = parse_u32(flag, &value)?;
                    if level > 9 {
                        return Err(invalid("--compression-level must be default or in 0..=9"));
                    }
                    compression_level = Some(u8::try_from(level).expect("level is at most nine"));
                }
            }
            "--simd-backend" => {
                let text = value
                    .to_str()
                    .ok_or_else(|| invalid("--simd-backend value is not UTF-8"))?;
                simd_backend = text.parse().map_err(|error| {
                    invalid(format!("invalid --simd-backend `{text}`: {error}"))
                })?;
            }
            "--output-contract" => {
                output_contract = parse_output_contract(flag, &value)?;
            }
            "--min-template-span" => minimum_template_span = parse_u64(flag, &value)?,
            "--max-template-span" => maximum_template_span = parse_u64(flag, &value)?,
            _ => unreachable!("value-bearing option was validated above"),
        }
    }
    if let Some(total_threads) = total_thread_budget {
        if seen_value_flags.contains("--threads")
            || seen_value_flags.contains("--compression-threads")
        {
            return Err(invalid(
                "--total-threads conflicts with --threads and --compression-threads",
            ));
        }
        if total_threads == 0 {
            return Err(invalid("--total-threads must be positive"));
        }
        if u32::try_from(total_threads).is_err() {
            return Err(invalid(
                "--total-threads exceeds the supported u32 worker domain",
            ));
        }
    }
    if threads == 0 {
        return Err(invalid("--threads must be positive"));
    }
    if u32::try_from(threads).is_err() {
        return Err(invalid("--threads exceeds the supported u32 worker domain"));
    }
    if batch_size == Some(0) {
        return Err(invalid("--batch-size must be positive"));
    }
    if i32::try_from(compression_threads).is_err() {
        return Err(invalid(
            "--compression-threads exceeds the native signed 32-bit worker domain",
        ));
    }
    if queue_batches == 0 {
        return Err(invalid("--queue-batches must be positive"));
    }
    if minimum_template_span > maximum_template_span {
        return Err(invalid(
            "--min-template-span must not exceed --max-template-span",
        ));
    }
    let output_policy = AlignmentOutputPolicy::new(
        adapter.as_deref(),
        adapter_minimum_overlap,
        adapter_maximum_clip_bases,
        soft_clip_mode,
        maximum_soft_clip_bases,
    )
    .map_err(|error| invalid(format!("invalid adapter/clipping policy: {error}")))?;
    let index = required(index, "--index")?;
    let (layout, read1, read2) = match (read1, read2) {
        (Some(read1), Some(read2)) => (ReadLayout::PairedEnd, read1, Some(read2)),
        (Some(read1), None) => (ReadLayout::SingleEnd, read1, None),
        (None, Some(_)) => return Err(invalid("--read2 requires --read1")),
        (None, None) => return Err(invalid("missing --read1")),
    };
    let output = required(output, "--output")?;
    let batch_size = batch_size.unwrap_or(match layout {
        ReadLayout::SingleEnd => 1_000,
        ReadLayout::PairedEnd => 16_384,
    });
    if matches!(layout, ReadLayout::SingleEnd) {
        let unsupported_flag = ["--min-template-span", "--max-template-span"]
            .into_iter()
            .find(|flag| seen_value_flags.contains(*flag));
        if let Some(flag) = unsupported_flag {
            return Err(invalid(format!("{flag} requires paired input via --read2")));
        }
    }
    Ok(Options {
        index,
        read1,
        read2,
        output,
        maximum_edit_distance,
        output_policy,
        batch_size,
        tie_break_seed,
        queue_batches,
        threads,
        total_thread_budget,
        compression_threads,
        compression_level,
        output_contract,
        library_profile,
        search_mode,
        read_output,
        minimum_template_span,
        maximum_template_span,
        emit_metrics,
        simd_backend,
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
