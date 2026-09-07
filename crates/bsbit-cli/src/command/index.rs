use std::ffi::OsString;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::mpsc::{RecvTimeoutError, sync_channel};
use std::thread;

use bsbit_core::reference::ReferenceSequenceMd5;
use bsbit_cpu::{BackendRequest, initialize};
use bsbit_hts::{
    AlignmentRecordAllocation, AlignmentRecordError, AlignmentRecordLimits, BsbitAlignmentMode,
    BsbitHeaderMetadata, DecodedFastaReader, FastaRecord, SamHeader, SamHeaderReference,
    TextRecordLimits,
};
use bsbit_index::build::combined::{
    CombinedIndexBuildError, CombinedIndexBuildOptions, DEFAULT_COMBINED_INDEX_MEMORY_MIB,
    build_combined_index_from_catalog_replace,
};
use bsbit_index::reference::ContigInput;
use bsbit_index::storage::combined::CombinedIndexSaStride;
use bsbit_index::storage::reference_catalog::write_reference_catalog_direct;
use bsbit_io::open_direct_output_distinct_from;

use crate::progress::{PROGRESS_INTERVAL, ProgressLog};
use crate::{CliError, INDEX_HELP, RunReport};

use super::{
    Action, internal_search_file_prefix, option_map_with_short_options, optional_text,
    optional_u64, required_path,
};

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct IndexOptions {
    pub(crate) reference: PathBuf,
    pub(crate) output: PathBuf,
    pub(crate) threads: u32,
    pub(crate) memory_mib: u64,
    pub(crate) speed: IndexSpeed,
    pub(crate) simd_backend: BackendRequest,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub(crate) enum IndexSpeed {
    #[default]
    Fast,
    Compact,
}

pub(super) fn parse_index(arguments: &[OsString]) -> Result<Action, CliError> {
    if matches!(arguments, [value] if value == "--help" || value == "-h") {
        return Ok(Action::Help(INDEX_HELP));
    }
    let (mut values, _) = option_map_with_short_options(
        arguments,
        &[
            "--reference",
            "--output",
            "--threads",
            "--memory-mib",
            "--index-speed",
            "--simd-backend",
        ],
        &[],
        &[
            ("-r", "--reference"),
            ("-o", "--output"),
            ("-t", "--threads"),
        ],
    )?;
    let reference = required_path(&mut values, "--reference")?;
    let output = required_path(&mut values, "--output")?;
    let threads = optional_u64(&mut values, "--threads")?.unwrap_or(1);
    let threads = u32::try_from(threads)
        .ok()
        .filter(|threads| *threads > 0 && i32::try_from(*threads).is_ok())
        .ok_or_else(|| {
            CliError::usage("--threads must fit a positive signed 32-bit worker count")
        })?;
    let memory_mib =
        optional_u64(&mut values, "--memory-mib")?.unwrap_or(DEFAULT_COMBINED_INDEX_MEMORY_MIB);
    if memory_mib == 0 || memory_mib.checked_mul(1 << 20).is_none() {
        return Err(CliError::usage(
            "--memory-mib must be positive and representable as bytes",
        ));
    }
    let speed = match optional_text(&mut values, "--index-speed")?.as_deref() {
        None | Some("fast") => IndexSpeed::Fast,
        Some("compact") => IndexSpeed::Compact,
        Some(value) => {
            return Err(CliError::usage(format!(
                "invalid value `{value}` for `--index-speed`; expected fast or compact"
            )));
        }
    };
    let simd_backend = optional_text(&mut values, "--simd-backend")?
        .map(|value| {
            value.parse().map_err(|error| {
                CliError::usage(format!("invalid --simd-backend `{value}`: {error}"))
            })
        })
        .transpose()?
        .unwrap_or_default();
    Ok(Action::Index(IndexOptions {
        reference,
        output,
        threads,
        memory_mib,
        speed,
        simd_backend,
    }))
}

pub(crate) fn run(options: &IndexOptions, output: &mut impl Write) -> Result<RunReport, CliError> {
    let output_file = open_direct_output_distinct_from(&options.output, &[&options.reference])
        .map_err(|error| operation_error("open output", &options.output, &error))?;
    let configuration = initialize(options.simd_backend)
        .map_err(|error| CliError::operation(format!("index: CPU backend selection: {error}")))?;
    let mut progress = ProgressLog::start(output, "index", configuration, true);
    progress.phase(
        "configuration",
        format_args!(
            "threads={} memory_mib={} reference={}",
            options.threads,
            options.memory_mib,
            options.reference.display()
        ),
    );
    let (contigs, total_reference_bases) = load_reference(options, &mut progress)?;
    let contig_count = u64::try_from(contigs.len()).expect("supported contig count fits u64");
    validate_alignment_output_reference(
        contigs
            .iter()
            .map(|contig| (contig.name(), contig.sequence().len())),
        AlignmentRecordLimits::default(),
    )
    .map_err(|error| {
        CliError::operation(format!(
            "index: reference {} cannot be represented in alignment BAM output: {error}; correct the FASTA and rerun `bsbit index`",
            options.reference.display()
        ))
    })?;
    progress.phase(
        "alignment-output-validated",
        format_args!("contigs={contig_count} bases={total_reference_bases}"),
    );
    let summary = write_reference_catalog_direct(&contigs, output_file)
        .map_err(|error| operation_error("write output", &options.output, &error))?;
    let semantic_digest = summary.semantic_digest();
    let internal_prefix = internal_search_file_prefix(&options.output);
    let sa_stride = match options.speed {
        IndexSpeed::Fast => CombinedIndexSaStride::Eight,
        IndexSpeed::Compact => CombinedIndexSaStride::Sixteen,
    };
    let build_options = CombinedIndexBuildOptions::new(options.threads)
        .expect("validated CLI thread count is accepted by the index builder")
        .with_sa_stride(sa_stride)
        .with_memory_mib(options.memory_mib)
        .expect("validated CLI memory budget is accepted by the index builder");
    progress.phase(
        "search-index-build",
        format_args!(
            "projected_bases={} threads={} memory_mib={} sa_stride={}",
            total_reference_bases.saturating_mul(2),
            build_options.threads(),
            build_options.memory_mib(),
            build_options.sa_stride().value(),
        ),
    );
    build_search_index_with_progress(
        contigs,
        semantic_digest,
        &internal_prefix,
        build_options,
        total_reference_bases,
        &mut progress,
    )
    .map_err(|error| operation_error("build internal search data for", &options.output, &error))?;
    progress.complete(
        total_reference_bases,
        "bases",
        format_args!("contigs={contig_count} output={}", options.output.display()),
    );
    Ok(RunReport::default())
}

fn build_search_index_with_progress(
    contigs: Vec<ContigInput>,
    semantic_digest: bsbit_core::reference::ReferenceSemanticDigest,
    internal_prefix: &Path,
    build_options: CombinedIndexBuildOptions,
    reference_bases: u64,
    progress: &mut ProgressLog<'_>,
) -> Result<(), CombinedIndexBuildError> {
    thread::scope(|scope| {
        let (sender, receiver) = sync_channel(1);
        let handle = thread::Builder::new()
            .name(String::from("bsbit-index-build"))
            .spawn_scoped(scope, move || {
                let result = build_combined_index_from_catalog_replace(
                    contigs,
                    semantic_digest,
                    internal_prefix,
                    build_options,
                );
                let _ = sender.send(result);
            })
            .map_err(|error| {
                CombinedIndexBuildError::Detail(format!("spawn search-index coordinator: {error}"))
            })?;
        let result = loop {
            match receiver.recv_timeout(PROGRESS_INTERVAL) {
                Ok(result) => break result,
                Err(RecvTimeoutError::Timeout) => progress.phase(
                    "search-index-build",
                    format_args!(
                        "status=running projected_bases={}",
                        reference_bases.saturating_mul(2)
                    ),
                ),
                Err(RecvTimeoutError::Disconnected) => {
                    break Err(CombinedIndexBuildError::Detail(String::from(
                        "search-index coordinator disconnected",
                    )));
                }
            }
        };
        if handle.join().is_err() {
            return Err(CombinedIndexBuildError::Detail(String::from(
                "search-index coordinator panicked",
            )));
        }
        result
    })
}

fn load_reference(
    options: &IndexOptions,
    progress: &mut ProgressLog<'_>,
) -> Result<(Vec<ContigInput>, u64), CliError> {
    let mut reader = DecodedFastaReader::open(&options.reference, reference_text_limits())
        .map_err(|error| operation_error("open reference", &options.reference, &error))?;
    let mut contigs = Vec::new();
    let mut total_reference_bases = 0_u64;
    loop {
        match reader.next_record() {
            Ok(Some(record)) => {
                contigs.try_reserve(1).map_err(|_| {
                    CliError::operation(format!(
                        "index: collect reference {}: allocation failed before record {}",
                        options.reference.display(),
                        record.ordinal().get()
                    ))
                })?;
                let record_bases = record.sequence().len();
                total_reference_bases = total_reference_bases
                    .checked_add(record_bases)
                    .ok_or_else(|| CliError::operation("index: reference base count overflow"))?;
                contigs.push(contig_input_from_fasta_record(&record));
                let loaded_contigs =
                    u64::try_from(contigs.len()).expect("supported contig count fits u64");
                progress.progress(
                    loaded_contigs,
                    "contigs",
                    format_args!("bases={total_reference_bases} phase=read-reference"),
                );
            }
            Ok(None) => break,
            Err(error) => {
                let _ = reader.close();
                return Err(operation_error(
                    "parse reference",
                    &options.reference,
                    &error,
                ));
            }
        }
    }
    reader
        .close()
        .map_err(|error| operation_error("close reference", &options.reference, &error))?;
    let contig_count = u64::try_from(contigs.len()).expect("supported contig count fits u64");
    progress.phase(
        "reference-loaded",
        format_args!("contigs={contig_count} bases={total_reference_bases}"),
    );
    Ok((contigs, total_reference_bases))
}

fn contig_input_from_fasta_record(record: &FastaRecord) -> ContigInput {
    ContigInput::new(
        record.record_name().name().to_vec(),
        record.sequence().clone(),
    )
}

fn validate_alignment_output_reference<'a>(
    references: impl ExactSizeIterator<Item = (&'a [u8], u64)>,
    limits: AlignmentRecordLimits,
) -> Result<(), AlignmentRecordError> {
    let requested = u64::try_from(references.len()).unwrap_or(u64::MAX);
    let mut entries = Vec::new();
    entries.try_reserve_exact(references.len()).map_err(|_| {
        AlignmentRecordError::AllocationFailed {
            allocation: AlignmentRecordAllocation::HeaderReferences,
            requested,
        }
    })?;
    for (ordinal, (name, length)) in references.enumerate() {
        entries.push(
            SamHeaderReference::new(u64::try_from(ordinal).unwrap_or(u64::MAX), name, length)?
                .with_md5(ReferenceSequenceMd5::from_bytes([0; 16])),
        );
    }
    SamHeader::new(entries, limits)?.with_bsbit_metadata(
        BsbitHeaderMetadata::new(BsbitAlignmentMode::NonDirectionalPairedEnd),
        limits,
    )?;
    Ok(())
}

fn operation_error(operation: &str, path: &Path, error: &impl std::fmt::Display) -> CliError {
    CliError::operation(format!("index: {operation} {}: {error}", path.display()))
}

const fn reference_text_limits() -> TextRecordLimits {
    TextRecordLimits::MAX
        .with_max_records(1_000_000)
        .with_max_name_bytes(64_000_000)
        .with_max_description_bytes(1_000_000)
        .with_max_quality_bytes(0)
}

#[cfg(test)]
mod tests {
    use std::io::{BufReader, Cursor};

    use bsbit_hts::FastaReader;

    use bsbit_hts::{AlignmentRecordError, AlignmentRecordLimits, AlignmentRecordResource};

    use super::{
        contig_input_from_fasta_record, reference_text_limits, validate_alignment_output_reference,
    };

    #[test]
    fn fasta_description_is_not_promoted_into_the_contig_name() {
        let input = b">chr1 retained human-readable description\nACGTN\n";
        let mut reader = FastaReader::new(
            BufReader::new(Cursor::new(input.as_slice())),
            reference_text_limits(),
        );
        let record = reader
            .next_record()
            .expect("FASTA parses")
            .expect("one FASTA record exists");
        assert_eq!(record.record_name().name(), b"chr1");
        assert_eq!(
            record.record_name().description(),
            b"retained human-readable description",
        );

        let contig = contig_input_from_fasta_record(&record);
        assert_eq!(contig.name(), b"chr1");
        assert_eq!(contig.sequence().to_ascii(), b"ACGTN");
    }

    #[test]
    fn alignment_output_preflight_rejects_invalid_names_and_lengths() {
        let limits = AlignmentRecordLimits::default();
        assert!(matches!(
            validate_alignment_output_reference([(b"chr,1".as_slice(), 10)].into_iter(), limits),
            Err(AlignmentRecordError::InvalidReferenceNameByte { .. })
        ));
        assert!(matches!(
            validate_alignment_output_reference(
                [(b"chr1".as_slice(), u64::from(i32::MAX as u32) + 1)].into_iter(),
                limits,
            ),
            Err(AlignmentRecordError::ReferenceLengthOutOfRange { .. })
        ));
    }

    #[test]
    fn alignment_output_preflight_applies_aggregate_header_limits() {
        let references = [(b"one".as_slice(), 10), (b"two".as_slice(), 20)];
        let count_error = validate_alignment_output_reference(
            references.into_iter(),
            AlignmentRecordLimits::default().with_max_header_references(1),
        );
        assert!(matches!(
            count_error,
            Err(AlignmentRecordError::LimitExceeded {
                resource: AlignmentRecordResource::HeaderReferences,
                ..
            })
        ));

        let name_error = validate_alignment_output_reference(
            references.into_iter(),
            AlignmentRecordLimits::default().with_max_header_name_bytes(5),
        );
        assert!(matches!(
            name_error,
            Err(AlignmentRecordError::LimitExceeded {
                resource: AlignmentRecordResource::HeaderNameBytes,
                ..
            })
        ));

        let header_error = validate_alignment_output_reference(
            references.into_iter(),
            AlignmentRecordLimits::default().with_max_header_bytes(1),
        );
        assert!(matches!(
            header_error,
            Err(AlignmentRecordError::LimitExceeded {
                resource: AlignmentRecordResource::HeaderBytes,
                ..
            })
        ));
    }
}
