//! FASTQ and BAM orchestration for canonical single-end alignment.

use std::path::{Path, PathBuf};
use std::sync::Mutex;
use std::sync::atomic::{AtomicBool, Ordering};

use bsbit_align::materialize::traceback_read_placement;
use bsbit_align::single_end::{
    SingleAlignmentResult, SingleBatchAligner, SingleMappingStatus, SingleSearchMode,
};
use bsbit_core::coordinate::{ReferenceInterval, ReferenceLength};
use bsbit_core::reference::ReferenceSemanticDigest;
use bsbit_hts::TextRecordLimits;
use bsbit_hts::{
    AlignmentRead, AlignmentRecord, AlignmentRecordLimits, BamStagingWriter, BsbitAlignmentMode,
    BsbitProgramProvenance, DecodedFastqReader, FastqRecord, RecordMappingQuality,
    SAM_MAX_QUERY_NAME_BYTES, SamHeader,
};
use bsbit_index::reference::ReferenceIndex;
use bsbit_index::storage::combined::load_combined_reference_catalog;
use bsbit_io::validate_replace_target;

use crate::parallel::{
    DispatchError, ProducerOutcome, WorkDispatcher, WorkerOutcome, run_ordered_parallel,
};
use crate::record_composition::{
    RecordBuildError, build_indexed_single_alignment_record, build_sam_header,
};
use crate::{CliError, CliWarning, RunReport};

use super::{internal_search_file_prefix, unused_staging_path};

const MAX_CLI_READ_BASES: u64 = 1_000_000;
const MAX_CLI_DESCRIPTION_BYTES: u64 = 1_000_000;
const SINGLE_SEARCH_BATCH_SIZE: usize = 64;
/// Validated inputs for canonical single-end alignment.
#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct SingleEndCommandOptions {
    pub(crate) index: PathBuf,
    pub(crate) read1: PathBuf,
    pub(crate) output_bam: PathBuf,
    pub(crate) search_mode: SingleSearchMode,
    pub(crate) max_edit_distance: u64,
    pub(crate) batch_records: u64,
    pub(crate) threads: u64,
    pub(crate) compression_threads: u32,
    pub(crate) compression_level: Option<u8>,
}

pub(crate) fn run_single_end(options: &SingleEndCommandOptions) -> Result<RunReport, CliError> {
    validate_output_target(&options.output_bam)?;
    let staging = unused_staging_path(&options.output_bam, "align", "output")?;
    let (reference, semantic_digest) = load_index(options)?;
    run_single_align(
        options,
        &reference,
        semantic_digest,
        &staging,
        BsbitAlignmentMode::CallerCompatibleDirectionalSingle,
    )
}

fn run_single_align(
    options: &SingleEndCommandOptions,
    reference: &ReferenceIndex,
    semantic_digest: ReferenceSemanticDigest,
    staging: &Path,
    alignment_mode: BsbitAlignmentMode,
) -> Result<RunReport, CliError> {
    if options.threads == 1 {
        run_single_align_scalar(options, reference, semantic_digest, staging, alignment_mode)
    } else {
        run_single_align_parallel(options, reference, semantic_digest, staging, alignment_mode)
    }
}

fn run_single_align_scalar(
    options: &SingleEndCommandOptions,
    reference: &ReferenceIndex,
    semantic_digest: ReferenceSemanticDigest,
    staging: &Path,
    alignment_mode: BsbitAlignmentMode,
) -> Result<RunReport, CliError> {
    let mut reader = DecodedFastqReader::open(&options.read1, read_text_limits())
        .map_err(|error| operation_error("align", "open reads", &options.read1, &error))?;
    let record_limits = AlignmentRecordLimits::default();
    let header = build_align_header(reference, semantic_digest, alignment_mode, record_limits)
        .map_err(|error| {
            operation_error("align", "validate output header", &options.index, &error)
        })?;
    let mut writer = OutputWriter::create(
        &options.output_bam,
        staging,
        &header,
        record_limits,
        options.compression_threads,
        options.compression_level,
    )?;
    let batch_size = physical_batch_size(options.batch_records)?;
    let mut batch = Vec::new();
    let mut aligner = SingleBatchAligner::with_capacity(SINGLE_SEARCH_BATCH_SIZE);
    batch.try_reserve_exact(batch_size).map_err(|_| {
        CliError::operation(format!(
            "align: reserve batch of {} records: allocation failed",
            options.batch_records
        ))
    })?;
    loop {
        match reader.next_record() {
            Ok(Some(record)) => {
                batch.push(record);
                if batch.len() == batch_size {
                    map_and_write_single_batch(
                        reference,
                        &batch,
                        options,
                        &mut aligner,
                        &mut writer,
                    )?;
                    batch.clear();
                }
            }
            Ok(None) => break,
            Err(error) => {
                let _ = reader.close();
                return Err(operation_error(
                    "align",
                    "parse reads",
                    &options.read1,
                    &error,
                ));
            }
        }
    }
    if !batch.is_empty() {
        map_and_write_single_batch(reference, &batch, options, &mut aligner, &mut writer)?;
    }
    reader
        .close()
        .map_err(|error| operation_error("align", "close reads", &options.read1, &error))?;
    writer.finish(&options.output_bam)
}

struct PreparedSingleInput {
    reader: DecodedFastqReader,
    batch_size: usize,
}

fn run_single_align_parallel(
    options: &SingleEndCommandOptions,
    reference: &ReferenceIndex,
    semantic_digest: ReferenceSemanticDigest,
    staging: &Path,
    alignment_mode: BsbitAlignmentMode,
) -> Result<RunReport, CliError> {
    let workers = physical_thread_count(options.threads)?;
    let aligners = (0..workers)
        .map(|_| Mutex::new(SingleBatchAligner::with_capacity(SINGLE_SEARCH_BATCH_SIZE)))
        .collect::<Vec<_>>();
    let record_limits = AlignmentRecordLimits::default();
    run_ordered_parallel(
        workers,
        || {
            let reader = DecodedFastqReader::open(&options.read1, read_text_limits())
                .map_err(|error| operation_error("align", "open reads", &options.read1, &error))?;
            Ok(PreparedSingleInput {
                reader,
                batch_size: physical_batch_size(options.batch_records)?,
            })
        },
        |prepared, dispatcher, cancellation| {
            produce_single_batches(prepared, dispatcher, cancellation, options)
        },
        |worker, input, cancellation| {
            let Ok(mut aligner) = aligners[worker].lock() else {
                return WorkerOutcome::Failed(CliError::operation(
                    "align: single mapping workspace lock was poisoned",
                ));
            };
            map_single_batch_parallel(reference, &input, options, cancellation, &mut aligner)
        },
        || {
            let header =
                build_align_header(reference, semantic_digest, alignment_mode, record_limits)
                    .map_err(|error| {
                        operation_error("align", "validate output header", &options.index, &error)
                    })?;
            OutputWriter::create(
                &options.output_bam,
                staging,
                &header,
                record_limits,
                options.compression_threads,
                options.compression_level,
            )
        },
        write_parallel_records,
        |writer| writer.finish(&options.output_bam),
    )
}

fn produce_single_batches(
    prepared: PreparedSingleInput,
    dispatcher: &mut WorkDispatcher<Vec<FastqRecord>>,
    cancellation: &AtomicBool,
    options: &SingleEndCommandOptions,
) -> ProducerOutcome {
    let PreparedSingleInput {
        mut reader,
        batch_size,
    } = prepared;
    let mut batch = match reserved_batch(batch_size, "single") {
        Ok(batch) => batch,
        Err(error) => return ProducerOutcome::Failed(error),
    };
    loop {
        if cancellation.load(Ordering::Relaxed) || dispatcher.is_cancelled() {
            let _ = reader.close();
            return ProducerOutcome::Cancelled;
        }
        match reader.next_record() {
            Ok(Some(record)) => {
                batch.push(record);
                if batch.len() == batch_size {
                    let next = match reserved_batch(batch_size, "single") {
                        Ok(next) => next,
                        Err(error) => {
                            let _ = reader.close();
                            return ProducerOutcome::Failed(error);
                        }
                    };
                    let complete = std::mem::replace(&mut batch, next);
                    if let Err(error) = dispatcher.send(complete) {
                        let _ = reader.close();
                        return dispatch_failure(error);
                    }
                }
            }
            Ok(None) => break,
            Err(error) => {
                let _ = reader.close();
                return ProducerOutcome::Failed(operation_error(
                    "align",
                    "parse reads",
                    &options.read1,
                    &error,
                ));
            }
        }
    }
    if !batch.is_empty()
        && let Err(error) = dispatcher.send(batch)
    {
        let _ = reader.close();
        return dispatch_failure(error);
    }
    match reader.close() {
        Ok(()) => ProducerOutcome::Completed,
        Err(error) => ProducerOutcome::Failed(operation_error(
            "align",
            "close reads",
            &options.read1,
            &error,
        )),
    }
}

fn reserved_batch<T>(batch_size: usize, label: &str) -> Result<Vec<T>, CliError> {
    let mut batch = Vec::new();
    batch.try_reserve_exact(batch_size).map_err(|_| {
        CliError::operation(format!(
            "align: reserve {label} batch of {batch_size} records: allocation failed"
        ))
    })?;
    Ok(batch)
}

fn dispatch_failure(error: DispatchError) -> ProducerOutcome {
    match error {
        DispatchError::Cancelled => ProducerOutcome::Cancelled,
        DispatchError::Disconnected { ordinal } => ProducerOutcome::Failed(CliError::operation(
            format!("align: parallel work queue disconnected before batch {ordinal}"),
        )),
        DispatchError::OrdinalOverflow => ProducerOutcome::Failed(CliError::operation(
            "align: parallel input batch ordinal overflow",
        )),
    }
}

fn map_single_batch_parallel(
    reference: &ReferenceIndex,
    input: &[FastqRecord],
    options: &SingleEndCommandOptions,
    cancellation: &AtomicBool,
    aligner: &mut SingleBatchAligner,
) -> WorkerOutcome<Vec<AlignmentRecord>> {
    match map_single_records(reference, input, options, Some(cancellation), aligner) {
        Ok(records) => WorkerOutcome::Completed(records),
        Err(_error) if cancellation.load(Ordering::Relaxed) => WorkerOutcome::Cancelled,
        Err(error) => WorkerOutcome::Failed(error),
    }
}

fn write_parallel_records(
    writer: &mut OutputWriter,
    records: Vec<AlignmentRecord>,
) -> Result<(), CliError> {
    for record in records {
        writer.write_record(&record)?;
    }
    Ok(())
}

fn map_and_write_single_batch(
    reference: &ReferenceIndex,
    input: &[FastqRecord],
    options: &SingleEndCommandOptions,
    aligner: &mut SingleBatchAligner,
    writer: &mut OutputWriter,
) -> Result<(), CliError> {
    for record in map_single_records(reference, input, options, None, aligner)? {
        writer.write_record(&record)?;
    }
    Ok(())
}

fn map_single_records(
    reference: &ReferenceIndex,
    input: &[FastqRecord],
    options: &SingleEndCommandOptions,
    cancellation: Option<&AtomicBool>,
    aligner: &mut SingleBatchAligner,
) -> Result<Vec<AlignmentRecord>, CliError> {
    let maximum_edit_distance = u8::try_from(options.max_edit_distance).map_err(|_| {
        CliError::operation(format!(
            "align: single edit distance {} is not representable",
            options.max_edit_distance
        ))
    })?;
    let mut output = Vec::new();
    output.try_reserve_exact(input.len()).map_err(|_| {
        CliError::operation(format!(
            "align: reserve {} single alignment records: allocation failed",
            input.len()
        ))
    })?;
    let mut reads = Vec::with_capacity(SINGLE_SEARCH_BATCH_SIZE);
    for chunk in input.chunks(SINGLE_SEARCH_BATCH_SIZE) {
        if cancellation.is_some_and(|flag| flag.load(Ordering::Relaxed)) {
            return Err(CliError::operation("align: single mapping cancelled"));
        }
        reads.clear();
        reads.extend(chunk.iter().map(|record| record.sequence().bases()));
        let mapped = aligner
            .map_reads_with_mode(
                reference,
                &reads,
                maximum_edit_distance,
                options.search_mode,
            )
            .map_err(|error| CliError::operation(format!("align: map single batch: {error}")))?;
        for (source, result) in chunk.iter().zip(mapped.iter().copied()) {
            if cancellation.is_some_and(|flag| flag.load(Ordering::Relaxed)) {
                return Err(CliError::operation("align: single mapping cancelled"));
            }
            output.push(materialize_single_record(
                reference, source, result, options,
            )?);
        }
    }
    Ok(output)
}

fn materialize_single_record(
    reference: &ReferenceIndex,
    source: &FastqRecord,
    result: SingleAlignmentResult,
    options: &SingleEndCommandOptions,
) -> Result<AlignmentRecord, CliError> {
    let alignment = result
        .placement()
        .map(|placement| {
            let contig_id = reference
                .contig_id(placement.contig_ordinal())
                .map_err(|error| error.to_string())?;
            let contig = reference
                .resolve_contig(&contig_id)
                .map_err(|error| error.to_string())?;
            let interval = ReferenceInterval::new(
                placement.start(),
                placement.end(),
                ReferenceLength::new(contig.sequence().len()),
            )
            .map_err(|error| error.to_string())?;
            traceback_read_placement(
                reference,
                source.sequence(),
                &contig_id,
                interval,
                placement.strand(),
                placement.distance(),
            )
            .map_err(|error| error.to_string())
        })
        .transpose()
        .map_err(|error| {
            CliError::operation(format!(
                "align: materialize record {} from {}: {error}",
                source.ordinal().get(),
                options.read1.display()
            ))
        })?;
    let mapping_quality = match result.status() {
        SingleMappingStatus::Unmapped => RecordMappingQuality::Unmapped,
        SingleMappingStatus::Unique => RecordMappingQuality::Calibrated(result.mapping_quality()),
        SingleMappingStatus::Ambiguous => RecordMappingQuality::Tied,
    };
    build_indexed_single_alignment_record(
        reference,
        source.record_name().name(),
        AlignmentRead::new(source.sequence(), Some(source.quality())),
        alignment.as_ref(),
        mapping_quality,
        AlignmentRecordLimits::default(),
    )
    .map_err(|error| {
        CliError::operation(format!(
            "align: construct record {} from {}: {error}",
            source.ordinal().get(),
            options.read1.display()
        ))
    })
}

fn load_index(
    options: &SingleEndCommandOptions,
) -> Result<(ReferenceIndex, ReferenceSemanticDigest), CliError> {
    let internal_prefix = internal_search_file_prefix(&options.index);
    let threads = physical_thread_count(options.threads)?;
    let loaded =
        load_combined_reference_catalog(&options.index, None, &internal_prefix, threads)
            .map_err(|error| operation_error("align", "validate index", &options.index, &error))?;
    let semantic_digest = loaded.summary().semantic_digest();
    Ok((loaded.into_index(), semantic_digest))
}

fn build_align_header(
    reference: &ReferenceIndex,
    semantic_digest: ReferenceSemanticDigest,
    alignment_mode: BsbitAlignmentMode,
    limits: AlignmentRecordLimits,
) -> Result<SamHeader, RecordBuildError> {
    build_sam_header(reference, limits)?
        .with_bsbit_provenance(
            BsbitProgramProvenance::new(semantic_digest.into_bytes(), alignment_mode),
            limits,
        )
        .map_err(Into::into)
}

struct OutputWriter(BamStagingWriter);

impl OutputWriter {
    fn create(
        target: &Path,
        staging: &Path,
        header: &SamHeader,
        limits: AlignmentRecordLimits,
        compression_threads: u32,
        compression_level: Option<u8>,
    ) -> Result<Self, CliError> {
        let writer = match compression_level {
            Some(level) => BamStagingWriter::create_new_with_threads_and_compression_level(
                staging,
                header,
                limits,
                compression_threads,
                level,
            ),
            None => BamStagingWriter::create_new_with_threads(
                staging,
                header,
                limits,
                compression_threads,
            ),
        };
        writer
            .map(Self)
            .map_err(|error| operation_error("align", "create BAM staging", target, &error))
    }

    fn write_record(&mut self, record: &AlignmentRecord) -> Result<(), CliError> {
        self.0
            .write_record_as_bam(record)
            .map_err(|error| CliError::operation(format!("align: write BAM record: {error}")))
    }

    fn finish(self, target: &Path) -> Result<RunReport, CliError> {
        let completed = self
            .0
            .finish()
            .map_err(|error| operation_error("align", "finalize BAM output", target, &error))?;
        let publication = completed
            .publish_replace(target)
            .map_err(|error| operation_error("align", "publish BAM output", target, &error))?;
        let warning = publication.cleanup_warning().map(|kind| {
            CliWarning::new(format!(
                "BAM output {} was published, but staging {} could not be removed: {kind:?}",
                target.display(),
                publication.staging_path().display()
            ))
        });
        Ok(RunReport::with_warning(warning))
    }
}

fn validate_output_target(path: &Path) -> Result<(), CliError> {
    match validate_replace_target(path) {
        Ok(()) => Ok(()),
        Err(error) if error.kind() == std::io::ErrorKind::NotADirectory => {
            let parent = output_parent(path);
            Err(CliError::operation(format!(
                "output parent {} is not a directory",
                parent.display()
            )))
        }
        Err(error) => Err(operation_error(
            "output",
            "inspect destination",
            path,
            &error,
        )),
    }
}

fn output_parent(path: &Path) -> &Path {
    path.parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."))
}

fn physical_batch_size(logical: u64) -> Result<usize, CliError> {
    usize::try_from(logical).map_err(|_| {
        CliError::operation(format!(
            "align: batch-record count {logical} is not addressable on this host"
        ))
    })
}

fn physical_thread_count(logical: u64) -> Result<usize, CliError> {
    usize::try_from(logical).map_err(|_| {
        CliError::operation(format!(
            "align: thread count {logical} is not addressable on this host"
        ))
    })
}

fn operation_error(
    operation: &str,
    context: &str,
    path: &Path,
    error: &impl std::fmt::Display,
) -> CliError {
    CliError::operation(format!(
        "{operation}: {context} {}: {error}",
        path.display()
    ))
}

const fn read_text_limits() -> TextRecordLimits {
    TextRecordLimits::new(
        MAX_CLI_READ_BASES,
        u64::MAX,
        SAM_MAX_QUERY_NAME_BYTES,
        MAX_CLI_DESCRIPTION_BYTES,
        MAX_CLI_READ_BASES,
        u64::MAX,
        MAX_CLI_READ_BASES,
    )
}
