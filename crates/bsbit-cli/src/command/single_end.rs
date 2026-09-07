//! FASTQ and BAM orchestration for canonical single-end alignment.

use std::fs::File;
use std::path::{Path, PathBuf};
use std::sync::Mutex;
use std::sync::atomic::{AtomicBool, Ordering};
use std::time::Instant;

use crate::parallel::{
    DispatchError, ProducerOutcome, WorkDispatcher, WorkerOutcome, run_ordered_parallel,
};
use crate::record_composition::{AlignmentRecordComposer, build_sam_header};
use crate::record_support::RecordBuildError;
use crate::{CliError, RunReport};
use bsbit_align::AlignmentOutputPolicy;
use bsbit_align::library::LibraryProfile;
use bsbit_align::materialize::traceback_read_placement;
use bsbit_align::single_end::{
    SINGLE_ALIGNMENT_BATCH_SIZE, SingleAlignmentResult, SingleBatchAligner, SingleMappingStatus,
    SingleSearchMode,
};
use bsbit_align::{ALIGNMENT_POLICY_ID, MAPQ_POLICY_ID};
use bsbit_core::coordinate::{ReferenceInterval, ReferenceLength};
use bsbit_core::sequence::NormalizedSequence;
use bsbit_hts::{
    AlignmentAuxiliaryMode, AlignmentPlacement, AlignmentRecordBatch, AlignmentRecordLimits,
    BamStagingWriter, BorrowedAlignmentRead, BorrowedFastqRecord, BsbitAlignmentMode,
    BsbitHeaderMetadata, DecodedFastqReader, FastqRecordBatch, SamHeader,
};
use bsbit_index::reference::{ContigId, ReferenceIndex};
use bsbit_index::storage::combined::load_combined_reference_catalog;

use super::alignment_input::alignment_fastq_limits;
use super::alignment_metrics::{
    MetricsRow, adapter_name, library_profile_name, output_contract_name, soft_clip_mode_name,
};
use super::internal_search_file_prefix;
use crate::progress::ProgressLog;

const SINGLE_END_METRICS_SCHEMA: &str = "bsbit-alignment-metrics-single-end-v1";

#[derive(Clone, Copy, Debug, Default)]
struct Observation {
    reads: u64,
    unique: u64,
    ambiguous: u64,
    unmapped: u64,
    located_rows: u64,
    verified_placements: u64,
    adapter_attempted: u64,
    adapter_unique: u64,
    adapter_ambiguous: u64,
    adapter_unmapped: u64,
    adapter_clipped_bases: u64,
    direct_records: u64,
    traceback_records: u64,
    mapping_worker_ns: u128,
    record_worker_ns: u128,
}

impl Observation {
    fn merge(&mut self, other: Self) {
        self.reads = self.reads.saturating_add(other.reads);
        self.unique = self.unique.saturating_add(other.unique);
        self.ambiguous = self.ambiguous.saturating_add(other.ambiguous);
        self.unmapped = self.unmapped.saturating_add(other.unmapped);
        self.located_rows = self.located_rows.saturating_add(other.located_rows);
        self.verified_placements = self
            .verified_placements
            .saturating_add(other.verified_placements);
        self.adapter_attempted = self
            .adapter_attempted
            .saturating_add(other.adapter_attempted);
        self.adapter_unique = self.adapter_unique.saturating_add(other.adapter_unique);
        self.adapter_ambiguous = self
            .adapter_ambiguous
            .saturating_add(other.adapter_ambiguous);
        self.adapter_unmapped = self.adapter_unmapped.saturating_add(other.adapter_unmapped);
        self.adapter_clipped_bases = self
            .adapter_clipped_bases
            .saturating_add(other.adapter_clipped_bases);
        self.direct_records = self.direct_records.saturating_add(other.direct_records);
        self.traceback_records = self
            .traceback_records
            .saturating_add(other.traceback_records);
        self.mapping_worker_ns = self
            .mapping_worker_ns
            .saturating_add(other.mapping_worker_ns);
        self.record_worker_ns = self.record_worker_ns.saturating_add(other.record_worker_ns);
    }

    fn observe_result(&mut self, result: SingleAlignmentResult) {
        self.reads = self.reads.saturating_add(1);
        match result.status() {
            SingleMappingStatus::Unique => self.unique = self.unique.saturating_add(1),
            SingleMappingStatus::Ambiguous => {
                self.ambiguous = self.ambiguous.saturating_add(1);
            }
            SingleMappingStatus::Unmapped => self.unmapped = self.unmapped.saturating_add(1),
        }
        self.located_rows = self.located_rows.saturating_add(result.located_rows());
        self.verified_placements = self
            .verified_placements
            .saturating_add(result.verified_placements());
        if result.adapter_attempted() {
            self.adapter_attempted = self.adapter_attempted.saturating_add(1);
            match result.adapter_status() {
                Some(SingleMappingStatus::Unique) => {
                    self.adapter_unique = self.adapter_unique.saturating_add(1);
                }
                Some(SingleMappingStatus::Ambiguous) => {
                    self.adapter_ambiguous = self.adapter_ambiguous.saturating_add(1);
                }
                Some(SingleMappingStatus::Unmapped) | None => {
                    self.adapter_unmapped = self.adapter_unmapped.saturating_add(1);
                }
            }
            self.adapter_clipped_bases = self
                .adapter_clipped_bases
                .saturating_add(u64::try_from(result.adapter_clipped_bases()).unwrap_or(u64::MAX));
        }
    }
}

struct BatchOutput {
    records: AlignmentRecordBatch,
    observation: Observation,
}

struct RunCompletion {
    report: RunReport,
    observation: Observation,
    records_written: u64,
}

#[derive(Clone, Copy)]
enum RecordPath {
    Omitted,
    Unmapped,
    Direct,
    Traceback,
}

/// Validated inputs for canonical single-end alignment.
#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct SingleEndCommandOptions {
    pub(crate) index: PathBuf,
    pub(crate) read1: PathBuf,
    pub(crate) output: PathBuf,
    pub(crate) search_mode: SingleSearchMode,
    pub(crate) library_profile: LibraryProfile,
    pub(crate) max_edit_distance: u8,
    pub(crate) output_policy: AlignmentOutputPolicy,
    pub(crate) batch_records: u64,
    pub(crate) tie_break_seed: u64,
    pub(crate) queue_batches: usize,
    pub(crate) threads: u64,
    pub(crate) compression_threads: u32,
    pub(crate) compression_level: Option<u8>,
    pub(crate) output_contract: AlignmentAuxiliaryMode,
    pub(crate) mapped_only: bool,
    pub(crate) configuration: bsbit_cpu::Configuration,
    pub(crate) emit_metrics: bool,
}

pub(crate) fn run_single_end(
    options: &SingleEndCommandOptions,
    output_file: File,
    progress: &mut ProgressLog<'_>,
) -> Result<RunReport, CliError> {
    let process_started = options.emit_metrics.then(Instant::now);
    progress.phase(
        "index-load",
        format_args!("path={}", options.index.display()),
    );
    let reference_started = options.emit_metrics.then(Instant::now);
    let reference = load_index(options)?;
    let reference_load_ns = elapsed_ns(reference_started);
    let reference_metrics = reference.metrics();
    progress.phase(
        "index-loaded",
        format_args!(
            "contigs={} bases={}",
            reference_metrics.contig_count(),
            reference_metrics.total_reference_bases(),
        ),
    );
    let alignment_mode = match options.library_profile {
        LibraryProfile::Directional => BsbitAlignmentMode::DirectionalSingleEnd,
        LibraryProfile::NonDirectional => BsbitAlignmentMode::NonDirectionalSingleEnd,
    };
    let completion = run_single_align(options, output_file, &reference, alignment_mode, progress)?;
    write_single_metrics(
        options,
        &completion.observation,
        completion.records_written,
        reference_load_ns,
        elapsed_ns(process_started),
        progress,
    )?;
    progress.complete(
        completion.observation.reads,
        "reads",
        format_args!(
            "unique={} ambiguous={} unmapped={} bam_records={} output={}",
            completion.observation.unique,
            completion.observation.ambiguous,
            completion.observation.unmapped,
            completion.records_written,
            options.output.display(),
        ),
    );
    Ok(completion.report)
}

fn run_single_align(
    options: &SingleEndCommandOptions,
    output_file: File,
    reference: &ReferenceIndex,
    alignment_mode: BsbitAlignmentMode,
    progress: &mut ProgressLog<'_>,
) -> Result<RunCompletion, CliError> {
    if options.threads == 1 {
        run_single_align_scalar(options, output_file, reference, alignment_mode, progress)
    } else {
        run_single_align_parallel(options, output_file, reference, alignment_mode, progress)
    }
}

fn run_single_align_scalar(
    options: &SingleEndCommandOptions,
    output_file: File,
    reference: &ReferenceIndex,
    alignment_mode: BsbitAlignmentMode,
    progress: &mut ProgressLog<'_>,
) -> Result<RunCompletion, CliError> {
    let mut reader = DecodedFastqReader::open(&options.read1, alignment_fastq_limits())
        .map_err(|error| operation_error("align", "open reads", &options.read1, &error))?;
    let record_limits = AlignmentRecordLimits::default();
    let header = build_align_header(reference, alignment_mode, record_limits).map_err(|error| {
        operation_error("align", "validate output header", &options.index, &error)
    })?;
    let mut writer = OutputWriter::create(
        output_file,
        &options.output,
        &header,
        record_limits,
        options.compression_threads,
        options.compression_level,
    )?;
    let batch_size = physical_batch_size(options.batch_records)?;
    let mut aligner = SingleBatchAligner::with_capacity(SINGLE_ALIGNMENT_BATCH_SIZE)
        .with_output_policy(options.output_policy.clone());
    let mut observation = Observation::default();
    loop {
        match reader.next_batch(batch_size) {
            Ok(batch) if batch.is_empty() => break,
            Ok(batch) => {
                observation.merge(map_and_write_single_batch(
                    reference,
                    &batch,
                    options,
                    &mut aligner,
                    &mut writer,
                )?);
                progress.progress(
                    observation.reads,
                    "reads",
                    format_args!(
                        "unique={} ambiguous={} unmapped={}",
                        observation.unique, observation.ambiguous, observation.unmapped,
                    ),
                );
            }
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
    reader
        .close()
        .map_err(|error| operation_error("align", "close reads", &options.read1, &error))?;
    let (report, records_written) = writer.finish(&options.output)?;
    Ok(RunCompletion {
        report,
        observation,
        records_written,
    })
}

struct PreparedSingleInput {
    reader: DecodedFastqReader,
    batch_size: usize,
}

struct SingleSink {
    writer: OutputWriter,
    observation: Observation,
}

fn run_single_align_parallel(
    options: &SingleEndCommandOptions,
    output_file: File,
    reference: &ReferenceIndex,
    alignment_mode: BsbitAlignmentMode,
    progress: &mut ProgressLog<'_>,
) -> Result<RunCompletion, CliError> {
    let workers = physical_thread_count(options.threads)?;
    let aligners = (0..workers)
        .map(|_| {
            Mutex::new(
                SingleBatchAligner::with_capacity(SINGLE_ALIGNMENT_BATCH_SIZE)
                    .with_output_policy(options.output_policy.clone()),
            )
        })
        .collect::<Vec<_>>();
    let record_limits = AlignmentRecordLimits::default();
    run_ordered_parallel(
        workers,
        options.queue_batches,
        || {
            let reader = DecodedFastqReader::open(&options.read1, alignment_fastq_limits())
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
                build_align_header(reference, alignment_mode, record_limits).map_err(|error| {
                    operation_error("align", "validate output header", &options.index, &error)
                })?;
            let writer = OutputWriter::create(
                output_file,
                &options.output,
                &header,
                record_limits,
                options.compression_threads,
                options.compression_level,
            )?;
            Ok(SingleSink {
                writer,
                observation: Observation::default(),
            })
        },
        |sink, output| {
            write_parallel_records(sink, output)?;
            progress.progress(
                sink.observation.reads,
                "reads",
                format_args!(
                    "unique={} ambiguous={} unmapped={}",
                    sink.observation.unique, sink.observation.ambiguous, sink.observation.unmapped,
                ),
            );
            Ok(())
        },
        |sink| {
            let (report, records_written) = sink.writer.finish(&options.output)?;
            Ok(RunCompletion {
                report,
                observation: sink.observation,
                records_written,
            })
        },
    )
}

fn produce_single_batches(
    prepared: PreparedSingleInput,
    dispatcher: &mut WorkDispatcher<FastqRecordBatch>,
    cancellation: &AtomicBool,
    options: &SingleEndCommandOptions,
) -> ProducerOutcome {
    let PreparedSingleInput {
        mut reader,
        batch_size,
    } = prepared;
    loop {
        if cancellation.load(Ordering::Relaxed) || dispatcher.is_cancelled() {
            let _ = reader.close();
            return ProducerOutcome::Cancelled;
        }
        match reader.next_batch(batch_size) {
            Ok(batch) if batch.is_empty() => break,
            Ok(batch) => {
                if let Err(error) = dispatcher.send(batch) {
                    let _ = reader.close();
                    return dispatch_failure(error);
                }
            }
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
    input: &FastqRecordBatch,
    options: &SingleEndCommandOptions,
    cancellation: &AtomicBool,
    aligner: &mut SingleBatchAligner,
) -> WorkerOutcome<BatchOutput> {
    match map_single_records(reference, input, options, Some(cancellation), aligner) {
        Ok(records) => WorkerOutcome::Completed(records),
        Err(_error) if cancellation.load(Ordering::Relaxed) => WorkerOutcome::Cancelled,
        Err(error) => WorkerOutcome::Failed(error),
    }
}

fn write_parallel_records(sink: &mut SingleSink, output: BatchOutput) -> Result<(), CliError> {
    let written = sink.writer.write_batch(&output.records);
    sink.observation.merge(output.observation);
    drop(output.records);
    written
}

fn map_and_write_single_batch(
    reference: &ReferenceIndex,
    input: &FastqRecordBatch,
    options: &SingleEndCommandOptions,
    aligner: &mut SingleBatchAligner,
    writer: &mut OutputWriter,
) -> Result<Observation, CliError> {
    let output = map_single_records(reference, input, options, None, aligner)?;
    writer.write_batch(&output.records)?;
    Ok(output.observation)
}

#[allow(clippy::too_many_lines)]
fn map_single_records(
    reference: &ReferenceIndex,
    input: &FastqRecordBatch,
    options: &SingleEndCommandOptions,
    cancellation: Option<&AtomicBool>,
    aligner: &mut SingleBatchAligner,
) -> Result<BatchOutput, CliError> {
    let maximum_edit_distance = options.max_edit_distance;
    let limits = AlignmentRecordLimits::default();
    let mut output = AlignmentRecordBatch::new();
    let mut composer = AlignmentRecordComposer::new();
    let mut reads = Vec::with_capacity(SINGLE_ALIGNMENT_BATCH_SIZE);
    let mut read_keys = Vec::with_capacity(SINGLE_ALIGNMENT_BATCH_SIZE);
    let mut observation = Observation::default();
    let mut expected_records = 0_usize;
    for chunk_start in (0..input.len()).step_by(SINGLE_ALIGNMENT_BATCH_SIZE) {
        let chunk_end = input
            .len()
            .min(chunk_start.saturating_add(SINGLE_ALIGNMENT_BATCH_SIZE));
        if cancellation.is_some_and(|flag| flag.load(Ordering::Relaxed)) {
            return Err(CliError::operation("align: single mapping cancelled"));
        }
        reads.clear();
        read_keys.clear();
        reads.extend((chunk_start..chunk_end).map(|index| {
            input
                .get(index)
                .expect("single FASTQ chunk index is bounded")
                .sequence()
        }));
        read_keys.extend((chunk_start..chunk_end).map(|index| {
            input
                .get(index)
                .expect("single FASTQ chunk index is bounded")
                .ordinal()
                .get()
        }));
        let mapping_started = options.emit_metrics.then(Instant::now);
        let mapped = aligner
            .map_reads_for_output_with_tie_break_keys(
                reference,
                &reads,
                maximum_edit_distance,
                options.library_profile,
                options.search_mode,
                options.tie_break_seed,
                &read_keys,
            )
            .map_err(|error| CliError::operation(format!("align: map single batch: {error}")))?;
        observation.mapping_worker_ns = observation
            .mapping_worker_ns
            .saturating_add(elapsed_ns(mapping_started));
        for result in mapped.iter().copied() {
            observation.observe_result(result);
        }
        let record_started = options.emit_metrics.then(Instant::now);
        for (offset, result) in mapped.iter().copied().enumerate() {
            if cancellation.is_some_and(|flag| flag.load(Ordering::Relaxed)) {
                return Err(CliError::operation("align: single mapping cancelled"));
            }
            let source = input
                .get(chunk_start + offset)
                .expect("mapped single result index is bounded");
            let path = materialize_single_record(
                reference,
                source,
                result,
                options,
                &mut composer,
                limits,
            )?;
            if !matches!(path, RecordPath::Omitted) {
                expected_records = expected_records.checked_add(1).ok_or_else(|| {
                    CliError::operation("align: single output record count overflow")
                })?;
            }
            if options.emit_metrics {
                match path {
                    RecordPath::Direct => {
                        observation.direct_records = observation.direct_records.saturating_add(1);
                    }
                    RecordPath::Traceback => {
                        observation.traceback_records =
                            observation.traceback_records.saturating_add(1);
                    }
                    RecordPath::Omitted | RecordPath::Unmapped => {}
                }
            }
        }
        composer.flush_into(&mut output, limits).map_err(|error| {
            CliError::operation(format!("align: store single alignment batch: {error}"))
        })?;
        observation.record_worker_ns = observation
            .record_worker_ns
            .saturating_add(elapsed_ns(record_started));
    }
    if output.len() != expected_records {
        return Err(CliError::operation(format!(
            "align: single record cardinality mismatch: expected {expected_records} output {}",
            output.len()
        )));
    }
    Ok(BatchOutput {
        records: output,
        observation,
    })
}

fn materialize_single_record(
    reference: &ReferenceIndex,
    source: BorrowedFastqRecord<'_>,
    result: SingleAlignmentResult,
    options: &SingleEndCommandOptions,
    composer: &mut AlignmentRecordComposer,
    limits: AlignmentRecordLimits,
) -> Result<RecordPath, CliError> {
    let full_read = BorrowedAlignmentRead::new(source.sequence(), source.quality());
    let Some(placement) = result.placement() else {
        return push_or_omit_unmapped_single(source, full_read, options, composer, limits);
    };
    let retained_range = result.retained_query_interval();
    let (contig_id, interval) = resolve_single_reference_interval(
        reference,
        placement.contig_ordinal(),
        placement.start(),
        placement.end(),
    )
    .map_err(|error| {
        CliError::operation(format!(
            "align: materialize record {} from {}: {error}",
            source.ordinal().get(),
            options.read1.display()
        ))
    })?;
    let mapping_quality = match result.status() {
        SingleMappingStatus::Unique => result.mapping_quality(),
        SingleMappingStatus::Ambiguous | SingleMappingStatus::Unmapped => 0,
    };
    let direct_placement = AlignmentPlacement::new(
        placement.contig_ordinal(),
        interval,
        placement.strand(),
        placement.distance(),
    );
    let pushed = if matches!(options.output_contract, AlignmentAuxiliaryMode::Minimal) {
        try_push_direct_single(
            composer,
            reference,
            source,
            full_read,
            retained_range.clone(),
            direct_placement,
            limits,
            mapping_quality,
        )
        .map_err(|error| {
            CliError::operation(format!(
                "align: construct record {} from {}: {error}",
                source.ordinal().get(),
                options.read1.display()
            ))
        })?
    } else {
        false
    };
    if pushed {
        return Ok(RecordPath::Direct);
    }

    let retained_sequence =
        NormalizedSequence::from_bases(source.sequence()[retained_range.clone()].iter().copied());
    let alignment = traceback_read_placement(
        reference,
        &retained_sequence,
        &contig_id,
        interval,
        placement.strand(),
        placement.distance(),
    )
    .map_err(|error| {
        CliError::operation(format!(
            "align: materialize record {} from {}: {error}",
            source.ordinal().get(),
            options.read1.display()
        ))
    })?;
    composer
        .push_retained_single_with_mapping_quality(
            reference,
            source.name(),
            full_read,
            retained_range,
            &retained_sequence,
            &alignment,
            limits,
            options.output_contract,
            mapping_quality,
        )
        .map_err(|error| {
            CliError::operation(format!(
                "align: construct record {} from {}: {error}",
                source.ordinal().get(),
                options.read1.display()
            ))
        })?;
    Ok(RecordPath::Traceback)
}

fn push_or_omit_unmapped_single(
    source: BorrowedFastqRecord<'_>,
    full_read: BorrowedAlignmentRead<'_>,
    options: &SingleEndCommandOptions,
    composer: &mut AlignmentRecordComposer,
    limits: AlignmentRecordLimits,
) -> Result<RecordPath, CliError> {
    if options.mapped_only {
        return Ok(RecordPath::Omitted);
    }
    composer
        .push_unmapped_single(source.name(), full_read, limits)
        .map_err(|error| {
            CliError::operation(format!(
                "align: construct unmapped record {} from {}: {error}",
                source.ordinal().get(),
                options.read1.display()
            ))
        })?;
    Ok(RecordPath::Unmapped)
}

#[allow(clippy::too_many_arguments)]
fn try_push_direct_single(
    composer: &mut AlignmentRecordComposer,
    reference: &ReferenceIndex,
    source: BorrowedFastqRecord<'_>,
    full_read: BorrowedAlignmentRead<'_>,
    retained_range: core::ops::Range<usize>,
    placement: AlignmentPlacement,
    limits: AlignmentRecordLimits,
    mapping_quality: u8,
) -> Result<bool, RecordBuildError> {
    if retained_range.start == 0 && retained_range.end == source.sequence().len() {
        composer.try_push_ungapped_single(
            reference,
            source.name(),
            full_read,
            placement,
            limits,
            mapping_quality,
        )
    } else {
        composer.try_push_soft_clipped_ungapped_single(
            reference,
            source.name(),
            full_read,
            retained_range,
            placement,
            limits,
            mapping_quality,
        )
    }
}

fn resolve_single_reference_interval(
    reference: &ReferenceIndex,
    contig_ordinal: u64,
    start: u64,
    end: u64,
) -> Result<(ContigId, ReferenceInterval), String> {
    let contig_id = reference
        .contig_id(contig_ordinal)
        .map_err(|error| error.to_string())?;
    let contig = reference
        .resolve_contig(&contig_id)
        .map_err(|error| error.to_string())?;
    let interval =
        ReferenceInterval::new(start, end, ReferenceLength::new(contig.sequence().len()))
            .map_err(|error| error.to_string())?;
    Ok((contig_id, interval))
}

fn load_index(options: &SingleEndCommandOptions) -> Result<ReferenceIndex, CliError> {
    let internal_prefix = internal_search_file_prefix(&options.index);
    let threads = physical_thread_count(options.threads)?;
    let loaded =
        load_combined_reference_catalog(&options.index, None, &internal_prefix, threads)
            .map_err(|error| operation_error("align", "validate index", &options.index, &error))?;
    Ok(loaded.into_index())
}

fn build_align_header(
    reference: &ReferenceIndex,
    alignment_mode: BsbitAlignmentMode,
    limits: AlignmentRecordLimits,
) -> Result<SamHeader, RecordBuildError> {
    build_sam_header(reference, limits)?
        .with_bsbit_metadata(BsbitHeaderMetadata::new(alignment_mode), limits)
        .map_err(Into::into)
}

struct OutputWriter {
    writer: BamStagingWriter,
    expected_records: u64,
}

impl OutputWriter {
    fn create(
        output_file: File,
        target: &Path,
        header: &SamHeader,
        limits: AlignmentRecordLimits,
        compression_threads: u32,
        compression_level: Option<u8>,
    ) -> Result<Self, CliError> {
        let writer = BamStagingWriter::create_direct_from_file(
            target,
            output_file,
            header,
            limits,
            compression_threads,
            compression_level,
        );
        writer
            .map(|writer| Self {
                writer,
                expected_records: 0,
            })
            .map_err(|error| operation_error("align", "open BAM output", target, &error))
    }

    fn write_batch(&mut self, batch: &AlignmentRecordBatch) -> Result<(), CliError> {
        for record in batch.records() {
            self.writer
                .write_borrowed_alignment_record(&record)
                .map_err(|error| {
                    CliError::operation(format!("align: write BAM record: {error}"))
                })?;
        }
        self.expected_records = self
            .expected_records
            .checked_add(u64::try_from(batch.len()).map_err(|_| {
                CliError::operation("align: BAM batch cardinality is not representable")
            })?)
            .ok_or_else(|| CliError::operation("align: expected BAM record count overflow"))?;
        Ok(())
    }

    fn finish(self, target: &Path) -> Result<(RunReport, u64), CliError> {
        let actual_records = self.writer.records_written();
        if actual_records != self.expected_records {
            return Err(CliError::operation(format!(
                "align: BAM writer cardinality mismatch: expected {} wrote {actual_records}",
                self.expected_records
            )));
        }
        let records_written = self
            .writer
            .finish_direct()
            .map_err(|error| operation_error("align", "finalize BAM output", target, &error))?;
        Ok((RunReport::default(), records_written))
    }
}

fn elapsed_ns(started: Option<Instant>) -> u128 {
    started.map_or(0, |started| started.elapsed().as_nanos())
}

// The body is a declarative schema table: keeping every field name adjacent
// to its value is safer than splitting the row across parallel helpers.
#[allow(clippy::too_many_lines)]
fn write_single_metrics(
    options: &SingleEndCommandOptions,
    observation: &Observation,
    records_written: u64,
    reference_load_ns: u128,
    process_total_ns: u128,
    progress: &mut ProgressLog<'_>,
) -> Result<(), CliError> {
    if !options.emit_metrics {
        return Ok(());
    }
    let search_mode = match options.search_mode {
        SingleSearchMode::Default => "default",
        SingleSearchMode::Sensitive => "sensitive",
    };
    let row = MetricsRow::new([
        ("schema", SINGLE_END_METRICS_SCHEMA.to_owned()),
        ("reads", observation.reads.to_string()),
        ("unique", observation.unique.to_string()),
        ("ambiguous", observation.ambiguous.to_string()),
        ("unmapped", observation.unmapped.to_string()),
        ("bam_records", records_written.to_string()),
        ("mapping_threads", options.threads.to_string()),
        (
            "compression_threads",
            options.compression_threads.to_string(),
        ),
        (
            "architecture",
            options.configuration.features().architecture().to_string(),
        ),
        ("backend", options.configuration.backend().to_string()),
        (
            "instruction_set",
            options.configuration.backend().instruction_set().to_owned(),
        ),
        (
            "output_contract",
            output_contract_name(options.output_contract).to_owned(),
        ),
        (
            "library_profile",
            library_profile_name(options.library_profile).to_owned(),
        ),
        ("search_mode", search_mode.to_owned()),
        ("alignment_policy", ALIGNMENT_POLICY_ID.to_owned()),
        ("mapq_policy", MAPQ_POLICY_ID.to_owned()),
        ("max_edit_distance", options.max_edit_distance.to_string()),
        ("adapter", adapter_name(&options.output_policy)),
        (
            "adapter_min_overlap",
            options.output_policy.adapter_minimum_overlap().to_string(),
        ),
        (
            "adapter_max_clip",
            options
                .output_policy
                .adapter_maximum_clip_bases()
                .to_string(),
        ),
        (
            "soft_clip_mode",
            soft_clip_mode_name(options.output_policy.soft_clip_mode()).to_owned(),
        ),
        (
            "max_soft_clip",
            options.output_policy.maximum_soft_clip_bases().to_string(),
        ),
        ("reference_load_ns", reference_load_ns.to_string()),
        ("process_total_ns", process_total_ns.to_string()),
        (
            "mapping_worker_total_ns",
            observation.mapping_worker_ns.to_string(),
        ),
        (
            "record_worker_total_ns",
            observation.record_worker_ns.to_string(),
        ),
        ("batch_size", options.batch_records.to_string()),
        ("queue_batches", options.queue_batches.to_string()),
        ("located_rows", observation.located_rows.to_string()),
        (
            "verified_placements",
            observation.verified_placements.to_string(),
        ),
        (
            "adapter_attempted_reads",
            observation.adapter_attempted.to_string(),
        ),
        (
            "adapter_unique_reads",
            observation.adapter_unique.to_string(),
        ),
        (
            "adapter_ambiguous_reads",
            observation.adapter_ambiguous.to_string(),
        ),
        (
            "adapter_unmapped_reads",
            observation.adapter_unmapped.to_string(),
        ),
        (
            "adapter_clipped_bases",
            observation.adapter_clipped_bases.to_string(),
        ),
        (
            "direct_ungapped_records",
            observation.direct_records.to_string(),
        ),
        (
            "traceback_records",
            observation.traceback_records.to_string(),
        ),
        (
            "read_output",
            if options.mapped_only {
                "mapped-only"
            } else {
                "complete"
            }
            .to_owned(),
        ),
        ("tie_break_seed", options.tie_break_seed.to_string()),
    ]);
    progress
        .machine_line(format_args!("{}", row.header()))
        .map_err(|error| CliError::operation(format!("align: write single metrics: {error}")))?;
    progress
        .machine_line(format_args!("{}", row.values()))
        .map_err(|error| CliError::operation(format!("align: write single metrics: {error}")))
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
