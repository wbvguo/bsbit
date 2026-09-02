//! FASTQ and BAM orchestration for canonical single-end alignment.

use std::path::{Path, PathBuf};
use std::sync::Mutex;
use std::sync::atomic::{AtomicBool, Ordering};
use std::time::Instant;

use bsbit_align::library::LibraryProfile;
use bsbit_align::materialize::traceback_read_placement;
use bsbit_align::single_end::{
    SingleAlignmentResult, SingleBatchAligner, SingleMappingStatus, SingleSearchMode,
};
use bsbit_core::coordinate::{ReferenceInterval, ReferenceLength};
use bsbit_core::reference::ReferenceSemanticDigest;
use bsbit_core::sequence::NormalizedSequence;
use bsbit_hts::TextRecordLimits;
use bsbit_hts::{
    AlignmentAuxiliaryMode, AlignmentPlacement, AlignmentRecordBatch, AlignmentRecordLimits,
    BamStagingWriter, BorrowedAlignmentRead, BorrowedFastqRecord, BsbitAlignmentMode,
    BsbitProgramProvenance, DecodedFastqReader, FastqRecordBatch, SAM_MAX_QUERY_NAME_BYTES,
    SamHeader,
};
use bsbit_index::reference::{ContigId, ReferenceIndex};
use bsbit_index::storage::combined::load_combined_reference_catalog;
use bsbit_io::validate_replace_target;

use super::{ReadOutputMode, library_profile_name, output_contract_name};
use crate::parallel::{
    DispatchError, ProducerOutcome, WorkDispatcher, WorkerOutcome, run_ordered_parallel,
};
use crate::record_composition::{RecordBuildError, SingleRecordComposer, build_sam_header};
use crate::{CliError, CliWarning, RunReport};

use crate::command::{
    ReadLayout, caller_compatible_alignment_mode, internal_search_file_prefix, unused_staging_path,
};

const MAX_CLI_READ_BASES: u64 = 1_000_000;
const MAX_CLI_DESCRIPTION_BYTES: u64 = 1_000_000;
const SEARCH_BATCH_SIZE: usize = 64;
const METRICS_SCHEMA: &str = "bsbit-single-alignment-metrics-v2";

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
pub(crate) struct Options {
    pub(super) index: PathBuf,
    pub(super) read1: PathBuf,
    pub(super) output_bam: PathBuf,
    pub(super) search_mode: SingleSearchMode,
    pub(super) library_profile: LibraryProfile,
    pub(super) max_edit_distance: u64,
    pub(super) batch_records: u64,
    pub(super) threads: u64,
    pub(super) bam_threads: u32,
    pub(super) bam_compression_level: Option<u8>,
    pub(super) output_contract: AlignmentAuxiliaryMode,
    pub(super) read_output: ReadOutputMode,
    pub(super) emit_metrics: bool,
}

pub(super) fn run(options: &Options) -> Result<RunReport, CliError> {
    let process_started = options.emit_metrics.then(Instant::now);
    validate_output_target(&options.output_bam)?;
    let staging = unused_staging_path(&options.output_bam, "align", "output")?;
    let reference_started = options.emit_metrics.then(Instant::now);
    let (reference, semantic_digest) = load_index(options)?;
    let reference_load_ns = elapsed_ns(reference_started);
    let alignment_mode =
        caller_compatible_alignment_mode(ReadLayout::SingleEnd, options.library_profile);
    let completion = run_single_align(
        options,
        &reference,
        semantic_digest,
        &staging,
        alignment_mode,
    )?;
    write_single_metrics(
        options,
        &completion.observation,
        completion.records_written,
        reference_load_ns,
        elapsed_ns(process_started),
    );
    Ok(completion.report)
}

fn run_single_align(
    options: &Options,
    reference: &ReferenceIndex,
    semantic_digest: ReferenceSemanticDigest,
    staging: &Path,
    alignment_mode: BsbitAlignmentMode,
) -> Result<RunCompletion, CliError> {
    if options.threads == 1 {
        run_single_align_scalar(options, reference, semantic_digest, staging, alignment_mode)
    } else {
        run_single_align_parallel(options, reference, semantic_digest, staging, alignment_mode)
    }
}

fn run_single_align_scalar(
    options: &Options,
    reference: &ReferenceIndex,
    semantic_digest: ReferenceSemanticDigest,
    staging: &Path,
    alignment_mode: BsbitAlignmentMode,
) -> Result<RunCompletion, CliError> {
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
        options.bam_threads,
        options.bam_compression_level,
    )?;
    let batch_size = physical_batch_size(options.batch_records)?;
    let mut aligner = SingleBatchAligner::with_capacity(SEARCH_BATCH_SIZE);
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
    let (report, records_written) = writer.finish(&options.output_bam)?;
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
    options: &Options,
    reference: &ReferenceIndex,
    semantic_digest: ReferenceSemanticDigest,
    staging: &Path,
    alignment_mode: BsbitAlignmentMode,
) -> Result<RunCompletion, CliError> {
    let workers = physical_thread_count(options.threads)?;
    let aligners = (0..workers)
        .map(|_| Mutex::new(SingleBatchAligner::with_capacity(SEARCH_BATCH_SIZE)))
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
            let writer = OutputWriter::create(
                &options.output_bam,
                staging,
                &header,
                record_limits,
                options.bam_threads,
                options.bam_compression_level,
            )?;
            Ok(SingleSink {
                writer,
                observation: Observation::default(),
            })
        },
        write_parallel_records,
        |sink| {
            let (report, records_written) = sink.writer.finish(&options.output_bam)?;
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
    options: &Options,
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
    options: &Options,
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
    options: &Options,
    aligner: &mut SingleBatchAligner,
    writer: &mut OutputWriter,
) -> Result<Observation, CliError> {
    let output = map_single_records(reference, input, options, None, aligner)?;
    writer.write_batch(&output.records)?;
    Ok(output.observation)
}

fn map_single_records(
    reference: &ReferenceIndex,
    input: &FastqRecordBatch,
    options: &Options,
    cancellation: Option<&AtomicBool>,
    aligner: &mut SingleBatchAligner,
) -> Result<BatchOutput, CliError> {
    let maximum_edit_distance = u8::try_from(options.max_edit_distance).map_err(|_| {
        CliError::operation(format!(
            "align: single edit distance {} is not representable",
            options.max_edit_distance
        ))
    })?;
    let limits = AlignmentRecordLimits::default();
    let mut output = AlignmentRecordBatch::new();
    let mut composer = SingleRecordComposer::new();
    let mut reads = Vec::with_capacity(SEARCH_BATCH_SIZE);
    let mut observation = Observation::default();
    let mut expected_records = 0_usize;
    for chunk_start in (0..input.len()).step_by(SEARCH_BATCH_SIZE) {
        let chunk_end = input
            .len()
            .min(chunk_start.saturating_add(SEARCH_BATCH_SIZE));
        if cancellation.is_some_and(|flag| flag.load(Ordering::Relaxed)) {
            return Err(CliError::operation("align: single mapping cancelled"));
        }
        reads.clear();
        reads.extend((chunk_start..chunk_end).map(|index| {
            input
                .get(index)
                .expect("single FASTQ chunk index is bounded")
                .sequence()
        }));
        let mapping_started = options.emit_metrics.then(Instant::now);
        let mapped = aligner
            .map_reads_for_output(
                reference,
                &reads,
                maximum_edit_distance,
                options.library_profile,
                options.search_mode,
            )
            .map_err(|error| CliError::operation(format!("align: map single batch: {error}")))?;
        observation.mapping_worker_ns = observation
            .mapping_worker_ns
            .saturating_add(elapsed_ns(mapping_started));
        if options.emit_metrics {
            for result in mapped.iter().copied() {
                observation.observe_result(result);
            }
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
    options: &Options,
    composer: &mut SingleRecordComposer,
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
    options: &Options,
    composer: &mut SingleRecordComposer,
    limits: AlignmentRecordLimits,
) -> Result<RecordPath, CliError> {
    if matches!(options.read_output, ReadOutputMode::MappedOnly) {
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
    composer: &mut SingleRecordComposer,
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

fn load_index(options: &Options) -> Result<(ReferenceIndex, ReferenceSemanticDigest), CliError> {
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

struct OutputWriter {
    writer: BamStagingWriter,
    expected_records: u64,
}

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
            .map(|writer| Self {
                writer,
                expected_records: 0,
            })
            .map_err(|error| operation_error("align", "create BAM staging", target, &error))
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
        let completed = self
            .writer
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
        Ok((
            RunReport::with_warning(warning),
            publication.records_written(),
        ))
    }
}

fn elapsed_ns(started: Option<Instant>) -> u128 {
    started.map_or(0, |started| started.elapsed().as_nanos())
}

fn write_single_metrics(
    options: &Options,
    observation: &Observation,
    records_written: u64,
    reference_load_ns: u128,
    process_total_ns: u128,
) {
    if !options.emit_metrics {
        return;
    }
    println!(
        "schema\treads\tunique\tambiguous\tunmapped\tbam_records\tmapping_threads\tbam_threads\toutput_contract\tlibrary_profile\tsearch_mode\treference_load_ns\tprocess_total_ns\tmapping_worker_total_ns\trecord_worker_total_ns\tlocated_rows\tverified_placements\tadapter_attempted_reads\tadapter_unique_reads\tadapter_ambiguous_reads\tadapter_unmapped_reads\tadapter_clipped_bases\tdirect_ungapped_records\ttraceback_records\tread_output"
    );
    let search_mode = match options.search_mode {
        SingleSearchMode::Default => "default",
        SingleSearchMode::Sensitive => "sensitive",
    };
    let fields = [
        METRICS_SCHEMA.to_owned(),
        observation.reads.to_string(),
        observation.unique.to_string(),
        observation.ambiguous.to_string(),
        observation.unmapped.to_string(),
        records_written.to_string(),
        options.threads.to_string(),
        options.bam_threads.to_string(),
        output_contract_name(options.output_contract).to_owned(),
        library_profile_name(options.library_profile).to_owned(),
        search_mode.to_owned(),
        reference_load_ns.to_string(),
        process_total_ns.to_string(),
        observation.mapping_worker_ns.to_string(),
        observation.record_worker_ns.to_string(),
        observation.located_rows.to_string(),
        observation.verified_placements.to_string(),
        observation.adapter_attempted.to_string(),
        observation.adapter_unique.to_string(),
        observation.adapter_ambiguous.to_string(),
        observation.adapter_unmapped.to_string(),
        observation.adapter_clipped_bases.to_string(),
        observation.direct_records.to_string(),
        observation.traceback_records.to_string(),
        options.read_output.name().to_owned(),
    ];
    println!("{}", fields.join("\t"));
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
