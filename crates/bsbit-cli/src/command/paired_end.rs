//! Paired-end execution pipeline for the standard alignment command.
//!
//! CLI parsing and layout dispatch remain in the parent module. This module
//! owns pair decoding, mapping orchestration, record composition, BAM writing,
//! and paired-only metrics.

use std::error::Error;
use std::fs::File;
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::mpsc::{Receiver, sync_channel};
use std::thread;

use bsbit_align::library::{LibraryProfile, TemplateSpan, TemplateSpanBounds};
use bsbit_align::materialize::traceback_read_placement;
use bsbit_align::paired_end::{
    PAIRED_ALIGNMENT_BATCH_SIZE, PairMappingStatus, PairedAlignmentOptions, PairedBatchAligner,
    PairedSearchMode,
};
use bsbit_align::{ALIGNMENT_POLICY_ID, AlignmentOutputPolicy, MAPQ_POLICY_ID};
use bsbit_core::coordinate::{ReferenceInterval, ReferenceLength};
use bsbit_core::sequence::NormalizedSequence;
use bsbit_hts::{
    AlignmentAuxiliaryMode, AlignmentPlacement, AlignmentRead, AlignmentRecordBatch,
    AlignmentRecordLimits, BamStagingWriter, BorrowedAlignmentRead, BorrowedFastqRecord,
    BsbitAlignmentMode, BsbitHeaderMetadata, DecodedFastqReader, FastqRecordBatch, SamHeader,
};
use bsbit_index::reference::ReferenceIndex;
use bsbit_index::storage::combined::load_combined_reference_catalog;

use super::align_options::{Options, ReadOutputMode, invalid};
use super::alignment_input::alignment_fastq_limits;
use super::alignment_metrics::{
    MetricsRow, MetricsTimer, adapter_name, library_profile_name, mate_rescue_name,
    output_contract_name, read_output_name, search_mode_name, soft_clip_fallback_name,
    soft_clip_mode_name, strategy_id,
};
use super::internal_search_file_prefix;
use crate::cpu_placement::CpuPlacement;
use crate::progress::ProgressLog;
use crate::record_composition::{AlignmentRecordComposer, build_sam_header};

const PAIRED_END_METRICS_SCHEMA: &str = "bsbit-alignment-metrics-paired-end-v1";

struct PairedInputBatch {
    first: FastqRecordBatch,
    second: FastqRecordBatch,
}

impl PairedInputBatch {
    fn len(&self) -> usize {
        self.first.len()
    }

    fn get(&self, index: usize) -> Option<PairedInputRecord<'_>> {
        let first = self.first.get(index)?;
        let second = self.second.get(index)?;
        Some(PairedInputRecord { first, second })
    }
}

#[derive(Clone, Copy)]
struct PairedInputRecord<'a> {
    first: BorrowedFastqRecord<'a>,
    second: BorrowedFastqRecord<'a>,
}

impl<'a> PairedInputRecord<'a> {
    fn shared_name(self) -> &'a [u8] {
        if self.first.name() == self.second.name() {
            self.first.name()
        } else {
            self.first
                .name()
                .strip_suffix(b"/1")
                .unwrap_or(self.first.name())
        }
    }
}

#[derive(Clone, Copy, Debug, Default)]
struct SoftClipObservation {
    attempted_pairs: u64,
    unique_pairs: u64,
    ambiguous_pairs: u64,
    unmapped_pairs: u64,
    clipped_mates: u64,
    clipped_bases: u64,
}

impl SoftClipObservation {
    fn observe(&mut self, class: PairMappingStatus) {
        match class {
            PairMappingStatus::Unique => {
                self.unique_pairs = self.unique_pairs.saturating_add(1);
            }
            PairMappingStatus::Ambiguous => {
                self.ambiguous_pairs = self.ambiguous_pairs.saturating_add(1);
            }
            PairMappingStatus::Unmapped => {
                self.unmapped_pairs = self.unmapped_pairs.saturating_add(1);
            }
        }
    }

    fn merge(&mut self, other: Self) {
        self.attempted_pairs = self.attempted_pairs.saturating_add(other.attempted_pairs);
        self.unique_pairs = self.unique_pairs.saturating_add(other.unique_pairs);
        self.ambiguous_pairs = self.ambiguous_pairs.saturating_add(other.ambiguous_pairs);
        self.unmapped_pairs = self.unmapped_pairs.saturating_add(other.unmapped_pairs);
        self.clipped_mates = self.clipped_mates.saturating_add(other.clipped_mates);
        self.clipped_bases = self.clipped_bases.saturating_add(other.clipped_bases);
    }
}

#[derive(Clone, Copy, Debug, Default)]
struct MateRescueObservation {
    attempted: u64,
    unique: u64,
    ambiguous: u64,
    unmapped: u64,
}

impl MateRescueObservation {
    fn observe(&mut self, class: PairMappingStatus) {
        self.attempted = self.attempted.saturating_add(1);
        match class {
            PairMappingStatus::Unique => {
                self.unique = self.unique.saturating_add(1);
            }
            PairMappingStatus::Ambiguous => {
                self.ambiguous = self.ambiguous.saturating_add(1);
            }
            PairMappingStatus::Unmapped => {
                self.unmapped = self.unmapped.saturating_add(1);
            }
        }
    }

    fn merge(&mut self, other: Self) {
        self.attempted = self.attempted.saturating_add(other.attempted);
        self.unique = self.unique.saturating_add(other.unique);
        self.ambiguous = self.ambiguous.saturating_add(other.ambiguous);
        self.unmapped = self.unmapped.saturating_add(other.unmapped);
    }
}

struct WriterObservation {
    records: u64,
    bam_write_ns: u128,
    finalize_ns: u128,
}

#[derive(Clone, Copy, Default)]
struct PairClassCounts {
    unique: u64,
    ambiguous: u64,
    unmapped: u64,
}

impl PairClassCounts {
    fn observe(&mut self, class: PairMappingStatus) {
        let count = match class {
            PairMappingStatus::Unique => &mut self.unique,
            PairMappingStatus::Ambiguous => &mut self.ambiguous,
            PairMappingStatus::Unmapped => &mut self.unmapped,
        };
        *count = count.saturating_add(1);
    }

    fn merge(&mut self, other: Self) {
        self.unique = self.unique.saturating_add(other.unique);
        self.ambiguous = self.ambiguous.saturating_add(other.ambiguous);
        self.unmapped = self.unmapped.saturating_add(other.unmapped);
    }

    const fn total(self) -> u64 {
        self.unique
            .saturating_add(self.ambiguous)
            .saturating_add(self.unmapped)
    }
}

#[derive(Default)]
struct Observation {
    classes: PairClassCounts,
    batch_processing_ns: u128,
    mapping_worker_total_ns: u128,
    record_worker_total_ns: u128,
    writer_queue_wait_ns: u128,
    writer_queue_sends: u64,
    soft_clip: SoftClipObservation,
    mate_rescue: MateRescueObservation,
}

pub(super) struct PairedRunSummary {
    observation: Observation,
    writer: WriterObservation,
    reference_load_ns: u128,
    decode_ns: u128,
}

pub(super) fn run_paired_alignment(
    options: &Options,
    output_file: File,
    progress: &mut ProgressLog<'_>,
) -> Result<PairedRunSummary, Box<dyn Error>> {
    let cpu_placement = CpuPlacement::detect(options.threads);
    let (producer, receiver) = spawn_paired_decoder(options, &cpu_placement);
    let limits = AlignmentRecordLimits::default();
    progress.phase(
        "index-load",
        format_args!("path={}", options.index.display()),
    );
    let (reference, header, reference_load_ns) = load_alignment_reference(options, limits)?;
    let reference_metrics = reference.metrics();
    progress.phase(
        "index-loaded",
        format_args!(
            "contigs={} bases={}",
            reference_metrics.contig_count(),
            reference_metrics.total_reference_bases(),
        ),
    );
    let (alignment_sender, alignment_receiver) = sync_channel(options.queue_batches);
    let writer = spawn_paired_writer(
        options,
        output_file,
        header,
        limits,
        alignment_receiver,
        &cpu_placement,
    );
    let bounds = TemplateSpanBounds::new(
        TemplateSpan::new(options.minimum_template_span),
        TemplateSpan::new(options.maximum_template_span),
    )?;
    let mut observation = Observation::default();
    let _coordinator_affinity = cpu_placement.pin_auxiliary_scoped();
    let consume_result = consume_batches(
        &reference,
        receiver,
        &alignment_sender,
        &mut observation,
        bounds,
        limits,
        options.output_contract,
        options.library_profile,
        options.search_mode.paired(),
        options.read_output,
        options.maximum_edit_distance,
        &options.output_policy,
        options.tie_break_seed,
        options.emit_metrics,
        options.threads,
        &cpu_placement,
        progress,
    );
    let producer_result = producer
        .join()
        .map_err(|_| invalid("FASTQ producer panicked"))?;
    drop(alignment_sender);
    let writer_result = writer
        .join()
        .map_err(|_| invalid("BAM writer worker panicked"))?;
    let decode_ns = producer_result.map_err(invalid)?;
    let writer_observation = writer_result.map_err(invalid)?;
    consume_result?;
    validate_paired_record_count(options, &observation, &writer_observation)?;
    let total_pairs = observation.classes.total();
    progress.complete(
        total_pairs,
        "pairs",
        format_args!(
            "reads={} unique={} ambiguous={} unmapped={} bam_records={} output={}",
            total_pairs.saturating_mul(2),
            observation.classes.unique,
            observation.classes.ambiguous,
            observation.classes.unmapped,
            writer_observation.records,
            options.output.display(),
        ),
    );
    Ok(PairedRunSummary {
        observation,
        writer: writer_observation,
        reference_load_ns,
        decode_ns,
    })
}

fn spawn_paired_decoder(
    options: &Options,
    cpu_placement: &CpuPlacement,
) -> (
    thread::JoinHandle<Result<u128, String>>,
    Receiver<PairedInputBatch>,
) {
    let (sender, receiver) = sync_channel(options.queue_batches);
    let read1 = options.read1.clone();
    let read2 = options
        .read2
        .clone()
        .expect("paired layout was validated with read 2");
    let batch_size = options.batch_size;
    let queue_batches = options.queue_batches;
    let emit_metrics = options.emit_metrics;
    let producer_cpu_placement = cpu_placement.clone();
    let producer = thread::spawn(move || {
        producer_cpu_placement.pin_auxiliary_worker();
        decode_batches(
            &read1,
            &read2,
            batch_size,
            queue_batches,
            &sender,
            emit_metrics,
        )
    });
    (producer, receiver)
}

fn spawn_paired_writer(
    options: &Options,
    output_file: File,
    header: SamHeader,
    limits: AlignmentRecordLimits,
    receiver: Receiver<Vec<AlignmentRecordBatch>>,
    cpu_placement: &CpuPlacement,
) -> thread::JoinHandle<Result<WriterObservation, String>> {
    let output = options.output.clone();
    let compression_threads = options.compression_threads;
    let compression_level = options.compression_level;
    let emit_metrics = options.emit_metrics;
    let writer_cpu_placement = cpu_placement.clone();
    thread::spawn(move || {
        writer_cpu_placement.pin_auxiliary_worker();
        write_batches(
            &output,
            output_file,
            &header,
            limits,
            compression_threads,
            compression_level,
            receiver,
            emit_metrics,
        )
    })
}

fn validate_paired_record_count(
    options: &Options,
    observation: &Observation,
    writer: &WriterObservation,
) -> Result<(), io::Error> {
    if !matches!(options.read_output, ReadOutputMode::Complete) {
        return Ok(());
    }
    let input_pairs = observation.classes.total();
    let expected_records = input_pairs
        .checked_mul(2)
        .ok_or_else(|| invalid("input primary-record count overflow"))?;
    if writer.records != expected_records {
        return Err(invalid(format!(
            "read-complete BAM wrote {} primary records for {input_pairs} input pairs; expected {expected_records}",
            writer.records
        )));
    }
    Ok(())
}

fn load_alignment_reference(
    options: &Options,
    limits: AlignmentRecordLimits,
) -> Result<(ReferenceIndex, SamHeader, u128), Box<dyn Error>> {
    let started = MetricsTimer::start(options.emit_metrics);
    let internal_prefix = internal_search_file_prefix(&options.index);
    let loaded =
        load_combined_reference_catalog(&options.index, None, &internal_prefix, options.threads)?;
    let reference = loaded.into_index();
    let alignment_mode = match options.library_profile {
        LibraryProfile::Directional => BsbitAlignmentMode::DirectionalPairedEnd,
        LibraryProfile::NonDirectional => BsbitAlignmentMode::NonDirectionalPairedEnd,
    };
    let header = build_sam_header(&reference, limits)?
        .with_bsbit_metadata(BsbitHeaderMetadata::new(alignment_mode), limits)?;
    Ok((reference, header, started.elapsed_ns()))
}

// The body is a declarative schema table: keeping every field name adjacent
// to its value is safer than splitting the row across parallel helpers.
#[allow(clippy::too_many_lines)]
pub(super) fn write_metrics(
    output: &mut impl Write,
    options: &Options,
    configuration: bsbit_cpu::Configuration,
    summary: &PairedRunSummary,
    process_total_ns: u128,
) -> io::Result<()> {
    if !options.emit_metrics {
        return Ok(());
    }
    let observation = &summary.observation;
    let writer_observation = &summary.writer;
    let row = MetricsRow::new([
        ("schema", PAIRED_END_METRICS_SCHEMA.to_owned()),
        ("pairs", observation.classes.total().to_string()),
        ("unique", observation.classes.unique.to_string()),
        ("ambiguous", observation.classes.ambiguous.to_string()),
        ("unmapped", observation.classes.unmapped.to_string()),
        ("bam_records", writer_observation.records.to_string()),
        ("mapping_threads", options.threads.to_string()),
        (
            "compression_threads",
            options.compression_threads.to_string(),
        ),
        (
            "architecture",
            configuration.features().architecture().to_string(),
        ),
        ("backend", configuration.backend().to_string()),
        (
            "instruction_set",
            configuration.backend().instruction_set().to_owned(),
        ),
        (
            "output_contract",
            output_contract_name(options.output_contract).to_owned(),
        ),
        (
            "library_profile",
            library_profile_name(options.library_profile).to_owned(),
        ),
        ("reference_mode", "indexed-reference".to_owned()),
        ("reference_load_ns", summary.reference_load_ns.to_string()),
        ("fastq_decode_ns", summary.decode_ns.to_string()),
        ("bam_write_ns", writer_observation.bam_write_ns.to_string()),
        (
            "bam_finalize_ns",
            writer_observation.finalize_ns.to_string(),
        ),
        ("process_total_ns", process_total_ns.to_string()),
        (
            "batch_processing_ns",
            observation.batch_processing_ns.to_string(),
        ),
        (
            "mapping_worker_total_ns",
            observation.mapping_worker_total_ns.to_string(),
        ),
        (
            "record_worker_total_ns",
            observation.record_worker_total_ns.to_string(),
        ),
        ("batch_size", options.batch_size.to_string()),
        ("queue_batches", options.queue_batches.to_string()),
        (
            "writer_queue_wait_ns",
            observation.writer_queue_wait_ns.to_string(),
        ),
        (
            "writer_queue_sends",
            observation.writer_queue_sends.to_string(),
        ),
        (
            "compression_level",
            options
                .compression_level
                .map_or_else(|| "default".to_owned(), |level| level.to_string()),
        ),
        (
            "max_edit_distance",
            options.maximum_edit_distance.to_string(),
        ),
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
        (
            "soft_clip_fallback",
            soft_clip_fallback_name(options.search_mode.paired(), &options.output_policy)
                .to_owned(),
        ),
        (
            "soft_clip_attempted_pairs",
            observation.soft_clip.attempted_pairs.to_string(),
        ),
        (
            "soft_clip_unique_pairs",
            observation.soft_clip.unique_pairs.to_string(),
        ),
        (
            "soft_clip_ambiguous_pairs",
            observation.soft_clip.ambiguous_pairs.to_string(),
        ),
        (
            "soft_clip_unmapped_pairs",
            observation.soft_clip.unmapped_pairs.to_string(),
        ),
        (
            "soft_clip_clipped_mates",
            observation.soft_clip.clipped_mates.to_string(),
        ),
        (
            "soft_clip_clipped_bases",
            observation.soft_clip.clipped_bases.to_string(),
        ),
        (
            "mate_rescue",
            mate_rescue_name(options.search_mode.paired()).to_owned(),
        ),
        (
            "mate_rescue_attempted_pairs",
            observation.mate_rescue.attempted.to_string(),
        ),
        (
            "mate_rescue_unique_pairs",
            observation.mate_rescue.unique.to_string(),
        ),
        (
            "mate_rescue_ambiguous_pairs",
            observation.mate_rescue.ambiguous.to_string(),
        ),
        (
            "mate_rescue_unmapped_pairs",
            observation.mate_rescue.unmapped.to_string(),
        ),
        (
            "search_mode",
            search_mode_name(options.search_mode).to_owned(),
        ),
        ("alignment_policy", ALIGNMENT_POLICY_ID.to_owned()),
        ("mapq_policy", MAPQ_POLICY_ID.to_owned()),
        ("mapq_zero_output", "all".to_owned()),
        ("strategy_id", strategy_id(options).to_owned()),
        (
            "read_output",
            read_output_name(options.read_output).to_owned(),
        ),
        ("tie_break_seed", options.tie_break_seed.to_string()),
    ]);
    writeln!(output, "{}", row.header())?;
    writeln!(output, "{}", row.values())?;
    output.flush()
}

#[allow(clippy::too_many_arguments)]
fn consume_batches(
    reference: &ReferenceIndex,
    receiver: Receiver<PairedInputBatch>,
    alignment_sender: &std::sync::mpsc::SyncSender<Vec<AlignmentRecordBatch>>,
    observation: &mut Observation,
    bounds: TemplateSpanBounds,
    limits: AlignmentRecordLimits,
    output_contract: AlignmentAuxiliaryMode,
    library_profile: LibraryProfile,
    search_mode: PairedSearchMode,
    read_output: ReadOutputMode,
    maximum_edit_distance: u8,
    output_policy: &AlignmentOutputPolicy,
    tie_break_seed: u64,
    emit_metrics: bool,
    threads: usize,
    cpu_placement: &CpuPlacement,
    progress: &mut ProgressLog<'_>,
) -> Result<(), Box<dyn Error>> {
    for batch in receiver {
        let started = MetricsTimer::start(emit_metrics);
        let processed = process_paired_batch(
            reference,
            &batch,
            bounds,
            limits,
            output_contract,
            library_profile,
            search_mode,
            read_output,
            maximum_edit_distance,
            output_policy,
            tie_break_seed,
            emit_metrics,
            threads,
            cpu_placement,
        )?;
        observation.batch_processing_ns = observation
            .batch_processing_ns
            .saturating_add(started.elapsed_ns());
        observation.mapping_worker_total_ns = observation
            .mapping_worker_total_ns
            .saturating_add(processed.mapping_worker_ns);
        observation.record_worker_total_ns = observation
            .record_worker_total_ns
            .saturating_add(processed.record_worker_ns);
        observation.classes.merge(processed.classes);
        observation.soft_clip.merge(processed.soft_clip);
        observation.mate_rescue.merge(processed.mate_rescue);
        if processed.records.iter().any(|chunk| !chunk.is_empty()) {
            let send_started = MetricsTimer::start(emit_metrics);
            alignment_sender
                .send(processed.records)
                .map_err(|_| invalid("BAM writer ended before mapping"))?;
            observation.writer_queue_wait_ns = observation
                .writer_queue_wait_ns
                .saturating_add(send_started.elapsed_ns());
            observation.writer_queue_sends = observation.writer_queue_sends.saturating_add(1);
        }
        let total_pairs = observation.classes.total();
        progress.progress(
            total_pairs,
            "pairs",
            format_args!(
                "reads={} unique={} ambiguous={} unmapped={}",
                total_pairs.saturating_mul(2),
                observation.classes.unique,
                observation.classes.ambiguous,
                observation.classes.unmapped,
            ),
        );
    }
    Ok(())
}

struct PairedBatchOutput {
    records: Vec<AlignmentRecordBatch>,
    classes: PairClassCounts,
    soft_clip: SoftClipObservation,
    mate_rescue: MateRescueObservation,
    mapping_worker_ns: u128,
    record_worker_ns: u128,
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn process_paired_batch(
    reference: &ReferenceIndex,
    records: &PairedInputBatch,
    bounds: TemplateSpanBounds,
    limits: AlignmentRecordLimits,
    output_contract: AlignmentAuxiliaryMode,
    library_profile: LibraryProfile,
    search_mode: PairedSearchMode,
    read_output: ReadOutputMode,
    maximum_edit_distance: u8,
    output_policy: &AlignmentOutputPolicy,
    tie_break_seed: u64,
    emit_metrics: bool,
    threads: usize,
    cpu_placement: &CpuPlacement,
) -> Result<PairedBatchOutput, Box<dyn Error>> {
    let next = AtomicUsize::new(0);
    let workers = threads.min(records.len().max(1));
    let mut chunks = thread::scope(|scope| {
        let mut handles = Vec::with_capacity(workers);
        for worker_ordinal in 0..workers {
            let next = &next;
            handles.push(scope.spawn(move || -> Result<_, String> {
                cpu_placement.pin_mapping_worker(worker_ordinal);
                let mut owned = Vec::new();
                let mut mapping_worker_ns = 0_u128;
                let mut record_worker_ns = 0_u128;
                let mut aligner = PairedBatchAligner::with_capacity(PAIRED_ALIGNMENT_BATCH_SIZE)
                    .with_output_policy(output_policy.clone());
                let mut composer = AlignmentRecordComposer::new();
                let mut pair_reads = Vec::with_capacity(PAIRED_ALIGNMENT_BATCH_SIZE);
                let mut read_keys = Vec::with_capacity(PAIRED_ALIGNMENT_BATCH_SIZE);
                loop {
                    let start = next.fetch_add(PAIRED_ALIGNMENT_BATCH_SIZE, Ordering::Relaxed);
                    if start >= records.len() {
                        break;
                    }
                    let end = start
                        .saturating_add(PAIRED_ALIGNMENT_BATCH_SIZE)
                        .min(records.len());
                    {
                        let mut classes = PairClassCounts::default();
                        let mut mate_rescue_observation = MateRescueObservation::default();
                        let mut processed = AlignmentRecordBatch::new();
                        pair_reads.clear();
                        read_keys.clear();
                        for index in start..end {
                            let pair = records.get(index).ok_or_else(|| {
                                String::from("paired FASTQ batch index is absent")
                            })?;
                            pair_reads.push([pair.first.sequence(), pair.second.sequence()]);
                            read_keys.push(pair.first.ordinal().get());
                        }
                        let mapping_started = MetricsTimer::start(emit_metrics);
                        let mapped = aligner
                            .map_pairs_for_output_with_tie_break_keys(
                                reference,
                                &pair_reads,
                                PairedAlignmentOptions::primary(
                                    library_profile,
                                    search_mode,
                                    bounds.minimum().get(),
                                    bounds.maximum().get(),
                                )
                                .with_maximum_edit_distance(maximum_edit_distance)
                                .map_err(|error| error.to_string())?,
                                tie_break_seed,
                                &read_keys,
                            )
                            .map_err(|error| error.to_string())?;
                        mapping_worker_ns =
                            mapping_worker_ns.saturating_add(mapping_started.elapsed_ns());
                        let mut soft_clip_observation = SoftClipObservation::default();
                        for (offset, output) in mapped.into_iter().enumerate() {
                            let index = start + offset;
                            let pair = records.get(index).ok_or_else(|| {
                                String::from("paired FASTQ batch index is absent")
                            })?;
                            if output.adapter_attempted() {
                                soft_clip_observation.attempted_pairs =
                                    soft_clip_observation.attempted_pairs.saturating_add(1);
                                soft_clip_observation
                                    .observe(output.adapter_class().unwrap_or(output.class()));
                                soft_clip_observation.clipped_mates = soft_clip_observation
                                    .clipped_mates
                                    .saturating_add(u64::from(output.adapter_clipped_mates()));
                                soft_clip_observation.clipped_bases =
                                    soft_clip_observation.clipped_bases.saturating_add(
                                        u64::try_from(output.adapter_clipped_bases())
                                            .unwrap_or(u64::MAX),
                                    );
                            }
                            if output.semi_global_attempted() {
                                soft_clip_observation.attempted_pairs =
                                    soft_clip_observation.attempted_pairs.saturating_add(1);
                                soft_clip_observation.observe(output.class());
                                soft_clip_observation.clipped_mates = soft_clip_observation
                                    .clipped_mates
                                    .saturating_add(u64::from(output.semi_global_clipped_mates()));
                                soft_clip_observation.clipped_bases =
                                    soft_clip_observation.clipped_bases.saturating_add(
                                        u64::try_from(output.semi_global_clipped_bases())
                                            .unwrap_or(u64::MAX),
                                    );
                            }
                            let class = output.class();
                            classes.observe(class);
                            if output.mate_rescue_attempted() {
                                mate_rescue_observation.observe(class);
                            }
                            let Some(selected) = output.placement() else {
                                if matches!(read_output, ReadOutputMode::Complete) {
                                    let record_started = MetricsTimer::start(emit_metrics);
                                    composer
                                        .push_unmapped_pair(
                                            pair.shared_name(),
                                            BorrowedAlignmentRead::new(
                                                pair.first.sequence(),
                                                pair.first.quality(),
                                            ),
                                            BorrowedAlignmentRead::new(
                                                pair.second.sequence(),
                                                pair.second.quality(),
                                            ),
                                            limits,
                                        )
                                        .map_err(|error| error.to_string())?;
                                    record_worker_ns = record_worker_ns
                                        .saturating_add(record_started.elapsed_ns());
                                }
                                continue;
                            };
                            let record_started = MetricsTimer::start(emit_metrics);
                            let mapping_quality = output.mapping_quality();
                            let retained_ranges = output.retained_query_intervals();
                            let first = selected.mate1();
                            let second = selected.mate2();
                            let soft_clipped = retained_ranges[0].start != 0
                                || retained_ranges[0].end != pair.first.sequence().len()
                                || retained_ranges[1].start != 0
                                || retained_ranges[1].end != pair.second.sequence().len();
                            let first_length = ReferenceLength::new(
                                reference
                                    .contig_by_ordinal(first.contig_ordinal())
                                    .ok_or_else(|| {
                                        String::from("first paired-end contig is absent")
                                    })?
                                    .sequence()
                                    .len(),
                            );
                            let second_length = ReferenceLength::new(
                                reference
                                    .contig_by_ordinal(second.contig_ordinal())
                                    .ok_or_else(|| {
                                        String::from("second paired-end contig is absent")
                                    })?
                                    .sequence()
                                    .len(),
                            );
                            let first_interval =
                                ReferenceInterval::new(first.start(), first.end(), first_length)
                                    .map_err(|error| error.to_string())?;
                            let second_interval =
                                ReferenceInterval::new(second.start(), second.end(), second_length)
                                    .map_err(|error| error.to_string())?;
                            let first_slab_read = BorrowedAlignmentRead::new(
                                pair.first.sequence(),
                                pair.first.quality(),
                            );
                            let second_slab_read = BorrowedAlignmentRead::new(
                                pair.second.sequence(),
                                pair.second.quality(),
                            );
                            if matches!(output_contract, AlignmentAuxiliaryMode::Minimal) {
                                let first_placement = AlignmentPlacement::new(
                                    first.contig_ordinal(),
                                    first_interval,
                                    first.strand(),
                                    first.distance(),
                                );
                                let second_placement = AlignmentPlacement::new(
                                    second.contig_ordinal(),
                                    second_interval,
                                    second.strand(),
                                    second.distance(),
                                );
                                let pushed = if soft_clipped {
                                    composer.try_push_soft_clipped_ungapped_pair(
                                        reference,
                                        pair.shared_name(),
                                        first_slab_read,
                                        second_slab_read,
                                        retained_ranges[0].clone(),
                                        retained_ranges[1].clone(),
                                        first_placement,
                                        second_placement,
                                        limits,
                                        mapping_quality,
                                    )
                                } else {
                                    composer.try_push_ungapped_pair(
                                        reference,
                                        pair.shared_name(),
                                        first_slab_read,
                                        second_slab_read,
                                        first_placement,
                                        second_placement,
                                        limits,
                                        mapping_quality,
                                    )
                                }
                                .map_err(|error| error.to_string())?;
                                if pushed {
                                    record_worker_ns = record_worker_ns
                                        .saturating_add(record_started.elapsed_ns());
                                    continue;
                                }
                            }
                            let first_sequence = NormalizedSequence::from_bases(
                                pair.first.sequence()[retained_ranges[0].clone()]
                                    .iter()
                                    .copied(),
                            );
                            let second_sequence = NormalizedSequence::from_bases(
                                pair.second.sequence()[retained_ranges[1].clone()]
                                    .iter()
                                    .copied(),
                            );
                            let first_contig = reference
                                .contig_id(first.contig_ordinal())
                                .map_err(|error| error.to_string())?;
                            let second_contig = reference
                                .contig_id(second.contig_ordinal())
                                .map_err(|error| error.to_string())?;
                            let first_alignment = traceback_read_placement(
                                reference,
                                &first_sequence,
                                &first_contig,
                                first_interval,
                                first.strand(),
                                first.distance(),
                            )
                            .map_err(|error| error.to_string())?;
                            let second_alignment = traceback_read_placement(
                                reference,
                                &second_sequence,
                                &second_contig,
                                second_interval,
                                second.strand(),
                                second.distance(),
                            )
                            .map_err(|error| error.to_string())?;
                            if soft_clipped {
                                composer
                                    .push_soft_clipped_retained_unique_pair(
                                        reference,
                                        pair.shared_name(),
                                        first_slab_read,
                                        second_slab_read,
                                        retained_ranges[0].clone(),
                                        retained_ranges[1].clone(),
                                        &first_sequence,
                                        &second_sequence,
                                        &first_alignment,
                                        &second_alignment,
                                        limits,
                                        output_contract,
                                        mapping_quality,
                                    )
                                    .map_err(|error| error.to_string())?;
                            } else {
                                composer
                                    .push_retained_unique_pair_with_mapping_quality(
                                        reference,
                                        pair.shared_name(),
                                        AlignmentRead::new(
                                            &first_sequence,
                                            Some(&pair.first.quality()[retained_ranges[0].clone()]),
                                        ),
                                        AlignmentRead::new(
                                            &second_sequence,
                                            Some(
                                                &pair.second.quality()[retained_ranges[1].clone()],
                                            ),
                                        ),
                                        &first_alignment,
                                        &second_alignment,
                                        limits,
                                        output_contract,
                                        mapping_quality,
                                    )
                                    .map_err(|error| error.to_string())?;
                            }
                            record_worker_ns =
                                record_worker_ns.saturating_add(record_started.elapsed_ns());
                        }
                        let flush_started = MetricsTimer::start(emit_metrics);
                        composer
                            .flush_into(&mut processed, limits)
                            .map_err(|error| error.to_string())?;
                        record_worker_ns =
                            record_worker_ns.saturating_add(flush_started.elapsed_ns());
                        owned.push((
                            start,
                            processed,
                            classes,
                            soft_clip_observation,
                            mate_rescue_observation,
                        ));
                    }
                }
                Ok((owned, mapping_worker_ns, record_worker_ns))
            }));
        }
        let mut chunks = Vec::new();
        let mut mapping_worker_ns = 0_u128;
        let mut record_worker_ns = 0_u128;
        for handle in handles {
            let (owned, mapping_ns, record_ns) = handle
                .join()
                .map_err(|_| invalid("paired mapping/record worker panicked"))?
                .map_err(invalid)?;
            chunks.extend(owned);
            mapping_worker_ns = mapping_worker_ns.saturating_add(mapping_ns);
            record_worker_ns = record_worker_ns.saturating_add(record_ns);
        }
        Ok::<_, io::Error>((chunks, mapping_worker_ns, record_worker_ns))
    })?;
    chunks.0.sort_unstable_by_key(|(start, _, _, _, _)| *start);
    let mut output = Vec::new();
    output.try_reserve_exact(chunks.0.len())?;
    let mut classes = PairClassCounts::default();
    let mut soft_clip = SoftClipObservation::default();
    let mut mate_rescue = MateRescueObservation::default();
    for (_, records, counts, clipped, rescued) in chunks.0 {
        classes.merge(counts);
        soft_clip.merge(clipped);
        mate_rescue.merge(rescued);
        output.push(records);
    }
    Ok(PairedBatchOutput {
        records: output,
        classes,
        soft_clip,
        mate_rescue,
        mapping_worker_ns: chunks.1,
        record_worker_ns: chunks.2,
    })
}

#[allow(clippy::too_many_arguments)]
fn write_batches(
    output: &PathBuf,
    output_file: File,
    header: &SamHeader,
    limits: AlignmentRecordLimits,
    compression_threads: u32,
    compression_level: Option<u8>,
    receiver: Receiver<Vec<AlignmentRecordBatch>>,
    emit_metrics: bool,
) -> Result<WriterObservation, String> {
    let mut writer = BamStagingWriter::create_direct_from_file(
        output,
        output_file,
        header,
        limits,
        compression_threads,
        compression_level,
    )
    .map_err(|error| error.to_string())?;
    let mut bam_write_ns = 0_u128;
    for batch in receiver {
        let write_started = MetricsTimer::start(emit_metrics);
        for chunk in batch {
            for record in chunk.records() {
                writer
                    .write_borrowed_alignment_record(&record)
                    .map_err(|error| error.to_string())?;
            }
        }
        bam_write_ns = bam_write_ns.saturating_add(write_started.elapsed_ns());
    }
    let finalize_started = MetricsTimer::start(emit_metrics);
    let records = writer.finish_direct().map_err(|error| error.to_string())?;
    Ok(WriterObservation {
        records,
        bam_write_ns,
        finalize_ns: finalize_started.elapsed_ns(),
    })
}

fn decode_batches(
    read1: &Path,
    read2: &Path,
    batch_size: usize,
    queue_batches: usize,
    sender: &std::sync::mpsc::SyncSender<PairedInputBatch>,
    emit_metrics: bool,
) -> Result<u128, String> {
    let (first_sender, first_receiver) = sync_channel(queue_batches);
    let (second_sender, second_receiver) = sync_channel(queue_batches);
    let first_path = read1.to_path_buf();
    let second_path = read2.to_path_buf();
    let first = thread::spawn(move || {
        decode_read_batches(first_path, batch_size, &first_sender, emit_metrics)
    });
    let second = thread::spawn(move || {
        decode_read_batches(second_path, batch_size, &second_sender, emit_metrics)
    });
    loop {
        let first_batch = first_receiver.recv();
        let second_batch = second_receiver.recv();
        let (first_batch, second_batch) = match (first_batch, second_batch) {
            (Ok(first), Ok(second)) => (first, second),
            (Err(_), Err(_)) => break,
            _ => {
                return Err(String::from(
                    "paired FASTQ inputs have different record counts",
                ));
            }
        };
        if first_batch.len() != second_batch.len() {
            return Err(String::from(
                "parallel FASTQ batches have different lengths",
            ));
        }
        sender
            .send(PairedInputBatch {
                first: first_batch,
                second: second_batch,
            })
            .map_err(|_| String::from("FASTQ consumer ended before the producer"))?;
    }
    let first_ns = first
        .join()
        .map_err(|_| String::from("R1 decoder panicked"))??;
    let second_ns = second
        .join()
        .map_err(|_| String::from("R2 decoder panicked"))??;
    Ok(first_ns.max(second_ns))
}

fn decode_read_batches(
    path: PathBuf,
    batch_records: usize,
    sender: &std::sync::mpsc::SyncSender<FastqRecordBatch>,
    emit_metrics: bool,
) -> Result<u128, String> {
    let started = MetricsTimer::start(emit_metrics);
    let mut reader = DecodedFastqReader::open(path, alignment_fastq_limits())
        .map_err(|error| error.to_string())?;
    loop {
        let batch = reader
            .next_batch(batch_records)
            .map_err(|error| error.to_string())?;
        if batch.is_empty() {
            break;
        }
        let reached_eof = batch.len() < batch_records;
        sender
            .send(batch)
            .map_err(|_| String::from("FASTQ pairing worker ended before its decoder"))?;
        if reached_eof {
            break;
        }
    }
    reader.close().map_err(|error| error.to_string())?;
    Ok(started.elapsed_ns())
}
