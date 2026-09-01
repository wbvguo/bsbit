//! Standard FASTQ-to-BAM alignment command.
//!
//! Single-end and paired-end input share the persisted combined index, bounded
//! d3/d5 verification core, canonical traceback, record construction, BAM
//! compression/finalization, and create-only publication path.

use std::collections::BTreeSet;
use std::error::Error;
use std::io;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::mpsc::{Receiver, sync_channel};
use std::thread;
use std::time::Instant;

use bsbit_align::library::{PairedLibraryProfile, TemplateSpan, TemplateSpanBounds};
use bsbit_align::materialize::traceback_read_placement;
use bsbit_align::paired_end::{
    PAIRED_ALIGNMENT_BATCH_SIZE, PAIRED_MAX_EDIT_DISTANCE, PairedAlignmentOptions, PairedSearchMode,
};
use bsbit_align::paired_end::{PairMappingStatus, PairedBatchAligner};
use bsbit_align::single_end::SingleSearchMode;
use bsbit_core::coordinate::{ReferenceInterval, ReferenceLength};
use bsbit_core::sequence::NormalizedSequence;
use bsbit_hts::{
    AlignmentAuxiliaryMode, AlignmentPlacement, AlignmentRead, AlignmentRecordBatch,
    AlignmentRecordLimits, BamStagingWriter, BorrowedAlignmentRead, BorrowedFastqRecord,
    BsbitAlignmentMode, BsbitProgramProvenance, DecodedFastqReader, FastqRecordBatch, SamHeader,
    TextRecordLimits,
};
use bsbit_index::reference::ReferenceIndex;
use bsbit_index::storage::combined::load_combined_reference_catalog;

use super::internal_search_file_prefix;
use crate::command::single_end::{SingleEndCommandOptions, run_single_end};
use crate::cpu_placement::CpuPlacement;
use crate::record_composition::{PairedRecordComposer, build_sam_header};

const SCHEMA: &str = "bsbit-alignment-metrics-v1";

#[derive(Clone, Copy)]
struct MetricsTimer(Option<Instant>);

impl MetricsTimer {
    fn start(enabled: bool) -> Self {
        Self(enabled.then(Instant::now))
    }

    fn elapsed_ns(self) -> u128 {
        self.0.map_or(0, |started| started.elapsed().as_nanos())
    }
}

pub(crate) const HELP: &str = r"bsbit align - standard bisulfite read alignment

USAGE:
  bsbit align --index PATH --read1 PATH [--read2 PATH] --output-bam PATH [OPTIONS]

REQUIRED:
  --index PATH                       complete index created by `bsbit index`
  -1, --read1 PATH                   single-end FASTQ, or R1 FASTQ when paired
  --output-bam PATH                  create-only published BAM path

INPUT LAYOUT:
  --read1 only                       directional single-end alignment
  --read1 and --read2                synchronized directional paired-end alignment
                                      (add --non-directional for four-strand PE)

OPTIONAL INPUT:
  -2, --read2 PATH                   R2 FASTQ; requires --read1

OPTIONS FOR BOTH LAYOUTS:
  --sensitive                        complete a wider bounded candidate frontier
  --threads N                        mapping workers; default: 1
  --bam-threads N                    BGZF workers; default: 1
  --bam-compression-level LEVEL      default|0..9; default: 1

PAIRED-END OPTIONS:
  --batch-pairs N                    default: 16384
  --alignment-queue-batches N        default: 2
  --output-contract CONTRACT         minimal|bismark; default: minimal
  --non-directional                  search all four bisulfite strands
  --mapped-only                      omit truly unmapped primary records
  --metrics                          write the full profiling TSV to stdout
  --min-template-span N              default: 0
  --max-template-span N              default: 1000

Single-end alignment uses the same persisted combined index and bounded d3/d5
verification core as paired-end alignment. Unique single reads receive numeric
MAPQ from their existing score-separation and repeat evidence; tied best
placements use MAPQ 0. Its caller-compatible directional-single BAM is accepted
by `bsbit call` after coordinate sorting, duplicate handling, and indexing.

Without --sensitive, default mode runs the low-latency d3 pass plus an
incremental d5 fallback. For single-end input, --sensitive completes the wider
bounded seed frontier before d5 verification, classification, and MAPQ. For
paired-end input it also enables the qualified pair-specific recovery policy.
Inputs may remain gzip-compressed; pre-decompression is not required or recommended.
";

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

type InputFastqBatch = PairedInputBatch;

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct Options {
    index: PathBuf,
    layout: ReadLayout,
    read1: PathBuf,
    read2: Option<PathBuf>,
    output_bam: PathBuf,
    batch_pairs: usize,
    alignment_queue_batches: usize,
    threads: usize,
    bam_threads: u32,
    bam_compression_level: Option<u8>,
    output_contract: AlignmentAuxiliaryMode,
    library_profile: PairedLibraryProfile,
    search_mode: PairedSearchMode,
    read_output: ReadOutputMode,
    minimum_template_span: u64,
    maximum_template_span: u64,
    emit_metrics: bool,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum ReadLayout {
    SingleEnd,
    PairedEnd,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
enum ReadOutputMode {
    #[default]
    Complete,
    MappedOnly,
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
    finalize_publish_ns: u128,
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

pub(super) fn parse(arguments: &[String]) -> Result<super::Action, crate::CliError> {
    if matches!(arguments, [value] if value == "--help" || value == "-h") {
        return Ok(super::Action::Help(HELP));
    }
    parse_options_from(arguments.iter().map(std::ffi::OsString::from))
        .map(super::Action::Align)
        .map_err(|error| crate::CliError::usage(error.to_string()))
}

pub(crate) fn run(options: Options) -> Result<(), Box<dyn Error>> {
    let process_started = MetricsTimer::start(options.emit_metrics);
    if matches!(options.layout, ReadLayout::SingleEnd) {
        return run_standard_single_from_options(options);
    }
    let cpu_placement = CpuPlacement::detect(options.threads);
    let (sender, receiver) = sync_channel(32);
    let read1 = options.read1.clone();
    let read2 = options
        .read2
        .clone()
        .expect("paired layout was validated with read 2");
    let batch_pairs = options.batch_pairs;
    let emit_metrics = options.emit_metrics;
    let producer_cpu_placement = cpu_placement.clone();
    let producer = thread::spawn(move || {
        producer_cpu_placement.pin_auxiliary_worker();
        decode_batches(&read1, &read2, batch_pairs, &sender, emit_metrics)
    });
    let limits = AlignmentRecordLimits::default();
    let (reference, header, reference_load_ns) = load_alignment_reference(&options, limits)?;
    let (alignment_sender, alignment_receiver) = sync_channel(options.alignment_queue_batches);
    let output_bam = options.output_bam.clone();
    let bam_threads = options.bam_threads;
    let bam_compression_level = options.bam_compression_level;
    let emit_metrics = options.emit_metrics;
    let writer_cpu_placement = cpu_placement.clone();
    let writer = thread::spawn(move || {
        writer_cpu_placement.pin_auxiliary_worker();
        write_batches(
            &output_bam,
            &header,
            limits,
            bam_threads,
            bam_compression_level,
            alignment_receiver,
            emit_metrics,
        )
    });
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
        options.search_mode,
        options.read_output,
        options.emit_metrics,
        options.threads,
        &cpu_placement,
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
    if matches!(options.read_output, ReadOutputMode::Complete) {
        let input_pairs = observation.classes.total();
        let expected_records = observation
            .classes
            .total()
            .checked_mul(2)
            .ok_or_else(|| invalid("input primary-record count overflow"))?;
        if writer_observation.records != expected_records {
            return Err(invalid(format!(
                "read-complete BAM wrote {} primary records for {} input pairs; expected {}",
                writer_observation.records, input_pairs, expected_records
            ))
            .into());
        }
    }
    write_metrics(
        &options,
        &observation,
        &writer_observation,
        reference_load_ns,
        decode_ns,
        process_started.elapsed_ns(),
    );
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
    let semantic_digest = loaded.summary().semantic_digest();
    let reference = loaded.into_index();
    let alignment_mode = match options.library_profile {
        PairedLibraryProfile::Directional => BsbitAlignmentMode::CallerCompatibleDirectionalPaired,
        PairedLibraryProfile::NonDirectional => {
            BsbitAlignmentMode::CallerCompatibleNondirectionalPaired
        }
    };
    let header = build_sam_header(&reference, limits)?.with_bsbit_provenance(
        BsbitProgramProvenance::new(semantic_digest.into_bytes(), alignment_mode),
        limits,
    )?;
    Ok((reference, header, started.elapsed_ns()))
}

fn write_metrics(
    options: &Options,
    observation: &Observation,
    writer_observation: &WriterObservation,
    reference_load_ns: u128,
    decode_ns: u128,
    process_total_ns: u128,
) {
    if !options.emit_metrics {
        return;
    }
    println!(
        "schema\tpairs\tunique\tambiguous\tunmapped\tbam_records\tmapping_threads\tbam_threads\toutput_contract\tlibrary_profile\treference_mode\treference_load_ns\tfastq_decode_ns\tbam_write_ns\tbam_finalize_publish_ns\tprocess_total_ns\tbatch_processing_ns\tmapping_worker_total_ns\trecord_worker_total_ns\talignment_queue_batches\twriter_queue_wait_ns\twriter_queue_sends\tbam_compression_level\tmax_edit_distance\tsoft_clip_fallback\tsoft_clip_attempted_pairs\tsoft_clip_unique_pairs\tsoft_clip_ambiguous_pairs\tsoft_clip_unmapped_pairs\tsoft_clip_clipped_mates\tsoft_clip_clipped_bases\tmate_rescue\tmate_rescue_attempted_pairs\tmate_rescue_unique_pairs\tmate_rescue_ambiguous_pairs\tmate_rescue_unmapped_pairs\tsearch_mode\tmapq_policy\tmapq_zero_output\tstrategy_id\tread_output"
    );
    let fields = [
        SCHEMA.to_owned(),
        observation.classes.total().to_string(),
        observation.classes.unique.to_string(),
        observation.classes.ambiguous.to_string(),
        observation.classes.unmapped.to_string(),
        writer_observation.records.to_string(),
        options.threads.to_string(),
        options.bam_threads.to_string(),
        output_contract_name(options.output_contract).to_owned(),
        library_profile_name(options.library_profile).to_owned(),
        "indexed-reference".to_owned(),
        reference_load_ns.to_string(),
        decode_ns.to_string(),
        writer_observation.bam_write_ns.to_string(),
        writer_observation.finalize_publish_ns.to_string(),
        process_total_ns.to_string(),
        observation.batch_processing_ns.to_string(),
        observation.mapping_worker_total_ns.to_string(),
        observation.record_worker_total_ns.to_string(),
        options.alignment_queue_batches.to_string(),
        observation.writer_queue_wait_ns.to_string(),
        observation.writer_queue_sends.to_string(),
        options
            .bam_compression_level
            .map_or_else(|| "default".to_owned(), |level| level.to_string()),
        PAIRED_MAX_EDIT_DISTANCE.to_string(),
        soft_clip_fallback_name(options.search_mode).to_owned(),
        observation.soft_clip.attempted_pairs.to_string(),
        observation.soft_clip.unique_pairs.to_string(),
        observation.soft_clip.ambiguous_pairs.to_string(),
        observation.soft_clip.unmapped_pairs.to_string(),
        observation.soft_clip.clipped_mates.to_string(),
        observation.soft_clip.clipped_bases.to_string(),
        mate_rescue_name(options.search_mode).to_owned(),
        observation.mate_rescue.attempted.to_string(),
        observation.mate_rescue.unique.to_string(),
        observation.mate_rescue.ambiguous.to_string(),
        observation.mate_rescue.unmapped.to_string(),
        search_mode_name(options.search_mode).to_owned(),
        "qualified".to_owned(),
        "all".to_owned(),
        strategy_id(options).to_owned(),
        read_output_name(options.read_output).to_owned(),
    ];
    println!("{}", fields.join("\t"));
}

fn run_standard_single_from_options(options: Options) -> Result<(), Box<dyn Error>> {
    let search_mode = match options.search_mode {
        PairedSearchMode::Default => SingleSearchMode::Default,
        PairedSearchMode::Sensitive => SingleSearchMode::Sensitive,
    };
    let align_options = SingleEndCommandOptions {
        index: options.index,
        read1: options.read1,
        output_bam: options.output_bam,
        search_mode,
        max_edit_distance: u64::from(PAIRED_MAX_EDIT_DISTANCE),
        batch_records: 1_000,
        threads: u64::try_from(options.threads).expect("validated thread count fits u64"),
        bam_threads: options.bam_threads,
        bam_compression_level: options.bam_compression_level,
        emit_metrics: options.emit_metrics,
    };
    run_single_end(&align_options)
        .map(|_| ())
        .map_err(|error| Box::new(error) as Box<dyn Error>)
}

#[allow(clippy::too_many_arguments)]
fn consume_batches(
    reference: &ReferenceIndex,
    receiver: Receiver<InputFastqBatch>,
    alignment_sender: &std::sync::mpsc::SyncSender<Vec<AlignmentRecordBatch>>,
    observation: &mut Observation,
    bounds: TemplateSpanBounds,
    limits: AlignmentRecordLimits,
    output_contract: AlignmentAuxiliaryMode,
    library_profile: PairedLibraryProfile,
    search_mode: PairedSearchMode,
    read_output: ReadOutputMode,
    emit_metrics: bool,
    threads: usize,
    cpu_placement: &CpuPlacement,
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
    records: &InputFastqBatch,
    bounds: TemplateSpanBounds,
    limits: AlignmentRecordLimits,
    output_contract: AlignmentAuxiliaryMode,
    library_profile: PairedLibraryProfile,
    search_mode: PairedSearchMode,
    read_output: ReadOutputMode,
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
                let mut aligner = PairedBatchAligner::with_capacity(PAIRED_ALIGNMENT_BATCH_SIZE);
                let mut composer = PairedRecordComposer::new();
                let mut pair_reads = Vec::with_capacity(PAIRED_ALIGNMENT_BATCH_SIZE);
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
                        for index in start..end {
                            let pair = records.get(index).ok_or_else(|| {
                                String::from("paired FASTQ batch index is absent")
                            })?;
                            pair_reads.push([pair.first.sequence(), pair.second.sequence()]);
                        }
                        let mapping_started = MetricsTimer::start(emit_metrics);
                        let mapped = aligner
                            .map_pairs_for_output(
                                reference,
                                &pair_reads,
                                PairedAlignmentOptions::primary(
                                    library_profile,
                                    search_mode,
                                    bounds.minimum().get(),
                                    bounds.maximum().get(),
                                ),
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

fn write_batches(
    output_bam: &PathBuf,
    header: &SamHeader,
    limits: AlignmentRecordLimits,
    bam_threads: u32,
    bam_compression_level: Option<u8>,
    receiver: Receiver<Vec<AlignmentRecordBatch>>,
    emit_metrics: bool,
) -> Result<WriterObservation, String> {
    let mut writer = match bam_compression_level {
        Some(level) => BamStagingWriter::create_sibling_with_threads_and_compression_level(
            output_bam,
            header,
            limits,
            bam_threads,
            level,
        ),
        None => {
            BamStagingWriter::create_sibling_with_threads(output_bam, header, limits, bam_threads)
        }
    }
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
    let publication = writer
        .finish()
        .map_err(|error| error.to_string())?
        .publish_create_new(output_bam)
        .map_err(|error| error.to_string())?;
    Ok(WriterObservation {
        records: publication.records_written(),
        bam_write_ns,
        finalize_publish_ns: finalize_started.elapsed_ns(),
    })
}

fn decode_batches(
    read1: &Path,
    read2: &Path,
    batch_pairs: usize,
    sender: &std::sync::mpsc::SyncSender<PairedInputBatch>,
    emit_metrics: bool,
) -> Result<u128, String> {
    let (first_sender, first_receiver) = sync_channel(2);
    let (second_sender, second_receiver) = sync_channel(2);
    let first_path = read1.to_path_buf();
    let second_path = read2.to_path_buf();
    let first = thread::spawn(move || {
        decode_read_batches(first_path, batch_pairs, &first_sender, emit_metrics)
    });
    let second = thread::spawn(move || {
        decode_read_batches(second_path, batch_pairs, &second_sender, emit_metrics)
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
    let mut reader =
        DecodedFastqReader::open(path, TextRecordLimits::MAX).map_err(|error| error.to_string())?;
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

fn option_takes_value(flag: &str) -> bool {
    matches!(
        flag,
        "--index"
            | "--read1"
            | "--read2"
            | "--output-bam"
            | "--batch-pairs"
            | "--alignment-queue-batches"
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
) -> Result<Options, io::Error> {
    let mut index = None;
    let mut read1 = None;
    let mut read2 = None;
    let mut output_bam = None;
    let mut batch_pairs = 16_384_usize;
    let mut alignment_queue_batches = 2_usize;
    let mut threads = 1_usize;
    // One BGZF worker lets record compression overlap mapping.  Zero remains
    // available for callers that require a strictly synchronous writer.
    let mut bam_threads = 1_u32;
    let mut bam_compression_level = Some(1_u8);
    let mut output_contract = AlignmentAuxiliaryMode::Minimal;
    let mut library_profile = PairedLibraryProfile::Directional;
    let mut search_mode = PairedSearchMode::Default;
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
            let requested = PairedSearchMode::Sensitive;
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
            if matches!(library_profile, PairedLibraryProfile::NonDirectional) {
                return Err(invalid("--non-directional may be specified only once"));
            }
            library_profile = PairedLibraryProfile::NonDirectional;
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
        let unsupported_flag = if matches!(library_profile, PairedLibraryProfile::NonDirectional) {
            Some("--non-directional")
        } else if read_output_explicit {
            Some("--mapped-only")
        } else {
            [
                "--batch-pairs",
                "--alignment-queue-batches",
                "--output-contract",
                "--min-template-span",
                "--max-template-span",
            ]
            .into_iter()
            .find(|flag| seen_value_flags.contains(*flag))
        };
        if let Some(flag) = unsupported_flag {
            return Err(invalid(format!("{flag} requires paired input via --read2")));
        }
    }
    Ok(Options {
        index,
        layout,
        read1,
        read2,
        output_bam,
        batch_pairs,
        alignment_queue_batches,
        threads,
        bam_threads,
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

const fn output_contract_name(mode: AlignmentAuxiliaryMode) -> &'static str {
    match mode {
        AlignmentAuxiliaryMode::Minimal => "minimal",
        AlignmentAuxiliaryMode::Bismark => "bismark",
    }
}

const fn library_profile_name(profile: PairedLibraryProfile) -> &'static str {
    match profile {
        PairedLibraryProfile::Directional => "directional",
        PairedLibraryProfile::NonDirectional => "non-directional",
    }
}

const fn soft_clip_fallback_name(mode: PairedSearchMode) -> &'static str {
    match mode {
        PairedSearchMode::Default => "adapter",
        PairedSearchMode::Sensitive => "semi-global",
    }
}

const fn mate_rescue_name(mode: PairedSearchMode) -> &'static str {
    match mode {
        PairedSearchMode::Default => "off",
        PairedSearchMode::Sensitive => "windowed",
    }
}

const fn search_mode_name(mode: PairedSearchMode) -> &'static str {
    match mode {
        PairedSearchMode::Default => "default",
        PairedSearchMode::Sensitive => "sensitive",
    }
}

const fn read_output_name(mode: ReadOutputMode) -> &'static str {
    match mode {
        ReadOutputMode::Complete => "complete",
        ReadOutputMode::MappedOnly => "mapped-only",
    }
}

const fn sensitive_mapq_zero_strategy_id() -> &'static str {
    "sensitive-bounded-integrated-mapq0-all-v1"
}

const fn sensitive_read_complete_strategy_id() -> &'static str {
    "sensitive-bounded-integrated-read-complete-v1"
}
fn strategy_id(options: &Options) -> &'static str {
    // The `balanced-d5` spellings are immutable identifiers recorded by prior
    // qualification reports. They now denote the stable `Default` policy; do
    // not rewrite persisted evidence merely to mirror the enum variant name.
    match (options.search_mode, options.read_output) {
        (PairedSearchMode::Sensitive, ReadOutputMode::Complete) => {
            sensitive_read_complete_strategy_id()
        }
        (PairedSearchMode::Sensitive, ReadOutputMode::MappedOnly) => {
            sensitive_mapq_zero_strategy_id()
        }
        (PairedSearchMode::Default, ReadOutputMode::Complete) => {
            "balanced-d5-adapter-recovery-read-complete-v2"
        }
        (PairedSearchMode::Default, ReadOutputMode::MappedOnly) => {
            "balanced-d5-adapter-recovery-mapq0-all-v2"
        }
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
#[path = "../../tests/whitebox/align.rs"]
mod tests;
