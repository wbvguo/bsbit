//! Regional memory planning, worker scheduling, and ordered aggregation.

use std::collections::BTreeMap;
use std::panic::{AssertUnwindSafe, catch_unwind};
use std::path::Path;
use std::sync::{Arc, Mutex, mpsc};
use std::thread;

use bsbit_hts::IndexedBamReader;

use crate::CallError;
use crate::CallErrorKind;
use crate::call_input::BamReference;
use crate::evidence::fragment::{
    EvidenceContext, EvidenceFilter, EvidenceWorkspace, for_each_region_fragment,
};
use crate::meth::Parameters as MethParameters;
use crate::meth::aggregation::{
    DenseMethRegion, accumulate_meth_fragment, meth_dense_bytes_per_block,
};
use crate::reference_context::CallReferenceReader;
use crate::region::CallRegion;
use crate::snp::candidate::{CandidateRegion, CandidateSite, snp_region_bytes_per_block};
use crate::snp::likelihood::{LikelihoodRegion, likelihood_site_bytes};
use crate::snp::result::{SnpConfig, VariantCall};

#[cfg(test)]
use crate::meth::aggregation::{CallKind, SiteCounts, SiteKey};

const MAX_REGION_BASES: u32 = 1 << 20;
const MIN_REGION_BASES: u32 = 1 << 12;
const REGION_REORDER_FACTOR: usize = 2;
const REGION_STATE_BUDGET_BYTES: usize = 256 << 20;
const LIKELIHOOD_STATE_BUDGET_BYTES: usize = 256 << 20;
const MAX_SNP_LIKELIHOOD_BATCH_SITES: usize = 4_096;
const MIN_SNP_LIKELIHOOD_BATCH_SITES: usize = 256;
const SNP_LIKELIHOOD_WINDOW_BASES: u32 = 1 << 16;

enum RegionWorkerMessage {
    Ready {
        worker: usize,
        result: Result<(), CallError>,
    },
    Region {
        worker: usize,
        ordinal: usize,
        result: Result<RegionAggregation, CallError>,
    },
    Closed {
        worker: usize,
        result: Result<(), CallError>,
    },
    Panicked {
        worker: usize,
        message: String,
    },
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) enum IndexedCallMode {
    Meth(MethParameters),
    Snp(SnpConfig),
    Joint(SnpConfig),
    #[cfg(test)]
    Panic,
}

impl IndexedCallMode {
    const fn region_bytes_per_block(self) -> usize {
        let meth = meth_dense_bytes_per_block();
        match self {
            Self::Meth(_) => meth,
            Self::Snp(_) => snp_region_bytes_per_block(),
            Self::Joint(_) => meth + snp_region_bytes_per_block(),
            #[cfg(test)]
            Self::Panic => meth,
        }
    }

    const fn evidence_filter(self) -> EvidenceFilter {
        match self {
            Self::Meth(parameters) => EvidenceFilter::new(
                parameters.minimum_base_quality,
                parameters.minimum_mapping_quality,
                true,
            ),
            Self::Snp(config) => EvidenceFilter::new(
                config.minimum_base_quality,
                config.minimum_mapping_quality,
                false,
            ),
            Self::Joint(config) => EvidenceFilter::new(
                config.minimum_base_quality,
                config.minimum_mapping_quality,
                true,
            ),
            #[cfg(test)]
            Self::Panic => EvidenceFilter::new(0, 0, true),
        }
    }
}

pub(super) fn region_bases_for(mode: IndexedCallMode, threads: usize) -> u32 {
    debug_assert!(threads > 0);
    let bounded_regions = threads.saturating_mul(REGION_REORDER_FACTOR).max(1);
    let bytes_per_region = REGION_STATE_BUDGET_BYTES / bounded_regions;
    let desired_blocks = bytes_per_region / mode.region_bytes_per_block();
    let minimum_blocks = usize::try_from(MIN_REGION_BASES / u64::BITS)
        .expect("minimum region block count fits usize");
    let maximum_blocks = usize::try_from(MAX_REGION_BASES / u64::BITS)
        .expect("maximum region block count fits usize");
    let blocks = desired_blocks.clamp(minimum_blocks, maximum_blocks);
    u32::try_from(blocks * u64::BITS as usize).expect("bounded calling region length fits u32")
}

fn likelihood_batch_sites_for(worker_count: usize) -> usize {
    debug_assert!(worker_count > 0);
    let bytes_per_worker = LIKELIHOOD_STATE_BUDGET_BYTES / worker_count;
    let lookup_bytes = usize::try_from(SNP_LIKELIHOOD_WINDOW_BASES)
        .expect("likelihood window fits usize")
        * std::mem::size_of::<u16>();
    let site_bytes = likelihood_site_bytes().max(1);
    bytes_per_worker
        .saturating_sub(lookup_bytes)
        .checked_div(site_bytes)
        .unwrap_or(0)
        .clamp(
            MIN_SNP_LIKELIHOOD_BATCH_SITES,
            MAX_SNP_LIKELIHOOD_BATCH_SITES,
        )
}

#[derive(Debug, Default)]
pub(super) struct RegionAggregation {
    pub(super) meth: Option<DenseMethRegion>,
    pub(super) variants: Vec<(u32, VariantCall)>,
}

#[cfg(test)]
#[derive(Debug, Default)]
struct CollectedRegionAggregation {
    meth_sites: Vec<(SiteKey, SiteCounts)>,
    variants: Vec<(u32, VariantCall)>,
}

#[cfg(test)]
fn run_indexed_region_workers(
    path: &Path,
    reference_path: &Path,
    references: &[BamReference],
    regions: &[CallRegion],
    worker_count: usize,
) -> Result<Vec<(SiteKey, SiteCounts)>, CallError> {
    let mut sites = Vec::new();
    stream_indexed_region_workers_mode(
        path,
        references,
        regions,
        worker_count,
        IndexedCallMode::Meth(MethParameters::default()),
        reference_path,
        |region| {
            let meth = region
                .meth
                .ok_or_else(|| CallError::operation("methylation region result is missing"))?;
            sites.extend(meth.into_sites()?);
            Ok(())
        },
    )?;
    Ok(sites)
}

#[cfg(test)]
fn collect_indexed_region_workers_mode(
    path: &Path,
    reference_path: &Path,
    references: &[BamReference],
    regions: &[CallRegion],
    worker_count: usize,
    mode: IndexedCallMode,
) -> Result<CollectedRegionAggregation, CallError> {
    let mut aggregation = CollectedRegionAggregation::default();
    stream_indexed_region_workers_mode(
        path,
        references,
        regions,
        worker_count,
        mode,
        reference_path,
        |mut region| {
            if let Some(meth) = region.meth.take() {
                aggregation.meth_sites.extend(meth.into_sites()?);
            }
            aggregation.variants.append(&mut region.variants);
            Ok(())
        },
    )?;
    Ok(aggregation)
}

pub(super) fn stream_indexed_region_workers_mode(
    path: &Path,
    references: &[BamReference],
    regions: &[CallRegion],
    worker_count: usize,
    mode: IndexedCallMode,
    reference_path: &Path,
    mut consume: impl FnMut(RegionAggregation) -> Result<(), CallError>,
) -> Result<(), CallError> {
    if regions.is_empty() {
        return Ok(());
    }
    debug_assert!(worker_count > 0);
    let lookahead = worker_count
        .saturating_mul(REGION_REORDER_FACTOR)
        .min(regions.len())
        .max(1);
    let likelihood_batch_sites = likelihood_batch_sites_for(worker_count);
    thread::scope(|scope| {
        let (task_sender, task_receiver) = mpsc::sync_channel::<CallRegion>(lookahead);
        let task_receiver = Arc::new(Mutex::new(task_receiver));
        let (result_sender, result_receiver) = mpsc::sync_channel::<RegionWorkerMessage>(lookahead);
        for worker in 0..worker_count {
            let task_receiver = Arc::clone(&task_receiver);
            let result_sender = result_sender.clone();
            let worker_plan = IndexedRegionWorkerPlan {
                worker,
                path,
                references,
                mode,
                reference_path,
                likelihood_batch_sites,
            };
            scope.spawn(move || {
                run_indexed_region_worker(worker_plan, &task_receiver, &result_sender);
            });
        }
        drop(result_sender);

        await_worker_readiness(&result_receiver, worker_count)?;
        dispatch_initial_regions(&task_sender, regions, lookahead)?;
        consume_region_results(
            &task_sender,
            &result_receiver,
            regions,
            lookahead,
            &mut consume,
        )?;
        drop(task_sender);
        await_worker_closure(&result_receiver, worker_count)
    })
}

#[derive(Clone, Copy)]
struct IndexedRegionWorkerPlan<'a> {
    worker: usize,
    path: &'a Path,
    references: &'a [BamReference],
    mode: IndexedCallMode,
    reference_path: &'a Path,
    likelihood_batch_sites: usize,
}

fn run_indexed_region_worker(
    plan: IndexedRegionWorkerPlan<'_>,
    task_receiver: &Mutex<mpsc::Receiver<CallRegion>>,
    result_sender: &mpsc::SyncSender<RegionWorkerMessage>,
) {
    let outcome = catch_unwind(AssertUnwindSafe(|| {
        indexed_region_worker_body(plan, task_receiver, result_sender);
    }));
    if let Err(payload) = outcome {
        let _ = result_sender.send(RegionWorkerMessage::Panicked {
            worker: plan.worker,
            message: panic_payload_message(&payload),
        });
    }
}

fn indexed_region_worker_body(
    plan: IndexedRegionWorkerPlan<'_>,
    task_receiver: &Mutex<mpsc::Receiver<CallRegion>>,
    result_sender: &mpsc::SyncSender<RegionWorkerMessage>,
) {
    let IndexedRegionWorkerPlan {
        worker,
        path,
        references,
        mode,
        reference_path,
        likelihood_batch_sites,
    } = plan;
    let mut reader = match IndexedBamReader::open(path) {
        Ok(reader) => reader,
        Err(error) => {
            let _ = result_sender.send(RegionWorkerMessage::Ready {
                worker,
                result: Err(CallError::with_source(
                    CallErrorKind::Input,
                    format!("worker {worker}: open indexed BAM {}", path.display()),
                    error,
                )),
            });
            return;
        }
    };
    let mut reference_reader = match CallReferenceReader::open(reference_path, references) {
        Ok(reader) => reader,
        Err(error) => {
            let _ = result_sender.send(RegionWorkerMessage::Ready {
                worker,
                result: Err(
                    error.with_context(format!("worker {worker}: open indexed reference FASTA"))
                ),
            });
            return;
        }
    };
    if result_sender
        .send(RegionWorkerMessage::Ready {
            worker,
            result: Ok(()),
        })
        .is_err()
    {
        return;
    }
    let mut workspace = EvidenceWorkspace::default();
    loop {
        let task = task_receiver
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .recv();
        let Ok(region) = task else {
            break;
        };
        let result = aggregate_indexed_region_mode(
            &mut reader,
            references,
            region,
            mode,
            likelihood_batch_sites,
            &mut workspace,
            &mut reference_reader,
        )
        .map_err(|error| {
            error.with_context(format!(
                "worker {worker}: call region {}:{}-{}",
                region.reference, region.start, region.end
            ))
        });
        if result_sender
            .send(RegionWorkerMessage::Region {
                worker,
                ordinal: region.ordinal,
                result,
            })
            .is_err()
        {
            return;
        }
    }
    let bam_close = reader.close().map_err(|error| {
        CallError::with_source(
            CallErrorKind::Input,
            format!("worker {worker}: close indexed BAM {}", path.display()),
            error,
        )
    });
    let reference_close = reference_reader.close();
    let close = bam_close.and(reference_close);
    let _ = result_sender.send(RegionWorkerMessage::Closed {
        worker,
        result: close,
    });
}

fn dispatch_initial_regions(
    sender: &mpsc::SyncSender<CallRegion>,
    regions: &[CallRegion],
    count: usize,
) -> Result<(), CallError> {
    for region in &regions[..count] {
        dispatch_region(sender, *region)?;
    }
    Ok(())
}

fn dispatch_region(
    sender: &mpsc::SyncSender<CallRegion>,
    region: CallRegion,
) -> Result<(), CallError> {
    sender
        .send(region)
        .map_err(|_| CallError::operation("indexed BAM task workers stopped during dispatch"))
}

fn consume_region_results(
    task_sender: &mpsc::SyncSender<CallRegion>,
    result_receiver: &mpsc::Receiver<RegionWorkerMessage>,
    regions: &[CallRegion],
    initial_sent: usize,
    consume: &mut impl FnMut(RegionAggregation) -> Result<(), CallError>,
) -> Result<(), CallError> {
    let mut next_to_send = initial_sent;
    let mut next_to_consume = 0_usize;
    let mut reorder = BTreeMap::<usize, Result<RegionAggregation, CallError>>::new();
    while next_to_consume < regions.len() {
        let message = result_receiver.recv().map_err(|_| {
            CallError::operation("indexed BAM workers stopped before returning all regions")
        })?;
        let RegionWorkerMessage::Region {
            worker,
            ordinal,
            result,
        } = message
        else {
            return Err(unexpected_region_message(message));
        };
        if reorder.insert(ordinal, result).is_some() {
            return Err(CallError::operation(format!(
                "worker {worker} returned duplicate region ordinal {ordinal}"
            )));
        }
        while next_to_consume < regions.len() {
            let expected = regions[next_to_consume].ordinal;
            let Some(result) = reorder.remove(&expected) else {
                break;
            };
            consume(result?)?;
            next_to_consume += 1;
            if let Some(region) = regions.get(next_to_send) {
                dispatch_region(task_sender, *region)?;
                next_to_send += 1;
            }
        }
    }
    if reorder.is_empty() {
        Ok(())
    } else {
        Err(CallError::operation(
            "indexed BAM workers returned unexpected region ordinals",
        ))
    }
}

fn unexpected_region_message(message: RegionWorkerMessage) -> CallError {
    match message {
        RegionWorkerMessage::Panicked { worker, message } => {
            CallError::operation(format!("indexed BAM worker {worker} panicked: {message}"))
        }
        RegionWorkerMessage::Ready { worker, .. } => CallError::operation(format!(
            "indexed BAM worker {worker} sent duplicate readiness"
        )),
        RegionWorkerMessage::Closed { worker, .. } => CallError::operation(format!(
            "indexed BAM worker {worker} closed before all regions completed"
        )),
        RegionWorkerMessage::Region { .. } => {
            CallError::operation("internal region-message classification failed")
        }
    }
}

fn await_worker_closure(
    receiver: &mpsc::Receiver<RegionWorkerMessage>,
    worker_count: usize,
) -> Result<(), CallError> {
    let mut close_errors = BTreeMap::new();
    for _ in 0..worker_count {
        match receiver
            .recv()
            .map_err(|_| CallError::operation("indexed BAM workers stopped before closing"))?
        {
            RegionWorkerMessage::Closed { worker, result } => {
                if let Err(error) = result {
                    close_errors.insert(worker, error);
                }
            }
            RegionWorkerMessage::Panicked { worker, message } => {
                return Err(CallError::operation(format!(
                    "indexed BAM worker {worker} panicked while closing: {message}"
                )));
            }
            RegionWorkerMessage::Ready { worker, .. }
            | RegionWorkerMessage::Region { worker, .. } => {
                return Err(CallError::operation(format!(
                    "indexed BAM worker {worker} returned an unexpected close message"
                )));
            }
        }
    }
    close_errors
        .into_iter()
        .next()
        .map_or(Ok(()), |(_, error)| Err(error))
}

fn panic_payload_message(payload: &Box<dyn std::any::Any + Send>) -> String {
    if let Some(message) = payload.downcast_ref::<&str>() {
        return (*message).to_owned();
    }
    if let Some(message) = payload.downcast_ref::<String>() {
        return message.clone();
    }
    String::from("non-string panic payload")
}

fn await_worker_readiness(
    receiver: &mpsc::Receiver<RegionWorkerMessage>,
    worker_count: usize,
) -> Result<(), CallError> {
    let mut ready = vec![false; worker_count];
    for _ in 0..worker_count {
        match receiver.recv().map_err(|_| {
            CallError::operation("indexed BAM workers stopped during initialization")
        })? {
            RegionWorkerMessage::Ready { worker, result } => {
                let slot = ready.get_mut(worker).ok_or_else(|| {
                    CallError::operation(format!(
                        "indexed BAM worker reported invalid identity {worker}"
                    ))
                })?;
                if std::mem::replace(slot, true) {
                    return Err(CallError::operation(format!(
                        "indexed BAM worker {worker} reported readiness twice"
                    )));
                }
                result?;
            }
            RegionWorkerMessage::Panicked { worker, message } => {
                return Err(CallError::operation(format!(
                    "indexed BAM worker {worker} panicked during initialization: {message}"
                )));
            }
            RegionWorkerMessage::Region { worker, .. }
            | RegionWorkerMessage::Closed { worker, .. } => {
                return Err(CallError::operation(format!(
                    "indexed BAM worker {worker} sent an out-of-order initialization message"
                )));
            }
        }
    }
    Ok(())
}

fn aggregate_indexed_region_mode(
    reader: &mut IndexedBamReader,
    references: &[BamReference],
    region: CallRegion,
    mode: IndexedCallMode,
    likelihood_batch_sites: usize,
    workspace: &mut EvidenceWorkspace,
    reference_reader: &mut CallReferenceReader,
) -> Result<RegionAggregation, CallError> {
    #[cfg(test)]
    if matches!(mode, IndexedCallMode::Panic) {
        panic!("injected regional worker panic");
    }
    let meth_parameters = match mode {
        IndexedCallMode::Meth(parameters) => Some(parameters),
        IndexedCallMode::Joint(config) => Some(MethParameters {
            minimum_base_quality: config.minimum_base_quality,
            minimum_mapping_quality: config.minimum_mapping_quality,
        }),
        IndexedCallMode::Snp(_) => None,
        #[cfg(test)]
        IndexedCallMode::Panic => None,
    };
    let mut meth = meth_parameters
        .is_some()
        .then(|| DenseMethRegion::new(region.reference, region.start, region.end))
        .transpose()?;
    let snp_config = match mode {
        IndexedCallMode::Meth(_) => None,
        IndexedCallMode::Snp(config) | IndexedCallMode::Joint(config) => Some(config),
        #[cfg(test)]
        IndexedCallMode::Panic => None,
    };
    let mut candidates = snp_config
        .map(|config| CandidateRegion::new(region.start, region.end, config))
        .transpose()?;
    let reference_window = reference_reader.fetch_context_window(region, references)?;

    let evidence_filter = mode.evidence_filter();
    let evidence_context = if meth_parameters.is_some() {
        EvidenceContext::WithCytosineContext(&reference_window)
    } else {
        EvidenceContext::WithoutCytosineContext(&reference_window)
    };
    let likelihood_evidence_context = EvidenceContext::WithoutCytosineContext(&reference_window);
    for_each_region_fragment(
        reader,
        references,
        region,
        evidence_filter,
        evidence_context,
        workspace,
        |observations| {
            if let (Some(meth), Some(parameters)) = (&mut meth, meth_parameters) {
                accumulate_meth_fragment(observations, meth, parameters)?;
            }
            if let Some(candidates) = &mut candidates {
                candidates.observe_fragment(observations)?;
            }
            Ok(())
        },
    )?;

    let mut variants = Vec::new();
    if let (Some(candidates), Some(config)) = (candidates, snp_config) {
        let candidates = candidates.candidates()?;
        let mut batch_start = 0_usize;
        while batch_start < candidates.len() {
            let first_position = candidates[batch_start].position;
            let maximum_position = first_position.saturating_add(SNP_LIKELIHOOD_WINDOW_BASES);
            let mut batch_end = batch_start + 1;
            while batch_end < candidates.len()
                && batch_end - batch_start < likelihood_batch_sites
                && candidates[batch_end].position < maximum_position
            {
                batch_end += 1;
            }
            let candidate_batch = &candidates[batch_start..batch_end];
            let mut likelihoods = LikelihoodRegion::new(candidate_batch, config)?;
            let query_region = candidate_query_region(region, candidate_batch)?;
            for_each_region_fragment(
                reader,
                references,
                query_region,
                evidence_filter,
                likelihood_evidence_context,
                workspace,
                |observations| likelihoods.observe_fragment(observations),
            )?;
            let calls = likelihoods.calls()?;
            variants.try_reserve(calls.len()).map_err(|error| {
                CallError::with_source(
                    CallErrorKind::Calling,
                    format!("reserve {} regional SNP calls", calls.len()),
                    error,
                )
            })?;
            variants.extend(calls.into_iter().map(|call| (region.reference, call)));
            batch_start = batch_end;
        }
    }

    Ok(RegionAggregation { meth, variants })
}

fn candidate_query_region(
    parent: CallRegion,
    candidates: &[CandidateSite],
) -> Result<CallRegion, CallError> {
    let first = candidates
        .first()
        .ok_or_else(|| CallError::operation("SNP likelihood batch is empty"))?
        .position;
    let last = candidates
        .last()
        .ok_or_else(|| CallError::operation("SNP likelihood batch is empty"))?
        .position;
    let end = last
        .checked_add(1)
        .ok_or_else(|| CallError::operation("SNP likelihood query end overflowed u32"))?;
    if first < parent.start || end > parent.end {
        return Err(CallError::operation(
            "SNP likelihood batch escaped its parent calling region",
        ));
    }
    Ok(CallRegion {
        ordinal: parent.ordinal,
        reference: parent.reference,
        start: first,
        end,
    })
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;
    use std::fs;
    use std::io::Read;
    use std::path::PathBuf;
    use std::time::{SystemTime, UNIX_EPOCH};

    use bsbit_core::bisulfite::BisulfiteStrand;
    use bsbit_hts::{
        AlignmentAuxiliaryMode, AlignmentCigarOp, AlignmentCigarRun, AlignmentRecordLimits,
        BamStagingWriter, BorrowedAlignmentRecord, DecodedReader, SamHeader, SamHeaderReference,
        SamSortOrder, TextOutputCompression, TextStagingWriter, build_bam_index_create_new,
    };

    use crate::evidence::fragment::merge_contexts;
    use crate::evidence::{ContextClass, CytosineContext, EvidenceObservation, EvidenceStrand};
    use crate::meth::OutputFormat as MethylationOutputFormat;
    use crate::meth::output::{
        UnresolvedContextSummary, render_bed, render_cgmap, render_region as render_meth_region,
    };
    use crate::snp::output::render_header as render_streaming_vcf_header;

    use super::*;

    fn unique_path(label: &str) -> PathBuf {
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        std::env::temp_dir().join(format!(
            "bsbit-call-meth-{label}-{}-{nonce}",
            std::process::id()
        ))
    }

    const INDEXED_FIXTURE_REFERENCE: &[u8] = b"ACGTTGCACTGATCGATGCTAGCTACGATCGTTCGAGTACCTGACGTA";

    struct IndexedXgFixture {
        directory: PathBuf,
        bam: PathBuf,
        fasta: PathBuf,
        references: Vec<BamReference>,
    }

    impl Drop for IndexedXgFixture {
        fn drop(&mut self) {
            let _ = fs::remove_dir_all(&self.directory);
        }
    }

    fn indexed_xg_fixture(label: &str) -> IndexedXgFixture {
        let directory = unique_path(label);
        fs::create_dir(&directory).expect("fixture directory is fresh");

        let mut observed = INDEXED_FIXTURE_REFERENCE.to_vec();
        observed[0] = b'G';

        let fasta = directory.join("reference.fa");
        let mut fasta_contents = b">chr1\n".to_vec();
        fasta_contents.extend_from_slice(INDEXED_FIXTURE_REFERENCE);
        fasta_contents.push(b'\n');
        fs::write(&fasta, fasta_contents).expect("fixture FASTA writes");
        fs::write(
            fasta.with_extension("fa.fai"),
            format!(
                "chr1\t{}\t6\t{}\t{}\n",
                INDEXED_FIXTURE_REFERENCE.len(),
                INDEXED_FIXTURE_REFERENCE.len(),
                INDEXED_FIXTURE_REFERENCE.len() + 1
            ),
        )
        .expect("fixture FAI writes");

        let limits = AlignmentRecordLimits::default();
        let reference_length =
            u64::try_from(INDEXED_FIXTURE_REFERENCE.len()).expect("fixture length fits u64");
        let mut digest = bsbit_core::reference::ReferenceSemanticDigestBuilder::new(1);
        digest
            .push_ascii_contig(b"chr1", INDEXED_FIXTURE_REFERENCE)
            .expect("fixture semantic digest input");
        let header = SamHeader::new(
            vec![
                SamHeaderReference::new(0, b"chr1", reference_length)
                    .expect("fixture dictionary entry"),
            ],
            limits,
        )
        .expect("fixture header builds")
        .with_bsbit_provenance(
            bsbit_hts::BsbitProgramProvenance::new(
                digest
                    .finish()
                    .expect("fixture semantic digest")
                    .into_bytes(),
                bsbit_hts::BsbitAlignmentMode::CallerCompatibleDirectionalPaired,
            ),
            limits,
        )
        .expect("fixture provenance fits")
        .with_sort_order(SamSortOrder::Coordinate);
        let staging = directory.join("fixture.bam.tmp");
        let bam = directory.join("fixture.bam");
        let mut writer =
            BamStagingWriter::create_new(&staging, &header, limits).expect("fixture BAM opens");
        let qualities = vec![b'I'; observed.len()];
        let cigar = [
            AlignmentCigarRun::new(AlignmentCigarOp::Match, reference_length)
                .expect("fixture CIGAR run"),
        ];
        for ordinal in 0..8 {
            let query_name = format!("read-{ordinal:02}");
            let record = BorrowedAlignmentRecord::new(
                query_name.as_bytes(),
                0,
                Some(0),
                1,
                60,
                &cigar,
                None,
                0,
                0,
                &observed,
                Some(&qualities),
                1,
                AlignmentAuxiliaryMode::Minimal,
                None,
                BisulfiteStrand::OT,
                None,
                limits,
            )
            .expect("fixture BAM record builds");
            writer
                .write_borrowed_alignment_record(&record)
                .expect("fixture BAM record writes");
        }
        writer
            .finish()
            .expect("fixture BAM finishes")
            .publish_create_new(&bam)
            .expect("fixture BAM publishes");
        build_bam_index_create_new(&bam, bam.with_extension("bam.bai"), 1)
            .expect("fixture BAI builds");

        IndexedXgFixture {
            directory,
            bam,
            fasta,
            references: vec![BamReference {
                name: b"chr1".to_vec(),
                length: u32::try_from(INDEXED_FIXTURE_REFERENCE.len())
                    .expect("fixture length fits u32"),
            }],
        }
    }

    #[test]
    fn renderers_use_cgmap_and_extended_bedmethyl_coordinates() {
        let key = SiteKey {
            reference: 0,
            position: 2,
            strand: EvidenceStrand::Top,
        };
        let counts = SiteCounts {
            context: Some(CytosineContext {
                class: ContextClass::Chg,
                second: b'A',
            }),
            methylated: 1,
            unmethylated: 1,
            deleted: 2,
            different: 3,
        };
        let mut cgmap = Vec::new();
        render_cgmap(
            &mut cgmap,
            b"chr1",
            key,
            &counts,
            counts.context.unwrap(),
            2,
        )
        .unwrap();
        assert_eq!(cgmap, b"chr1\tC\t3\tCHG\tCA\t0.500000\t1\t2\n");

        let mut bed = Vec::new();
        render_bed(&mut bed, b"chr1", key, &counts, counts.context.unwrap(), 2).unwrap();
        assert_eq!(
            bed,
            b"chr1\t2\t3\tm,CHG,0\t2\t+\t2\t3\t255,0,0\t2\t50.00\t1\t1\t0\t2\t0\t3\t0\n"
        );
    }

    #[test]
    fn region_and_likelihood_sizes_adapt_to_parallel_memory_budget() {
        let config = SnpConfig::default();
        let meth_parameters = MethParameters::default();
        assert_eq!(
            region_bases_for(IndexedCallMode::Meth(meth_parameters), 1),
            MAX_REGION_BASES
        );
        assert_eq!(
            region_bases_for(IndexedCallMode::Joint(config), 1),
            MAX_REGION_BASES
        );

        let meth = region_bases_for(IndexedCallMode::Meth(meth_parameters), 64);
        let snp = region_bases_for(IndexedCallMode::Snp(config), 64);
        let joint = region_bases_for(IndexedCallMode::Joint(config), 64);
        assert!((MIN_REGION_BASES..=MAX_REGION_BASES).contains(&joint));
        assert!(joint < snp && snp < meth);
        for (mode, bases) in [
            (IndexedCallMode::Meth(meth_parameters), meth),
            (IndexedCallMode::Snp(config), snp),
            (IndexedCallMode::Joint(config), joint),
        ] {
            let blocks = usize::try_from(bases / u64::BITS).unwrap();
            assert!(
                blocks * mode.region_bytes_per_block() * 64 * REGION_REORDER_FACTOR
                    <= REGION_STATE_BUDGET_BYTES
            );
        }

        assert_eq!(
            likelihood_batch_sites_for(1),
            MAX_SNP_LIKELIHOOD_BATCH_SITES
        );
        assert!(likelihood_batch_sites_for(64) < MAX_SNP_LIKELIHOOD_BATCH_SITES);
        assert!(likelihood_batch_sites_for(64) >= MIN_SNP_LIKELIHOOD_BATCH_SITES);
    }

    #[test]
    fn dense_bit_sliced_region_matches_scalar_counts_past_u8_coverage() {
        let context = Some(CytosineContext {
            class: ContextClass::Cg,
            second: b'G',
        });
        let mut dense = DenseMethRegion::new(0, 100, 164).unwrap();
        let mut scalar = HashMap::new();
        for round in 0..700_u16 {
            let call = match round {
                0..300 => CallKind::Methylated,
                300..500 => CallKind::Unmethylated,
                500..600 => CallKind::Deleted,
                _ => CallKind::Different,
            };
            for position in [100, 101, 163] {
                let key = SiteKey {
                    reference: 0,
                    position,
                    strand: EvidenceStrand::Top,
                };
                dense.add_observation(key, context, call).unwrap();
                let counts: &mut SiteCounts = scalar.entry(key).or_default();
                counts.context =
                    merge_contexts(key.reference, key.position, counts.context, context).unwrap();
                match call {
                    CallKind::Methylated => counts.methylated += 1,
                    CallKind::Unmethylated => counts.unmethylated += 1,
                    CallKind::Deleted => counts.deleted += 1,
                    CallKind::Different => counts.different += 1,
                }
            }
        }
        let dense = dense.into_sites().unwrap();
        let mut scalar = scalar.into_iter().collect::<Vec<_>>();
        scalar.sort_unstable_by_key(|(key, _)| *key);
        assert_eq!(dense, scalar);
        assert_eq!(dense[0].1.methylated, 300);
        assert_eq!(dense[0].1.different, 100);
    }

    #[test]
    fn classified_word_updates_match_scalar_across_region_blocks() {
        let context = Some(CytosineContext {
            class: ContextClass::Cg,
            second: b'G',
        });
        let observations = (100..180_u32)
            .map(|position| EvidenceObservation {
                reference: 0,
                position,
                reference_base: b'C',
                query_base: match position % 4 {
                    0 => Some(b'C'),
                    1 => Some(b'T'),
                    2 => Some(b'A'),
                    _ => None,
                },
                base_quality: Some(30),
                mapping_quality: 60,
                strand: EvidenceStrand::Top,
                context,
            })
            .collect::<Vec<_>>();
        let mut dense = DenseMethRegion::new(0, 100, 180).unwrap();
        accumulate_meth_fragment(&observations, &mut dense, MethParameters::default()).unwrap();
        let dense = dense.into_sites().unwrap();
        assert_eq!(dense.len(), observations.len());
        for (key, counts) in dense {
            let expected = match key.position % 4 {
                0 => (1, 0, 0, 0),
                1 => (0, 1, 0, 0),
                2 => (0, 0, 0, 1),
                _ => (0, 0, 1, 0),
            };
            assert_eq!(
                (
                    counts.methylated,
                    counts.unmethylated,
                    counts.deleted,
                    counts.different,
                ),
                expected
            );
        }
    }

    #[test]
    fn methylation_quality_thresholds_filter_before_counting() {
        let context = Some(CytosineContext {
            class: ContextClass::Cg,
            second: b'G',
        });
        let observation =
            |position, query_base, base_quality, mapping_quality| EvidenceObservation {
                reference: 0,
                position,
                reference_base: b'C',
                query_base,
                base_quality,
                mapping_quality,
                strand: EvidenceStrand::Top,
                context,
            };
        let observations = vec![
            observation(100, Some(b'C'), Some(30), 60),
            observation(101, Some(b'T'), Some(14), 60),
            observation(102, Some(b'C'), Some(30), 19),
            observation(103, Some(b'C'), Some(30), u8::MAX),
            observation(104, None, None, 60),
            observation(105, Some(b'C'), None, 60),
        ];
        let mut dense = DenseMethRegion::new(0, 100, 106).unwrap();
        accumulate_meth_fragment(&observations, &mut dense, MethParameters::default()).unwrap();
        let sites = dense.into_sites().unwrap();

        assert_eq!(sites.len(), 2);
        assert_eq!(sites[0].0.position, 100);
        assert_eq!(sites[0].1.methylated, 1);
        assert_eq!(sites[1].0.position, 104);
        assert_eq!(sites[1].1.deleted, 1);
    }

    #[test]
    fn indexed_region_workers_are_thread_count_and_boundary_invariant() {
        let fixture = indexed_xg_fixture("indexed-workers");
        let boundary = fixture.references[0].length / 2;
        let regions = [
            CallRegion {
                ordinal: 0,
                reference: 0,
                start: 0,
                end: boundary,
            },
            CallRegion {
                ordinal: 1,
                reference: 0,
                start: boundary,
                end: fixture.references[0].length,
            },
        ];
        let single = run_indexed_region_workers(
            &fixture.bam,
            &fixture.fasta,
            &fixture.references,
            &regions,
            1,
        )
        .unwrap();
        let parallel = run_indexed_region_workers(
            &fixture.bam,
            &fixture.fasta,
            &fixture.references,
            &regions,
            2,
        )
        .unwrap();
        assert_eq!(parallel, single);
        assert!(parallel.windows(2).all(|rows| rows[0].0 < rows[1].0));
        assert!(
            parallel
                .iter()
                .all(|(key, _)| key.position < fixture.references[0].length)
        );
    }

    #[test]
    fn indexed_region_worker_panic_is_returned_as_an_error() {
        let path =
            Path::new(env!("CARGO_MANIFEST_DIR")).join("../../external/htslib/test/range.bam");
        let reference_path =
            Path::new(env!("CARGO_MANIFEST_DIR")).join("../../external/htslib/test/ce.fa");
        let reader = IndexedBamReader::open(&path).unwrap();
        let references = reader
            .header()
            .references()
            .iter()
            .map(|reference| BamReference {
                name: reference.name().to_vec(),
                length: u32::try_from(reference.length()).unwrap(),
            })
            .collect::<Vec<_>>();
        reader.close().unwrap();
        let regions = [CallRegion {
            ordinal: 0,
            reference: 0,
            start: 900,
            end: 1_000,
        }];
        let error = collect_indexed_region_workers_mode(
            &path,
            &reference_path,
            &references,
            &regions,
            1,
            IndexedCallMode::Panic,
        )
        .unwrap_err();
        assert!(error.to_string().contains("worker 0 panicked"));
        assert!(error.to_string().contains("injected regional worker panic"));
    }

    #[test]
    fn joint_regions_equal_separate_meth_and_snp_modules() {
        let fixture = indexed_xg_fixture("joint-workers");
        let boundary = fixture.references[0].length / 2;
        let regions = [
            CallRegion {
                ordinal: 0,
                reference: 0,
                start: 0,
                end: boundary,
            },
            CallRegion {
                ordinal: 1,
                reference: 0,
                start: boundary,
                end: fixture.references[0].length,
            },
        ];
        let config = SnpConfig {
            minimum_base_quality: 0,
            minimum_mapping_quality: 0,
            minimum_depth: 1,
            minimum_alternate_count: 1,
            minimum_genotype_quality: 0,
            ..SnpConfig::default()
        };
        let meth = collect_indexed_region_workers_mode(
            &fixture.bam,
            &fixture.fasta,
            &fixture.references,
            &regions,
            2,
            IndexedCallMode::Meth(MethParameters {
                minimum_base_quality: config.minimum_base_quality,
                minimum_mapping_quality: config.minimum_mapping_quality,
            }),
        )
        .unwrap();
        let snp = collect_indexed_region_workers_mode(
            &fixture.bam,
            &fixture.fasta,
            &fixture.references,
            &regions,
            2,
            IndexedCallMode::Snp(config),
        )
        .unwrap();
        let joint = collect_indexed_region_workers_mode(
            &fixture.bam,
            &fixture.fasta,
            &fixture.references,
            &regions,
            2,
            IndexedCallMode::Joint(config),
        )
        .unwrap();
        assert_eq!(joint.meth_sites, meth.meth_sites);
        assert_eq!(joint.variants, snp.variants);
    }

    #[test]
    fn streaming_methylation_output_is_bgzf_and_create_only() {
        let output = unique_path("output.cgmap.gz");
        let references = vec![BamReference {
            name: b"chr1".to_vec(),
            length: 100,
        }];
        let key = SiteKey {
            reference: 0,
            position: 2,
            strand: EvidenceStrand::Top,
        };
        let context = Some(CytosineContext {
            class: ContextClass::Chg,
            second: b'A',
        });
        let mut region = DenseMethRegion::new(0, 0, 100).unwrap();
        for call in [
            CallKind::Methylated,
            CallKind::Unmethylated,
            CallKind::Deleted,
            CallKind::Deleted,
            CallKind::Different,
            CallKind::Different,
            CallKind::Different,
        ] {
            region.add_observation(key, context, call).unwrap();
        }
        let mut writer =
            TextStagingWriter::create_sibling(&output, TextOutputCompression::Bgzf, 1).unwrap();
        render_meth_region(
            &mut writer,
            MethylationOutputFormat::Cgmap,
            &references,
            &region,
            &mut UnresolvedContextSummary::default(),
        )
        .unwrap();
        writer.finish().unwrap().publish_create_new().unwrap();

        let mut reader = DecodedReader::open(&output).unwrap();
        assert_eq!(reader.compression(), bsbit_hts::Compression::Bgzf);
        let mut decoded = Vec::new();
        reader.read_to_end(&mut decoded).unwrap();
        reader.close().unwrap();
        assert_eq!(decoded, b"chr1\tC\t3\tCHG\tCA\t0.500000\t1\t2\n");
        let published = fs::read(&output).unwrap();
        assert!(
            TextStagingWriter::create_sibling(&output, TextOutputCompression::Bgzf, 1).is_err()
        );
        assert_eq!(fs::read(&output).unwrap(), published);
        fs::remove_file(output).unwrap();
    }

    #[test]
    fn compressed_vcf_is_bgzf_and_decodes_as_vcf() {
        let output = unique_path("output.vcf.gz");
        let references = vec![BamReference {
            name: b"chr1".to_vec(),
            length: 100,
        }];
        let mut writer =
            TextStagingWriter::create_sibling(&output, TextOutputCompression::Bgzf, 1).unwrap();
        render_streaming_vcf_header(&mut writer, &references, SnpConfig::default(), b"sample")
            .unwrap();
        writer.finish().unwrap().publish_create_new().unwrap();

        let mut reader = DecodedReader::open(&output).unwrap();
        assert_eq!(reader.compression(), bsbit_hts::Compression::Bgzf);
        let mut decoded = String::new();
        reader.read_to_string(&mut decoded).unwrap();
        reader.close().unwrap();
        assert!(decoded.starts_with("##fileformat=VCFv4.3\n"));
        assert!(
            decoded.ends_with("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n")
        );
        fs::remove_file(output).unwrap();
    }

    #[test]
    fn candidate_batches_query_only_their_covered_span() {
        let parent = CallRegion {
            ordinal: 3,
            reference: 2,
            start: 100,
            end: 1_000,
        };
        let candidates = [
            crate::snp::candidate::CandidateSite::for_test(120, b'C'),
            crate::snp::candidate::CandidateSite::for_test(450, b'G'),
        ];
        let query = candidate_query_region(parent, &candidates).unwrap();
        assert_eq!(query.reference, 2);
        assert_eq!((query.start, query.end), (120, 451));
    }
}
