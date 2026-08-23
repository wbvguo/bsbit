//! Bounded streaming merge implementation for sorted extended bedMethyl inputs.
//!
//! The implementation performs a bounded-memory hierarchical k-way merge.
//! Input workers retain one record per sample and the ordered coordinator
//! retains one row per worker, so resident methylation state is proportional
//! to the number of samples rather than the number of genomic sites.

#![forbid(unsafe_code)]

use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashSet};
use std::sync::mpsc::{self, Receiver, SyncSender};
use std::thread;

use bsbit_io::validate_distinct_paths;

use crate::input::{BedCursor, ContigCatalog, preflight_catalog};
use crate::output::{
    MatrixOutput, OutputSpec, create_outputs, finish_outputs, output_error, output_specs,
    publish_outputs, write_header, write_matrix_row,
};
use crate::request::{Input, MAX_THREADS, Options};
use crate::result::{CombineError, CombineErrorKind, CombineReport};
use crate::site::{Counts, SiteKey};

const PROPORTION_SCALE: u64 = 1_000_000_000;

/// Combines sorted bsbit extended bedMethyl files into wide matrices.
///
/// Inputs are decoded by content, so ordinary text, gzip, and BGZF are
/// accepted regardless of filename suffix. The output is staged beside its
/// absent destination(s), finalized, synchronized, and published create-only.
///
/// # Errors
///
/// Returns an error for invalid options, malformed or inconsistently ordered
/// input, worker failure, output encoding failure, or publication failure.
pub fn combine(options: &Options) -> Result<CombineReport, CombineError> {
    let output_specs = validate_options(options)?;
    let mut outputs = create_outputs(options, &output_specs)?;
    let contigs = preflight_catalog(&options.inputs, options.threads)?;
    for output in &mut outputs {
        write_header(&mut output.writer, options, output.spec.kind)
            .map_err(|error| output_error(&output.spec.path, error))?;
    }

    let mut report = run_parallel_merge(options, &contigs, &mut outputs)?;
    let completed = finish_outputs(outputs)?;
    publish_outputs(completed, &mut report)?;
    Ok(report)
}

fn validate_options(options: &Options) -> Result<Vec<OutputSpec>, CombineError> {
    if options.inputs.is_empty() {
        return Err(CombineError::configuration(
            "combine: at least one methylation input is required",
        ));
    }
    if !(1..=MAX_THREADS).contains(&options.threads) {
        return Err(CombineError::configuration(format!(
            "combine: thread count must be within 1..={MAX_THREADS}"
        )));
    }
    if u64::from(
        options
            .parameters
            .minimum_sample_proportion_parts_per_billion,
    ) > PROPORTION_SCALE
    {
        return Err(CombineError::configuration(
            "combine: minimum sample proportion must be within 0..=1",
        ));
    }
    if options.output.as_os_str().is_empty() {
        return Err(CombineError::configuration(
            "combine: output path must not be empty",
        ));
    }

    let output_specs = output_specs(options)?;
    for pair in output_specs.windows(2) {
        validate_distinct_paths(&pair[0].path, &pair[1].path).map_err(|error| {
            CombineError::with_source(
                CombineErrorKind::Configuration,
                "combine: derived output paths must differ",
                error,
            )
        })?;
    }

    let mut samples = HashSet::with_capacity(options.inputs.len());
    let mut paths = HashSet::with_capacity(options.inputs.len());
    for input in &options.inputs {
        if input.sample.is_empty() || input.sample.bytes().any(|byte| byte.is_ascii_control()) {
            return Err(CombineError::configuration(format!(
                "combine: invalid sample label `{}`; labels must be nonempty and contain no control bytes",
                input.sample
            )));
        }
        if !samples.insert(input.sample.as_str()) {
            return Err(CombineError::configuration(format!(
                "combine: duplicate sample label `{}`",
                input.sample
            )));
        }
        if input.path.as_os_str().is_empty() {
            return Err(CombineError::configuration(format!(
                "combine: input path for sample `{}` must not be empty",
                input.sample
            )));
        }
        if !paths.insert(input.path.as_path()) {
            return Err(CombineError::configuration(format!(
                "combine: duplicate input path {}",
                input.path.display()
            )));
        }
        for output in &output_specs {
            validate_distinct_paths(&input.path, &output.path).map_err(|error| {
                CombineError::with_source(
                    CombineErrorKind::Configuration,
                    format!(
                        "combine: input {} and output {} must be different paths",
                        input.path.display(),
                        output.path.display()
                    ),
                    error,
                )
            })?;
        }
    }
    Ok(output_specs)
}

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
struct StreamHeapItem {
    key: SiteKey,
    index: usize,
}

#[derive(Debug)]
struct GroupRow {
    key: SiteKey,
    modification: Vec<u8>,
    values: Vec<(usize, Counts)>,
}

enum WorkerMessage {
    Row(GroupRow),
    Error(CombineError),
    Done,
}

fn merge_group(
    inputs: &[Input],
    global_start: usize,
    contigs: &ContigCatalog,
    sender: &SyncSender<WorkerMessage>,
) -> Result<(), CombineError> {
    let mut cursors = inputs
        .iter()
        .map(|input| BedCursor::open(input, contigs))
        .collect::<Result<Vec<_>, _>>()?;
    let mut heads = (0..inputs.len()).map(|_| None).collect::<Vec<_>>();
    let mut heap = BinaryHeap::new();
    for (index, cursor) in cursors.iter_mut().enumerate() {
        if let Some(record) = cursor.next_record()? {
            heap.push(Reverse(StreamHeapItem {
                key: record.key,
                index,
            }));
            heads[index] = Some(record);
        }
    }

    while let Some(Reverse(first)) = heap.pop() {
        let key = first.key;
        let mut matching = vec![first.index];
        while heap.peek().is_some_and(|Reverse(item)| item.key == key) {
            let Reverse(item) = heap.pop().expect("peeked group heap item");
            matching.push(item.index);
        }

        let first_record = heads[matching[0]]
            .as_ref()
            .expect("group heap points to a record");
        let modification = first_record.modification.clone();
        let mut values = Vec::with_capacity(matching.len());
        for index in matching {
            let record = heads[index]
                .take()
                .expect("group heap points to an owned record");
            if record.modification != modification {
                return Err(metadata_mismatch(contigs, key));
            }
            values.push((global_start + index, record.counts));
            if let Some(next) = cursors[index].next_record()? {
                heap.push(Reverse(StreamHeapItem {
                    key: next.key,
                    index,
                }));
                heads[index] = Some(next);
            }
        }
        if sender
            .send(WorkerMessage::Row(GroupRow {
                key,
                modification,
                values,
            }))
            .is_err()
        {
            return Ok(());
        }
    }
    Ok(())
}

fn metadata_mismatch(contigs: &ContigCatalog, key: SiteKey) -> CombineError {
    let name = contigs
        .names
        .get(usize::try_from(key.contig).expect("u32 fits usize"))
        .map_or_else(|| "<unknown>".into(), |name| String::from_utf8_lossy(name));
    CombineError::input(format!(
        "combine: inputs disagree on modification/context at {name}:{}-{} strand {}",
        key.start,
        key.end,
        if key.strand == 0 { '+' } else { '-' }
    ))
}

fn run_parallel_merge(
    options: &Options,
    contigs: &ContigCatalog,
    outputs: &mut [MatrixOutput],
) -> Result<CombineReport, CombineError> {
    let worker_count = usize::try_from(options.threads)
        .expect("validated thread count fits usize")
        .min(options.inputs.len());
    thread::scope(|scope| {
        let mut handles = Vec::with_capacity(worker_count);
        let mut receivers = Vec::with_capacity(worker_count);
        for worker in 0..worker_count {
            let start = options.inputs.len() * worker / worker_count;
            let end = options.inputs.len() * (worker + 1) / worker_count;
            let (sender, receiver) = mpsc::sync_channel(1);
            receivers.push(receiver);
            handles.push(scope.spawn(move || {
                let result = merge_group(&options.inputs[start..end], start, contigs, &sender);
                let terminal = match result {
                    Ok(()) => WorkerMessage::Done,
                    Err(error) => WorkerMessage::Error(error),
                };
                let _ = sender.send(terminal);
            }));
        }

        let result = coordinate_groups(options, contigs, outputs, &receivers);
        drop(receivers);
        let mut worker_panic = false;
        for handle in handles {
            worker_panic |= handle.join().is_err();
        }
        if worker_panic {
            Err(CombineError::worker("combine: input merge worker panicked"))
        } else {
            result
        }
    })
}

fn coordinate_groups(
    options: &Options,
    contigs: &ContigCatalog,
    outputs: &mut [MatrixOutput],
    receivers: &[Receiver<WorkerMessage>],
) -> Result<CombineReport, CombineError> {
    let mut heads = (0..receivers.len()).map(|_| None).collect::<Vec<_>>();
    let mut heap = BinaryHeap::new();
    for index in 0..receivers.len() {
        receive_group(index, receivers, &mut heads, &mut heap)?;
    }
    let required_samples = required_sample_count(options)?;
    let mut report = CombineReport::default();

    while let Some(Reverse(first)) = heap.pop() {
        let key = first.key;
        let mut matching = vec![first.index];
        while heap.peek().is_some_and(|Reverse(item)| item.key == key) {
            let Reverse(item) = heap.pop().expect("peeked coordinator heap item");
            matching.push(item.index);
        }
        report.sites_seen = report
            .sites_seen
            .checked_add(1)
            .ok_or_else(|| CombineError::input("combine: distinct site count overflowed u64"))?;

        let first_row = heads[matching[0]]
            .as_ref()
            .expect("coordinator heap points to a row");
        let modification = first_row.modification.clone();
        let present_count = matching
            .iter()
            .map(|index| {
                heads[*index]
                    .as_ref()
                    .expect("coordinator heap points to a row")
                    .values
                    .len()
            })
            .sum();
        let mut values = Vec::with_capacity(present_count);
        for index in &matching {
            let row = heads[*index]
                .take()
                .expect("coordinator heap points to an owned row");
            if row.modification != modification {
                return Err(metadata_mismatch(contigs, key));
            }
            values.extend(row.values);
        }
        let valid_samples = values
            .iter()
            .filter(|(_, counts)| sample_is_valid(*counts, options.parameters.minimum_count))
            .count();
        if valid_samples >= required_samples {
            for output in outputs.iter_mut() {
                write_matrix_row(
                    &mut output.writer,
                    output.spec.kind,
                    options,
                    contigs,
                    key,
                    &modification,
                    &values,
                )
                .map_err(|error| output_error(&output.spec.path, error))?;
            }
            report.sites_written = report
                .sites_written
                .checked_add(1)
                .ok_or_else(|| CombineError::input("combine: written site count overflowed u64"))?;
        }
        for index in matching {
            receive_group(index, receivers, &mut heads, &mut heap)?;
        }
    }
    Ok(report)
}

fn receive_group(
    index: usize,
    receivers: &[Receiver<WorkerMessage>],
    heads: &mut [Option<GroupRow>],
    heap: &mut BinaryHeap<Reverse<StreamHeapItem>>,
) -> Result<(), CombineError> {
    match receivers[index].recv() {
        Ok(WorkerMessage::Row(row)) => {
            heap.push(Reverse(StreamHeapItem {
                key: row.key,
                index,
            }));
            heads[index] = Some(row);
            Ok(())
        }
        Ok(WorkerMessage::Done) => Ok(()),
        Ok(WorkerMessage::Error(error)) => Err(error),
        Err(_) => Err(CombineError::worker(format!(
            "combine: input merge worker {} disconnected without a result",
            index + 1
        ))),
    }
}

fn required_sample_count(options: &Options) -> Result<usize, CombineError> {
    let numerator = (options.inputs.len() as u128)
        * u128::from(
            options
                .parameters
                .minimum_sample_proportion_parts_per_billion,
        );
    let required = numerator.div_ceil(u128::from(PROPORTION_SCALE)).max(1);
    usize::try_from(required)
        .map_err(|_| CombineError::configuration("combine: required sample count exceeds usize"))
}

const fn sample_is_valid(counts: Counts, minimum_count: u64) -> bool {
    counts.total != 0 && counts.total >= minimum_count
}
