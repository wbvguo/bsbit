//! Input decoding, validation, and ordered bedMethyl cursors.

#![forbid(unsafe_code)]

use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashMap, HashSet};
use std::io::{BufRead, BufReader};
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering as AtomicOrdering};
use std::thread;

use bsbit_hts::{BedMethylRecord, BedMethylStrand, DecodedReader};

use crate::request::Input;
use crate::result::{CombineError, CombineErrorKind};
use crate::site::{Counts, SiteKey};

const INPUT_BUFFER_BYTES: usize = 64 * 1024;
const MAX_INPUT_LINE_BYTES: usize = 1024 * 1024;

#[derive(Debug)]
struct InputPreflight {
    contigs: Vec<Vec<u8>>,
}

pub(crate) fn preflight_catalog(
    inputs: &[Input],
    threads: u64,
) -> Result<ContigCatalog, CombineError> {
    let preflights = preflight_inputs(inputs, threads)?;
    ContigCatalog::build(&preflights)
}

fn preflight_inputs(inputs: &[Input], threads: u64) -> Result<Vec<InputPreflight>, CombineError> {
    let worker_count = usize::try_from(threads)
        .expect("validated thread count fits usize")
        .min(inputs.len());
    let next = AtomicUsize::new(0);
    let worker_results =
        thread::scope(|scope| {
            let mut handles = Vec::with_capacity(worker_count);
            for _ in 0..worker_count {
                handles.push(scope.spawn(|| {
                    let mut results = Vec::new();
                    loop {
                        let index = next.fetch_add(1, AtomicOrdering::Relaxed);
                        let Some(input) = inputs.get(index) else {
                            break;
                        };
                        results.push((index, preflight_input(input)));
                    }
                    results
                }));
            }
            let mut results = Vec::with_capacity(worker_count);
            for handle in handles {
                results.push(handle.join().map_err(|_| {
                    CombineError::worker("combine: input preflight worker panicked")
                })?);
            }
            Ok::<_, CombineError>(results)
        })?;

    let mut ordered = (0..inputs.len()).map(|_| None).collect::<Vec<_>>();
    for worker in worker_results {
        for (index, result) in worker {
            ordered[index] = Some(result);
        }
    }
    ordered
        .into_iter()
        .enumerate()
        .map(|(index, result)| {
            result.ok_or_else(|| {
                CombineError::worker(format!(
                    "combine: preflight worker omitted input {}",
                    inputs[index].path.display()
                ))
            })?
        })
        .collect()
}

fn preflight_input(input: &Input) -> Result<InputPreflight, CombineError> {
    let mut lines = DecodedLines::open(input)?;
    let mut contigs = Vec::new();
    let mut seen_contigs = HashSet::new();
    let mut active_contig = None::<Vec<u8>>;
    let mut previous_coordinate = None;

    while let Some(line_number) = lines.next_data_line()? {
        let parsed = BedMethylRecord::parse(lines.current_line())
            .map(SampleSite::from)
            .map_err(|error| input_line_error(input, line_number, error.to_string()))?;
        match active_contig.as_deref() {
            Some(contig) if contig == parsed.contig => {
                let coordinate = parsed.coordinate();
                if previous_coordinate.is_some_and(|previous| previous >= coordinate) {
                    return Err(input_line_error(
                        input,
                        line_number,
                        "site coordinates are duplicated or not strictly increasing within the contig",
                    ));
                }
                previous_coordinate = Some(coordinate);
            }
            _ => {
                if seen_contigs.contains(parsed.contig) {
                    return Err(input_line_error(
                        input,
                        line_number,
                        "a contig appears in more than one non-contiguous block",
                    ));
                }
                let contig = parsed.contig.to_vec();
                seen_contigs.insert(contig.clone());
                contigs.push(contig.clone());
                active_contig = Some(contig);
                previous_coordinate = Some(parsed.coordinate());
            }
        }
    }
    Ok(InputPreflight { contigs })
}

fn input_line_error(
    input: &Input,
    line_number: u64,
    message: impl std::fmt::Display,
) -> CombineError {
    CombineError::input(format!(
        "combine: sample `{}` input {} line {line_number}: {message}",
        input.sample,
        input.path.display()
    ))
}

struct DecodedLines {
    sample: String,
    path: PathBuf,
    reader: Option<BufReader<DecodedReader>>,
    buffer: Vec<u8>,
    line_number: u64,
}

impl DecodedLines {
    fn open(input: &Input) -> Result<Self, CombineError> {
        let reader = DecodedReader::open(&input.path).map_err(|error| {
            CombineError::with_source(
                CombineErrorKind::Input,
                format!(
                    "combine: open sample `{}` input {}",
                    input.sample,
                    input.path.display()
                ),
                error,
            )
        })?;
        Ok(Self {
            sample: input.sample.clone(),
            path: input.path.clone(),
            reader: Some(BufReader::with_capacity(INPUT_BUFFER_BYTES, reader)),
            buffer: Vec::with_capacity(256),
            line_number: 0,
        })
    }

    fn next_data_line(&mut self) -> Result<Option<u64>, CombineError> {
        loop {
            self.buffer.clear();
            if !self.read_physical_line()? {
                return Ok(None);
            }
            if self.buffer.last() == Some(&b'\n') {
                self.buffer.pop();
            }
            if self.buffer.last() == Some(&b'\r') {
                self.buffer.pop();
            }
            if self.buffer.is_empty() || self.buffer.first() == Some(&b'#') {
                continue;
            }
            return Ok(Some(self.line_number));
        }
    }

    fn current_line(&self) -> &[u8] {
        &self.buffer
    }

    fn read_physical_line(&mut self) -> Result<bool, CombineError> {
        if self.reader.is_none() {
            return Ok(false);
        }
        let mut read_any = false;
        loop {
            let (available_length, take, complete) = {
                let reader = self.reader.as_mut().ok_or_else(|| {
                    CombineError::input(format!(
                        "combine: sample `{}` input {} was read after end of stream",
                        self.sample,
                        self.path.display()
                    ))
                })?;
                let available = reader.fill_buf().map_err(|error| {
                    CombineError::with_source(
                        CombineErrorKind::Input,
                        format!(
                            "combine: decode sample `{}` input {} near line {}",
                            self.sample,
                            self.path.display(),
                            self.line_number.saturating_add(1)
                        ),
                        error,
                    )
                })?;
                if available.is_empty() {
                    (0, 0, false)
                } else if let Some(index) = available.iter().position(|byte| *byte == b'\n') {
                    (available.len(), index + 1, true)
                } else {
                    (available.len(), available.len(), false)
                }
            };

            if available_length == 0 {
                self.close()?;
                if read_any {
                    self.line_number = self.line_number.checked_add(1).ok_or_else(|| {
                        CombineError::input(format!(
                            "combine: line count overflowed for {}",
                            self.path.display()
                        ))
                    })?;
                }
                return Ok(read_any);
            }
            if self.buffer.len().saturating_add(take) > MAX_INPUT_LINE_BYTES {
                return Err(CombineError::input(format!(
                    "combine: sample `{}` input {} line {} exceeds {MAX_INPUT_LINE_BYTES} bytes",
                    self.sample,
                    self.path.display(),
                    self.line_number.saturating_add(1)
                )));
            }
            {
                let reader = self.reader.as_mut().expect("nonempty decoded reader");
                let available = reader.fill_buf().map_err(|error| {
                    CombineError::with_source(
                        CombineErrorKind::Input,
                        format!("combine: decode input {}", self.path.display()),
                        error,
                    )
                })?;
                self.buffer.extend_from_slice(&available[..take]);
                reader.consume(take);
            }
            read_any = true;
            if complete {
                self.line_number = self.line_number.checked_add(1).ok_or_else(|| {
                    CombineError::input(format!(
                        "combine: line count overflowed for {}",
                        self.path.display()
                    ))
                })?;
                return Ok(true);
            }
        }
    }

    fn close(&mut self) -> Result<(), CombineError> {
        let Some(reader) = self.reader.take() else {
            return Ok(());
        };
        reader.into_inner().close().map_err(|error| {
            CombineError::with_source(
                CombineErrorKind::Input,
                format!(
                    "combine: close sample `{}` input {}",
                    self.sample,
                    self.path.display()
                ),
                error,
            )
        })
    }
}

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
struct Coordinate {
    start: u64,
    end: u64,
    strand: u8,
}

#[derive(Clone, Copy, Debug)]
struct SampleSite<'a> {
    contig: &'a [u8],
    start: u64,
    end: u64,
    modification: &'a [u8],
    strand: u8,
    methylated: u64,
    total: u64,
}

impl SampleSite<'_> {
    const fn coordinate(self) -> Coordinate {
        Coordinate {
            start: self.start,
            end: self.end,
            strand: self.strand,
        }
    }
}

impl<'a> From<BedMethylRecord<'a>> for SampleSite<'a> {
    fn from(record: BedMethylRecord<'a>) -> Self {
        let strand = match record.strand() {
            BedMethylStrand::Forward => 0,
            BedMethylStrand::Reverse => 1,
        };
        Self {
            contig: record.contig(),
            start: record.start(),
            end: record.end(),
            modification: record.context().modification(),
            strand,
            methylated: record.methylated(),
            total: record.coverage(),
        }
    }
}

#[derive(Debug)]
pub(crate) struct ContigCatalog {
    pub(crate) names: Vec<Vec<u8>>,
    ranks: HashMap<Vec<u8>, u32>,
}

impl ContigCatalog {
    fn build(preflights: &[InputPreflight]) -> Result<Self, CombineError> {
        let mut node_names = Vec::<Vec<u8>>::new();
        let mut nodes = HashMap::<Vec<u8>, usize>::new();
        for preflight in preflights {
            for contig in &preflight.contigs {
                if !nodes.contains_key(contig.as_slice()) {
                    let index = node_names.len();
                    nodes.insert(contig.clone(), index);
                    node_names.push(contig.clone());
                }
            }
        }

        let mut outgoing = (0..node_names.len())
            .map(|_| HashSet::<usize>::new())
            .collect::<Vec<_>>();
        let mut indegree = vec![0_usize; node_names.len()];
        for preflight in preflights {
            for pair in preflight.contigs.windows(2) {
                let first = nodes[pair[0].as_slice()];
                let second = nodes[pair[1].as_slice()];
                if outgoing[first].insert(second) {
                    indegree[second] = indegree[second].checked_add(1).ok_or_else(|| {
                        CombineError::input("combine: contig-order indegree overflowed")
                    })?;
                }
            }
        }

        let mut ready = BinaryHeap::new();
        for (index, degree) in indegree.iter().enumerate() {
            if *degree == 0 {
                ready.push(Reverse(index));
            }
        }
        let mut ordered_nodes = Vec::with_capacity(node_names.len());
        while let Some(Reverse(index)) = ready.pop() {
            ordered_nodes.push(index);
            let mut successors = outgoing[index].iter().copied().collect::<Vec<_>>();
            successors.sort_unstable();
            for successor in successors {
                indegree[successor] -= 1;
                if indegree[successor] == 0 {
                    ready.push(Reverse(successor));
                }
            }
        }
        if ordered_nodes.len() != node_names.len() {
            let conflicting = indegree
                .iter()
                .enumerate()
                .filter(|(_, degree)| **degree != 0)
                .take(8)
                .map(|(index, _)| String::from_utf8_lossy(&node_names[index]).into_owned())
                .collect::<Vec<_>>()
                .join(", ");
            return Err(CombineError::input(format!(
                "combine: inputs disagree on contig order; ordering cycle includes {conflicting}"
            )));
        }

        let mut names = Vec::with_capacity(node_names.len());
        let mut ranks = HashMap::with_capacity(node_names.len());
        for (rank, node) in ordered_nodes.into_iter().enumerate() {
            let rank = u32::try_from(rank)
                .map_err(|_| CombineError::input("combine: more than u32::MAX contigs"))?;
            let name = node_names[node].clone();
            ranks.insert(name.clone(), rank);
            names.push(name);
        }
        Ok(Self { names, ranks })
    }
}

#[derive(Debug)]
pub(crate) struct BedRecord {
    pub(crate) key: SiteKey,
    pub(crate) modification: Vec<u8>,
    pub(crate) counts: Counts,
}

pub(crate) struct BedCursor<'a> {
    input: &'a Input,
    contigs: &'a ContigCatalog,
    lines: DecodedLines,
    previous: Option<SiteKey>,
}

impl<'a> BedCursor<'a> {
    pub(crate) fn open(input: &'a Input, contigs: &'a ContigCatalog) -> Result<Self, CombineError> {
        Ok(Self {
            input,
            contigs,
            lines: DecodedLines::open(input)?,
            previous: None,
        })
    }

    pub(crate) fn next_record(&mut self) -> Result<Option<BedRecord>, CombineError> {
        let Some(line_number) = self.lines.next_data_line()? else {
            return Ok(None);
        };
        let parsed = BedMethylRecord::parse(self.lines.current_line())
            .map(SampleSite::from)
            .map_err(|error| input_line_error(self.input, line_number, error.to_string()))?;
        let contig = self
            .contigs
            .ranks
            .get(parsed.contig)
            .copied()
            .ok_or_else(|| {
                input_line_error(
                    self.input,
                    line_number,
                    "contig was not present during input preflight; input changed while combine was running",
                )
            })?;
        let key = SiteKey {
            contig,
            start: parsed.start,
            end: parsed.end,
            strand: parsed.strand,
        };
        if self.previous.is_some_and(|previous| previous >= key) {
            return Err(input_line_error(
                self.input,
                line_number,
                "site coordinates are duplicated or not strictly increasing in the shared contig order",
            ));
        }
        self.previous = Some(key);
        Ok(Some(BedRecord {
            key,
            modification: parsed.modification.to_vec(),
            counts: Counts {
                methylated: parsed.methylated,
                total: parsed.total,
            },
        }))
    }
}
