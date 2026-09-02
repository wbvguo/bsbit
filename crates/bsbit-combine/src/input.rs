//! Input decoding, validation, and ordered methylation-record cursors.

#![forbid(unsafe_code)]

use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashMap, HashSet};
use std::io::{BufRead, BufReader};
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering as AtomicOrdering};
use std::thread;

use bsbit_hts::{BedMethylContext, BedMethylRecord, BedMethylStrand, DecodedReader};

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
        let parsed = parse_sample_site(lines.current_line())
            .map_err(|error| input_line_error(input, line_number, error))?;
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

fn parse_sample_site(line: &[u8]) -> Result<SampleSite<'_>, String> {
    let columns = line.split(|byte| *byte == b'\t').count();
    match columns {
        8 => parse_cgmap_site(line),
        18 => BedMethylRecord::parse(line)
            .map(SampleSite::from)
            .map_err(|error| error.to_string()),
        _ => Err(format!(
            "expected exactly 8 CGmap or 18 extended bedMethyl columns, observed {columns}"
        )),
    }
}

fn parse_cgmap_site(line: &[u8]) -> Result<SampleSite<'_>, String> {
    let mut columns = [&[][..]; 8];
    let mut fields = line.split(|byte| *byte == b'\t');
    for column in &mut columns {
        *column = fields
            .next()
            .expect("CGmap column count was validated before parsing");
    }
    debug_assert!(fields.next().is_none());
    if columns[0].is_empty() {
        return Err("CGmap contig in column 1 must not be empty".to_owned());
    }
    let strand = match columns[1] {
        b"C" => 0,
        b"G" => 1,
        _ => return Err("CGmap column 2 must be `C` or `G`".to_owned()),
    };
    let position = parse_cgmap_u64(columns[2], 3)?;
    let start = position
        .checked_sub(1)
        .ok_or_else(|| "CGmap column 3 must be a positive 1-based position".to_owned())?;
    let modification = match columns[3] {
        b"CG" => BedMethylContext::Cg.modification(),
        b"CHG" => BedMethylContext::Chg.modification(),
        b"CHH" => BedMethylContext::Chh.modification(),
        _ => return Err("CGmap column 4 must be `CG`, `CHG`, or `CHH`".to_owned()),
    };
    let valid_dinucleotide = match columns[3] {
        b"CG" => columns[4] == b"CG",
        b"CHG" | b"CHH" => matches!(columns[4], b"CA" | b"CC" | b"CT"),
        _ => false,
    };
    if !valid_dinucleotide {
        return Err("CGmap column 5 is inconsistent with the context in column 4".to_owned());
    }
    let methylated = parse_cgmap_u64(columns[6], 7)?;
    let total = parse_cgmap_u64(columns[7], 8)?;
    if methylated > total {
        return Err(
            "CGmap methylated count in column 7 exceeds total count in column 8".to_owned(),
        );
    }
    validate_cgmap_level(columns[5], total)?;
    Ok(SampleSite {
        contig: columns[0],
        start,
        end: position,
        modification,
        strand,
        methylated,
        total,
    })
}

fn parse_cgmap_u64(value: &[u8], column: u8) -> Result<u64, String> {
    if value.is_empty() {
        return Err(format!(
            "CGmap column {column} must be a nonnegative integer"
        ));
    }
    let mut parsed = 0_u64;
    for &byte in value {
        if !byte.is_ascii_digit() {
            return Err(format!(
                "CGmap column {column} must be a nonnegative integer"
            ));
        }
        parsed = parsed
            .checked_mul(10)
            .and_then(|current| current.checked_add(u64::from(byte - b'0')))
            .ok_or_else(|| format!("CGmap column {column} overflows u64"))?;
    }
    Ok(parsed)
}

fn validate_cgmap_level(value: &[u8], total: u64) -> Result<(), String> {
    if value == b"na" {
        return (total == 0)
            .then_some(())
            .ok_or_else(|| "CGmap column 6 may be `na` only when total count is zero".to_owned());
    }
    if total == 0 {
        return Err("CGmap column 6 must be `na` when total count is zero".to_owned());
    }
    let mut parts = value.split(|byte| *byte == b'.');
    let whole = parts.next().unwrap_or_default();
    let fraction = parts.next();
    if parts.next().is_some()
        || whole.is_empty()
        || !whole.iter().all(u8::is_ascii_digit)
        || fraction
            .is_some_and(|digits| digits.is_empty() || !digits.iter().all(u8::is_ascii_digit))
        || !matches!(whole, b"0" | b"1")
        || (whole == b"1"
            && fraction.is_some_and(|digits| digits.iter().any(|digit| *digit != b'0')))
    {
        return Err("CGmap column 6 must be a decimal within 0..=1 or `na`".to_owned());
    }
    Ok(())
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
pub(crate) struct MethylationRecord {
    pub(crate) key: SiteKey,
    pub(crate) modification: Vec<u8>,
    pub(crate) counts: Counts,
}

pub(crate) struct MethylationCursor<'a> {
    input: &'a Input,
    contigs: &'a ContigCatalog,
    lines: DecodedLines,
    previous: Option<SiteKey>,
}

impl<'a> MethylationCursor<'a> {
    pub(crate) fn open(input: &'a Input, contigs: &'a ContigCatalog) -> Result<Self, CombineError> {
        Ok(Self {
            input,
            contigs,
            lines: DecodedLines::open(input)?,
            previous: None,
        })
    }

    pub(crate) fn next_record(&mut self) -> Result<Option<MethylationRecord>, CombineError> {
        let Some(line_number) = self.lines.next_data_line()? else {
            return Ok(None);
        };
        let parsed = parse_sample_site(self.lines.current_line())
            .map_err(|error| input_line_error(self.input, line_number, error))?;
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
        Ok(Some(MethylationRecord {
            key,
            modification: parsed.modification.to_vec(),
            counts: Counts {
                methylated: parsed.methylated,
                total: parsed.total,
            },
        }))
    }
}
