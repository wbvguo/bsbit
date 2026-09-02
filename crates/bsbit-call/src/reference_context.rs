//! Bounded FASTA context lookup for methylation and variant calling.

use std::ffi::OsString;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use bsbit_core::reference::{ReferenceSemanticDigest, ReferenceSemanticDigestBuilder};
use bsbit_hts::{Compression, DecodedReader, IndexedFastaReader};

use crate::call_input::BamReference;
use crate::evidence::{ContextClass, CytosineContext, EvidenceStrand};
use crate::region::CallRegion;
use crate::{CallError, CallErrorKind};

pub(crate) struct CallReferenceSource {
    path: PathBuf,
    backend: CallReferenceBackend,
    fasta_reference_ids: Arc<[u32]>,
}

enum CallReferenceBackend {
    Indexed,
    Scanned(Arc<ScannedFastaIndex>),
}

impl CallReferenceSource {
    pub(crate) fn prepare(path: &Path, references: &[BamReference]) -> Result<Self, CallError> {
        let path = path.to_path_buf();
        let compression = detect_reference_compression(&path)?;
        let (backend, fasta_references) = match compression {
            Compression::Gzip => return Err(unsupported_gzip_reference(&path)),
            Compression::Bgzf => {
                require_bgzf_fasta_indexes(&path)?;
                indexed_reference_metadata(&path)?
            }
            Compression::Plain => prepare_plain_reference(&path)?,
        };
        let fasta_reference_ids = map_reference_dictionary(&path, references, &fasta_references)?;
        Ok(Self {
            path,
            backend,
            fasta_reference_ids: fasta_reference_ids.into(),
        })
    }

    pub(crate) fn open(&self) -> Result<CallReferenceReader, CallError> {
        let reader = match &self.backend {
            CallReferenceBackend::Indexed => {
                ReferenceReaderBackend::Indexed(open_indexed_reader(&self.path)?)
            }
            CallReferenceBackend::Scanned(index) => ReferenceReaderBackend::Scanned(
                ScannedFastaReader::open(&self.path, Arc::clone(index))?,
            ),
        };
        Ok(CallReferenceReader {
            reader,
            fasta_reference_ids: Arc::clone(&self.fasta_reference_ids),
        })
    }
}

pub(crate) struct CallReferenceReader {
    reader: ReferenceReaderBackend,
    fasta_reference_ids: Arc<[u32]>,
}

enum ReferenceReaderBackend {
    Indexed(IndexedFastaReader),
    Scanned(ScannedFastaReader),
}

impl CallReferenceReader {
    fn fetch(&mut self, reference_id: u32, start: u64, end: u64) -> Result<Vec<u8>, CallError> {
        match &mut self.reader {
            ReferenceReaderBackend::Indexed(reader) => {
                reader.fetch(reference_id, start, end).map_err(|error| {
                    CallError::with_source(
                        CallErrorKind::Input,
                        "fetch indexed reference FASTA",
                        error,
                    )
                })
            }
            ReferenceReaderBackend::Scanned(reader) => reader.fetch(reference_id, start, end),
        }
    }

    pub(crate) fn fetch_context_window(
        &mut self,
        region: CallRegion,
        references: &[BamReference],
    ) -> Result<ReferenceWindow, CallError> {
        let reference = references
            .get(usize::try_from(region.reference).expect("u32 fits usize"))
            .ok_or_else(|| CallError::operation("calling region references a missing contig"))?;
        let fasta_reference_id = *self
            .fasta_reference_ids
            .get(usize::try_from(region.reference).expect("u32 fits usize"))
            .ok_or_else(|| CallError::operation("reference dictionary mapping is incomplete"))?;
        let start = region.start.saturating_sub(2);
        let end = region.end.saturating_add(2).min(reference.length);
        let mut bases = self
            .fetch(fasta_reference_id, u64::from(start), u64::from(end))
            .map_err(|error| {
                error.with_context(format!(
                    "fetch reference FASTA contig `{}` interval {}-{}",
                    String::from_utf8_lossy(&reference.name),
                    start,
                    end
                ))
            })?;
        bases.make_ascii_uppercase();
        Ok(ReferenceWindow {
            reference: region.reference,
            start,
            bases,
        })
    }

    pub(crate) fn validate_semantic_digest(
        &mut self,
        references: &[BamReference],
        expected: ReferenceSemanticDigest,
    ) -> Result<(), CallError> {
        let contig_count = u64::try_from(references.len())
            .map_err(|_| CallError::input("BAM reference count exceeds u64"))?;
        let mut builder = ReferenceSemanticDigestBuilder::new(contig_count);
        for (ordinal, reference) in references.iter().enumerate() {
            let fasta_reference_id = *self.fasta_reference_ids.get(ordinal).ok_or_else(|| {
                CallError::operation("reference dictionary mapping is incomplete")
            })?;
            builder
                .begin_ascii_contig(&reference.name, u64::from(reference.length))
                .map_err(|error| {
                    CallError::with_source(
                        CallErrorKind::Input,
                        format!(
                            "normalize reference FASTA contig `{}` for provenance validation",
                            String::from_utf8_lossy(&reference.name)
                        ),
                        error,
                    )
                })?;
            let mut start = 0_u64;
            let length = u64::from(reference.length);
            while start < length {
                let end = start.saturating_add(8 * 1024 * 1024).min(length);
                let bases = self.fetch(fasta_reference_id, start, end).map_err(|error| {
                    error.with_context(format!(
                        "read reference FASTA contig `{}` interval {start}-{end} for provenance validation",
                        String::from_utf8_lossy(&reference.name)
                    ))
                })?;
                builder.push_ascii_bases(&bases).map_err(|error| {
                    CallError::with_source(
                        CallErrorKind::Input,
                        format!(
                            "normalize reference FASTA contig `{}` for provenance validation",
                            String::from_utf8_lossy(&reference.name)
                        ),
                        error,
                    )
                })?;
                start = end;
            }
            builder.end_ascii_contig().map_err(|error| {
                CallError::with_source(
                    CallErrorKind::Input,
                    format!(
                        "finish reference FASTA contig `{}` provenance validation",
                        String::from_utf8_lossy(&reference.name)
                    ),
                    error,
                )
            })?;
        }
        let observed = builder.finish().map_err(|error| {
            CallError::with_source(
                CallErrorKind::Input,
                "finish reference FASTA provenance validation",
                error,
            )
        })?;
        if observed != expected {
            return Err(CallError::input(format!(
                "reference FASTA semantic digest {observed} differs from BAM provenance {expected}"
            )));
        }
        Ok(())
    }

    pub(crate) fn close(self) -> Result<(), CallError> {
        match self.reader {
            ReferenceReaderBackend::Indexed(reader) => reader.close().map_err(|error| {
                CallError::with_source(CallErrorKind::Input, "close indexed reference FASTA", error)
            }),
            ReferenceReaderBackend::Scanned(_) => Ok(()),
        }
    }
}

#[derive(Clone)]
struct FastaReferenceMetadata {
    name: Vec<u8>,
    length: u64,
}

fn open_indexed_reader(path: &Path) -> Result<IndexedFastaReader, CallError> {
    IndexedFastaReader::open(path).map_err(|error| {
        CallError::with_source(
            CallErrorKind::Input,
            format!("open indexed reference FASTA {}", path.display()),
            error,
        )
    })
}

fn detect_reference_compression(path: &Path) -> Result<Compression, CallError> {
    let reader = DecodedReader::open(path).map_err(|error| {
        CallError::with_source(
            CallErrorKind::Input,
            format!("inspect reference FASTA {} compression", path.display()),
            error,
        )
    })?;
    let compression = reader.compression();
    reader.close().map_err(|error| {
        CallError::with_source(
            CallErrorKind::Input,
            format!(
                "close reference FASTA {} after compression detection",
                path.display()
            ),
            error,
        )
    })?;
    Ok(compression)
}

fn prepare_plain_reference(
    path: &Path,
) -> Result<(CallReferenceBackend, Vec<FastaReferenceMetadata>), CallError> {
    let fai_path = adjacent_fai_path(path);
    match fs::metadata(&fai_path) {
        Ok(_) => indexed_reference_metadata(path),
        Err(error) if error.kind() == io::ErrorKind::NotFound => {
            let index = Arc::new(ScannedFastaIndex::build(path)?);
            let fasta_references = index
                .references
                .iter()
                .map(|reference| FastaReferenceMetadata {
                    name: reference.name.clone(),
                    length: reference.length,
                })
                .collect::<Vec<_>>();
            Ok((CallReferenceBackend::Scanned(index), fasta_references))
        }
        Err(error) => Err(CallError::with_source(
            CallErrorKind::Input,
            format!("inspect FASTA index {}", fai_path.display()),
            error,
        )),
    }
}

fn indexed_reference_metadata(
    path: &Path,
) -> Result<(CallReferenceBackend, Vec<FastaReferenceMetadata>), CallError> {
    let reader = open_indexed_reader(path)?;
    let fasta_references = reader
        .references()
        .iter()
        .map(|reference| FastaReferenceMetadata {
            name: reference.name().to_vec(),
            length: reference.length(),
        })
        .collect::<Vec<_>>();
    reader.close().map_err(|error| {
        CallError::with_source(
            CallErrorKind::Input,
            format!("close indexed reference FASTA {}", path.display()),
            error,
        )
    })?;
    Ok((CallReferenceBackend::Indexed, fasta_references))
}

fn require_bgzf_fasta_indexes(path: &Path) -> Result<(), CallError> {
    for sidecar in [adjacent_fai_path(path), adjacent_gzi_path(path)] {
        match fs::metadata(&sidecar) {
            Ok(_) => {}
            Err(error) if error.kind() == io::ErrorKind::NotFound => {
                return Err(CallError::input(format!(
                    "BGZF-compressed reference FASTA {} requires adjacent .fai and .gzi indexes; missing {}. Create both with `samtools faidx REFERENCE.bgzf.fa.gz`",
                    path.display(),
                    sidecar.display()
                )));
            }
            Err(error) => {
                return Err(CallError::with_source(
                    CallErrorKind::Input,
                    format!("inspect BGZF FASTA index {}", sidecar.display()),
                    error,
                ));
            }
        }
    }
    Ok(())
}

fn unsupported_gzip_reference(path: &Path) -> CallError {
    CallError::input(format!(
        "reference FASTA {} uses ordinary gzip compression, which is unsupported because it cannot provide random access; use plain FASTA or BGZF-compressed FASTA. To convert, run `gzip -cd INPUT.fa.gz | bgzip -c > REFERENCE.bgzf.fa.gz`, then `samtools faidx REFERENCE.bgzf.fa.gz`",
        path.display()
    ))
}

fn adjacent_fai_path(path: &Path) -> PathBuf {
    let mut value = OsString::from(path.as_os_str());
    value.push(".fai");
    PathBuf::from(value)
}

fn adjacent_gzi_path(path: &Path) -> PathBuf {
    let mut value = OsString::from(path.as_os_str());
    value.push(".gzi");
    PathBuf::from(value)
}

fn map_reference_dictionary(
    path: &Path,
    references: &[BamReference],
    fasta_references: &[FastaReferenceMetadata],
) -> Result<Vec<u32>, CallError> {
    let mut fasta_reference_ids = Vec::new();
    fasta_reference_ids
        .try_reserve_exact(references.len())
        .map_err(|error| {
            CallError::with_source(
                CallErrorKind::Input,
                "reserve reference dictionary mapping",
                error,
            )
        })?;
    for (bam_id, bam_reference) in references.iter().enumerate() {
        let matches = fasta_references
            .iter()
            .enumerate()
            .filter(|(_, fasta_reference)| fasta_reference.name == bam_reference.name)
            .collect::<Vec<_>>();
        let [(fasta_id, fasta_reference)] = matches.as_slice() else {
            return Err(CallError::input(if matches.is_empty() {
                format!(
                    "reference FASTA {} is missing BAM contig {} (`{}`)",
                    path.display(),
                    bam_id,
                    String::from_utf8_lossy(&bam_reference.name)
                )
            } else {
                format!(
                    "reference FASTA {} contains duplicate contig `{}`",
                    path.display(),
                    String::from_utf8_lossy(&bam_reference.name)
                )
            }));
        };
        if fasta_reference.length != u64::from(bam_reference.length) {
            return Err(CallError::input(format!(
                "reference FASTA {} contig `{}` length {} differs from BAM dictionary length {}",
                path.display(),
                String::from_utf8_lossy(&bam_reference.name),
                fasta_reference.length,
                bam_reference.length
            )));
        }
        fasta_reference_ids.push(
            u32::try_from(*fasta_id)
                .map_err(|_| CallError::input("reference FASTA dictionary ordinal exceeds u32"))?,
        );
    }
    Ok(fasta_reference_ids)
}

struct ScannedFastaIndex {
    references: Vec<ScannedFastaReference>,
}

struct ScannedFastaReference {
    name: Vec<u8>,
    length: u64,
    runs: Vec<ScannedFastaRun>,
}

#[derive(Clone, Copy)]
struct ScannedFastaRun {
    start_base: u64,
    end_base: u64,
    file_offset: u64,
    line_bases: u64,
    line_bytes: u64,
}

struct ScannedReferenceBuilder {
    name: Vec<u8>,
    length: u64,
    runs: Vec<ScannedFastaRun>,
}

impl ScannedReferenceBuilder {
    fn new(name: Vec<u8>) -> Self {
        Self {
            name,
            length: 0,
            runs: Vec::new(),
        }
    }

    fn push_line(
        &mut self,
        path: &Path,
        line_number: u64,
        file_offset: u64,
        line_bases: u64,
        line_bytes: u64,
    ) -> Result<(), CallError> {
        if line_bases == 0 {
            return Err(CallError::input(format!(
                "reference FASTA {} line {line_number} contains an empty sequence line",
                path.display()
            )));
        }
        let end_base = self.length.checked_add(line_bases).ok_or_else(|| {
            CallError::input(format!(
                "reference FASTA {} contig `{}` length exceeds u64",
                path.display(),
                String::from_utf8_lossy(&self.name)
            ))
        })?;
        if let Some(run) = self.runs.last_mut()
            && run.line_bases == line_bases
            && run.line_bytes == line_bytes
        {
            run.end_base = end_base;
        } else {
            self.runs.try_reserve(1).map_err(|error| {
                CallError::with_source(
                    CallErrorKind::Input,
                    format!(
                        "reserve FASTA line-layout run for contig `{}`",
                        String::from_utf8_lossy(&self.name)
                    ),
                    error,
                )
            })?;
            self.runs.push(ScannedFastaRun {
                start_base: self.length,
                end_base,
                file_offset,
                line_bases,
                line_bytes,
            });
        }
        self.length = end_base;
        Ok(())
    }

    fn finish(self, path: &Path) -> Result<ScannedFastaReference, CallError> {
        if self.length == 0 {
            return Err(CallError::input(format!(
                "reference FASTA {} contig `{}` has no sequence",
                path.display(),
                String::from_utf8_lossy(&self.name)
            )));
        }
        Ok(ScannedFastaReference {
            name: self.name,
            length: self.length,
            runs: self.runs,
        })
    }
}

impl ScannedFastaIndex {
    fn build(path: &Path) -> Result<Self, CallError> {
        const MAX_FASTA_LINE_BASES: usize = 1_000_000;
        const MAX_FASTA_PHYSICAL_LINE_BYTES: u64 = 1_000_003;

        let file = File::open(path).map_err(|error| {
            CallError::with_source(
                CallErrorKind::Input,
                format!("open reference FASTA {} for direct scan", path.display()),
                error,
            )
        })?;
        let mut reader = BufReader::new(file);
        let mut references = Vec::new();
        let mut current: Option<ScannedReferenceBuilder> = None;
        let mut line = Vec::new();
        let mut line_number = 0_u64;
        let mut file_offset = 0_u64;
        loop {
            line.clear();
            let read = Read::by_ref(&mut reader)
                .take(MAX_FASTA_PHYSICAL_LINE_BYTES)
                .read_until(b'\n', &mut line)
                .map_err(|error| {
                    CallError::with_source(
                        CallErrorKind::Input,
                        format!("scan reference FASTA {}", path.display()),
                        error,
                    )
                })?;
            if read == 0 {
                break;
            }
            line_number = line_number.checked_add(1).ok_or_else(|| {
                CallError::input(format!(
                    "reference FASTA {} physical line count exceeds u64",
                    path.display()
                ))
            })?;
            let line_bytes = u64::try_from(read).map_err(|_| {
                CallError::input(format!(
                    "reference FASTA {} line {line_number} length exceeds u64",
                    path.display()
                ))
            })?;
            let content = fasta_line_content(&line);
            if content.len() > MAX_FASTA_LINE_BASES {
                return Err(CallError::input(format!(
                    "reference FASTA {} line {line_number} exceeds {MAX_FASTA_LINE_BASES} bytes",
                    path.display()
                )));
            }
            if content.first() == Some(&b'>') {
                if let Some(builder) = current.take() {
                    push_scanned_reference(&mut references, builder.finish(path)?)?;
                }
                current = Some(ScannedReferenceBuilder::new(parse_fasta_name(
                    path,
                    line_number,
                    content,
                )?));
            } else {
                let builder = current.as_mut().ok_or_else(|| {
                    CallError::input(format!(
                        "reference FASTA {} line {line_number} appears before the first header",
                        path.display()
                    ))
                })?;
                let line_bases = u64::try_from(content.len()).map_err(|_| {
                    CallError::input(format!(
                        "reference FASTA {} line {line_number} base count exceeds u64",
                        path.display()
                    ))
                })?;
                builder.push_line(path, line_number, file_offset, line_bases, line_bytes)?;
            }
            file_offset = file_offset.checked_add(line_bytes).ok_or_else(|| {
                CallError::input(format!(
                    "reference FASTA {} byte offset exceeds u64",
                    path.display()
                ))
            })?;
        }
        if let Some(builder) = current {
            push_scanned_reference(&mut references, builder.finish(path)?)?;
        }
        if references.is_empty() {
            return Err(CallError::input(format!(
                "reference FASTA {} contains no records",
                path.display()
            )));
        }
        Ok(Self { references })
    }
}

fn push_scanned_reference(
    references: &mut Vec<ScannedFastaReference>,
    reference: ScannedFastaReference,
) -> Result<(), CallError> {
    references.try_reserve(1).map_err(|error| {
        CallError::with_source(
            CallErrorKind::Input,
            "reserve scanned FASTA reference dictionary",
            error,
        )
    })?;
    references.push(reference);
    Ok(())
}

fn fasta_line_content(line: &[u8]) -> &[u8] {
    let Some(without_newline) = line.strip_suffix(b"\n") else {
        return line;
    };
    without_newline
        .strip_suffix(b"\r")
        .unwrap_or(without_newline)
}

fn parse_fasta_name(path: &Path, line_number: u64, line: &[u8]) -> Result<Vec<u8>, CallError> {
    let tail = line.strip_prefix(b">").ok_or_else(|| {
        CallError::input(format!(
            "reference FASTA {} line {line_number} has an invalid header marker",
            path.display()
        ))
    })?;
    let separator = tail
        .iter()
        .position(|byte| matches!(*byte, b' ' | b'\t'))
        .unwrap_or(tail.len());
    let name = &tail[..separator];
    if name.is_empty()
        || name
            .iter()
            .any(|byte| !byte.is_ascii_graphic() || byte.is_ascii_whitespace())
    {
        return Err(CallError::input(format!(
            "reference FASTA {} line {line_number} has an invalid or empty contig name",
            path.display()
        )));
    }
    if tail[separator..]
        .iter()
        .any(|byte| !matches!(*byte, b' ' | b'\t') && !byte.is_ascii_graphic())
    {
        return Err(CallError::input(format!(
            "reference FASTA {} line {line_number} has invalid header text",
            path.display()
        )));
    }
    Ok(name.to_vec())
}

struct ScannedFastaReader {
    path: PathBuf,
    reader: BufReader<File>,
    index: Arc<ScannedFastaIndex>,
    scratch: Vec<u8>,
}

impl ScannedFastaReader {
    fn open(path: &Path, index: Arc<ScannedFastaIndex>) -> Result<Self, CallError> {
        let file = File::open(path).map_err(|error| {
            CallError::with_source(
                CallErrorKind::Input,
                format!("open scanned reference FASTA {}", path.display()),
                error,
            )
        })?;
        Ok(Self {
            path: path.to_path_buf(),
            reader: BufReader::new(file),
            index,
            scratch: Vec::new(),
        })
    }

    fn fetch(&mut self, reference_id: u32, start: u64, end: u64) -> Result<Vec<u8>, CallError> {
        let reference_ordinal = usize::try_from(reference_id)
            .map_err(|_| CallError::input("reference FASTA ordinal exceeds usize"))?;
        let index = Arc::clone(&self.index);
        let reference = index
            .references
            .get(reference_ordinal)
            .ok_or_else(|| CallError::input("reference FASTA ordinal is out of range"))?;
        if start >= end || end > reference.length {
            return Err(CallError::input(format!(
                "reference FASTA {} interval {start}-{end} is outside contig `{}` length {}",
                self.path.display(),
                String::from_utf8_lossy(&reference.name),
                reference.length
            )));
        }
        let requested = usize::try_from(end - start).map_err(|_| {
            CallError::input(format!(
                "reference FASTA {} requested interval exceeds usize",
                self.path.display()
            ))
        })?;
        let mut bases = Vec::new();
        bases.try_reserve_exact(requested).map_err(|error| {
            CallError::with_source(
                CallErrorKind::Input,
                format!("reserve {requested} fetched reference bases"),
                error,
            )
        })?;
        let mut position = start;
        while position < end {
            let run_index = reference
                .runs
                .partition_point(|run| run.end_base <= position);
            let run = reference.runs.get(run_index).ok_or_else(|| {
                CallError::input(format!(
                    "reference FASTA {} scanned layout does not cover contig `{}` position {position}",
                    self.path.display(),
                    String::from_utf8_lossy(&reference.name)
                ))
            })?;
            if position < run.start_base {
                return Err(CallError::input(format!(
                    "reference FASTA {} scanned layout has a gap before contig `{}` position {position}",
                    self.path.display(),
                    String::from_utf8_lossy(&reference.name)
                )));
            }
            let segment_end = end.min(run.end_base);
            self.fetch_run_segment(run, position, segment_end, &mut bases)?;
            position = segment_end;
        }
        if bases.len() != requested {
            return Err(CallError::input(format!(
                "reference FASTA {} changed after its in-memory position table was built",
                self.path.display()
            )));
        }
        bases.make_ascii_uppercase();
        Ok(bases)
    }

    fn fetch_run_segment(
        &mut self,
        run: &ScannedFastaRun,
        start: u64,
        end: u64,
        bases: &mut Vec<u8>,
    ) -> Result<(), CallError> {
        let relative_start = start
            .checked_sub(run.start_base)
            .ok_or_else(|| CallError::input("reference FASTA run-relative start underflow"))?;
        let relative_last = end
            .checked_sub(run.start_base)
            .and_then(|value| value.checked_sub(1))
            .ok_or_else(|| CallError::input("reference FASTA run-relative end underflow"))?;
        let start_byte = run
            .file_offset
            .checked_add(
                (relative_start / run.line_bases)
                    .checked_mul(run.line_bytes)
                    .ok_or_else(|| CallError::input("reference FASTA byte offset overflow"))?,
            )
            .and_then(|value| value.checked_add(relative_start % run.line_bases))
            .ok_or_else(|| CallError::input("reference FASTA byte offset overflow"))?;
        let last_byte = run
            .file_offset
            .checked_add(
                (relative_last / run.line_bases)
                    .checked_mul(run.line_bytes)
                    .ok_or_else(|| CallError::input("reference FASTA byte offset overflow"))?,
            )
            .and_then(|value| value.checked_add(relative_last % run.line_bases))
            .ok_or_else(|| CallError::input("reference FASTA byte offset overflow"))?;
        let physical_length = last_byte
            .checked_sub(start_byte)
            .and_then(|value| value.checked_add(1))
            .ok_or_else(|| CallError::input("reference FASTA physical span overflow"))?;
        let physical_length = usize::try_from(physical_length)
            .map_err(|_| CallError::input("reference FASTA physical span exceeds usize"))?;
        self.scratch.clear();
        self.scratch
            .try_reserve_exact(physical_length)
            .map_err(|error| {
                CallError::with_source(
                    CallErrorKind::Input,
                    format!("reserve {physical_length} reference input bytes"),
                    error,
                )
            })?;
        self.scratch.resize(physical_length, 0);
        self.reader
            .seek(SeekFrom::Start(start_byte))
            .map_err(|error| {
                CallError::with_source(
                    CallErrorKind::Input,
                    format!(
                        "seek reference FASTA {} to byte {start_byte}",
                        self.path.display()
                    ),
                    error,
                )
            })?;
        self.reader.read_exact(&mut self.scratch).map_err(|error| {
            CallError::with_source(
                CallErrorKind::Input,
                format!(
                    "read reference FASTA {} byte interval {start_byte}-{}",
                    self.path.display(),
                    last_byte.saturating_add(1)
                ),
                error,
            )
        })?;
        bases.extend(
            self.scratch
                .iter()
                .copied()
                .filter(|byte| !matches!(*byte, b'\n' | b'\r')),
        );
        Ok(())
    }
}

pub(crate) struct ReferenceWindow {
    reference: u32,
    start: u32,
    bases: Vec<u8>,
}

impl ReferenceWindow {
    pub(crate) fn base(&self, reference: u32, position: u32) -> Option<u8> {
        if reference != self.reference {
            return None;
        }
        let offset = usize::try_from(position.checked_sub(self.start)?).ok()?;
        self.bases.get(offset).copied()
    }

    pub(crate) fn context(
        &self,
        reference: u32,
        position: u32,
        strand: EvidenceStrand,
    ) -> Option<CytosineContext> {
        let (first, second) = match strand {
            EvidenceStrand::Top => {
                let first = canonical(self.base(reference, position.checked_add(1)?)?)?;
                if first == b'G' {
                    return Some(CytosineContext {
                        class: ContextClass::Cg,
                        second: first,
                    });
                }
                let second = canonical(self.base(reference, position.checked_add(2)?)?)?;
                (first, second)
            }
            EvidenceStrand::Bottom => {
                let first = complement(self.base(reference, position.checked_sub(1)?)?)?;
                if first == b'G' {
                    return Some(CytosineContext {
                        class: ContextClass::Cg,
                        second: first,
                    });
                }
                let second = complement(self.base(reference, position.checked_sub(2)?)?)?;
                (first, second)
            }
        };
        Some(CytosineContext {
            class: if second == b'G' {
                ContextClass::Chg
            } else {
                ContextClass::Chh
            },
            second: first,
        })
    }
}

const fn canonical(base: u8) -> Option<u8> {
    match base {
        b'A' | b'C' | b'G' | b'T' => Some(base),
        _ => None,
    }
}

const fn complement(base: u8) -> Option<u8> {
    match base {
        b'A' => Some(b'T'),
        b'C' => Some(b'G'),
        b'G' => Some(b'C'),
        b'T' => Some(b'A'),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use std::fs;
    use std::io::Write;
    use std::path::PathBuf;
    use std::time::{SystemTime, UNIX_EPOCH};

    use bsbit_core::reference::ReferenceSemanticDigestBuilder;
    use bsbit_hts::BgzfWriter;

    use super::{CallReferenceSource, ReferenceWindow, adjacent_fai_path, adjacent_gzi_path};
    use crate::call_input::BamReference;
    use crate::evidence::{ContextClass, CytosineContext, EvidenceStrand};

    fn unique_directory(label: &str) -> PathBuf {
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock is after epoch")
            .as_nanos();
        std::env::temp_dir().join(format!(
            "bsbit-call-reference-{label}-{}-{nonce}",
            std::process::id()
        ))
    }

    #[test]
    fn missing_fai_builds_a_shared_in_memory_layout_without_a_sidecar() {
        let directory = unique_directory("scanned-layout");
        fs::create_dir(&directory).unwrap();
        let fasta = directory.join("reference.fa");
        fs::write(&fasta, b">chr1 description\r\nACGT\r\ntgca\r\nCT\nGATC").unwrap();
        let references = [BamReference {
            name: b"chr1".to_vec(),
            length: 14,
        }];

        let source = CallReferenceSource::prepare(&fasta, &references).unwrap();
        assert!(!adjacent_fai_path(&fasta).exists());
        let mut reader = source.open().unwrap();
        assert_eq!(reader.fetch(0, 0, 14).unwrap(), b"ACGTTGCACTGATC");
        assert_eq!(reader.fetch(0, 2, 13).unwrap(), b"GTTGCACTGAT");
        let mut digest = ReferenceSemanticDigestBuilder::new(1);
        digest
            .push_ascii_contig(b"chr1", b"ACGTTGCACTGATC")
            .unwrap();
        reader
            .validate_semantic_digest(&references, digest.finish().unwrap())
            .unwrap();
        reader.close().unwrap();
        assert!(!adjacent_fai_path(&fasta).exists());
        fs::remove_dir_all(directory).unwrap();
    }

    #[test]
    fn ordinary_gzip_fasta_is_rejected_with_conversion_guidance() {
        const GZIP_FASTA: &[u8] = &[
            0x1f, 0x8b, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x02, 0x03, 0xb3, 0x4b, 0xce, 0x28,
            0x32, 0xe4, 0x72, 0xe4, 0x02, 0x00, 0x8e, 0x7b, 0x39, 0x64, 0x08, 0x00, 0x00, 0x00,
        ];
        let directory = unique_directory("ordinary-gzip");
        fs::create_dir(&directory).unwrap();
        let fasta = directory.join("reference.fa.gz");
        fs::write(&fasta, GZIP_FASTA).unwrap();
        let references = [BamReference {
            name: b"chr1".to_vec(),
            length: 1,
        }];

        let Err(error) = CallReferenceSource::prepare(&fasta, &references) else {
            panic!("ordinary gzip FASTA must fail");
        };
        assert!(error.to_string().contains("ordinary gzip compression"));
        assert!(error.to_string().contains("bgzip -c"));
        assert!(error.to_string().contains("samtools faidx"));
        assert!(!adjacent_fai_path(&fasta).exists());
        fs::remove_dir_all(directory).unwrap();
    }

    #[test]
    fn bgzf_fasta_without_indexes_reports_both_required_sidecars() {
        let directory = unique_directory("bgzf-no-index");
        fs::create_dir(&directory).unwrap();
        let fasta = directory.join("reference.data");
        let file = fs::File::create(&fasta).unwrap();
        let mut writer = BgzfWriter::from_file(file, 0).unwrap();
        writer.write_all(b">chr1\nA\n").unwrap();
        writer.finish().unwrap().sync_all().unwrap();
        let references = [BamReference {
            name: b"chr1".to_vec(),
            length: 1,
        }];

        let Err(error) = CallReferenceSource::prepare(&fasta, &references) else {
            panic!("BGZF FASTA without indexes must fail");
        };
        assert!(
            error
                .to_string()
                .contains("BGZF-compressed reference FASTA")
        );
        assert!(error.to_string().contains(".fai and .gzi"));
        assert!(error.to_string().contains("samtools faidx"));
        assert!(!adjacent_fai_path(&fasta).exists());
        assert!(!adjacent_gzi_path(&fasta).exists());
        fs::remove_dir_all(directory).unwrap();
    }

    #[test]
    fn reference_window_resolves_context_on_both_strands_and_edges() {
        let window = ReferenceWindow {
            reference: 0,
            start: 8,
            bases: b"ACGCTGG".to_vec(),
        };
        assert_eq!(
            window.context(0, 9, EvidenceStrand::Top),
            Some(CytosineContext {
                class: ContextClass::Cg,
                second: b'G',
            })
        );
        assert_eq!(
            window.context(0, 13, EvidenceStrand::Bottom),
            Some(CytosineContext {
                class: ContextClass::Chg,
                second: b'A',
            })
        );
        assert_eq!(window.context(0, 8, EvidenceStrand::Bottom), None);
        assert_eq!(window.context(1, 9, EvidenceStrand::Top), None);
    }
}
