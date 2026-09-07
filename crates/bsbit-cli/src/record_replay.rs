//! Literal CIGAR replay and calling-oriented alignment annotations.
//!
//! This module derives authoritative NM/MD values, Bismark methylation calls,
//! and orientation-correct read buffers without owning record assembly.

use bsbit_align::extension::VerifiedAlignment;
use bsbit_core::alphabet::Base;
use bsbit_core::bisulfite::{AlignmentOrientation, CytosineStrand};
use bsbit_core::cigar::CoreCigarOp;
use bsbit_core::coordinate::ReferenceInterval;
use bsbit_core::sequence::NormalizedSequence;
#[cfg(test)]
use bsbit_hts::{AlignmentAuxiliaryMode, AlignmentRecordLimits};
use bsbit_hts::{AlignmentRecordAllocation, AlignmentRecordResource};

#[cfg(test)]
use crate::record_support::check_limit;
use crate::record_support::{
    RecordBuildError, append_u64, checked_add_resource, decimal_digits, storage_count, storage_len,
};

#[cfg(test)]
pub(super) struct LiteralReplay {
    pub(super) literal_nm: u64,
    pub(super) md: Option<Box<[u8]>>,
    pub(super) bismark_xm: Option<Box<[u8]>>,
}

#[cfg(test)]
const INITIAL_MD_CAPACITY: u64 = 64;

#[cfg(test)]
pub(super) fn literal_replay(
    reference: &NormalizedSequence,
    interval: ReferenceInterval,
    read: &NormalizedSequence,
    alignment: &VerifiedAlignment,
    limits: AlignmentRecordLimits,
    auxiliary_mode: AlignmentAuxiliaryMode,
) -> Result<LiteralReplay, RecordBuildError> {
    let emit_bismark = matches!(auxiliary_mode, AlignmentAuxiliaryMode::Bismark);
    let mut md = Vec::new();
    if emit_bismark {
        let initial_capacity = limits.max_md_bytes().min(INITIAL_MD_CAPACITY);
        let initial_capacity = storage_count(initial_capacity, AlignmentRecordAllocation::Md)?;
        md.try_reserve_exact(initial_capacity)
            .map_err(|_| RecordBuildError::AllocationFailed {
                allocation: AlignmentRecordAllocation::Md,
                requested: storage_len(initial_capacity),
            })?;
    }
    let mut bismark_xm = Vec::new();
    let summary = replay_pass(
        reference,
        interval,
        read,
        alignment,
        emit_bismark.then_some(&mut md),
        emit_bismark.then_some(&mut bismark_xm),
        limits.max_md_bytes(),
    )?;
    check_limit(
        summary.md_bytes,
        limits.max_md_bytes(),
        AlignmentRecordResource::MdBytes,
    )?;
    debug_assert_eq!(storage_len(md.len()), summary.md_bytes);
    debug_assert_eq!(storage_len(bismark_xm.len()), summary.bismark_xm_bytes);
    Ok(LiteralReplay {
        literal_nm: summary.literal_nm,
        md: emit_bismark.then(|| md.into_boxed_slice()),
        bismark_xm: emit_bismark.then(|| bismark_xm.into_boxed_slice()),
    })
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) struct ReplaySummary {
    pub(super) literal_nm: u64,
    pub(super) md_bytes: u64,
    pub(super) bismark_xm_bytes: u64,
}

// This is one bounds-checked pass over the CIGAR. Keeping its cursor and output
// state together avoids subtly divergent NM, MD, and XM traversal semantics.
#[allow(clippy::too_many_lines)]
pub(super) fn replay_pass(
    reference: &NormalizedSequence,
    interval: ReferenceInterval,
    read: &NormalizedSequence,
    alignment: &VerifiedAlignment,
    mut md: Option<&mut Vec<u8>>,
    mut bismark_xm: Option<&mut Vec<u8>>,
    max_md_bytes: u64,
) -> Result<ReplaySummary, RecordBuildError> {
    let md_start = md.as_ref().map_or(0, |output| storage_len(output.len()));
    let bismark_xm_start = bismark_xm.as_ref().map_or(0, |output| output.len());
    if let Some(output) = bismark_xm.as_deref_mut() {
        let requested = read.len();
        output
            .try_reserve(read.bases().len())
            .map_err(|_| RecordBuildError::AllocationFailed {
                allocation: AlignmentRecordAllocation::MethylationCall,
                requested,
            })?;
        output.resize(output.len().saturating_add(read.bases().len()), b'.');
    }
    let mut reference_index = storage_count(interval.start(), AlignmentRecordAllocation::Md)?;
    let mut query_index = 0_usize;
    let mut matches = 0_u64;
    let mut literal_nm = 0_u64;
    let mut md_bytes = 0_u64;

    for run in alignment.cigar().runs() {
        let length = storage_count(run.length(), AlignmentRecordAllocation::Md)?;
        match run.operation() {
            CoreCigarOp::M => {
                for _ in 0..length {
                    let reference_base = reference.bases()[reference_index];
                    let query_base = oriented_base(read, alignment.orientation(), query_index);
                    if let Some(output) = bismark_xm.as_deref_mut() {
                        output[bismark_xm_start + query_index] = bismark_methylation_call(
                            reference.bases(),
                            reference_index,
                            reference_base,
                            query_base,
                            alignment.cytosine_strand(),
                        );
                    }
                    let literal = is_literal_acgt_match(reference_base, query_base);
                    if literal {
                        matches =
                            checked_add_resource(matches, 1, AlignmentRecordResource::MdBytes)?;
                    } else {
                        if let Some(output) = md.as_deref_mut() {
                            let next_md_bytes = checked_add_resource(
                                md_bytes,
                                decimal_digits(matches) + 1,
                                AlignmentRecordResource::MdBytes,
                            )?;
                            if next_md_bytes <= max_md_bytes {
                                reserve_md_total(
                                    output,
                                    checked_add_resource(
                                        md_start,
                                        next_md_bytes,
                                        AlignmentRecordResource::MdBytes,
                                    )?,
                                )?;
                                append_u64(output, matches);
                                output.push(reference_base.as_ascii());
                            }
                            md_bytes = next_md_bytes;
                        }
                        matches = 0;
                        literal_nm =
                            checked_add_resource(literal_nm, 1, AlignmentRecordResource::MdBytes)?;
                    }
                    reference_index += 1;
                    query_index += 1;
                }
            }
            CoreCigarOp::I => {
                let run_length = run.length();
                literal_nm =
                    checked_add_resource(literal_nm, run_length, AlignmentRecordResource::MdBytes)?;
                query_index += length;
            }
            CoreCigarOp::D => {
                let run_length = run.length();
                if let Some(output) = md.as_deref_mut() {
                    let next_md_bytes = checked_add_resource(
                        md_bytes,
                        decimal_digits(matches) + 1 + run_length,
                        AlignmentRecordResource::MdBytes,
                    )?;
                    if next_md_bytes <= max_md_bytes {
                        reserve_md_total(
                            output,
                            checked_add_resource(
                                md_start,
                                next_md_bytes,
                                AlignmentRecordResource::MdBytes,
                            )?,
                        )?;
                        append_u64(output, matches);
                        output.push(b'^');
                        for base in &reference.bases()[reference_index..reference_index + length] {
                            output.push(base.as_ascii());
                        }
                    }
                    md_bytes = next_md_bytes;
                }
                matches = 0;
                literal_nm =
                    checked_add_resource(literal_nm, run_length, AlignmentRecordResource::MdBytes)?;
                reference_index += length;
            }
        }
    }
    debug_assert_eq!(query_index, read.bases().len());
    if let Some(output) = md {
        let next_md_bytes = checked_add_resource(
            md_bytes,
            decimal_digits(matches),
            AlignmentRecordResource::MdBytes,
        )?;
        if next_md_bytes <= max_md_bytes {
            reserve_md_total(
                output,
                checked_add_resource(md_start, next_md_bytes, AlignmentRecordResource::MdBytes)?,
            )?;
            append_u64(output, matches);
        }
        md_bytes = next_md_bytes;
    }
    Ok(ReplaySummary {
        literal_nm,
        md_bytes,
        bismark_xm_bytes: bismark_xm.map_or(0, |_| storage_len(read.bases().len())),
    })
}

#[derive(Clone, Copy)]
enum BismarkMethylationContext {
    CpG,
    Chg,
    Chh,
    Unknown,
}

pub(super) fn bismark_methylation_call(
    reference: &[Base],
    reference_index: usize,
    reference_base: Base,
    query_base: Base,
    cytosine_strand: CytosineStrand,
) -> u8 {
    let (methylated, context) = match cytosine_strand {
        CytosineStrand::Top if reference_base == Base::C => match query_base {
            Base::C => (true, bismark_top_context(reference, reference_index)),
            Base::T => (false, bismark_top_context(reference, reference_index)),
            _ => return b'.',
        },
        CytosineStrand::Bottom if reference_base == Base::G => match query_base {
            Base::G => (true, bismark_bottom_context(reference, reference_index)),
            Base::A => (false, bismark_bottom_context(reference, reference_index)),
            _ => return b'.',
        },
        CytosineStrand::Top | CytosineStrand::Bottom => return b'.',
    };
    match (context, methylated) {
        (BismarkMethylationContext::CpG, false) => b'z',
        (BismarkMethylationContext::CpG, true) => b'Z',
        (BismarkMethylationContext::Chg, false) => b'x',
        (BismarkMethylationContext::Chg, true) => b'X',
        (BismarkMethylationContext::Chh, false) => b'h',
        (BismarkMethylationContext::Chh, true) => b'H',
        (BismarkMethylationContext::Unknown, false) => b'u',
        (BismarkMethylationContext::Unknown, true) => b'U',
    }
}

fn bismark_top_context(reference: &[Base], index: usize) -> BismarkMethylationContext {
    match reference.get(index.saturating_add(1)).copied() {
        Some(Base::G) => BismarkMethylationContext::CpG,
        Some(Base::A | Base::C | Base::T) => {
            match reference.get(index.saturating_add(2)).copied() {
                Some(Base::G) => BismarkMethylationContext::Chg,
                Some(Base::A | Base::C | Base::T) => BismarkMethylationContext::Chh,
                _ => BismarkMethylationContext::Unknown,
            }
        }
        _ => BismarkMethylationContext::Unknown,
    }
}

fn bismark_bottom_context(reference: &[Base], index: usize) -> BismarkMethylationContext {
    let Some(first_index) = index.checked_sub(1) else {
        return BismarkMethylationContext::Unknown;
    };
    match reference.get(first_index).copied() {
        Some(Base::C) => BismarkMethylationContext::CpG,
        Some(Base::A | Base::G | Base::T) => {
            let Some(second_index) = index.checked_sub(2) else {
                return BismarkMethylationContext::Unknown;
            };
            match reference.get(second_index).copied() {
                Some(Base::C) => BismarkMethylationContext::Chg,
                Some(Base::A | Base::G | Base::T) => BismarkMethylationContext::Chh,
                _ => BismarkMethylationContext::Unknown,
            }
        }
        _ => BismarkMethylationContext::Unknown,
    }
}

fn reserve_md_total(output: &mut Vec<u8>, requested: u64) -> Result<(), RecordBuildError> {
    let requested_storage = storage_count(requested, AlignmentRecordAllocation::Md)?;
    if output.capacity() < requested_storage {
        output
            .try_reserve_exact(requested_storage.saturating_sub(output.len()))
            .map_err(|_| RecordBuildError::AllocationFailed {
                allocation: AlignmentRecordAllocation::Md,
                requested,
            })?;
    }
    Ok(())
}

pub(super) fn oriented_base(
    read: &NormalizedSequence,
    orientation: AlignmentOrientation,
    index: usize,
) -> Base {
    match orientation {
        AlignmentOrientation::Forward => read.bases()[index],
        AlignmentOrientation::Reverse => read.bases()[read.bases().len() - 1 - index].complement(),
    }
}

fn is_literal_acgt_match(reference: Base, query: Base) -> bool {
    reference == query && !reference.is_unknown()
}

#[cfg(test)]
pub(super) fn oriented_sequence(
    sequence: &NormalizedSequence,
    orientation: AlignmentOrientation,
) -> Result<Box<[u8]>, RecordBuildError> {
    let requested = sequence.len();
    let capacity = storage_count(requested, AlignmentRecordAllocation::Sequence)?;
    let mut output = Vec::new();
    output
        .try_reserve_exact(capacity)
        .map_err(|_| RecordBuildError::AllocationFailed {
            allocation: AlignmentRecordAllocation::Sequence,
            requested,
        })?;
    for index in 0..capacity {
        output.push(oriented_base(sequence, orientation, index).as_ascii());
    }
    Ok(output.into_boxed_slice())
}

#[cfg(test)]
pub(super) fn oriented_quality(
    quality: Option<&[u8]>,
    orientation: AlignmentOrientation,
) -> Result<Option<Box<[u8]>>, RecordBuildError> {
    quality
        .map(|quality| {
            let requested = storage_len(quality.len());
            let mut output = Vec::new();
            output.try_reserve_exact(quality.len()).map_err(|_| {
                RecordBuildError::AllocationFailed {
                    allocation: AlignmentRecordAllocation::Quality,
                    requested,
                }
            })?;
            match orientation {
                AlignmentOrientation::Forward => output.extend_from_slice(quality),
                AlignmentOrientation::Reverse => output.extend(quality.iter().rev().copied()),
            }
            Ok(output.into_boxed_slice())
        })
        .transpose()
}
