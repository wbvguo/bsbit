//! Owning record materialization used as a white-box qualification oracle.

use super::auxiliary::{oriented_base, replay_pass};
use super::direct::authoritative_literal_nm;
use super::*;

/// Constructs one deterministic single-read record.
///
/// # Errors
///
/// Returns [`RecordBuildError`] for invalid read metadata, limits,
/// owner/coordinate/CIGAR failures, replay disagreement, or allocation failure.
#[cfg(test)]
pub(crate) fn build_single_alignment_record(
    reference: &ReferenceIndex,
    query_name: &[u8],
    read: AlignmentRead<'_>,
    alignment: Option<&VerifiedAlignment>,
    mapping_quality: RecordMappingQuality,
    limits: AlignmentRecordLimits,
) -> Result<AlignmentRecord, RecordBuildError> {
    build_single_alignment_record_with_auxiliary_mode(
        reference,
        query_name,
        read,
        alignment,
        mapping_quality,
        limits,
        AlignmentAuxiliaryMode::Minimal,
    )
}

/// Constructs one deterministic single-read record with an explicit optional
/// field materialization policy.
///
/// # Errors
///
/// Returns the same failures as [`build_single_alignment_record`].
#[doc(hidden)]
#[cfg(test)]
pub(crate) fn build_single_alignment_record_with_auxiliary_mode(
    reference: &ReferenceIndex,
    query_name: &[u8],
    read: AlignmentRead<'_>,
    alignment: Option<&VerifiedAlignment>,
    mapping_quality: RecordMappingQuality,
    limits: AlignmentRecordLimits,
    auxiliary_mode: AlignmentAuxiliaryMode,
) -> Result<AlignmentRecord, RecordBuildError> {
    prepare_record(
        reference,
        query_name,
        read,
        alignment,
        RecordSegment::Unpaired,
        false,
        mapping_quality,
        None,
        0,
        limits,
        auxiliary_mode,
    )
}

#[cfg(test)]
struct PreparedRecord {
    mapping: Option<MappedAlignmentRecord>,
    sequence: Box<[u8]>,
    quality: Option<Box<[u8]>>,
}

#[allow(clippy::too_many_arguments)]
#[cfg(test)]
fn prepare_record(
    reference: &ReferenceIndex,
    query_name: &[u8],
    read: AlignmentRead<'_>,
    alignment: Option<&VerifiedAlignment>,
    segment: RecordSegment,
    proper_pair: bool,
    mapping_quality: RecordMappingQuality,
    mate: Option<RecordMateLocation>,
    template_length: i32,
    limits: AlignmentRecordLimits,
    auxiliary_mode: AlignmentAuxiliaryMode,
) -> Result<AlignmentRecord, RecordBuildError> {
    let prepared = prepare_mapping(reference, read, alignment, limits, auxiliary_mode)?;
    finish_prepared_record(
        query_name,
        segment,
        proper_pair,
        mapping_quality,
        prepared,
        mate,
        template_length,
        limits,
    )
}

#[allow(clippy::too_many_arguments)]
#[cfg(test)]
fn finish_prepared_record(
    query_name: &[u8],
    segment: RecordSegment,
    proper_pair: bool,
    mapping_quality: RecordMappingQuality,
    prepared: PreparedRecord,
    mate: Option<RecordMateLocation>,
    template_length: i32,
    limits: AlignmentRecordLimits,
) -> Result<AlignmentRecord, RecordBuildError> {
    AlignmentRecord::new(
        query_name,
        segment,
        proper_pair,
        mapping_quality,
        prepared.mapping,
        mate,
        template_length,
        &prepared.sequence,
        prepared.quality.as_deref(),
        limits,
    )
    .map_err(Into::into)
}

#[cfg(test)]
fn prepare_mapping(
    reference: &ReferenceIndex,
    read: AlignmentRead<'_>,
    alignment: Option<&VerifiedAlignment>,
    limits: AlignmentRecordLimits,
    auxiliary_mode: AlignmentAuxiliaryMode,
) -> Result<PreparedRecord, RecordBuildError> {
    let orientation = alignment.map_or(
        AlignmentOrientation::Forward,
        VerifiedAlignment::orientation,
    );
    let sequence = oriented_sequence(read.sequence(), orientation)?;
    let quality = oriented_quality(read.quality(), orientation)?;
    let mapping = alignment
        .map(|alignment| {
            build_mapped_record(
                reference,
                read.sequence(),
                alignment,
                limits,
                auxiliary_mode,
            )
        })
        .transpose()?;
    Ok(PreparedRecord {
        mapping,
        sequence,
        quality,
    })
}

#[cfg(test)]
fn build_mapped_record(
    reference: &ReferenceIndex,
    read: &NormalizedSequence,
    alignment: &VerifiedAlignment,
    limits: AlignmentRecordLimits,
    auxiliary_mode: AlignmentAuxiliaryMode,
) -> Result<MappedAlignmentRecord, RecordBuildError> {
    let contig = reference
        .resolve_contig(alignment.contig())
        .map_err(|source| RecordBuildError::ReferenceAccess { source })?;
    let interval = alignment.interval();

    let authoritative_literal_nm = authoritative_literal_nm(reference, read, alignment)?;
    let replay = match auxiliary_mode {
        AlignmentAuxiliaryMode::Minimal => LiteralReplay {
            literal_nm: authoritative_literal_nm,
            md: None,
            bismark_xm: None,
        },
        AlignmentAuxiliaryMode::Bismark => literal_replay(
            contig.sequence(),
            interval,
            read,
            alignment,
            limits,
            auxiliary_mode,
        )?,
    };
    if replay.literal_nm != authoritative_literal_nm {
        return Err(RecordBuildError::AlignmentLiteralNmMismatch {
            expected: authoritative_literal_nm,
            observed: replay.literal_nm,
        });
    }
    let literal_nm =
        u32::try_from(replay.literal_nm).map_err(|_| RecordBuildError::FieldOutOfRange {
            field: AlignmentRecordField::Nm,
            value: replay.literal_nm,
        })?;
    let record_reference = RecordReference::new(
        contig.ordinal(),
        contig.name(),
        contig.sequence().len(),
        interval,
    )?;
    MappedAlignmentRecord::new(
        record_reference,
        alignment.orientation(),
        alignment.strand(),
        alignment.cytosine_strand(),
        alignment.cigar().clone(),
        read.len(),
        literal_nm,
        replay.md.as_deref(),
        replay.bismark_xm.as_deref(),
        limits,
    )
    .map_err(Into::into)
}

#[cfg(test)]
struct LiteralReplay {
    literal_nm: u64,
    md: Option<Box<[u8]>>,
    bismark_xm: Option<Box<[u8]>>,
}

#[cfg(test)]
const INITIAL_MD_CAPACITY: u64 = 64;

#[cfg(test)]
fn literal_replay(
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

#[cfg(test)]
fn oriented_sequence(
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
fn oriented_quality(
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
