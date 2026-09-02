//! Slab-backed direct SAM/BAM record composition.

use super::auxiliary::{ReplaySummary, oriented_base, replay_pass};
use super::{
    AlignmentAuxiliaryMode, AlignmentCigarOp, AlignmentCigarRun, AlignmentOrientation,
    AlignmentPlacement, AlignmentRead, AlignmentRecordAllocation, AlignmentRecordBatch,
    AlignmentRecordField, AlignmentRecordLimits, Base, BisulfiteStrand, BorrowedAlignmentRead,
    BorrowedAlignmentRecord, CoreCigarOp, CoreCigarRun, NormalizedSequence, RecordBuildError,
    RecordSegment, ReferenceIndex, ReferenceInterval, SAM_MAX_REFERENCE_LENGTH, VerifiedAlignment,
    evaluate_certified_ungapped_alignment, evaluate_ungapped_alignment,
    evaluate_verified_alignment, storage_count, storage_len, strand_semantics,
};

trait AlignmentCigarRunAdapter {
    fn from_core(run: CoreCigarRun) -> Self;
    fn all_match(length: u64) -> Self;
    fn soft_clip(length: u64) -> Self;
}

impl AlignmentCigarRunAdapter for AlignmentCigarRun {
    fn from_core(run: CoreCigarRun) -> Self {
        let operation = match run.operation() {
            CoreCigarOp::M => AlignmentCigarOp::Match,
            CoreCigarOp::I => AlignmentCigarOp::Insertion,
            CoreCigarOp::D => AlignmentCigarOp::Deletion,
        };
        Self::new(operation, run.length()).expect("validated core CIGAR run is positive")
    }

    fn all_match(length: u64) -> Self {
        Self::new(AlignmentCigarOp::Match, length).expect("validated aligned span is positive")
    }

    fn soft_clip(length: u64) -> Self {
        Self::new(AlignmentCigarOp::SoftClip, length).expect("validated soft clip is positive")
    }
}

/// Sequencing-order identity of one emitted record.
#[derive(Debug)]
struct DirectRecordOffsets {
    query_name_start: usize,
    query_name_len: usize,
    flag: u16,
    reference_ordinal: Option<u64>,
    position: u32,
    mapping_quality: u8,
    cigar_start: usize,
    cigar_len: usize,
    mate_reference_ordinal: Option<u64>,
    mate_position: u32,
    template_length: i32,
    sequence_start: usize,
    sequence_len: usize,
    quality_start: usize,
    quality_len: usize,
    literal_nm: u32,
    md_start: usize,
    md_len: usize,
    has_md: bool,
    strand: BisulfiteStrand,
    bismark_xm_start: usize,
    bismark_xm_len: usize,
    has_bismark_xm: bool,
}

/// Worker-local compact storage for retained direct-BAM records.
#[derive(Debug, Default)]
pub(crate) struct DirectRecordComposer {
    records: Vec<DirectRecordOffsets>,
    bytes: Vec<u8>,
    cigar_runs: Vec<AlignmentCigarRun>,
    md: Vec<u8>,
    bismark_xm: Vec<u8>,
}

pub(crate) type PairedRecordComposer = DirectRecordComposer;
pub(crate) type SingleRecordComposer = DirectRecordComposer;

impl DirectRecordComposer {
    /// Creates an empty worker-local batch.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            records: Vec::new(),
            bytes: Vec::new(),
            cigar_runs: Vec::new(),
            md: Vec::new(),
            bismark_xm: Vec::new(),
        }
    }

    /// Appends one unmapped unpaired record while preserving sequencing
    /// orientation and the complete input sequence and qualities.
    pub(crate) fn push_unmapped_single(
        &mut self,
        query_name: &[u8],
        read: BorrowedAlignmentRead<'_>,
        _limits: AlignmentRecordLimits,
    ) -> Result<(), RecordBuildError> {
        let query_name_start = append_pool_bytes(
            &mut self.bytes,
            query_name,
            AlignmentRecordAllocation::QueryName,
        )?;
        let sequence_start = append_oriented_read_sequence(
            &mut self.bytes,
            read.sequence(),
            AlignmentOrientation::Forward,
        )?;
        let (quality_start, quality_len) = append_oriented_quality(
            &mut self.bytes,
            Some(read.quality()),
            AlignmentOrientation::Forward,
        )?;
        self.records
            .try_reserve_exact(1)
            .map_err(|_| RecordBuildError::AllocationFailed {
                allocation: AlignmentRecordAllocation::Sequence,
                requested: 1,
            })?;
        self.records.push(DirectRecordOffsets {
            query_name_start,
            query_name_len: query_name.len(),
            flag: 0x4,
            reference_ordinal: None,
            position: 0,
            mapping_quality: 0,
            cigar_start: self.cigar_runs.len(),
            cigar_len: 0,
            mate_reference_ordinal: None,
            mate_position: 0,
            template_length: 0,
            sequence_start,
            sequence_len: read.sequence().len(),
            quality_start,
            quality_len,
            literal_nm: 0,
            md_start: 0,
            md_len: 0,
            has_md: false,
            strand: BisulfiteStrand::OT,
            bismark_xm_start: 0,
            bismark_xm_len: 0,
            has_bismark_xm: false,
        });
        Ok(())
    }

    /// Appends one mapped unpaired record from a retained contiguous query
    /// interval. The complete read is serialized and omitted terminal bases
    /// are represented as strand-correct soft clips.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn push_retained_single_with_mapping_quality(
        &mut self,
        reference: &ReferenceIndex,
        query_name: &[u8],
        full_read: BorrowedAlignmentRead<'_>,
        retained_range: core::ops::Range<usize>,
        retained_sequence: &NormalizedSequence,
        alignment: &VerifiedAlignment,
        limits: AlignmentRecordLimits,
        auxiliary_mode: AlignmentAuxiliaryMode,
        mapping_quality: u8,
    ) -> Result<(), RecordBuildError> {
        validate_soft_clipped_subsequence(full_read, retained_range.clone(), retained_sequence, 1)?;
        let query_name_start = append_pool_bytes(
            &mut self.bytes,
            query_name,
            AlignmentRecordAllocation::QueryName,
        )?;
        let prepared = prepare_direct_soft_clipped_mapping(
            &mut self.bytes,
            &mut self.cigar_runs,
            &mut self.md,
            &mut self.bismark_xm,
            reference,
            full_read,
            retained_range,
            retained_sequence,
            alignment,
            limits,
            auxiliary_mode,
        )?;
        self.records
            .try_reserve_exact(1)
            .map_err(|_| RecordBuildError::AllocationFailed {
                allocation: AlignmentRecordAllocation::Sequence,
                requested: 1,
            })?;
        self.records.push(prepared.finish_single(
            query_name_start,
            query_name.len(),
            mapping_quality,
        ));
        Ok(())
    }

    /// Appends one unpaired ungapped alignment directly from its selected
    /// placement. `false` leaves the batch unchanged and requests the
    /// traceback path.
    #[doc(hidden)]
    pub(crate) fn try_push_ungapped_single(
        &mut self,
        reference: &ReferenceIndex,
        query_name: &[u8],
        read: BorrowedAlignmentRead<'_>,
        placement: AlignmentPlacement,
        _limits: AlignmentRecordLimits,
        mapping_quality: u8,
    ) -> Result<bool, RecordBuildError> {
        let Some(inspected) = inspect_certified_ungapped_single(reference, read, placement)? else {
            return Ok(false);
        };
        let query_name_start = append_pool_bytes(
            &mut self.bytes,
            query_name,
            AlignmentRecordAllocation::QueryName,
        )?;
        let prepared =
            append_ungapped_mapping(&mut self.bytes, &mut self.cigar_runs, read, inspected)?;
        self.records
            .try_reserve_exact(1)
            .map_err(|_| RecordBuildError::AllocationFailed {
                allocation: AlignmentRecordAllocation::Sequence,
                requested: 1,
            })?;
        self.records.push(prepared.finish_single(
            query_name_start,
            query_name.len(),
            mapping_quality,
        ));
        Ok(true)
    }

    /// Soft-clipped counterpart of the unpaired ungapped fast path. The full
    /// read is retained in SEQ/QUAL and omitted terminal bases are represented
    /// as strand-correct soft clips.
    #[doc(hidden)]
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn try_push_soft_clipped_ungapped_single(
        &mut self,
        reference: &ReferenceIndex,
        query_name: &[u8],
        full_read: BorrowedAlignmentRead<'_>,
        retained_range: core::ops::Range<usize>,
        placement: AlignmentPlacement,
        _limits: AlignmentRecordLimits,
        mapping_quality: u8,
    ) -> Result<bool, RecordBuildError> {
        validate_soft_clipped_range(full_read, &retained_range, 1)?;
        let retained = BorrowedAlignmentRead::new(
            &full_read.sequence()[retained_range.clone()],
            &full_read.quality()[retained_range.clone()],
        );
        let Some(inspected) = inspect_certified_ungapped_single(reference, retained, placement)?
        else {
            return Ok(false);
        };
        let query_name_start = append_pool_bytes(
            &mut self.bytes,
            query_name,
            AlignmentRecordAllocation::QueryName,
        )?;
        let prepared = append_soft_clipped_ungapped_mapping(
            &mut self.bytes,
            &mut self.cigar_runs,
            full_read,
            retained_range,
            inspected,
        )?;
        self.records
            .try_reserve_exact(1)
            .map_err(|_| RecordBuildError::AllocationFailed {
                allocation: AlignmentRecordAllocation::Sequence,
                requested: 1,
            })?;
        self.records.push(prepared.finish_single(
            query_name_start,
            query_name.len(),
            mapping_quality,
        ));
        Ok(true)
    }

    /// Appends the two primary records for a paired template with no accepted
    /// placement. The input sequence and qualities remain in sequencing
    /// orientation, both records carry MAPQ 0, and neither record has a
    /// reference, mate coordinate, CIGAR, or alignment auxiliary fields.
    ///
    /// # Errors
    ///
    /// Returns the same validation and allocation failures as mapped direct
    /// record construction.
    #[doc(hidden)]
    pub(crate) fn push_unmapped_pair(
        &mut self,
        query_name: &[u8],
        first_read: BorrowedAlignmentRead<'_>,
        second_read: BorrowedAlignmentRead<'_>,
        _limits: AlignmentRecordLimits,
    ) -> Result<(), RecordBuildError> {
        let query_name_start = append_pool_bytes(
            &mut self.bytes,
            query_name,
            AlignmentRecordAllocation::QueryName,
        )?;
        let first_sequence_start = append_oriented_read_sequence(
            &mut self.bytes,
            first_read.sequence(),
            AlignmentOrientation::Forward,
        )?;
        let (first_quality_start, first_quality_len) = append_oriented_quality(
            &mut self.bytes,
            Some(first_read.quality()),
            AlignmentOrientation::Forward,
        )?;
        let second_sequence_start = append_oriented_read_sequence(
            &mut self.bytes,
            second_read.sequence(),
            AlignmentOrientation::Forward,
        )?;
        let (second_quality_start, second_quality_len) = append_oriented_quality(
            &mut self.bytes,
            Some(second_read.quality()),
            AlignmentOrientation::Forward,
        )?;
        self.records
            .try_reserve_exact(2)
            .map_err(|_| RecordBuildError::AllocationFailed {
                allocation: AlignmentRecordAllocation::Sequence,
                requested: 2,
            })?;
        let cigar_start = self.cigar_runs.len();
        let unmapped = |segment, sequence_start, sequence_len, quality_start, quality_len| {
            let segment_flag = match segment {
                RecordSegment::First => 0x40,
                RecordSegment::Last => 0x80,
                RecordSegment::Unpaired => 0,
            };
            DirectRecordOffsets {
                query_name_start,
                query_name_len: query_name.len(),
                flag: 0x1 | 0x4 | 0x8 | segment_flag,
                reference_ordinal: None,
                position: 0,
                mapping_quality: 0,
                cigar_start,
                cigar_len: 0,
                mate_reference_ordinal: None,
                mate_position: 0,
                template_length: 0,
                sequence_start,
                sequence_len,
                quality_start,
                quality_len,
                literal_nm: 0,
                md_start: 0,
                md_len: 0,
                has_md: false,
                strand: BisulfiteStrand::OT,
                bismark_xm_start: 0,
                bismark_xm_len: 0,
                has_bismark_xm: false,
            }
        };
        self.records.push(unmapped(
            RecordSegment::First,
            first_sequence_start,
            first_read.sequence().len(),
            first_quality_start,
            first_quality_len,
        ));
        self.records.push(unmapped(
            RecordSegment::Last,
            second_sequence_start,
            second_read.sequence().len(),
            second_quality_start,
            second_quality_len,
        ));
        Ok(())
    }

    /// Mapping-quality-aware variant used by the paired-end selector.
    #[doc(hidden)]
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn push_retained_unique_pair_with_mapping_quality(
        &mut self,
        reference: &ReferenceIndex,
        query_name: &[u8],
        first_read: AlignmentRead<'_>,
        second_read: AlignmentRead<'_>,
        first_alignment: &VerifiedAlignment,
        second_alignment: &VerifiedAlignment,
        limits: AlignmentRecordLimits,
        auxiliary_mode: AlignmentAuxiliaryMode,
        mapping_quality: u8,
    ) -> Result<(), RecordBuildError> {
        let query_name_start = append_pool_bytes(
            &mut self.bytes,
            query_name,
            AlignmentRecordAllocation::QueryName,
        )?;
        let first = prepare_direct_mapping(
            DirectMappingBuffers {
                bytes: &mut self.bytes,
                cigar_runs: &mut self.cigar_runs,
                md: &mut self.md,
                bismark_xm: &mut self.bismark_xm,
            },
            reference,
            first_read,
            first_alignment,
            limits,
            auxiliary_mode,
        )?;
        let second = prepare_direct_mapping(
            DirectMappingBuffers {
                bytes: &mut self.bytes,
                cigar_runs: &mut self.cigar_runs,
                md: &mut self.md,
                bismark_xm: &mut self.bismark_xm,
            },
            reference,
            second_read,
            second_alignment,
            limits,
            auxiliary_mode,
        )?;
        self.finish_direct_pair(
            query_name_start,
            query_name.len(),
            first,
            second,
            mapping_quality,
        )
    }

    /// Appends one unique pair aligned from arbitrary contiguous retained
    /// subsequences of the full reads. Omitted sequencing-orientation prefixes
    /// and suffixes are emitted at the strand-correct CIGAR ends.
    #[doc(hidden)]
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn push_soft_clipped_retained_unique_pair(
        &mut self,
        reference: &ReferenceIndex,
        query_name: &[u8],
        first_full_read: BorrowedAlignmentRead<'_>,
        second_full_read: BorrowedAlignmentRead<'_>,
        first_retained_range: core::ops::Range<usize>,
        second_retained_range: core::ops::Range<usize>,
        first_retained_sequence: &NormalizedSequence,
        second_retained_sequence: &NormalizedSequence,
        first_alignment: &VerifiedAlignment,
        second_alignment: &VerifiedAlignment,
        limits: AlignmentRecordLimits,
        auxiliary_mode: AlignmentAuxiliaryMode,
        mapping_quality: u8,
    ) -> Result<(), RecordBuildError> {
        validate_soft_clipped_subsequence(
            first_full_read,
            first_retained_range.clone(),
            first_retained_sequence,
            1,
        )?;
        validate_soft_clipped_subsequence(
            second_full_read,
            second_retained_range.clone(),
            second_retained_sequence,
            2,
        )?;
        let query_name_start = append_pool_bytes(
            &mut self.bytes,
            query_name,
            AlignmentRecordAllocation::QueryName,
        )?;
        let first = prepare_direct_soft_clipped_mapping(
            &mut self.bytes,
            &mut self.cigar_runs,
            &mut self.md,
            &mut self.bismark_xm,
            reference,
            first_full_read,
            first_retained_range,
            first_retained_sequence,
            first_alignment,
            limits,
            auxiliary_mode,
        )?;
        let second = prepare_direct_soft_clipped_mapping(
            &mut self.bytes,
            &mut self.cigar_runs,
            &mut self.md,
            &mut self.bismark_xm,
            reference,
            second_full_read,
            second_retained_range,
            second_retained_sequence,
            second_alignment,
            limits,
            auxiliary_mode,
        )?;
        self.finish_direct_pair(
            query_name_start,
            query_name.len(),
            first,
            second,
            mapping_quality,
        )
    }

    /// Mapping-quality-aware slab-backed ungapped path.
    #[doc(hidden)]
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn try_push_ungapped_pair(
        &mut self,
        reference: &ReferenceIndex,
        query_name: &[u8],
        first_read: BorrowedAlignmentRead<'_>,
        second_read: BorrowedAlignmentRead<'_>,
        first_placement: AlignmentPlacement,
        second_placement: AlignmentPlacement,
        _limits: AlignmentRecordLimits,
        mapping_quality: u8,
    ) -> Result<bool, RecordBuildError> {
        let Some(first) = inspect_ungapped_mapping(reference, first_read, first_placement)? else {
            return Ok(false);
        };
        let Some(second) = inspect_ungapped_mapping(reference, second_read, second_placement)?
        else {
            return Ok(false);
        };
        let query_name_start = append_pool_bytes(
            &mut self.bytes,
            query_name,
            AlignmentRecordAllocation::QueryName,
        )?;
        let first =
            append_ungapped_mapping(&mut self.bytes, &mut self.cigar_runs, first_read, first)?;
        let second =
            append_ungapped_mapping(&mut self.bytes, &mut self.cigar_runs, second_read, second)?;
        self.finish_direct_pair(
            query_name_start,
            query_name.len(),
            first,
            second,
            mapping_quality,
        )?;
        Ok(true)
    }

    /// Soft-clipped counterpart of the slab-backed ungapped fast path. It
    /// validates both retained subsequences against their selected placements
    /// and emits full-read `S/M/S` records without constructing traceback
    /// objects. `false` leaves the batch unchanged and requests the traceback
    /// path.
    #[doc(hidden)]
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn try_push_soft_clipped_ungapped_pair(
        &mut self,
        reference: &ReferenceIndex,
        query_name: &[u8],
        first_read: BorrowedAlignmentRead<'_>,
        second_read: BorrowedAlignmentRead<'_>,
        first_retained_range: core::ops::Range<usize>,
        second_retained_range: core::ops::Range<usize>,
        first_placement: AlignmentPlacement,
        second_placement: AlignmentPlacement,
        _limits: AlignmentRecordLimits,
        mapping_quality: u8,
    ) -> Result<bool, RecordBuildError> {
        validate_soft_clipped_range(first_read, &first_retained_range, 1)?;
        validate_soft_clipped_range(second_read, &second_retained_range, 2)?;
        let first_retained = BorrowedAlignmentRead::new(
            &first_read.sequence()[first_retained_range.clone()],
            &first_read.quality()[first_retained_range.clone()],
        );
        let second_retained = BorrowedAlignmentRead::new(
            &second_read.sequence()[second_retained_range.clone()],
            &second_read.quality()[second_retained_range.clone()],
        );
        let Some(first) = inspect_ungapped_mapping(reference, first_retained, first_placement)?
        else {
            return Ok(false);
        };
        let Some(second) = inspect_ungapped_mapping(reference, second_retained, second_placement)?
        else {
            return Ok(false);
        };
        let query_name_start = append_pool_bytes(
            &mut self.bytes,
            query_name,
            AlignmentRecordAllocation::QueryName,
        )?;
        let first = append_soft_clipped_ungapped_mapping(
            &mut self.bytes,
            &mut self.cigar_runs,
            first_read,
            first_retained_range,
            first,
        )?;
        let second = append_soft_clipped_ungapped_mapping(
            &mut self.bytes,
            &mut self.cigar_runs,
            second_read,
            second_retained_range,
            second,
        )?;
        self.finish_direct_pair(
            query_name_start,
            query_name.len(),
            first,
            second,
            mapping_quality,
        )?;
        Ok(true)
    }

    fn finish_direct_pair(
        &mut self,
        query_name_start: usize,
        query_name_len: usize,
        first: DirectPreparedMapping,
        second: DirectPreparedMapping,
        mapping_quality: u8,
    ) -> Result<(), RecordBuildError> {
        if first.reference_ordinal != second.reference_ordinal {
            return Err(RecordBuildError::ConcordantReferenceMismatch);
        }
        let (first_tlen, second_tlen) = direct_template_lengths(first.interval, second.interval)?;
        self.records
            .try_reserve_exact(2)
            .map_err(|_| RecordBuildError::AllocationFailed {
                allocation: AlignmentRecordAllocation::Sequence,
                requested: 2,
            })?;
        self.records.push(first.finish(
            query_name_start,
            query_name_len,
            RecordSegment::First,
            second.reference_ordinal,
            second.position,
            second.orientation,
            first_tlen,
            mapping_quality,
        ));
        self.records.push(second.finish(
            query_name_start,
            query_name_len,
            RecordSegment::Last,
            first.reference_ordinal,
            first.position,
            first.orientation,
            second_tlen,
            mapping_quality,
        ));
        Ok(())
    }

    pub(crate) fn flush_into(
        &mut self,
        batch: &mut AlignmentRecordBatch,
        limits: AlignmentRecordLimits,
    ) -> Result<(), RecordBuildError> {
        for record in &self.records {
            let query_name = &self.bytes
                [record.query_name_start..record.query_name_start + record.query_name_len];
            let cigar = &self.cigar_runs[record.cigar_start..record.cigar_start + record.cigar_len];
            let sequence =
                &self.bytes[record.sequence_start..record.sequence_start + record.sequence_len];
            let quality = (record.quality_len != 0).then(|| {
                &self.bytes[record.quality_start..record.quality_start + record.quality_len]
            });
            let md = record
                .has_md
                .then(|| &self.md[record.md_start..record.md_start + record.md_len]);
            let bismark_xm = record.has_bismark_xm.then(|| {
                &self.bismark_xm
                    [record.bismark_xm_start..record.bismark_xm_start + record.bismark_xm_len]
            });
            let record = BorrowedAlignmentRecord::new(
                query_name,
                record.flag,
                record.reference_ordinal,
                record.position,
                record.mapping_quality,
                cigar,
                record.mate_reference_ordinal,
                record.mate_position,
                record.template_length,
                sequence,
                quality,
                record.literal_nm,
                if record.has_md {
                    AlignmentAuxiliaryMode::Bismark
                } else {
                    AlignmentAuxiliaryMode::Minimal
                },
                md,
                record.strand,
                bismark_xm,
                limits,
            )?;
            batch.push(&record)?;
        }
        self.records.clear();
        self.bytes.clear();
        self.cigar_runs.clear();
        self.md.clear();
        self.bismark_xm.clear();
        Ok(())
    }
}

#[derive(Clone, Copy)]
struct InspectedUngappedMapping {
    reference_ordinal: u64,
    interval: ReferenceInterval,
    position: u32,
    orientation: AlignmentOrientation,
    strand: BisulfiteStrand,
    literal_nm: u32,
}

fn inspect_certified_ungapped_single(
    reference: &ReferenceIndex,
    read: BorrowedAlignmentRead<'_>,
    placement: AlignmentPlacement,
) -> Result<Option<InspectedUngappedMapping>, RecordBuildError> {
    inspect_ungapped_mapping_with_policy(reference, read, placement, true)
}

fn inspect_ungapped_mapping(
    reference: &ReferenceIndex,
    read: BorrowedAlignmentRead<'_>,
    placement: AlignmentPlacement,
) -> Result<Option<InspectedUngappedMapping>, RecordBuildError> {
    inspect_ungapped_mapping_with_policy(reference, read, placement, false)
}

fn inspect_ungapped_mapping_with_policy(
    reference: &ReferenceIndex,
    read: BorrowedAlignmentRead<'_>,
    placement: AlignmentPlacement,
    require_canonical_certificate: bool,
) -> Result<Option<InspectedUngappedMapping>, RecordBuildError> {
    if placement.interval().len() != storage_len(read.sequence().len()) {
        return Ok(None);
    }
    let Some(contig) = reference.contig_by_ordinal(placement.reference_ordinal()) else {
        let source = reference
            .contig_id(placement.reference_ordinal())
            .expect_err("absent paired-end ordinal must fail canonical lookup");
        return Err(RecordBuildError::ReferenceAccess { source });
    };
    let start = usize::try_from(placement.interval().start()).ok();
    let end = usize::try_from(placement.interval().end()).ok();
    let Some(reference_bases) = start
        .zip(end)
        .and_then(|(start, end)| contig.sequence().bases().get(start..end))
    else {
        return Ok(None);
    };
    let evaluation = if require_canonical_certificate {
        let Some(evaluation) = evaluate_certified_ungapped_alignment(
            reference_bases,
            read.sequence(),
            placement.strand(),
        )
        .map_err(|source| RecordBuildError::AlignmentEvaluation { source })?
        else {
            return Ok(None);
        };
        evaluation
    } else {
        evaluate_ungapped_alignment(reference_bases, read.sequence(), placement.strand())
            .map_err(|source| RecordBuildError::AlignmentEvaluation { source })?
    };
    if evaluation.distance().get() != u64::from(placement.distance()) {
        return Ok(None);
    }
    let literal_nm =
        u32::try_from(evaluation.literal_nm()).map_err(|_| RecordBuildError::FieldOutOfRange {
            field: AlignmentRecordField::Nm,
            value: evaluation.literal_nm(),
        })?;
    let semantics = strand_semantics(placement.strand());
    let position_u64 =
        placement
            .interval()
            .start()
            .checked_add(1)
            .ok_or(RecordBuildError::FieldOutOfRange {
                field: AlignmentRecordField::Position,
                value: placement.interval().start(),
            })?;
    let position = u32::try_from(position_u64).map_err(|_| RecordBuildError::FieldOutOfRange {
        field: AlignmentRecordField::Position,
        value: position_u64,
    })?;
    Ok(Some(InspectedUngappedMapping {
        reference_ordinal: placement.reference_ordinal(),
        interval: placement.interval(),
        position,
        orientation: semantics.orientation(),
        strand: placement.strand(),
        literal_nm,
    }))
}

fn append_ungapped_mapping(
    bytes: &mut Vec<u8>,
    cigar_runs: &mut Vec<AlignmentCigarRun>,
    read: BorrowedAlignmentRead<'_>,
    inspected: InspectedUngappedMapping,
) -> Result<DirectPreparedMapping, RecordBuildError> {
    let cigar_start = cigar_runs.len();
    cigar_runs
        .try_reserve_exact(1)
        .map_err(|_| RecordBuildError::AllocationFailed {
            allocation: AlignmentRecordAllocation::Cigar,
            requested: 1,
        })?;
    cigar_runs.push(AlignmentCigarRun::all_match(storage_len(
        read.sequence().len(),
    )));
    let sequence_start = bytes.len();
    bytes
        .try_reserve(read.sequence().len())
        .map_err(|_| RecordBuildError::AllocationFailed {
            allocation: AlignmentRecordAllocation::Sequence,
            requested: storage_len(read.sequence().len()),
        })?;
    match inspected.orientation {
        AlignmentOrientation::Forward => {
            bytes.extend(read.sequence().iter().map(|base| base.as_ascii()));
        }
        AlignmentOrientation::Reverse => {
            bytes.extend(
                read.sequence()
                    .iter()
                    .rev()
                    .map(|base| base.complement().as_ascii()),
            );
        }
    }
    let quality_start = bytes.len();
    bytes
        .try_reserve(read.quality().len())
        .map_err(|_| RecordBuildError::AllocationFailed {
            allocation: AlignmentRecordAllocation::Quality,
            requested: storage_len(read.quality().len()),
        })?;
    match inspected.orientation {
        AlignmentOrientation::Forward => bytes.extend_from_slice(read.quality()),
        AlignmentOrientation::Reverse => bytes.extend(read.quality().iter().rev().copied()),
    }
    Ok(DirectPreparedMapping {
        reference_ordinal: inspected.reference_ordinal,
        interval: inspected.interval,
        position: inspected.position,
        orientation: inspected.orientation,
        cigar_start,
        cigar_len: 1,
        sequence_start,
        sequence_len: read.sequence().len(),
        quality_start,
        quality_len: read.quality().len(),
        literal_nm: inspected.literal_nm,
        md_start: 0,
        md_len: 0,
        has_md: false,
        strand: inspected.strand,
        bismark_xm_start: 0,
        bismark_xm_len: 0,
        has_bismark_xm: false,
    })
}

fn append_soft_clipped_ungapped_mapping(
    bytes: &mut Vec<u8>,
    cigar_runs: &mut Vec<AlignmentCigarRun>,
    read: BorrowedAlignmentRead<'_>,
    retained_range: core::ops::Range<usize>,
    inspected: InspectedUngappedMapping,
) -> Result<DirectPreparedMapping, RecordBuildError> {
    let five_prime_clip = retained_range.start;
    let three_prime_clip = read.sequence().len().saturating_sub(retained_range.end);
    let (leading_clip, trailing_clip) = match inspected.orientation {
        AlignmentOrientation::Forward => (five_prime_clip, three_prime_clip),
        AlignmentOrientation::Reverse => (three_prime_clip, five_prime_clip),
    };
    let cigar_len = 1 + usize::from(leading_clip != 0) + usize::from(trailing_clip != 0);
    let cigar_start = cigar_runs.len();
    cigar_runs
        .try_reserve_exact(cigar_len)
        .map_err(|_| RecordBuildError::AllocationFailed {
            allocation: AlignmentRecordAllocation::Cigar,
            requested: storage_len(cigar_len),
        })?;
    if leading_clip != 0 {
        cigar_runs.push(AlignmentCigarRun::soft_clip(storage_len(leading_clip)));
    }
    cigar_runs.push(AlignmentCigarRun::all_match(storage_len(
        retained_range.end - retained_range.start,
    )));
    if trailing_clip != 0 {
        cigar_runs.push(AlignmentCigarRun::soft_clip(storage_len(trailing_clip)));
    }

    let sequence_start =
        append_oriented_read_sequence(bytes, read.sequence(), inspected.orientation)?;
    let (quality_start, quality_len) =
        append_oriented_quality(bytes, Some(read.quality()), inspected.orientation)?;
    Ok(DirectPreparedMapping {
        reference_ordinal: inspected.reference_ordinal,
        interval: inspected.interval,
        position: inspected.position,
        orientation: inspected.orientation,
        cigar_start,
        cigar_len,
        sequence_start,
        sequence_len: read.sequence().len(),
        quality_start,
        quality_len,
        literal_nm: inspected.literal_nm,
        md_start: 0,
        md_len: 0,
        has_md: false,
        strand: inspected.strand,
        bismark_xm_start: 0,
        bismark_xm_len: 0,
        has_bismark_xm: false,
    })
}

#[derive(Clone, Copy, Debug)]
struct DirectPreparedMapping {
    reference_ordinal: u64,
    interval: ReferenceInterval,
    position: u32,
    orientation: AlignmentOrientation,
    cigar_start: usize,
    cigar_len: usize,
    sequence_start: usize,
    sequence_len: usize,
    quality_start: usize,
    quality_len: usize,
    literal_nm: u32,
    md_start: usize,
    md_len: usize,
    has_md: bool,
    strand: BisulfiteStrand,
    bismark_xm_start: usize,
    bismark_xm_len: usize,
    has_bismark_xm: bool,
}

impl DirectPreparedMapping {
    fn finish_single(
        self,
        query_name_start: usize,
        query_name_len: usize,
        mapping_quality: u8,
    ) -> DirectRecordOffsets {
        let flag = if matches!(self.orientation, AlignmentOrientation::Reverse) {
            0x10
        } else {
            0
        };
        DirectRecordOffsets {
            query_name_start,
            query_name_len,
            flag,
            reference_ordinal: Some(self.reference_ordinal),
            position: self.position,
            mapping_quality,
            cigar_start: self.cigar_start,
            cigar_len: self.cigar_len,
            mate_reference_ordinal: None,
            mate_position: 0,
            template_length: 0,
            sequence_start: self.sequence_start,
            sequence_len: self.sequence_len,
            quality_start: self.quality_start,
            quality_len: self.quality_len,
            literal_nm: self.literal_nm,
            md_start: self.md_start,
            md_len: self.md_len,
            has_md: self.has_md,
            strand: self.strand,
            bismark_xm_start: self.bismark_xm_start,
            bismark_xm_len: self.bismark_xm_len,
            has_bismark_xm: self.has_bismark_xm,
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn finish(
        self,
        query_name_start: usize,
        query_name_len: usize,
        segment: RecordSegment,
        mate_reference_ordinal: u64,
        mate_position: u32,
        mate_orientation: AlignmentOrientation,
        template_length: i32,
        mapping_quality: u8,
    ) -> DirectRecordOffsets {
        let mut flag = 0x1 | 0x2;
        if matches!(self.orientation, AlignmentOrientation::Reverse) {
            flag |= 0x10;
        }
        if matches!(mate_orientation, AlignmentOrientation::Reverse) {
            flag |= 0x20;
        }
        flag |= match segment {
            RecordSegment::First => 0x40,
            RecordSegment::Last => 0x80,
            RecordSegment::Unpaired => 0,
        };
        DirectRecordOffsets {
            query_name_start,
            query_name_len,
            flag,
            reference_ordinal: Some(self.reference_ordinal),
            position: self.position,
            mapping_quality,
            cigar_start: self.cigar_start,
            cigar_len: self.cigar_len,
            mate_reference_ordinal: Some(mate_reference_ordinal),
            mate_position,
            template_length,
            sequence_start: self.sequence_start,
            sequence_len: self.sequence_len,
            quality_start: self.quality_start,
            quality_len: self.quality_len,
            literal_nm: self.literal_nm,
            md_start: self.md_start,
            md_len: self.md_len,
            has_md: self.has_md,
            strand: self.strand,
            bismark_xm_start: self.bismark_xm_start,
            bismark_xm_len: self.bismark_xm_len,
            has_bismark_xm: self.has_bismark_xm,
        }
    }
}

fn append_pool_bytes(
    output: &mut Vec<u8>,
    bytes: &[u8],
    allocation: AlignmentRecordAllocation,
) -> Result<usize, RecordBuildError> {
    let start = output.len();
    output
        .try_reserve(bytes.len())
        .map_err(|_| RecordBuildError::AllocationFailed {
            allocation,
            requested: storage_len(bytes.len()),
        })?;
    output.extend_from_slice(bytes);
    Ok(start)
}

fn append_oriented_sequence(
    output: &mut Vec<u8>,
    sequence: &NormalizedSequence,
    orientation: AlignmentOrientation,
) -> Result<usize, RecordBuildError> {
    let start = output.len();
    let capacity = storage_count(sequence.len(), AlignmentRecordAllocation::Sequence)?;
    output
        .try_reserve(capacity)
        .map_err(|_| RecordBuildError::AllocationFailed {
            allocation: AlignmentRecordAllocation::Sequence,
            requested: sequence.len(),
        })?;
    for index in 0..capacity {
        output.push(oriented_base(sequence, orientation, index).as_ascii());
    }
    Ok(start)
}

fn append_oriented_quality(
    output: &mut Vec<u8>,
    quality: Option<&[u8]>,
    orientation: AlignmentOrientation,
) -> Result<(usize, usize), RecordBuildError> {
    let Some(quality) = quality else {
        return Ok((0, 0));
    };
    let start = output.len();
    output
        .try_reserve(quality.len())
        .map_err(|_| RecordBuildError::AllocationFailed {
            allocation: AlignmentRecordAllocation::Quality,
            requested: storage_len(quality.len()),
        })?;
    match orientation {
        AlignmentOrientation::Forward => output.extend_from_slice(quality),
        AlignmentOrientation::Reverse => output.extend(quality.iter().rev().copied()),
    }
    Ok((start, quality.len()))
}

#[allow(clippy::too_many_arguments)]
fn validate_soft_clipped_subsequence(
    full_read: BorrowedAlignmentRead<'_>,
    retained_range: core::ops::Range<usize>,
    retained: &NormalizedSequence,
    mate: u8,
) -> Result<(), RecordBuildError> {
    validate_soft_clipped_range(full_read, &retained_range, mate)?;
    if retained.bases().len() != retained_range.end - retained_range.start
        || retained.bases() != &full_read.sequence()[retained_range]
    {
        return Err(RecordBuildError::SoftClippedSequenceMismatch { mate });
    }
    Ok(())
}

fn validate_soft_clipped_range(
    full_read: BorrowedAlignmentRead<'_>,
    retained_range: &core::ops::Range<usize>,
    mate: u8,
) -> Result<(), RecordBuildError> {
    if retained_range.start >= retained_range.end
        || retained_range.end > full_read.sequence().len()
        || retained_range.end > full_read.quality().len()
    {
        return Err(RecordBuildError::SoftClippedSequenceMismatch { mate });
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn prepare_direct_soft_clipped_mapping(
    bytes: &mut Vec<u8>,
    cigar_runs: &mut Vec<AlignmentCigarRun>,
    md: &mut Vec<u8>,
    bismark_xm: &mut Vec<u8>,
    reference: &ReferenceIndex,
    full_read: BorrowedAlignmentRead<'_>,
    retained_range: core::ops::Range<usize>,
    retained_sequence: &NormalizedSequence,
    alignment: &VerifiedAlignment,
    limits: AlignmentRecordLimits,
    auxiliary_mode: AlignmentAuxiliaryMode,
) -> Result<DirectPreparedMapping, RecordBuildError> {
    let five_prime_clip = retained_range.start;
    let three_prime_clip = full_read
        .sequence()
        .len()
        .saturating_sub(retained_range.end);
    let clip_runs = usize::from(five_prime_clip != 0) + usize::from(three_prime_clip != 0);
    if clip_runs != 0 {
        cigar_runs
            .try_reserve(clip_runs)
            .map_err(|_| RecordBuildError::AllocationFailed {
                allocation: AlignmentRecordAllocation::Cigar,
                requested: storage_len(clip_runs),
            })?;
    }
    let retained_quality = &full_read.quality()[retained_range.clone()];
    let mut prepared = prepare_direct_mapping(
        DirectMappingBuffers {
            bytes,
            cigar_runs,
            md,
            bismark_xm,
        },
        reference,
        AlignmentRead::new(retained_sequence, Some(retained_quality)),
        alignment,
        limits,
        auxiliary_mode,
    )?;

    debug_assert_eq!(
        bytes.len(),
        prepared.quality_start.saturating_add(prepared.quality_len)
    );
    bytes.truncate(prepared.sequence_start);
    prepared.sequence_start =
        append_oriented_read_sequence(bytes, full_read.sequence(), prepared.orientation)?;
    prepared.sequence_len = storage_count(
        storage_len(full_read.sequence().len()),
        AlignmentRecordAllocation::Sequence,
    )?;
    (prepared.quality_start, prepared.quality_len) =
        append_oriented_quality(bytes, Some(full_read.quality()), prepared.orientation)?;

    let (leading_clip, trailing_clip) = match prepared.orientation {
        AlignmentOrientation::Forward => (five_prime_clip, three_prime_clip),
        AlignmentOrientation::Reverse => (three_prime_clip, five_prime_clip),
    };
    if prepared.has_bismark_xm {
        debug_assert_eq!(
            bismark_xm.len(),
            prepared
                .bismark_xm_start
                .saturating_add(prepared.bismark_xm_len)
        );
        let additional = leading_clip.saturating_add(trailing_clip);
        bismark_xm
            .try_reserve(additional)
            .map_err(|_| RecordBuildError::AllocationFailed {
                allocation: AlignmentRecordAllocation::MethylationCall,
                requested: storage_len(full_read.sequence().len()),
            })?;
        bismark_xm.resize(
            prepared
                .bismark_xm_start
                .saturating_add(full_read.sequence().len()),
            b'.',
        );
        bismark_xm.copy_within(
            prepared.bismark_xm_start..prepared.bismark_xm_start + prepared.bismark_xm_len,
            prepared.bismark_xm_start + leading_clip,
        );
        bismark_xm[prepared.bismark_xm_start..prepared.bismark_xm_start + leading_clip].fill(b'.');
        bismark_xm[prepared.bismark_xm_start + leading_clip + prepared.bismark_xm_len
            ..prepared.bismark_xm_start + full_read.sequence().len()]
            .fill(b'.');
        prepared.bismark_xm_len = full_read.sequence().len();
    }
    if trailing_clip != 0 {
        cigar_runs.insert(
            prepared.cigar_start + prepared.cigar_len,
            AlignmentCigarRun::soft_clip(storage_len(trailing_clip)),
        );
        prepared.cigar_len += 1;
    }
    if leading_clip != 0 {
        cigar_runs.insert(
            prepared.cigar_start,
            AlignmentCigarRun::soft_clip(storage_len(leading_clip)),
        );
        prepared.cigar_len += 1;
    }
    Ok(prepared)
}

fn append_oriented_read_sequence(
    output: &mut Vec<u8>,
    sequence: &[Base],
    orientation: AlignmentOrientation,
) -> Result<usize, RecordBuildError> {
    let start = output.len();
    output
        .try_reserve(sequence.len())
        .map_err(|_| RecordBuildError::AllocationFailed {
            allocation: AlignmentRecordAllocation::Sequence,
            requested: storage_len(sequence.len()),
        })?;
    match orientation {
        AlignmentOrientation::Forward => {
            output.extend(sequence.iter().map(|base| base.as_ascii()));
        }
        AlignmentOrientation::Reverse => {
            output.extend(
                sequence
                    .iter()
                    .rev()
                    .map(|base| base.complement().as_ascii()),
            );
        }
    }
    Ok(start)
}

pub(super) fn authoritative_literal_nm(
    reference: &ReferenceIndex,
    read: &NormalizedSequence,
    alignment: &VerifiedAlignment,
) -> Result<u64, RecordBuildError> {
    if let Some(literal_nm) = alignment.cached_literal_nm() {
        return Ok(literal_nm);
    }
    evaluate_verified_alignment(reference, read, alignment)
        .map(bsbit_align::materialize::AlignmentEvaluation::literal_nm)
        .map_err(|source| RecordBuildError::AlignmentEvaluation { source })
}

struct DirectMappingBuffers<'a> {
    bytes: &'a mut Vec<u8>,
    cigar_runs: &'a mut Vec<AlignmentCigarRun>,
    md: &'a mut Vec<u8>,
    bismark_xm: &'a mut Vec<u8>,
}

fn sam_position(interval: &ReferenceInterval) -> Result<u32, RecordBuildError> {
    let position_u64 =
        interval
            .start()
            .checked_add(1)
            .ok_or(RecordBuildError::FieldOutOfRange {
                field: AlignmentRecordField::Position,
                value: interval.start(),
            })?;
    let position = u32::try_from(position_u64).map_err(|_| RecordBuildError::FieldOutOfRange {
        field: AlignmentRecordField::Position,
        value: position_u64,
    })?;
    if position_u64 > SAM_MAX_REFERENCE_LENGTH {
        return Err(RecordBuildError::FieldOutOfRange {
            field: AlignmentRecordField::Position,
            value: position_u64,
        });
    }
    Ok(position)
}

fn prepare_direct_mapping(
    buffers: DirectMappingBuffers<'_>,
    reference: &ReferenceIndex,
    read: AlignmentRead<'_>,
    alignment: &VerifiedAlignment,
    limits: AlignmentRecordLimits,
    auxiliary_mode: AlignmentAuxiliaryMode,
) -> Result<DirectPreparedMapping, RecordBuildError> {
    let DirectMappingBuffers {
        bytes,
        cigar_runs,
        md,
        bismark_xm,
    } = buffers;
    let contig = reference
        .resolve_contig(alignment.contig())
        .map_err(|source| RecordBuildError::ReferenceAccess { source })?;
    let interval = alignment.interval();
    let position = sam_position(&interval)?;
    let md_start = md.len();
    let bismark_xm_start = bismark_xm.len();
    let emit_md = matches!(auxiliary_mode, AlignmentAuxiliaryMode::Bismark);
    let emit_bismark = matches!(auxiliary_mode, AlignmentAuxiliaryMode::Bismark);
    let authoritative_literal_nm = authoritative_literal_nm(reference, read.sequence(), alignment)?;
    let summary = match auxiliary_mode {
        AlignmentAuxiliaryMode::Minimal => ReplaySummary {
            literal_nm: authoritative_literal_nm,
            md_bytes: 0,
            bismark_xm_bytes: 0,
        },
        AlignmentAuxiliaryMode::Bismark => replay_pass(
            contig.sequence(),
            interval,
            read.sequence(),
            alignment,
            emit_md.then_some(&mut *md),
            emit_bismark.then_some(&mut *bismark_xm),
            limits.max_md_bytes(),
        )?,
    };
    if summary.literal_nm != authoritative_literal_nm {
        return Err(RecordBuildError::AlignmentLiteralNmMismatch {
            expected: authoritative_literal_nm,
            observed: summary.literal_nm,
        });
    }
    let literal_nm =
        u32::try_from(summary.literal_nm).map_err(|_| RecordBuildError::FieldOutOfRange {
            field: AlignmentRecordField::Nm,
            value: summary.literal_nm,
        })?;
    let cigar_start = cigar_runs.len();
    cigar_runs
        .try_reserve(alignment.cigar().run_count())
        .map_err(|_| RecordBuildError::AllocationFailed {
            allocation: AlignmentRecordAllocation::Cigar,
            requested: storage_len(alignment.cigar().run_count()),
        })?;
    cigar_runs.extend(
        alignment
            .cigar()
            .runs()
            .iter()
            .copied()
            .map(AlignmentCigarRun::from_core),
    );
    let sequence_start = append_oriented_sequence(bytes, read.sequence(), alignment.orientation())?;
    let (quality_start, quality_len) =
        append_oriented_quality(bytes, read.quality(), alignment.orientation())?;
    let md_len = md.len().saturating_sub(md_start);
    debug_assert_eq!(storage_len(md_len), summary.md_bytes);
    let bismark_xm_len = bismark_xm.len().saturating_sub(bismark_xm_start);
    debug_assert_eq!(storage_len(bismark_xm_len), summary.bismark_xm_bytes);
    Ok(DirectPreparedMapping {
        reference_ordinal: contig.ordinal(),
        interval,
        position,
        orientation: alignment.orientation(),
        cigar_start,
        cigar_len: alignment.cigar().run_count(),
        sequence_start,
        sequence_len: storage_count(read.sequence().len(), AlignmentRecordAllocation::Sequence)?,
        quality_start,
        quality_len,
        literal_nm,
        md_start,
        md_len,
        has_md: emit_md,
        strand: alignment.strand(),
        bismark_xm_start,
        bismark_xm_len,
        has_bismark_xm: emit_bismark,
    })
}

fn direct_template_lengths(
    first: ReferenceInterval,
    second: ReferenceInterval,
) -> Result<(i32, i32), RecordBuildError> {
    let start = first.start().min(second.start());
    let end = first.end().max(second.end());
    let span = end
        .checked_sub(start)
        .ok_or(RecordBuildError::FieldOutOfRange {
            field: AlignmentRecordField::TemplateLength,
            value: end,
        })?;
    let magnitude = i32::try_from(span).map_err(|_| RecordBuildError::FieldOutOfRange {
        field: AlignmentRecordField::TemplateLength,
        value: span,
    })?;
    if first.start() <= second.start() {
        Ok((magnitude, -magnitude))
    } else {
        Ok((-magnitude, magnitude))
    }
}
