//! Application-layer mapping-result to SAM/BAM record composition.
//!
//! The aligner and index stay independent of output formats, while bsbit-hts
//! owns the shared SAM/BAM format model. This module is the explicit composition
//! boundary between those domains.

use core::fmt;
use core::mem::size_of;

use bsbit_align::extension::VerifiedAlignment;
use bsbit_align::materialize::{
    AlignmentEvaluationError, evaluate_ungapped_alignment, evaluate_verified_alignment,
};
use bsbit_core::alphabet::Base;
use bsbit_core::bisulfite::{
    AlignmentOrientation, BisulfiteStrand, CytosineStrand, strand_semantics,
};
use bsbit_core::cigar::{CoreCigarOp, CoreCigarRun};
use bsbit_core::coordinate::ReferenceInterval;
use bsbit_core::sequence::NormalizedSequence;
use bsbit_hts::{
    AlignmentAuxiliaryMode, AlignmentCigarOp, AlignmentCigarRun, AlignmentPlacement, AlignmentRead,
    AlignmentRecordAllocation, AlignmentRecordBatch,
    AlignmentRecordError as HtsAlignmentRecordError, AlignmentRecordField, AlignmentRecordLimits,
    AlignmentRecordResource, BorrowedAlignmentRead, BorrowedAlignmentRecord, RecordSegment,
    SAM_MAX_REFERENCE_LENGTH, SamHeader, SamHeaderReference,
};
#[cfg(test)]
use bsbit_hts::{
    AlignmentRecord, MappedAlignmentRecord, RecordMappingQuality, RecordMateLocation,
    RecordReference,
};
use bsbit_index::reference::{ReferenceAccessError, ReferenceIndex};

#[derive(Debug)]
pub(crate) enum RecordBuildError {
    /// A configured cap was exceeded.
    LimitExceeded {
        /// Controlled resource.
        resource: AlignmentRecordResource,
        /// First known complete value outside the cap.
        observed: u64,
        /// Configured cap.
        limit: u64,
    },
    /// A checked logical sum overflowed.
    ArithmeticOverflow {
        /// Resource being counted.
        resource: AlignmentRecordResource,
        /// Value before addition.
        current: u64,
        /// Attempted increment.
        increment: u64,
    },
    /// A bounded allocation could not be reserved.
    AllocationFailed {
        /// Allocation site.
        allocation: AlignmentRecordAllocation,
        /// Requested bytes or elements.
        requested: u64,
    },
    /// A retained soft-clipped query was not an exact prefix of its full read.
    SoftClippedSequenceMismatch {
        /// One-based mate number.
        mate: u8,
    },
    /// An owner-bound alignment did not resolve through the supplied reference.
    ReferenceAccess {
        /// Underlying exact-owner access failure.
        source: ReferenceAccessError,
    },
    /// A record integer cannot be represented in its SAM field.
    FieldOutOfRange {
        /// Failed field.
        field: AlignmentRecordField,
        /// Observed unsigned magnitude.
        value: u64,
    },
    /// Authoritative alignment-fact evaluation failed.
    AlignmentEvaluation {
        /// Alignment-layer failure.
        source: AlignmentEvaluationError,
    },
    /// Literal SAM replay disagreed with the alignment-layer NM fact.
    AlignmentLiteralNmMismatch {
        /// Alignment-layer literal NM.
        expected: u64,
        /// Format replay literal NM.
        observed: u64,
    },
    /// A concordant pair did not satisfy the same-reference record invariant.
    ConcordantReferenceMismatch,
    /// The shared HTS alignment-record model rejected composed format fields.
    Format { source: HtsAlignmentRecordError },
}

impl fmt::Display for RecordBuildError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::LimitExceeded {
                resource,
                observed,
                limit,
            } => write!(
                formatter,
                "alignment-record resource {resource:?} observed {observed}, exceeding {limit}"
            ),
            Self::ArithmeticOverflow {
                resource,
                current,
                increment,
            } => write!(
                formatter,
                "alignment-record resource {resource:?} overflowed: {current} + {increment}"
            ),
            Self::AllocationFailed {
                allocation,
                requested,
            } => write!(
                formatter,
                "failed to reserve {requested} bytes/elements for {allocation:?}"
            ),
            Self::SoftClippedSequenceMismatch { mate } => write!(
                formatter,
                "mate {mate} soft-clipped query is not an exact nonempty prefix of the full read"
            ),
            Self::ReferenceAccess { source } => {
                write!(formatter, "alignment reference access failed: {source}")
            }
            Self::FieldOutOfRange { field, value } => {
                write!(
                    formatter,
                    "alignment-record field {field:?} cannot represent {value}"
                )
            }
            Self::AlignmentEvaluation { source } => {
                write!(formatter, "alignment evaluation failed: {source}")
            }
            Self::AlignmentLiteralNmMismatch { expected, observed } => write!(
                formatter,
                "alignment literal NM {expected} differs from format replay NM {observed}"
            ),
            Self::ConcordantReferenceMismatch => {
                formatter.write_str("concordant pair resolved to different reference sequences")
            }
            Self::Format { source } => {
                write!(formatter, "alignment record validation failed: {source}")
            }
        }
    }
}

impl std::error::Error for RecordBuildError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::ReferenceAccess { source } => Some(source),
            Self::AlignmentEvaluation { source } => Some(source),
            Self::Format { source } => Some(source),
            _ => None,
        }
    }
}

impl From<HtsAlignmentRecordError> for RecordBuildError {
    fn from(source: HtsAlignmentRecordError) -> Self {
        match source {
            HtsAlignmentRecordError::LimitExceeded {
                resource,
                observed,
                limit,
            } => Self::LimitExceeded {
                resource,
                observed,
                limit,
            },
            HtsAlignmentRecordError::ArithmeticOverflow {
                resource,
                current,
                increment,
            } => Self::ArithmeticOverflow {
                resource,
                current,
                increment,
            },
            HtsAlignmentRecordError::AllocationFailed {
                allocation,
                requested,
            } => Self::AllocationFailed {
                allocation,
                requested,
            },
            HtsAlignmentRecordError::FieldOutOfRange { field, value } => {
                Self::FieldOutOfRange { field, value }
            }
            source => Self::Format { source },
        }
    }
}

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
        // With an equal-length query/reference span, a gapped traceback needs
        // at least one insertion and one deletion. Distances below two are
        // therefore the subset whose canonical traceback is necessarily
        // ungapped; higher distances retain the established tie-breaking path.
        if placement.distance() >= 2 {
            return Ok(false);
        }
        let Some(inspected) = inspect_ungapped_mapping(reference, read, placement)? else {
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
        if placement.distance() >= 2 {
            return Ok(false);
        }
        validate_soft_clipped_range(full_read, &retained_range, 1)?;
        let retained = BorrowedAlignmentRead::new(
            &full_read.sequence()[retained_range.clone()],
            &full_read.quality()[retained_range.clone()],
        );
        let Some(inspected) = inspect_ungapped_mapping(reference, retained, placement)? else {
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

fn inspect_ungapped_mapping(
    reference: &ReferenceIndex,
    read: BorrowedAlignmentRead<'_>,
    placement: AlignmentPlacement,
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
    let Ok(evaluation) =
        evaluate_ungapped_alignment(reference_bases, read.sequence(), placement.strand())
    else {
        return Ok(None);
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

fn authoritative_literal_nm(
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

pub(crate) fn build_sam_header(
    reference: &ReferenceIndex,
    limits: AlignmentRecordLimits,
) -> Result<SamHeader, RecordBuildError> {
    let capacity = usize::try_from(reference.contig_count()).map_err(|_| {
        RecordBuildError::AllocationFailed {
            allocation: AlignmentRecordAllocation::HeaderReferences,
            requested: reference.contig_count(),
        }
    })?;
    let mut entries = Vec::new();
    entries
        .try_reserve_exact(capacity)
        .map_err(|_| RecordBuildError::AllocationFailed {
            allocation: AlignmentRecordAllocation::HeaderReferences,
            requested: reference.contig_count(),
        })?;
    for ordinal in 0..reference.contig_count() {
        let id = reference
            .contig_id(ordinal)
            .map_err(|source| RecordBuildError::ReferenceAccess { source })?;
        let contig = reference
            .resolve_contig(&id)
            .map_err(|source| RecordBuildError::ReferenceAccess { source })?;
        entries.push(SamHeaderReference::new(
            ordinal,
            contig.name(),
            contig.sequence().len(),
        )?);
    }
    SamHeader::new(entries, limits).map_err(Into::into)
}

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
fn build_single_alignment_record_with_auxiliary_mode(
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

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct ReplaySummary {
    literal_nm: u64,
    md_bytes: u64,
    bismark_xm_bytes: u64,
}

// This is one bounds-checked pass over the CIGAR. Keeping its cursor and output
// state together avoids subtly divergent NM, MD, and XM traversal semantics.
#[allow(clippy::too_many_lines)]
fn replay_pass(
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

fn bismark_methylation_call(
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

fn oriented_base(
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

#[cfg(test)]
pub(crate) fn check_limit(
    observed: u64,
    limit: u64,
    resource: AlignmentRecordResource,
) -> Result<(), RecordBuildError> {
    if observed > limit {
        Err(RecordBuildError::LimitExceeded {
            resource,
            observed,
            limit,
        })
    } else {
        Ok(())
    }
}

pub(crate) fn checked_add_resource(
    current: u64,
    increment: u64,
    resource: AlignmentRecordResource,
) -> Result<u64, RecordBuildError> {
    current
        .checked_add(increment)
        .ok_or(RecordBuildError::ArithmeticOverflow {
            resource,
            current,
            increment,
        })
}

pub(crate) const fn decimal_digits(mut value: u64) -> u64 {
    let mut digits = 1;
    while value >= 10 {
        value /= 10;
        digits += 1;
    }
    digits
}

pub(crate) fn append_u64(output: &mut Vec<u8>, mut value: u64) {
    let mut digits = [0_u8; 20];
    let mut start = digits.len();
    loop {
        start -= 1;
        digits[start] = b'0' + u8::try_from(value % 10).expect("single decimal digit");
        value /= 10;
        if value == 0 {
            break;
        }
    }
    output.extend_from_slice(&digits[start..]);
}

pub(crate) const fn storage_len(length: usize) -> u64 {
    length as u64
}

fn storage_count(
    value: u64,
    allocation: AlignmentRecordAllocation,
) -> Result<usize, RecordBuildError> {
    usize::try_from(value).map_err(|_| RecordBuildError::AllocationFailed {
        allocation,
        requested: value.saturating_mul(size_of::<u8>() as u64),
    })
}

#[cfg(test)]
#[path = "../tests/whitebox/alignment_record.rs"]
mod tests;

#[cfg(test)]
#[path = "../tests/whitebox/alignment_record_fuzz_smoke.rs"]
mod fuzz_tests;

#[cfg(test)]
#[path = "../tests/whitebox/bam_fields.rs"]
mod bam_tests;

#[cfg(test)]
#[path = "../tests/whitebox/record_fixture.rs"]
mod record_fixture;
