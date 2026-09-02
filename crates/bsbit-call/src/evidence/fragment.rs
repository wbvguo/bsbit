//! BAM alignment reconstruction and fragment-level evidence normalization.

use std::collections::HashMap;

use bsbit_hts::{BamAlignmentColumn, BamRecordDecodeWorkspace, IndexedBamReader, IndexedBamRecord};

use super::{BaseCode, CytosineContext, EvidenceObservation, EvidenceStrand};
use crate::call_input::BamReference;
use crate::reference_context::ReferenceWindow;
use crate::region::CallRegion;
use crate::{CallError, CallErrorKind};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum MateSegment {
    First,
    Second,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum ReconstructedRecord {
    Ignored,
    Unpaired,
    Paired {
        segment: MateSegment,
        alignment_end: u32,
    },
}

#[derive(Clone, Copy)]
pub(crate) enum EvidenceContext<'a> {
    WithCytosineContext(&'a ReferenceWindow),
    WithoutCytosineContext(&'a ReferenceWindow),
}

#[derive(Clone, Copy)]
struct EvidenceRecordFields {
    reference: u32,
    core_region: (u32, u32),
    mapping_quality: u8,
    strand: EvidenceStrand,
}

#[derive(Debug)]
struct PendingMate {
    read_group: Option<Vec<u8>>,
    segment: MateSegment,
    reference_id: i32,
    position: i64,
    mate_reference_id: i32,
    mate_position: i64,
    observations: Vec<EvidenceObservation>,
}

#[derive(Debug, Default)]
pub(crate) struct EvidenceWorkspace {
    record: IndexedBamRecord,
    evidence: EvidenceDecodeWorkspace,
}

#[derive(Debug, Default)]
struct EvidenceDecodeWorkspace {
    bam: BamRecordDecodeWorkspace,
    observations: Vec<EvidenceObservation>,
    merged: Vec<EvidenceObservation>,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct EvidenceFilter {
    minimum_base_quality: u8,
    minimum_mapping_quality: u8,
    retain_deletions: bool,
    ignore_orphans: bool,
}

impl EvidenceFilter {
    pub(crate) const fn new(
        minimum_base_quality: u8,
        minimum_mapping_quality: u8,
        retain_deletions: bool,
        ignore_orphans: bool,
    ) -> Self {
        Self {
            minimum_base_quality,
            minimum_mapping_quality,
            retain_deletions,
            ignore_orphans,
        }
    }
}

pub(crate) fn for_each_region_fragment(
    reader: &mut IndexedBamReader,
    references: &[BamReference],
    region: CallRegion,
    evidence_filter: EvidenceFilter,
    evidence_context: EvidenceContext<'_>,
    workspace: &mut EvidenceWorkspace,
    consume: impl FnMut(&[EvidenceObservation]) -> Result<(), CallError>,
) -> Result<(), CallError> {
    reader
        .query(
            region.reference,
            u64::from(region.start),
            u64::from(region.end),
        )
        .map_err(|error| {
            CallError::with_source(
                CallErrorKind::Input,
                format!(
                    "query indexed BAM reference {}:{}-{}",
                    region.reference, region.start, region.end
                ),
                error,
            )
        })?;
    let EvidenceWorkspace { record, evidence } = workspace;
    let mut collector = RegionFragmentCollector::new(region, evidence_filter, evidence, consume);
    let mut ordinal = 0_u64;
    loop {
        let has_record = reader.next_record_into(record).map_err(|error| {
            CallError::with_source(
                CallErrorKind::Input,
                "read indexed BAM region record",
                error,
            )
        })?;
        if !has_record {
            break;
        }
        if should_ignore_orphan(record.flag(), evidence_filter.ignore_orphans) {
            continue;
        }
        ordinal = ordinal
            .checked_add(1)
            .ok_or_else(|| CallError::operation("region BAM record ordinal overflowed u64"))?;
        let reconstructed = reconstruct_indexed_evidence(
            record,
            references,
            (region.start, region.end),
            evidence_context,
            &mut *collector.workspace,
        )
        .map_err(|error| error.with_context(format!("region record {ordinal}")))?;
        collector.accept(record, reconstructed, ordinal)?;
    }
    collector.finish()
}

const fn should_ignore_orphan(flag: u16, ignore_orphans: bool) -> bool {
    ignore_orphans && flag & 0x1 != 0 && flag & 0x2 == 0
}

fn reconstruct_indexed_evidence(
    record: &IndexedBamRecord,
    references: &[BamReference],
    core_region: (u32, u32),
    evidence_context: EvidenceContext<'_>,
    workspace: &mut EvidenceDecodeWorkspace,
) -> Result<ReconstructedRecord, CallError> {
    workspace.observations.clear();
    let flag = record.flag();
    if flag & (0x100 | 0x200 | 0x400 | 0x800) != 0 || flag & 0x04 != 0 {
        return Ok(ReconstructedRecord::Ignored);
    }
    let segment = mate_segment(flag)?;
    let reference_ordinal = usize::try_from(record.reference_id())
        .ok()
        .filter(|ordinal| *ordinal < references.len())
        .ok_or_else(|| {
            CallError::input(format!(
                "mapped reference id {} is outside the BAM dictionary",
                record.reference_id()
            ))
        })?;
    let genome_conversion = record
        .string_auxiliary(*b"XG")
        .map_err(|error| {
            CallError::with_source(CallErrorKind::Input, "decode BAM XG auxiliary field", error)
        })?
        .ok_or_else(|| CallError::input("mapped primary record has no required XG:Z tag"))?;
    let strand = evidence_strand_from_xg(genome_conversion)?;
    let columns = record
        .project_alignment_into(references[reference_ordinal].length, &mut workspace.bam)
        .map_err(|error| {
            CallError::with_source(
                CallErrorKind::Input,
                "reconstruct BAM alignment evidence",
                error,
            )
        })?;
    let alignment_end = columns
        .last()
        .map(|column| {
            column
                .position
                .checked_add(1)
                .ok_or_else(|| CallError::input("BAM alignment end overflowed u32"))
        })
        .transpose()?
        .unwrap_or_else(|| u32::try_from(record.position()).unwrap_or(0));
    let reference = u32::try_from(reference_ordinal)
        .map_err(|_| CallError::input("reference ordinal does not fit u32"))?;
    workspace
        .observations
        .try_reserve(columns.len())
        .map_err(|error| {
            CallError::with_source(
                CallErrorKind::Calling,
                format!("reserve {} evidence columns", columns.len()),
                error,
            )
        })?;
    let fields = EvidenceRecordFields {
        reference,
        core_region,
        mapping_quality: record.mapping_quality(),
        strand,
    };
    append_evidence_for_context(
        columns,
        fields,
        evidence_context,
        &mut workspace.observations,
    )?;
    Ok(segment.map_or(ReconstructedRecord::Unpaired, |segment| {
        ReconstructedRecord::Paired {
            segment,
            alignment_end,
        }
    }))
}

fn append_evidence_for_context(
    columns: &[BamAlignmentColumn],
    fields: EvidenceRecordFields,
    evidence_context: EvidenceContext<'_>,
    observations: &mut Vec<EvidenceObservation>,
) -> Result<(), CallError> {
    match evidence_context {
        EvidenceContext::WithCytosineContext(window) => {
            append_reference_evidence::<true>(columns, fields, window, observations)?;
        }
        EvidenceContext::WithoutCytosineContext(window) => {
            append_reference_evidence::<false>(columns, fields, window, observations)?;
        }
    }
    Ok(())
}

fn append_reference_evidence<const INCLUDE_CYTOSINE_CONTEXT: bool>(
    columns: &[BamAlignmentColumn],
    fields: EvidenceRecordFields,
    reference_window: &ReferenceWindow,
    observations: &mut Vec<EvidenceObservation>,
) -> Result<(), CallError> {
    for column in columns {
        if column.position < fields.core_region.0 || column.position >= fields.core_region.1 {
            continue;
        }
        let fasta_base = reference_window
            .base(fields.reference, column.position)
            .ok_or_else(|| {
                CallError::operation(format!(
                    "reference context window does not cover position {}:{}",
                    fields.reference, column.position
                ))
            })?;
        observations.push(EvidenceObservation {
            reference: fields.reference,
            position: column.position,
            reference_base: fasta_base,
            query_base: column.query_base,
            base_quality: column.query_quality,
            mapping_quality: fields.mapping_quality,
            strand: fields.strand,
            context: if INCLUDE_CYTOSINE_CONTEXT {
                reference_window.context(fields.reference, column.position, fields.strand)
            } else {
                None
            },
        });
    }
    Ok(())
}

struct RegionFragmentCollector<'a, F> {
    region: CallRegion,
    evidence_filter: EvidenceFilter,
    workspace: &'a mut EvidenceDecodeWorkspace,
    pending_mates: HashMap<Vec<u8>, Vec<PendingMate>>,
    consume: F,
}

impl<'a, F> RegionFragmentCollector<'a, F>
where
    F: FnMut(&[EvidenceObservation]) -> Result<(), CallError>,
{
    fn new(
        region: CallRegion,
        evidence_filter: EvidenceFilter,
        workspace: &'a mut EvidenceDecodeWorkspace,
        consume: F,
    ) -> Self {
        Self {
            region,
            evidence_filter,
            workspace,
            pending_mates: HashMap::new(),
            consume,
        }
    }

    fn accept(
        &mut self,
        record: &IndexedBamRecord,
        reconstructed: ReconstructedRecord,
        ordinal: u64,
    ) -> Result<(), CallError> {
        match reconstructed {
            ReconstructedRecord::Ignored => Ok(()),
            ReconstructedRecord::Unpaired => self.consume_current(ordinal),
            ReconstructedRecord::Paired {
                segment,
                alignment_end,
            } => self.accept_paired(record, segment, alignment_end, ordinal),
        }
    }

    fn accept_paired(
        &mut self,
        record: &IndexedBamRecord,
        segment: MateSegment,
        alignment_end: u32,
        ordinal: u64,
    ) -> Result<(), CallError> {
        let read_group = record.string_auxiliary(*b"RG").map_err(|error| {
            CallError::with_source(
                CallErrorKind::Input,
                format!("region record {ordinal}: decode RG auxiliary field"),
                error,
            )
        })?;
        if let Some(previous) =
            take_matching_pending_mate(&mut self.pending_mates, record, read_group)
        {
            return self.consume_pair(record, previous, segment, ordinal);
        }
        if may_overlap_future_mate(record, alignment_end, self.region)
            && !self.workspace.observations.is_empty()
        {
            self.stash_pending(record, read_group, segment)
        } else {
            self.consume_current(ordinal)
        }
    }

    fn consume_pair(
        &mut self,
        record: &IndexedBamRecord,
        previous: PendingMate,
        segment: MateSegment,
        ordinal: u64,
    ) -> Result<(), CallError> {
        if previous.segment == segment {
            return Err(CallError::input(format!(
                "region record {ordinal}: paired QNAME `{}` repeats the same segment",
                String::from_utf8_lossy(record.query_name())
            )));
        }
        let current = std::mem::take(&mut self.workspace.observations);
        merge_evidence_pair_into(
            previous.segment,
            &previous.observations,
            segment,
            &current,
            self.evidence_filter,
            &mut self.workspace.merged,
        )?;
        (self.consume)(&self.workspace.merged)
            .map_err(|error| error.with_context(format!("region record {ordinal}")))?;
        self.workspace.merged.clear();
        self.workspace.observations = if previous.observations.capacity() >= current.capacity() {
            previous.observations
        } else {
            current
        };
        self.workspace.observations.clear();
        Ok(())
    }

    fn stash_pending(
        &mut self,
        record: &IndexedBamRecord,
        read_group: Option<&[u8]>,
        segment: MateSegment,
    ) -> Result<(), CallError> {
        self.pending_mates.try_reserve(1).map_err(|error| {
            CallError::with_source(
                CallErrorKind::Calling,
                "reserve pending overlapping-mate map",
                error,
            )
        })?;
        let query_name = copy_evidence_key(record.query_name(), "overlapping-mate query name")?;
        let read_group = read_group
            .map(|value| copy_evidence_key(value, "overlapping-mate read group"))
            .transpose()?;
        let bucket = self.pending_mates.entry(query_name).or_default();
        bucket.try_reserve(1).map_err(|error| {
            CallError::with_source(
                CallErrorKind::Calling,
                "reserve pending overlapping-mate bucket",
                error,
            )
        })?;
        bucket.push(PendingMate {
            read_group,
            segment,
            reference_id: record.reference_id(),
            position: record.position(),
            mate_reference_id: record.mate_reference_id(),
            mate_position: record.mate_position(),
            observations: std::mem::take(&mut self.workspace.observations),
        });
        Ok(())
    }

    fn consume_current(&mut self, ordinal: u64) -> Result<(), CallError> {
        (self.consume)(&self.workspace.observations)
            .map_err(|error| error.with_context(format!("region record {ordinal}")))
    }

    fn finish(mut self) -> Result<(), CallError> {
        for bucket in self.pending_mates.into_values() {
            for pending in bucket {
                (self.consume)(&pending.observations).map_err(|error| {
                    error.with_context("aggregate unmatched paired region record")
                })?;
            }
        }
        Ok(())
    }
}

fn take_matching_pending_mate(
    pending_mates: &mut HashMap<Vec<u8>, Vec<PendingMate>>,
    record: &IndexedBamRecord,
    read_group: Option<&[u8]>,
) -> Option<PendingMate> {
    let bucket = pending_mates.get_mut(record.query_name())?;
    let index = bucket
        .iter()
        .position(|pending| pending_matches_record(pending, record, read_group))?;
    let pending = bucket.swap_remove(index);
    let remove_bucket = bucket.is_empty();
    if remove_bucket {
        pending_mates.remove(record.query_name());
    }
    Some(pending)
}

fn copy_evidence_key(bytes: &[u8], label: &str) -> Result<Vec<u8>, CallError> {
    let mut copy = Vec::new();
    copy.try_reserve_exact(bytes.len()).map_err(|error| {
        CallError::with_source(
            CallErrorKind::Calling,
            format!("reserve {label} bytes"),
            error,
        )
    })?;
    copy.extend_from_slice(bytes);
    Ok(copy)
}

fn pending_matches_record(
    pending: &PendingMate,
    record: &IndexedBamRecord,
    read_group: Option<&[u8]>,
) -> bool {
    pending.read_group.as_deref() == read_group
        && pending.reference_id == record.mate_reference_id()
        && pending.position == record.mate_position()
        && pending.mate_reference_id == record.reference_id()
        && pending.mate_position == record.position()
}

fn may_overlap_future_mate(
    record: &IndexedBamRecord,
    alignment_end: u32,
    region: CallRegion,
) -> bool {
    if record.flag() & 0x08 != 0 || record.reference_id() != record.mate_reference_id() {
        return false;
    }
    let Ok(alignment_start) = u32::try_from(record.position()) else {
        return false;
    };
    let Ok(mate_start) = u32::try_from(record.mate_position()) else {
        return false;
    };
    mate_start >= alignment_start
        && mate_start < alignment_end
        && mate_start < region.end
        && alignment_end > region.start
}

fn merge_evidence_pair_into(
    left_segment: MateSegment,
    left: &[EvidenceObservation],
    right_segment: MateSegment,
    right: &[EvidenceObservation],
    evidence_filter: EvidenceFilter,
    merged: &mut Vec<EvidenceObservation>,
) -> Result<(), CallError> {
    debug_assert_ne!(left_segment, right_segment);
    merged.clear();
    merged
        .try_reserve(left.len().saturating_add(right.len()))
        .map_err(|error| {
            CallError::with_source(
                CallErrorKind::Calling,
                "reserve overlapping-pair evidence",
                error,
            )
        })?;
    let mut left_index = 0_usize;
    let mut right_index = 0_usize;
    while left_index < left.len() || right_index < right.len() {
        match (left.get(left_index), right.get(right_index)) {
            (Some(left_observation), Some(right_observation)) => {
                match left_observation.key().cmp(&right_observation.key()) {
                    std::cmp::Ordering::Less => {
                        merged.push(*left_observation);
                        left_index += 1;
                    }
                    std::cmp::Ordering::Greater => {
                        merged.push(*right_observation);
                        right_index += 1;
                    }
                    std::cmp::Ordering::Equal => {
                        if left_observation.reference_base != right_observation.reference_base
                            || left_observation.strand != right_observation.strand
                        {
                            return Err(CallError::input(format!(
                                "inconsistent overlapping-pair reference evidence at reference {}, position {}",
                                left_observation.reference, left_observation.position
                            )));
                        }
                        let context = merge_contexts(
                            left_observation.reference,
                            left_observation.position,
                            left_observation.context,
                            right_observation.context,
                        )?;
                        let mut selected = if prefer_left_evidence(
                            left_segment,
                            *left_observation,
                            right_segment,
                            *right_observation,
                            evidence_filter,
                        ) {
                            *left_observation
                        } else {
                            *right_observation
                        };
                        selected.context = context;
                        merged.push(selected);
                        left_index += 1;
                        right_index += 1;
                    }
                }
            }
            (Some(observation), None) => {
                merged.push(*observation);
                left_index += 1;
            }
            (None, Some(observation)) => {
                merged.push(*observation);
                right_index += 1;
            }
            (None, None) => break,
        }
    }
    Ok(())
}

fn prefer_left_evidence(
    left_segment: MateSegment,
    left: EvidenceObservation,
    right_segment: MateSegment,
    right: EvidenceObservation,
    evidence_filter: EvidenceFilter,
) -> bool {
    let structural_priority = |observation: EvidenceObservation| {
        (
            observation_passes_filter(observation, evidence_filter),
            observation
                .query_base
                .and_then(BaseCode::from_ascii)
                .is_some(),
            observation.query_base.is_some(),
            observation.base_quality.is_some(),
        )
    };
    let left_priority = structural_priority(left);
    let right_priority = structural_priority(right);
    if left_priority != right_priority {
        return left_priority > right_priority;
    }
    let prefer_first_segment =
        left_segment == MateSegment::First && right_segment == MateSegment::Second;
    match (
        effective_observation_error(left),
        effective_observation_error(right),
    ) {
        (Some(left_error), Some(right_error)) => match left_error.total_cmp(&right_error) {
            std::cmp::Ordering::Less => true,
            std::cmp::Ordering::Greater => false,
            std::cmp::Ordering::Equal => prefer_first_segment,
        },
        (Some(_), None) => true,
        (None, Some(_)) => false,
        (None, None) => prefer_first_segment,
    }
}

fn observation_passes_filter(
    observation: EvidenceObservation,
    evidence_filter: EvidenceFilter,
) -> bool {
    if observation.mapping_quality == u8::MAX
        || observation.mapping_quality < evidence_filter.minimum_mapping_quality
    {
        return false;
    }
    match observation.query_base {
        Some(_) => observation
            .base_quality
            .is_some_and(|quality| quality >= evidence_filter.minimum_base_quality),
        None => evidence_filter.retain_deletions,
    }
}

fn effective_observation_error(observation: EvidenceObservation) -> Option<f64> {
    let mapping_quality =
        (observation.mapping_quality != u8::MAX).then_some(observation.mapping_quality)?;
    observation.base_quality.map_or_else(
        || Some(10_f64.powf(-f64::from(mapping_quality.min(60)) / 10.0)),
        |base_quality| Some(combined_observation_error(base_quality, mapping_quality)),
    )
}

pub(crate) fn combined_observation_error(base_quality: u8, mapping_quality: u8) -> f64 {
    let base_error = 10_f64.powf(-f64::from(base_quality.min(60)) / 10.0);
    let mapping_error = 10_f64.powf(-f64::from(mapping_quality.min(60)) / 10.0);
    (1.0 - (1.0 - base_error) * (1.0 - mapping_error)).min(0.75)
}

fn mate_segment(flag: u16) -> Result<Option<MateSegment>, CallError> {
    let paired = flag & 0x1 != 0;
    let first = flag & 0x40 != 0;
    let second = flag & 0x80 != 0;
    if paired && first == second {
        return Err(CallError::input(
            "paired primary record must set exactly one of first/last segment FLAG",
        ));
    }
    if !paired && (first || second) {
        return Err(CallError::input(
            "unpaired primary record sets a paired segment FLAG",
        ));
    }
    Ok(if first {
        Some(MateSegment::First)
    } else if second {
        Some(MateSegment::Second)
    } else {
        None
    })
}

fn evidence_strand_from_xg(xg: &[u8]) -> Result<EvidenceStrand, CallError> {
    match xg {
        b"CT" => Ok(EvidenceStrand::Top),
        b"GA" => Ok(EvidenceStrand::Bottom),
        _ => Err(CallError::input(format!(
            "XG:Z must be CT or GA, observed `{}`",
            String::from_utf8_lossy(xg)
        ))),
    }
}

pub(crate) fn merge_contexts(
    reference: u32,
    position: u32,
    left: Option<CytosineContext>,
    right: Option<CytosineContext>,
) -> Result<Option<CytosineContext>, CallError> {
    match (left, right) {
        (Some(left), Some(right)) if left != right => Err(CallError::input(format!(
            "inconsistent reconstructed context at reference {reference}, position {position}: {left:?} versus {right:?}"
        ))),
        (Some(context), _) | (_, Some(context)) => Ok(Some(context)),
        (None, None) => Ok(None),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn overlap_evidence(query_base: Option<u8>, quality: Option<u8>) -> EvidenceObservation {
        EvidenceObservation {
            reference: 0,
            position: 10,
            reference_base: b'C',
            query_base,
            base_quality: quality,
            mapping_quality: 60,
            strand: EvidenceStrand::Top,
            context: None,
        }
    }

    const fn permissive_evidence_filter() -> EvidenceFilter {
        EvidenceFilter::new(0, 0, true, false)
    }

    fn merge_evidence_pair(
        left_segment: MateSegment,
        left: &[EvidenceObservation],
        right_segment: MateSegment,
        right: &[EvidenceObservation],
        evidence_filter: EvidenceFilter,
    ) -> Result<Vec<EvidenceObservation>, CallError> {
        let mut merged = Vec::new();
        merge_evidence_pair_into(
            left_segment,
            left,
            right_segment,
            right,
            evidence_filter,
            &mut merged,
        )?;
        Ok(merged)
    }

    #[test]
    fn xg_defines_bisulfite_evidence_strand() {
        assert_eq!(evidence_strand_from_xg(b"CT").unwrap(), EvidenceStrand::Top);
        assert_eq!(
            evidence_strand_from_xg(b"GA").unwrap(),
            EvidenceStrand::Bottom
        );
        assert!(evidence_strand_from_xg(b"CA").is_err());
    }

    #[test]
    fn overlap_collapse_uses_eligibility_combined_quality_then_r1_once() {
        let r1 = overlap_evidence(Some(b'C'), Some(20));
        let r2 = overlap_evidence(Some(b'T'), Some(30));
        let merged = merge_evidence_pair(
            MateSegment::First,
            &[r1],
            MateSegment::Second,
            &[r2],
            permissive_evidence_filter(),
        )
        .unwrap();
        assert_eq!(merged.len(), 1);
        assert_eq!(merged[0].query_base, Some(b'T'));

        let tied_r2 = overlap_evidence(Some(b'T'), Some(20));
        let merged = merge_evidence_pair(
            MateSegment::First,
            &[r1],
            MateSegment::Second,
            &[tied_r2],
            permissive_evidence_filter(),
        )
        .unwrap();
        assert_eq!(merged, vec![r1]);

        let deletion = overlap_evidence(None, None);
        let missing_quality = overlap_evidence(Some(b'T'), None);
        let merged = merge_evidence_pair(
            MateSegment::First,
            &[deletion],
            MateSegment::Second,
            &[missing_quality],
            permissive_evidence_filter(),
        )
        .unwrap();
        assert_eq!(merged, vec![deletion]);

        let ambiguous = overlap_evidence(Some(b'N'), Some(40));
        let canonical = overlap_evidence(Some(b'C'), Some(20));
        let merged = merge_evidence_pair(
            MateSegment::First,
            &[ambiguous],
            MateSegment::Second,
            &[canonical],
            permissive_evidence_filter(),
        )
        .unwrap();
        assert_eq!(merged, vec![canonical]);

        let mut low_mapping_quality = overlap_evidence(Some(b'C'), Some(40));
        low_mapping_quality.mapping_quality = 0;
        let high_confidence = overlap_evidence(Some(b'T'), Some(30));
        let merged = merge_evidence_pair(
            MateSegment::First,
            &[low_mapping_quality],
            MateSegment::Second,
            &[high_confidence],
            EvidenceFilter::new(15, 20, true, false),
        )
        .unwrap();
        assert_eq!(merged, vec![high_confidence]);
    }

    #[test]
    fn orphan_filter_retains_single_reads_and_proper_pairs() {
        assert!(!should_ignore_orphan(0, true));
        assert!(!should_ignore_orphan(0x1 | 0x2 | 0x40, true));
        assert!(should_ignore_orphan(0x1 | 0x40, true));
        assert!(!should_ignore_orphan(0x1 | 0x40, false));
    }
}
