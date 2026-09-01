//! Layout-neutral candidate verification for one read.
//!
//! Single-end and paired-end orchestration both own these same candidate,
//! verifier, metrics, and reusable-workspace primitives. Pair geometry and
//! mate rescue remain in the `paired_end` module.

use bsbit_core::alphabet::Base;
use bsbit_core::bisulfite::{AlignmentOrientation, BisulfiteStrand, strand_semantics};
use bsbit_index::reference::ReferenceIndex;

use crate::AlignmentError;
use crate::placement::ReadPlacement;
use crate::read_mapping_limits::{
    INITIAL_EDIT_DISTANCE, MAX_EDIT_DISTANCE, MAX_READ_BASES, VERIFICATION_BATCH,
};
use crate::search::combined_adaptive::{DIRECT_SINGLETON_PROOF, FLEXIBLE_NOMINAL_PROOF};
use crate::verification::ungapped::{bounded_complete_distance, reference_masks_by_query};
use crate::verification::{
    NarrowEndpointDistances, NarrowPlacementDistances, narrow_banded_fixed_start_batch,
    narrow_banded_placement_distances, narrow_banded_placement_distances_batch,
    narrow_banded_placement_distances_batch_d3, narrow_banded_placement_distances_batch_d5,
    narrow_banded_placement_distances_d3, narrow_banded_placement_distances_d5,
};

const MINIMUM_READ_BASES: usize = 3;
const EDIT_BUDGET: u64 = INITIAL_EDIT_DISTANCE as u64;
pub(crate) const LOCAL_FILTER_BLOCKS: usize = 8;
const LOCAL_FILTER_SUPPORT: usize = LOCAL_FILTER_BLOCKS - INITIAL_EDIT_DISTANCE as usize;
const VERIFICATION_PATTERN_BASES: usize = MAX_READ_BASES + 2 * MAX_EDIT_DISTANCE as usize;
const MAX_FLEXIBLE_PLACEMENTS: usize = (2 * MAX_EDIT_DISTANCE as usize + 1).pow(2);
// Collisions are ordinary misses. Generation tags invalidate the table
// between reads without clearing it on the hot path.
const VERIFICATION_CACHE_SLOTS: usize = 256;
const DIRECT_SINGLETON_DISTANCE_MASK: u8 = 0b11;

#[derive(Clone, Copy)]
struct LocalFilterBlock {
    query_start: u16,
    length: u8,
    query_masks: [u64; 5],
}

#[derive(Clone, Copy)]
pub(crate) struct LocalCandidateFilter {
    strand: BisulfiteStrand,
    blocks: [LocalFilterBlock; LOCAL_FILTER_BLOCKS],
}

impl LocalCandidateFilter {
    pub(crate) fn new(read: &[Base], strand: BisulfiteStrand) -> Self {
        let semantics = strand_semantics(strand);
        let blocks = core::array::from_fn(|ordinal| {
            let query_start = ordinal * read.len() / LOCAL_FILTER_BLOCKS;
            let query_end = (ordinal + 1) * read.len() / LOCAL_FILTER_BLOCKS;
            let mut query_masks = [0_u64; 5];
            for oriented_position in query_start..query_end {
                let query = match semantics.orientation() {
                    AlignmentOrientation::Forward => read[oriented_position],
                    AlignmentOrientation::Reverse => {
                        read[read.len() - oriented_position - 1].complement()
                    }
                };
                query_masks[usize::from(base_code(query))] |=
                    1_u64 << (oriented_position - query_start);
            }
            LocalFilterBlock {
                query_start: u16::try_from(query_start)
                    .expect("bounded local block start fits u16"),
                length: u8::try_from(query_end - query_start)
                    .expect("bounded local block length fits u8"),
                query_masks,
            }
        });
        Self { strand, blocks }
    }

    pub(crate) fn supports(self, reference: &ReferenceIndex, candidate: ReadCandidate) -> bool {
        let cytosine = strand_semantics(self.strand).cytosine_strand();
        let mut support = 0_usize;
        for (ordinal, block) in self.blocks.into_iter().enumerate() {
            let required = (1_u64 << block.length) - 1;
            let mut block_supported = false;
            for displacement in -i64::from(INITIAL_EDIT_DISTANCE)..=i64::from(INITIAL_EDIT_DISTANCE)
            {
                let Some(start) = candidate
                    .start()
                    .checked_add(u64::from(block.query_start))
                    .and_then(|start| start.checked_add_signed(displacement))
                    .and_then(|start| usize::try_from(start).ok())
                else {
                    continue;
                };
                let Some(reference_word) =
                    reference.packed_reference_word(candidate.contig_ordinal(), start)
                else {
                    continue;
                };
                let matched =
                    projected_word_match_mask(reference_word, &block.query_masks, cytosine);
                if matched & required == required {
                    block_supported = true;
                    break;
                }
            }
            support += usize::from(block_supported);
            if support + LOCAL_FILTER_BLOCKS - ordinal - 1 < LOCAL_FILTER_SUPPORT {
                return false;
            }
        }
        support >= LOCAL_FILTER_SUPPORT
    }
}

fn projected_word_match_mask(
    reference: bsbit_index::reference::PackedReferenceWord,
    query_masks: &[u64; 5],
    strand: bsbit_core::bisulfite::CytosineStrand,
) -> u64 {
    let low = reference.low();
    let high = reference.high();
    let equal = match strand {
        bsbit_core::bisulfite::CytosineStrand::Top => [!low & !high, low, !low & high, low, 0],
        bsbit_core::bisulfite::CytosineStrand::Bottom => [!low, low & !high, !low, low & high, 0],
    };
    equal
        .iter()
        .zip(query_masks)
        .fold(0_u64, |matched, (&reference, &query)| {
            matched | (reference & query)
        })
}

/// One integer-only whole-read candidate in forward contig coordinates.
#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct ReadCandidate {
    pub(crate) contig_ordinal: u64,
    pub(crate) start: u64,
    pub(crate) strand: BisulfiteStrand,
    pub(crate) proof_mask: u8,
}

#[derive(Clone, Copy, Default)]
pub(crate) struct VerificationCacheEntry {
    contig_ordinal: u64,
    start: u64,
    placement_start: usize,
    generation: u32,
    placement_len: u16,
    maximum_edit_distance: u8,
    mode: u8,
}

impl ReadCandidate {
    #[must_use]
    pub(crate) const fn contig_ordinal(self) -> u64 {
        self.contig_ordinal
    }

    #[must_use]
    pub(crate) const fn start(self) -> u64 {
        self.start
    }

    #[must_use]
    pub(crate) const fn strand(self) -> BisulfiteStrand {
        self.strand
    }
}

/// Worker-local, allocation-free verifier for one read.
///
/// The narrow-band kernel performs only five byte lanes of work per query
/// position at edit distance two and uses AVX2 when available. Query
/// orientation and bisulfite masks are prepared once per strand, not once per
/// candidate.
pub(crate) struct PlacementVerifier {
    read_len: usize,
    query_codes: [[u8; MAX_READ_BASES]; 2],
    reference_masks: [[u8; 5]; 2],
    pattern_codes: [u8; VERIFICATION_PATTERN_BASES * VERIFICATION_BATCH],
    endpoint_distances: [NarrowEndpointDistances; VERIFICATION_BATCH],
    placement_distances: [NarrowPlacementDistances; 4],
}

impl PlacementVerifier {
    pub(crate) fn new(read: &[Base]) -> Result<Self, AlignmentError> {
        if !(MINIMUM_READ_BASES..=MAX_READ_BASES).contains(&read.len()) {
            return Err(AlignmentError::UnsupportedReadLength { length: read.len() });
        }
        let mut query_codes = [[4_u8; MAX_READ_BASES]; 2];
        for (position, &base) in read.iter().enumerate() {
            query_codes[0][position] = base_code(base);
            query_codes[1][read.len() - position - 1] = base_code(base.complement());
        }
        Ok(Self {
            read_len: read.len(),
            query_codes,
            reference_masks: [
                reference_masks(BisulfiteStrand::OT),
                reference_masks(BisulfiteStrand::CTOB),
            ],
            pattern_codes: [4_u8; VERIFICATION_PATTERN_BASES * VERIFICATION_BATCH],
            endpoint_distances: [NarrowEndpointDistances::EMPTY; VERIFICATION_BATCH],
            placement_distances: [NarrowPlacementDistances::EMPTY; 4],
        })
    }

    /// Returns the best score and all tied in-budget reference lengths for one
    /// fixed candidate start. Bit `d` represents length `read_len - 2 + d`.
    #[cfg(test)]
    pub(crate) fn verify(
        &mut self,
        reference: &ReferenceIndex,
        candidate: ReadCandidate,
    ) -> Result<Option<PlacementVerification>, AlignmentError> {
        let mut output = [None];
        self.verify_batch(reference, core::slice::from_ref(&candidate), &mut output)?;
        Ok(output[0])
    }

    fn verify_flexible_nominal_batch(
        &mut self,
        reference: &ReferenceIndex,
        candidates: &[ReadCandidate],
        output: &mut [[Option<FlexibleVerification>; MAX_FLEXIBLE_PLACEMENTS]; 4],
        retained: &mut [usize; 4],
    ) -> Result<(), AlignmentError> {
        debug_assert!(!candidates.is_empty() && candidates.len() <= 4);
        // `retained` is the authoritative initialized prefix for every row.
        // A lane with no hit remains empty, so stale slots need no full-row
        // clearing between candidate batches.
        retained.fill(0);
        let first = candidates[0];
        let contig = reference.contig_by_ordinal(first.contig_ordinal()).ok_or(
            AlignmentError::InvalidContigOrdinal {
                ordinal: first.contig_ordinal(),
            },
        )?;
        let budget = usize::from(INITIAL_EDIT_DISTANCE);
        let pattern_len = self.read_len + 2 * budget;
        let mut window_starts = [0_usize; 4];
        for (ordinal, &candidate) in candidates.iter().enumerate() {
            debug_assert_eq!(candidate.contig_ordinal(), first.contig_ordinal());
            debug_assert_eq!(candidate.strand(), first.strand());
            let nominal = usize::try_from(candidate.start()).map_err(|_| {
                AlignmentError::CandidateCoordinateOverflow {
                    start: candidate.start(),
                }
            })?;
            let window_start = nominal.saturating_sub(budget);
            window_starts[ordinal] = window_start;
            let pattern = &mut self.pattern_codes
                [ordinal * pattern_len..ordinal.saturating_add(1) * pattern_len];
            pattern.fill(4);
            let available = contig.sequence().bases().len().saturating_sub(window_start);
            let copied = available.min(pattern_len);
            for (destination, &base) in pattern[..copied]
                .iter_mut()
                .zip(&contig.sequence().bases()[window_start..window_start + copied])
            {
                *destination = base_code(base);
            }
        }
        let semantics = strand_semantics(first.strand());
        let axis = usize::from(matches!(
            semantics.orientation(),
            AlignmentOrientation::Reverse
        ));
        let cytosine_axis = usize::from(matches!(
            semantics.cytosine_strand(),
            bsbit_core::bisulfite::CytosineStrand::Bottom
        ));
        if candidates.len() == 1 {
            debug_assert_eq!(budget, usize::from(INITIAL_EDIT_DISTANCE));
            self.placement_distances[0] = narrow_banded_placement_distances_d3(
                &self.reference_masks[cytosine_axis],
                &self.query_codes[axis][..self.read_len],
                &self.pattern_codes[..pattern_len],
            )?;
        } else {
            debug_assert_eq!(budget, usize::from(INITIAL_EDIT_DISTANCE));
            narrow_banded_placement_distances_batch_d3(
                &self.reference_masks[cytosine_axis],
                &self.query_codes[axis][..self.read_len],
                &self.pattern_codes[..pattern_len * candidates.len()],
                &mut self.placement_distances[..candidates.len()],
            )?;
        }
        for ordinal in 0..candidates.len() {
            let mut best = u32::MAX;
            for start_delta in 0..=2 * budget {
                for endpoint_delta in 0..=2 * budget {
                    let Some(distance) =
                        self.placement_distances[ordinal].distance(start_delta, endpoint_delta)
                    else {
                        continue;
                    };
                    if window_starts[ordinal] + self.read_len + endpoint_delta
                        > contig.sequence().bases().len()
                        || self.read_len + endpoint_delta < start_delta
                    {
                        continue;
                    }
                    if distance < best {
                        best = distance;
                        // `retained` is the authoritative initialized prefix;
                        // resetting its length makes the previous best row
                        // unreachable without clearing the whole tie buffer.
                        retained[ordinal] = 0;
                    }
                    if distance == best {
                        debug_assert!(retained[ordinal] < MAX_FLEXIBLE_PLACEMENTS);
                        output[ordinal][retained[ordinal]] = Some(FlexibleVerification {
                            start: u64::try_from(window_starts[ordinal] + start_delta)
                                .expect("bounded flexible start fits u64"),
                            end: u64::try_from(
                                window_starts[ordinal] + self.read_len + endpoint_delta,
                            )
                            .expect("bounded flexible end fits u64"),
                            distance: u8::try_from(distance).expect("in-budget distance fits u8"),
                        });
                        retained[ordinal] += 1;
                    }
                }
            }
        }
        Ok(())
    }

    // Pattern materialization, SIMD dispatch, and tie extraction share the
    // same fixed buffers and must remain one ordered verification pass.
    #[allow(clippy::too_many_lines)]
    fn verify_batch(
        &mut self,
        reference: &ReferenceIndex,
        candidates: &[ReadCandidate],
        output: &mut [Option<PlacementVerification>],
    ) -> Result<(), AlignmentError> {
        if candidates.is_empty() || candidates.len() > VERIFICATION_BATCH {
            return Err(AlignmentError::VerificationBatchSize {
                observed: candidates.len(),
            });
        }
        if output.len() != candidates.len() {
            return Err(AlignmentError::VerificationOutputSize {
                candidates: candidates.len(),
                output: output.len(),
            });
        }
        let strand = candidates[0].strand();
        if candidates
            .iter()
            .any(|candidate| candidate.strand() != strand)
        {
            return Err(AlignmentError::MixedVerificationStrands);
        }
        output.fill(None);
        let budget = usize::from(INITIAL_EDIT_DISTANCE);
        let pattern_len = self.read_len + 2 * budget;
        let mut available = [0_usize; VERIFICATION_BATCH];
        let mut starts = [0_usize; VERIFICATION_BATCH];
        for (ordinal, &candidate) in candidates.iter().enumerate() {
            let contig = reference
                .contig_by_ordinal(candidate.contig_ordinal())
                .ok_or(AlignmentError::InvalidContigOrdinal {
                    ordinal: candidate.contig_ordinal(),
                })?;
            let start = usize::try_from(candidate.start()).map_err(|_| {
                AlignmentError::CandidateCoordinateOverflow {
                    start: candidate.start(),
                }
            })?;
            starts[ordinal] = start;
            if start >= contig.sequence().bases().len() {
                continue;
            }
            available[ordinal] = contig.sequence().bases().len() - start;
        }

        let orientation = strand_semantics(strand).orientation();
        let axis = usize::from(matches!(orientation, AlignmentOrientation::Reverse));
        let cytosine_axis = usize::from(matches!(
            strand_semantics(strand).cytosine_strand(),
            bsbit_core::bisulfite::CytosineStrand::Bottom
        ));
        for (ordinal, (&candidate, &start)) in candidates.iter().zip(&starts).enumerate() {
            let contig = reference
                .contig_by_ordinal(candidate.contig_ordinal())
                .ok_or(AlignmentError::InvalidContigOrdinal {
                    ordinal: candidate.contig_ordinal(),
                })?;
            let pattern = &mut self.pattern_codes
                [ordinal * pattern_len..ordinal.saturating_add(1) * pattern_len];
            let copied = available[ordinal].min(self.read_len + budget);
            pattern.fill(4);
            if copied != 0 {
                for (destination, &base) in pattern[budget..budget + copied]
                    .iter_mut()
                    .zip(&contig.sequence().bases()[start..start + copied])
                {
                    *destination = base_code(base);
                }
            }
        }
        narrow_banded_fixed_start_batch(
            &self.reference_masks[cytosine_axis],
            &self.query_codes[axis][..self.read_len],
            &self.pattern_codes[..pattern_len * candidates.len()],
            budget,
            &mut self.endpoint_distances[..candidates.len()],
        )?;
        for (ordinal, result) in output.iter_mut().enumerate() {
            let mut best = u32::MAX;
            let mut tied_lengths = 0_u8;
            for delta in 0..=2 * budget {
                let reference_len = self.read_len - budget + delta;
                if reference_len > available[ordinal] {
                    continue;
                }
                let Some(distance) = self.endpoint_distances[ordinal].distance(delta) else {
                    continue;
                };
                if distance < best {
                    best = distance;
                    tied_lengths = 1_u8 << delta;
                } else if distance == best {
                    tied_lengths |= 1_u8 << delta;
                }
            }
            *result = (best <= u32::from(INITIAL_EDIT_DISTANCE)).then(|| PlacementVerification {
                distance: u8::try_from(best).expect("in-budget distance fits u8"),
                tied_lengths,
            });
        }
        Ok(())
    }
}

/// Flexible-only verifier used by the combined paired-end backend.
///
/// Keeping this state in a separate non-inlined call avoids reserving and
/// initializing the unused 32-candidate fixed-start buffers on workloads whose
/// candidates all carry `FLEXIBLE_NOMINAL_PROOF`.
struct FlexibleVerifier {
    read_len: usize,
    query_codes: [[u8; MAX_READ_BASES]; 2],
    reference_masks: [[u8; 5]; 2],
    pattern_codes: [u8; VERIFICATION_PATTERN_BASES * 4],
    placement_distances: [NarrowPlacementDistances; 4],
}

impl FlexibleVerifier {
    fn new(read: &[Base]) -> Result<Self, AlignmentError> {
        if !(MINIMUM_READ_BASES..=MAX_READ_BASES).contains(&read.len()) {
            return Err(AlignmentError::UnsupportedReadLength { length: read.len() });
        }
        let mut query_codes = [[4_u8; MAX_READ_BASES]; 2];
        for (position, &base) in read.iter().enumerate() {
            query_codes[0][position] = base_code(base);
            query_codes[1][read.len() - position - 1] = base_code(base.complement());
        }
        Ok(Self {
            read_len: read.len(),
            query_codes,
            reference_masks: [
                reference_masks(BisulfiteStrand::OT),
                reference_masks(BisulfiteStrand::CTOB),
            ],
            pattern_codes: [4_u8; VERIFICATION_PATTERN_BASES * 4],
            placement_distances: [NarrowPlacementDistances::EMPTY; 4],
        })
    }

    // Flexible start/end enumeration, SIMD dispatch, and tie retention share
    // the same fixed buffers and must remain one ordered verification pass.
    #[allow(clippy::too_many_lines)]
    fn verify_batch(
        &mut self,
        reference: &ReferenceIndex,
        candidates: &[ReadCandidate],
        maximum_edit_distance: u8,
        output: &mut [[Option<FlexibleVerification>; MAX_FLEXIBLE_PLACEMENTS]; 4],
        retained: &mut [usize; 4],
    ) -> Result<(), AlignmentError> {
        let budget = usize::from(maximum_edit_distance);
        let maximum_batch = 32 / (2 * budget + 1);
        debug_assert!(!candidates.is_empty() && candidates.len() <= maximum_batch);
        // `retained` is the authoritative initialized prefix for every row.
        // A lane with no hit remains empty, so stale slots need no full-row
        // clearing between candidate batches.
        retained.fill(0);
        let first = candidates[0];
        let contig = reference.contig_by_ordinal(first.contig_ordinal()).ok_or(
            AlignmentError::InvalidContigOrdinal {
                ordinal: first.contig_ordinal(),
            },
        )?;
        let pattern_len = self.read_len + 2 * budget;
        let mut window_starts = [0_usize; 4];
        for (ordinal, &candidate) in candidates.iter().enumerate() {
            debug_assert_eq!(candidate.contig_ordinal(), first.contig_ordinal());
            debug_assert_eq!(candidate.strand(), first.strand());
            let nominal = usize::try_from(candidate.start()).map_err(|_| {
                AlignmentError::CandidateCoordinateOverflow {
                    start: candidate.start(),
                }
            })?;
            let window_start = nominal.saturating_sub(budget);
            window_starts[ordinal] = window_start;
            let pattern =
                &mut self.pattern_codes[ordinal * pattern_len..(ordinal + 1) * pattern_len];
            pattern.fill(4);
            let available = contig.sequence().bases().len().saturating_sub(window_start);
            let copied = available.min(pattern_len);
            for (destination, &base) in pattern[..copied]
                .iter_mut()
                .zip(&contig.sequence().bases()[window_start..window_start + copied])
            {
                *destination = base_code(base);
            }
        }
        let semantics = strand_semantics(first.strand());
        let axis = usize::from(matches!(
            semantics.orientation(),
            AlignmentOrientation::Reverse
        ));
        let cytosine_axis = usize::from(matches!(
            semantics.cytosine_strand(),
            bsbit_core::bisulfite::CytosineStrand::Bottom
        ));
        if candidates.len() == 1 {
            self.placement_distances[0] = if maximum_edit_distance == INITIAL_EDIT_DISTANCE {
                narrow_banded_placement_distances_d3(
                    &self.reference_masks[cytosine_axis],
                    &self.query_codes[axis][..self.read_len],
                    &self.pattern_codes[..pattern_len],
                )?
            } else if maximum_edit_distance == MAX_EDIT_DISTANCE {
                narrow_banded_placement_distances_d5(
                    &self.reference_masks[cytosine_axis],
                    &self.query_codes[axis][..self.read_len],
                    &self.pattern_codes[..pattern_len],
                )?
            } else {
                narrow_banded_placement_distances(
                    &self.reference_masks[cytosine_axis],
                    &self.query_codes[axis][..self.read_len],
                    &self.pattern_codes[..pattern_len],
                    budget,
                )?
            };
        } else if maximum_edit_distance == INITIAL_EDIT_DISTANCE {
            narrow_banded_placement_distances_batch_d3(
                &self.reference_masks[cytosine_axis],
                &self.query_codes[axis][..self.read_len],
                &self.pattern_codes[..pattern_len * candidates.len()],
                &mut self.placement_distances[..candidates.len()],
            )?;
        } else if maximum_edit_distance == MAX_EDIT_DISTANCE {
            narrow_banded_placement_distances_batch_d5(
                &self.reference_masks[cytosine_axis],
                &self.query_codes[axis][..self.read_len],
                &self.pattern_codes[..pattern_len * candidates.len()],
                &mut self.placement_distances[..candidates.len()],
            )?;
        } else {
            narrow_banded_placement_distances_batch(
                &self.reference_masks[cytosine_axis],
                &self.query_codes[axis][..self.read_len],
                &self.pattern_codes[..pattern_len * candidates.len()],
                budget,
                &mut self.placement_distances[..candidates.len()],
            )?;
        }
        for ordinal in 0..candidates.len() {
            let mut best = u32::MAX;
            for start_delta in 0..=2 * budget {
                for endpoint_delta in 0..=2 * budget {
                    let Some(distance) =
                        self.placement_distances[ordinal].distance(start_delta, endpoint_delta)
                    else {
                        continue;
                    };
                    if window_starts[ordinal] + self.read_len + endpoint_delta
                        > contig.sequence().bases().len()
                        || self.read_len + endpoint_delta < start_delta
                    {
                        continue;
                    }
                    if distance < best {
                        best = distance;
                        // Callers consume only `[..retained[ordinal]]`, so a
                        // new best invalidates the old row by resetting its
                        // initialized prefix rather than clearing every slot.
                        retained[ordinal] = 0;
                    }
                    if distance == best {
                        debug_assert!(retained[ordinal] < MAX_FLEXIBLE_PLACEMENTS);
                        output[ordinal][retained[ordinal]] = Some(FlexibleVerification {
                            start: u64::try_from(window_starts[ordinal] + start_delta)
                                .expect("bounded flexible start fits u64"),
                            end: u64::try_from(
                                window_starts[ordinal] + self.read_len + endpoint_delta,
                            )
                            .expect("bounded flexible end fits u64"),
                            distance: u8::try_from(distance).expect("in-budget distance fits u8"),
                        });
                        retained[ordinal] += 1;
                    }
                }
            }
        }
        Ok(())
    }
}

#[derive(Clone, Copy)]
struct FlexibleVerification {
    start: u64,
    end: u64,
    distance: u8,
}

/// Complete best endpoint frontier for one fixed candidate start.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct PlacementVerification {
    distance: u8,
    tied_lengths: u8,
}

impl PlacementVerification {
    #[must_use]
    pub(crate) const fn distance(self) -> u8 {
        self.distance
    }

    #[must_use]
    pub(crate) const fn tied_lengths(self) -> u8 {
        self.tied_lengths
    }
}

/// Observable work counters for the lean candidate path.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub(crate) struct ReadAlignmentMetrics {
    pub(crate) located_rows: u64,
    pub(crate) emitted_candidate_starts: u64,
    pub(crate) distinct_candidate_starts: u64,
    pub(crate) verified_placements: u64,
}

/// Reusable storage owned by one mapping worker.
pub(crate) struct ReadWorkspace {
    pub(crate) candidate_nominals: Vec<ReadCandidate>,
    pub(crate) candidates: Vec<ReadCandidate>,
    pub(crate) placements: Vec<ReadPlacement>,
    pub(crate) verification_cache: Vec<VerificationCacheEntry>,
    pub(crate) verification_cache_placements: Vec<ReadPlacement>,
    pub(crate) verification_cache_generation: u32,
    pub(crate) verification_cache_population: usize,
}

impl ReadWorkspace {
    #[must_use]
    pub(crate) fn with_capacity(candidate_capacity: usize, placement_capacity: usize) -> Self {
        Self {
            candidate_nominals: Vec::with_capacity(candidate_capacity / 5 + 1),
            candidates: Vec::with_capacity(candidate_capacity),
            placements: Vec::with_capacity(placement_capacity),
            verification_cache: vec![VerificationCacheEntry::default(); VERIFICATION_CACHE_SLOTS],
            verification_cache_placements: Vec::with_capacity(placement_capacity),
            verification_cache_generation: 0,
            verification_cache_population: 0,
        }
    }

    pub(crate) fn begin_verification_cache_read(&mut self) {
        self.verification_cache_generation = self.verification_cache_generation.wrapping_add(1);
        if self.verification_cache_generation == 0 {
            self.verification_cache
                .fill(VerificationCacheEntry::default());
            self.verification_cache_generation = 1;
        }
        self.verification_cache_population = 0;
        self.verification_cache_placements.clear();
    }

    #[inline]
    fn verification_cache_slot(
        candidate: ReadCandidate,
        maximum_edit_distance: u8,
        flexible: bool,
    ) -> usize {
        let mut mixed = candidate.start()
            ^ candidate.contig_ordinal().rotate_left(23)
            ^ (u64::from(candidate.strand() as u8) << 57)
            ^ (u64::from(maximum_edit_distance) << 49)
            ^ (u64::from(flexible) << 48);
        mixed ^= mixed >> 33;
        mixed = mixed.wrapping_mul(0xff51_afd7_ed55_8ccd);
        mixed ^= mixed >> 33;
        usize::try_from(mixed & (VERIFICATION_CACHE_SLOTS as u64 - 1))
            .expect("verification cache slot fits usize")
    }

    #[inline]
    fn verification_cache_mode(candidate: ReadCandidate, flexible: bool) -> u8 {
        candidate.strand() as u8 | (u8::from(flexible) << 2)
    }

    #[inline]
    fn cached_verification_range(
        &self,
        candidate: ReadCandidate,
        maximum_edit_distance: u8,
        flexible: bool,
    ) -> Option<(usize, usize)> {
        let entry = self.verification_cache
            [Self::verification_cache_slot(candidate, maximum_edit_distance, flexible)];
        (entry.generation == self.verification_cache_generation
            && entry.contig_ordinal == candidate.contig_ordinal()
            && entry.start == candidate.start()
            && entry.maximum_edit_distance == maximum_edit_distance
            && entry.mode == Self::verification_cache_mode(candidate, flexible))
        .then_some((
            entry.placement_start,
            entry.placement_start + usize::from(entry.placement_len),
        ))
    }

    pub(crate) fn retain_uncached_candidates(
        &mut self,
        maximum_edit_distance: u8,
        all_flexible: bool,
    ) -> (u64, u64) {
        if self.verification_cache_population == 0 {
            return (0, 0);
        }
        let mut retained = 0_usize;
        let mut hits = 0_u64;
        let mut misses = 0_u64;
        for index in 0..self.candidates.len() {
            let candidate = self.candidates[index];
            let flexible = candidate.proof_mask & FLEXIBLE_NOMINAL_PROOF != 0;
            let effective_budget = if flexible && all_flexible {
                maximum_edit_distance
            } else {
                INITIAL_EDIT_DISTANCE
            };
            if let Some((start, end)) =
                self.cached_verification_range(candidate, effective_budget, flexible)
            {
                hits = hits.saturating_add(1);
                self.placements
                    .extend_from_slice(&self.verification_cache_placements[start..end]);
            } else {
                misses = misses.saturating_add(1);
                self.candidates[retained] = candidate;
                retained += 1;
            }
        }
        self.candidates.truncate(retained);
        (hits, misses)
    }

    // Cache storage is split into parallel worker-owned buffers; keeping the
    // insertion seam explicit prevents hidden allocation or aliasing.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn cache_candidate_verification(
        verification_cache: &mut [VerificationCacheEntry],
        verification_cache_placements: &mut Vec<ReadPlacement>,
        verification_cache_generation: u32,
        verification_cache_population: &mut usize,
        placements: &[ReadPlacement],
        candidate: ReadCandidate,
        maximum_edit_distance: u8,
        flexible: bool,
        placement_start: usize,
    ) {
        let cached_placement_start = verification_cache_placements.len();
        verification_cache_placements.extend_from_slice(&placements[placement_start..]);
        let slot = Self::verification_cache_slot(candidate, maximum_edit_distance, flexible);
        if verification_cache[slot].generation != verification_cache_generation {
            *verification_cache_population += 1;
        }
        verification_cache[slot] = VerificationCacheEntry {
            contig_ordinal: candidate.contig_ordinal(),
            start: candidate.start(),
            placement_start: cached_placement_start,
            generation: verification_cache_generation,
            placement_len: u16::try_from(placements.len() - placement_start)
                .expect("bounded candidate verification placements fit u16"),
            maximum_edit_distance,
            mode: Self::verification_cache_mode(candidate, flexible),
        };
    }
}

impl ReadWorkspace {
    pub(crate) fn verify_candidates_with_budget(
        &mut self,
        reference: &ReferenceIndex,
        read: &[Base],
        metrics: ReadAlignmentMetrics,
        maximum_edit_distance: u8,
    ) -> Result<(&[ReadPlacement], ReadAlignmentMetrics), AlignmentError> {
        self.verify_candidates_with_budget_and_order(
            reference,
            read,
            metrics,
            maximum_edit_distance,
            false,
        )
    }

    pub(crate) fn verify_sorted_candidates_with_budget(
        &mut self,
        reference: &ReferenceIndex,
        read: &[Base],
        metrics: ReadAlignmentMetrics,
        maximum_edit_distance: u8,
    ) -> Result<(&[ReadPlacement], ReadAlignmentMetrics), AlignmentError> {
        self.verify_candidates_with_budget_and_order(
            reference,
            read,
            metrics,
            maximum_edit_distance,
            true,
        )
    }

    // Candidate deduplication, cache lookup, and verification form one ordered
    // pass over the same reusable workspace buffers.
    #[allow(clippy::too_many_lines)]
    fn verify_candidates_with_budget_and_order(
        &mut self,
        reference: &ReferenceIndex,
        read: &[Base],
        mut metrics: ReadAlignmentMetrics,
        maximum_edit_distance: u8,
        nominal_candidates_sorted: bool,
    ) -> Result<(&[ReadPlacement], ReadAlignmentMetrics), AlignmentError> {
        if !nominal_candidates_sorted {
            sort_nominal_candidates(&mut self.candidate_nominals);
        }
        let mut nominal_retained = 0_usize;
        for index in 0..self.candidate_nominals.len() {
            let incoming = self.candidate_nominals[index];
            if nominal_retained != 0
                && candidate_key(self.candidate_nominals[nominal_retained - 1])
                    == candidate_key(incoming)
            {
                self.candidate_nominals[nominal_retained - 1].proof_mask |= incoming.proof_mask;
            } else {
                self.candidate_nominals[nominal_retained] = incoming;
                nominal_retained += 1;
            }
        }
        self.candidate_nominals.truncate(nominal_retained);

        self.candidates.clear();
        let mut index = 0_usize;
        while index < self.candidate_nominals.len() {
            let first = self.candidate_nominals[index];
            if first.proof_mask & DIRECT_SINGLETON_PROOF != 0 {
                self.placements.push(ReadPlacement::strict(
                    first.contig_ordinal(),
                    first.start(),
                    first.start().saturating_add(
                        u64::try_from(read.len()).expect("bounded paired-end read length fits u64"),
                    ),
                    first.strand(),
                    first.proof_mask & DIRECT_SINGLETON_DISTANCE_MASK,
                ));
                index += 1;
                continue;
            }
            if first.proof_mask & FLEXIBLE_NOMINAL_PROOF != 0 {
                self.candidates.push(first);
                index += 1;
                continue;
            }
            let mut range_end = first.start().saturating_add(EDIT_BUDGET);
            let range_start = first.start().saturating_sub(EDIT_BUDGET);
            index += 1;
            while index < self.candidate_nominals.len() {
                let incoming = self.candidate_nominals[index];
                if incoming.strand() != first.strand()
                    || incoming.contig_ordinal() != first.contig_ordinal()
                    || incoming.proof_mask & DIRECT_SINGLETON_PROOF != 0
                    || incoming.proof_mask & FLEXIBLE_NOMINAL_PROOF != 0
                    || incoming.start().saturating_sub(EDIT_BUDGET) > range_end.saturating_add(1)
                {
                    break;
                }
                range_end = range_end.max(incoming.start().saturating_add(EDIT_BUDGET));
                index += 1;
            }
            for start in range_start..=range_end {
                self.candidates.push(ReadCandidate {
                    contig_ordinal: first.contig_ordinal(),
                    start,
                    strand: first.strand(),
                    proof_mask: 1,
                });
            }
        }
        metrics.emitted_candidate_starts = u64::try_from(self.candidates.len()).unwrap_or(u64::MAX);
        if self
            .candidates
            .iter()
            .any(|candidate| candidate.proof_mask & FLEXIBLE_NOMINAL_PROOF == 0)
        {
            let local_filters =
                BisulfiteStrand::ALL.map(|strand| LocalCandidateFilter::new(read, strand));
            let mut locally_retained = 0_usize;
            for index in 0..self.candidates.len() {
                let candidate = self.candidates[index];
                if candidate.proof_mask & FLEXIBLE_NOMINAL_PROOF != 0
                    || local_filters[strand_index(candidate.strand())]
                        .supports(reference, candidate)
                {
                    self.candidates[locally_retained] = candidate;
                    locally_retained += 1;
                }
            }
            self.candidates.truncate(locally_retained);
        }
        metrics.distinct_candidate_starts =
            u64::try_from(self.candidates.len()).unwrap_or(u64::MAX);
        let all_flexible = self
            .candidates
            .iter()
            .all(|candidate| candidate.proof_mask & FLEXIBLE_NOMINAL_PROOF != 0);
        let _ = self.retain_uncached_candidates(maximum_edit_distance, all_flexible);
        if self.candidates.is_empty() {
            // Direct singleton placements, or no evidence at all, require no
            // verifier state. In particular, do not construct the large
            // general verifier merely to iterate an empty candidate slice.
        } else if all_flexible {
            self.verify_flexible_candidates(reference, read, maximum_edit_distance)?;
        } else {
            self.verify_general_candidates(reference, read)?;
        }
        self.placements.sort_unstable_by_key(|placement| {
            (
                placement.contig_ordinal,
                placement.strand,
                placement.start,
                placement.end,
                placement.distance,
            )
        });
        self.placements.dedup();
        metrics.verified_placements = u64::try_from(self.placements.len()).unwrap_or(u64::MAX);
        Ok((&self.placements, metrics))
    }

    #[inline(never)]
    fn verify_flexible_candidates(
        &mut self,
        reference: &ReferenceIndex,
        read: &[Base],
        maximum_edit_distance: u8,
    ) -> Result<(), AlignmentError> {
        let mut verifier = FlexibleVerifier::new(read)?;
        let batch_size = (32 / (2 * usize::from(maximum_edit_distance) + 1)).min(4);
        let mut run_start = 0_usize;
        while run_start < self.candidates.len() {
            let strand = self.candidates[run_start].strand();
            let contig = self.candidates[run_start].contig_ordinal();
            let run_end = self.candidates[run_start..].partition_point(|candidate| {
                candidate.strand() == strand && candidate.contig_ordinal() == contig
            }) + run_start;
            let mut flexible = [[None; MAX_FLEXIBLE_PLACEMENTS]; 4];
            let mut retained = [0_usize; 4];
            for batch in self.candidates[run_start..run_end].chunks(batch_size) {
                verifier.verify_batch(
                    reference,
                    batch,
                    maximum_edit_distance,
                    &mut flexible,
                    &mut retained,
                )?;
                for (ordinal, &candidate) in batch.iter().enumerate() {
                    let placement_start = self.placements.len();
                    for result in flexible[ordinal][..retained[ordinal]].iter().flatten() {
                        self.placements.push(ReadPlacement::strict(
                            candidate.contig_ordinal(),
                            result.start,
                            result.end,
                            candidate.strand(),
                            result.distance,
                        ));
                    }
                    Self::cache_candidate_verification(
                        &mut self.verification_cache,
                        &mut self.verification_cache_placements,
                        self.verification_cache_generation,
                        &mut self.verification_cache_population,
                        &self.placements,
                        candidate,
                        maximum_edit_distance,
                        true,
                        placement_start,
                    );
                }
            }
            run_start = run_end;
        }
        Ok(())
    }

    #[inline(never)]
    fn verify_general_candidates(
        &mut self,
        reference: &ReferenceIndex,
        read: &[Base],
    ) -> Result<(), AlignmentError> {
        let mut verifier = PlacementVerifier::new(read)?;
        let minimum_length = read.len() - usize::from(INITIAL_EDIT_DISTANCE);
        let mut run_start = 0_usize;
        while run_start < self.candidates.len() {
            let strand = self.candidates[run_start].strand();
            let contig = self.candidates[run_start].contig_ordinal();
            let run_end = self.candidates[run_start..].partition_point(|candidate| {
                candidate.strand() == strand && candidate.contig_ordinal() == contig
            }) + run_start;
            if self.candidates[run_start].proof_mask & FLEXIBLE_NOMINAL_PROOF != 0 {
                let mut flexible = [[None; MAX_FLEXIBLE_PLACEMENTS]; 4];
                let mut retained = [0_usize; 4];
                for batch in self.candidates[run_start..run_end].chunks(4) {
                    verifier.verify_flexible_nominal_batch(
                        reference,
                        batch,
                        &mut flexible,
                        &mut retained,
                    )?;
                    for (ordinal, &candidate) in batch.iter().enumerate() {
                        let placement_start = self.placements.len();
                        for result in flexible[ordinal][..retained[ordinal]].iter().flatten() {
                            self.placements.push(ReadPlacement::strict(
                                candidate.contig_ordinal(),
                                result.start,
                                result.end,
                                candidate.strand(),
                                result.distance,
                            ));
                        }
                        Self::cache_candidate_verification(
                            &mut self.verification_cache,
                            &mut self.verification_cache_placements,
                            self.verification_cache_generation,
                            &mut self.verification_cache_population,
                            &self.placements,
                            candidate,
                            INITIAL_EDIT_DISTANCE,
                            true,
                            placement_start,
                        );
                    }
                }
                run_start = run_end;
                continue;
            }
            for batch in self.candidates[run_start..run_end].chunks(VERIFICATION_BATCH) {
                let mut batch_results = [None; VERIFICATION_BATCH];
                verifier.verify_batch(reference, batch, &mut batch_results[..batch.len()])?;
                for (&candidate, result) in batch.iter().zip(batch_results) {
                    let placement_start = self.placements.len();
                    if let Some(result) = result {
                        let delta = canonical_tied_delta(result.tied_lengths());
                        let reference_len = minimum_length + delta;
                        let end = candidate
                            .start()
                            .checked_add(
                                u64::try_from(reference_len).expect("bounded length fits u64"),
                            )
                            .ok_or(AlignmentError::CandidateEndpointOverflow {
                                start: candidate.start(),
                                length: reference_len,
                            })?;
                        self.placements.push(ReadPlacement::strict(
                            candidate.contig_ordinal(),
                            candidate.start(),
                            end,
                            candidate.strand(),
                            result.distance(),
                        ));
                    }
                    Self::cache_candidate_verification(
                        &mut self.verification_cache,
                        &mut self.verification_cache_placements,
                        self.verification_cache_generation,
                        &mut self.verification_cache_population,
                        &self.placements,
                        candidate,
                        INITIAL_EDIT_DISTANCE,
                        false,
                        placement_start,
                    );
                }
            }
            run_start = run_end;
        }
        Ok(())
    }
}

const fn candidate_key(candidate: ReadCandidate) -> (BisulfiteStrand, u64, u64) {
    (candidate.strand, candidate.contig_ordinal, candidate.start)
}

fn canonical_tied_delta(tied: u8) -> usize {
    let preferred = usize::from(INITIAL_EDIT_DISTANCE);
    if tied & (1_u8 << preferred) != 0 {
        return preferred;
    }
    for distance in 1..=usize::from(INITIAL_EDIT_DISTANCE) {
        if preferred >= distance && tied & (1_u8 << (preferred - distance)) != 0 {
            return preferred - distance;
        }
        let longer = preferred + distance;
        if tied & (1_u8 << longer) != 0 {
            return longer;
        }
    }
    tied.trailing_zeros() as usize
}

pub(crate) fn sort_nominal_candidates(candidates: &mut [ReadCandidate]) {
    candidates.sort_unstable_by_key(|candidate| {
        (
            candidate.strand(),
            candidate.contig_ordinal(),
            candidate.start(),
        )
    });
}

pub(crate) fn ungapped_distance(
    reference: &ReferenceIndex,
    read: &[Base],
    candidate: ReadCandidate,
) -> Option<u8> {
    let contig = reference.contig_by_ordinal(candidate.contig_ordinal())?;
    let start = usize::try_from(candidate.start()).ok()?;
    bounded_complete_distance(
        contig.sequence().bases(),
        start,
        read,
        candidate.strand(),
        INITIAL_EDIT_DISTANCE,
    )
}

const fn base_code(base: Base) -> u8 {
    base.storage_code()
}

pub(crate) const fn strand_index(strand: BisulfiteStrand) -> usize {
    match strand {
        BisulfiteStrand::OT => 0,
        BisulfiteStrand::OB => 1,
        BisulfiteStrand::CTOT => 2,
        BisulfiteStrand::CTOB => 3,
    }
}

fn reference_masks(strand: BisulfiteStrand) -> [u8; 5] {
    reference_masks_by_query(strand_semantics(strand).cytosine_strand())
}
