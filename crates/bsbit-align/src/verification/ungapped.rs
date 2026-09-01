//! Allocation-free ungapped distance and bounded semi-global endpoint search.

use bsbit_core::alphabet::Base;
use bsbit_core::bisulfite::{
    AlignmentOrientation, BisulfiteStrand, CytosineStrand, classify_bases, strand_semantics,
};

/// Largest query represented by the allocation-free ungapped profile.
pub const MAX_UNGAPPED_QUERY_BASES: usize = 192;
type UngappedSelection = ((u8, usize, u8, usize, usize), UngappedEndpoint);

/// Conversion-aware mismatch and barrier prefixes for one selected reference
/// origin and query.
///
/// This value deliberately knows nothing about reference indexes, candidates,
/// mates, or mapping confidence. Callers select a reference origin first and
/// then use the profile to compare possible ungapped endpoints.
#[derive(Clone, Debug)]
pub struct UngappedProfile {
    read_length: usize,
    nominal_start: usize,
    orientation: AlignmentOrientation,
    complete_reference: bool,
    mismatch_prefix: [u16; MAX_UNGAPPED_QUERY_BASES + 1],
    barrier_prefix: [u16; MAX_UNGAPPED_QUERY_BASES + 1],
}

/// Returns the conversion-aware whole-query mismatch distance when a complete
/// ungapped reference span is present and the distance stays within the bound.
///
/// Unlike [`UngappedProfile`], this direct scan does not build prefix arrays;
/// callers that only need the complete distance can also stop at the first
/// mismatch beyond their bound.
#[must_use]
pub(crate) fn bounded_complete_distance(
    reference: &[Base],
    nominal_start: usize,
    read: &[Base],
    strand: BisulfiteStrand,
    maximum_distance: u8,
) -> Option<u8> {
    if read.is_empty() || read.len() > MAX_UNGAPPED_QUERY_BASES {
        return None;
    }
    let reference_end = nominal_start.checked_add(read.len())?;
    let reference = reference.get(nominal_start..reference_end)?;
    let semantics = strand_semantics(strand);
    let mut distance = 0_u8;
    match semantics.orientation() {
        AlignmentOrientation::Forward => {
            for (&reference_base, &query_base) in reference.iter().zip(read) {
                distance += u8::from(
                    !classify_bases(reference_base, query_base, semantics.cytosine_strand())
                        .is_zero_cost(),
                );
                if distance > maximum_distance {
                    return None;
                }
            }
        }
        AlignmentOrientation::Reverse => {
            for (&reference_base, query_base) in reference
                .iter()
                .zip(read.iter().rev().map(|base| base.complement()))
            {
                distance += u8::from(
                    !classify_bases(reference_base, query_base, semantics.cytosine_strand())
                        .is_zero_cost(),
                );
                if distance > maximum_distance {
                    return None;
                }
            }
        }
    }
    Some(distance)
}

impl UngappedProfile {
    /// Builds the allocation-free profile for an already selected reference
    /// origin.
    ///
    /// Returns `None` when the query exceeds the supported fixed capacity or
    /// when its nominal origin is not an existing reference base. A query may
    /// extend beyond the right edge; such positions become endpoint barriers
    /// so terminal clipping can still retain an in-bounds prefix.
    #[must_use]
    pub fn new(
        reference: &[Base],
        nominal_start: usize,
        read: &[Base],
        strand: BisulfiteStrand,
    ) -> Option<Self> {
        if read.is_empty()
            || read.len() > MAX_UNGAPPED_QUERY_BASES
            || nominal_start >= reference.len()
        {
            return None;
        }

        let semantics = strand_semantics(strand);
        let mut profile = Self {
            read_length: read.len(),
            nominal_start,
            orientation: semantics.orientation(),
            complete_reference: nominal_start
                .checked_add(read.len())
                .is_some_and(|end| end <= reference.len()),
            mismatch_prefix: [0; MAX_UNGAPPED_QUERY_BASES + 1],
            barrier_prefix: [0; MAX_UNGAPPED_QUERY_BASES + 1],
        };

        for position in 0..read.len() {
            let query_base = match semantics.orientation() {
                AlignmentOrientation::Forward => read[position],
                AlignmentOrientation::Reverse => read[read.len() - 1 - position].complement(),
            };
            if let Some(reference_base) = nominal_start
                .checked_add(position)
                .and_then(|offset| reference.get(offset))
                .copied()
            {
                let mismatch =
                    !classify_bases(reference_base, query_base, semantics.cytosine_strand())
                        .is_zero_cost();
                profile.mismatch_prefix[position + 1] =
                    profile.mismatch_prefix[position] + u16::from(mismatch);
                profile.barrier_prefix[position + 1] =
                    profile.barrier_prefix[position] + u16::from(reference_base == Base::N);
            } else {
                profile.mismatch_prefix[position + 1] = profile.mismatch_prefix[position] + 1;
                profile.barrier_prefix[position + 1] = profile.barrier_prefix[position] + 1;
            }
        }
        Some(profile)
    }

    /// Returns the whole-query distance when the complete reference span is
    /// present and its distance does not exceed `maximum_distance`.
    ///
    /// An `N` base is a unit mismatch here, matching global alignment
    /// semantics; bounded endpoint search treats it as a barrier instead.
    #[must_use]
    pub fn complete_distance(&self, maximum_distance: u8) -> Option<u8> {
        if !self.complete_reference {
            return None;
        }
        let distance = u8::try_from(self.mismatch_prefix[self.read_length]).ok()?;
        (distance <= maximum_distance).then_some(distance)
    }

    /// Returns one barrier-free ungapped endpoint after removing the supplied
    /// oriented terminal lengths.
    #[must_use]
    pub fn endpoint(
        &self,
        oriented_left_clip: usize,
        oriented_right_clip: usize,
    ) -> Option<UngappedEndpoint> {
        let clipped = oriented_left_clip.checked_add(oriented_right_clip)?;
        if clipped > self.read_length {
            return None;
        }
        let aligned_end = self.read_length - oriented_right_clip;
        if self.barrier_prefix[aligned_end] != self.barrier_prefix[oriented_left_clip] {
            return None;
        }
        let distance = u8::try_from(
            self.mismatch_prefix[aligned_end] - self.mismatch_prefix[oriented_left_clip],
        )
        .ok()?;
        let (query_start, query_end) = match self.orientation {
            AlignmentOrientation::Forward => (oriented_left_clip, aligned_end),
            AlignmentOrientation::Reverse => (
                oriented_right_clip,
                self.read_length.saturating_sub(oriented_left_clip),
            ),
        };
        Some(UngappedEndpoint {
            reference_start: self.nominal_start.checked_add(oriented_left_clip)?,
            reference_end: self.nominal_start.checked_add(aligned_end)?,
            query_start,
            query_end,
            distance,
            oriented_left_clip,
            oriented_right_clip,
        })
    }

    /// Chooses the best bounded terminally clipped endpoint under an explicit
    /// caller-supplied score and admission policy.
    #[must_use]
    // Admission, tie-breaking, and barrier handling share one bounded scan;
    // splitting it would obscure the frozen endpoint policy.
    #[allow(clippy::too_many_lines)]
    pub fn best_bounded_semiglobal(
        &self,
        config: BoundedSemiglobalConfig,
    ) -> Option<BoundedSemiglobalAlignment> {
        if self.read_length < config.minimum_aligned_bases {
            return None;
        }
        let maximum_clip = config.maximum_clip_bases.min(
            self.read_length
                .saturating_sub(config.minimum_aligned_bases),
        );
        let mut best: Option<UngappedSelection> = None;
        let mut consider = |oriented_left_clip: usize, oriented_right_clip: usize| {
            let clipped = oriented_left_clip.saturating_add(oriented_right_clip);
            if clipped == 0
                || self.read_length.saturating_sub(clipped) < config.minimum_aligned_bases
            {
                return false;
            }
            let Some(endpoint) = self.endpoint(oriented_left_clip, oriented_right_clip) else {
                return false;
            };
            if endpoint.distance > config.maximum_edit_distance {
                return false;
            }
            let clipped_u8 = u8::try_from(clipped).unwrap_or(u8::MAX);
            let score = endpoint
                .distance
                .saturating_mul(config.edit_penalty)
                .saturating_add(clipped_u8.saturating_mul(config.clip_penalty));
            let admission_score = endpoint
                .distance
                .saturating_mul(config.admission_edit_penalty)
                .saturating_add(clipped_u8.saturating_mul(config.admission_clip_penalty));
            if admission_score > config.maximum_admission_score {
                return false;
            }
            let key = (
                score,
                clipped,
                endpoint.distance,
                oriented_left_clip,
                oriented_right_clip,
            );
            if best.as_ref().is_none_or(|(current, _)| key < *current) {
                best = Some((key, endpoint));
            }
            true
        };

        if self.barrier_prefix[self.read_length] == 0
            && maximum_clip.saturating_mul(2)
                <= self
                    .read_length
                    .saturating_sub(config.minimum_aligned_bases)
        {
            let mut exact_error_right = [None; MAX_UNGAPPED_QUERY_BASES + 1];
            let total_mismatches = usize::from(self.mismatch_prefix[self.read_length]);
            for right_clip in 0..=maximum_clip {
                let removed = usize::from(
                    self.mismatch_prefix[self.read_length]
                        - self.mismatch_prefix[self.read_length.saturating_sub(right_clip)],
                );
                if exact_error_right[removed].is_none() {
                    exact_error_right[removed] = Some(right_clip);
                }
            }
            let mut at_least_error_right = [None; MAX_UNGAPPED_QUERY_BASES + 1];
            let mut suffix_best: Option<usize> = None;
            for removed in (0..=maximum_clip).rev() {
                if let Some(incoming) = exact_error_right[removed] {
                    let incoming_key = (
                        i32::try_from(incoming).unwrap_or(i32::MAX)
                            - i32::from(config.edit_penalty)
                                * i32::try_from(removed).unwrap_or(i32::MAX),
                        incoming,
                    );
                    if suffix_best.is_none_or(|current| {
                        let current_removed = usize::from(
                            self.mismatch_prefix[self.read_length]
                                - self.mismatch_prefix[self.read_length.saturating_sub(current)],
                        );
                        incoming_key
                            < (
                                i32::try_from(current).unwrap_or(i32::MAX)
                                    - i32::from(config.edit_penalty)
                                        * i32::try_from(current_removed).unwrap_or(i32::MAX),
                                current,
                            )
                    }) {
                        suffix_best = Some(incoming);
                    }
                }
                at_least_error_right[removed] = suffix_best;
            }
            for right_clip in 1..=maximum_clip {
                let _ = consider(0, right_clip);
            }
            for left_clip in 1..=maximum_clip {
                let left_mismatches = usize::from(self.mismatch_prefix[left_clip]);
                let after_left = total_mismatches.saturating_sub(left_mismatches);
                let required = after_left.saturating_sub(usize::from(config.maximum_edit_distance));
                if required <= maximum_clip
                    && let Some(right_clip) = at_least_error_right[required]
                    && consider(left_clip, right_clip)
                {
                    continue;
                }
                for distance in 0..=usize::from(config.maximum_edit_distance) {
                    let Some(right_mismatches) = after_left.checked_sub(distance) else {
                        continue;
                    };
                    if right_mismatches <= maximum_clip
                        && let Some(right_clip) = exact_error_right[right_mismatches]
                    {
                        let _ = consider(left_clip, right_clip);
                    }
                }
            }
        } else {
            for oriented_left_clip in 0..=maximum_clip {
                for oriented_right_clip in 0..=maximum_clip {
                    let _ = consider(oriented_left_clip, oriented_right_clip);
                }
            }
        }

        best.map(|((score, ..), endpoint)| BoundedSemiglobalAlignment { endpoint, score })
    }
}

/// One barrier-free ungapped endpoint in forward-reference coordinates.
#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub struct UngappedEndpoint {
    reference_start: usize,
    reference_end: usize,
    query_start: usize,
    query_end: usize,
    distance: u8,
    oriented_left_clip: usize,
    oriented_right_clip: usize,
}

impl UngappedEndpoint {
    /// Inclusive forward-reference start.
    #[must_use]
    pub const fn reference_start(self) -> usize {
        self.reference_start
    }

    /// Exclusive forward-reference end.
    #[must_use]
    pub const fn reference_end(self) -> usize {
        self.reference_end
    }

    /// Inclusive query start in original sequencing orientation.
    #[must_use]
    pub const fn query_start(self) -> usize {
        self.query_start
    }

    /// Exclusive query end in original sequencing orientation.
    #[must_use]
    pub const fn query_end(self) -> usize {
        self.query_end
    }

    /// Conversion-aware unit mismatch count.
    #[must_use]
    pub const fn distance(self) -> u8 {
        self.distance
    }

    /// Left clip after orienting the query to forward-reference order.
    #[must_use]
    pub const fn oriented_left_clip(self) -> usize {
        self.oriented_left_clip
    }

    /// Right clip after orienting the query to forward-reference order.
    #[must_use]
    pub const fn oriented_right_clip(self) -> usize {
        self.oriented_right_clip
    }
}

/// Explicit scoring and admission parameters for bounded ungapped endpoint
/// selection.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct BoundedSemiglobalConfig {
    maximum_edit_distance: u8,
    maximum_clip_bases: usize,
    minimum_aligned_bases: usize,
    edit_penalty: u8,
    clip_penalty: u8,
    admission_edit_penalty: u8,
    admission_clip_penalty: u8,
    maximum_admission_score: u8,
}

impl BoundedSemiglobalConfig {
    /// Constructs an explicit endpoint scoring policy.
    #[must_use]
    #[allow(clippy::too_many_arguments)]
    pub const fn new(
        maximum_edit_distance: u8,
        maximum_clip_bases: usize,
        minimum_aligned_bases: usize,
        edit_penalty: u8,
        clip_penalty: u8,
        admission_edit_penalty: u8,
        admission_clip_penalty: u8,
        maximum_admission_score: u8,
    ) -> Self {
        Self {
            maximum_edit_distance,
            maximum_clip_bases,
            minimum_aligned_bases,
            edit_penalty,
            clip_penalty,
            admission_edit_penalty,
            admission_clip_penalty,
            maximum_admission_score,
        }
    }
}

/// Selected bounded semi-global endpoint and its caller-defined score.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct BoundedSemiglobalAlignment {
    endpoint: UngappedEndpoint,
    score: u8,
}

impl BoundedSemiglobalAlignment {
    /// Returns the selected endpoint.
    #[must_use]
    pub const fn endpoint(self) -> UngappedEndpoint {
        self.endpoint
    }

    /// Returns the saturated endpoint score.
    #[must_use]
    pub const fn score(self) -> u8 {
        self.score
    }
}

/// Builds the reference-code masks accepted by each query base under one
/// bisulfite cytosine strand.
#[must_use]
pub fn reference_masks_by_query(strand: CytosineStrand) -> [u8; 5] {
    let mut masks = [0_u8; 5];
    for query in Base::ALL {
        for reference in Base::ALL {
            if classify_bases(reference, query, strand).is_zero_cost() {
                masks[usize::from(query.storage_code())] |= 1 << reference.storage_code();
            }
        }
    }
    masks
}

#[cfg(test)]
mod tests {
    use super::*;
    use bsbit_core::sequence::{NormalizedSequence, normalize_dna};

    fn sequence(input: &str) -> NormalizedSequence {
        normalize_dna(input.as_bytes()).expect("valid test sequence")
    }

    #[test]
    fn reverse_endpoints_are_reported_in_sequencing_coordinates() {
        let reference = sequence(&format!("{}{}", "A".repeat(70), "C".repeat(80)));
        let read = sequence(&format!("{}{}", "G".repeat(75), "T".repeat(5)));
        let profile =
            UngappedProfile::new(reference.bases(), 70, read.bases(), BisulfiteStrand::OB)
                .expect("profile");
        let endpoint = profile.endpoint(5, 0).expect("endpoint");
        assert_eq!((endpoint.query_start(), endpoint.query_end()), (0, 75));
        assert_eq!(
            (endpoint.reference_start(), endpoint.reference_end()),
            (75, 150)
        );
    }

    #[test]
    fn complete_distance_and_endpoint_barriers_have_distinct_semantics() {
        let reference = sequence("ACNGT");
        let read = sequence("ACAGT");
        let profile = UngappedProfile::new(reference.bases(), 0, read.bases(), BisulfiteStrand::OT)
            .expect("profile");
        assert_eq!(profile.complete_distance(1), Some(1));
        assert_eq!(profile.endpoint(0, 0), None);
    }

    #[test]
    fn masks_match_the_complete_relation_table() {
        for strand in [CytosineStrand::Top, CytosineStrand::Bottom] {
            let masks = reference_masks_by_query(strand);
            for query in Base::ALL {
                for reference in Base::ALL {
                    let actual = masks[usize::from(query.storage_code())]
                        & (1 << reference.storage_code())
                        != 0;
                    assert_eq!(
                        actual,
                        classify_bases(reference, query, strand).is_zero_cost()
                    );
                }
            }
        }
    }

    #[test]
    fn bounded_complete_distance_matches_profile_for_every_short_sequence() {
        let strands = [
            BisulfiteStrand::OT,
            BisulfiteStrand::OB,
            BisulfiteStrand::CTOT,
            BisulfiteStrand::CTOB,
        ];
        let all_pairs: Vec<(Base, Base)> = Base::ALL
            .into_iter()
            .flat_map(|reference| Base::ALL.into_iter().map(move |query| (reference, query)))
            .collect();
        for &strand in &strands {
            for maximum_distance in 0..=3 {
                for &(reference_first, query_first) in &all_pairs {
                    for &(reference_second, query_second) in &all_pairs {
                        let reference = [Base::A, reference_first, reference_second, Base::T];
                        let read = [query_first, query_second];
                        let expected = UngappedProfile::new(&reference, 1, &read, strand)
                            .and_then(|profile| profile.complete_distance(maximum_distance));
                        assert_eq!(
                            bounded_complete_distance(
                                &reference,
                                1,
                                &read,
                                strand,
                                maximum_distance,
                            ),
                            expected,
                            "strand={strand:?}, maximum_distance={maximum_distance}, reference={reference:?}, read={read:?}",
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn bounded_complete_distance_preserves_profile_boundaries() {
        let reference = [Base::A; MAX_UNGAPPED_QUERY_BASES];
        assert_eq!(
            bounded_complete_distance(&reference, 0, &[], BisulfiteStrand::OT, 3),
            None
        );
        assert_eq!(
            bounded_complete_distance(
                &reference,
                0,
                &[Base::A; MAX_UNGAPPED_QUERY_BASES + 1],
                BisulfiteStrand::OT,
                3,
            ),
            None
        );
        assert_eq!(
            bounded_complete_distance(
                &reference,
                reference.len(),
                &[Base::A],
                BisulfiteStrand::OT,
                3
            ),
            None
        );
        assert_eq!(
            bounded_complete_distance(
                &reference,
                reference.len() - 1,
                &[Base::A, Base::A],
                BisulfiteStrand::OT,
                3,
            ),
            None
        );
    }
}
