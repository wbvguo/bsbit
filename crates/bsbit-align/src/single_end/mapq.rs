//! Single-end mapping-quality evidence and confidence policy.
//!
//! The mapper supplies evidence already observed while selecting one read;
//! MAPQ calculation performs no additional index search or verification.

/// Already-observed evidence used to score one unique single-read origin.
///
/// The policy intentionally consumes only the candidate and verification
/// frontier that selected the alignment.  Producing MAPQ therefore adds no
/// second FM-index search and no second verification pass.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct SingleMapqEvidence {
    pub(crate) best_distance: u8,
    pub(crate) second_best_distance: Option<u8>,
    pub(crate) verified_distance_limit: u8,
    pub(crate) located_rows: u64,
    pub(crate) distinct_candidate_starts: u64,
    pub(crate) verified_placements: u64,
    pub(crate) first_seed_hits: u64,
    pub(crate) first_seed_bases: u64,
    pub(crate) direct_singleton: bool,
}

/// Returns an integer-only, conservative single-read MAPQ from one retained
/// search frontier.  Caps are evidence gates: edit separation supplies the
/// raw confidence, while seed length, repeat pressure, and incomplete search
/// can only lower it.
#[must_use]
pub(crate) fn single_mapping_quality_from_evidence(evidence: SingleMapqEvidence) -> u8 {
    let separation = evidence.second_best_distance.map_or_else(
        || {
            evidence
                .verified_distance_limit
                .saturating_add(1)
                .saturating_sub(evidence.best_distance)
        },
        |second| second.saturating_sub(evidence.best_distance),
    );

    // A unique bounded result is the Q10 floor.  Q20 additionally excludes
    // weak edit separation, high edit burden, and a repeat-heavy frontier.
    // These are caps rather than bonuses, so adding adverse evidence can never
    // raise confidence.
    let mut mapping_quality = match (evidence.best_distance, separation) {
        (_, 0) => 1,
        (5.., _) => 10,
        (4, _) | (_, 1) => 15,
        _ => 20,
    };
    if evidence.first_seed_hits > 64
        || evidence.located_rows > 256
        || evidence.distinct_candidate_starts > 64
        || evidence.verified_placements > 64
    {
        mapping_quality = mapping_quality.min(10);
    }

    // A singleton interval reached after a short projected suffix is stronger
    // uniqueness evidence than one that needs a long suffix: the FM interval
    // becomes unique earlier.  Reuse that already-observed fact for Q30/Q40;
    // do not continue the search solely to manufacture MAPQ evidence.
    let short_singleton = evidence.direct_singleton
        && evidence.first_seed_hits == 1
        && (16..=46).contains(&evidence.first_seed_bases);
    if mapping_quality >= 20 && short_singleton {
        mapping_quality = 30;
        if evidence.best_distance == 0 && evidence.first_seed_bases <= 41 {
            mapping_quality = 40;
        }
    }
    mapping_quality
}

#[cfg(test)]
#[path = "../../tests/whitebox/single_mapq.rs"]
mod whitebox;
