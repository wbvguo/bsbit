//! Alignment-owned seed policy over neutral index query primitives.

use bsbit_core::bisulfite::BisulfiteStrand;
use bsbit_index::reference::{
    COMBINED_EXACT_LOOKUP_BASES, CombinedIndexBackendError, CombinedIndexCoordinate,
    CombinedIndexInterval, CombinedIndexQuery, CombinedIndexQueryError, ReferenceAllocation,
    ReferenceIndex, ReferenceLocateError, ReferenceLocateInvariant, ReferenceLocateMetrics,
    ReferenceQueryError,
};
use bsbit_index::storage::fm::{ProjectedBase, SearchBase};

/// One alignment-selected interval from the combined image.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct CombinedSeedMatches {
    interval: CombinedIndexInterval,
    matched_bases: u64,
}

impl CombinedSeedMatches {
    pub(crate) const fn exact_hit_count(self) -> u64 {
        self.interval.len()
    }

    pub(crate) const fn matched_bases(self) -> u64 {
        self.matched_bases
    }
}

/// One combined-index coordinate interpreted as an alignment seed hit.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct CombinedSeedHit {
    contig_ordinal: u64,
    strand: BisulfiteStrand,
    start: u64,
}

impl CombinedSeedHit {
    pub(crate) const fn contig_ordinal(self) -> u64 {
        self.contig_ordinal
    }

    pub(crate) const fn strand(self) -> BisulfiteStrand {
        self.strand
    }

    pub(crate) const fn start(self) -> u64 {
        self.start
    }
}

#[derive(Clone, Copy)]
struct CombinedSuffixState {
    interval: CombinedIndexInterval,
    matched_bases: usize,
    remaining_prefix_bases: usize,
    finished: bool,
}

impl CombinedSuffixState {
    const fn new(interval: CombinedIndexInterval, remaining_prefix_bases: usize) -> Self {
        Self {
            interval,
            matched_bases: COMBINED_EXACT_LOOKUP_BASES,
            remaining_prefix_bases,
            finished: interval.len() == 1 || remaining_prefix_bases == 0,
        }
    }

    fn accept(&mut self, extended: CombinedIndexInterval) {
        if extended.is_empty() {
            self.finished = true;
            return;
        }
        self.interval = extended;
        self.matched_bases += 1;
        self.remaining_prefix_bases -= 1;
        self.finished = extended.len() == 1 || self.remaining_prefix_bases == 0;
    }

    fn finish(self) -> Result<CombinedSeedMatches, ReferenceQueryError> {
        Ok(CombinedSeedMatches {
            interval: self.interval,
            matched_bases: u64::try_from(self.matched_bases).map_err(|_| {
                ReferenceQueryError::PatternLengthNotRepresentable {
                    pattern_len: self.matched_bases,
                }
            })?,
        })
    }
}

trait CombinedQuerySymbol: Copy {
    fn lookup(
        query: CombinedIndexQuery<'_>,
        symbols: &[Self],
    ) -> Result<Option<CombinedIndexInterval>, CombinedIndexQueryError>;

    fn extend(
        query: CombinedIndexQuery<'_>,
        interval: CombinedIndexInterval,
        symbol: Self,
    ) -> Result<CombinedIndexInterval, CombinedIndexQueryError>;

    fn extend_batch(
        query: CombinedIndexQuery<'_>,
        intervals: &[CombinedIndexInterval],
        symbols: &[Self],
        output: &mut [CombinedIndexInterval],
    ) -> Result<(), CombinedIndexQueryError>;
}

impl CombinedQuerySymbol for SearchBase {
    fn lookup(
        query: CombinedIndexQuery<'_>,
        symbols: &[Self],
    ) -> Result<Option<CombinedIndexInterval>, CombinedIndexQueryError> {
        query.lookup_interval(symbols)
    }

    fn extend(
        query: CombinedIndexQuery<'_>,
        interval: CombinedIndexInterval,
        symbol: Self,
    ) -> Result<CombinedIndexInterval, CombinedIndexQueryError> {
        query.backward_extend(interval, symbol)
    }

    fn extend_batch(
        query: CombinedIndexQuery<'_>,
        intervals: &[CombinedIndexInterval],
        symbols: &[Self],
        output: &mut [CombinedIndexInterval],
    ) -> Result<(), CombinedIndexQueryError> {
        query.backward_extend_intervals(intervals, symbols, output)
    }
}

impl CombinedQuerySymbol for ProjectedBase {
    fn lookup(
        query: CombinedIndexQuery<'_>,
        symbols: &[Self],
    ) -> Result<Option<CombinedIndexInterval>, CombinedIndexQueryError> {
        query.lookup_projected_interval(symbols)
    }

    fn extend(
        query: CombinedIndexQuery<'_>,
        interval: CombinedIndexInterval,
        symbol: Self,
    ) -> Result<CombinedIndexInterval, CombinedIndexQueryError> {
        query.backward_extend_projected(interval, symbol)
    }

    fn extend_batch(
        query: CombinedIndexQuery<'_>,
        intervals: &[CombinedIndexInterval],
        symbols: &[Self],
        output: &mut [CombinedIndexInterval],
    ) -> Result<(), CombinedIndexQueryError> {
        query.backward_extend_projected_intervals(intervals, symbols, output)
    }
}

/// Alignment-policy methods kept beside candidate discovery rather than the
/// reference owner. The receiver syntax keeps call sites readable while this
/// trait remains crate-private.
pub(crate) trait CombinedSearchReferenceExt {
    fn combined_exact_seed(
        &self,
        pattern: &[SearchBase],
    ) -> Result<Option<CombinedSeedMatches>, ReferenceQueryError>;

    fn combined_maximal_suffix_projected(
        &self,
        pattern: &[ProjectedBase],
        minimum_suffix_bases: usize,
    ) -> Result<Option<CombinedSeedMatches>, ReferenceQueryError>;

    fn combined_maximal_suffix_projected_two_lanes(
        &self,
        patterns: [&[ProjectedBase]; 2],
        minimum_suffix_bases: usize,
    ) -> Result<[Option<CombinedSeedMatches>; 2], ReferenceQueryError>;

    fn combined_maximal_suffix_projected_wavefront(
        &self,
        patterns: &[&[ProjectedBase]],
        minimum_suffix_bases: usize,
    ) -> Result<Vec<Option<CombinedSeedMatches>>, ReferenceQueryError>;

    fn combined_maximal_suffix_projected_wavefront_into(
        &self,
        patterns: &[&[ProjectedBase]],
        minimum_suffix_bases: usize,
        output: &mut [Option<CombinedSeedMatches>],
    ) -> Result<(), ReferenceQueryError>;

    fn visit_combined_seed(
        &self,
        matches: CombinedSeedMatches,
        seed_offset: u64,
        query_len: u64,
        visitor: &mut dyn FnMut(CombinedSeedHit) -> bool,
    ) -> Result<ReferenceLocateMetrics, ReferenceLocateError>;

    fn visit_combined_seed_two_lanes_complete(
        &self,
        matches: [CombinedSeedMatches; 2],
        seed_offsets: [u64; 2],
        query_lens: [u64; 2],
        visitor: &mut dyn FnMut(usize, CombinedSeedHit),
    ) -> Result<[ReferenceLocateMetrics; 2], ReferenceLocateError>;
}

impl CombinedSearchReferenceExt for ReferenceIndex {
    fn combined_exact_seed(
        &self,
        pattern: &[SearchBase],
    ) -> Result<Option<CombinedSeedMatches>, ReferenceQueryError> {
        if pattern.len() < COMBINED_EXACT_LOOKUP_BASES {
            return Ok(None);
        }
        let Some(query) = self.combined_index_query() else {
            return Ok(None);
        };
        query
            .exact_interval(pattern)
            .map_err(|error| combined_query_error(&error))?
            .map(|interval| {
                Ok(CombinedSeedMatches {
                    interval,
                    matched_bases: u64::try_from(pattern.len()).map_err(|_| {
                        ReferenceQueryError::PatternLengthNotRepresentable {
                            pattern_len: pattern.len(),
                        }
                    })?,
                })
            })
            .transpose()
    }

    fn combined_maximal_suffix_projected(
        &self,
        pattern: &[ProjectedBase],
        minimum_suffix_bases: usize,
    ) -> Result<Option<CombinedSeedMatches>, ReferenceQueryError> {
        combined_maximal_suffix(self, pattern, minimum_suffix_bases)
    }

    fn combined_maximal_suffix_projected_two_lanes(
        &self,
        patterns: [&[ProjectedBase]; 2],
        minimum_suffix_bases: usize,
    ) -> Result<[Option<CombinedSeedMatches>; 2], ReferenceQueryError> {
        let seeds =
            combined_maximal_suffix_projected_wavefront(self, &patterns, minimum_suffix_bases)?;
        Ok([seeds[0], seeds[1]])
    }

    fn combined_maximal_suffix_projected_wavefront(
        &self,
        patterns: &[&[ProjectedBase]],
        minimum_suffix_bases: usize,
    ) -> Result<Vec<Option<CombinedSeedMatches>>, ReferenceQueryError> {
        combined_maximal_suffix_projected_wavefront(self, patterns, minimum_suffix_bases)
    }

    fn combined_maximal_suffix_projected_wavefront_into(
        &self,
        patterns: &[&[ProjectedBase]],
        minimum_suffix_bases: usize,
        output: &mut [Option<CombinedSeedMatches>],
    ) -> Result<(), ReferenceQueryError> {
        combined_maximal_suffix_projected_wavefront_into(
            self,
            patterns,
            minimum_suffix_bases,
            output,
        )
    }

    fn visit_combined_seed(
        &self,
        matches: CombinedSeedMatches,
        seed_offset: u64,
        query_len: u64,
        visitor: &mut dyn FnMut(CombinedSeedHit) -> bool,
    ) -> Result<ReferenceLocateMetrics, ReferenceLocateError> {
        let query = combined_query_for_locate(self)?;
        query
            .visit_interval(
                matches.interval,
                matches.matched_bases,
                seed_offset,
                query_len,
                &mut |coordinate| visitor(combined_hit(coordinate)),
            )
            .map_err(combined_locate_error)
    }

    fn visit_combined_seed_two_lanes_complete(
        &self,
        matches: [CombinedSeedMatches; 2],
        seed_offsets: [u64; 2],
        query_lens: [u64; 2],
        visitor: &mut dyn FnMut(usize, CombinedSeedHit),
    ) -> Result<[ReferenceLocateMetrics; 2], ReferenceLocateError> {
        let query = combined_query_for_locate(self)?;
        query
            .visit_raw_intervals_two_lanes_complete(
                [matches[0].interval, matches[1].interval],
                &mut |lane, position| {
                    if let Some(coordinate) = query.recover_coordinate(
                        position,
                        matches[lane].matched_bases,
                        seed_offsets[lane],
                        query_lens[lane],
                    ) {
                        visitor(lane, combined_hit(coordinate));
                    }
                },
            )
            .map_err(combined_locate_error)
    }
}

fn combined_maximal_suffix<T: CombinedQuerySymbol>(
    reference: &ReferenceIndex,
    pattern: &[T],
    minimum_suffix_bases: usize,
) -> Result<Option<CombinedSeedMatches>, ReferenceQueryError> {
    if minimum_suffix_bases < COMBINED_EXACT_LOOKUP_BASES || pattern.len() < minimum_suffix_bases {
        return Ok(None);
    }
    let Some(query) = reference.combined_index_query() else {
        return Ok(None);
    };
    let suffix_start = pattern.len() - COMBINED_EXACT_LOOKUP_BASES;
    let Some(interval) =
        T::lookup(query, &pattern[suffix_start..]).map_err(|error| combined_query_error(&error))?
    else {
        return Ok(None);
    };
    let mut state = CombinedSuffixState::new(interval, suffix_start);
    while !state.finished {
        let symbol = pattern[state.remaining_prefix_bases - 1];
        let extended = T::extend(query, state.interval, symbol)
            .map_err(|error| combined_query_error(&error))?;
        state.accept(extended);
    }
    state.finish().map(Some)
}

fn combined_maximal_suffix_projected_wavefront(
    reference: &ReferenceIndex,
    patterns: &[&[ProjectedBase]],
    minimum_suffix_bases: usize,
) -> Result<Vec<Option<CombinedSeedMatches>>, ReferenceQueryError> {
    const MAX_LANES: usize = 64;
    if patterns.len() > MAX_LANES {
        return combined_maximal_suffix_wavefront(reference, patterns, minimum_suffix_bases);
    }
    let mut output = Vec::new();
    output
        .try_reserve_exact(patterns.len())
        .map_err(|_| query_allocation_error(patterns.len()))?;
    output.resize(patterns.len(), None);
    combined_maximal_suffix_projected_wavefront_into(
        reference,
        patterns,
        minimum_suffix_bases,
        &mut output,
    )?;
    Ok(output)
}

fn combined_maximal_suffix_projected_wavefront_into(
    reference: &ReferenceIndex,
    patterns: &[&[ProjectedBase]],
    minimum_suffix_bases: usize,
    output: &mut [Option<CombinedSeedMatches>],
) -> Result<(), ReferenceQueryError> {
    const MAX_LANES: usize = 64;
    if patterns.len() != output.len() || patterns.len() > MAX_LANES {
        return Err(ReferenceQueryError::CombinedIndex {
            source: CombinedIndexBackendError::Structure,
        });
    }
    output.fill(None);
    let Some(query) = reference.combined_index_query() else {
        return Ok(());
    };
    // Alignment chooses singleton completion as its stopping rule. The index
    // backend receives that rule explicitly and owns only the memory schedule
    // for the dense lookup and adjacent rank rounds.
    let mut intervals = [None; MAX_LANES];
    query
        .resolve_projected_suffix_intervals(
            patterns,
            minimum_suffix_bases,
            1,
            &mut intervals[..patterns.len()],
        )
        .map_err(|error| combined_query_error(&error))?;
    for (destination, result) in output.iter_mut().zip(intervals) {
        *destination = result.map(|(interval, matched_bases)| CombinedSeedMatches {
            interval,
            matched_bases,
        });
    }
    Ok(())
}

fn combined_maximal_suffix_wavefront<T: CombinedQuerySymbol>(
    reference: &ReferenceIndex,
    patterns: &[&[T]],
    minimum_suffix_bases: usize,
) -> Result<Vec<Option<CombinedSeedMatches>>, ReferenceQueryError> {
    const MAX_LANES: usize = 64;
    let mut output = Vec::new();
    output
        .try_reserve_exact(patterns.len())
        .map_err(|_| query_allocation_error(patterns.len()))?;
    output.resize(patterns.len(), None);
    if patterns.len() > MAX_LANES {
        for (output, pattern) in output.iter_mut().zip(patterns) {
            *output = combined_maximal_suffix(reference, pattern, minimum_suffix_bases)?;
        }
        return Ok(output);
    }
    let Some(query) = reference.combined_index_query() else {
        return Ok(output);
    };
    let mut states = [None; MAX_LANES];
    for (lane, pattern) in patterns.iter().enumerate() {
        if minimum_suffix_bases < COMBINED_EXACT_LOOKUP_BASES
            || pattern.len() < minimum_suffix_bases
        {
            continue;
        }
        let suffix_start = pattern.len() - COMBINED_EXACT_LOOKUP_BASES;
        let Some(interval) = T::lookup(query, &pattern[suffix_start..])
            .map_err(|error| combined_query_error(&error))?
        else {
            continue;
        };
        states[lane] = Some(CombinedSuffixState::new(interval, suffix_start));
    }
    loop {
        let mut active_lanes = [0_usize; MAX_LANES];
        let mut active_count = 0_usize;
        for (lane, state) in states.iter().enumerate().take(patterns.len()) {
            if state.is_some_and(|state| !state.finished) {
                active_lanes[active_count] = lane;
                active_count += 1;
            }
        }
        if active_count == 0 {
            break;
        }
        let first_state = states[active_lanes[0]].expect("active suffix state exists");
        let mut intervals = [first_state.interval; MAX_LANES];
        let first_symbol = patterns[active_lanes[0]][first_state.remaining_prefix_bases - 1];
        let mut symbols = [first_symbol; MAX_LANES];
        let mut extended = [first_state.interval; MAX_LANES];
        for active in 0..active_count {
            let lane = active_lanes[active];
            let state = states[lane].expect("active suffix state exists");
            intervals[active] = state.interval;
            symbols[active] = patterns[lane][state.remaining_prefix_bases - 1];
        }
        T::extend_batch(
            query,
            &intervals[..active_count],
            &symbols[..active_count],
            &mut extended[..active_count],
        )
        .map_err(|error| combined_query_error(&error))?;
        for active in 0..active_count {
            states[active_lanes[active]]
                .as_mut()
                .expect("active suffix state exists")
                .accept(extended[active]);
        }
    }
    for lane in 0..patterns.len() {
        output[lane] = states[lane].map(CombinedSuffixState::finish).transpose()?;
    }
    Ok(output)
}

fn combined_query_for_locate(
    reference: &ReferenceIndex,
) -> Result<CombinedIndexQuery<'_>, ReferenceLocateError> {
    reference
        .combined_index_query()
        .ok_or_else(missing_run_error)
}

const fn combined_hit(coordinate: CombinedIndexCoordinate) -> CombinedSeedHit {
    CombinedSeedHit {
        contig_ordinal: coordinate.contig_ordinal(),
        strand: coordinate.strand(),
        start: coordinate.start(),
    }
}

fn combined_query_error(error: &CombinedIndexQueryError) -> ReferenceQueryError {
    let source = match *error {
        CombinedIndexQueryError::Backend(source) => source,
        CombinedIndexQueryError::ForeignInterval | CombinedIndexQueryError::Coordinate(_) => {
            CombinedIndexBackendError::Structure
        }
        _ => CombinedIndexBackendError::Structure,
    };
    ReferenceQueryError::CombinedIndex { source }
}

fn combined_locate_error(error: CombinedIndexQueryError) -> ReferenceLocateError {
    match error {
        CombinedIndexQueryError::ForeignInterval => ReferenceLocateError::ForeignMatches,
        CombinedIndexQueryError::Backend(source) => combined_backend_locate_error(source),
        CombinedIndexQueryError::Coordinate(source) => source,
        _ => combined_backend_locate_error(CombinedIndexBackendError::Structure),
    }
}

const fn combined_backend_locate_error(source: CombinedIndexBackendError) -> ReferenceLocateError {
    ReferenceLocateError::CombinedIndex { source }
}

const fn missing_run_error() -> ReferenceLocateError {
    ReferenceLocateError::Invariant {
        invariant: ReferenceLocateInvariant::MissingRun,
        expected: 1,
        observed: 0,
    }
}

fn query_allocation_error(elements: usize) -> ReferenceQueryError {
    ReferenceQueryError::AllocationFailed {
        allocation: ReferenceAllocation::ProjectedPattern,
        elements: u64::try_from(elements).unwrap_or(u64::MAX),
    }
}
