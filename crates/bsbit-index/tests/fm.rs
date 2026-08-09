//! Independent ground-truth and exhaustive tests for the exact FM reference backend.

#[path = "support/fm_oracle.rs"]
mod fm_oracle;

use std::collections::HashSet;
use std::sync::Arc;
use std::thread;

use bsbit_core::alphabet::Base;
use bsbit_index::storage::fm::{
    FmAllocation, FmBuildLimit, FmError, FmIndex, FmInterval, SearchBase, TextOffset,
};
use fm_oracle::{
    CANONICAL, CT_PROJECTED, GA_PROJECTED, OracleBwt, direct_hits, enumerate_strings, naive_bwt,
    naive_interval, naive_rank, naive_suffix_array, to_search_bases, to_u64,
};

fn build(raw: &[u8]) -> FmIndex {
    FmIndex::build_reference(&to_search_bases(raw), FmBuildLimit::MAX)
        .expect("small oracle index should build")
}

fn located_values(index: &FmIndex, interval: FmInterval) -> Vec<u64> {
    index
        .locate(interval)
        .expect("valid interval should locate")
        .iter()
        .map(|offset| offset.get())
        .collect()
}

fn search_values(index: &FmIndex, raw_pattern: &[u8]) -> (FmInterval, Vec<u64>) {
    let interval = index.exact_search(&to_search_bases(raw_pattern));
    let locations = located_values(index, interval);
    (interval, locations)
}

fn observed_bwt(index: &FmIndex) -> Vec<OracleBwt> {
    let mut observed = Vec::new();
    for row in 0..index.suffix_count() {
        let mut changed = Vec::new();
        for base in SearchBase::ALL {
            let before = index.rank(base, row).expect("rank boundary is valid");
            let after = index
                .rank(base, row + 1)
                .expect("next rank boundary is valid");
            match after - before {
                0 => {}
                1 => changed.push(base.as_ascii()),
                delta => panic!("rank changed by {delta} at row {row}"),
            }
        }
        match changed.as_slice() {
            [] => observed.push(OracleBwt::Sentinel),
            [base] => observed.push(OracleBwt::Base(*base)),
            _ => panic!("multiple canonical symbols occupy BWT row {row}"),
        }
    }
    observed
}

fn assert_named_search(
    index: &FmIndex,
    pattern: &[u8],
    expected_interval: (u64, u64),
    expected_locations: &[u64],
) {
    let (interval, locations) = search_values(index, pattern);
    assert_eq!(
        (interval.lower(), interval.upper()),
        expected_interval,
        "pattern {}",
        String::from_utf8_lossy(pattern)
    );
    assert_eq!(
        locations,
        expected_locations,
        "pattern {}",
        String::from_utf8_lossy(pattern)
    );
}

fn assert_singleton(base: u8) {
    let index = build(&[base]);
    assert_eq!(
        observed_bwt(&index),
        vec![OracleBwt::Base(base), OracleBwt::Sentinel]
    );
    assert_named_search(&index, b"", (0, 2), &[1, 0]);
    assert_named_search(&index, &[base], (1, 2), &[0]);
}

#[test]
fn named_ground_truth_has_exact_public_semantics() {
    let empty = build(b"");
    assert_eq!(empty.text_len(), 0);
    assert_eq!(empty.suffix_count(), 1);
    assert_eq!(observed_bwt(&empty), vec![OracleBwt::Sentinel]);
    assert_named_search(&empty, b"", (0, 1), &[0]);
    assert_named_search(&empty, b"A", (1, 1), &[]);

    for base in CANONICAL {
        assert_singleton(base);
    }

    let singleton = build(b"A");
    assert_named_search(&singleton, b"AA", (2, 2), &[]);
    assert_named_search(&singleton, b"C", (2, 2), &[]);

    let one_c = build(b"C");
    assert_named_search(&one_c, b"A", (1, 1), &[]);
    assert_named_search(&one_c, b"CA", (2, 2), &[]);
}

#[test]
fn repetitive_and_all_symbol_ground_truth_is_exact() {
    let repeated = build(b"AAAA");
    assert_eq!(
        observed_bwt(&repeated),
        vec![
            OracleBwt::Base(b'A'),
            OracleBwt::Base(b'A'),
            OracleBwt::Base(b'A'),
            OracleBwt::Base(b'A'),
            OracleBwt::Sentinel,
        ]
    );
    assert_named_search(&repeated, b"A", (1, 5), &[3, 2, 1, 0]);
    assert_named_search(&repeated, b"AA", (2, 5), &[2, 1, 0]);
    assert_named_search(&repeated, b"AAAA", (4, 5), &[0]);
    assert_named_search(&repeated, b"AAAAA", (5, 5), &[]);

    let alternating = build(b"ACAC");
    assert_named_search(&alternating, b"AC", (1, 3), &[2, 0]);
    assert_named_search(&alternating, b"ACA", (2, 3), &[0]);
    assert_named_search(&alternating, b"C", (3, 5), &[3, 1]);
    assert_named_search(&alternating, b"CA", (4, 5), &[1]);
    assert_named_search(&alternating, b"G", (5, 5), &[]);

    let all_symbols = build(b"ACGT");
    assert_eq!(
        observed_bwt(&all_symbols),
        vec![
            OracleBwt::Base(b'T'),
            OracleBwt::Sentinel,
            OracleBwt::Base(b'A'),
            OracleBwt::Base(b'C'),
            OracleBwt::Base(b'G'),
        ]
    );

    let non_palindrome = build(b"GACTA");
    assert_eq!(
        located_values(&non_palindrome, non_palindrome.full_interval()),
        vec![5, 4, 1, 2, 0, 3]
    );
    assert_eq!(
        observed_bwt(&non_palindrome),
        vec![
            OracleBwt::Base(b'A'),
            OracleBwt::Base(b'T'),
            OracleBwt::Base(b'G'),
            OracleBwt::Base(b'A'),
            OracleBwt::Sentinel,
            OracleBwt::Base(b'C'),
        ]
    );
    assert_named_search(&non_palindrome, b"A", (1, 3), &[4, 1]);
    assert_named_search(&non_palindrome, b"ACT", (2, 3), &[1]);
    assert_named_search(&non_palindrome, b"CTA", (3, 4), &[2]);
    assert_named_search(&non_palindrome, b"GACTA", (4, 5), &[0]);
}

#[test]
fn public_values_and_every_error_field_are_structured() {
    for (search, base, ascii) in [
        (SearchBase::A, Base::A, b'A'),
        (SearchBase::C, Base::C, b'C'),
        (SearchBase::G, Base::G, b'G'),
        (SearchBase::T, Base::T, b'T'),
    ] {
        assert_eq!(SearchBase::from_base(base), Some(search));
        assert_eq!(search.as_base(), base);
        assert_eq!(search.as_ascii(), ascii);
        assert_eq!(search.to_string(), char::from(ascii).to_string());
    }
    assert_eq!(SearchBase::from_base(Base::N), None);
    assert_eq!(FmBuildLimit::new(7).get(), 7);

    let empty_source = Vec::<SearchBase>::new();
    let Err(rejected) = FmIndex::build_reference(&empty_source, FmBuildLimit::new(0)) else {
        panic!("zero suffix limit unexpectedly built an index");
    };
    assert_eq!(
        rejected,
        FmError::BuildLimitExceeded {
            suffix_count: 1,
            max_suffix_count: 0,
        }
    );
    let index = FmIndex::build_reference(&empty_source, FmBuildLimit::new(1))
        .expect("one suffix row admits empty text");
    let one_source = [SearchBase::A];
    let Err(one_rejected) = FmIndex::build_reference(&one_source, FmBuildLimit::new(1)) else {
        panic!("one-row limit unexpectedly admitted two suffixes");
    };
    assert_eq!(
        one_rejected,
        FmError::BuildLimitExceeded {
            suffix_count: 2,
            max_suffix_count: 1,
        }
    );
    FmIndex::build_reference(&one_source, FmBuildLimit::new(2))
        .expect("exact suffix limit should admit nonempty text");
    assert_eq!(
        index.rank(SearchBase::A, 2),
        Err(FmError::RankBoundaryOutOfBounds {
            boundary: 2,
            suffix_count: 1,
        })
    );
    assert_eq!(
        index.interval(3, 2),
        Err(FmError::InvertedInterval { lower: 3, upper: 2 })
    );
    assert_eq!(
        index.interval(0, 2),
        Err(FmError::IntervalOutOfBounds {
            lower: 0,
            upper: 2,
            suffix_count: 1,
        })
    );

    assert_all_error_diagnostics();
}

fn assert_all_error_diagnostics() {
    let diagnostics = [
        (
            FmError::TextLengthNotRepresentable { text_len: 7 },
            "text length 7 is not representable as u64".to_owned(),
        ),
        (
            FmError::SuffixCountOverflow { text_len: u64::MAX },
            format!("text length {} cannot include a terminal suffix", u64::MAX),
        ),
        (
            FmError::RankRowCountOverflow {
                suffix_count: u64::MAX,
            },
            format!(
                "suffix count {} cannot include the final rank boundary",
                u64::MAX
            ),
        ),
        (
            FmError::BuildLimitExceeded {
                suffix_count: 9,
                max_suffix_count: 8,
            },
            "requested 9 suffix rows exceeds build limit 8".to_owned(),
        ),
        (
            FmError::RankBoundaryOutOfBounds {
                boundary: 10,
                suffix_count: 9,
            },
            "rank boundary 10 exceeds BWT length 9".to_owned(),
        ),
        (
            FmError::InvertedInterval { lower: 4, upper: 3 },
            "FM interval [4, 3) is inverted".to_owned(),
        ),
        (
            FmError::IntervalDomainMismatch {
                interval_suffix_count: 8,
                index_suffix_count: 9,
            },
            "FM interval domain has 8 suffix rows; index has 9".to_owned(),
        ),
        (
            FmError::IntervalOutOfBounds {
                lower: 0,
                upper: 10,
                suffix_count: 9,
            },
            "FM interval [0, 10) exceeds suffix count 9".to_owned(),
        ),
    ];
    for (error, expected) in diagnostics {
        assert_eq!(error.to_string(), expected);
        let error_ref: &dyn std::error::Error = &error;
        assert!(error_ref.source().is_none());
    }

    for component in [
        FmAllocation::SuffixArray,
        FmAllocation::RankClasses,
        FmAllocation::NextRankClasses,
        FmAllocation::Bwt,
        FmAllocation::RankPrefixes,
        FmAllocation::LocateResults,
    ] {
        let overflow = FmError::AllocationSizeOverflow {
            component,
            elements: 11,
            element_size: 13,
        };
        let failure = FmError::AllocationFailed {
            component,
            elements: 17,
        };
        assert_eq!(
            overflow.to_string(),
            format!("cannot size {component:?}: 11 elements of 13 bytes")
        );
        assert_eq!(
            failure.to_string(),
            format!("failed to reserve 17 elements for {component:?}")
        );
    }
}

#[test]
fn interval_suffix_domains_and_bounds_are_revalidated() {
    let larger = build(b"ACGT");
    let interval = larger
        .interval(0, larger.suffix_count())
        .expect("full larger interval is valid");
    let smaller = build(b"A");
    assert_eq!(
        smaller.locate(interval),
        Err(FmError::IntervalDomainMismatch {
            interval_suffix_count: 5,
            index_suffix_count: 2,
        })
    );

    let empty = smaller.interval(1, 1).expect("empty interval is valid");
    assert!(empty.is_empty());
    assert_eq!(empty.len(), 0);
    assert_eq!(empty.to_string(), "[1, 1)");
    assert_eq!(smaller.locate(empty), Ok(Vec::new()));

    let same_size_other_text = build(b"T");
    let row_coordinate = smaller.interval(0, 1).expect("row coordinate is valid");
    assert!(same_size_other_text.locate(row_coordinate).is_ok());
}

fn verify_text_representation(text: &[u8], index: &FmIndex, suffix_array: &[usize]) -> (u64, u64) {
    let oracle_bwt = naive_bwt(text, suffix_array);
    assert_eq!(observed_bwt(index), oracle_bwt, "{text:?}");

    let full = located_values(index, index.full_interval());
    let expected_sa = suffix_array.iter().copied().map(to_u64).collect::<Vec<_>>();
    assert_eq!(full, expected_sa, "{text:?}");

    let mut rank_cases = 0_u64;
    for boundary in 0..=index.suffix_count() {
        let boundary_storage = usize::try_from(boundary).expect("tiny boundary fits usize");
        for (search_base, raw_base) in SearchBase::ALL.into_iter().zip(CANONICAL) {
            assert_eq!(
                index
                    .rank(search_base, boundary)
                    .expect("exhaustive rank boundary is valid"),
                naive_rank(&oracle_bwt, raw_base, boundary_storage),
                "text={:?}, base={}, boundary={boundary}",
                text,
                char::from(raw_base)
            );
            rank_cases += 1;
        }
    }

    let mut interval_cases = 0_u64;
    for lower in 0..=index.suffix_count() {
        for upper in lower..=index.suffix_count() {
            let interval = index
                .interval(lower, upper)
                .expect("enumerated interval is valid");
            let observed = located_values(index, interval);
            let lower_storage = usize::try_from(lower).expect("tiny lower fits usize");
            let upper_storage = usize::try_from(upper).expect("tiny upper fits usize");
            let expected = suffix_array[lower_storage..upper_storage]
                .iter()
                .copied()
                .map(to_u64)
                .collect::<Vec<_>>();
            assert_eq!(observed, expected, "text={text:?}, interval={interval}");
            interval_cases += 1;
        }
    }
    (rank_cases, interval_cases)
}

fn assert_enumeration_complete(values: &[Vec<u8>], alphabet: &[u8], maximum_length: usize) {
    let unique = values.iter().map(Vec::as_slice).collect::<HashSet<_>>();
    assert_eq!(unique.len(), values.len());
    for length in 0..=maximum_length {
        let observed = values.iter().filter(|value| value.len() == length).count();
        let exponent = u32::try_from(length).expect("test length fits u32");
        assert_eq!(observed, alphabet.len().pow(exponent));
    }
    assert!(values.iter().flatten().all(|base| alphabet.contains(base)));
}

fn verify_search_case(text: &[u8], index: &FmIndex, suffix_array: &[usize], pattern: &[u8]) {
    let original_pattern = pattern.to_vec();
    let search_pattern = to_search_bases(pattern);
    let original_search_pattern = search_pattern.clone();
    let interval = index.exact_search(&search_pattern);
    assert_eq!(search_pattern, original_search_pattern);

    let expected_interval = naive_interval(text, suffix_array, pattern);
    assert_eq!(
        (interval.lower(), interval.upper()),
        expected_interval,
        "text={text:?}, pattern={pattern:?}"
    );
    let observed = located_values(index, interval);
    let lower = usize::try_from(expected_interval.0).expect("tiny lower fits usize");
    let upper = usize::try_from(expected_interval.1).expect("tiny upper fits usize");
    let expected_order = suffix_array[lower..upper]
        .iter()
        .copied()
        .map(to_u64)
        .collect::<Vec<_>>();
    assert_eq!(
        observed, expected_order,
        "text={text:?}, pattern={pattern:?}"
    );

    let mut observed_set = observed;
    observed_set.sort_unstable();
    let mut expected_set = direct_hits(text, pattern);
    expected_set.sort_unstable();
    assert_eq!(
        observed_set, expected_set,
        "text={text:?}, pattern={pattern:?}"
    );
    if pattern.is_empty() {
        assert!(expected_set.contains(&to_u64(text.len())));
    } else {
        assert!(!expected_set.contains(&to_u64(text.len())));
    }
    assert_eq!(pattern, original_pattern);
}

fn verify_all_substrings(text: &[u8], index: &FmIndex, suffix_array: &[usize]) -> u64 {
    let mut cases = 0_u64;
    for start in 0..text.len() {
        for end in start + 1..=text.len() {
            verify_search_case(text, index, suffix_array, &text[start..end]);
            cases += 1;
        }
    }
    cases
}

#[test]
fn canonical_exhaustive_oracle_covers_sa_bwt_rank_search_and_locate() {
    let texts = enumerate_strings(&CANONICAL, 5);
    let patterns = enumerate_strings(&CANONICAL, 4);
    assert_eq!(texts.len(), 1_365);
    assert_eq!(patterns.len(), 341);
    assert_enumeration_complete(&texts, &CANONICAL, 5);
    assert_enumeration_complete(&patterns, &CANONICAL, 4);

    let mut search_cases = 0_u64;
    let mut rank_cases = 0_u64;
    let mut interval_cases = 0_u64;
    let mut substring_cases = 0_u64;
    for text in &texts {
        let source = to_search_bases(text);
        let original_source = source.clone();
        let index = FmIndex::build_reference(&source, FmBuildLimit::MAX)
            .expect("exhaustive tiny index should build");
        assert_eq!(source, original_source);

        let suffix_array = naive_suffix_array(text);
        let (text_rank_cases, text_interval_cases) =
            verify_text_representation(text, &index, &suffix_array);
        rank_cases += text_rank_cases;
        interval_cases += text_interval_cases;

        for pattern in &patterns {
            verify_search_case(text, &index, &suffix_array, pattern);
            search_cases += 1;
        }
        substring_cases += verify_all_substrings(text, &index, &suffix_array);
    }

    assert_eq!(search_cases, 465_465);
    assert_eq!(rank_cases, 36_408);
    assert_eq!(interval_cases, 35_195);
    assert_eq!(substring_cases, 18_356);
}

fn verify_projected_alphabet(alphabet: &[u8], missing_base: SearchBase, expected_cases: u64) {
    let texts = enumerate_strings(alphabet, 6);
    let patterns = enumerate_strings(alphabet, 4);
    assert_eq!(texts.len(), 1_093);
    assert_eq!(patterns.len(), 121);
    assert_enumeration_complete(&texts, alphabet, 6);
    assert_enumeration_complete(&patterns, alphabet, 4);

    let mut cases = 0_u64;
    let mut substring_cases = 0_u64;
    for text in &texts {
        let index = build(text);
        assert_eq!(
            index
                .rank(missing_base, index.suffix_count())
                .expect("complete rank boundary is valid"),
            0
        );
        let suffix_array = naive_suffix_array(text);
        for pattern in &patterns {
            verify_search_case(text, &index, &suffix_array, pattern);
            cases += 1;
        }
        substring_cases += verify_all_substrings(text, &index, &suffix_array);
    }
    assert_eq!(cases, expected_cases);
    assert_eq!(substring_cases, 19_956);
}

#[test]
fn projected_ct_and_ga_alphabets_match_independent_oracles() {
    verify_projected_alphabet(&CT_PROJECTED, SearchBase::C, 132_253);
    verify_projected_alphabet(&GA_PROJECTED, SearchBase::G, 132_253);
}

fn query_snapshot(index: &FmIndex) -> Vec<u64> {
    let mut snapshot = Vec::new();
    for boundary in 0..=index.suffix_count() {
        for base in SearchBase::ALL {
            snapshot.push(
                index
                    .rank(base, boundary)
                    .expect("snapshot rank boundary is valid"),
            );
        }
    }
    for pattern in [b"".as_slice(), b"A", b"AC", b"GT", b"TTTT"] {
        let interval = index.exact_search(&to_search_bases(pattern));
        snapshot.push(interval.lower());
        snapshot.push(interval.upper());
        snapshot.extend(located_values(index, interval));
    }
    snapshot
}

fn assert_send_sync<T: Send + Sync>() {}

#[test]
fn ownership_rebuild_and_shared_queries_are_deterministic() {
    assert_send_sync::<FmIndex>();

    let raw = b"ACGTACACGTAAAACT";
    let mut source = to_search_bases(raw);
    let index = Arc::new(
        FmIndex::build_reference(&source, FmBuildLimit::MAX).expect("owned index should build"),
    );
    let expected = query_snapshot(&index);
    source.fill(SearchBase::T);
    assert_eq!(query_snapshot(&index), expected);

    let mut disposable = index
        .locate(index.full_interval())
        .expect("full interval should locate");
    disposable.reverse();
    disposable.clear();
    assert_eq!(query_snapshot(&index), expected);

    let rebuilt = build(raw);
    assert_eq!(query_snapshot(&rebuilt), expected);

    let mut workers = Vec::new();
    for _ in 0..8 {
        let shared = Arc::clone(&index);
        let expected = expected.clone();
        workers.push(thread::spawn(move || {
            for _ in 0..50 {
                assert_eq!(query_snapshot(&shared), expected);
            }
        }));
    }
    for worker in workers {
        worker.join().expect("query worker should not panic");
    }
    assert_eq!(query_snapshot(&index), expected);
}

fn verify_medium_text(raw: &[u8]) {
    let index = build(raw);
    let full = located_values(&index, index.full_interval());
    let suffix_array = naive_suffix_array(raw);
    let expected_order = suffix_array.iter().copied().map(to_u64).collect::<Vec<_>>();
    assert_eq!(full, expected_order);

    let mut sorted = full;
    sorted.sort_unstable();
    let expected_set = (0..=raw.len()).map(to_u64).collect::<Vec<_>>();
    assert_eq!(sorted, expected_set);
    for pattern in [b"A".as_slice(), b"ACGT", b"AAAA", b"TACG", b"CCCC"] {
        verify_search_case(raw, &index, &suffix_array, pattern);
    }
}

#[test]
fn medium_balanced_and_repetitive_indexes_have_complete_locate() {
    let balanced = (0..4_096)
        .map(|index| CANONICAL[index % CANONICAL.len()])
        .collect::<Vec<_>>();
    let repetitive = vec![b'A'; 4_096];
    verify_medium_text(&balanced);
    verify_medium_text(&repetitive);
}

#[test]
fn text_offset_is_a_boundary_capable_logical_value() {
    let index = build(b"AC");
    let offsets = index
        .locate(index.full_interval())
        .expect("full interval should locate");
    let values = offsets
        .iter()
        .copied()
        .map(TextOffset::get)
        .collect::<Vec<_>>();
    assert!(values.contains(&index.text_len()));
    assert_eq!(values.len(), 3);
}
