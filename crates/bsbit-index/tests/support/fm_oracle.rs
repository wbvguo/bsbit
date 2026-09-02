use core::cmp::Ordering;

use bsbit_index::storage::fm::SearchBase;

pub(crate) const CANONICAL: [u8; 4] = *b"ACGT";
pub(crate) const CT_PROJECTED: [u8; 3] = *b"AGT";
pub(crate) const GA_PROJECTED: [u8; 3] = *b"ACT";

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum OracleBwt {
    Sentinel,
    Base(u8),
}

pub(crate) fn to_search_bases(raw: &[u8]) -> Vec<SearchBase> {
    raw.iter()
        .map(|byte| match byte {
            b'A' => SearchBase::A,
            b'C' => SearchBase::C,
            b'G' => SearchBase::G,
            b'T' => SearchBase::T,
            _ => panic!("oracle input contains a noncanonical byte"),
        })
        .collect()
}

pub(crate) fn enumerate_strings(alphabet: &[u8], maximum_length: usize) -> Vec<Vec<u8>> {
    let mut all = vec![Vec::new()];
    let mut current = vec![Vec::new()];
    for _ in 0..maximum_length {
        let mut next = Vec::new();
        for prefix in &current {
            for &base in alphabet {
                let mut value = prefix.clone();
                value.push(base);
                next.push(value);
            }
        }
        all.extend(next.iter().cloned());
        current = next;
    }
    all
}

pub(crate) fn naive_suffix_array(text: &[u8]) -> Vec<usize> {
    let mut suffixes = (0..=text.len()).collect::<Vec<_>>();
    suffixes.sort_by(|left, right| compare_suffixes(text, *left, *right));
    suffixes
}

pub(crate) fn naive_bwt(text: &[u8], suffix_array: &[usize]) -> Vec<OracleBwt> {
    suffix_array
        .iter()
        .map(|&start| {
            if start == 0 {
                OracleBwt::Sentinel
            } else {
                OracleBwt::Base(text[start - 1])
            }
        })
        .collect()
}

pub(crate) fn naive_rank(bwt: &[OracleBwt], base: u8, boundary: usize) -> u64 {
    let count = bwt[..boundary]
        .iter()
        .filter(|&&symbol| symbol == OracleBwt::Base(base))
        .count();
    to_u64(count)
}

pub(crate) fn naive_interval(text: &[u8], suffix_array: &[usize], pattern: &[u8]) -> (u64, u64) {
    let mut lower = 0;
    while lower < suffix_array.len()
        && compare_suffix_pattern(text, suffix_array[lower], pattern) == Ordering::Less
    {
        lower += 1;
    }
    let mut upper = lower;
    while upper < suffix_array.len()
        && compare_suffix_pattern(text, suffix_array[upper], pattern) == Ordering::Equal
    {
        upper += 1;
    }
    (to_u64(lower), to_u64(upper))
}

pub(crate) fn direct_hits(text: &[u8], pattern: &[u8]) -> Vec<u64> {
    (0..=text.len())
        .filter(|&offset| text[offset..].starts_with(pattern))
        .map(to_u64)
        .collect()
}

pub(crate) fn to_u64(value: usize) -> u64 {
    u64::try_from(value).expect("test sizes fit u64")
}

fn raw_rank(base: u8) -> u8 {
    match base {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => panic!("oracle comparison contains a noncanonical byte"),
    }
}

fn compare_suffixes(text: &[u8], left: usize, right: usize) -> Ordering {
    let mut left_index = left;
    let mut right_index = right;
    loop {
        match (text.get(left_index), text.get(right_index)) {
            (None, None) => return Ordering::Equal,
            (None, Some(_)) => return Ordering::Less,
            (Some(_), None) => return Ordering::Greater,
            (Some(&left_base), Some(&right_base)) => {
                let order = raw_rank(left_base).cmp(&raw_rank(right_base));
                if order != Ordering::Equal {
                    return order;
                }
                left_index += 1;
                right_index += 1;
            }
        }
    }
}

fn compare_suffix_pattern(text: &[u8], start: usize, pattern: &[u8]) -> Ordering {
    for (pattern_index, &pattern_base) in pattern.iter().enumerate() {
        let Some(&text_base) = text.get(start + pattern_index) else {
            return Ordering::Less;
        };
        let order = raw_rank(text_base).cmp(&raw_rank(pattern_base));
        if order != Ordering::Equal {
            return order;
        }
    }
    Ordering::Equal
}
