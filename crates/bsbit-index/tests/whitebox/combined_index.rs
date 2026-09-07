//! Hermetic white-box contract tests for the combined-index runtime.
//!
//! The file is loaded as a `#[cfg(test)]` child module so private runtime
//! invariants remain testable without widening the production API. One tiny
//! index is constructed per suite and reused across all assertions; this keeps
//! the dense 16-mer table out of Git while making CI independent of local data.

use super::*;
use crate::build::combined::{
    CombinedIndexBuildOptions, build_combined_index_from_catalog_create_new,
};
use crate::reference::ContigInput;
use bsbit_core::reference::ReferenceSemanticDigestBuilder;
use bsbit_core::sequence::normalize_dna;

struct TestDirectory(PathBuf);

impl TestDirectory {
    fn new(label: &str) -> Self {
        Self(unique_test_directory(label))
    }
}

impl Drop for TestDirectory {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

struct TinyFixture {
    _directory: TestDirectory,
    prefix: PathBuf,
    reference: Vec<u8>,
}

impl TinyFixture {
    fn build() -> Self {
        let directory = TestDirectory::new("hermetic");
        let prefix = directory.0.join("genome.index.bs.index");
        let mut state = 0x6a09_e667_f3bc_c909_u64;
        let reference = (0..257)
            .map(|_| {
                state = state
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);
                b"ACGT"[usize::try_from((state >> 32) & 3).expect("two bits fit usize")]
            })
            .collect::<Vec<_>>();
        let sequence = normalize_dna(&reference).expect("canonical tiny reference");
        let contig = ContigInput::new(b"tiny".to_vec(), sequence);
        let mut digest = ReferenceSemanticDigestBuilder::new(1);
        digest
            .push_normalized_contig(contig.name(), contig.sequence().bases())
            .expect("tiny reference contributes to its semantic digest");
        let digest = digest.finish().expect("tiny semantic digest completes");
        let options = CombinedIndexBuildOptions::new(1)
            .expect("one build thread")
            .with_memory_mib(512)
            .expect("bounded test memory budget");
        build_combined_index_from_catalog_create_new(vec![contig], digest, &prefix, options)
            .expect("hermetic tiny combined index builds");
        Self {
            _directory: directory,
            prefix,
            reference,
        }
    }
}

fn unique_test_directory(label: &str) -> PathBuf {
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .expect("system time follows epoch")
        .as_nanos();
    let directory = std::env::temp_dir().join(format!(
        "bsbit-combined-index-{label}-{}-{nonce}",
        std::process::id()
    ));
    std::fs::create_dir(&directory).expect("unique test directory is created");
    directory
}

fn copied_tiny_prefix(fixture: &TinyFixture, label: &str) -> (PathBuf, PathBuf) {
    let directory = unique_test_directory(label);
    let target = directory.join("genome.index.bs.index");
    for suffix in ["", ".bwt", ".sa", ".occ"] {
        std::fs::copy(
            suffixed_path(&fixture.prefix, suffix),
            suffixed_path(&target, suffix),
        )
        .expect("tiny combined-index component is copied");
    }
    (directory, target)
}

fn overwrite_u64(path: &Path, offset: u64, value: u64) {
    use std::io::{Seek, SeekFrom, Write};

    let mut file = std::fs::OpenOptions::new()
        .write(true)
        .open(path)
        .expect("copied combined-index component opens for mutation");
    file.seek(SeekFrom::Start(offset))
        .expect("combined-index field offset is reachable");
    file.write_all(&value.to_le_bytes())
        .expect("mutated combined-index value is written");
    file.sync_all()
        .expect("mutated combined-index component syncs");
}

fn search_bases(pattern: &[u8]) -> Vec<SearchBase> {
    pattern
        .iter()
        .map(|base| match base {
            b'A' => SearchBase::A,
            b'G' => SearchBase::G,
            b'T' => SearchBase::T,
            _ => panic!("projected tiny pattern"),
        })
        .collect()
}

fn public_backward_extend_rejects_foreign_domains_and_unsupported_symbols(fixture: &TinyFixture) {
    let index = CombinedIndex::open(&fixture.prefix).expect("tiny combined index opens");
    let foreign = FmInterval::private_checked(0, 1, index.suffix_count() + 1)
        .expect("foreign interval is internally valid");
    assert_eq!(index.backward_extend(foreign, SearchBase::A), None);

    let local =
        FmInterval::private_checked(0, 1, index.suffix_count()).expect("local interval is valid");
    assert_eq!(index.backward_extend(local, SearchBase::C), None);
}

fn combined_index_open_rejects_cumulative_counts_beyond_the_suffix_domain(fixture: &TinyFixture) {
    let suffix_count = CombinedIndex::open(&fixture.prefix)
        .expect("tiny combined index opens")
        .suffix_count();
    let (directory, target) = copied_tiny_prefix(fixture, "metadata-domain");
    overwrite_u64(&target, 40, suffix_count + 1);

    assert!(matches!(
        CombinedIndex::open(&target),
        Err(CombinedIndexError::Structure(
            "metadata suffix or cumulative-count domain is invalid"
        ))
    ));
    std::fs::remove_dir_all(directory).expect("unique test directory is removed");
}

fn combined_index_open_rejects_high_occurrence_values_beyond_the_suffix_domain(
    fixture: &TinyFixture,
) {
    let suffix_count = CombinedIndex::open(&fixture.prefix)
        .expect("tiny combined index opens")
        .suffix_count();
    let (directory, target) = copied_tiny_prefix(fixture, "occ-domain");
    overwrite_u64(&suffixed_path(&target, ".occ"), 8, suffix_count + 1);

    assert!(matches!(
        CombinedIndex::open(&target),
        Err(CombinedIndexError::Structure(
            "high-occurrence value exceeds suffix domain"
        ))
    ));
    std::fs::remove_dir_all(directory).expect("unique test directory is removed");
}

fn wavefront_boundary_rank_matches_scalar_for_every_batch_lane(fixture: &TinyFixture) {
    let index = CombinedIndex::open(&fixture.prefix).expect("tiny combined index opens");
    let suffix_count = index.suffix_count();
    for round in 0_u64..257 {
        let lane_count = usize::try_from(round % MAX_WAVEFRONT_LANES as u64 + 1).unwrap();
        let mut intervals = [None; MAX_WAVEFRONT_LANES];
        let mut digits = [0_u8; MAX_WAVEFRONT_LANES];
        for lane in 0..lane_count {
            let lane = u64::try_from(lane).unwrap();
            let lower = (round * 131 + lane * 67) % (suffix_count + 1);
            let width =
                [0_u64, 1, 2, 7, 31, 63, 64, 65][usize::try_from((round + lane) % 8).unwrap()];
            let upper = lower.saturating_add(width).min(suffix_count);
            intervals[usize::try_from(lane).unwrap()] = Some(
                FmInterval::private_checked(lower, upper, suffix_count)
                    .expect("sampled wavefront interval is valid"),
            );
            digits[usize::try_from(lane).unwrap()] = u8::try_from((round + lane) % 3).unwrap();
        }
        let mut observed = [None; MAX_WAVEFRONT_LANES];
        let mut workspace = BackwardExtendWavefrontWorkspace::default();
        index.backward_extend_wavefront_validated_with_workspace(
            &intervals,
            &digits,
            lane_count,
            &mut observed,
            &mut workspace,
        );
        for lane in 0..lane_count {
            let interval = intervals[lane].unwrap();
            assert_eq!(
                observed[lane],
                index.backward_extend_validated(interval, digits[lane]),
                "round={round}, lane={lane}"
            );
        }
    }
}

fn two_lane_complete_direct_locate_matches_two_scalar_intervals(fixture: &TinyFixture) {
    let index = CombinedIndex::open(&fixture.prefix).expect("tiny combined index opens");
    for first_row in 0..index.suffix_count() {
        for second_row in 0..index.suffix_count() {
            assert_eq!(
                index
                    .locate_rows_two_lanes([first_row, second_row])
                    .expect("two-lane rows locate"),
                [
                    index.locate_row(first_row).expect("first row locates"),
                    index.locate_row(second_row).expect("second row locates"),
                ],
            );
        }
    }

    let starts = [
        0_u64,
        1,
        index.sentinel_row.saturating_sub(1),
        index.sentinel_row,
        index.sentinel_row.saturating_add(1),
        index.suffix_count().saturating_sub(15),
        index.suffix_count().saturating_sub(31),
    ];
    for &first_start in &starts {
        for &second_start in starts.iter().rev() {
            for first_width in [0_u64, 1, 2, 7, 15, 17, 31] {
                for second_width in [0_u64, 1, 3, 8, 15, 17, 31] {
                    let first_upper = first_start
                        .saturating_add(first_width)
                        .min(index.suffix_count());
                    let second_upper = second_start
                        .saturating_add(second_width)
                        .min(index.suffix_count());
                    let intervals = [
                        FmInterval::private_checked(
                            first_start.min(first_upper),
                            first_upper,
                            index.suffix_count(),
                        )
                        .expect("first interval"),
                        FmInterval::private_checked(
                            second_start.min(second_upper),
                            second_upper,
                            index.suffix_count(),
                        )
                        .expect("second interval"),
                    ];
                    let mut expected = [Vec::new(), Vec::new()];
                    let expected_metrics = [
                        index
                            .visit_interval(intervals[0], &mut |position| {
                                expected[0].push(position);
                                true
                            })
                            .expect("first scalar interval locates"),
                        index
                            .visit_interval(intervals[1], &mut |position| {
                                expected[1].push(position);
                                true
                            })
                            .expect("second scalar interval locates"),
                    ];
                    let mut observed = [Vec::new(), Vec::new()];
                    let observed_metrics = index
                        .visit_interval_two_lanes_complete(intervals, &mut |lane, position| {
                            observed[lane].push(position);
                        })
                        .expect("two-lane complete intervals locate");
                    assert_eq!(observed, expected);
                    assert_eq!(observed_metrics, expected_metrics);
                }
            }
        }
    }
}

fn combined_index_lookup_rank_and_sa16_match_tiny_directional_text(fixture: &TinyFixture) {
    let index = CombinedIndex::open(&fixture.prefix).expect("tiny combined index opens");
    let reference = &fixture.reference;
    assert_eq!(
        index.reference_length(),
        u64::try_from(reference.len()).unwrap()
    );
    let mut directional = Vec::with_capacity(reference.len() * 2);
    directional.extend(reference.iter().map(|base| match base {
        b'A' | b'G' => b'T',
        b'C' => b'G',
        b'T' => b'A',
        _ => panic!("canonical tiny reference"),
    }));
    directional.extend(reference.iter().rev().map(|base| match base {
        b'C' => b'T',
        base => *base,
    }));

    for start in (0..directional.len() - 24).step_by(37) {
        let pattern = search_bases(&directional[start..start + 24]);
        let interval = index
            .exact_search(&pattern)
            .expect("combined-index exact search");
        let mut observed = Vec::new();
        let metrics = index
            .visit_interval(interval, &mut |position| {
                observed.push(usize::try_from(position).expect("tiny position fits usize"));
                true
            })
            .expect("combined-index interval locates");
        let mut expected = directional
            .windows(pattern.len())
            .enumerate()
            .filter_map(|(position, window)| {
                (window == &directional[start..start + pattern.len()]).then_some(position)
            })
            .collect::<Vec<_>>();
        observed.sort_unstable();
        expected.sort_unstable();
        assert_eq!(observed, expected, "pattern start {start}");
        assert_eq!(metrics.located_rows, interval.len());
        assert!(metrics.lf_steps < interval.len() * SA_STRIDE);
    }
}

fn boundary_pair_matches_two_scalar_ranks(fixture: &TinyFixture) {
    let index = CombinedIndex::open(&fixture.prefix).expect("tiny combined index opens");
    for lower in 0..=index.suffix_count() {
        for width in [0_u64, 1, 2, 7, 31, 63, 64, 65, 127] {
            let upper = lower.saturating_add(width).min(index.suffix_count());
            for digit in 0..3 {
                assert_eq!(
                    index
                        .lf_boundary_pair(lower, upper, digit)
                        .expect("fused boundary pair"),
                    [
                        index.lf_boundary(lower, digit).expect("scalar lower"),
                        index.lf_boundary(upper, digit).expect("scalar upper"),
                    ],
                    "lower={lower}, upper={upper}, digit={digit}"
                );
            }
        }
    }
}

fn two_lane_backward_extension_matches_two_scalar_extensions(fixture: &TinyFixture) {
    let index = CombinedIndex::open(&fixture.prefix).expect("tiny combined index opens");
    let suffix_count = index.suffix_count();
    for lower in (0..=suffix_count).step_by(7) {
        for width in [0_u64, 1, 2, 7, 31, 63, 64, 65, 127] {
            let upper = lower.saturating_add(width).min(suffix_count);
            let second_lower = suffix_count.saturating_sub(upper);
            let second_upper = second_lower.saturating_add(width).min(suffix_count);
            let first = FmInterval::private_checked(lower, upper, suffix_count)
                .expect("first two-lane interval");
            let second = FmInterval::private_checked(second_lower, second_upper, suffix_count)
                .expect("second two-lane interval");
            for first_digit in 0..3 {
                for second_digit in 0..3 {
                    let expected = [
                        index.backward_extend_validated(first, first_digit),
                        index.backward_extend_validated(second, second_digit),
                    ];
                    let mut observed = [first, second];
                    index
                        .backward_extend_interval_round(
                            &[first, second],
                            &[first_digit, second_digit],
                            &mut observed,
                        )
                        .expect("batched backward extension succeeds");
                    assert_eq!(
                        observed.map(Some),
                        expected,
                        "first={lower}..{upper}, second={second_lower}..{second_upper}, digits={first_digit}/{second_digit}"
                    );
                }
            }
        }
    }
}

#[test]
fn hermetic_combined_index_runtime_contract_suite() {
    let fixture = TinyFixture::build();
    public_backward_extend_rejects_foreign_domains_and_unsupported_symbols(&fixture);
    combined_index_open_rejects_cumulative_counts_beyond_the_suffix_domain(&fixture);
    combined_index_open_rejects_high_occurrence_values_beyond_the_suffix_domain(&fixture);
    wavefront_boundary_rank_matches_scalar_for_every_batch_lane(&fixture);
    two_lane_complete_direct_locate_matches_two_scalar_intervals(&fixture);
    combined_index_lookup_rank_and_sa16_match_tiny_directional_text(&fixture);
    boundary_pair_matches_two_scalar_ranks(&fixture);
    two_lane_backward_extension_matches_two_scalar_extensions(&fixture);
}
