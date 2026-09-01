//! White-box tests for the combined sparse-SA builder implementation.
//!
//! Kept outside production `src/` while remaining a child module so private
//! invariants can be tested without widening the crate API.

use std::fs;
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use crate::reference::ContigInput;
use crate::storage::combined::{META_EXTENSION_MINOR, META_EXTENSION_MINOR_SA8};
use bsbit_core::sequence::normalize_dna;

use super::*;

struct TestDirectory(PathBuf);

impl TestDirectory {
    fn new(label: &str) -> Self {
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("test clock follows epoch")
            .as_nanos();
        let path = std::env::temp_dir().join(format!(
            "bsbit-combined-build-test-{}-{label}-{nonce}",
            std::process::id()
        ));
        fs::create_dir(&path).expect("create unique combined-build test directory");
        Self(path)
    }

    fn path(&self, name: &str) -> PathBuf {
        self.0.join(name)
    }
}

impl Drop for TestDirectory {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn contig(name: &[u8], bases: &[u8]) -> ContigInput {
    ContigInput::new(name.to_vec(), normalize_dna(bases).expect("valid test DNA"))
}

fn read_u64_at(bytes: &[u8], offset: usize) -> u64 {
    u64::from_le_bytes(bytes[offset..offset + 8].try_into().expect("u64 bytes"))
}

fn read_u32_at(bytes: &[u8], offset: usize) -> u32 {
    u32::from_le_bytes(bytes[offset..offset + 4].try_into().expect("u32 bytes"))
}

#[test]
fn metadata_binds_the_combined_image_to_the_reference_digest() {
    let directory = TestDirectory::new("bound-metadata");
    let path = directory.path("index");
    let file = create_new_file(&path).expect("create metadata file");
    let digest = ReferenceSemanticDigest::from_bytes([0xa5; 32]);
    write_metadata(
        &file,
        BwtDimensions {
            suffix_count: 17,
            sentinel_row: 3,
            first_occurrence: [1, 5, 9, 17],
            bwt_words: 3,
            high_occ_entries: 2,
        },
        CombinedIndexSaStride::Sixteen,
        digest,
    )
    .expect("write bound metadata");
    let bytes = fs::read(path).expect("read bound metadata");
    assert_eq!(bytes.len(), META_BYTES);
    assert_eq!(
        &bytes[META_EXTENSION_OFFSET..META_EXTENSION_OFFSET + 8],
        META_EXTENSION_MAGIC
    );
    assert_eq!(
        read_u32_at(&bytes, 80),
        u32::try_from(META_BYTES).expect("metadata length fits u32")
    );
    assert_eq!(
        &bytes[META_DIGEST_OFFSET..META_DIGEST_OFFSET + 32],
        digest.as_bytes()
    );
    assert_eq!(&bytes[116..120], &[0; 4]);
    assert_eq!(read_u32_at(&bytes, 56), 16);
    assert_eq!(
        u16::from_le_bytes(bytes[78..80].try_into().expect("minor bytes")),
        META_EXTENSION_MINOR
    );
}

#[test]
fn metadata_distinguishes_balanced_and_fast_sparse_sa_layouts() {
    let directory = TestDirectory::new("stride-metadata");
    let dimensions = BwtDimensions {
        suffix_count: 17,
        sentinel_row: 3,
        first_occurrence: [1, 5, 9, 17],
        bwt_words: 3,
        high_occ_entries: 2,
    };
    let digest = ReferenceSemanticDigest::from_bytes([0x5a; 32]);
    for (name, stride, expected_value, expected_minor) in [
        (
            "balanced",
            CombinedIndexSaStride::Sixteen,
            16,
            META_EXTENSION_MINOR,
        ),
        (
            "fast",
            CombinedIndexSaStride::Eight,
            8,
            META_EXTENSION_MINOR_SA8,
        ),
    ] {
        let path = directory.path(name);
        let file = create_new_file(&path).expect("create metadata file");
        write_metadata(&file, dimensions, stride, digest).expect("write stride metadata");
        let bytes = fs::read(path).expect("read stride metadata");
        assert_eq!(read_u32_at(&bytes, 56), expected_value);
        assert_eq!(
            u16::from_le_bytes(bytes[78..80].try_into().expect("minor bytes")),
            expected_minor
        );
    }
}

#[test]
fn builder_rejects_a_digest_that_does_not_describe_its_source_catalog() {
    let directory = TestDirectory::new("foreign-source-digest");
    let path = directory.path("index");
    let error = build_combined_index_from_catalog_create_new(
        vec![contig(b"chr1", b"ACGTN")],
        ReferenceSemanticDigest::from_bytes([0; 32]),
        &path,
        CombinedIndexBuildOptions::new(1).expect("builder options"),
    )
    .expect_err("foreign digest must fail before construction");
    assert!(matches!(
        error,
        CombinedIndexBuildError::ReferenceDigestMismatch { .. }
    ));
    assert!(!path.exists());
}

#[test]
fn projection_crosses_contig_boundaries_without_separators() {
    let projected = project_combined_text(&[contig(b"one", b"AC"), contig(b"two", b"GA")], 7)
        .expect("project reference");
    assert_eq!(projected.reference_bases, 4);
    assert_eq!(projected.replaced_unknown_bases, 0);
    assert_eq!(projected.digits, [1, 0, 1, 1, 2, 0, 1, 2]);
}

#[test]
fn unknown_projection_is_repeatable_and_bound_to_reference_identity() {
    let reference = [contig(b"unknowns", b"NNNNNNNNNNNNNNNN")];
    let first = project_combined_text(&reference, 11).expect("first projection");
    let repeated = project_combined_text(&reference, 11).expect("repeat projection");
    let parallel =
        project_combined_text_with_threads(&reference, 11, 8).expect("parallel projection");
    let different = project_combined_text(&reference, 12).expect("different projection");
    assert_eq!(first.digits, repeated.digits);
    assert_eq!(first.digits, parallel.digits);
    assert_eq!(
        first.replaced_unknown_bases,
        parallel.replaced_unknown_bases
    );
    assert_ne!(first.digits, different.digits);
    assert_eq!(first.replaced_unknown_bases, 16);
    assert!(first.digits.iter().all(|&digit| digit <= 2));
}

#[test]
fn parallel_projection_matches_scalar_across_contig_segments() {
    let contigs = (0..19_u8)
        .map(|ordinal| {
            let name = format!("contig-{ordinal}").into_bytes();
            let bases = (0..(257 + usize::from(ordinal)))
                .map(|offset| b"ACGTN"[(offset + usize::from(ordinal)) % 5])
                .collect::<Vec<_>>();
            contig(&name, &bases)
        })
        .collect::<Vec<_>>();
    let scalar = project_combined_text(&contigs, 0x1234_5678).expect("scalar projection");
    let parallel =
        project_combined_text_with_threads(&contigs, 0x1234_5678, 8).expect("parallel projection");
    assert_eq!(scalar.digits, parallel.digits);
    assert_eq!(scalar.reference_bases, parallel.reference_bases);
    assert_eq!(
        scalar.replaced_unknown_bases,
        parallel.replaced_unknown_bases
    );
    let packed =
        project_combined_packed_text(&contigs, 0x1234_5678, 8).expect("packed parallel projection");
    assert_eq!(packed.reference_bases(), scalar.reference_bases);
    assert_eq!(
        (0..packed.len())
            .map(|position| packed.get(position))
            .collect::<Vec<_>>(),
        scalar.digits
    );
}

#[test]
fn direct_32_and_64_bit_libsais_outputs_are_identical() {
    let text = (0..1_037_usize)
        .map(|offset| u8::try_from((offset * 11 + offset / 7) % 3).expect("combined digit"))
        .collect::<Vec<_>>();
    let narrow = build_direct_bwt32(text.clone(), 3, DEFAULT_COMBINED_INDEX_SA_STRIDE)
        .expect("libsais32 direct BWT");
    let wide = build_direct_bwt64(text, 3, DEFAULT_COMBINED_INDEX_SA_STRIDE)
        .expect("libsais64 direct BWT");
    assert_eq!(narrow.transformed, wide.transformed);
    assert_eq!(narrow.sentinel_row, wide.sentinel_row);
    assert_eq!(narrow.sampled_rows, wide.sampled_rows);
}

#[test]
fn bounded_backend_matches_direct_bwt_occ_and_sa16_bytes() {
    let digits = (0..2_074_usize)
        .map(|offset| u8::try_from((offset * 11 + offset / 7) % 3).expect("combined digit"))
        .collect::<Vec<_>>();
    let mut direct = build_direct_bwt(digits.clone(), 4, DEFAULT_COMBINED_INDEX_SA_STRIDE)
        .expect("direct combined-index BWT");
    let packed =
        PackedProjectedText::from_projected_digits(digits, 1_037).expect("pack projected text");
    let config = BoundedBwtConfig::new(64, 4)
        .expect("bounded config")
        .with_block_bases(127)
        .expect("test block size");
    let bounded = build_bounded_bwt(packed, config).expect("bounded combined-index BWT");
    assert_eq!(
        u64::try_from(bounded.sentinel_row()).unwrap(),
        direct.sentinel_row
    );
    assert_eq!(
        (0..bounded.text_len())
            .map(|line| bounded.transformed_digit(line))
            .collect::<Vec<_>>(),
        direct.transformed
    );

    let directory = TestDirectory::new("bounded-components");
    let direct_bwt = directory.path("direct.bwt");
    let direct_occ = directory.path("direct.occ");
    let bounded_bwt = directory.path("bounded.bwt");
    let bounded_occ = directory.path("bounded.occ");
    let direct_dimensions = write_bwt_and_occ(
        &direct.transformed,
        direct.sentinel_row,
        &direct_bwt,
        &direct_occ,
    )
    .expect("write direct BWT/Occ");
    let bounded_bwt_file = create_new_file(&bounded_bwt).expect("reserve bounded BWT");
    let bounded_occ_file = create_new_file(&bounded_occ).expect("reserve bounded Occ");
    let bounded_dimensions =
        write_bounded_bwt_and_occ(&bounded, &bounded_bwt_file, &bounded_occ_file)
            .expect("write bounded BWT/Occ");
    assert_eq!(
        direct_dimensions.suffix_count,
        bounded_dimensions.suffix_count
    );
    assert_eq!(
        fs::read(direct_bwt).unwrap(),
        fs::read(bounded_bwt).unwrap()
    );
    assert_eq!(
        fs::read(direct_occ).unwrap(),
        fs::read(bounded_occ).unwrap()
    );

    let direct_sa = directory.path("direct.sa");
    let bounded_sa = directory.path("bounded.sa");
    write_sa16(
        &mut direct.sampled_rows,
        direct_dimensions.suffix_count,
        4,
        DEFAULT_COMBINED_INDEX_SA_STRIDE,
        &direct_sa,
    )
    .expect("write direct SA16");
    let direct_audit = SaAudit::prepare(
        direct.sampled_rows.clone(),
        37,
        DEFAULT_COMBINED_INDEX_SA_STRIDE,
    )
    .expect("prepare direct audit");
    let bounded_sa_file = create_new_file(&bounded_sa).expect("reserve bounded SA16");
    write_bounded_sa16(&bounded, DEFAULT_COMBINED_INDEX_SA_STRIDE, &bounded_sa_file)
        .expect("write bounded SA16");
    let bounded_audit = SaAudit::prepare(
        bounded
            .row_ordered_samples()
            .map(|(row, quotient)| (row << SA_VALUE_BITS) | u64::from(quotient))
            .collect(),
        37,
        DEFAULT_COMBINED_INDEX_SA_STRIDE,
    )
    .expect("prepare bounded audit");
    assert_eq!(fs::read(direct_sa).unwrap(), fs::read(bounded_sa).unwrap());
    assert_eq!(bounded_audit, direct_audit);
}

#[test]
fn packed_bwt_rank_matches_naive_prefixes_at_layout_boundaries() {
    let directory = TestDirectory::new("bwt-rank");
    for length in [
        2_usize, 63, 64, 65, 127, 128, 129, 255, 256, 257, 65_535, 65_536, 65_537,
    ] {
        let text = (0..length)
            .map(|offset| u8::try_from((offset * 7 + offset / 5) % 3).expect("combined digit"))
            .collect::<Vec<_>>();
        let direct = build_direct_bwt(text.clone(), 2, DEFAULT_COMBINED_INDEX_SA_STRIDE)
            .expect("build direct BWT");
        let bwt_path = directory.path(&format!("{length}.bwt"));
        let occ_path = directory.path(&format!("{length}.occ"));
        let dimensions = write_bwt_and_occ(
            &direct.transformed,
            direct.sentinel_row,
            &bwt_path,
            &occ_path,
        )
        .expect("pack combined rank");
        assert_eq!(
            fs::metadata(&bwt_path).expect("BWT metadata").len(),
            8 + dimensions.bwt_words * 8
        );
        let bwt_file = File::open(&bwt_path).expect("open BWT descriptor");
        let occ_file = File::open(&occ_path).expect("open Occ descriptor");
        let rank = BuildRank::open(&bwt_file, &occ_file, dimensions).expect("open packed rank");
        let mut prefixes = Vec::with_capacity(direct.transformed.len() + 1);
        prefixes.push([0_u64; 3]);
        for &digit in &direct.transformed {
            let mut next = *prefixes.last().expect("prefix exists");
            next[usize::from(digit)] += 1;
            prefixes.push(next);
        }
        for boundary in 0..=dimensions.suffix_count {
            let line = boundary - u64::from(boundary > dimensions.sentinel_row);
            let counts = prefixes[usize::try_from(line).expect("test line fits usize")];
            let expected = [
                dimensions.first_occurrence[0] + counts[0],
                dimensions.first_occurrence[1] + counts[1],
                dimensions.first_occurrence[2] + counts[2],
            ];
            assert_eq!(rank.all_boundaries(boundary), Some(expected));
        }
    }
}

fn assert_sa16_tail_contract(suffix_count: u64, expected_flag_entries: u64) {
    let directory = TestDirectory::new("sa-tail");
    let sample_count = (suffix_count - 1) / u64::from(DEFAULT_COMBINED_INDEX_SA_STRIDE) + 1;
    let mut rows = (0..sample_count)
        .map(|quotient| quotient * u64::from(DEFAULT_COMBINED_INDEX_SA_STRIDE))
        .collect::<Vec<_>>();
    let sa_path = directory.path("index.sa");
    let dimensions = write_sa16(
        &mut rows,
        suffix_count,
        2,
        DEFAULT_COMBINED_INDEX_SA_STRIDE,
        &sa_path,
    )
    .expect("write SA16");
    assert_eq!(dimensions.sparse_entries, sample_count);
    assert_eq!(dimensions.flag_entries, expected_flag_entries);

    let bytes = fs::read(sa_path).expect("read SA16 image");
    assert_eq!(read_u64_at(&bytes, 0), sample_count);
    for quotient in 0..sample_count {
        assert_eq!(
            read_u32_at(
                &bytes,
                8 + usize::try_from(quotient).expect("test quotient fits usize") * 4,
            ),
            u32::try_from(quotient).expect("test quotient fits u32")
        );
    }
    let flag_header = 8 + usize::try_from(sample_count).expect("sample count fits usize") * 4;
    assert_eq!(read_u64_at(&bytes, flag_header), expected_flag_entries);
    let flags_offset = flag_header + 8;
    let flag = |ordinal: usize| read_u64_at(&bytes, flags_offset + ordinal * 8);
    assert_eq!(flag(0), 0);
    assert_eq!(flag(5), 16);
    assert_eq!(flag(10), 32);
    assert_eq!(flag(15), 48);
    assert_eq!(flag(20), 64);
    if suffix_count == 1_025 {
        assert_eq!(flag(21), 1_u64 << 63);
        assert_eq!(flag(22), 0);
    } else {
        assert_eq!(flag(21), 0);
    }
    assert_eq!(
        bytes.len(),
        flags_offset + usize::try_from(expected_flag_entries).expect("flag entries fit usize") * 8
    );
}

#[test]
fn sa16_writer_retains_combined_boundary_and_guard_words() {
    assert_sa16_tail_contract(1_024, 22);
    assert_sa16_tail_contract(1_025, 23);
}

#[test]
fn parallel_radix_sort_orders_packed_samples_by_row() {
    let length = PARALLEL_RADIX_MIN_ENTRIES + 13;
    let mut values = (0..length)
        .map(|ordinal| {
            let row = u64::try_from(length - ordinal - 1).expect("test row fits u64");
            let quotient = u64::try_from(ordinal).expect("test ordinal fits u64") & SA_VALUE_MASK;
            (row << SA_VALUE_BITS) | quotient
        })
        .collect::<Vec<_>>();
    sort_packed_samples_by_row(
        &mut values,
        u64::try_from(length - 1).expect("test maximum row fits u64"),
        8,
    )
    .expect("parallel radix sort");
    for (row, &packed) in values.iter().enumerate() {
        let row = u64::try_from(row).expect("test row fits u64");
        assert_eq!(packed >> SA_VALUE_BITS, row);
        assert_eq!(
            packed & SA_VALUE_MASK,
            (u64::try_from(length).expect("test length fits u64") - row - 1) & SA_VALUE_MASK
        );
    }
}

#[test]
fn create_only_staging_rejects_an_existing_component() {
    let directory = TestDirectory::new("create-only");
    let prefix = directory.path("index");
    fs::write(suffixed_path(&prefix, ".occ"), b"occupied").expect("write occupied component");
    let error = StagedCombinedIndex::create(Path::new(&prefix)).expect_err("must reject target");
    assert!(matches!(error, CombinedIndexBuildError::Argument(_)));
}

#[test]
fn create_only_staging_rejects_a_dangling_target_symlink() {
    use std::os::unix::fs::symlink;

    let directory = TestDirectory::new("dangling-target");
    let prefix = directory.path("index");
    symlink(directory.path("missing"), suffixed_path(&prefix, ".sa"))
        .expect("create dangling target symlink");
    let error = StagedCombinedIndex::create(&prefix).expect_err("must reject dangling target");
    assert!(matches!(error, CombinedIndexBuildError::Argument(_)));
}

#[test]
fn staging_name_does_not_inherit_a_long_target_name() {
    let directory = TestDirectory::new("long-target");
    let prefix = directory.path(&"x".repeat(220));
    let staging = StagedCombinedIndex::create(&prefix).expect("create short staging name");
    assert!(
        staging
            .stage
            .meta
            .file_name()
            .expect("staging file name")
            .len()
            < 100
    );
}

#[test]
fn unsealed_staging_drop_preserves_a_replacement_and_cleans_owned_components() {
    let directory = TestDirectory::new("unsealed-replacement");
    let prefix = directory.path("index");
    let staging = StagedCombinedIndex::create(&prefix).expect("create staging files");
    let stage_paths = staging.stage.all().map(Path::to_path_buf);
    let replaced_path = staging.stage.bwt.clone();
    fs::remove_file(&replaced_path).expect("unlink owned BWT staging name");
    fs::write(&replaced_path, b"concurrent replacement").expect("replace staging name");

    drop(staging);

    assert_eq!(
        fs::read(&replaced_path).expect("replacement survives identity-safe drop"),
        b"concurrent replacement"
    );
    for path in stage_paths {
        if path != replaced_path {
            assert!(!path.exists(), "owned staging component was cleaned");
        }
    }
}

fn populate_staged_components(staging: &StagedCombinedIndex) {
    for (ordinal, component) in [&staging.meta, &staging.bwt, &staging.sa, &staging.occ]
        .into_iter()
        .enumerate()
    {
        let mut file = reset_component_file(&component.file).expect("reset staged component");
        file.write_all(&[u8::try_from(ordinal).expect("component ordinal fits u8")])
            .expect("write staged component");
        file.sync_all().expect("sync staged component");
    }
}

#[test]
fn sealed_staging_publishes_held_descriptors_create_only() {
    let directory = TestDirectory::new("sealed-publication");
    let prefix = directory.path("index");
    let staging = StagedCombinedIndex::create(&prefix).expect("create staging paths");
    populate_staged_components(&staging);
    let stage_paths = staging.stage.all().map(std::path::Path::to_path_buf);
    let targets = staging.target.all().map(std::path::Path::to_path_buf);
    let completed = staging.seal().expect("seal staged descriptors");
    completed
        .publish()
        .expect("publish descriptor-bound components");
    for (ordinal, path) in targets.iter().enumerate() {
        assert_eq!(
            fs::read(path).expect("read published component"),
            [u8::try_from(ordinal).expect("component ordinal fits u8")]
        );
    }
    assert!(stage_paths.iter().all(|path| !path.exists()));
}

#[test]
fn sealed_staging_rejects_and_preserves_a_replacement() {
    let directory = TestDirectory::new("sealed-replacement");
    let prefix = directory.path("index");
    let staging = StagedCombinedIndex::create(&prefix).expect("create staging paths");
    populate_staged_components(&staging);
    let completed = staging.seal().expect("seal staged descriptors");
    let replaced_path = completed.stage.bwt.clone();
    fs::remove_file(&replaced_path).expect("unlink sealed BWT path");
    fs::write(&replaced_path, b"replacement").expect("replace sealed BWT path");
    let error = completed
        .publish()
        .expect_err("replacement must fail publication");
    assert!(matches!(error, CombinedIndexBuildError::Io(_)));
    assert!(
        IndexComponentPaths::from_prefix(&prefix)
            .all()
            .into_iter()
            .all(|path| !path.exists())
    );
    assert_eq!(
        fs::read(replaced_path).expect("replacement survives identity-safe cleanup"),
        b"replacement"
    );
}

#[test]
fn late_component_owner_rolls_back_earlier_links_and_cleans_all_staging() {
    let directory = TestDirectory::new("publication-rollback");
    let prefix = directory.path("index");
    let staging = StagedCombinedIndex::create(&prefix).expect("create staging files");
    populate_staged_components(&staging);
    let stage_paths = staging.stage.all().map(Path::to_path_buf);
    let targets = staging.target.all().map(Path::to_path_buf);
    let completed = staging.seal().expect("seal staged descriptors");
    fs::write(&targets[3], b"late owner").expect("occupy Occ target after staging");

    let error = completed
        .publish()
        .expect_err("late owner must fail multi-file publication");
    assert!(matches!(error, CombinedIndexBuildError::Io(_)));
    assert!(!targets[0].exists(), "metadata commit marker stays absent");
    assert!(!targets[1].exists(), "published BWT is rolled back");
    assert!(!targets[2].exists(), "published SA is rolled back");
    assert_eq!(
        fs::read(&targets[3]).expect("late owner survives"),
        b"late owner"
    );
    assert!(stage_paths.iter().all(|path| !path.exists()));
}
