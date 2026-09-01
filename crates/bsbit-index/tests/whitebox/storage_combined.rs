//! Hermetic white-box tests for combined-index layout invariants.
//!
//! This file is loaded as a `#[cfg(test)]` child module so it can verify
//! private format arithmetic without widening the production API.

use super::*;

struct TestDirectory(std::path::PathBuf);

impl TestDirectory {
    fn new(label: &str) -> Self {
        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .expect("test clock follows epoch")
            .as_nanos();
        let path = std::env::temp_dir().join(format!(
            "bsbit-combined-storage-test-{}-{label}-{nonce}",
            std::process::id()
        ));
        std::fs::create_dir(&path).expect("create unique storage test directory");
        Self(path)
    }

    fn prefix(&self) -> std::path::PathBuf {
        self.0.join("index")
    }
}

impl Drop for TestDirectory {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

fn put_u16(bytes: &mut [u8], offset: usize, value: u16) {
    bytes[offset..offset + 2].copy_from_slice(&value.to_le_bytes());
}

fn put_u32(bytes: &mut [u8], offset: usize, value: u32) {
    bytes[offset..offset + 4].copy_from_slice(&value.to_le_bytes());
}

fn put_u64(bytes: &mut [u8], offset: usize, value: u64) {
    bytes[offset..offset + 8].copy_from_slice(&value.to_le_bytes());
}

#[test]
fn reader_rejects_retired_unbound_and_mismatched_stride_metadata() {
    let directory = TestDirectory::new("retired-metadata");
    let prefix = directory.prefix();

    std::fs::write(&prefix, [0_u8; META_EXTENSION_OFFSET]).expect("write retired unbound metadata");
    assert!(matches!(
        CombinedIndex::open(&prefix),
        Err(CombinedIndexError::Structure(
            "metadata file is not the current 120-byte bound format"
        ))
    ));

    let mut metadata = [0_u8; META_BYTES];
    put_u64(&mut metadata, 0, 17);
    put_u64(&mut metadata, 8, 0);
    for (ordinal, value) in [1_u64, 1, 1, 17].into_iter().enumerate() {
        put_u64(&mut metadata, 16 + ordinal * 8, value);
    }
    put_u64(&mut metadata, 48, 17);
    put_u32(&mut metadata, 56, 8);
    put_u32(&mut metadata, 60, 64);
    put_u32(&mut metadata, 64, 128);
    metadata[META_EXTENSION_OFFSET..META_EXTENSION_OFFSET + META_EXTENSION_MAGIC.len()]
        .copy_from_slice(META_EXTENSION_MAGIC);
    put_u16(&mut metadata, 76, META_EXTENSION_MAJOR);
    put_u16(&mut metadata, 78, META_EXTENSION_MINOR);
    put_u32(&mut metadata, 80, META_BYTES_U32);
    std::fs::write(&prefix, metadata).expect("write retired SA8 metadata");
    assert!(matches!(
        CombinedIndex::open(&prefix),
        Err(CombinedIndexError::Structure(
            "metadata sparse-SA stride and format minor are unsupported"
        ))
    ));
}

#[test]
fn sparse_sa_stride_is_bound_to_its_format_minor() {
    assert_eq!(
        CombinedIndexSaStride::from_metadata(16, META_EXTENSION_MINOR),
        Some(CombinedIndexSaStride::Sixteen)
    );
    assert_eq!(
        CombinedIndexSaStride::from_metadata(8, META_EXTENSION_MINOR_SA8),
        Some(CombinedIndexSaStride::Eight)
    );
    assert_eq!(
        CombinedIndexSaStride::from_metadata(8, META_EXTENSION_MINOR),
        None
    );
    assert_eq!(
        CombinedIndexSaStride::from_metadata(16, META_EXTENSION_MINOR_SA8),
        None
    );
    assert_eq!(CombinedIndexSaStride::from_metadata(4, 2), None);
}

#[test]
fn validated_runtime_dimensions_cover_every_hot_access() {
    for suffix_count in 3_u64..=2049 {
        let text_rows = suffix_count - 1;
        let minimum_bwt = minimum_bwt_words_for_suffix_count(suffix_count)
            .expect("small suffix domain has BWT dimensions");
        let minimum_flags = minimum_sa_flag_entries_for_suffix_count(suffix_count)
            .expect("small suffix domain has flag dimensions");
        let minimum_high_occ = minimum_high_occ_entries_for_suffix_count(suffix_count)
            .expect("small suffix domain has Occ dimensions");

        for line in 0..=text_rows {
            let high_word = (line >> 7) * BWT_WORDS_PER_128_ROWS;
            assert!(high_word < minimum_bwt, "counter line {line}");
            if line < text_rows || line & 63 != 0 {
                let low_block = (line & 127) >> 6;
                let plane_start = high_word + 1 + (low_block << 1);
                assert!(plane_start + 1 < minimum_bwt, "planes line {line}");
            }
            let high_occ_block = (line >> 16) * 2;
            assert!(
                high_occ_block + 1 < minimum_high_occ,
                "high Occ line {line}"
            );
        }

        for row in 0..suffix_count {
            let block_word = (row >> 8) * SA_FLAG_WORDS_PER_256_ROWS;
            let flag_word = block_word + 1 + ((row & 255) >> 6);
            assert!(block_word < minimum_flags, "flag counter row {row}");
            assert!(flag_word < minimum_flags, "flag word row {row}");
        }
    }
}

#[test]
fn dense_lookup_key_domain_covers_every_projected_16_mer_boundary() {
    assert_eq!(
        3_u64.pow(u32::try_from(LOOKUP_BASES).unwrap()) + 1,
        LOOKUP_ENTRIES
    );
    let maximum_key = (0..LOOKUP_BASES).fold(0_u64, |key, _| key * 3 + 2);
    assert_eq!(maximum_key + 1, LOOKUP_ENTRIES - 1);
    assert!(usize::try_from(LOOKUP_ENTRIES).is_ok());
}
