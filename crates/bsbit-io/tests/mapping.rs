//! Contract tests for immutable file mappings.

use std::fs::{self, OpenOptions};
use std::io::Write;
use std::path::PathBuf;

use bsbit_io::ReadOnlyMapping;

fn temporary_path(label: &str) -> PathBuf {
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .expect("test clock follows epoch")
        .as_nanos();
    std::env::temp_dir().join(format!(
        "bsbit-read-only-mapping-{}-{label}-{nonce}",
        std::process::id()
    ))
}

#[test]
fn complete_mapping_exposes_bounded_little_endian_reads() {
    let path = temporary_path("bytes");
    let mut file = OpenOptions::new()
        .read(true)
        .write(true)
        .create_new(true)
        .open(&path)
        .expect("create fixture");
    file.write_all(&[1, 2, 3, 4, 5, 6, 7, 8])
        .expect("write fixture");
    file.sync_all().expect("sync fixture");

    let mapping = ReadOnlyMapping::map(&file).expect("map fixture");
    assert_eq!(mapping.len(), 8);
    assert!(!mapping.is_empty());
    assert_eq!(mapping.as_slice(), &[1, 2, 3, 4, 5, 6, 7, 8]);
    assert_eq!(mapping.read_u8(3), 4);
    assert_eq!(mapping.read_u32(1), 0x0504_0302);
    assert_eq!(mapping.read_u64(0), 0x0807_0605_0403_0201);

    drop(mapping);
    drop(file);
    fs::remove_file(path).expect("remove fixture");
}

#[test]
fn empty_file_is_rejected_before_mmap() {
    let path = temporary_path("empty");
    let file = OpenOptions::new()
        .read(true)
        .write(true)
        .create_new(true)
        .open(&path)
        .expect("create fixture");
    let error = ReadOnlyMapping::map(&file).expect_err("empty file must fail");
    assert_eq!(error.kind(), std::io::ErrorKind::InvalidData);

    drop(file);
    fs::remove_file(path).expect("remove fixture");
}
