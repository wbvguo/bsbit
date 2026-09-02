//! Shared constants and component naming for the combined-index wire format.

use std::ffi::OsString;
use std::path::{Path, PathBuf};

pub(crate) const BWT_WORDS_PER_128_ROWS: u64 = 5;
pub(crate) const SA_FLAG_WORDS_PER_256_ROWS: u64 = 5;

pub(crate) const META_BYTES: usize = 120;
pub(crate) const META_BYTES_U32: u32 = 120;
pub(crate) const META_EXTENSION_MAGIC: &[u8; 8] = b"BSBICMB1";
pub(crate) const META_EXTENSION_MAJOR: u16 = 1;
pub(crate) const META_EXTENSION_MINOR: u16 = 0;
pub(crate) const META_EXTENSION_MINOR_SA8: u16 = 1;
pub(crate) const META_EXTENSION_OFFSET: usize = 68;
pub(crate) const META_DIGEST_OFFSET: usize = 84;

pub(crate) const OCC_STRIDE: u32 = 64;
pub(crate) const HIGH_OCC_STRIDE: u32 = 128;
pub(crate) const LOOKUP_BASES: usize = 16;
pub(crate) const LOOKUP_KEYS: u64 = 43_046_721;
pub(crate) const LOOKUP_KEYS_USIZE: usize = 43_046_721;
pub(crate) const LOOKUP_ENTRIES: u64 = LOOKUP_KEYS + 1;

pub(crate) fn suffixed_path(prefix: &Path, suffix: &str) -> PathBuf {
    let mut name: OsString = prefix.as_os_str().to_owned();
    name.push(suffix);
    PathBuf::from(name)
}
