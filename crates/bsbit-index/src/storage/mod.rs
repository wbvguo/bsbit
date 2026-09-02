//! Immutable in-memory and memory-mapped index representations.

#[cfg(feature = "combined-index")]
#[allow(unsafe_code)]
pub mod combined;
#[cfg(feature = "combined-index")]
pub(crate) mod combined_format;
pub mod fm;
#[cfg(feature = "combined-index")]
#[allow(unsafe_code)]
pub mod reference_catalog;
