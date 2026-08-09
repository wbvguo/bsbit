//! Immutable in-memory and memory-mapped index representations.

#[cfg(feature = "combined-index")]
#[allow(unsafe_code)]
pub mod combined;
pub mod fm;
#[cfg(feature = "combined-index")]
#[allow(unsafe_code)]
pub mod reference_catalog;
