//! Immutable in-memory and memory-mapped index representations.

#[cfg(feature = "combined-index")]
// The combined image validates every mapped component extent once at open,
// then uses three tiny unchecked integer readers in its rank hot path.
#[allow(unsafe_code)]
pub mod combined;
#[cfg(feature = "combined-index")]
pub(crate) mod combined_layout;
pub mod fm;
#[cfg(feature = "combined-index")]
pub mod reference_catalog;
