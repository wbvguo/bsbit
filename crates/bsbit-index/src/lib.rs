//! Reference catalogs, FM-index primitives, storage layouts, and builders.
//!
//! The crate separates construction from the immutable representations queried
//! at runtime.  Higher-level seed selection and candidate policy belong to
//! `bsbit-align`.

#![deny(unsafe_op_in_unsafe_fn)]
#![deny(unsafe_code)]

#[cfg(any(test, feature = "index-construction"))]
pub mod build;
pub mod reference;
// The only unsafe operations in this module are private, target-feature-
// isolated POPCNT/NEON calls installed after bsbit-cpu validation.
#[cfg(feature = "combined-index")]
#[allow(unsafe_code)]
mod simd;
pub mod storage;
