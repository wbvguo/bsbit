//! CPU capability detection and process-wide implementation selection.
//!
//! This leaf crate contains no alignment, indexing, or biological semantics.
//! It publishes one immutable backend decision that algorithm crates use to
//! construct their own private dispatch tables. CPU-specific instructions stay
//! in those algorithm crates behind explicit `target_feature` boundaries.

#![deny(unsafe_code)]

#[cfg(target_arch = "aarch64")]
mod aarch64;
mod backend;
mod features;
#[cfg(target_arch = "x86_64")]
mod x86_64;

pub use backend::{
    Backend, BackendParseError, BackendRequest, BackendUnavailable, Configuration,
    InitializationError, configuration, initialize, selected_backend,
};
pub use features::{Architecture, CpuFeatures};
