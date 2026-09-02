//! Index-driven seed discovery and candidate construction for read mapping.
//!
//! Search owns alignment policy that queries `bsbit-index`. Sequence-only
//! scoring and verification algorithms remain in `crate::verification`.

pub mod candidate;
pub mod fixed_seed;
pub mod seed;

pub(crate) mod combined_adaptive;
pub(crate) mod combined_query;
