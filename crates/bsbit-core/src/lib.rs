//! Stable domain values and invariants shared by bsbit crates.
//!
//! This crate deliberately contains no reference index, candidate search,
//! alignment algorithm, pairing policy, file I/O, FFI, unsafe code, or mutable
//! global state.

#![forbid(unsafe_code)]

pub mod alphabet;
pub mod bisulfite;
pub mod cigar;
pub mod coordinate;
pub mod reference;
pub mod sequence;
