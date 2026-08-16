//! Fixed bounds shared by single-read verification and paired-read selection.

/// Edit-distance budget used by the initial low-latency verification pass.
pub(crate) const INITIAL_EDIT_DISTANCE: u8 = 3;

/// Largest per-read edit-distance budget supported by the mapping core.
pub(crate) const MAX_EDIT_DISTANCE: u8 = 5;

/// Largest read accepted by the fixed verification buffers.
pub(crate) const MAX_READ_BASES: usize = 3 * 64;

/// Maximum number of candidates processed by one vectorized verifier call.
pub(crate) const VERIFICATION_BATCH: usize = 32;

/// Dense combined-index suffix width used to begin candidate discovery.
pub(crate) const MIN_SUFFIX_BASES: usize = 16;
