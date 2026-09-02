//! Shared exact 3' adapter evidence for read-layout output policies.

use bsbit_core::alphabet::Base;

pub(crate) const ADAPTER_STABILITY_DELTA: usize = 8;
pub(crate) const ILLUMINA_UNIVERSAL_ADAPTER: &[u8] = b"AGATCGGAAGAGC";
pub(crate) const MIN_ADAPTER_RETAINED_BASES: usize = 50;
pub(crate) const MIN_ADAPTER_SUPPORT_BASES: usize = 8;
pub(crate) const THREE_PRIME_ADAPTER_MAX_CLIP_BASES: usize = 30;

pub(crate) fn sequencing_three_prime_adapter_supported(read: &[Base], retained_end: usize) -> bool {
    let clipped = read.get(retained_end..).unwrap_or_default();
    let supported = clipped.len().min(ILLUMINA_UNIVERSAL_ADAPTER.len());
    supported >= MIN_ADAPTER_SUPPORT_BASES
        && clipped
            .iter()
            .take(supported)
            .zip(ILLUMINA_UNIVERSAL_ADAPTER.iter().take(supported))
            .all(|(observed, expected)| observed.as_ascii() == *expected)
}

#[must_use]
pub(crate) fn supported_three_prime_adapter_start(read: &[Base]) -> Option<usize> {
    let earliest = read
        .len()
        .saturating_sub(THREE_PRIME_ADAPTER_MAX_CLIP_BASES);
    let latest = read.len().checked_sub(MIN_ADAPTER_SUPPORT_BASES)?;
    (earliest..=latest).find(|&start| sequencing_three_prime_adapter_supported(read, start))
}

pub(crate) fn read_has_supported_three_prime_adapter(read: &[Base]) -> bool {
    supported_three_prime_adapter_start(read).is_some()
}
