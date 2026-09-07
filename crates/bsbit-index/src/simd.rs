//! Combined-index population-count implementations.
//!
//! Process-wide CPU detection and backend policy belong to `bsbit-cpu`. This
//! module owns only the index algorithm's backend-specific function table.

#![deny(unsafe_op_in_unsafe_fn)]

use bsbit_cpu::{Backend, Configuration};

/// Backend-bound 64-bit population count used by the combined FM index.
///
/// The implementation is selected once when an index is opened. The POPCNT
/// trampoline is never installed until `bsbit-cpu` has validated the complete
/// backend feature contract.
#[derive(Clone, Copy, Debug)]
pub(crate) struct PopcountDispatch(fn(u64) -> u32);

impl PopcountDispatch {
    /// Creates a function table from a process configuration that has already
    /// validated every backend feature requirement.
    pub(crate) fn from_configuration(configuration: Configuration) -> Self {
        Self::for_backend(configuration.backend())
    }

    fn for_backend(backend: Backend) -> Self {
        match backend {
            Backend::Scalar | Backend::Sse2 => Self(scalar_popcount),
            Backend::Sse42Popcnt | Backend::Avx2Popcnt | Backend::Avx512BwPopcnt => {
                Self(popcnt_trampoline)
            }
            Backend::Neon => Self(neon_popcount_trampoline),
        }
    }

    #[inline]
    pub(crate) fn count(self, value: u64) -> u32 {
        (self.0)(value)
    }
}

#[inline]
fn scalar_popcount(value: u64) -> u32 {
    value.count_ones()
}

#[cfg(target_arch = "x86_64")]
fn popcnt_trampoline(value: u64) -> u32 {
    // SAFETY: this private trampoline is installed only for a backend whose
    // bsbit-cpu validation proved POPCNT support for the current process.
    unsafe { popcnt_u64(value) }
}

#[cfg(not(target_arch = "x86_64"))]
fn popcnt_trampoline(value: u64) -> u32 {
    scalar_popcount(value)
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "popcnt")]
unsafe fn popcnt_u64(value: u64) -> u32 {
    value.count_ones()
}

#[cfg(target_arch = "aarch64")]
fn neon_popcount_trampoline(value: u64) -> u32 {
    // SAFETY: the validated NEON backend guarantees Advanced SIMD before this
    // private trampoline is installed.
    unsafe { neon_popcount_u64(value) }
}

#[cfg(not(target_arch = "aarch64"))]
fn neon_popcount_trampoline(value: u64) -> u32 {
    scalar_popcount(value)
}

#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn neon_popcount_u64(value: u64) -> u32 {
    use std::arch::aarch64::{vaddv_u8, vcnt_u8, vcreate_u64, vreinterpret_u8_u64};

    u32::from(vaddv_u8(vcnt_u8(vreinterpret_u8_u64(vcreate_u64(value)))))
}

#[cfg(test)]
mod tests {
    use super::*;
    use bsbit_cpu::CpuFeatures;

    #[test]
    fn installed_implementations_are_identical_on_supported_hardware() {
        let features = CpuFeatures::detect();
        let mut values = vec![
            0,
            1,
            u64::MAX,
            0x0123_4567_89ab_cdef,
            0xaaaa_aaaa_aaaa_aaaa,
            0x5555_5555_5555_5555,
            1_u64 << 63,
        ];
        let mut packed = 0xd1b5_4a32_d192_ed03_u64;
        for _ in 0..4_096 {
            packed = packed
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            values.push(packed);
        }
        for backend in [
            Backend::Scalar,
            Backend::Sse2,
            Backend::Sse42Popcnt,
            Backend::Avx2Popcnt,
            Backend::Avx512BwPopcnt,
            Backend::Neon,
        ] {
            if !backend.is_supported_by(features) {
                continue;
            }
            let implementation = PopcountDispatch::for_backend(backend);
            for &value in &values {
                assert_eq!(implementation.count(value), value.count_ones());
            }
        }
    }

    #[test]
    fn scalar_popcount_covers_machine_word_boundaries() {
        let cases = [
            (0, 0),
            (1, 1),
            (u64::MAX, 64),
            (1_u64 << 63, 1),
            (0x8000_0000_0000_0001, 2),
            (0xaaaa_aaaa_aaaa_aaaa, 32),
        ];
        for (value, expected) in cases {
            assert_eq!(scalar_popcount(value), expected);
        }
    }
}
