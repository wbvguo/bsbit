//! Baseline-safe `AArch64` feature detection.

use crate::{Architecture, CpuFeatures};

pub(crate) fn detect() -> CpuFeatures {
    CpuFeatures {
        architecture: Architecture::Aarch64,
        sse2: false,
        sse41: false,
        sse42: false,
        popcnt: false,
        avx2: false,
        avx512f: false,
        avx512bw: false,
        neon: std::arch::is_aarch64_feature_detected!("neon"),
    }
}
