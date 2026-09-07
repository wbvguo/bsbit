//! Baseline-safe x86-64 feature detection.

use crate::{Architecture, CpuFeatures};

pub(crate) fn detect() -> CpuFeatures {
    CpuFeatures {
        architecture: Architecture::X86_64,
        sse2: std::arch::is_x86_feature_detected!("sse2"),
        sse41: std::arch::is_x86_feature_detected!("sse4.1"),
        sse42: std::arch::is_x86_feature_detected!("sse4.2"),
        popcnt: std::arch::is_x86_feature_detected!("popcnt"),
        avx2: std::arch::is_x86_feature_detected!("avx2"),
        avx512f: std::arch::is_x86_feature_detected!("avx512f"),
        avx512bw: std::arch::is_x86_feature_detected!("avx512bw"),
        neon: false,
    }
}
