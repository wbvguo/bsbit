//! Architecture capability values used by backend policy.

use core::fmt;

/// Architecture of the running process and its executable.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum Architecture {
    /// 64-bit x86.
    X86_64,
    /// 64-bit Arm, with NEON detected independently.
    Aarch64,
    /// Any other architecture supported by the Rust implementation.
    Other,
}

impl Architecture {
    /// Stable diagnostic spelling.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::X86_64 => "x86_64",
            Self::Aarch64 => "aarch64",
            Self::Other => std::env::consts::ARCH,
        }
    }
}

impl fmt::Display for Architecture {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.as_str())
    }
}

/// CPU features relevant to bsbit's supported runtime backends.
///
/// The x86 bits are intentionally independent. In particular, vector features
/// do not imply POPCNT or one another in bsbit's feature contract. NEON is
/// detected independently so custom `AArch64` targets can fall back safely.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[allow(clippy::struct_excessive_bools)]
pub struct CpuFeatures {
    pub(crate) architecture: Architecture,
    pub(crate) sse2: bool,
    pub(crate) sse41: bool,
    pub(crate) sse42: bool,
    pub(crate) popcnt: bool,
    pub(crate) avx2: bool,
    pub(crate) avx512f: bool,
    pub(crate) avx512bw: bool,
    pub(crate) neon: bool,
}

impl CpuFeatures {
    /// Creates a synthetic x86-64 capability set for deterministic policy tests.
    #[must_use]
    pub const fn new_x86_64(avx2: bool, sse41: bool, popcnt: bool) -> Self {
        Self {
            architecture: Architecture::X86_64,
            sse2: true,
            sse41,
            sse42: false,
            popcnt,
            avx2,
            avx512f: false,
            avx512bw: false,
            neon: false,
        }
    }

    /// Overrides the synthetic x86-64 SSE2 bit for deterministic fallback
    /// tests. Real x86-64 targets provide SSE2 as an architectural baseline.
    #[must_use]
    pub const fn with_sse2(mut self, sse2: bool) -> Self {
        self.sse2 = sse2;
        self
    }

    /// Adds an explicit SSE4.2 feature bit to a synthetic x86-64 capability
    /// set for deterministic policy tests.
    #[must_use]
    pub const fn with_sse42(mut self, sse42: bool) -> Self {
        self.sse42 = sse42;
        self
    }

    /// Adds explicit AVX-512 feature bits to a synthetic x86-64 capability
    /// set for deterministic policy tests.
    #[must_use]
    pub const fn with_avx512(mut self, avx512f: bool, avx512bw: bool) -> Self {
        self.avx512f = avx512f;
        self.avx512bw = avx512bw;
        self
    }

    /// Creates the standard NEON-capable `AArch64` set for policy tests.
    #[must_use]
    pub const fn new_aarch64() -> Self {
        Self::new_aarch64_with_neon(true)
    }

    /// Creates an `AArch64` capability set with an explicit NEON value for
    /// deterministic fallback-policy tests.
    #[must_use]
    pub const fn new_aarch64_with_neon(neon: bool) -> Self {
        Self {
            architecture: Architecture::Aarch64,
            sse2: false,
            sse41: false,
            sse42: false,
            popcnt: false,
            avx2: false,
            avx512f: false,
            avx512bw: false,
            neon,
        }
    }

    /// Creates a capability set for a non-x86, non-Arm architecture.
    #[must_use]
    pub const fn new_other() -> Self {
        Self {
            architecture: Architecture::Other,
            sse2: false,
            sse41: false,
            sse42: false,
            popcnt: false,
            avx2: false,
            avx512f: false,
            avx512bw: false,
            neon: false,
        }
    }

    /// Detects the current process-visible CPU features.
    ///
    /// On x86-64, hypervisor-masked CPUID bits are reflected as exposed to the
    /// process. Rust's AVX2 detector also verifies operating-system AVX state
    /// support before reporting AVX2.
    #[must_use]
    #[cfg(target_arch = "x86_64")]
    pub fn detect() -> Self {
        crate::x86_64::detect()
    }

    /// Detects process-visible `AArch64` NEON support.
    #[must_use]
    #[cfg(target_arch = "aarch64")]
    pub fn detect() -> Self {
        crate::aarch64::detect()
    }

    /// Selects the portable feature set on architectures without a specialized
    /// bsbit backend.
    #[must_use]
    #[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64")))]
    pub const fn detect() -> Self {
        Self::new_other()
    }

    /// Architecture represented by this capability set.
    #[must_use]
    pub const fn architecture(self) -> Architecture {
        self.architecture
    }

    /// Whether SSE2 is exposed to this process.
    #[must_use]
    pub const fn sse2(self) -> bool {
        self.sse2
    }

    /// Whether SSE4.1 is exposed to this process.
    #[must_use]
    pub const fn sse41(self) -> bool {
        self.sse41
    }

    /// Whether SSE4.2 is exposed to this process.
    #[must_use]
    pub const fn sse42(self) -> bool {
        self.sse42
    }

    /// Whether POPCNT is exposed to this process.
    #[must_use]
    pub const fn popcnt(self) -> bool {
        self.popcnt
    }

    /// Whether AVX2, including operating-system AVX state support, is exposed.
    #[must_use]
    pub const fn avx2(self) -> bool {
        self.avx2
    }

    /// Whether AVX-512 Foundation, including operating-system state support,
    /// is exposed to this process.
    #[must_use]
    pub const fn avx512f(self) -> bool {
        self.avx512f
    }

    /// Whether AVX-512 byte/word operations are exposed to this process.
    #[must_use]
    pub const fn avx512bw(self) -> bool {
        self.avx512bw
    }

    /// Whether `AArch64` Advanced SIMD/NEON is available.
    #[must_use]
    pub const fn neon(self) -> bool {
        self.neon
    }
}
