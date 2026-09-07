//! Process-wide backend selection policy.

use core::fmt;
use core::str::FromStr;
use std::sync::OnceLock;

use crate::{Architecture, CpuFeatures};

static CONFIGURATION: OnceLock<Configuration> = OnceLock::new();
static DETECTED_CPU_FEATURES: OnceLock<CpuFeatures> = OnceLock::new();

/// Concrete implementation contract selected for the current process.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum Backend {
    /// Portable scalar kernels with no explicit SIMD or CPU-feature requirement.
    Scalar,
    /// 128-bit SSE2 vector kernels without an explicit POPCNT requirement.
    Sse2,
    /// SSE4.2-qualified vector kernels plus hardware POPCNT.
    Sse42Popcnt,
    /// AVX2 vector kernels plus hardware POPCNT.
    Avx2Popcnt,
    /// AVX-512 byte/word kernels plus the AVX2 helpers and hardware POPCNT.
    Avx512BwPopcnt,
    /// `AArch64` Advanced SIMD/NEON kernels.
    Neon,
}

impl Backend {
    /// Complete instruction-set requirement used by diagnostics and provenance.
    #[must_use]
    pub const fn instruction_set(self) -> &'static str {
        match self {
            Self::Scalar => "portable",
            Self::Sse2 => "sse2",
            Self::Sse42Popcnt => "sse4.2+popcnt",
            Self::Avx2Popcnt => "avx2+popcnt",
            Self::Avx512BwPopcnt => "avx512f+avx512bw+avx2+popcnt",
            Self::Neon => "neon",
        }
    }

    /// Whether every independent feature required by this backend is present.
    #[must_use]
    pub const fn is_supported_by(self, features: CpuFeatures) -> bool {
        match self {
            Self::Scalar => true,
            Self::Sse2 => {
                matches!(features.architecture(), Architecture::X86_64) && features.sse2()
            }
            Self::Sse42Popcnt => {
                matches!(features.architecture(), Architecture::X86_64)
                    && features.sse41()
                    && features.sse42()
                    && features.popcnt()
            }
            Self::Avx2Popcnt => {
                matches!(features.architecture(), Architecture::X86_64)
                    && features.avx2()
                    && features.popcnt()
            }
            Self::Avx512BwPopcnt => {
                matches!(features.architecture(), Architecture::X86_64)
                    && features.avx512f()
                    && features.avx512bw()
                    && features.avx2()
                    && features.popcnt()
            }
            Self::Neon => {
                matches!(features.architecture(), Architecture::Aarch64) && features.neon()
            }
        }
    }

    /// Validates this backend against the process-visible CPU.
    ///
    /// # Errors
    ///
    /// Returns the exact missing feature set instead of allowing a caller to
    /// enter an unsupported target-feature function.
    pub fn validate_current(self) -> Result<Self, BackendUnavailable> {
        BackendRequest::from(self).select(detected_features())
    }
}

impl fmt::Display for Backend {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::Scalar => "scalar",
            Self::Sse2 => "sse2",
            Self::Sse42Popcnt => "sse4.2",
            Self::Avx2Popcnt => "avx2",
            Self::Avx512BwPopcnt => "avx512",
            Self::Neon => "neon",
        })
    }
}

/// Process-visible CPU features and the immutable selected backend.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Configuration {
    features: CpuFeatures,
    backend: Backend,
}

impl Configuration {
    /// Process-visible feature set used for selection.
    #[must_use]
    pub const fn features(self) -> CpuFeatures {
        self.features
    }

    /// Concrete backend selected for this process.
    #[must_use]
    pub const fn backend(self) -> Backend {
        self.backend
    }
}

/// Requested runtime backend, including automatic selection.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum BackendRequest {
    /// Best backend for the current architecture.
    #[default]
    Auto,
    /// Force the portable scalar backend.
    Scalar,
    /// Force the 128-bit SSE2 backend.
    Sse2,
    /// Force SSE4.2+POPCNT.
    Sse42Popcnt,
    /// Force AVX2+POPCNT.
    Avx2Popcnt,
    /// Force AVX-512F+AVX-512BW with the AVX2 helpers and POPCNT.
    Avx512BwPopcnt,
    /// Force the `AArch64` NEON backend.
    Neon,
}

impl BackendRequest {
    /// Resolves this request against an explicit feature set.
    ///
    /// # Errors
    ///
    /// A forced backend fails when any independent requirement is absent.
    /// Automatic selection follows the architecture-specific baseline policy.
    pub const fn select(self, features: CpuFeatures) -> Result<Backend, BackendUnavailable> {
        if matches!(self, Self::Auto) {
            return Ok(automatic_backend(features));
        }
        let selected = match self {
            Self::Auto | Self::Scalar => Backend::Scalar,
            Self::Sse2 => Backend::Sse2,
            Self::Sse42Popcnt => Backend::Sse42Popcnt,
            Self::Avx2Popcnt => Backend::Avx2Popcnt,
            Self::Avx512BwPopcnt => Backend::Avx512BwPopcnt,
            Self::Neon => Backend::Neon,
        };
        if selected.is_supported_by(features) {
            Ok(selected)
        } else {
            Err(BackendUnavailable {
                requested: selected,
                features,
            })
        }
    }
}

const fn automatic_backend(features: CpuFeatures) -> Backend {
    match features.architecture() {
        Architecture::X86_64 => {
            if Backend::Avx512BwPopcnt.is_supported_by(features) {
                Backend::Avx512BwPopcnt
            } else if Backend::Avx2Popcnt.is_supported_by(features) {
                Backend::Avx2Popcnt
            } else if Backend::Sse42Popcnt.is_supported_by(features) {
                Backend::Sse42Popcnt
            } else if Backend::Sse2.is_supported_by(features) {
                Backend::Sse2
            } else {
                Backend::Scalar
            }
        }
        Architecture::Aarch64 => {
            if Backend::Neon.is_supported_by(features) {
                Backend::Neon
            } else {
                Backend::Scalar
            }
        }
        Architecture::Other => Backend::Scalar,
    }
}

impl From<Backend> for BackendRequest {
    fn from(backend: Backend) -> Self {
        match backend {
            Backend::Scalar => Self::Scalar,
            Backend::Sse2 => Self::Sse2,
            Backend::Sse42Popcnt => Self::Sse42Popcnt,
            Backend::Avx2Popcnt => Self::Avx2Popcnt,
            Backend::Avx512BwPopcnt => Self::Avx512BwPopcnt,
            Backend::Neon => Self::Neon,
        }
    }
}

impl fmt::Display for BackendRequest {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Auto => formatter.write_str("auto"),
            Self::Scalar => Backend::Scalar.fmt(formatter),
            Self::Sse2 => Backend::Sse2.fmt(formatter),
            Self::Sse42Popcnt => Backend::Sse42Popcnt.fmt(formatter),
            Self::Avx2Popcnt => Backend::Avx2Popcnt.fmt(formatter),
            Self::Avx512BwPopcnt => Backend::Avx512BwPopcnt.fmt(formatter),
            Self::Neon => Backend::Neon.fmt(formatter),
        }
    }
}

impl FromStr for BackendRequest {
    type Err = BackendParseError;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        match value {
            "auto" => Ok(Self::Auto),
            "scalar" => Ok(Self::Scalar),
            "sse2" => Ok(Self::Sse2),
            "sse4.2" => Ok(Self::Sse42Popcnt),
            "avx2" => Ok(Self::Avx2Popcnt),
            "avx512" => Ok(Self::Avx512BwPopcnt),
            "neon" => Ok(Self::Neon),
            _ => Err(BackendParseError),
        }
    }
}

/// Invalid textual backend request.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct BackendParseError;

impl fmt::Display for BackendParseError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str("expected auto, scalar, sse2, sse4.2, avx2, avx512, or neon")
    }
}

impl std::error::Error for BackendParseError {}

/// A forced backend is not safe on the supplied CPU.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct BackendUnavailable {
    requested: Backend,
    features: CpuFeatures,
}

impl BackendUnavailable {
    /// Requested concrete backend.
    #[must_use]
    pub const fn requested(self) -> Backend {
        self.requested
    }

    /// Feature set used for validation.
    #[must_use]
    pub const fn features(self) -> CpuFeatures {
        self.features
    }
}

impl fmt::Display for BackendUnavailable {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        if matches!(self.requested, Backend::Scalar) {
            return formatter.write_str("scalar backend is available on every architecture");
        }
        write!(formatter, "SIMD backend {} requires ", self.requested)?;
        let required_architecture = match self.requested {
            Backend::Sse2
            | Backend::Sse42Popcnt
            | Backend::Avx2Popcnt
            | Backend::Avx512BwPopcnt => Architecture::X86_64,
            Backend::Neon => Architecture::Aarch64,
            Backend::Scalar => unreachable!("portable scalar backend has no requirements"),
        };
        if self.features.architecture() != required_architecture {
            return write!(
                formatter,
                "{required_architecture}; current architecture is {}",
                self.features.architecture()
            );
        }
        match self.requested {
            Backend::Scalar => unreachable!("portable scalar backend has no requirements"),
            Backend::Sse2 => formatter.write_str("SSE2; CPU reports SSE2=0"),
            Backend::Sse42Popcnt => {
                formatter.write_str("SSE4.2+POPCNT; missing")?;
                write_missing(
                    formatter,
                    self.features,
                    &[
                        RequiredFeature::Sse41,
                        RequiredFeature::Sse42,
                        RequiredFeature::Popcnt,
                    ],
                )
            }
            Backend::Avx2Popcnt => {
                formatter.write_str("AVX2+POPCNT; missing")?;
                write_missing(
                    formatter,
                    self.features,
                    &[RequiredFeature::Avx2, RequiredFeature::Popcnt],
                )
            }
            Backend::Avx512BwPopcnt => {
                formatter.write_str("AVX-512F+AVX-512BW+AVX2+POPCNT; missing")?;
                write_missing(
                    formatter,
                    self.features,
                    &[
                        RequiredFeature::Avx2,
                        RequiredFeature::Avx512F,
                        RequiredFeature::Avx512Bw,
                        RequiredFeature::Popcnt,
                    ],
                )
            }
            Backend::Neon => formatter.write_str("NEON; CPU reports NEON=0"),
        }
    }
}

#[derive(Clone, Copy)]
enum RequiredFeature {
    Sse41,
    Sse42,
    Avx2,
    Avx512F,
    Avx512Bw,
    Popcnt,
}

fn write_missing(
    formatter: &mut fmt::Formatter<'_>,
    features: CpuFeatures,
    required: &[RequiredFeature],
) -> fmt::Result {
    for feature in required {
        let (available, label) = match feature {
            RequiredFeature::Sse41 => (features.sse41(), "SSE4.1"),
            RequiredFeature::Sse42 => (features.sse42(), "SSE4.2"),
            RequiredFeature::Avx2 => (features.avx2(), "AVX2"),
            RequiredFeature::Avx512F => (features.avx512f(), "AVX-512F"),
            RequiredFeature::Avx512Bw => (features.avx512bw(), "AVX-512BW"),
            RequiredFeature::Popcnt => (features.popcnt(), "POPCNT"),
        };
        if !available {
            write!(formatter, " {label}")?;
        }
    }
    Ok(())
}

impl std::error::Error for BackendUnavailable {}

/// Failure to establish one immutable process CPU configuration.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum InitializationError {
    /// The requested backend lacks at least one CPU feature.
    Unavailable(BackendUnavailable),
    /// A different backend was already published for this process.
    AlreadyInitialized {
        /// Existing immutable backend.
        selected: Backend,
        /// Conflicting requested backend.
        requested: Backend,
    },
}

impl fmt::Display for InitializationError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Unavailable(error) => error.fmt(formatter),
            Self::AlreadyInitialized {
                selected,
                requested,
            } => write!(
                formatter,
                "SIMD backend is already initialized as {selected}; cannot switch to {requested}"
            ),
        }
    }
}

impl std::error::Error for InitializationError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Unavailable(error) => Some(error),
            Self::AlreadyInitialized { .. } => None,
        }
    }
}

/// Initializes CPU detection and backend selection exactly once.
///
/// Repeating `auto` or the installed concrete backend is idempotent. A
/// conflicting request fails instead of mutating dispatch while workers run.
///
/// # Errors
///
/// Returns an unsupported-feature or conflicting-initialization error.
pub fn initialize(request: BackendRequest) -> Result<Configuration, InitializationError> {
    if let Some(configuration) = CONFIGURATION.get() {
        return compatible_existing(*configuration, request);
    }
    let features = detected_features();
    let backend = request
        .select(features)
        .map_err(InitializationError::Unavailable)?;
    let proposed = Configuration { features, backend };
    let installed = *CONFIGURATION.get_or_init(|| proposed);
    compatible_existing(installed, request)
}

fn compatible_existing(
    configuration: Configuration,
    request: BackendRequest,
) -> Result<Configuration, InitializationError> {
    if matches!(request, BackendRequest::Auto)
        || BackendRequest::from(configuration.backend) == request
    {
        Ok(configuration)
    } else {
        Err(InitializationError::AlreadyInitialized {
            selected: configuration.backend,
            requested: match request {
                BackendRequest::Scalar => Backend::Scalar,
                BackendRequest::Sse2 => Backend::Sse2,
                BackendRequest::Sse42Popcnt => Backend::Sse42Popcnt,
                BackendRequest::Avx2Popcnt => Backend::Avx2Popcnt,
                BackendRequest::Avx512BwPopcnt => Backend::Avx512BwPopcnt,
                BackendRequest::Neon => Backend::Neon,
                BackendRequest::Auto => unreachable!("auto is compatible with existing state"),
            },
        })
    }
}

/// Returns the process configuration, selecting automatically on first use.
#[must_use]
pub fn configuration() -> Configuration {
    *CONFIGURATION.get_or_init(|| {
        let features = detected_features();
        let backend = automatic_backend(features);
        Configuration { features, backend }
    })
}

fn detected_features() -> CpuFeatures {
    *DETECTED_CPU_FEATURES.get_or_init(CpuFeatures::detect)
}

/// Returns the concrete process backend without repeating feature detection.
#[must_use]
pub fn selected_backend() -> Backend {
    configuration().backend
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn existing_configuration_accepts_only_compatible_requests() {
        let configuration = Configuration {
            features: CpuFeatures::new_x86_64(true, true, true),
            backend: Backend::Avx2Popcnt,
        };
        assert_eq!(
            compatible_existing(configuration, BackendRequest::Auto),
            Ok(configuration)
        );
        assert_eq!(
            compatible_existing(configuration, BackendRequest::Avx2Popcnt),
            Ok(configuration)
        );
        assert_eq!(
            compatible_existing(configuration, BackendRequest::Scalar),
            Err(InitializationError::AlreadyInitialized {
                selected: Backend::Avx2Popcnt,
                requested: Backend::Scalar,
            })
        );
    }

    #[test]
    fn scalar_configuration_is_idempotent() {
        let configuration = Configuration {
            features: CpuFeatures::new_other(),
            backend: Backend::Scalar,
        };
        assert_eq!(
            compatible_existing(configuration, BackendRequest::Auto),
            Ok(configuration)
        );
        assert_eq!(
            compatible_existing(configuration, BackendRequest::Scalar),
            Ok(configuration)
        );
    }
}
