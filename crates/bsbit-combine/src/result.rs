//! Public failures, warnings, and completion report for matrix assembly.

use core::fmt;

/// Stable high-level class for one combine failure.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum CombineErrorKind {
    /// Options violate the public contract.
    Configuration,
    /// An input path, stream, row, or ordering contract is invalid.
    Input,
    /// A worker stopped without returning a normal result.
    Worker,
    /// Output staging or encoding failed.
    Output,
    /// Replacement publication or rollback failed.
    Publication,
}

/// One methylation-matrix assembly failure.
#[derive(Debug)]
pub struct CombineError {
    kind: CombineErrorKind,
    context: String,
    source: Option<Box<dyn std::error::Error + Send + Sync>>,
}

impl CombineError {
    pub(crate) fn configuration(message: impl Into<String>) -> Self {
        Self {
            kind: CombineErrorKind::Configuration,
            context: message.into(),
            source: None,
        }
    }

    pub(crate) fn input(message: impl Into<String>) -> Self {
        Self {
            kind: CombineErrorKind::Input,
            context: message.into(),
            source: None,
        }
    }

    pub(crate) fn worker(message: impl Into<String>) -> Self {
        Self {
            kind: CombineErrorKind::Worker,
            context: message.into(),
            source: None,
        }
    }

    pub(crate) fn with_source(
        kind: CombineErrorKind,
        context: impl Into<String>,
        source: impl std::error::Error + Send + Sync + 'static,
    ) -> Self {
        Self {
            kind,
            context: context.into(),
            source: Some(Box::new(source)),
        }
    }

    /// Returns the stable high-level failure class.
    #[must_use]
    pub const fn kind(&self) -> CombineErrorKind {
        self.kind
    }
}

impl fmt::Display for CombineError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(&self.context)?;
        if let Some(source) = &self.source {
            write!(formatter, ": {source}")?;
        }
        Ok(())
    }
}

impl std::error::Error for CombineError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        self.source
            .as_deref()
            .map(|source| source as &(dyn std::error::Error + 'static))
    }
}

/// One non-fatal warning returned after successful publication.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CombineWarning {
    pub(crate) message: String,
}

impl CombineWarning {
    /// Returns the warning text.
    #[must_use]
    pub fn message(&self) -> &str {
        &self.message
    }
}

/// Summary of one successfully published combine operation.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct CombineReport {
    pub(crate) sites_seen: u64,
    pub(crate) sites_written: u64,
    pub(crate) warnings: Vec<CombineWarning>,
}

impl CombineReport {
    /// Returns the number of distinct input sites considered.
    #[must_use]
    pub const fn sites_seen(&self) -> u64 {
        self.sites_seen
    }

    /// Returns the number of sites retained after filtering.
    #[must_use]
    pub const fn sites_written(&self) -> u64 {
        self.sites_written
    }

    /// Returns post-publication warnings.
    #[must_use]
    pub fn warnings(&self) -> &[CombineWarning] {
        &self.warnings
    }
}
