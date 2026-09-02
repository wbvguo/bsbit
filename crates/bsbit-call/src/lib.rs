//! Bit-parallel methylation, SNP, and joint calling from canonical bsbit BAM.
//!
//! The library owns evidence reconstruction, overlapping-mate collapse,
//! regional parallelism, likelihood evaluation, and biological output
//! rendering. BAM field access, BGZF transport, and atomic publication are
//! delegated to `bsbit-hts`; command-line parsing remains in `bsbit-cli`.

#![forbid(unsafe_code)]

mod call_input;
mod evidence;
mod publication;
mod reference_context;

pub mod joint;
pub mod meth;
pub mod region;
pub mod snp;

use core::fmt;

/// Maximum supported regional calling worker count.
pub const MAX_THREADS: u64 = 64;

/// Stable high-level class for one calling failure.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum CallErrorKind {
    /// Command options violate the public contract.
    Configuration,
    /// BAM opening, indexing, header, or record evidence failed.
    Input,
    /// Regional aggregation or likelihood evaluation failed.
    Calling,
    /// Staging creation, encoding, or finalization failed.
    Output,
    /// Replacement publication or transactional rollback failed.
    Publication,
}

/// One operational calling failure.
#[derive(Debug)]
pub struct CallError {
    kind: CallErrorKind,
    context: String,
    source: Option<Box<dyn std::error::Error + Send + Sync>>,
}

impl CallError {
    pub(crate) fn operation(message: impl Into<String>) -> Self {
        Self {
            kind: CallErrorKind::Calling,
            context: message.into(),
            source: None,
        }
    }

    pub(crate) fn configuration(message: impl Into<String>) -> Self {
        Self {
            kind: CallErrorKind::Configuration,
            context: message.into(),
            source: None,
        }
    }

    pub(crate) fn input(message: impl Into<String>) -> Self {
        Self {
            kind: CallErrorKind::Input,
            context: message.into(),
            source: None,
        }
    }

    pub(crate) fn with_source(
        kind: CallErrorKind,
        context: impl Into<String>,
        source: impl std::error::Error + Send + Sync + 'static,
    ) -> Self {
        Self {
            kind,
            context: context.into(),
            source: Some(Box::new(source)),
        }
    }

    pub(crate) fn with_context(self, context: impl Into<String>) -> Self {
        let kind = self.kind;
        Self {
            kind,
            context: context.into(),
            source: Some(Box::new(self)),
        }
    }

    /// Returns the stable high-level failure class.
    #[must_use]
    pub const fn kind(&self) -> CallErrorKind {
        self.kind
    }
}

impl fmt::Display for CallError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(&self.context)?;
        if let Some(source) = &self.source {
            write!(formatter, ": {source}")?;
        }
        Ok(())
    }
}

impl std::error::Error for CallError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        self.source
            .as_deref()
            .map(|source| source as &(dyn std::error::Error + 'static))
    }
}

/// One non-fatal warning produced after a successful calling run.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CallWarning {
    message: String,
}

impl CallWarning {
    pub(crate) fn new(message: impl Into<String>) -> Self {
        Self {
            message: message.into(),
        }
    }

    /// Returns the stable warning text.
    #[must_use]
    pub fn message(&self) -> &str {
        &self.message
    }
}

/// Complete successful calling report.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct CallReport {
    warnings: Vec<CallWarning>,
}

impl CallReport {
    /// Returns warnings emitted after or alongside successful publication.
    #[must_use]
    pub fn warnings(&self) -> &[CallWarning] {
        &self.warnings
    }

    pub(crate) fn with_warning(warning: Option<CallWarning>) -> Self {
        Self {
            warnings: warning.into_iter().collect(),
        }
    }

    pub(crate) fn with_prior_warning(mut self, warning: Option<CallWarning>) -> Self {
        if let Some(warning) = warning {
            self.warnings.insert(0, warning);
        }
        self
    }
}

pub(crate) fn validate_threads(command: &str, threads: u64) -> Result<(), CallError> {
    if (1..=MAX_THREADS).contains(&threads) {
        Ok(())
    } else {
        Err(CallError::configuration(format!(
            "{command}: thread count must be within 1..={MAX_THREADS}"
        )))
    }
}
