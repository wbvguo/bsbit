//! Stable command errors, warnings, and successful-run reports.

use core::fmt;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum ErrorClass {
    Usage,
    Operation,
}

/// Stable command failure rendered by the `bsbit` executable.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CliError {
    class: ErrorClass,
    message: String,
}

impl CliError {
    pub(crate) fn usage(message: impl Into<String>) -> Self {
        Self {
            class: ErrorClass::Usage,
            message: message.into(),
        }
    }

    pub(crate) fn operation(message: impl Into<String>) -> Self {
        Self {
            class: ErrorClass::Operation,
            message: message.into(),
        }
    }

    /// Returns 2 for command syntax/unsupported-mode errors and 1 for
    /// operational failures.
    #[must_use]
    pub const fn exit_code(&self) -> u8 {
        match self.class {
            ErrorClass::Usage => 2,
            ErrorClass::Operation => 1,
        }
    }
}

impl fmt::Display for CliError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(&self.message)
    }
}

impl std::error::Error for CliError {}

/// One non-fatal post-publication warning.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CliWarning {
    message: String,
}

impl CliWarning {
    /// Returns the stable warning text.
    #[must_use]
    pub fn message(&self) -> &str {
        &self.message
    }

    pub(crate) fn new(message: impl Into<String>) -> Self {
        Self {
            message: message.into(),
        }
    }
}

/// Complete successful command report.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct RunReport {
    pub(crate) warnings: Vec<CliWarning>,
}

impl RunReport {
    /// Returns post-publication warnings. A warning never changes exit status
    /// because the final target is already visible and complete.
    #[must_use]
    pub fn warnings(&self) -> &[CliWarning] {
        &self.warnings
    }

    pub(crate) fn with_warning(warning: Option<CliWarning>) -> Self {
        Self {
            warnings: warning.into_iter().collect(),
        }
    }
}
