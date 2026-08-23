//! Streaming methylation-matrix assembly from sorted extended bedMethyl files.

#![forbid(unsafe_code)]

mod input;
mod merge;
mod output;
mod request;
mod result;
mod site;

pub use merge::combine;
pub use request::{Input, MatrixFormat, Options, Parameters};
pub use result::{CombineError, CombineErrorKind, CombineReport, CombineWarning};
