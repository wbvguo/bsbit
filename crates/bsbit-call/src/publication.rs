//! Replaceable text staging and publication shared by all call modes.

use std::io;
use std::path::Path;

use bsbit_hts::{TextOutputCompression, TextPublication, TextStagingWriter};

use crate::{CallError, CallErrorKind, CallWarning};

pub(crate) fn create_text_staging(
    command: &str,
    target: &Path,
    compress: bool,
    threads: u64,
) -> Result<TextStagingWriter, CallError> {
    let compression = if compress {
        TextOutputCompression::Bgzf
    } else {
        TextOutputCompression::Plain
    };
    // Region calculation and encoding overlap. One private BGZF worker keeps
    // compression off the ordered writer without multiplying `-t` threads.
    let compression_threads = u32::from(compress && threads > 1);
    TextStagingWriter::create_sibling_replace(target, compression, compression_threads).map_err(
        |error| {
            CallError::with_source(
                CallErrorKind::Output,
                format!("{command}: create output staging for {}", target.display()),
                error,
            )
        },
    )
}

pub(crate) fn finish_and_publish(
    command: &str,
    output: TextStagingWriter,
) -> Result<TextPublication, CallError> {
    let completed = output.finish().map_err(|error| {
        CallError::with_source(
            CallErrorKind::Output,
            format!("{command}: finalize output"),
            error,
        )
    })?;
    completed.publish_replace().map_err(|error| {
        CallError::with_source(
            CallErrorKind::Publication,
            format!("{command}: publish output"),
            error,
        )
    })
}

pub(crate) fn publication_warning(
    publication: &TextPublication,
    artifact_label: &str,
) -> Option<CallWarning> {
    publication.cleanup_warning().map(|warning| {
        CallWarning::new(format!(
            "{artifact_label} {} was published, but staging {} could not be cleaned: {warning:?}",
            publication.target_path().display(),
            publication.staging_path().display()
        ))
    })
}

pub(crate) fn output_write_error(command: &str, target: &Path, error: io::Error) -> CallError {
    CallError::with_source(
        CallErrorKind::Output,
        format!("{command}: write output staging for {}", target.display()),
        error,
    )
}
