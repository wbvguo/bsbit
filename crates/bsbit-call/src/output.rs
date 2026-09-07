//! Direct text-output creation and finalization shared by all call modes.

use std::io;
use std::path::{Path, PathBuf};

use bsbit_hts::{TextOutputCompression, TextStagingWriter};

use crate::{CallError, CallErrorKind};

pub(crate) fn create_text_output(
    command: &str,
    target: &Path,
    input: &Path,
    reference: &Path,
    additional_protected_paths: &[&Path],
    compress: bool,
    compression_threads: u32,
) -> Result<TextStagingWriter, CallError> {
    let compression = if compress {
        TextOutputCompression::Bgzf
    } else {
        TextOutputCompression::Plain
    };
    let sidecars = [
        append_path_suffix(input, ".bai"),
        append_path_suffix(input, ".csi"),
        input.with_extension("bai"),
        input.with_extension("csi"),
        append_path_suffix(reference, ".fai"),
        append_path_suffix(reference, ".gzi"),
    ];
    let mut protected_paths = vec![input, reference];
    protected_paths.extend(sidecars.iter().map(PathBuf::as_path));
    protected_paths.extend_from_slice(additional_protected_paths);
    TextStagingWriter::create_direct_distinct_from(
        target,
        &protected_paths,
        compression,
        compression_threads,
    )
    .map_err(|error| {
        CallError::with_source(
            CallErrorKind::Output,
            format!("{command}: open output {}", target.display()),
            error,
        )
    })
}

fn append_path_suffix(path: &Path, suffix: &str) -> PathBuf {
    let mut value = path.as_os_str().to_owned();
    value.push(suffix);
    PathBuf::from(value)
}

pub(crate) fn finish_output(command: &str, output: TextStagingWriter) -> Result<(), CallError> {
    output.finish_direct().map_err(|error| {
        CallError::with_source(
            CallErrorKind::Output,
            format!("{command}: finalize output"),
            error,
        )
    })
}

pub(crate) fn output_write_error(command: &str, target: &Path, error: io::Error) -> CallError {
    CallError::with_source(
        CallErrorKind::Output,
        format!("{command}: write output {}", target.display()),
        error,
    )
}
