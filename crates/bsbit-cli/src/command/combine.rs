use std::collections::BTreeSet;
use std::ffi::OsString;

use bsbit_combine::{
    Input as CombineInput, MatrixFormat as CombineMatrixFormat, Options as CombineOptions,
    Parameters as CombineParameters,
};

use crate::{COMBINE_HELP, CliError};

use super::{
    Action, option_map, optional_text, optional_u64, parse_bool, parse_compression_threads,
    parse_threads, path, probability_parts_per_billion, required_path,
};

pub(super) fn parse_combine(arguments: &[OsString]) -> Result<Action, CliError> {
    if matches!(arguments, [value] if value == "--help" || value == "-h" || value == "help") {
        return Ok(Action::Help(COMBINE_HELP));
    }
    let normalized = normalize_combine_options(arguments)?;
    let RawCombineOptions {
        input_spec,
        sample_name_spec,
        scalar_options,
        cg_only,
    } = partition_combine_options(&normalized)?;
    let (mut values, _) = option_map(
        &scalar_options,
        &[
            "--output",
            "--matrix",
            "--compress",
            "--threads",
            "--compression-threads",
            "--min-count",
            "--min-prop",
        ],
        &[],
    )?;
    let matrix_format = match optional_text(&mut values, "--matrix")?
        .as_deref()
        .unwrap_or("level")
    {
        "level" => CombineMatrixFormat::Level,
        "count" => CombineMatrixFormat::Count,
        "both" => CombineMatrixFormat::Both,
        value => {
            return Err(CliError::usage(format!(
                "unsupported --matrix `{value}`; expected `level`, `count`, or `both`"
            )));
        }
    };
    let compress = optional_text(&mut values, "--compress")?
        .map(|value| parse_bool("--compress", &value))
        .transpose()?
        .unwrap_or(false);
    let output = required_path(&mut values, "--output")?;
    let threads = parse_threads(&mut values)?;
    let compression_threads = parse_compression_threads(&mut values, compress, threads)?;
    let minimum_count = optional_u64(&mut values, "--min-count")?.unwrap_or(1);
    let minimum_sample_proportion_parts_per_billion =
        probability_parts_per_billion(&mut values, "--min-prop", 0)?;

    let inputs = parse_combine_inputs(input_spec, sample_name_spec.as_ref())?;

    Ok(Action::Combine(CombineOptions {
        inputs,
        output,
        matrix_format,
        compress,
        threads,
        compression_threads,
        parameters: CombineParameters {
            minimum_count,
            minimum_sample_proportion_parts_per_billion,
            cg_only,
        },
    }))
}

struct RawCombineOptions {
    input_spec: OsString,
    sample_name_spec: Option<OsString>,
    scalar_options: Vec<OsString>,
    cg_only: bool,
}

fn partition_combine_options(normalized: &[OsString]) -> Result<RawCombineOptions, CliError> {
    let mut input_spec = None;
    let mut sample_name_spec = None;
    let mut scalar_options = Vec::new();
    let mut cg_only = false;
    let mut cursor = 0;
    while cursor < normalized.len() {
        let option = normalized[cursor]
            .to_str()
            .ok_or_else(|| CliError::usage("option name must be valid UTF-8"))?;
        if !option.starts_with("--") {
            return Err(CliError::usage(format!("unknown option `{option}`")));
        }
        if option == "--cg-only" {
            if cg_only {
                return Err(CliError::usage("duplicate option `--cg-only`"));
            }
            cg_only = true;
            cursor += 1;
            continue;
        }
        let Some(value) = normalized.get(cursor + 1) else {
            return Err(CliError::usage(format!("missing value for `{option}`")));
        };
        if value.to_str().is_some_and(|value| value.starts_with("--")) {
            return Err(CliError::usage(format!("missing value for `{option}`")));
        }
        match option {
            "--input" => set_unique_list(&mut input_spec, value, option)?,
            "--sample-name" => set_unique_list(&mut sample_name_spec, value, option)?,
            _ => scalar_options.extend([normalized[cursor].clone(), value.clone()]),
        }
        cursor += 2;
    }
    let input_spec = input_spec
        .ok_or_else(|| CliError::usage("missing required option `--input`; use PATH[,PATH...]"))?;
    Ok(RawCombineOptions {
        input_spec,
        sample_name_spec,
        scalar_options,
        cg_only,
    })
}

fn set_unique_list(
    target: &mut Option<OsString>,
    value: &OsString,
    option: &str,
) -> Result<(), CliError> {
    if target.replace(value.clone()).is_some() {
        return Err(CliError::usage(format!(
            "duplicate option `{option}`; provide one comma-separated list"
        )));
    }
    Ok(())
}

fn parse_combine_inputs(
    input_spec: OsString,
    sample_name_spec: Option<&OsString>,
) -> Result<Vec<CombineInput>, CliError> {
    let input_paths = comma_separated_combine_paths(input_spec)?;
    let sample_names_were_supplied = sample_name_spec.is_some();
    let supplied_sample_names = sample_name_spec
        .map(|specification| comma_separated_combine_text("--sample-name", specification))
        .transpose()?
        .unwrap_or_default();
    if sample_names_were_supplied && supplied_sample_names.len() != input_paths.len() {
        return Err(CliError::usage(format!(
            "--sample-name supplies {} name(s), but --input supplies {} path(s)",
            supplied_sample_names.len(),
            input_paths.len()
        )));
    }

    let named_paths = if sample_names_were_supplied {
        supplied_sample_names
            .into_iter()
            .zip(input_paths)
            .collect::<Vec<_>>()
    } else {
        input_paths
            .into_iter()
            .map(|path| {
                let sample = path.to_str().ok_or_else(|| {
                    CliError::usage("a non-UTF-8 --input path requires an explicit --sample-name")
                })?;
                Ok((sample.to_owned(), path))
            })
            .collect::<Result<Vec<_>, CliError>>()?
    };

    let mut sample_names = BTreeSet::new();
    named_paths
        .into_iter()
        .map(|(sample, input_path)| {
            if sample.is_empty() || sample.bytes().any(|byte| byte.is_ascii_control()) {
                return Err(CliError::usage(format!(
                    "invalid sample label `{sample}`; labels must be nonempty and contain no control bytes"
                )));
            }
            if !sample_names.insert(sample.clone()) {
                return Err(CliError::usage(format!(
                    "duplicate sample label `{sample}`"
                )));
            }
            Ok(CombineInput {
                sample,
                path: path("--input", input_path)?,
            })
        })
        .collect()
}

fn comma_separated_combine_text(
    option: &str,
    specification: &OsString,
) -> Result<Vec<String>, CliError> {
    let mut values = Vec::new();
    let specification = specification
        .to_str()
        .ok_or_else(|| CliError::usage(format!("value for `{option}` must be valid UTF-8")))?;
    for value in specification.split(',') {
        if value.is_empty() {
            return Err(CliError::usage(format!(
                "empty item in comma-separated `{option}` value `{specification}`"
            )));
        }
        values.push(value.to_owned());
    }
    Ok(values)
}

fn comma_separated_combine_paths(specification: OsString) -> Result<Vec<OsString>, CliError> {
    let mut values = Vec::new();
    if let Some(text) = specification.to_str() {
        for value in text.split(',') {
            if value.is_empty() {
                return Err(CliError::usage(format!(
                    "empty item in comma-separated `--input` value `{text}`"
                )));
            }
            values.push(OsString::from(value));
        }
    } else {
        values.push(specification);
    }
    Ok(values)
}

fn normalize_combine_options(arguments: &[OsString]) -> Result<Vec<OsString>, CliError> {
    arguments
        .iter()
        .map(|argument| {
            Ok(match argument.to_str() {
                Some("-i") => OsString::from("--input"),
                Some("-o") => OsString::from("--output"),
                Some("-c") => OsString::from("--compress"),
                Some("-t") => OsString::from("--threads"),
                Some(value) if value.starts_with('-') && !value.starts_with("--") => {
                    return Err(CliError::usage(format!("unknown option `{value}`")));
                }
                _ => argument.clone(),
            })
        })
        .collect()
}
