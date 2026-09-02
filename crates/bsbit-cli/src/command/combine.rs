use std::collections::BTreeSet;

use bsbit_combine::{
    Input as CombineInput, MatrixFormat as CombineMatrixFormat, Options as CombineOptions,
    Parameters as CombineParameters,
};

use crate::{COMBINE_HELP, CliError};

use super::{
    Action, option_map, optional_u64, parse_bool, parse_threads, path,
    probability_parts_per_billion, required_path,
};

pub(super) fn parse_combine(arguments: &[String]) -> Result<Action, CliError> {
    if matches!(arguments, [value] if value == "--help" || value == "-h" || value == "help") {
        return Ok(Action::Help(COMBINE_HELP));
    }
    let normalized = normalize_combine_options(arguments)?;
    let mut input_specs = Vec::new();
    let mut sample_name_specs = Vec::new();
    let mut scalar_options = Vec::new();
    let mut cg_only = false;
    let mut cursor = 0;
    while cursor < normalized.len() {
        let option = &normalized[cursor];
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
        if value.starts_with("--") {
            return Err(CliError::usage(format!("missing value for `{option}`")));
        }
        match option.as_str() {
            "--input" => input_specs.push(value.clone()),
            "--sample-name" => sample_name_specs.push(value.clone()),
            _ => {
                scalar_options.push(option.clone());
                scalar_options.push(value.clone());
            }
        }
        cursor += 2;
    }
    if input_specs.is_empty() {
        return Err(CliError::usage(
            "missing required option `--input`; use PATH[,PATH...]",
        ));
    }
    let (mut values, _) = option_map(
        &scalar_options,
        &[
            "--prefix",
            "--matrix",
            "--compress",
            "--threads",
            "--min-count",
            "--min-prop",
        ],
        &[],
    )?;
    let output_prefix = required_path(&mut values, "--prefix")?;
    let matrix_format = match values.remove("--matrix").as_deref().unwrap_or("level") {
        "level" => CombineMatrixFormat::Level,
        "count" => CombineMatrixFormat::Count,
        "both" => CombineMatrixFormat::Both,
        value => {
            return Err(CliError::usage(format!(
                "unsupported --matrix `{value}`; expected `level`, `count`, or `both`"
            )));
        }
    };
    let compress = values
        .remove("--compress")
        .map(|value| parse_bool("--compress", &value))
        .transpose()?
        .unwrap_or(true);
    let threads = parse_threads(&mut values)?;
    let minimum_count = optional_u64(&mut values, "--min-count")?.unwrap_or(1);
    let minimum_sample_proportion_parts_per_billion =
        probability_parts_per_billion(&mut values, "--min-prop", 0)?;

    let inputs = parse_combine_inputs(input_specs, sample_name_specs)?;

    Ok(Action::Combine(CombineOptions {
        inputs,
        output_prefix,
        matrix_format,
        compress,
        threads,
        parameters: CombineParameters {
            minimum_count,
            minimum_sample_proportion_parts_per_billion,
            cg_only,
        },
    }))
}

fn parse_combine_inputs(
    input_specs: Vec<String>,
    sample_name_specs: Vec<String>,
) -> Result<Vec<CombineInput>, CliError> {
    let input_paths = comma_separated_combine_values("--input", input_specs)?;
    let sample_names_were_supplied = !sample_name_specs.is_empty();
    let supplied_sample_names = comma_separated_combine_values("--sample-name", sample_name_specs)?;
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
            .map(|path| (path.clone(), path))
            .collect::<Vec<_>>()
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

fn comma_separated_combine_values(
    option: &str,
    specifications: Vec<String>,
) -> Result<Vec<String>, CliError> {
    let mut values = Vec::new();
    for specification in specifications {
        for value in specification.split(',') {
            if value.is_empty() {
                return Err(CliError::usage(format!(
                    "empty item in comma-separated `{option}` value `{specification}`"
                )));
            }
            values.push(value.to_owned());
        }
    }
    Ok(values)
}

fn normalize_combine_options(arguments: &[String]) -> Result<Vec<String>, CliError> {
    arguments
        .iter()
        .map(|argument| {
            Ok(match argument.as_str() {
                "-i" => String::from("--input"),
                "-p" => String::from("--prefix"),
                "-m" => String::from("--matrix"),
                "-c" => String::from("--compress"),
                "-t" => String::from("--threads"),
                value if value.starts_with('-') && !value.starts_with("--") => {
                    return Err(CliError::usage(format!("unknown option `{value}`")));
                }
                _ => argument.clone(),
            })
        })
        .collect()
}
