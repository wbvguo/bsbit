pub(crate) mod align;
mod align_options;
mod alignment_input;
mod alignment_metrics;
mod call;
mod combine;
pub(crate) mod cpu;
pub(crate) mod index;
mod paired_end;
pub(crate) mod single_end;

use std::collections::{BTreeMap, BTreeSet};
use std::ffi::{OsStr, OsString};
use std::path::{Path, PathBuf};

use crate::{CliError, GENERAL_HELP};
use bsbit_call::joint::Options as CallJointOptions;
use bsbit_call::meth::Options as CallMethOptions;
use bsbit_call::snp::Options as CallSnpOptions;
use bsbit_combine::Options as CombineOptions;

use call::parse_call;
use combine::parse_combine;
use index::{IndexOptions, parse_index};

/// Derives the hidden search-image prefix associated with one opaque index
/// handle. Index construction and both alignment layouts must use this one
/// physical-layout rule without depending on one another's command modules.
fn internal_search_file_prefix(index: &Path) -> PathBuf {
    const FNV_OFFSET_BASIS: u64 = 0xcbf2_9ce4_8422_2325;
    const FNV_PRIME: u64 = 0x0000_0100_0000_01b3;
    let mut hash = FNV_OFFSET_BASIS;
    for byte in index
        .file_name()
        .unwrap_or(index.as_os_str())
        .as_encoded_bytes()
    {
        hash ^= u64::from(*byte);
        hash = hash.wrapping_mul(FNV_PRIME);
    }
    index.with_file_name(format!(".bsbit-index-{hash:016x}"))
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) enum Action {
    Help(&'static str),
    Version,
    Cpu(cpu::Options),
    Index(IndexOptions),
    Align(align_options::Options),
    CallMeth(CallMethOptions),
    CallSnp(CallSnpOptions),
    CallJoint(CallJointOptions),
    Combine(CombineOptions),
}

pub(crate) fn parse(arguments: impl IntoIterator<Item = OsString>) -> Result<Action, CliError> {
    let arguments = arguments.into_iter().collect::<Vec<_>>();
    let Some(command) = arguments.first() else {
        return Err(CliError::usage("missing command; run `bsbit --help`"));
    };
    let command = command
        .to_str()
        .ok_or_else(|| CliError::usage("command name must be valid UTF-8"))?;
    match command {
        "--help" | "-h" | "help" if arguments.len() == 1 => Ok(Action::Help(GENERAL_HELP)),
        "--version" | "-v" if arguments.len() == 1 => Ok(Action::Version),
        "index" => parse_index(&arguments[1..]),
        "align" => align::parse(&arguments[1..]),
        "cpu" => cpu::parse(&arguments[1..]),
        "call" => parse_call(&arguments[1..]),
        "combine" => parse_combine(&arguments[1..]),
        value => Err(CliError::usage(format!(
            "unknown command `{value}`; run `bsbit --help`"
        ))),
    }
}

type OptionValues = BTreeMap<String, OsString>;

fn parse_threads(values: &mut OptionValues) -> Result<u64, CliError> {
    let threads = optional_u64(values, "--threads")?.unwrap_or(1);
    if threads == 0 {
        return Err(CliError::usage("--threads must be positive"));
    }
    if u32::try_from(threads).is_err() {
        return Err(CliError::usage(
            "--threads exceeds the supported u32 worker domain",
        ));
    }
    Ok(threads)
}

fn parse_compression_threads(
    values: &mut OptionValues,
    compress: bool,
    worker_threads: u64,
) -> Result<u32, CliError> {
    let default = u64::from(compress && worker_threads > 1);
    let value = optional_u64(values, "--compression-threads")?.unwrap_or(default);
    let threads = u32::try_from(value)
        .ok()
        .filter(|threads| i32::try_from(*threads).is_ok())
        .ok_or_else(|| {
            CliError::usage(
                "--compression-threads must fit the native nonnegative signed 32-bit worker domain",
            )
        })?;
    if !compress && threads != 0 {
        return Err(CliError::usage(
            "--compression-threads must be 0 when --compress is false",
        ));
    }
    Ok(threads)
}

fn probability_parts_per_billion(
    values: &mut OptionValues,
    option: &str,
    default: u32,
) -> Result<u32, CliError> {
    let Some(value) = values.remove(option) else {
        return Ok(default);
    };
    let value = text_value(option, &value)?;
    let mut parts = value.split('.');
    let whole = parts.next().unwrap_or_default();
    let fraction = parts.next();
    let invalid_fraction = fraction.is_some_and(|digits| {
        digits.is_empty() || digits.len() > 9 || !digits.bytes().all(|byte| byte.is_ascii_digit())
    });
    let nonzero_fraction_at_one =
        whole == "1" && fraction.is_some_and(|digits| !digits.bytes().all(|byte| byte == b'0'));
    if parts.next().is_some()
        || !matches!(whole, "0" | "1")
        || invalid_fraction
        || nonzero_fraction_at_one
    {
        return Err(CliError::usage(format!(
            "invalid value `{value}` for `{option}`; expected a decimal probability in 0..=1 with at most 9 fractional digits"
        )));
    }
    if whole == "1" {
        return Ok(1_000_000_000);
    }
    let Some(fraction) = fraction else {
        return Ok(0);
    };
    let fractional = fraction.parse::<u32>().map_err(|_| {
        CliError::usage(format!(
            "invalid decimal probability `{value}` for `{option}`"
        ))
    })?;
    let exponent = u32::try_from(9 - fraction.len()).expect("fraction length is at most 9");
    let scale = 10_u32
        .checked_pow(exponent)
        .expect("probability scale fits u32");
    fractional
        .checked_mul(scale)
        .ok_or_else(|| CliError::usage(format!("probability `{value}` overflows for `{option}`")))
}

fn parse_bool(option: &str, value: &str) -> Result<bool, CliError> {
    match value {
        "true" => Ok(true),
        "false" => Ok(false),
        _ => Err(CliError::usage(format!(
            "invalid value `{value}` for `{option}`; expected `true` or `false`"
        ))),
    }
}

fn option_map(
    arguments: &[OsString],
    accepted: &[&str],
    accepted_flags: &[&str],
) -> Result<(OptionValues, BTreeSet<String>), CliError> {
    option_map_with_short_options(arguments, accepted, accepted_flags, &[])
}

fn option_map_with_short_options(
    arguments: &[OsString],
    accepted: &[&str],
    accepted_flags: &[&str],
    short_options: &[(&str, &str)],
) -> Result<(OptionValues, BTreeSet<String>), CliError> {
    let mut values = BTreeMap::new();
    let mut flags = BTreeSet::new();
    let mut cursor = 0;
    while cursor < arguments.len() {
        let supplied = arguments[cursor]
            .to_str()
            .ok_or_else(|| CliError::usage("option name must be valid UTF-8"))?;
        let option = short_options
            .iter()
            .find_map(|&(alias, canonical)| (supplied == alias).then_some(canonical))
            .unwrap_or(supplied);
        if accepted_flags.contains(&option) {
            if !flags.insert(option.to_owned()) {
                return Err(CliError::usage(format!("duplicate option `{option}`")));
            }
            cursor += 1;
            continue;
        }
        if !option.starts_with("--") || !accepted.contains(&option) {
            return Err(CliError::usage(format!("unknown option `{option}`")));
        }
        let Some(value) = arguments.get(cursor + 1) else {
            return Err(CliError::usage(format!("missing value for `{option}`")));
        };
        let value_is_option = value.to_str().is_some_and(|text| {
            text.starts_with("--") || short_options.iter().any(|(short, _)| text == *short)
        });
        if value_is_option {
            return Err(CliError::usage(format!("missing value for `{option}`")));
        }
        if values.insert(option.to_owned(), value.clone()).is_some() {
            return Err(CliError::usage(format!("duplicate option `{option}`")));
        }
        cursor += 2;
    }
    Ok((values, flags))
}

fn required(values: &mut OptionValues, option: &str) -> Result<OsString, CliError> {
    values
        .remove(option)
        .ok_or_else(|| CliError::usage(format!("missing required option `{option}`")))
}

fn required_text(values: &mut OptionValues, option: &str) -> Result<String, CliError> {
    let value = required(values, option)?;
    Ok(text_value(option, &value)?.to_owned())
}

fn optional_text(values: &mut OptionValues, option: &str) -> Result<Option<String>, CliError> {
    values
        .remove(option)
        .map(|value| text_value(option, &value).map(str::to_owned))
        .transpose()
}

fn optional_u64(values: &mut OptionValues, option: &str) -> Result<Option<u64>, CliError> {
    values
        .remove(option)
        .map(|value| text_value(option, &value).and_then(|value| parse_u64(option, value)))
        .transpose()
}

fn parse_u64(option: &str, value: &str) -> Result<u64, CliError> {
    value.parse::<u64>().map_err(|_| {
        CliError::usage(format!(
            "invalid value `{value}` for `{option}`; expected a nonnegative integer"
        ))
    })
}

fn required_path(values: &mut OptionValues, option: &str) -> Result<PathBuf, CliError> {
    path(option, required(values, option)?)
}

fn optional_path(values: &mut OptionValues, option: &str) -> Result<Option<PathBuf>, CliError> {
    values
        .remove(option)
        .map(|value| path(option, value))
        .transpose()
}

fn path(option: &str, value: OsString) -> Result<PathBuf, CliError> {
    if value.is_empty() {
        return Err(CliError::usage(format!("empty path for `{option}`")));
    }
    if value == OsStr::new("-") || value.to_str().is_some_and(|text| text.contains("://")) {
        return Err(CliError::usage(format!(
            "unsupported non-local path `{}` for `{option}`",
            value.to_string_lossy()
        )));
    }
    let path = PathBuf::from(value);
    if path.file_name().is_none() {
        return Err(CliError::usage(format!(
            "path for `{option}` must name a file"
        )));
    }
    Ok(path)
}

fn text_value<'a>(option: &str, value: &'a OsStr) -> Result<&'a str, CliError> {
    value.to_str().ok_or_else(|| {
        CliError::usage(format!(
            "value for `{option}` must be valid UTF-8; filesystem paths may use arbitrary platform bytes"
        ))
    })
}

#[cfg(test)]
mod tests;
