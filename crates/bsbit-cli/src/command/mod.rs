pub(crate) mod align;
mod call;
mod combine;
pub(crate) mod index;
pub(crate) mod single_end;

use std::collections::{BTreeMap, BTreeSet};
use std::ffi::OsString;
use std::path::{Path, PathBuf};

use bsbit_align::library::LibraryProfile;
use bsbit_call::joint::Options as CallJointOptions;
use bsbit_call::meth::Options as CallMethOptions;
use bsbit_call::snp::Options as CallSnpOptions;
use bsbit_combine::Options as CombineOptions;
use bsbit_hts::BsbitAlignmentMode;
use bsbit_io::select_sibling_staging_path;

use crate::{CliError, GENERAL_HELP};

pub(crate) use index::IndexOptions;

use call::parse_call;
use combine::parse_combine;
use index::parse_index;

pub(crate) const MAX_CLI_THREADS: u64 = 64;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum ReadLayout {
    SingleEnd,
    PairedEnd,
}

const fn caller_compatible_alignment_mode(
    layout: ReadLayout,
    profile: LibraryProfile,
) -> BsbitAlignmentMode {
    match (layout, profile) {
        (ReadLayout::SingleEnd, LibraryProfile::Directional) => {
            BsbitAlignmentMode::CallerCompatibleDirectionalSingle
        }
        (ReadLayout::SingleEnd, LibraryProfile::NonDirectional) => {
            BsbitAlignmentMode::CallerCompatibleNondirectionalSingle
        }
        (ReadLayout::PairedEnd, LibraryProfile::Directional) => {
            BsbitAlignmentMode::CallerCompatibleDirectionalPaired
        }
        (ReadLayout::PairedEnd, LibraryProfile::NonDirectional) => {
            BsbitAlignmentMode::CallerCompatibleNondirectionalPaired
        }
    }
}

fn unused_staging_path(
    target: &Path,
    label: &str,
    error_subject: &str,
) -> Result<PathBuf, CliError> {
    select_sibling_staging_path(target, label).map_err(|error| {
        CliError::operation(format!(
            "{error_subject}: inspect staging path {}: {error}",
            target.display()
        ))
    })
}

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
    Index(IndexOptions),
    Align(align::Options),
    CallMeth(CallMethOptions),
    CallSnp(CallSnpOptions),
    CallJoint(CallJointOptions),
    Combine(CombineOptions),
}

pub(crate) fn parse(arguments: impl IntoIterator<Item = OsString>) -> Result<Action, CliError> {
    let arguments = arguments
        .into_iter()
        .map(|value| {
            value
                .into_string()
                .map_err(|_| CliError::usage("arguments and paths must be valid UTF-8"))
        })
        .collect::<Result<Vec<_>, _>>()?;
    let Some(command) = arguments.first().map(String::as_str) else {
        return Err(CliError::usage("missing command; run `bsbit --help`"));
    };
    match command {
        "--help" | "-h" | "help" if arguments.len() == 1 => Ok(Action::Help(GENERAL_HELP)),
        "--version" | "-V" if arguments.len() == 1 => Ok(Action::Version),
        "index" => parse_index(&arguments[1..]),
        "align" => align::parse(&arguments[1..]),
        "call" => parse_call(&arguments[1..]),
        "combine" => parse_combine(&arguments[1..]),
        value => Err(CliError::usage(format!(
            "unknown command `{value}`; run `bsbit --help`"
        ))),
    }
}

fn parse_threads(values: &mut BTreeMap<String, String>) -> Result<u64, CliError> {
    let threads = optional_u64(values, "--threads")?.unwrap_or(1);
    if !(1..=MAX_CLI_THREADS).contains(&threads) {
        return Err(CliError::usage(format!(
            "--threads must be in 1..={MAX_CLI_THREADS}"
        )));
    }
    Ok(threads)
}

fn probability_parts_per_billion(
    values: &mut BTreeMap<String, String>,
    option: &str,
    default: u32,
) -> Result<u32, CliError> {
    let Some(value) = values.remove(option) else {
        return Ok(default);
    };
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
    arguments: &[String],
    accepted: &[&str],
    accepted_flags: &[&str],
) -> Result<(BTreeMap<String, String>, BTreeSet<String>), CliError> {
    option_map_with_aliases(arguments, accepted, accepted_flags, &[])
}

fn option_map_with_aliases(
    arguments: &[String],
    accepted: &[&str],
    accepted_flags: &[&str],
    aliases: &[(&str, &str)],
) -> Result<(BTreeMap<String, String>, BTreeSet<String>), CliError> {
    let mut values = BTreeMap::new();
    let mut flags = BTreeSet::new();
    let mut cursor = 0;
    while cursor < arguments.len() {
        let supplied = &arguments[cursor];
        let option = aliases
            .iter()
            .find_map(|&(alias, canonical)| (supplied == alias).then_some(canonical))
            .unwrap_or(supplied.as_str());
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
        if value.starts_with("--") || aliases.iter().any(|(alias, _)| value == alias) {
            return Err(CliError::usage(format!("missing value for `{option}`")));
        }
        if values.insert(option.to_owned(), value.clone()).is_some() {
            return Err(CliError::usage(format!("duplicate option `{option}`")));
        }
        cursor += 2;
    }
    Ok((values, flags))
}

fn required(values: &mut BTreeMap<String, String>, option: &str) -> Result<String, CliError> {
    values
        .remove(option)
        .ok_or_else(|| CliError::usage(format!("missing required option `{option}`")))
}

fn optional_u64(
    values: &mut BTreeMap<String, String>,
    option: &str,
) -> Result<Option<u64>, CliError> {
    values
        .remove(option)
        .map(|value| parse_u64(option, &value))
        .transpose()
}

fn parse_u64(option: &str, value: &str) -> Result<u64, CliError> {
    value.parse::<u64>().map_err(|_| {
        CliError::usage(format!(
            "invalid value `{value}` for `{option}`; expected a nonnegative integer"
        ))
    })
}

fn required_path(values: &mut BTreeMap<String, String>, option: &str) -> Result<PathBuf, CliError> {
    path(option, required(values, option)?)
}

fn optional_path(
    values: &mut BTreeMap<String, String>,
    option: &str,
) -> Result<Option<PathBuf>, CliError> {
    values
        .remove(option)
        .map(|value| path(option, value))
        .transpose()
}

fn path(option: &str, value: String) -> Result<PathBuf, CliError> {
    if value.is_empty() {
        return Err(CliError::usage(format!("empty path for `{option}`")));
    }
    if value == "-" || value.contains("://") {
        return Err(CliError::usage(format!(
            "unsupported non-local path `{value}` for `{option}`"
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

#[cfg(test)]
mod tests;
