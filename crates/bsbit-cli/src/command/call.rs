use std::ffi::OsString;
use std::path::PathBuf;

use bsbit_call::joint::Options as CallJointOptions;
use bsbit_call::meth::{
    Options as CallMethOptions, OutputFormat as MethylationOutputFormat,
    Parameters as MethylationCallParameters,
};
use bsbit_call::region::{GenomicInterval, RegionSelection};
use bsbit_call::snp::{Options as CallSnpOptions, Parameters as SnpCallParameters};

use crate::{CALL_HELP, CALL_JOINT_HELP, CALL_METH_HELP, CALL_SNP_HELP, CliError};

use super::{
    Action, OptionValues, option_map, optional_path, optional_text, optional_u64, parse_bool,
    parse_compression_threads, parse_threads, probability_parts_per_billion, required_path,
    required_text, text_value,
};

const CALL_METH_VALUE_OPTIONS: &[&str] = &[
    "--input",
    "--reference",
    "--regions-bed",
    "--output",
    "--format",
    "--compress",
    "--threads",
    "--compression-threads",
    "--min-bq",
    "--min-mapq",
    "--min-depth",
];
const CALL_METH_FLAGS: &[&str] = &["--cg-only", "--ignore-orphan"];

const CALL_SNP_VALUE_OPTIONS: &[&str] = &[
    "--input",
    "--reference",
    "--output",
    "--sample-name",
    "--regions-bed",
    "--compress",
    "--threads",
    "--compression-threads",
    "--min-bq",
    "--min-mapq",
    "--min-depth",
    "--min-alt-count",
    "--min-alt-fraction",
    "--min-gq",
    "--min-aq",
    "--heterozygosity",
    "--underconversion-rate",
    "--overconversion-rate",
];
const CALL_SNP_FLAGS: &[&str] = &["--ignore-orphan"];

const CALL_JOINT_VALUE_OPTIONS: &[&str] = &[
    "--input",
    "--reference",
    "--sample-name",
    "--regions-bed",
    "--prefix",
    "--meth-format",
    "--compress",
    "--threads",
    "--compression-threads",
    "--min-bq",
    "--min-mapq",
    "--min-depth",
    "--min-alt-count",
    "--min-alt-fraction",
    "--min-gq",
    "--min-aq",
    "--heterozygosity",
    "--underconversion-rate",
    "--overconversion-rate",
];
const CALL_JOINT_FLAGS: &[&str] = &["--cg-only", "--ignore-orphan"];

pub(super) fn parse_call(arguments: &[OsString]) -> Result<Action, CliError> {
    if matches!(arguments, [value] if value == "--help" || value == "-h" || value == "help") {
        return Ok(Action::Help(CALL_HELP));
    }
    let Some(module) = arguments.first() else {
        return Err(CliError::usage(
            "missing call module; run `bsbit call --help`",
        ));
    };
    let module = module
        .to_str()
        .ok_or_else(|| CliError::usage("call module name must be valid UTF-8"))?;
    match module {
        "meth" => parse_call_meth(&arguments[1..]),
        "snp" => parse_call_snp(&arguments[1..]),
        "joint" => parse_call_joint(&arguments[1..]),
        value => Err(CliError::usage(format!(
            "unknown call module `{value}`; run `bsbit call --help`"
        ))),
    }
}

fn parse_call_snp(arguments: &[OsString]) -> Result<Action, CliError> {
    if matches!(arguments, [value] if value == "--help" || value == "-h") {
        return Ok(Action::Help(CALL_SNP_HELP));
    }
    let normalized = normalize_call_snp_options(arguments)?;
    let (normalized, region_specs) =
        extract_repeatable_option(&normalized, "--region", CALL_SNP_FLAGS)?;
    let (mut values, flags) = option_map(&normalized, CALL_SNP_VALUE_OPTIONS, CALL_SNP_FLAGS)?;
    let input = required_path(&mut values, "--input")?;
    let reference = required_path(&mut values, "--reference")?;
    let output = required_path(&mut values, "--output")?;
    let sample_name = optional_text(&mut values, "--sample-name")?;
    let regions = parse_call_regions(region_specs, optional_path(&mut values, "--regions-bed")?)?;
    let compress = optional_text(&mut values, "--compress")?
        .map(|value| parse_bool("--compress", &value))
        .transpose()?
        .unwrap_or(false);
    let threads = parse_threads(&mut values)?;
    let compression_threads = parse_compression_threads(&mut values, compress, threads)?;
    let parameters = parse_snp_parameters(&mut values, flags.contains("--ignore-orphan"))?;
    Ok(Action::CallSnp(CallSnpOptions {
        input,
        reference,
        sample_name,
        regions,
        output,
        compress,
        threads,
        compression_threads,
        parameters,
    }))
}

fn parse_call_joint(arguments: &[OsString]) -> Result<Action, CliError> {
    if matches!(arguments, [value] if value == "--help" || value == "-h") {
        return Ok(Action::Help(CALL_JOINT_HELP));
    }
    let normalized = normalize_call_joint_options(arguments)?;
    let (normalized, region_specs) =
        extract_repeatable_option(&normalized, "--region", CALL_JOINT_FLAGS)?;
    let (mut values, flags) = option_map(&normalized, CALL_JOINT_VALUE_OPTIONS, CALL_JOINT_FLAGS)?;
    let input = required_path(&mut values, "--input")?;
    let reference = required_path(&mut values, "--reference")?;
    let sample_name = optional_text(&mut values, "--sample-name")?;
    let regions = parse_call_regions(region_specs, optional_path(&mut values, "--regions-bed")?)?;
    let prefix = required_path(&mut values, "--prefix")?;
    let meth_format = match required_text(&mut values, "--meth-format")?.as_str() {
        "cgmap" => MethylationOutputFormat::Cgmap,
        "bed" => MethylationOutputFormat::Bed,
        value => {
            return Err(CliError::usage(format!(
                "unsupported --meth-format `{value}`; expected `cgmap` or `bed`"
            )));
        }
    };
    let compress = optional_text(&mut values, "--compress")?
        .map(|value| parse_bool("--compress", &value))
        .transpose()?
        .unwrap_or(false);
    let meth_output = joint_output_path(
        &prefix,
        match meth_format {
            MethylationOutputFormat::Cgmap => ".CGmap",
            MethylationOutputFormat::Bed => ".bed",
        },
        compress,
    );
    let vcf_output = joint_output_path(&prefix, ".vcf", compress);
    let threads = parse_threads(&mut values)?;
    let compression_threads = parse_compression_threads(&mut values, compress, threads)?;
    let parameters = parse_snp_parameters(&mut values, flags.contains("--ignore-orphan"))?;
    Ok(Action::CallJoint(CallJointOptions {
        input,
        reference,
        sample_name,
        regions,
        meth_output,
        meth_format,
        vcf_output,
        compress,
        threads,
        compression_threads,
        cg_only: flags.contains("--cg-only"),
        parameters,
    }))
}

fn parse_call_meth(arguments: &[OsString]) -> Result<Action, CliError> {
    if matches!(arguments, [value] if value == "--help" || value == "-h") {
        return Ok(Action::Help(CALL_METH_HELP));
    }
    let normalized = normalize_call_meth_options(arguments)?;
    let (normalized, region_specs) =
        extract_repeatable_option(&normalized, "--region", CALL_METH_FLAGS)?;
    let (mut values, flags) = option_map(&normalized, CALL_METH_VALUE_OPTIONS, CALL_METH_FLAGS)?;
    let input = required_path(&mut values, "--input")?;
    let reference = required_path(&mut values, "--reference")?;
    let regions = parse_call_regions(region_specs, optional_path(&mut values, "--regions-bed")?)?;
    let output = required_path(&mut values, "--output")?;
    let format = match required_text(&mut values, "--format")?.as_str() {
        "cgmap" => MethylationOutputFormat::Cgmap,
        "bed" => MethylationOutputFormat::Bed,
        value => {
            return Err(CliError::usage(format!(
                "unsupported --format `{value}`; expected `cgmap` or `bed`"
            )));
        }
    };
    let compress = optional_text(&mut values, "--compress")?
        .map(|value| parse_bool("--compress", &value))
        .transpose()?
        .unwrap_or(false);
    let threads = parse_threads(&mut values)?;
    let compression_threads = parse_compression_threads(&mut values, compress, threads)?;
    let parameters = parse_meth_parameters(
        &mut values,
        flags.contains("--cg-only"),
        flags.contains("--ignore-orphan"),
    )?;
    Ok(Action::CallMeth(CallMethOptions {
        input,
        reference,
        regions,
        output,
        format,
        compress,
        threads,
        compression_threads,
        parameters,
    }))
}

fn normalize_call_meth_options(arguments: &[OsString]) -> Result<Vec<OsString>, CliError> {
    arguments
        .iter()
        .map(|argument| {
            Ok(match argument.to_str() {
                Some("-i") => OsString::from("--input"),
                Some("-r") => OsString::from("--reference"),
                Some("-o") => OsString::from("--output"),
                Some("-f") => OsString::from("--format"),
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

fn normalize_call_snp_options(arguments: &[OsString]) -> Result<Vec<OsString>, CliError> {
    normalize_call_options(arguments, false)
}

fn normalize_call_joint_options(arguments: &[OsString]) -> Result<Vec<OsString>, CliError> {
    normalize_call_options(arguments, true)
}

fn normalize_call_options(arguments: &[OsString], joint: bool) -> Result<Vec<OsString>, CliError> {
    arguments
        .iter()
        .map(|argument| {
            Ok(match argument.to_str() {
                Some("-i") => OsString::from("--input"),
                Some("-r") => OsString::from("--reference"),
                Some("-o") if !joint => OsString::from("--output"),
                Some("-p") if joint => OsString::from("--prefix"),
                Some("-f") if joint => OsString::from("--meth-format"),
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

fn joint_output_path(prefix: &std::path::Path, suffix: &str, compress: bool) -> PathBuf {
    let mut output = prefix.as_os_str().to_os_string();
    output.push(suffix);
    if compress {
        output.push(".gz");
    }
    output.into()
}

fn extract_repeatable_option(
    arguments: &[OsString],
    repeatable: &str,
    accepted_flags: &[&str],
) -> Result<(Vec<OsString>, Vec<OsString>), CliError> {
    let mut remaining = Vec::with_capacity(arguments.len());
    let mut repeated = Vec::new();
    let mut cursor = 0;
    while cursor < arguments.len() {
        let option = arguments[cursor]
            .to_str()
            .ok_or_else(|| CliError::usage("option name must be valid UTF-8"))?;
        if accepted_flags.contains(&option) {
            remaining.push(arguments[cursor].clone());
            cursor += 1;
            continue;
        }
        let Some(value) = arguments.get(cursor + 1) else {
            return Err(CliError::usage(format!("missing value for `{option}`")));
        };
        if !option.starts_with("--") {
            return Err(CliError::usage(format!("unknown option `{option}`")));
        }
        if value.to_str().is_some_and(|value| value.starts_with("--")) {
            return Err(CliError::usage(format!("missing value for `{option}`")));
        }
        if option == repeatable {
            repeated.push(value.clone());
        } else {
            remaining.push(arguments[cursor].clone());
            remaining.push(value.clone());
        }
        cursor += 2;
    }
    Ok((remaining, repeated))
}

fn parse_call_regions(
    specifications: Vec<OsString>,
    regions_file: Option<PathBuf>,
) -> Result<RegionSelection, CliError> {
    let intervals = specifications
        .into_iter()
        .map(|specification| text_value("--region", &specification).and_then(parse_call_region))
        .collect::<Result<Vec<_>, _>>()?;
    Ok(RegionSelection {
        intervals,
        regions_file,
    })
}

fn parse_call_region(specification: &str) -> Result<GenomicInterval, CliError> {
    let Some((contig, coordinates)) = specification.rsplit_once(':') else {
        return Err(invalid_call_region(specification));
    };
    let Some((start, end)) = coordinates.split_once('-') else {
        return Err(invalid_call_region(specification));
    };
    if contig.is_empty()
        || contig
            .bytes()
            .any(|byte| byte.is_ascii_control() || byte.is_ascii_whitespace())
    {
        return Err(invalid_call_region(specification));
    }
    let start = parse_one_based_region_coordinate(specification, start)?;
    let end = parse_one_based_region_coordinate(specification, end)?;
    if start == 0 || start > end {
        return Err(invalid_call_region(specification));
    }
    Ok(GenomicInterval {
        contig: contig.to_owned(),
        start: start - 1,
        end,
    })
}

fn parse_one_based_region_coordinate(specification: &str, value: &str) -> Result<u64, CliError> {
    let groups = value.split(',').collect::<Vec<_>>();
    let all_digits = groups
        .iter()
        .all(|group| !group.is_empty() && group.bytes().all(|byte| byte.is_ascii_digit()));
    let valid_grouping = all_digits
        && (groups.len() == 1
            || (groups.first().is_some_and(|group| group.len() <= 3)
                && groups.iter().skip(1).all(|group| group.len() == 3)));
    if !valid_grouping {
        return Err(invalid_call_region(specification));
    }
    value
        .bytes()
        .filter(|byte| *byte != b',')
        .try_fold(0_u64, |coordinate, digit| {
            coordinate
                .checked_mul(10)?
                .checked_add(u64::from(digit - b'0'))
        })
        .ok_or_else(|| invalid_call_region(specification))
}

fn invalid_call_region(specification: &str) -> CliError {
    CliError::usage(format!(
        "invalid --region `{specification}`; expected CONTIG:START-END with 1-based inclusive coordinates"
    ))
}

fn parse_meth_parameters(
    values: &mut OptionValues,
    cg_only: bool,
    ignore_orphans: bool,
) -> Result<MethylationCallParameters, CliError> {
    let defaults = MethylationCallParameters::default();
    Ok(MethylationCallParameters {
        minimum_base_quality: bounded_u8(values, "--min-bq", defaults.minimum_base_quality, 93)?,
        minimum_mapping_quality: bounded_u8(
            values,
            "--min-mapq",
            defaults.minimum_mapping_quality,
            254,
        )?,
        minimum_depth: nonzero_u32(values, "--min-depth", defaults.minimum_depth)?,
        cg_only,
        ignore_orphans,
    })
}

fn parse_snp_parameters(
    values: &mut OptionValues,
    ignore_orphans: bool,
) -> Result<SnpCallParameters, CliError> {
    let defaults = SnpCallParameters::default();
    let minimum_base_quality = bounded_u8(values, "--min-bq", defaults.minimum_base_quality, 93)?;
    let minimum_mapping_quality =
        bounded_u8(values, "--min-mapq", defaults.minimum_mapping_quality, 254)?;
    let minimum_depth = nonzero_u32(values, "--min-depth", defaults.minimum_depth)?;
    let minimum_alternate_count =
        nonzero_u32(values, "--min-alt-count", defaults.minimum_alternate_count)?;
    let minimum_alternate_fraction_parts_per_billion = probability_parts_per_billion(
        values,
        "--min-alt-fraction",
        defaults.minimum_alternate_fraction_parts_per_billion,
    )?;
    let minimum_genotype_quality =
        bounded_u8(values, "--min-gq", defaults.minimum_genotype_quality, 99)?;
    let minimum_allele_quality =
        bounded_u8(values, "--min-aq", defaults.minimum_allele_quality, 99)?;
    let heterozygosity_parts_per_billion = probability_parts_per_billion(
        values,
        "--heterozygosity",
        defaults.heterozygosity_parts_per_billion,
    )?;
    if heterozygosity_parts_per_billion == 0 || heterozygosity_parts_per_billion >= 1_000_000_000 {
        return Err(CliError::usage(
            "--heterozygosity must be strictly between 0 and 1",
        ));
    }
    let underconversion_parts_per_billion = probability_parts_per_billion(
        values,
        "--underconversion-rate",
        defaults.underconversion_parts_per_billion,
    )?;
    let overconversion_parts_per_billion = probability_parts_per_billion(
        values,
        "--overconversion-rate",
        defaults.overconversion_parts_per_billion,
    )?;
    Ok(SnpCallParameters {
        minimum_base_quality,
        minimum_mapping_quality,
        ignore_orphans,
        minimum_depth,
        minimum_alternate_count,
        minimum_alternate_fraction_parts_per_billion,
        minimum_genotype_quality,
        minimum_allele_quality,
        heterozygosity_parts_per_billion,
        underconversion_parts_per_billion,
        overconversion_parts_per_billion,
    })
}

fn bounded_u8(
    values: &mut OptionValues,
    option: &str,
    default: u8,
    maximum: u8,
) -> Result<u8, CliError> {
    let value = optional_u64(values, option)?.unwrap_or(u64::from(default));
    u8::try_from(value)
        .ok()
        .filter(|value| *value <= maximum)
        .ok_or_else(|| CliError::usage(format!("{option} must be in 0..={maximum}")))
}

fn nonzero_u32(values: &mut OptionValues, option: &str, default: u32) -> Result<u32, CliError> {
    let value = optional_u64(values, option)?.unwrap_or(u64::from(default));
    u32::try_from(value)
        .ok()
        .filter(|value| *value != 0)
        .ok_or_else(|| CliError::usage(format!("{option} must be in 1..=4294967295")))
}

#[cfg(test)]
mod tests {
    use crate::{CALL_JOINT_HELP, CALL_METH_HELP, CALL_SNP_HELP};

    use super::{
        CALL_JOINT_FLAGS, CALL_JOINT_VALUE_OPTIONS, CALL_METH_FLAGS, CALL_METH_VALUE_OPTIONS,
        CALL_SNP_FLAGS, CALL_SNP_VALUE_OPTIONS,
    };

    #[test]
    fn module_help_names_every_accepted_long_option() {
        for (help, values, flags) in [
            (CALL_METH_HELP, CALL_METH_VALUE_OPTIONS, CALL_METH_FLAGS),
            (CALL_SNP_HELP, CALL_SNP_VALUE_OPTIONS, CALL_SNP_FLAGS),
            (CALL_JOINT_HELP, CALL_JOINT_VALUE_OPTIONS, CALL_JOINT_FLAGS),
        ] {
            for option in values.iter().chain(flags).chain([&"--region"]) {
                assert!(help.contains(option), "module help omits {option}");
            }
        }
    }
}
