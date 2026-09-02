use std::collections::BTreeMap;
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
    Action, option_map, optional_path, optional_u64, parse_bool, parse_threads,
    probability_parts_per_billion, required, required_path,
};

pub(super) fn parse_call(arguments: &[String]) -> Result<Action, CliError> {
    if matches!(arguments, [value] if value == "--help" || value == "-h" || value == "help") {
        return Ok(Action::Help(CALL_HELP));
    }
    let Some(module) = arguments.first().map(String::as_str) else {
        return Err(CliError::usage(
            "missing call module; run `bsbit call --help`",
        ));
    };
    match module {
        "meth" => parse_call_meth(&arguments[1..]),
        "snp" => parse_call_snp(&arguments[1..]),
        "joint" => parse_call_joint(&arguments[1..]),
        value => Err(CliError::usage(format!(
            "unknown call module `{value}`; run `bsbit call --help`"
        ))),
    }
}

fn parse_call_snp(arguments: &[String]) -> Result<Action, CliError> {
    if matches!(arguments, [value] if value == "--help" || value == "-h") {
        return Ok(Action::Help(CALL_SNP_HELP));
    }
    let normalized = normalize_call_snp_options(arguments)?;
    let (normalized, region_specs) = extract_repeatable_option(&normalized, "--region")?;
    let (mut values, _) = option_map(
        &normalized,
        &[
            "--input",
            "--reference",
            "--output",
            "--sample-name",
            "--regions-file",
            "--compress",
            "--threads",
            "--min-base-quality",
            "--min-mapq",
            "--min-depth",
            "--min-alt-count",
            "--min-alt-fraction",
            "--min-gq",
            "--min-aq",
            "--heterozygosity",
            "--underconversion-rate",
            "--overconversion-rate",
        ],
        &[],
    )?;
    let input = required_path(&mut values, "--input")?;
    let reference = required_path(&mut values, "--reference")?;
    let output = required_path(&mut values, "--output")?;
    let sample_name = values.remove("--sample-name");
    let regions = parse_call_regions(region_specs, optional_path(&mut values, "--regions-file")?)?;
    let compress = values
        .remove("--compress")
        .map(|value| parse_bool("--compress", &value))
        .transpose()?
        .unwrap_or(false);
    let threads = parse_threads(&mut values)?;
    let parameters = parse_snp_parameters(&mut values)?;
    Ok(Action::CallSnp(CallSnpOptions {
        input,
        reference,
        sample_name,
        regions,
        output,
        compress,
        threads,
        parameters,
    }))
}

fn parse_call_joint(arguments: &[String]) -> Result<Action, CliError> {
    if matches!(arguments, [value] if value == "--help" || value == "-h") {
        return Ok(Action::Help(CALL_JOINT_HELP));
    }
    let normalized = normalize_call_joint_options(arguments)?;
    let (normalized, region_specs) = extract_repeatable_option(&normalized, "--region")?;
    let (mut values, _) = option_map(
        &normalized,
        &[
            "--input",
            "--reference",
            "--sample-name",
            "--regions-file",
            "--meth",
            "--meth-format",
            "--vcf",
            "--compress",
            "--threads",
            "--min-base-quality",
            "--min-mapq",
            "--min-depth",
            "--min-alt-count",
            "--min-alt-fraction",
            "--min-gq",
            "--min-aq",
            "--heterozygosity",
            "--underconversion-rate",
            "--overconversion-rate",
        ],
        &[],
    )?;
    let input = required_path(&mut values, "--input")?;
    let reference = required_path(&mut values, "--reference")?;
    let sample_name = values.remove("--sample-name");
    let regions = parse_call_regions(region_specs, optional_path(&mut values, "--regions-file")?)?;
    let meth_output = required_path(&mut values, "--meth")?;
    let vcf_output = required_path(&mut values, "--vcf")?;
    if meth_output == vcf_output {
        return Err(CliError::usage(
            "--meth and --vcf must name different output files",
        ));
    }
    let meth_format = match required(&mut values, "--meth-format")?.as_str() {
        "cgmap" => MethylationOutputFormat::Cgmap,
        "bed" => MethylationOutputFormat::Bed,
        value => {
            return Err(CliError::usage(format!(
                "unsupported --meth-format `{value}`; expected `cgmap` or `bed`"
            )));
        }
    };
    let compress = values
        .remove("--compress")
        .map(|value| parse_bool("--compress", &value))
        .transpose()?
        .unwrap_or(false);
    let threads = parse_threads(&mut values)?;
    let parameters = parse_snp_parameters(&mut values)?;
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
        parameters,
    }))
}

fn parse_call_meth(arguments: &[String]) -> Result<Action, CliError> {
    if matches!(arguments, [value] if value == "--help" || value == "-h") {
        return Ok(Action::Help(CALL_METH_HELP));
    }
    let normalized = normalize_call_meth_options(arguments)?;
    let (normalized, region_specs) = extract_repeatable_option(&normalized, "--region")?;
    let (mut values, _) = option_map(
        &normalized,
        &[
            "--input",
            "--reference",
            "--regions-file",
            "--output",
            "--format",
            "--compress",
            "--threads",
            "--min-base-quality",
            "--min-mapq",
        ],
        &[],
    )?;
    let input = required_path(&mut values, "--input")?;
    let reference = required_path(&mut values, "--reference")?;
    let regions = parse_call_regions(region_specs, optional_path(&mut values, "--regions-file")?)?;
    let output = required_path(&mut values, "--output")?;
    let format = match required(&mut values, "--format")?.as_str() {
        "cgmap" => MethylationOutputFormat::Cgmap,
        "bed" => MethylationOutputFormat::Bed,
        value => {
            return Err(CliError::usage(format!(
                "unsupported --format `{value}`; expected `cgmap` or `bed`"
            )));
        }
    };
    let compress = values
        .remove("--compress")
        .map(|value| parse_bool("--compress", &value))
        .transpose()?
        .unwrap_or(false);
    let threads = parse_threads(&mut values)?;
    let parameters = parse_meth_parameters(&mut values)?;
    Ok(Action::CallMeth(CallMethOptions {
        input,
        reference,
        regions,
        output,
        format,
        compress,
        threads,
        parameters,
    }))
}

fn normalize_call_meth_options(arguments: &[String]) -> Result<Vec<String>, CliError> {
    arguments
        .iter()
        .map(|argument| {
            Ok(match argument.as_str() {
                "-i" => String::from("--input"),
                "-r" => String::from("--reference"),
                "-o" => String::from("--output"),
                "-f" => String::from("--format"),
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

fn normalize_call_snp_options(arguments: &[String]) -> Result<Vec<String>, CliError> {
    normalize_call_options(arguments, false)
}

fn normalize_call_joint_options(arguments: &[String]) -> Result<Vec<String>, CliError> {
    normalize_call_options(arguments, true)
}

fn normalize_call_options(arguments: &[String], joint: bool) -> Result<Vec<String>, CliError> {
    arguments
        .iter()
        .map(|argument| {
            Ok(match argument.as_str() {
                "-i" => String::from("--input"),
                "-r" => String::from("--reference"),
                "-o" if !joint => String::from("--output"),
                "-m" if joint => String::from("--meth"),
                "-v" if joint => String::from("--vcf"),
                "-f" if joint => String::from("--meth-format"),
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

fn extract_repeatable_option(
    arguments: &[String],
    repeatable: &str,
) -> Result<(Vec<String>, Vec<String>), CliError> {
    let mut remaining = Vec::with_capacity(arguments.len());
    let mut repeated = Vec::new();
    let mut cursor = 0;
    while cursor < arguments.len() {
        let option = &arguments[cursor];
        let Some(value) = arguments.get(cursor + 1) else {
            return Err(CliError::usage(format!("missing value for `{option}`")));
        };
        if !option.starts_with("--") {
            return Err(CliError::usage(format!("unknown option `{option}`")));
        }
        if value.starts_with("--") {
            return Err(CliError::usage(format!("missing value for `{option}`")));
        }
        if option == repeatable {
            repeated.push(value.clone());
        } else {
            remaining.push(option.clone());
            remaining.push(value.clone());
        }
        cursor += 2;
    }
    Ok((remaining, repeated))
}

fn parse_call_regions(
    specifications: Vec<String>,
    regions_file: Option<PathBuf>,
) -> Result<RegionSelection, CliError> {
    let intervals = specifications
        .into_iter()
        .map(|specification| parse_call_region(&specification))
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
    values: &mut BTreeMap<String, String>,
) -> Result<MethylationCallParameters, CliError> {
    let defaults = MethylationCallParameters::default();
    Ok(MethylationCallParameters {
        minimum_base_quality: bounded_u8(
            values,
            "--min-base-quality",
            defaults.minimum_base_quality,
            93,
        )?,
        minimum_mapping_quality: bounded_u8(
            values,
            "--min-mapq",
            defaults.minimum_mapping_quality,
            254,
        )?,
    })
}

fn parse_snp_parameters(
    values: &mut BTreeMap<String, String>,
) -> Result<SnpCallParameters, CliError> {
    let defaults = SnpCallParameters::default();
    let minimum_base_quality = bounded_u8(
        values,
        "--min-base-quality",
        defaults.minimum_base_quality,
        93,
    )?;
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
    values: &mut BTreeMap<String, String>,
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

fn nonzero_u32(
    values: &mut BTreeMap<String, String>,
    option: &str,
    default: u32,
) -> Result<u32, CliError> {
    let value = optional_u64(values, option)?.unwrap_or(u64::from(default));
    u32::try_from(value)
        .ok()
        .filter(|value| *value != 0)
        .ok_or_else(|| CliError::usage(format!("{option} must be in 1..=4294967295")))
}
