//! CPU capability and SIMD backend diagnostics.

use std::ffi::OsString;
use std::io::Write;

use bsbit_cpu::{BackendRequest, initialize};

use crate::{CPU_HELP, CliError};

/// Validated CPU diagnostic options.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub(crate) struct Options {
    request: BackendRequest,
}

pub(crate) fn parse(arguments: &[OsString]) -> Result<super::Action, CliError> {
    if arguments.is_empty() {
        return Ok(super::Action::Cpu(Options::default()));
    }
    if matches!(arguments, [value] if value == "--help" || value == "-h") {
        return Ok(super::Action::Help(CPU_HELP));
    }
    let [flag, value] = arguments else {
        return Err(CliError::usage(
            "usage: bsbit cpu [--simd-backend auto|scalar|sse2|sse4.2|avx2|avx512|neon]",
        ));
    };
    let flag = flag
        .to_str()
        .ok_or_else(|| CliError::usage("option name must be valid UTF-8"))?;
    if flag != "--simd-backend" {
        return Err(CliError::usage(format!("unknown option `{flag}`")));
    }
    let value = value
        .to_str()
        .ok_or_else(|| CliError::usage("--simd-backend value must be valid UTF-8"))?;
    let request = value
        .parse()
        .map_err(|error| CliError::usage(format!("invalid --simd-backend `{value}`: {error}")))?;
    Ok(super::Action::Cpu(Options { request }))
}

pub(crate) fn run(options: Options, output: &mut impl Write) -> Result<(), CliError> {
    let configuration = initialize(options.request)
        .map_err(|error| CliError::operation(format!("CPU backend selection: {error}")))?;
    let features = configuration.features();
    writeln!(output, "architecture={}", features.architecture())
        .and_then(|()| writeln!(output, "sse2={}", u8::from(features.sse2())))
        .and_then(|()| writeln!(output, "sse4.1={}", u8::from(features.sse41())))
        .and_then(|()| writeln!(output, "sse4.2={}", u8::from(features.sse42())))
        .and_then(|()| writeln!(output, "popcnt={}", u8::from(features.popcnt())))
        .and_then(|()| writeln!(output, "avx2={}", u8::from(features.avx2())))
        .and_then(|()| writeln!(output, "avx512f={}", u8::from(features.avx512f())))
        .and_then(|()| writeln!(output, "avx512bw={}", u8::from(features.avx512bw())))
        .and_then(|()| writeln!(output, "neon={}", u8::from(features.neon())))
        .and_then(|()| writeln!(output, "backend={}", configuration.backend()))
        .and_then(|()| {
            writeln!(
                output,
                "instruction_set={}",
                configuration.backend().instruction_set()
            )
        })
        .map_err(|error| CliError::operation(format!("write CPU diagnostics: {error}")))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parser_accepts_default_and_stable_backend_spellings() {
        assert!(matches!(parse(&[]), Ok(super::super::Action::Cpu(_))));
        for value in ["auto", "scalar", "sse2", "sse4.2", "avx2", "avx512", "neon"] {
            assert!(matches!(
                parse(&[OsString::from("--simd-backend"), OsString::from(value)]),
                Ok(super::super::Action::Cpu(_))
            ));
        }
        assert!(parse(&[OsString::from("--simd-backend"), OsString::from("native"),]).is_err());
    }

    #[test]
    fn report_contains_independent_features_and_selected_backend() {
        let mut output = Vec::new();
        run(Options::default(), &mut output).expect("automatic CPU report");
        let report = String::from_utf8(output).expect("ASCII report");
        let architecture = format!("architecture={}\n", std::env::consts::ARCH);
        assert!(
            report.contains(&architecture),
            "missing {architecture:?} from {report:?}"
        );
        for key in [
            "sse2=",
            "sse4.1=",
            "sse4.2=",
            "popcnt=",
            "avx2=",
            "avx512f=",
            "avx512bw=",
            "neon=",
            "backend=",
            "instruction_set=",
        ] {
            assert!(report.contains(key), "missing {key:?} from {report:?}");
        }
    }
}
