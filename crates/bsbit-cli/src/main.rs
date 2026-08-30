//! Process entry point for the `bsbit` CLI.

use std::io;
use std::process::ExitCode;

fn main() -> ExitCode {
    let mut output = io::stdout().lock();
    match bsbit_cli::run(std::env::args_os().skip(1), &mut output) {
        Ok(report) => {
            for warning in report.warnings() {
                eprintln!("bsbit: warning: {}", warning.message());
            }
            ExitCode::SUCCESS
        }
        Err(error) => {
            eprintln!("bsbit: {error}");
            ExitCode::from(error.exit_code())
        }
    }
}
