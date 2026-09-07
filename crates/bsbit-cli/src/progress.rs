//! Human-readable command progress written by the CLI coordinator.

use core::fmt;
use std::io::Write;
use std::time::{Duration, Instant};

use bsbit_cpu::Configuration;

pub(crate) const PROGRESS_INTERVAL: Duration = Duration::from_secs(5);

/// Coordinator-owned, rate-limited stdout progress for one command.
pub(crate) struct ProgressLog<'a> {
    output: &'a mut dyn Write,
    command: &'static str,
    started: Instant,
    last_progress: Instant,
    enabled: bool,
}

impl<'a> ProgressLog<'a> {
    pub(crate) fn start(
        output: &'a mut dyn Write,
        command: &'static str,
        configuration: Configuration,
        enabled: bool,
    ) -> Self {
        let started = Instant::now();
        let mut log = Self {
            output,
            command,
            started,
            last_progress: started,
            enabled,
        };
        if enabled
            && writeln!(
                log.output,
                "[bsbit::{command}] start architecture={} backend={} instruction_set={}",
                configuration.features().architecture(),
                configuration.backend(),
                configuration.backend().instruction_set(),
            )
            .and_then(|()| log.output.flush())
            .is_err()
        {
            log.enabled = false;
        }
        log
    }

    pub(crate) fn phase(&mut self, phase: &'static str, details: fmt::Arguments<'_>) {
        if !self.enabled {
            return;
        }
        if writeln!(
            self.output,
            "[bsbit::{}] phase={phase} {details} elapsed_seconds={:.3}",
            self.command,
            self.started.elapsed().as_secs_f64(),
        )
        .and_then(|()| self.output.flush())
        .is_err()
        {
            self.enabled = false;
        }
    }

    pub(crate) fn progress(
        &mut self,
        processed: u64,
        unit: &'static str,
        details: fmt::Arguments<'_>,
    ) {
        if !self.enabled || self.last_progress.elapsed() < PROGRESS_INTERVAL {
            return;
        }
        self.last_progress = Instant::now();
        self.write_count("progress", processed, unit, details);
    }

    pub(crate) fn complete(
        &mut self,
        processed: u64,
        unit: &'static str,
        details: fmt::Arguments<'_>,
    ) {
        if !self.enabled {
            return;
        }
        self.write_count("completed", processed, unit, details);
    }

    /// Writes one machine-readable line even when human progress is disabled.
    ///
    /// Metrics mode deliberately disables the human log stream, but it still
    /// shares this coordinator-owned writer so tests and embedding callers do
    /// not get bypassed through process-global stdout.
    pub(crate) fn machine_line(&mut self, line: fmt::Arguments<'_>) -> std::io::Result<()> {
        writeln!(self.output, "{line}")?;
        self.output.flush()
    }

    fn write_count(
        &mut self,
        event: &'static str,
        processed: u64,
        unit: &'static str,
        details: fmt::Arguments<'_>,
    ) {
        let elapsed = self.started.elapsed();
        let elapsed_millis = elapsed.as_millis().max(1);
        let rate = u128::from(processed).saturating_mul(1_000) / elapsed_millis;
        if writeln!(
            self.output,
            "[bsbit::{}] {event} {unit}={processed} {unit}_per_second={rate} {details} elapsed_seconds={:.3}",
            self.command,
            elapsed.as_secs_f64(),
        )
        .and_then(|()| self.output.flush())
        .is_err()
        {
            self.enabled = false;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bsbit_cpu::{BackendRequest, initialize};

    #[test]
    fn progress_is_rate_limited_but_completion_is_always_visible() {
        let configuration = initialize(BackendRequest::Auto).expect("CPU initializes");
        let mut output = Vec::new();
        {
            let mut log = ProgressLog::start(&mut output, "align", configuration, true);
            log.progress(10, "reads", format_args!("mapped=8"));
            log.complete(10, "reads", format_args!("mapped=8"));
        }

        let output = String::from_utf8(output).expect("progress is UTF-8");
        assert!(output.contains("[bsbit::align] start architecture="));
        assert!(output.contains(" instruction_set="));
        assert!(!output.contains("] progress "));
        assert!(output.contains("] completed reads=10 "));
        assert!(output.contains(" mapped=8 elapsed_seconds="));
    }
}
