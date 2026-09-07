//! Methylation call orchestration.

use super::Options;
use super::output::{UnresolvedContextSummary, render_region};
use crate::call_input::prepare_call_input;
use crate::output::{create_text_output, finish_output};
use crate::region_workers::{IndexedCallMode, stream_indexed_region_workers_mode};
use crate::{CallError, CallReport};

pub(super) fn run(options: &Options) -> Result<CallReport, CallError> {
    let mode = IndexedCallMode::Meth(options.parameters);
    let mut output = create_text_output(
        "call meth",
        &options.output,
        &options.input,
        &options.reference,
        &[],
        options.compress,
        options.compression_threads,
    )?;
    let input = prepare_call_input(
        "call meth",
        &options.input,
        &options.reference,
        &options.regions,
        usize::try_from(options.threads).expect("validated thread count fits usize"),
        mode,
    )?;
    let mut summary = UnresolvedContextSummary::default();
    stream_indexed_region_workers_mode(
        &options.input,
        &input.references,
        &input.regions,
        input.worker_count,
        mode,
        &input.reference,
        |region| {
            let meth = region
                .meth
                .as_ref()
                .ok_or_else(|| CallError::operation("methylation region result is missing"))?;
            render_region(
                &mut output,
                options.format,
                options.parameters,
                &input.references,
                meth,
                &mut summary,
            )
            .map_err(|error| error.with_context("call meth: render methylation region"))
        },
    )?;
    let unresolved_warning = summary.into_warning("call meth");
    finish_output("call meth", output)?;
    Ok(CallReport::with_warning(unresolved_warning).with_prior_warning(input.reference_warning))
}
