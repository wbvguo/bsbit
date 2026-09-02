//! Methylation call orchestration.

use super::Options;
use super::output::{UnresolvedContextSummary, render_region};
use crate::call_input::prepare_call_input;
use crate::publication::{create_text_staging, finish_and_publish, publication_warning};
use crate::region::workers::{IndexedCallMode, stream_indexed_region_workers_mode};
use crate::{CallError, CallReport};

pub(super) fn run(options: &Options) -> Result<CallReport, CallError> {
    let mode = IndexedCallMode::Meth(options.parameters);
    let input = prepare_call_input(
        "call meth",
        &options.input,
        &options.reference,
        &options.regions,
        usize::try_from(options.threads).expect("validated thread count fits usize"),
        mode,
    )?;
    let mut output = create_text_staging(
        "call meth",
        &options.output,
        options.compress,
        options.threads,
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
                &input.references,
                meth,
                &mut summary,
            )
            .map_err(|error| error.with_context("call meth: render methylation region"))
        },
    )?;
    let unresolved_warning = summary.into_warning("call meth");
    let publication = finish_and_publish("call meth", output)?;
    let cleanup_warning = publication_warning(&publication, "methylation output");
    Ok(CallReport::with_warning(cleanup_warning).with_prior_warning(unresolved_warning))
}
