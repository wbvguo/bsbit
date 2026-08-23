//! SNP call orchestration.

use super::Options;
use super::output::{render_header, render_region};
use super::result::SnpConfig;
use crate::call_input::{prepare_call_input, resolve_sample_name, validate_explicit_sample_name};
use crate::publication::{
    create_text_staging, finish_and_publish, output_write_error, publication_warning,
};
use crate::region_workers::{IndexedCallMode, stream_indexed_region_workers_mode};
use crate::{CallError, CallReport};

pub(super) fn run(options: &Options) -> Result<CallReport, CallError> {
    let config = SnpConfig::from(options.parameters);
    let mode = IndexedCallMode::Snp(config);
    validate_explicit_sample_name("call snp", options.sample_name.as_deref())?;
    let input = prepare_call_input(
        "call snp",
        &options.input,
        &options.reference,
        &options.regions,
        usize::try_from(options.threads).expect("validated thread count fits usize"),
        mode,
    )?;
    let sample_name = resolve_sample_name(
        "call snp",
        &options.input,
        options.sample_name.as_deref(),
        input.bam_sample_name.as_deref(),
    )?;
    let mut output = create_text_staging(
        "call snp",
        &options.output,
        options.compress,
        options.threads,
    )?;
    render_header(&mut output, &input.references, config, &sample_name)
        .map_err(|error| output_write_error("call snp", &options.output, error))?;
    stream_indexed_region_workers_mode(
        &options.input,
        &input.references,
        &input.regions,
        input.worker_count,
        mode,
        &input.reference,
        |region| render_region(&mut output, &input.references, &region.variants),
    )?;
    let publication = finish_and_publish("call snp", output)?;
    let warning = publication_warning(&publication, "SNP output");
    Ok(CallReport::with_warning(warning))
}
