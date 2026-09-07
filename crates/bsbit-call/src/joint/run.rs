//! Joint methylation and SNP call orchestration.

use bsbit_io::validate_distinct_paths;

use super::Options;
use crate::call_input::{prepare_call_input, resolve_sample_name, validate_explicit_sample_name};
use crate::meth::Parameters as MethParameters;
use crate::meth::output::{UnresolvedContextSummary, render_region as render_meth_region};
use crate::output::{create_text_output, finish_output, output_write_error};
use crate::region_workers::{IndexedCallMode, stream_indexed_region_workers_mode};
use crate::snp::output::{render_header as render_vcf_header, render_region as render_vcf_region};
use crate::snp::result::SnpConfig;
use crate::{CallError, CallErrorKind, CallReport};

#[allow(clippy::too_many_lines)]
pub(super) fn run(options: &Options) -> Result<CallReport, CallError> {
    validate_distinct_paths(&options.meth_output, &options.vcf_output).map_err(|error| {
        CallError::with_source(
            CallErrorKind::Configuration,
            "call joint: output paths must differ",
            error,
        )
    })?;
    let config = SnpConfig::from(options.parameters);
    let meth_parameters = MethParameters {
        minimum_base_quality: options.parameters.minimum_base_quality,
        minimum_mapping_quality: options.parameters.minimum_mapping_quality,
        minimum_depth: options.parameters.minimum_depth,
        cg_only: options.cg_only,
        ignore_orphans: options.parameters.ignore_orphans,
    };
    let mode = IndexedCallMode::Joint(config);
    validate_explicit_sample_name("call joint", options.sample_name.as_deref())?;
    let mut meth_output = create_text_output(
        "call joint",
        &options.meth_output,
        &options.input,
        &options.reference,
        &[&options.vcf_output],
        options.compress,
        options.compression_threads,
    )?;
    let mut vcf_output = create_text_output(
        "call joint",
        &options.vcf_output,
        &options.input,
        &options.reference,
        &[&options.meth_output],
        options.compress,
        options.compression_threads,
    )?;
    let input = prepare_call_input(
        "call joint",
        &options.input,
        &options.reference,
        &options.regions,
        usize::try_from(options.threads).expect("validated thread count fits usize"),
        mode,
    )?;
    let sample_name = resolve_sample_name(
        "call joint",
        &options.input,
        options.sample_name.as_deref(),
        input.bam_sample_name.as_deref(),
    )?;
    render_vcf_header(&mut vcf_output, &input.references, config, &sample_name)
        .map_err(|error| output_write_error("call joint", &options.vcf_output, error))?;
    let mut summary = UnresolvedContextSummary::default();
    stream_indexed_region_workers_mode(
        &options.input,
        &input.references,
        &input.regions,
        input.worker_count,
        mode,
        &input.reference,
        |region| {
            let meth = region.meth.as_ref().ok_or_else(|| {
                CallError::operation("joint methylation region result is missing")
            })?;
            render_meth_region(
                &mut meth_output,
                options.meth_format,
                meth_parameters,
                &input.references,
                meth,
                &mut summary,
            )
            .map_err(|error| error.with_context("call joint: render methylation region"))?;
            render_vcf_region(&mut vcf_output, &input.references, &region.variants)
        },
    )?;
    let unresolved_warning = summary.into_warning("call joint");
    finish_output("call joint methylation", meth_output)?;
    finish_output("call joint SNP", vcf_output)?;
    Ok(CallReport::with_warning(unresolved_warning).with_prior_warning(input.reference_warning))
}
