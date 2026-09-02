//! Joint methylation and SNP call orchestration.

use bsbit_io::validate_distinct_paths;

use super::Options;
use crate::call_input::{prepare_call_input, resolve_sample_name, validate_explicit_sample_name};
use crate::meth::output::{UnresolvedContextSummary, render_region as render_meth_region};
use crate::publication::{create_text_staging, output_write_error, publication_warning};
use crate::region::workers::{IndexedCallMode, stream_indexed_region_workers_mode};
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
    let mode = IndexedCallMode::Joint(config);
    validate_explicit_sample_name("call joint", options.sample_name.as_deref())?;
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
    let mut meth_output = create_text_staging(
        "call joint",
        &options.meth_output,
        options.compress,
        options.threads,
    )?;
    let mut vcf_output = create_text_staging(
        "call joint",
        &options.vcf_output,
        options.compress,
        options.threads,
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
                &input.references,
                meth,
                &mut summary,
            )
            .map_err(|error| error.with_context("call joint: render methylation region"))?;
            render_vcf_region(&mut vcf_output, &input.references, &region.variants)
        },
    )?;
    let unresolved_warning = summary.into_warning("call joint");
    let meth_completed = meth_output.finish().map_err(|error| {
        CallError::with_source(
            CallErrorKind::Output,
            "call joint: finalize methylation output",
            error,
        )
    })?;
    let vcf_completed = vcf_output.finish().map_err(|error| {
        CallError::with_source(
            CallErrorKind::Output,
            "call joint: finalize SNP output",
            error,
        )
    })?;
    let meth_publication = meth_completed.publish_create_new().map_err(|error| {
        CallError::with_source(
            CallErrorKind::Publication,
            "call joint: publish methylation output",
            error,
        )
    })?;
    let vcf_publication = match vcf_completed.publish_create_new() {
        Ok(publication) => publication,
        Err(error) => {
            meth_publication.rollback().map_err(|rollback| {
                CallError::with_source(
                    CallErrorKind::Publication,
                    format!(
                        "call joint: publish SNP output failed ({error}); methylation rollback failed"
                    ),
                    rollback,
                )
            })?;
            return Err(CallError::with_source(
                CallErrorKind::Publication,
                "call joint: publish SNP output",
                error,
            ));
        }
    };
    let meth_warning = publication_warning(&meth_publication, "methylation output");
    let vcf_warning = publication_warning(&vcf_publication, "SNP output");
    Ok(CallReport::with_warning(vcf_warning)
        .with_prior_warning(meth_warning)
        .with_prior_warning(unresolved_warning))
}
