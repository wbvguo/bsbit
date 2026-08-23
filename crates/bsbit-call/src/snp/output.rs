//! Ordered VCF rendering for SNP results.

use std::io::{self, Write};

use super::result::{SnpConfig, VariantCall};
use super::vcf::{render_vcf_call, render_vcf_header};
use crate::call_input::BamReference;
use crate::{CallError, CallErrorKind};

pub(crate) fn render_header(
    writer: &mut (impl Write + ?Sized),
    references: &[BamReference],
    config: SnpConfig,
    sample_name: &[u8],
) -> io::Result<()> {
    let dictionary = references
        .iter()
        .map(|reference| (reference.name.as_slice(), reference.length))
        .collect::<Vec<_>>();
    render_vcf_header(writer, &dictionary, config, sample_name)
}

pub(crate) fn render_region(
    writer: &mut (impl Write + ?Sized),
    references: &[BamReference],
    variants: &[(u32, VariantCall)],
) -> Result<(), CallError> {
    for (reference_id, call) in variants {
        let reference = references
            .get(usize::try_from(*reference_id).expect("u32 fits usize"))
            .ok_or_else(|| {
                CallError::operation("variant references a missing BAM dictionary entry")
            })?;
        render_vcf_call(writer, &reference.name, call).map_err(|error| {
            CallError::with_source(CallErrorKind::Output, "write VCF region", error)
        })?;
    }
    Ok(())
}
