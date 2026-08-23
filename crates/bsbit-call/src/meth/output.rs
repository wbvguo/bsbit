//! Ordered `CGmap` and `bedMethyl` rendering.

use std::io::{self, Write};

use bsbit_hts::{BedMethylContext, BedMethylRecord, BedMethylStrand};

use super::OutputFormat;
use super::aggregation::{DenseMethRegion, SiteCounts, SiteKey};
use crate::call_input::BamReference;
use crate::evidence::{ContextClass, CytosineContext, EvidenceStrand};
use crate::{CallError, CallErrorKind, CallWarning};

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub(crate) struct UnresolvedContextSummary {
    sites: u64,
    observations: u64,
}

impl UnresolvedContextSummary {
    fn observe(&mut self, counts: &SiteCounts, coverage: u64) -> Result<(), CallError> {
        if counts.context.is_some() || coverage == 0 {
            return Ok(());
        }
        self.sites = self
            .sites
            .checked_add(1)
            .ok_or_else(|| CallError::operation("unresolved-context site count overflowed"))?;
        self.observations = self.observations.checked_add(coverage).ok_or_else(|| {
            CallError::operation("unresolved-context observation count overflowed")
        })?;
        Ok(())
    }

    pub(crate) fn into_warning(self, command: &str) -> Option<CallWarning> {
        (self.sites != 0).then(|| {
            CallWarning::new(format!(
                "{command} omitted {} observations at {} covered sites because two canonical neighboring reference bases were unavailable for context classification",
                self.observations, self.sites
            ))
        })
    }
}

pub(crate) fn render_region(
    writer: &mut (impl Write + ?Sized),
    format: OutputFormat,
    references: &[BamReference],
    region: &DenseMethRegion,
    unresolved: &mut UnresolvedContextSummary,
) -> Result<(), CallError> {
    region.for_each_site(|key, counts| {
        let coverage = counts.valid_coverage()?;
        unresolved.observe(&counts, coverage)?;
        let Some(context) = counts.context else {
            return Ok(());
        };
        if coverage == 0 {
            return Ok(());
        }
        let reference = references
            .get(usize::try_from(key.reference).expect("u32 fits usize"))
            .ok_or_else(|| {
                CallError::operation("site references a missing BAM dictionary entry")
            })?;
        match format {
            OutputFormat::Cgmap => {
                render_cgmap(writer, &reference.name, key, &counts, context, coverage).map_err(
                    |error| {
                        CallError::with_source(CallErrorKind::Output, "write CGmap region", error)
                    },
                )?;
            }
            OutputFormat::Bed => {
                render_bed(writer, &reference.name, key, &counts, context, coverage).map_err(
                    |error| {
                        CallError::with_source(
                            CallErrorKind::Output,
                            "write bedMethyl region",
                            error,
                        )
                    },
                )?;
            }
        }
        Ok(())
    })
}

pub(crate) fn render_cgmap(
    writer: &mut (impl Write + ?Sized),
    reference: &[u8],
    key: SiteKey,
    counts: &SiteCounts,
    context: CytosineContext,
    coverage: u64,
) -> io::Result<()> {
    writer.write_all(reference)?;
    let nucleotide = match key.strand {
        EvidenceStrand::Top => b'C',
        EvidenceStrand::Bottom => b'G',
    };
    let level = rounded_scaled_ratio(counts.methylated, coverage, 1_000_000);
    let level_whole = level / 1_000_000;
    let level_fraction = level % 1_000_000;
    write!(
        writer,
        "\t{}\t{}\t",
        char::from(nucleotide),
        u64::from(key.position) + 1
    )?;
    writer.write_all(context.class.name())?;
    writer.write_all(b"\tC")?;
    writer.write_all(&[context.second])?;
    writeln!(
        writer,
        "\t{level_whole}.{level_fraction:06}\t{}\t{coverage}",
        counts.methylated
    )
}

pub(crate) fn render_bed(
    writer: &mut (impl Write + ?Sized),
    reference: &[u8],
    key: SiteKey,
    counts: &SiteCounts,
    context: CytosineContext,
    coverage: u64,
) -> io::Result<()> {
    let strand = match key.strand {
        EvidenceStrand::Top => BedMethylStrand::Forward,
        EvidenceStrand::Bottom => BedMethylStrand::Reverse,
    };
    let context = match context.class {
        ContextClass::Cg => BedMethylContext::Cg,
        ContextClass::Chg => BedMethylContext::Chg,
        ContextClass::Chh => BedMethylContext::Chh,
    };
    let record = BedMethylRecord::new(
        reference,
        u64::from(key.position),
        context,
        strand,
        b"255,0,0",
        counts.methylated,
        counts.unmethylated,
        0,
        counts.deleted,
        0,
        counts.different,
        0,
    )
    .map_err(io::Error::other)?;
    debug_assert_eq!(record.coverage(), coverage);
    record.encode(writer)
}

fn rounded_scaled_ratio(numerator: u64, denominator: u64, scale: u64) -> u64 {
    debug_assert!(denominator != 0);
    debug_assert!(numerator <= denominator);
    let denominator = u128::from(denominator);
    let scaled = (u128::from(numerator) * u128::from(scale) + denominator / 2) / denominator;
    u64::try_from(scaled).expect("bounded methylation ratio fits u64")
}
