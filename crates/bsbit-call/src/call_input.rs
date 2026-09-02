//! Shared BAM, reference, region, and sample preflight for call modes.

use std::collections::BTreeSet;
use std::path::Path;

use bsbit_core::reference::ReferenceSemanticDigest;
use bsbit_hts::IndexedBamReader;

use crate::reference_context::CallReferenceSource;
use crate::region::{CallRegion, RegionSelection, plan_call_regions};
use crate::region_workers::{IndexedCallMode, region_bases_for};
use crate::{CallError, CallErrorKind};

/// One validated entry from the input BAM reference dictionary.
#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct BamReference {
    pub(crate) name: Vec<u8>,
    pub(crate) length: u32,
}

pub(crate) fn resolve_sample_name(
    command: &str,
    input: &Path,
    explicit: Option<&str>,
    bam_sample: Option<&[u8]>,
) -> Result<Vec<u8>, CallError> {
    if let Some(sample) = explicit {
        return validate_sample_name(command, sample.as_bytes(), "--sample-name");
    }
    if let Some(sample) = bam_sample {
        return validate_sample_name(command, sample, "BAM @RG SM");
    }
    let sample = input
        .file_stem()
        .and_then(|value| value.to_str())
        .ok_or_else(|| {
            CallError::configuration(format!(
                "{command}: cannot derive a UTF-8 sample name from BAM {}; supply --sample-name",
                input.display()
            ))
        })?;
    validate_sample_name(command, sample.as_bytes(), "BAM filename stem")
}

pub(crate) fn validate_explicit_sample_name(
    command: &str,
    explicit: Option<&str>,
) -> Result<(), CallError> {
    if let Some(sample) = explicit {
        validate_sample_name(command, sample.as_bytes(), "--sample-name")?;
    }
    Ok(())
}

fn validate_sample_name(command: &str, sample: &[u8], source: &str) -> Result<Vec<u8>, CallError> {
    let sample_text = std::str::from_utf8(sample).map_err(|_| {
        CallError::configuration(format!(
            "{command}: {source} is not valid UTF-8; supply --sample-name"
        ))
    })?;
    if sample.is_empty()
        || sample_text
            .chars()
            .any(|character| character.is_control() || character.is_whitespace())
    {
        return Err(CallError::configuration(format!(
            "{command}: {source} sample name must be nonempty and contain no whitespace or control characters"
        )));
    }
    Ok(sample.to_vec())
}

pub(crate) struct PreparedCallInput {
    pub(crate) references: Vec<BamReference>,
    pub(crate) regions: Vec<CallRegion>,
    pub(crate) worker_count: usize,
    pub(crate) reference: CallReferenceSource,
    pub(crate) bam_sample_name: Option<Vec<u8>>,
}

pub(crate) fn prepare_call_input(
    command: &str,
    path: &Path,
    reference_path: &Path,
    region_selection: &RegionSelection,
    threads: usize,
    mode: IndexedCallMode,
) -> Result<PreparedCallInput, CallError> {
    if threads == 0 {
        return Err(CallError::configuration(format!(
            "{command}: indexed calling requires at least one thread"
        )));
    }
    let reader = IndexedBamReader::open(path).map_err(|error| {
        CallError::with_source(
            CallErrorKind::Input,
            format!("{command}: open indexed BAM {}", path.display()),
            error,
        )
    })?;
    if !reader.header().has_program(b"bsbit", b"bsbit") {
        return Err(CallError::input(format!(
            "{command}: BAM {} header does not contain `@PG ID:bsbit PN:bsbit`",
            path.display()
        )));
    }
    let provenance = reader
        .header()
        .bsbit_program_provenance()
        .map_err(|error| {
            CallError::with_source(
                CallErrorKind::Input,
                format!("{command}: validate BAM {} provenance", path.display()),
                error,
            )
        })?
        .ok_or_else(|| {
            CallError::input(format!(
                "{command}: BAM {} lacks structured bsbit reference/alignment provenance",
                path.display()
            ))
        })?;
    if !provenance.alignment_mode().is_caller_compatible() {
        return Err(CallError::input(format!(
            "{command}: BAM {} was produced by the {:?} alignment path, which does not provide caller-calibrated mapping quality",
            path.display(),
            provenance.alignment_mode()
        )));
    }
    if !reader.header().is_coordinate_sorted() {
        return Err(CallError::input(format!(
            "{command}: BAM {} requires `@HD SO:coordinate`; run name sort, fixmate -m, coordinate sort, markdup, and index first",
            path.display()
        )));
    }
    let bam_sample_name = resolve_single_bam_sample(command, path, reader.header())?;
    let references = reader
        .header()
        .references()
        .iter()
        .enumerate()
        .map(|(ordinal, reference)| {
            let length = u32::try_from(reference.length()).map_err(|_| {
                CallError::input(format!(
                    "reference {ordinal} length {} exceeds the current u32 calling coordinate contract",
                    reference.length()
                ))
            })?;
            Ok(BamReference {
                name: reference.name().to_vec(),
                length,
            })
        })
        .collect::<Result<Vec<_>, CallError>>()
        .map_err(|error| error.with_context(format!("{command}: read BAM dictionary")))?;
    reader.close().map_err(|error| {
        CallError::with_source(
            CallErrorKind::Input,
            format!("{command}: close indexed BAM {}", path.display()),
            error,
        )
    })?;
    let reference = CallReferenceSource::prepare(reference_path, &references)?;
    let mut reference_reader = reference.open()?;
    reference_reader.validate_semantic_digest(
        &references,
        ReferenceSemanticDigest::from_bytes(provenance.reference_semantic_digest()),
    )?;
    reference_reader.close().map_err(|error| {
        error.with_context(format!(
            "{command}: validate reference FASTA {}",
            reference_path.display()
        ))
    })?;
    let region_bases = region_bases_for(mode, threads);
    let regions = plan_call_regions(&references, region_bases, region_selection)
        .map_err(|error| error.with_context(format!("{command}: plan calling regions")))?;
    let worker_count = threads.min(regions.len()).max(1);
    Ok(PreparedCallInput {
        references,
        regions,
        worker_count,
        reference,
        bam_sample_name,
    })
}

fn resolve_single_bam_sample(
    command: &str,
    path: &Path,
    header: &bsbit_hts::IndexedBamHeader,
) -> Result<Option<Vec<u8>>, CallError> {
    resolve_single_sample_names(command, path, header.read_group_sample_names())
}

fn resolve_single_sample_names<'a>(
    command: &str,
    path: &Path,
    sample_names: impl IntoIterator<Item = &'a [u8]>,
) -> Result<Option<Vec<u8>>, CallError> {
    let mut samples = BTreeSet::new();
    for sample in sample_names {
        if sample.is_empty() {
            return Err(CallError::input(format!(
                "{command}: BAM {} contains an empty @RG SM field",
                path.display()
            )));
        }
        samples.insert(sample.to_vec());
    }
    if samples.len() > 1 {
        let names = samples
            .iter()
            .map(|sample| String::from_utf8_lossy(sample).into_owned())
            .collect::<Vec<_>>()
            .join(", ");
        return Err(CallError::input(format!(
            "{command}: BAM {} contains multiple biological samples in @RG SM ({names}); bsbit call accepts one sample per BAM",
            path.display()
        )));
    }
    Ok(samples.into_iter().next())
}

#[cfg(test)]
mod tests {
    use std::path::Path;

    use super::{resolve_sample_name, resolve_single_sample_names};
    use crate::CallErrorKind;

    #[test]
    fn bam_read_groups_must_resolve_to_one_biological_sample() {
        let path = Path::new("sample.bam");
        assert_eq!(
            resolve_single_sample_names(
                "call snp",
                path,
                [b"donor-A".as_slice(), b"donor-A".as_slice()]
            )
            .unwrap(),
            Some(b"donor-A".to_vec())
        );
        assert_eq!(
            resolve_single_sample_names("call snp", path, std::iter::empty()).unwrap(),
            None
        );
        assert_eq!(
            resolve_sample_name("call snp", path, None, Some(b"donor-A")).unwrap(),
            b"donor-A"
        );
        assert_eq!(
            resolve_sample_name("call snp", path, Some("renamed"), Some(b"donor-A")).unwrap(),
            b"renamed"
        );

        let error = resolve_single_sample_names(
            "call snp",
            path,
            [b"donor-A".as_slice(), b"donor-B".as_slice()],
        )
        .unwrap_err();
        assert_eq!(error.kind(), CallErrorKind::Input);
        assert!(error.to_string().contains("one sample per BAM"));
    }
}
