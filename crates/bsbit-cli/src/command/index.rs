use std::path::{Path, PathBuf};

use bsbit_hts::{Compression, DecodedFastaReader, FastaRecord, TextRecordLimits};
use bsbit_index::build::combined::{
    CombinedIndexBuildOptions, build_combined_index_from_catalog_replace,
};
use bsbit_index::reference::ContigInput;
use bsbit_index::storage::reference_catalog::publish_reference_catalog_replace;
use bsbit_io::validate_replace_target;

use crate::{CliError, CliWarning, INDEX_HELP, RunReport};

use super::{
    Action, internal_search_file_prefix, option_map_with_aliases, required_path,
    unused_staging_path,
};
use super::{MAX_CLI_THREADS, optional_u64};

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct IndexOptions {
    pub(crate) reference: PathBuf,
    pub(crate) output: PathBuf,
    pub(crate) threads: u64,
}

pub(super) fn parse_index(arguments: &[String]) -> Result<Action, CliError> {
    if matches!(arguments, [value] if value == "--help" || value == "-h") {
        return Ok(Action::Help(INDEX_HELP));
    }
    let (mut values, _) = option_map_with_aliases(
        arguments,
        &["--reference", "--output", "--threads"],
        &[],
        &[
            ("-r", "--reference"),
            ("-o", "--output"),
            ("-t", "--threads"),
        ],
    )?;
    let reference = required_path(&mut values, "--reference")?;
    let output = required_path(&mut values, "--output")?;
    let threads = optional_u64(&mut values, "--threads")?.unwrap_or(1);
    if !(1..=MAX_CLI_THREADS).contains(&threads) {
        return Err(CliError::usage(format!(
            "--threads must be in 1..={MAX_CLI_THREADS}"
        )));
    }
    Ok(Action::Index(IndexOptions {
        reference,
        output,
        threads,
    }))
}

pub(crate) fn run(options: &IndexOptions) -> Result<RunReport, CliError> {
    validate_output_target(&options.output)?;
    let staging = unused_staging_path(&options.output, "index", "index")?;
    let mut reader = DecodedFastaReader::open(&options.reference, reference_text_limits())
        .map_err(|error| operation_error("open reference", &options.reference, &error))?;
    if reader.compression() == Compression::Gzip {
        let _ = reader.close();
        return Err(unsupported_gzip_reference(&options.reference));
    }
    let mut contigs = Vec::new();
    loop {
        match reader.next_record() {
            Ok(Some(record)) => {
                contigs.try_reserve(1).map_err(|_| {
                    CliError::operation(format!(
                        "index: collect reference {}: allocation failed before record {}",
                        options.reference.display(),
                        record.ordinal().get()
                    ))
                })?;
                contigs.push(contig_input_from_fasta_record(&record));
            }
            Ok(None) => break,
            Err(error) => {
                let _ = reader.close();
                return Err(operation_error(
                    "parse reference",
                    &options.reference,
                    &error,
                ));
            }
        }
    }
    reader
        .close()
        .map_err(|error| operation_error("close reference", &options.reference, &error))?;
    let publication = publish_reference_catalog_replace(&contigs, &options.output, &staging)
        .map_err(|error| operation_error("publish output", &options.output, &error))?;
    let semantic_digest = publication.summary().semantic_digest();
    let internal_prefix = internal_search_file_prefix(&options.output);
    let threads = u32::try_from(options.threads).expect("validated CLI thread count fits u32");
    let build_options = CombinedIndexBuildOptions::new(threads)
        .expect("validated CLI thread count is accepted by the index builder");
    if let Err(error) = build_combined_index_from_catalog_replace(
        contigs,
        semantic_digest,
        &internal_prefix,
        build_options,
    ) {
        let build_error =
            operation_error("build internal search data for", &options.output, &error);
        if let Err(rollback_error) = publication.rollback() {
            return Err(CliError::operation(format!(
                "{build_error}; rollback published index {}: {rollback_error}",
                options.output.display()
            )));
        }
        return Err(build_error);
    }
    let warning = publication.cleanup_error().map(|kind| {
        CliWarning::new(format!(
            "index output {} was published, but staging {} could not be removed: {kind:?}",
            options.output.display(),
            publication.staging_path().display()
        ))
    });
    Ok(RunReport::with_warning(warning))
}

fn contig_input_from_fasta_record(record: &FastaRecord) -> ContigInput {
    ContigInput::new(
        record.record_name().name().to_vec(),
        record.sequence().clone(),
    )
}

fn validate_output_target(path: &Path) -> Result<(), CliError> {
    match validate_replace_target(path) {
        Ok(()) => Ok(()),
        Err(error) if error.kind() == std::io::ErrorKind::NotADirectory => {
            let parent = output_parent(path);
            Err(CliError::operation(format!(
                "output parent {} is not a directory",
                parent.display()
            )))
        }
        Err(error) => Err(operation_error("inspect destination", path, &error)),
    }
}

fn output_parent(path: &Path) -> &Path {
    path.parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."))
}

fn operation_error(operation: &str, path: &Path, error: &impl std::fmt::Display) -> CliError {
    CliError::operation(format!("index: {operation} {}: {error}", path.display()))
}

fn unsupported_gzip_reference(path: &Path) -> CliError {
    CliError::operation(format!(
        "index: reference FASTA {} uses ordinary gzip compression, which is unsupported because it cannot provide random access; use plain FASTA or BGZF-compressed FASTA. To convert, run `gzip -cd INPUT.fa.gz | bgzip -c > REFERENCE.bgzf.fa.gz`, then `samtools faidx REFERENCE.bgzf.fa.gz` before calling",
        path.display()
    ))
}

const fn reference_text_limits() -> TextRecordLimits {
    TextRecordLimits::new(
        u64::MAX,
        1_000_000,
        64_000_000,
        1_000_000,
        u64::MAX,
        u64::MAX,
        0,
    )
}

#[cfg(test)]
mod tests {
    use std::io::{BufReader, Cursor};

    use bsbit_hts::FastaReader;

    use super::{contig_input_from_fasta_record, reference_text_limits};

    #[test]
    fn fasta_description_is_not_promoted_into_the_contig_name() {
        let input = b">chr1 retained human-readable description\nACGTN\n";
        let mut reader = FastaReader::new(
            BufReader::new(Cursor::new(input.as_slice())),
            reference_text_limits(),
        );
        let record = reader
            .next_record()
            .expect("FASTA parses")
            .expect("one FASTA record exists");
        assert_eq!(record.record_name().name(), b"chr1");
        assert_eq!(
            record.record_name().description(),
            b"retained human-readable description",
        );

        let contig = contig_input_from_fasta_record(&record);
        assert_eq!(contig.name(), b"chr1");
        assert_eq!(contig.sequence().to_ascii(), b"ACGTN");
    }
}
