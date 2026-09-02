//! Matrix output planning, encoding, and replaceable publication.

use std::io::{self, Write};
use std::path::{Path, PathBuf};

use bsbit_hts::{CompletedTextOutput, TextOutputCompression, TextPublication, TextStagingWriter};

use crate::input::ContigCatalog;
use crate::request::{MatrixFormat, Options};
use crate::result::{CombineError, CombineErrorKind, CombineReport, CombineWarning};
use crate::site::{Counts, SiteKey};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum MatrixKind {
    Level,
    Count,
}

impl MatrixKind {
    const fn name(self) -> &'static str {
        match self {
            Self::Level => "level",
            Self::Count => "count",
        }
    }

    const fn suffix(self, compressed: bool) -> &'static str {
        match (self, compressed) {
            (Self::Level, true) => ".level.bed.gz",
            (Self::Level, false) => ".level.bed",
            (Self::Count, true) => ".count.bed.gz",
            (Self::Count, false) => ".count.bed",
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct OutputSpec {
    pub(crate) kind: MatrixKind,
    pub(crate) path: PathBuf,
}

pub(crate) struct MatrixOutput {
    pub(crate) spec: OutputSpec,
    pub(crate) writer: TextStagingWriter,
}

pub(crate) struct CompletedMatrixOutput {
    spec: OutputSpec,
    output: CompletedTextOutput,
}

pub(crate) fn output_specs(options: &Options) -> Result<Vec<OutputSpec>, CombineError> {
    match options.matrix_format {
        MatrixFormat::Level => Ok(vec![OutputSpec {
            kind: MatrixKind::Level,
            path: matrix_path(options, MatrixKind::Level)?,
        }]),
        MatrixFormat::Count => Ok(vec![OutputSpec {
            kind: MatrixKind::Count,
            path: matrix_path(options, MatrixKind::Count)?,
        }]),
        MatrixFormat::Both => Ok(vec![
            OutputSpec {
                kind: MatrixKind::Level,
                path: matrix_path(options, MatrixKind::Level)?,
            },
            OutputSpec {
                kind: MatrixKind::Count,
                path: matrix_path(options, MatrixKind::Count)?,
            },
        ]),
    }
}

fn matrix_path(options: &Options, kind: MatrixKind) -> Result<PathBuf, CombineError> {
    let file_name = options
        .output_prefix
        .file_name()
        .filter(|name| !name.is_empty())
        .ok_or_else(|| {
            CombineError::configuration("combine: output prefix must contain a filename")
        })?;
    let mut output_name = file_name.to_os_string();
    output_name.push(kind.suffix(options.compress));
    Ok(options.output_prefix.with_file_name(output_name))
}

pub(crate) fn create_outputs(
    options: &Options,
    specs: &[OutputSpec],
) -> Result<Vec<MatrixOutput>, CombineError> {
    specs
        .iter()
        .map(|spec| {
            Ok(MatrixOutput {
                spec: spec.clone(),
                writer: create_output(options, &spec.path)?,
            })
        })
        .collect()
}

fn create_output(options: &Options, path: &Path) -> Result<TextStagingWriter, CombineError> {
    let compression = if options.compress {
        TextOutputCompression::Bgzf
    } else {
        TextOutputCompression::Plain
    };
    let compression_threads = u32::from(options.compress && options.threads > 1);
    TextStagingWriter::create_sibling_replace(path, compression, compression_threads).map_err(
        |error| {
            CombineError::with_source(
                CombineErrorKind::Output,
                format!("combine: create output staging for {}", path.display()),
                error,
            )
        },
    )
}

pub(crate) fn output_error(path: &Path, error: io::Error) -> CombineError {
    CombineError::with_source(
        CombineErrorKind::Output,
        format!("combine: write output staging for {}", path.display()),
        error,
    )
}

pub(crate) fn finish_outputs(
    outputs: Vec<MatrixOutput>,
) -> Result<Vec<CompletedMatrixOutput>, CombineError> {
    outputs
        .into_iter()
        .map(|output| {
            let completed = output.writer.finish().map_err(|error| {
                CombineError::with_source(
                    CombineErrorKind::Output,
                    format!(
                        "combine: finalize output for {}",
                        output.spec.path.display()
                    ),
                    error,
                )
            })?;
            Ok(CompletedMatrixOutput {
                spec: output.spec,
                output: completed,
            })
        })
        .collect()
}

pub(crate) fn publish_outputs(
    outputs: Vec<CompletedMatrixOutput>,
    report: &mut CombineReport,
) -> Result<(), CombineError> {
    let mut publications = Vec::<TextPublication>::with_capacity(outputs.len());
    for output in outputs {
        let path = output.spec.path;
        match output.output.publish_replace() {
            Ok(publication) => publications.push(publication),
            Err(error) => {
                for publication in publications.into_iter().rev() {
                    if let Err(rollback) = publication.rollback() {
                        return Err(CombineError::with_source(
                            CombineErrorKind::Publication,
                            format!(
                                "combine: publish output {} failed ({error}); rollback of an earlier matrix failed",
                                path.display()
                            ),
                            rollback,
                        ));
                    }
                }
                return Err(CombineError::with_source(
                    CombineErrorKind::Publication,
                    format!("combine: publish output {}", path.display()),
                    error,
                ));
            }
        }
    }
    for publication in publications {
        if let Some(warning) = publication.cleanup_warning() {
            report.warnings.push(CombineWarning {
                message: format!(
                    "combined matrix {} was published, but staging {} could not be cleaned: {warning:?}",
                    publication.target_path().display(),
                    publication.staging_path().display()
                ),
            });
        }
    }
    Ok(())
}

pub(crate) fn write_header(
    writer: &mut impl Write,
    options: &Options,
    matrix_kind: MatrixKind,
) -> io::Result<()> {
    writeln!(writer, "##bsbit_matrix_format={}", matrix_kind.name())?;
    writeln!(
        writer,
        "##bsbit_min_count={}",
        options.parameters.minimum_count
    )?;
    let proportion = options
        .parameters
        .minimum_sample_proportion_parts_per_billion;
    writeln!(
        writer,
        "##bsbit_min_prop={}.{:09}",
        proportion / 1_000_000_000,
        proportion % 1_000_000_000
    )?;
    writeln!(writer, "##bsbit_cg_only={}", options.parameters.cg_only)?;
    writer.write_all(b"#chrom\tstart\tend\tmodification\tscore\tstrand")?;
    for input in &options.inputs {
        match matrix_kind {
            MatrixKind::Level => write!(writer, "\t{}", input.sample)?,
            MatrixKind::Count => write!(
                writer,
                "\t{}_meth_count\t{}_total_count",
                input.sample, input.sample
            )?,
        }
    }
    writer.write_all(b"\n")
}

pub(crate) fn write_matrix_row(
    writer: &mut impl Write,
    matrix_kind: MatrixKind,
    options: &Options,
    contigs: &ContigCatalog,
    key: SiteKey,
    modification: &[u8],
    values: &[(usize, Counts)],
) -> io::Result<()> {
    let contig = &contigs.names[usize::try_from(key.contig).expect("u32 fits usize")];
    writer.write_all(contig)?;
    write!(writer, "\t{}\t{}\t", key.start, key.end)?;
    writer.write_all(modification)?;
    write!(writer, "\t0\t{}", if key.strand == 0 { '+' } else { '-' })?;

    let mut value_index = 0;
    for sample_index in 0..options.inputs.len() {
        let counts = values
            .get(value_index)
            .filter(|(index, _)| *index == sample_index)
            .map(|(_, counts)| *counts);
        if counts.is_some() {
            value_index += 1;
        }
        let valid = counts
            .filter(|counts| counts.total != 0 && counts.total >= options.parameters.minimum_count);
        match (matrix_kind, valid) {
            (MatrixKind::Level, Some(counts)) => {
                writer.write_all(b"\t")?;
                write_level(writer, counts)?;
            }
            (MatrixKind::Level, None) => writer.write_all(b"\t.")?,
            (MatrixKind::Count, Some(counts)) => {
                write!(writer, "\t{}\t{}", counts.methylated, counts.total)?;
            }
            (MatrixKind::Count, None) => writer.write_all(b"\t.\t.")?,
        }
    }
    debug_assert_eq!(value_index, values.len());
    writer.write_all(b"\n")
}

fn write_level(writer: &mut impl Write, counts: Counts) -> io::Result<()> {
    debug_assert!(counts.total != 0);
    debug_assert!(counts.methylated <= counts.total);
    let denominator = u128::from(counts.total);
    let scaled = (u128::from(counts.methylated) * 1_000_000 + denominator / 2) / denominator;
    let scaled = u64::try_from(scaled).expect("bounded methylation level fits u64");
    write!(writer, "{}.{:06}", scaled / 1_000_000, scaled % 1_000_000)
}
