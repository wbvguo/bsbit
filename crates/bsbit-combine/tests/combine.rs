//! Public integration tests for streaming methylation-matrix assembly.

use std::fs;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use bsbit_combine::{CombineErrorKind, Input, MatrixFormat, Options, Parameters, combine};
use bsbit_hts::{DecodedReader, TextOutputCompression, TextStagingWriter};

fn unique_directory(label: &str) -> PathBuf {
    let nonce = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("clock is after epoch")
        .as_nanos();
    std::env::temp_dir().join(format!(
        "bsbit-combine-{label}-{}-{nonce}",
        std::process::id()
    ))
}

fn bed_row(
    contig: &str,
    start: u64,
    modification: &str,
    strand: char,
    methylated: u64,
    unmethylated: u64,
) -> String {
    let total = methylated + unmethylated;
    let percent = if total == 0 {
        0
    } else {
        (u128::from(methylated) * 10_000 + u128::from(total) / 2) / u128::from(total)
    };
    format!(
        "{contig}\t{start}\t{}\t{modification}\t{total}\t{strand}\t{start}\t{}\t255,0,0\t{total}\t{}.{:02}\t{methylated}\t{unmethylated}\t0\t0\t0\t0\t0\n",
        start + 1,
        start + 1,
        percent / 100,
        percent % 100
    )
}

fn cgmap_row(
    contig: &str,
    nucleotide: char,
    position: u64,
    context: &str,
    dinucleotide: &str,
    methylated: u64,
    total: u64,
) -> String {
    let level = if total == 0 {
        "na".to_owned()
    } else {
        let scaled =
            (u128::from(methylated) * 1_000_000 + u128::from(total) / 2) / u128::from(total);
        format!("{}.{:06}", scaled / 1_000_000, scaled % 1_000_000)
    };
    format!(
        "{contig}\t{nucleotide}\t{position}\t{context}\t{dinucleotide}\t{level}\t{methylated}\t{total}\n"
    )
}

fn write_plain(path: &Path, rows: &[String]) {
    fs::write(path, rows.concat()).expect("plain bedMethyl writes");
}

fn write_bgzf(path: &Path, rows: &[String]) {
    let mut writer = TextStagingWriter::create_sibling(path, TextOutputCompression::Bgzf, 1)
        .expect("BGZF staging opens");
    writer
        .write_all(rows.concat().as_bytes())
        .expect("BGZF rows write");
    writer
        .finish()
        .expect("BGZF finishes")
        .publish_create_new()
        .expect("BGZF publishes");
}

fn decoded(path: &Path) -> String {
    let mut reader = DecodedReader::open(path).expect("combined output opens");
    let mut text = String::new();
    reader
        .read_to_string(&mut text)
        .expect("combined output decodes");
    reader.close().expect("combined output closes");
    text
}

fn input(sample: &str, path: &Path) -> Input {
    Input {
        sample: sample.to_owned(),
        path: path.to_path_buf(),
    }
}

fn options(inputs: Vec<Input>, output: PathBuf) -> Options {
    Options {
        inputs,
        output,
        matrix_format: MatrixFormat::Level,
        compress: false,
        threads: 1,
        parameters: Parameters::default(),
    }
}

#[test]
fn public_configuration_errors_have_stable_classification() {
    let directory = unique_directory("configuration-errors");
    fs::create_dir(&directory).expect("fixture directory");
    let first = directory.join("first.bed");
    let second = directory.join("second.bed");
    write_plain(&first, &[bed_row("chr1", 0, "m,CG,0", '+', 1, 1)]);
    write_plain(&second, &[bed_row("chr1", 1, "m,CG,0", '+', 1, 1)]);

    let empty = combine(&options(Vec::new(), directory.join("empty.bed")))
        .expect_err("empty input list fails");
    assert_eq!(empty.kind(), CombineErrorKind::Configuration);

    for threads in [0, 65] {
        let mut invalid_threads = options(
            vec![input("sample", &first)],
            directory.join(format!("threads-{threads}.bed")),
        );
        invalid_threads.threads = threads;
        let error = combine(&invalid_threads).expect_err("invalid thread count fails");
        assert_eq!(error.kind(), CombineErrorKind::Configuration);
    }

    let duplicate_sample = combine(&options(
        vec![input("sample", &first), input("sample", &second)],
        directory.join("duplicate-sample.bed"),
    ))
    .expect_err("duplicate sample label fails");
    assert_eq!(duplicate_sample.kind(), CombineErrorKind::Configuration);

    let duplicate_path = combine(&options(
        vec![input("first", &first), input("second", &first)],
        directory.join("duplicate-path.bed"),
    ))
    .expect_err("duplicate input path fails");
    assert_eq!(duplicate_path.kind(), CombineErrorKind::Configuration);

    let same_path = combine(&options(vec![input("sample", &first)], first.clone()))
        .expect_err("input and output path collision fails");
    assert_eq!(same_path.kind(), CombineErrorKind::Configuration);
    assert_eq!(decoded(&first), bed_row("chr1", 0, "m,CG,0", '+', 1, 1));

    fs::remove_dir_all(directory).expect("fixture cleanup");
}

#[test]
fn invalid_methylation_rows_have_stable_input_classification() {
    let directory = unique_directory("input-errors");
    fs::create_dir(&directory).expect("fixture directory");
    let malformed = directory.join("malformed.bed");
    let duplicate = directory.join("duplicate.bed");
    let non_contiguous = directory.join("non-contiguous.bed");
    write_plain(&malformed, &["chr1\t0\t1\n".to_owned()]);
    let repeated = bed_row("chr1", 0, "m,CG,0", '+', 1, 1);
    write_plain(&duplicate, &[repeated.clone(), repeated]);
    write_plain(
        &non_contiguous,
        &[
            bed_row("chr1", 0, "m,CG,0", '+', 1, 1),
            bed_row("chr2", 0, "m,CG,0", '+', 1, 1),
            bed_row("chr1", 1, "m,CG,0", '+', 1, 1),
        ],
    );

    for (label, path) in [
        ("malformed", malformed),
        ("duplicate", duplicate),
        ("non-contiguous", non_contiguous),
    ] {
        let output = directory.join(format!("{label}-matrix.bed"));
        let error = combine(&options(vec![input("sample", &path)], output.clone()))
            .expect_err("invalid methylation input fails");
        assert_eq!(error.kind(), CombineErrorKind::Input, "case {label}");
        assert!(!output.exists(), "case {label} must not publish output");
    }

    fs::remove_dir_all(directory).expect("fixture cleanup");
}

#[test]
fn cgmap_and_bed_methyl_inputs_share_one_matrix_coordinate_model() {
    let directory = unique_directory("cgmap-bed-parity");
    fs::create_dir(&directory).expect("fixture directory");
    let bed = directory.join("first.bed.gz");
    let cgmap = directory.join("second.cgmap");
    write_bgzf(
        &bed,
        &[
            bed_row("chr1", 0, "m,CG,0", '+', 3, 2),
            bed_row("chr1", 1, "m,CG,0", '-', 2, 3),
            bed_row("chr1", 3, "m,CHH,0", '+', 1, 1),
        ],
    );
    write_plain(
        &cgmap,
        &[
            cgmap_row("chr1", 'C', 1, "CG", "CG", 4, 5),
            cgmap_row("chr1", 'G', 2, "CG", "CG", 1, 5),
            cgmap_row("chr1", 'C', 4, "CHH", "CA", 0, 2),
        ],
    );
    let output = directory.join("matrix.bed");
    let report = combine(&Options {
        inputs: vec![input("bed", &bed), input("cgmap", &cgmap)],
        output: output.clone(),
        matrix_format: MatrixFormat::Count,
        compress: false,
        threads: 2,
        parameters: Parameters::default(),
    })
    .expect("CGmap and bedMethyl combine together");

    assert_eq!(report.sites_seen(), 3);
    assert_eq!(report.sites_written(), 3);
    assert_eq!(
        decoded(&output),
        concat!(
            "##bsbit_matrix_format=count\n",
            "##bsbit_min_count=1\n",
            "##bsbit_min_prop=0.000000000\n",
            "#chrom\tstart\tend\tmodification\tscore\tstrand",
            "\tbed_meth_count\tbed_total_count",
            "\tcgmap_meth_count\tcgmap_total_count\n",
            "chr1\t0\t1\tm,CG,0\t0\t+\t3\t5\t4\t5\n",
            "chr1\t1\t2\tm,CG,0\t0\t-\t2\t5\t1\t5\n",
            "chr1\t3\t4\tm,CHH,0\t0\t+\t1\t2\t0\t2\n",
        )
    );
    fs::remove_dir_all(directory).expect("fixture cleanup");
}

#[test]
fn malformed_cgmap_rows_fail_before_publication() {
    let directory = unique_directory("cgmap-errors");
    fs::create_dir(&directory).expect("fixture directory");
    let invalid_rows = [
        ("zero-position", "chr1\tC\t0\tCG\tCG\t0.5\t1\t2\n"),
        ("wrong-context", "chr1\tC\t1\tCG\tCA\t0.5\t1\t2\n"),
        ("count-overflow", "chr1\tC\t1\tCG\tCG\t1.0\t3\t2\n"),
    ];
    for (label, row) in invalid_rows {
        let input_path = directory.join(format!("{label}.cgmap"));
        fs::write(&input_path, row).expect("invalid CGmap fixture writes");
        let output = directory.join(format!("{label}.bed"));
        let error = combine(&options(vec![input("sample", &input_path)], output.clone()))
            .expect_err("invalid CGmap fails");
        assert_eq!(error.kind(), CombineErrorKind::Input);
        assert!(!output.exists());
    }
    fs::remove_dir_all(directory).expect("fixture cleanup");
}

#[test]
fn union_filter_preserves_missing_cells_and_raw_counts() {
    let directory = unique_directory("union-filter");
    fs::create_dir(&directory).expect("fixture directory");
    let first = directory.join("first.bed");
    let second = directory.join("second.bed");
    let third = directory.join("third.bed");
    write_plain(
        &first,
        &[
            bed_row("chr1", 0, "m,CG,0", '+', 8, 2),
            bed_row("chr1", 2, "m,CHH,0", '+', 1, 1),
            bed_row("chr3", 0, "m,CG,0", '-', 4, 1),
        ],
    );
    write_plain(
        &second,
        &[
            bed_row("chr1", 0, "m,CG,0", '+', 2, 3),
            bed_row("chr1", 1, "m,CG,0", '+', 0, 6),
            bed_row("chr2", 0, "m,CG,0", '+', 3, 2),
        ],
    );
    write_plain(
        &third,
        &[
            bed_row("chr1", 0, "m,CG,0", '+', 9, 1),
            bed_row("chr1", 1, "m,CG,0", '+', 1, 1),
            bed_row("chr3", 0, "m,CG,0", '-', 1, 4),
        ],
    );
    let output = directory.join("matrix.bed");
    let report = combine(&Options {
        inputs: vec![
            input("tumor", &first),
            input("normal", &second),
            input("control", &third),
        ],
        output: output.clone(),
        matrix_format: MatrixFormat::Both,
        compress: false,
        threads: 3,
        parameters: Parameters {
            minimum_count: 5,
            minimum_sample_proportion_parts_per_billion: 666_666_666,
        },
    })
    .expect("matrix combines");
    assert_eq!(report.sites_seen(), 5);
    assert_eq!(report.sites_written(), 2);
    assert!(report.warnings().is_empty());
    assert!(!output.exists());
    let level_output = directory.join("matrix.level.bed");
    let count_output = directory.join("matrix.count.bed");
    assert_eq!(
        decoded(&level_output),
        concat!(
            "##bsbit_matrix_format=level\n",
            "##bsbit_min_count=5\n",
            "##bsbit_min_prop=0.666666666\n",
            "#chrom\tstart\tend\tmodification\tscore\tstrand",
            "\ttumor\tnormal\tcontrol\n",
            "chr1\t0\t1\tm,CG,0\t0\t+\t0.800000\t0.400000\t0.900000\n",
            "chr3\t0\t1\tm,CG,0\t0\t-\t0.800000\t.\t0.200000\n",
        )
    );
    assert_eq!(
        decoded(&count_output),
        concat!(
            "##bsbit_matrix_format=count\n",
            "##bsbit_min_count=5\n",
            "##bsbit_min_prop=0.666666666\n",
            "#chrom\tstart\tend\tmodification\tscore\tstrand",
            "\ttumor_meth_count\ttumor_total_count",
            "\tnormal_meth_count\tnormal_total_count",
            "\tcontrol_meth_count\tcontrol_total_count\n",
            "chr1\t0\t1\tm,CG,0\t0\t+\t8\t10\t2\t5\t9\t10\n",
            "chr3\t0\t1\tm,CG,0\t0\t-\t4\t5\t.\t.\t1\t5\n",
        )
    );
    fs::remove_dir_all(directory).expect("fixture cleanup");
}

#[test]
fn bgzf_and_thread_counts_are_byte_deterministic() {
    let directory = unique_directory("bgzf-parity");
    fs::create_dir(&directory).expect("fixture directory");
    let first = directory.join("first.bed.gz");
    let second = directory.join("second.bed");
    let rows = [
        bed_row("chr1", 3, "m,CG,0", '+', 1, 2),
        bed_row("chr1", 9, "m,CHG,0", '+', 2, 1),
    ];
    write_bgzf(&first, &rows);
    write_plain(&second, &rows);

    let one = directory.join("one.bed.gz");
    let many = directory.join("many.bed.gz");
    for (output, threads) in [(&one, 1), (&many, 4)] {
        let report = combine(&Options {
            inputs: vec![input("a", &first), input("b", &second)],
            output: output.clone(),
            matrix_format: MatrixFormat::Level,
            compress: true,
            threads,
            parameters: Parameters::default(),
        })
        .expect("compressed matrix combines");
        assert_eq!(report.sites_seen(), 2);
        assert_eq!(report.sites_written(), 2);
    }
    assert_eq!(
        fs::read(&one).expect("one-thread bytes"),
        fs::read(&many).expect("many-thread bytes")
    );
    assert!(decoded(&one).contains("chr1\t3\t4\tm,CG,0\t0\t+\t0.333333\t0.333333\n"));

    let both = directory.join("both.bed.gz");
    combine(&Options {
        inputs: vec![input("a", &first), input("b", &second)],
        output: both.clone(),
        matrix_format: MatrixFormat::Both,
        compress: true,
        threads: 4,
        parameters: Parameters::default(),
    })
    .expect("both compressed matrices combine");
    assert!(!both.exists());
    assert!(
        decoded(&directory.join("both.level.bed.gz"))
            .contains("chr1\t3\t4\tm,CG,0\t0\t+\t0.333333\t0.333333\n")
    );
    assert!(
        decoded(&directory.join("both.count.bed.gz"))
            .contains("\ta_meth_count\ta_total_count\tb_meth_count\tb_total_count\n")
    );
    fs::remove_dir_all(directory).expect("fixture cleanup");
}

#[test]
fn incompatible_contig_order_fails_before_publication() {
    let directory = unique_directory("contig-cycle");
    fs::create_dir(&directory).expect("fixture directory");
    let first = directory.join("first.bed");
    let second = directory.join("second.bed");
    write_plain(
        &first,
        &[
            bed_row("chr1", 0, "m,CG,0", '+', 1, 1),
            bed_row("chr2", 0, "m,CG,0", '+', 1, 1),
        ],
    );
    write_plain(
        &second,
        &[
            bed_row("chr2", 0, "m,CG,0", '+', 1, 1),
            bed_row("chr1", 0, "m,CG,0", '+', 1, 1),
        ],
    );
    let output = directory.join("matrix.bed");
    let error = combine(&Options {
        inputs: vec![input("a", &first), input("b", &second)],
        output: output.clone(),
        matrix_format: MatrixFormat::Count,
        compress: false,
        threads: 2,
        parameters: Parameters::default(),
    })
    .expect_err("contig-order cycle fails");
    assert_eq!(error.kind(), CombineErrorKind::Input);
    assert!(error.to_string().contains("contig order"));
    assert!(!output.exists());
    fs::remove_dir_all(directory).expect("fixture cleanup");
}

#[test]
fn metadata_mismatch_fails_closed_and_existing_targets_are_replaced() {
    let directory = unique_directory("fail-closed");
    fs::create_dir(&directory).expect("fixture directory");
    let first = directory.join("first.bed");
    let second = directory.join("second.bed");
    write_plain(&first, &[bed_row("chr1", 0, "m,CG,0", '+', 1, 1)]);
    write_plain(&second, &[bed_row("chr1", 0, "m,CHH,0", '+', 1, 1)]);
    let mismatch_output = directory.join("mismatch.bed");
    let mismatch = combine(&Options {
        inputs: vec![input("a", &first), input("b", &second)],
        output: mismatch_output.clone(),
        matrix_format: MatrixFormat::Level,
        compress: false,
        threads: 2,
        parameters: Parameters::default(),
    })
    .expect_err("context mismatch fails");
    assert_eq!(mismatch.kind(), CombineErrorKind::Input);
    assert!(mismatch.to_string().contains("modification/context"));
    assert!(!mismatch_output.exists());

    let existing = directory.join("existing.bed");
    fs::write(&existing, b"owned\n").expect("existing target");
    combine(&Options {
        inputs: vec![input("a", &first)],
        output: existing.clone(),
        matrix_format: MatrixFormat::Level,
        compress: false,
        threads: 1,
        parameters: Parameters::default(),
    })
    .expect("existing target is replaced");
    assert!(decoded(&existing).contains("#chrom"));
    assert_ne!(fs::read(&existing).expect("replacement bytes"), b"owned\n");

    let both_template = directory.join("cohort.bed.gz");
    let existing_count = directory.join("cohort.count.bed.gz");
    let absent_level = directory.join("cohort.level.bed.gz");
    fs::write(&existing_count, b"owned count\n").expect("existing count target");
    combine(&Options {
        inputs: vec![input("a", &first)],
        output: both_template.clone(),
        matrix_format: MatrixFormat::Both,
        compress: true,
        threads: 1,
        parameters: Parameters::default(),
    })
    .expect("both outputs replace existing destinations");
    assert!(decoded(&existing_count).contains("##bsbit_matrix_format=count"));
    assert!(decoded(&absent_level).contains("##bsbit_matrix_format=level"));
    assert!(!both_template.exists());
    fs::remove_dir_all(directory).expect("fixture cleanup");
}
