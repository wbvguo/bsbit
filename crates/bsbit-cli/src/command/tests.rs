use super::index::IndexSpeed;
use super::{Action, parse};
use bsbit_call::meth::OutputFormat as MethylationOutputFormat;
use bsbit_combine::MatrixFormat as CombineMatrixFormat;
use bsbit_cpu::BackendRequest;
use bsbit_index::build::combined::DEFAULT_COMBINED_INDEX_MEMORY_MIB;

fn arguments(values: &[&str]) -> Vec<std::ffi::OsString> {
    values.iter().map(std::ffi::OsString::from).collect()
}

#[cfg(unix)]
#[test]
fn filesystem_paths_preserve_non_utf8_bytes_across_every_command_parser() {
    use std::ffi::OsString;
    use std::os::unix::ffi::{OsStrExt, OsStringExt};

    let raw = |stem: &[u8]| {
        let mut bytes = stem.to_vec();
        bytes.push(0xff);
        OsString::from_vec(bytes)
    };
    let index_path = raw(b"reference-");
    let reference_path = raw(b"source-");
    let output_path = raw(b"result-");
    let read_path = raw(b"reads-");

    let Action::Index(index) = parse(vec![
        OsString::from("index"),
        OsString::from("--reference"),
        reference_path.clone(),
        OsString::from("--output"),
        index_path.clone(),
    ])
    .expect("non-UTF-8 index paths parse") else {
        panic!("expected index action");
    };
    assert_eq!(
        index.reference.as_os_str().as_bytes(),
        reference_path.as_bytes()
    );
    assert_eq!(index.output.as_os_str().as_bytes(), index_path.as_bytes());

    let Action::Align(_) = parse(vec![
        OsString::from("align"),
        OsString::from("--index"),
        index_path.clone(),
        OsString::from("--read1"),
        read_path.clone(),
        OsString::from("--output"),
        output_path.clone(),
    ])
    .expect("non-UTF-8 alignment paths parse") else {
        panic!("expected align action");
    };
    let Action::CallMeth(call) = parse(vec![
        OsString::from("call"),
        OsString::from("meth"),
        OsString::from("--input"),
        read_path.clone(),
        OsString::from("--reference"),
        reference_path.clone(),
        OsString::from("--output"),
        output_path.clone(),
        OsString::from("--format"),
        OsString::from("cgmap"),
    ])
    .expect("non-UTF-8 calling paths parse") else {
        panic!("expected methylation-call action");
    };
    assert_eq!(call.input.as_os_str().as_bytes(), read_path.as_bytes());
    assert_eq!(
        call.reference.as_os_str().as_bytes(),
        reference_path.as_bytes()
    );
    assert_eq!(call.output.as_os_str().as_bytes(), output_path.as_bytes());

    let Action::Combine(combine) = parse(vec![
        OsString::from("combine"),
        OsString::from("--input"),
        read_path.clone(),
        OsString::from("--sample-name"),
        OsString::from("sample"),
        OsString::from("--output"),
        output_path.clone(),
    ])
    .expect("non-UTF-8 combine paths parse with an explicit sample name") else {
        panic!("expected combine action");
    };
    assert_eq!(
        combine.inputs[0].path.as_os_str().as_bytes(),
        read_path.as_bytes()
    );
    assert_eq!(
        combine.output.as_os_str().as_bytes(),
        output_path.as_bytes()
    );
}

#[test]
fn exact_index_options_parse() {
    let Action::Index(index) = parse(arguments(&[
        "index",
        "--reference",
        "ref.fa",
        "--output",
        "ref.bsbit",
        "--threads",
        "4",
    ]))
    .expect("index parses") else {
        panic!("expected index action");
    };
    assert_eq!(index.reference, std::path::Path::new("ref.fa"));
    assert_eq!(index.output, std::path::Path::new("ref.bsbit"));
    assert_eq!(index.threads, 4);
    assert_eq!(index.memory_mib, DEFAULT_COMBINED_INDEX_MEMORY_MIB);
    assert_eq!(index.speed, IndexSpeed::Fast);
    assert_eq!(index.simd_backend, BackendRequest::Auto);
}

#[test]
fn short_index_options_and_fast_layout_parse() {
    let Action::Index(index) = parse(arguments(&[
        "index",
        "-r",
        "ref.fa",
        "-o",
        "ref.bsbit",
        "-t",
        "4",
        "--index-speed",
        "fast",
        "--memory-mib",
        "512",
        "--simd-backend",
        "scalar",
    ]))
    .expect("short index options parse") else {
        panic!("expected index action");
    };
    assert_eq!(index.reference, std::path::Path::new("ref.fa"));
    assert_eq!(index.output, std::path::Path::new("ref.bsbit"));
    assert_eq!(index.threads, 4);
    assert_eq!(index.memory_mib, 512);
    assert_eq!(index.speed, IndexSpeed::Fast);
    assert_eq!(index.simd_backend, BackendRequest::Scalar);

    let error = parse(arguments(&[
        "index",
        "-r",
        "first.fa",
        "--reference",
        "second.fa",
        "-o",
        "ref.bsbit",
    ]))
    .expect_err("short and long reference forms are one option");
    assert!(error.to_string().contains("duplicate option `--reference`"));
}

#[test]
fn compact_index_layout_is_explicit_and_retired_balanced_spelling_is_rejected() {
    let Action::Index(index) = parse(arguments(&[
        "index",
        "-r",
        "ref.fa",
        "-o",
        "ref.bsbit",
        "--index-speed",
        "compact",
    ]))
    .expect("compact index layout parses") else {
        panic!("expected index action");
    };
    assert_eq!(index.speed, IndexSpeed::Compact);

    let error = parse(arguments(&[
        "index",
        "-r",
        "ref.fa",
        "-o",
        "ref.bsbit",
        "--index-speed",
        "balanced",
    ]))
    .expect_err("retired balanced spelling must fail closed");
    assert!(error.to_string().contains("expected fast or compact"));
}

#[test]
fn internal_index_construction_is_not_a_public_subcommand() {
    assert!(parse(arguments(&["index", "combined", "--snapshot", "ref.bsbit"])).is_err());
}

#[test]
fn index_resource_and_backend_controls_fail_closed() {
    for extra in [
        ["--memory-mib", "0"],
        ["--threads", "0"],
        ["--threads", "2147483648"],
        ["--simd-backend", "native"],
    ] {
        let supplied = [
            vec!["index", "-r", "ref.fa", "-o", "ref.bsbit"],
            extra.to_vec(),
        ]
        .concat();
        assert!(parse(arguments(&supplied)).is_err(), "accepted {extra:?}");
    }
}

#[test]
fn public_help_lists_supported_index_and_alignment_entry_points() {
    let help = crate::GENERAL_HELP;
    assert!(help.contains("bsbit index"));
    assert!(help.contains("bsbit align -x PATH -1 PATH"));
    assert!(crate::ALIGN_HELP.contains("-x, --index PATH"));
    assert!(crate::ALIGN_HELP.contains("-1, --read1 PATH"));
    assert!(crate::ALIGN_HELP.contains("-2, --read2 PATH"));
}

#[test]
fn standard_alignment_layout_and_official_short_forms_parse() {
    let Action::Align(_) = parse(arguments(&[
        "align",
        "--index",
        "ref.bsbit",
        "-1",
        "reads.fq",
        "--output",
        "out.bam",
    ]))
    .expect("single parses") else {
        panic!("expected align action");
    };

    let Action::Align(_) = parse(arguments(&[
        "align",
        "--index",
        "ref.bsbit",
        "--read1",
        "r1.fq",
        "--read2",
        "r2.fq",
        "--output",
        "out.bam",
    ]))
    .expect("paired input parses") else {
        panic!("expected align action");
    };
    let duplicate = parse(arguments(&[
        "align",
        "--index",
        "ref.bsbit",
        "--read1",
        "first.fq",
        "-1",
        "second.fq",
        "--output",
        "out.bam",
    ]))
    .expect_err("short and long read-1 forms are one option");
    assert!(
        duplicate
            .to_string()
            .contains("--read1 may be specified only once")
    );
}

#[test]
fn nested_methylation_call_accepts_short_and_long_options() {
    let Action::CallMeth(short) = parse(arguments(&[
        "call",
        "meth",
        "-i",
        "reads.bam",
        "-r",
        "reference.fa",
        "-o",
        "calls.cgmap.gz",
        "-f",
        "cgmap",
        "-c",
        "true",
        "--region",
        "chr1:1-10",
        "--region",
        "chr1:21-1,000",
        "--regions-bed",
        "targets.bed.gz",
    ]))
    .expect("short call options parse") else {
        panic!("expected methylation-call action");
    };
    assert_eq!(short.format, MethylationOutputFormat::Cgmap);
    assert!(short.compress);
    assert_eq!(short.threads, 1);
    assert_eq!(short.compression_threads, 0);
    assert_eq!(short.parameters.minimum_base_quality, 20);
    assert_eq!(short.parameters.minimum_depth, 10);
    assert_eq!(short.parameters.minimum_mapping_quality, 20);
    assert_eq!(short.regions.intervals.len(), 2);
    assert_eq!(short.regions.intervals[0].start, 0);
    assert_eq!(short.regions.intervals[0].end, 10);
    assert_eq!(short.regions.intervals[1].start, 20);
    assert_eq!(short.regions.intervals[1].end, 1_000);
    assert_eq!(
        short.regions.regions_file.as_deref(),
        Some(std::path::Path::new("targets.bed.gz"))
    );

    let Action::CallMeth(long) = parse(arguments(&[
        "call",
        "meth",
        "--input",
        "reads.bam",
        "--reference",
        "reference.fa",
        "--output",
        "calls.bed",
        "--format",
        "bed",
        "--compress",
        "true",
        "--threads",
        "4",
        "--compression-threads",
        "3",
        "--min-bq",
        "25",
        "--min-mapq",
        "30",
        "--min-depth",
        "12",
        "--cg-only",
        "--ignore-orphan",
    ]))
    .expect("long call options parse") else {
        panic!("expected methylation-call action");
    };
    assert_eq!(long.format, MethylationOutputFormat::Bed);
    assert_eq!(long.reference, std::path::Path::new("reference.fa"));
    assert!(long.compress);
    assert_eq!(long.threads, 4);
    assert_eq!(long.compression_threads, 3);
    assert_eq!(long.parameters.minimum_base_quality, 25);
    assert_eq!(long.parameters.minimum_mapping_quality, 30);
    assert_eq!(long.parameters.minimum_depth, 12);
    assert!(long.parameters.cg_only);
    assert!(long.parameters.ignore_orphans);

    let duplicate_reference = parse(arguments(&[
        "call",
        "meth",
        "-i",
        "reads.bam",
        "-r",
        "first.fa",
        "--reference",
        "second.fa",
        "-o",
        "calls.cgmap",
        "-f",
        "cgmap",
    ]))
    .expect_err("short and long reference forms are one option");
    assert!(
        duplicate_reference
            .to_string()
            .contains("duplicate option `--reference`")
    );
}

#[test]
fn nested_methylation_call_rejects_invalid_options() {
    for invalid in [
        arguments(&["call"]),
        arguments(&["call", "unknown"]),
        arguments(&["call", "meth", "-i", "reads.bam"]),
        arguments(&[
            "call",
            "meth",
            "-i",
            "reads.bam",
            "-o",
            "calls",
            "-f",
            "wig",
        ]),
        arguments(&[
            "call",
            "meth",
            "-i",
            "reads.bam",
            "--input",
            "other.bam",
            "-o",
            "calls",
            "-f",
            "cgmap",
        ]),
        arguments(&[
            "call",
            "meth",
            "-i",
            "reads.bam",
            "-o",
            "calls",
            "-f",
            "cgmap",
            "-c",
            "yes",
        ]),
        arguments(&[
            "call",
            "meth",
            "-i",
            "reads.bam",
            "-o",
            "calls",
            "-f",
            "cgmap",
            "-t",
            "0",
        ]),
        arguments(&[
            "call",
            "meth",
            "-i",
            "reads.bam",
            "-o",
            "calls",
            "-f",
            "cgmap",
            "--min-bq",
            "94",
        ]),
        arguments(&[
            "call",
            "meth",
            "-i",
            "reads.bam",
            "-o",
            "calls",
            "-f",
            "cgmap",
            "--min-mapq",
            "255",
        ]),
    ] {
        assert!(parse(invalid).is_err());
    }
}

#[test]
fn call_region_coordinates_fail_closed() {
    for region in ["chr1:0-10", "chr1:20-10", "chr1:1,00-200", "chr1"] {
        assert!(
            parse(arguments(&[
                "call",
                "meth",
                "-i",
                "reads.bam",
                "--reference",
                "reference.fa",
                "-o",
                "calls",
                "-f",
                "cgmap",
                "--region",
                region,
            ]))
            .is_err()
        );
    }
}

#[test]
fn all_call_modules_require_a_reference_argument() {
    assert!(
        parse(arguments(&[
            "call",
            "snp",
            "-i",
            "reads.bam",
            "-o",
            "calls.vcf",
        ]))
        .is_err()
    );
    assert!(
        parse(arguments(&[
            "call",
            "meth",
            "-i",
            "reads.bam",
            "-o",
            "calls.cgmap",
            "-f",
            "cgmap",
        ]))
        .is_err()
    );
    assert!(
        parse(arguments(&[
            "call",
            "joint",
            "-i",
            "reads.bam",
            "-p",
            "calls",
            "-f",
            "cgmap",
        ]))
        .is_err()
    );
}

#[test]
fn methylation_combine_parses_comma_separated_inputs_names_and_filters() {
    let Action::Combine(options) = parse(arguments(&[
        "combine",
        "-i",
        "tumor.bed.gz,normal.bed",
        "--sample-name",
        "tumor,normal",
        "-o",
        "matrix",
        "--matrix",
        "both",
        "--min-count",
        "10",
        "--min-prop",
        "0.8",
        "--cg-only",
        "-c",
        "true",
        "-t",
        "8",
        "--compression-threads",
        "3",
    ]))
    .expect("combine options parse") else {
        panic!("expected combine action");
    };
    assert_eq!(options.inputs.len(), 2);
    assert_eq!(options.inputs[0].sample, "tumor");
    assert_eq!(options.inputs[1].sample, "normal");
    assert_eq!(options.matrix_format, CombineMatrixFormat::Both);
    assert_eq!(options.output, std::path::Path::new("matrix"));
    assert_eq!(options.parameters.minimum_count, 10);
    assert_eq!(
        options
            .parameters
            .minimum_sample_proportion_parts_per_billion,
        800_000_000
    );
    assert!(options.compress);
    assert!(options.parameters.cg_only);
    assert_eq!(options.threads, 8);
    assert_eq!(options.compression_threads, 3);
}

#[test]
fn methylation_combine_defaults_names_to_exact_paths() {
    let Action::Combine(defaulted) = parse(arguments(&[
        "combine",
        "--input",
        "cohort/tumor.bed.gz,cohort/normal sample.bed",
        "--output",
        "matrix",
    ]))
    .expect("path-derived sample names parse") else {
        panic!("expected combine action");
    };
    assert_eq!(defaulted.inputs.len(), 2);
    assert_eq!(defaulted.output, std::path::Path::new("matrix"));
    assert!(!defaulted.compress);
    assert!(!defaulted.parameters.cg_only);
    assert_eq!(defaulted.inputs[0].sample, "cohort/tumor.bed.gz");
    assert_eq!(defaulted.inputs[1].sample, "cohort/normal sample.bed");
    assert_eq!(
        defaulted.inputs[1].path,
        std::path::Path::new("cohort/normal sample.bed")
    );

    let Action::Combine(equals_path) = parse(arguments(&[
        "combine",
        "--input",
        "tumor=one.bed",
        "--output",
        "matrix.bed",
    ]))
    .expect("equals sign remains part of the path") else {
        panic!("expected combine action");
    };
    assert_eq!(equals_path.inputs[0].sample, "tumor=one.bed");
    assert_eq!(
        equals_path.inputs[0].path,
        std::path::Path::new("tumor=one.bed")
    );
}

#[test]
fn methylation_combine_rejects_ambiguous_or_invalid_options() {
    for invalid in [
        arguments(&["combine", "-o", "matrix.bed"]),
        arguments(&[
            "combine",
            "-i",
            "one.bed,two.bed",
            "--sample-name",
            "only-one",
            "-o",
            "matrix.bed",
        ]),
        arguments(&["combine", "-i", "one.bed,,two.bed", "-o", "matrix.bed"]),
        arguments(&[
            "combine",
            "-i",
            "one.bed,two.bed",
            "--sample-name",
            "sample,,control",
            "-o",
            "matrix.bed",
        ]),
        arguments(&[
            "combine",
            "-i",
            "one.bed",
            "-i",
            "one.bed",
            "-o",
            "matrix.bed",
        ]),
        arguments(&[
            "combine",
            "-i",
            "sample=one.bed",
            "-o",
            "matrix.bed",
            "--matrix",
            "wide",
        ]),
        arguments(&[
            "combine",
            "-i",
            "sample=one.bed",
            "-o",
            "matrix.bed",
            "--min-prop",
            "1.1",
        ]),
        arguments(&[
            "combine",
            "-i",
            "sample=one.bed",
            "-o",
            "matrix.bed",
            "--threads",
            "0",
        ]),
        arguments(&[
            "combine",
            "-i",
            "sample=one.bed",
            "--prefix",
            "matrix",
            "--cg-only",
            "--cg-only",
        ]),
        arguments(&[
            "combine",
            "-i",
            "sample=one.bed",
            "--output",
            "matrix.bed",
            "--output",
            "matrix",
        ]),
    ] {
        assert!(parse(invalid).is_err());
    }
}

#[test]
fn snp_and_joint_call_modules_parse_exact_quality_controls() {
    let Action::CallSnp(snp) = parse(arguments(&[
        "call",
        "snp",
        "-i",
        "reads.bam",
        "-r",
        "reference.fa",
        "-o",
        "calls.vcf.gz",
        "--sample-name",
        "tumor-A",
        "-c",
        "true",
        "-t",
        "8",
        "--compression-threads",
        "2",
        "--min-bq",
        "25",
        "--min-mapq",
        "30",
        "--min-depth",
        "6",
        "--min-alt-count",
        "3",
        "--min-alt-fraction",
        "0.05",
        "--min-gq",
        "40",
        "--min-aq",
        "30",
        "--heterozygosity",
        "0.0005",
        "--underconversion-rate",
        "0.0025",
        "--overconversion-rate",
        "0.000001",
        "--ignore-orphan",
    ]))
    .expect("SNP call parses") else {
        panic!("expected SNP-call action");
    };
    assert!(snp.compress);
    assert_eq!(snp.reference, std::path::Path::new("reference.fa"));
    assert_eq!(snp.sample_name.as_deref(), Some("tumor-A"));
    assert_eq!(snp.threads, 8);
    assert_eq!(snp.compression_threads, 2);
    assert_eq!(snp.parameters.minimum_base_quality, 25);
    assert_eq!(snp.parameters.minimum_mapping_quality, 30);
    assert_eq!(snp.parameters.minimum_depth, 6);
    assert_eq!(snp.parameters.minimum_alternate_count, 3);
    assert_eq!(
        snp.parameters.minimum_alternate_fraction_parts_per_billion,
        50_000_000
    );
    assert_eq!(snp.parameters.minimum_genotype_quality, 40);
    assert_eq!(snp.parameters.minimum_allele_quality, 30);
    assert_eq!(snp.parameters.heterozygosity_parts_per_billion, 500_000);
    assert_eq!(snp.parameters.underconversion_parts_per_billion, 2_500_000);
    assert_eq!(snp.parameters.overconversion_parts_per_billion, 1_000);
    assert!(snp.parameters.ignore_orphans);

    let Action::CallJoint(joint) = parse(arguments(&[
        "call",
        "joint",
        "-i",
        "reads.bam",
        "-r",
        "reference.fa",
        "--sample-name",
        "tumor-A",
        "-p",
        "calls",
        "-f",
        "cgmap",
        "--heterozygosity",
        "0.002",
        "--cg-only",
        "--ignore-orphan",
    ]))
    .expect("joint call parses") else {
        panic!("expected joint-call action");
    };
    assert_eq!(joint.meth_format, MethylationOutputFormat::Cgmap);
    assert_eq!(joint.meth_output, std::path::Path::new("calls.CGmap"));
    assert_eq!(joint.vcf_output, std::path::Path::new("calls.vcf"));
    assert_eq!(joint.reference, std::path::Path::new("reference.fa"));
    assert_eq!(joint.sample_name.as_deref(), Some("tumor-A"));
    assert_eq!(joint.parameters.minimum_base_quality, 20);
    assert_eq!(joint.parameters.minimum_depth, 10);
    assert_eq!(joint.parameters.minimum_mapping_quality, 20);
    assert_eq!(joint.parameters.heterozygosity_parts_per_billion, 2_000_000);
    assert!(joint.cg_only);
    assert!(joint.parameters.ignore_orphans);
    assert!(!joint.compress);
    assert_eq!(joint.compression_threads, 0);
}

#[test]
fn call_and_combine_compression_workers_are_explicit_and_bounded() {
    for invalid in [
        arguments(&[
            "call",
            "meth",
            "--input",
            "reads.bam",
            "--reference",
            "reference.fa",
            "--output",
            "calls.cgmap",
            "--format",
            "cgmap",
            "--compress",
            "false",
            "--compression-threads",
            "1",
        ]),
        arguments(&[
            "call",
            "snp",
            "--input",
            "reads.bam",
            "--reference",
            "reference.fa",
            "--output",
            "calls.vcf.gz",
            "--compression-threads",
            "2147483648",
        ]),
        arguments(&[
            "combine",
            "--input",
            "calls.cgmap",
            "--output",
            "matrix.bed",
            "--compress",
            "false",
            "--compression-threads",
            "1",
        ]),
    ] {
        assert!(parse(invalid).is_err());
    }
}

#[test]
fn call_and_combine_default_to_plain_text() {
    let Action::CallMeth(meth) = parse(arguments(&[
        "call",
        "meth",
        "--input",
        "reads.bam",
        "--reference",
        "reference.fa",
        "--output",
        "calls.cgmap",
        "--format",
        "cgmap",
    ]))
    .expect("methylation call defaults parse") else {
        panic!("expected methylation-call action");
    };
    assert!(!meth.compress);
    assert_eq!(meth.compression_threads, 0);

    let Action::CallSnp(snp) = parse(arguments(&[
        "call",
        "snp",
        "--input",
        "reads.bam",
        "--reference",
        "reference.fa",
        "--output",
        "calls.vcf",
    ]))
    .expect("SNP call defaults parse") else {
        panic!("expected SNP-call action");
    };
    assert!(!snp.compress);
    assert_eq!(snp.compression_threads, 0);

    let Action::CallJoint(joint) = parse(arguments(&[
        "call",
        "joint",
        "--input",
        "reads.bam",
        "--reference",
        "reference.fa",
        "--prefix",
        "calls",
        "--meth-format",
        "cgmap",
    ]))
    .expect("joint call defaults parse") else {
        panic!("expected joint-call action");
    };
    assert!(!joint.compress);
    assert_eq!(joint.compression_threads, 0);
    assert_eq!(joint.meth_output, std::path::Path::new("calls.CGmap"));
    assert_eq!(joint.vcf_output, std::path::Path::new("calls.vcf"));

    let Action::Combine(combine) = parse(arguments(&[
        "combine",
        "--input",
        "calls.cgmap",
        "--output",
        "matrix.bed",
    ]))
    .expect("combine defaults parse") else {
        panic!("expected combine action");
    };
    assert!(!combine.compress);
    assert_eq!(combine.compression_threads, 0);
}

#[test]
fn joint_prefix_derives_compressed_bed_and_vcf_paths() {
    let Action::CallJoint(joint) = parse(arguments(&[
        "call",
        "joint",
        "--input",
        "reads.bam",
        "--reference",
        "reference.fa",
        "--prefix",
        "cohort/sample",
        "--meth-format",
        "bed",
        "--compress",
        "true",
    ]))
    .expect("compressed joint outputs parse") else {
        panic!("expected joint-call action");
    };
    assert_eq!(
        joint.meth_output,
        std::path::Path::new("cohort/sample.bed.gz")
    );
    assert_eq!(
        joint.vcf_output,
        std::path::Path::new("cohort/sample.vcf.gz")
    );
}

#[test]
fn invalid_snp_and_joint_quality_controls_are_rejected() {
    for invalid in [
        arguments(&[
            "call",
            "snp",
            "-i",
            "reads.bam",
            "--reference",
            "reference.fa",
            "-o",
            "calls.vcf",
            "-t",
            "0",
        ]),
        arguments(&[
            "call",
            "snp",
            "-i",
            "reads.bam",
            "--reference",
            "reference.fa",
            "-o",
            "calls.vcf",
            "--min-bq",
            "94",
        ]),
        arguments(&[
            "call",
            "snp",
            "-i",
            "reads.bam",
            "--reference",
            "reference.fa",
            "-o",
            "calls.vcf",
            "--underconversion-rate",
            "1.000000001",
        ]),
        arguments(&[
            "call",
            "snp",
            "-i",
            "reads.bam",
            "--reference",
            "reference.fa",
            "-o",
            "calls.vcf",
            "--heterozygosity",
            "0",
        ]),
        arguments(&[
            "call",
            "joint",
            "-i",
            "reads.bam",
            "-p",
            "calls",
            "-f",
            "cgmap",
            "--heterozygosity",
            "1",
        ]),
    ] {
        assert!(parse(invalid).is_err());
    }
}

#[test]
fn retired_min_base_quality_spelling_is_rejected() {
    for invalid in [
        arguments(&[
            "call",
            "meth",
            "-i",
            "reads.bam",
            "-r",
            "reference.fa",
            "-o",
            "calls.CGmap",
            "-f",
            "cgmap",
            "--min-base-quality",
            "20",
        ]),
        arguments(&[
            "call",
            "snp",
            "-i",
            "reads.bam",
            "-r",
            "reference.fa",
            "-o",
            "calls.vcf",
            "--min-base-quality",
            "20",
        ]),
        arguments(&[
            "call",
            "joint",
            "-i",
            "reads.bam",
            "-r",
            "reference.fa",
            "-p",
            "calls",
            "-f",
            "cgmap",
            "--min-base-quality",
            "20",
        ]),
    ] {
        assert!(parse(invalid).is_err());
    }
}

#[test]
fn alignment_entry_points_are_explicit() {
    let parsed = parse(arguments(&[
        "align",
        "--index",
        "reference.bsbit",
        "--read1",
        "r1.fq",
        "-2",
        "r2.fq",
        "--output",
        "out.bam",
    ]))
    .expect("canonical alignment parses");
    assert!(matches!(parsed, Action::Align(_)));

    let single = parse(arguments(&[
        "align",
        "--index",
        "reference.bsbit",
        "--read1",
        "single.fq",
        "--output",
        "single.bam",
    ]))
    .expect("canonical single-end alignment parses");
    assert!(matches!(single, Action::Align(_)));
}

#[test]
fn paired_span_and_fail_closed_rules_are_exact() {
    let paired = arguments(&[
        "align",
        "--index",
        "ref.bsbit",
        "-1",
        "r1.fq",
        "-2",
        "r2.fq",
        "--output",
        "out.bam",
        "--min-template-span",
        "10",
        "--max-template-span",
        "500",
        "--batch-size",
        "17",
    ]);
    assert!(matches!(parse(paired), Ok(Action::Align(_))));

    for invalid in [
        arguments(&[
            "align",
            "--index",
            "ref.bsbit",
            "--read2",
            "r2",
            "--output",
            "o.bam",
        ]),
        arguments(&[
            "align",
            "--index",
            "x",
            "--read1",
            "r1",
            "--read2",
            "r2",
            "--output",
            "o.bam",
            "--min-template-span",
            "501",
            "--max-template-span",
            "500",
        ]),
    ] {
        assert!(parse(invalid).is_err());
    }
}

#[test]
fn thread_domain_is_exact_and_independent_of_host_cpu_count() {
    let base = [
        "align",
        "--index",
        "ref.bsbit",
        "--read1",
        "reads.fq",
        "--output",
        "out.bam",
        "--threads",
    ];
    for value in ["1", "2", "64", "65", "1024"] {
        let mut supplied = base.to_vec();
        supplied.push(value);
        let Action::Align(_) = parse(arguments(&supplied)).expect("threads parse") else {
            panic!("expected align action");
        };
    }
    for value in ["0", "4294967296", "18446744073709551615", "many"] {
        let mut supplied = base.to_vec();
        supplied.push(value);
        assert!(parse(arguments(&supplied)).is_err(), "accepted {value}");
    }
    assert!(
        parse(arguments(&base)).is_err(),
        "accepted a missing thread value"
    );
    let mut duplicate = base.to_vec();
    duplicate.extend(["1", "--threads", "2"]);
    assert!(
        parse(arguments(&duplicate)).is_err(),
        "accepted duplicate thread values"
    );

    assert!(matches!(
        parse(arguments(&[
            "index",
            "--reference",
            "ref.fa",
            "--output",
            "ref.bsbit",
            "--threads",
            "65",
        ])),
        Ok(Action::Index(_))
    ));
    assert!(matches!(
        parse(arguments(&[
            "call",
            "meth",
            "--input",
            "reads.bam",
            "--reference",
            "ref.fa",
            "--output",
            "calls.cgmap",
            "--format",
            "cgmap",
            "--threads",
            "65",
        ])),
        Ok(Action::CallMeth(_))
    ));
    assert!(matches!(
        parse(arguments(&[
            "combine",
            "--input",
            "calls.cgmap",
            "--output",
            "matrix.bed",
            "--threads",
            "65",
        ])),
        Ok(Action::Combine(_))
    ));
}
