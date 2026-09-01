use super::index::IndexSpeed;
use super::{Action, parse};
use bsbit_call::meth::OutputFormat as MethylationOutputFormat;
use bsbit_combine::MatrixFormat as CombineMatrixFormat;

fn arguments(values: &[&str]) -> Vec<std::ffi::OsString> {
    values.iter().map(std::ffi::OsString::from).collect()
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
    assert_eq!(index.speed, IndexSpeed::Balanced);
}

#[test]
fn internal_index_construction_is_not_a_public_subcommand() {
    assert!(parse(arguments(&["index", "combined", "--snapshot", "ref.bsbit"])).is_err());
}

#[test]
fn public_help_exposes_only_supported_index_and_alignment_entry_points() {
    let help = crate::GENERAL_HELP;
    assert!(help.contains("bsbit index"));
    assert!(help.contains("bsbit align --index"));
    for hidden in [
        "align-general",
        "index combined",
        "build-bsbit",
        "bsbit-align",
        "bsbit-call",
        "bsbit cache",
        "cache-gc",
    ] {
        assert!(!help.contains(hidden), "public help leaked {hidden}");
    }
}

#[test]
fn standard_alignment_layout_and_read_aliases_parse() {
    let Action::Align(_) = parse(arguments(&[
        "align",
        "--index",
        "ref.bsbit",
        "-1",
        "reads.fq",
        "--output-bam",
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
        "--output-bam",
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
        "--output-bam",
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
        "--reference",
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
        "--regions-file",
        "targets.bed.gz",
    ]))
    .expect("short call options parse") else {
        panic!("expected methylation-call action");
    };
    assert_eq!(short.format, MethylationOutputFormat::Cgmap);
    assert!(short.compress);
    assert_eq!(short.threads, 1);
    assert_eq!(short.parameters.minimum_base_quality, 15);
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
        "--threads",
        "4",
        "--min-base-quality",
        "25",
        "--min-mapq",
        "30",
    ]))
    .expect("long call options parse") else {
        panic!("expected methylation-call action");
    };
    assert_eq!(long.format, MethylationOutputFormat::Bed);
    assert_eq!(long.reference, std::path::Path::new("reference.fa"));
    assert!(!long.compress);
    assert_eq!(long.threads, 4);
    assert_eq!(long.parameters.minimum_base_quality, 25);
    assert_eq!(long.parameters.minimum_mapping_quality, 30);
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
            "--min-base-quality",
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
fn all_call_modules_require_an_indexed_reference_argument() {
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
            "-m",
            "calls.cgmap",
            "-f",
            "cgmap",
            "-v",
            "calls.vcf",
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
        "matrix.bed.gz",
        "-m",
        "both",
        "--min-count",
        "10",
        "--min-prop",
        "0.8",
        "-c",
        "true",
        "-t",
        "8",
    ]))
    .expect("combine options parse") else {
        panic!("expected combine action");
    };
    assert_eq!(options.inputs.len(), 2);
    assert_eq!(options.inputs[0].sample, "tumor");
    assert_eq!(options.inputs[1].sample, "normal");
    assert_eq!(options.matrix_format, CombineMatrixFormat::Both);
    assert_eq!(options.parameters.minimum_count, 10);
    assert_eq!(
        options
            .parameters
            .minimum_sample_proportion_parts_per_billion,
        800_000_000
    );
    assert!(options.compress);
    assert_eq!(options.threads, 8);
}

#[test]
fn methylation_combine_defaults_names_to_exact_paths() {
    let Action::Combine(defaulted) = parse(arguments(&[
        "combine",
        "--input",
        "cohort/tumor.bed.gz,cohort/normal sample.bed",
        "--output",
        "matrix.bed",
    ]))
    .expect("path-derived sample names parse") else {
        panic!("expected combine action");
    };
    assert_eq!(defaulted.inputs.len(), 2);
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
        "--reference",
        "reference.fa",
        "-o",
        "calls.vcf.gz",
        "--sample-name",
        "tumor-A",
        "-c",
        "true",
        "-t",
        "8",
        "--min-base-quality",
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
    ]))
    .expect("SNP call parses") else {
        panic!("expected SNP-call action");
    };
    assert!(snp.compress);
    assert_eq!(snp.reference, std::path::Path::new("reference.fa"));
    assert_eq!(snp.sample_name.as_deref(), Some("tumor-A"));
    assert_eq!(snp.threads, 8);
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

    let Action::CallJoint(joint) = parse(arguments(&[
        "call",
        "joint",
        "-i",
        "reads.bam",
        "--reference",
        "reference.fa",
        "--sample-name",
        "tumor-A",
        "-m",
        "calls.cgmap.gz",
        "-f",
        "cgmap",
        "-v",
        "calls.vcf.gz",
        "--heterozygosity",
        "0.002",
    ]))
    .expect("joint call parses") else {
        panic!("expected joint-call action");
    };
    assert_eq!(joint.meth_format, MethylationOutputFormat::Cgmap);
    assert_eq!(joint.reference, std::path::Path::new("reference.fa"));
    assert_eq!(joint.sample_name.as_deref(), Some("tumor-A"));
    assert_eq!(joint.parameters.minimum_base_quality, 15);
    assert_eq!(joint.parameters.minimum_mapping_quality, 20);
    assert_eq!(joint.parameters.heterozygosity_parts_per_billion, 2_000_000);
    assert!(!joint.compress);
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
            "--min-base-quality",
            "94",
        ]),
        arguments(&[
            "call",
            "joint",
            "-i",
            "reads.bam",
            "-m",
            "same.gz",
            "-f",
            "cgmap",
            "-v",
            "same.gz",
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
            "-m",
            "calls.cgmap",
            "-f",
            "cgmap",
            "-v",
            "calls.vcf",
            "--heterozygosity",
            "1",
        ]),
    ] {
        assert!(parse(invalid).is_err());
    }
}

#[test]
fn retired_joint_output_aliases_are_rejected() {
    for retired in ["--meth-output", "--vcf-output"] {
        let supplied = [
            "call",
            "joint",
            "--input",
            "reads.bam",
            "--reference",
            "reference.fa",
            "--meth",
            "calls.cgmap",
            "--meth-format",
            "cgmap",
            "--vcf",
            "calls.vcf",
            retired,
            "retired-output",
        ];
        assert!(parse(arguments(&supplied)).is_err());
    }
}

#[test]
fn internal_cache_and_historical_alignment_options_are_absent() {
    for removed in ["cache", "cache-gc"] {
        assert!(parse(arguments(&[removed, "--index", "ref.bsbit"])).is_err());
    }
    assert!(parse(arguments(&["align-general", "--help"])).is_err());

    for retired in [
        "--reads",
        "--output",
        "--output-format",
        "--max-edit-distance",
        "--reference-backend",
        "--packed-cache",
        "--debug",
    ] {
        let supplied = [
            "align",
            "--index",
            "ref.bsbit",
            "--read1",
            "reads.fq",
            "--output-bam",
            "out.bam",
            retired,
            "value",
        ];
        assert!(parse(arguments(&supplied)).is_err(), "accepted {retired}");
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
        "--output-bam",
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
        "--output-bam",
        "single.bam",
    ]))
    .expect("canonical single-end alignment parses");
    assert!(matches!(single, Action::Align(_)));

    let implicit = parse(arguments(&[
        "align",
        "--index",
        "reference.bsbit",
        "--reads",
        "reads.fq",
        "--output",
        "out.sam",
        "--output-format",
        "sam",
        "--max-edit-distance",
        "1",
    ]));
    assert!(
        implicit.is_err(),
        "standard align must reject historical option names"
    );
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
        "--output-bam",
        "out.bam",
        "--min-template-span",
        "10",
        "--max-template-span",
        "500",
        "--batch-pairs",
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
            "--output-bam",
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
            "--output-bam",
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
        "--output-bam",
        "out.bam",
        "--threads",
    ];
    for value in ["1", "2", "64"] {
        let mut supplied = base.to_vec();
        supplied.push(value);
        let Action::Align(_) = parse(arguments(&supplied)).expect("threads parse") else {
            panic!("expected align action");
        };
    }
    for value in ["0", "65", "18446744073709551615", "many"] {
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
}
