//! White-box tests for the standard alignment command orchestration.
//!
//! Kept outside implementation `src/` while remaining a child module so private
//! invariants can be tested without widening the crate API.

use super::{
    AlignmentAuxiliaryMode, HELP, MetricsTimer, PairedLibraryProfile, PairedSearchMode, ReadLayout,
    ReadOutputMode, parse_options_from, sensitive_mapq_zero_strategy_id,
    sensitive_read_complete_strategy_id, strategy_id, throughput_thread_split,
};
use std::ffi::OsString;

#[test]
fn disabled_metrics_timer_never_starts_a_clock() {
    let timer = MetricsTimer::start(false);
    assert!(timer.0.is_none());
    assert_eq!(timer.elapsed_ns(), 0);
}

#[test]
fn help_exposes_only_default_and_sensitive_modes() {
    assert!(HELP.contains("Without --sensitive, default mode"));
    assert!(HELP.contains("--sensitive"));
    assert!(HELP.contains("--mapped-only"));
    assert!(HELP.contains("minimal|bismark"));
    assert!(HELP.contains("--non-directional"));
    assert!(HELP.contains("--read1 only"));
    assert!(HELP.contains("-1, --read1 PATH"));
    assert!(HELP.contains("-2, --read2 PATH"));
    assert!(HELP.contains("same persisted combined index and bounded d3/d5"));
    for hidden in [
        "align-general",
        "--fast",
        "--repeat-frontier-v1",
        "--repeat-frontier-v2",
        "--sensitive-audited",
        "--mapq-zero-output",
        "--mapq-policy",
        "--soft-clip-fallback",
        "--mate-rescue",
        "--insert-prior",
        "--records",
        "--skip-records",
        "--staging-bam",
        "--packed-reference-catalog",
        "--expected-reference-digest",
        "--combined-index-prefix",
        "--reads PATH",
        "--reads1",
        "--reads2",
    ] {
        assert!(!HELP.contains(hidden), "public help leaked {hidden}");
    }
}

#[test]
fn standard_single_input_is_first_class_and_pair_only_options_fail_closed() {
    let single = [
        "--index",
        "reference.bsbit",
        "--read1",
        "reads.fastq.gz",
        "--output-bam",
        "output.bam",
    ]
    .map(OsString::from)
    .to_vec();
    let parsed = parse_options_from(single.clone()).expect("single-end input parses");
    assert_eq!(parsed.layout, ReadLayout::SingleEnd);
    assert!(parsed.read2.is_none());
    assert_eq!(parsed.threads, 1);
    assert_eq!(parsed.bam_threads, 1);
    assert_eq!(parsed.auxiliary_core_budget, None);
    assert_eq!(parsed.bam_compression_level, Some(1));

    let mut sensitive = single.clone();
    sensitive.push(OsString::from("--sensitive"));
    let parsed = parse_options_from(sensitive).expect("single-end sensitive input parses");
    assert_eq!(parsed.layout, ReadLayout::SingleEnd);
    assert_eq!(parsed.search_mode, PairedSearchMode::Sensitive);

    let mut threaded = single.clone();
    threaded.extend(["--threads", "4"].map(OsString::from));
    assert_eq!(
        parse_options_from(threaded)
            .expect("single-end threads parse")
            .threads,
        4
    );

    let mut explicit_bam = single.clone();
    explicit_bam
        .extend(["--bam-threads", "2", "--bam-compression-level", "default"].map(OsString::from));
    let parsed = parse_options_from(explicit_bam).expect("single-end BAM controls parse");
    assert_eq!(parsed.bam_threads, 2);
    assert_eq!(parsed.bam_compression_level, None);

    for arguments in [
        [
            single.clone(),
            vec![
                OsString::from("--output-contract"),
                OsString::from("bismark"),
            ],
        ]
        .concat(),
        [single, vec![OsString::from("--metrics")]].concat(),
    ] {
        assert!(
            parse_options_from(arguments)
                .expect_err("paired-only option must reject single input")
                .to_string()
                .contains("requires paired input via --read2")
        );
    }
}

#[test]
fn paired_total_thread_budget_selects_qualified_mapping_output_splits() {
    assert_eq!(throughput_thread_split(1), (1, 0));
    assert_eq!(throughput_thread_split(2), (1, 1));
    assert_eq!(throughput_thread_split(10), (8, 2));
    assert_eq!(throughput_thread_split(14), (11, 3));
    assert_eq!(throughput_thread_split(64), (60, 4));

    let paired = || {
        [
            "--index",
            "reference.bsbit",
            "--read1",
            "r1.fastq.gz",
            "--read2",
            "r2.fastq.gz",
            "--output-bam",
            "output.bam",
        ]
        .map(OsString::from)
        .to_vec()
    };
    let mut automatic = paired();
    automatic.extend(["--total-threads", "14"].map(OsString::from));
    let automatic = parse_options_from(automatic).expect("total budget parses");
    assert_eq!(automatic.threads, 11);
    assert_eq!(automatic.bam_threads, 3);
    assert_eq!(automatic.auxiliary_core_budget, Some(3));

    for explicit in ["--threads", "--bam-threads"] {
        let mut conflicting = paired();
        conflicting.extend(["--total-threads", "14", explicit, "2"].map(OsString::from));
        assert_eq!(
            parse_options_from(conflicting)
                .expect_err("automatic and explicit thread controls conflict")
                .to_string(),
            "--total-threads conflicts with --threads and --bam-threads"
        );
    }
}

#[test]
fn read_layout_accepts_canonical_short_forms_and_rejects_retired_spellings() {
    let short_single = [
        "--index",
        "reference.bsbit",
        "-1",
        "single.fastq.gz",
        "--output-bam",
        "single.bam",
    ]
    .map(OsString::from);
    let parsed = parse_options_from(short_single).expect("-1 selects single-end input");
    assert_eq!(parsed.layout, ReadLayout::SingleEnd);
    assert_eq!(parsed.read1, std::path::PathBuf::from("single.fastq.gz"));
    assert!(parsed.read2.is_none());

    let short_pair = [
        "--index",
        "reference.bsbit",
        "-1",
        "r1.fastq.gz",
        "-2",
        "r2.fastq.gz",
        "--output-bam",
        "paired.bam",
    ]
    .map(OsString::from);
    let parsed = parse_options_from(short_pair).expect("-1 and -2 select paired input");
    assert_eq!(parsed.layout, ReadLayout::PairedEnd);
    assert_eq!(parsed.read1, std::path::PathBuf::from("r1.fastq.gz"));
    assert_eq!(
        parsed.read2.as_deref(),
        Some(std::path::Path::new("r2.fastq.gz"))
    );

    for retired in ["--reads1", "--reads2", "--reads"] {
        let arguments = [
            "--index",
            "reference.bsbit",
            retired,
            "reads.fastq.gz",
            "--output-bam",
            "output.bam",
        ]
        .map(OsString::from);
        assert_eq!(
            parse_options_from(arguments)
                .expect_err("retired read flag must be rejected")
                .to_string(),
            format!("unknown option {retired}")
        );
    }

    let duplicate_alias = [
        "--index",
        "reference.bsbit",
        "--read1",
        "first.fastq.gz",
        "-1",
        "second.fastq.gz",
        "--output-bam",
        "output.bam",
    ]
    .map(OsString::from);
    assert_eq!(
        parse_options_from(duplicate_alias)
            .expect_err("short and long forms identify one option")
            .to_string(),
        "--read1 may be specified only once"
    );

    let read2_only = [
        "--index",
        "reference.bsbit",
        "-2",
        "r2.fastq.gz",
        "--output-bam",
        "output.bam",
    ]
    .map(OsString::from);
    assert_eq!(
        parse_options_from(read2_only)
            .expect_err("read 2 alone is not a valid layout")
            .to_string(),
        "--read2 requires --read1"
    );
}

#[test]
fn minimal_is_default_and_bismark_output_is_explicit() {
    let required = || {
        [
            "--index",
            "reference.bsbit",
            "--read1",
            "reads.R1.fastq.gz",
            "--read2",
            "reads.R2.fastq.gz",
            "--output-bam",
            "output.bam",
        ]
        .map(OsString::from)
        .to_vec()
    };

    let defaults = parse_options_from(required()).expect("default output contract");
    assert_eq!(defaults.output_contract, AlignmentAuxiliaryMode::Minimal);
    assert_eq!(defaults.library_profile, PairedLibraryProfile::Directional);
    let mut compatible = required();
    compatible.extend(["--output-contract", "bismark"].map(OsString::from));
    assert_eq!(
        parse_options_from(compatible)
            .expect("Bismark output contract")
            .output_contract,
        AlignmentAuxiliaryMode::Bismark
    );

    let mut non_directional = required();
    non_directional.push(OsString::from("--non-directional"));
    let non_directional = parse_options_from(non_directional).expect("non-directional defaults");
    assert_eq!(
        non_directional.library_profile,
        PairedLibraryProfile::NonDirectional
    );
    assert_eq!(
        non_directional.output_contract,
        AlignmentAuxiliaryMode::Minimal
    );

    let mut retired_alias = required();
    retired_alias.push(OsString::from("--non_directional"));
    assert!(
        parse_options_from(retired_alias)
            .expect_err("retired Bismark spelling must be rejected")
            .to_string()
            .contains("unknown option")
    );

    let mut explicit_minimal = required();
    explicit_minimal
        .extend(["--non-directional", "--output-contract", "minimal"].map(OsString::from));
    assert_eq!(
        parse_options_from(explicit_minimal)
            .expect("non-directional minimal contract")
            .output_contract,
        AlignmentAuxiliaryMode::Minimal
    );

    let mut legacy = required();
    legacy.extend(["--output-contract", "nm-md"].map(OsString::from));
    assert_eq!(
        parse_options_from(legacy)
            .expect_err("legacy contract spelling is rejected")
            .to_string(),
        "invalid --output-contract; expected minimal or bismark"
    );
}

#[test]
fn parser_accepts_only_the_opaque_index_handle_and_rejects_duplicate_value_flags() {
    let snapshot = [
        "--index",
        "reference.bsbit",
        "--read1",
        "r1.fq",
        "--read2",
        "r2.fq",
        "--output-bam",
        "out.bam",
    ]
    .map(OsString::from)
    .to_vec();
    let options = parse_options_from(snapshot.clone()).expect("index form");
    assert_eq!(options.index, std::path::PathBuf::from("reference.bsbit"));
    assert!(!options.emit_metrics);

    for hidden in [
        "--packed-reference-catalog",
        "--expected-reference-digest",
        "--combined-index-prefix",
    ] {
        let mut arguments = snapshot.clone();
        arguments.extend([hidden, "internal"].map(OsString::from));
        assert_eq!(
            parse_options_from(arguments)
                .expect_err("internal option must stay hidden")
                .to_string(),
            format!("unknown option {hidden}")
        );
    }

    let mut metrics = snapshot.clone();
    metrics.push(OsString::from("--metrics"));
    assert!(
        parse_options_from(metrics)
            .expect("metrics opt in")
            .emit_metrics
    );

    let mut duplicate = snapshot.clone();
    duplicate.extend(["--threads", "2", "--threads", "3"].map(OsString::from));
    assert_eq!(
        parse_options_from(duplicate)
            .expect_err("duplicate value flag")
            .to_string(),
        "--threads may be specified only once"
    );

    let mut inverted = snapshot;
    inverted.extend(["--min-template-span", "20", "--max-template-span", "10"].map(OsString::from));
    assert_eq!(
        parse_options_from(inverted)
            .expect_err("inverted template bounds")
            .to_string(),
        "--min-template-span must not exceed --max-template-span"
    );
}

#[test]
fn sensitive_strategy_ids_describe_only_supported_output_modes() {
    let mapped_only = sensitive_mapq_zero_strategy_id();
    let complete = sensitive_read_complete_strategy_id();
    assert!(mapped_only.ends_with("-mapq0-all-v1"));
    assert!(complete.ends_with("-read-complete-v1"));
    assert_ne!(mapped_only, complete);
    assert_eq!(mapped_only, "sensitive-bounded-integrated-mapq0-all-v1");
    assert_eq!(complete, "sensitive-bounded-integrated-read-complete-v1");
}
#[test]
fn public_modes_select_fixed_default_and_sensitive_strategies() {
    let required = || {
        [
            "--index",
            "reference.bsbit",
            "--read1",
            "reads.R1.fastq.gz",
            "--read2",
            "reads.R2.fastq.gz",
            "--output-bam",
            "output.bam",
        ]
        .map(OsString::from)
        .to_vec()
    };

    let default = parse_options_from(required()).expect("default options");
    assert_eq!(
        strategy_id(&default),
        "balanced-d5-adapter-recovery-read-complete-v2"
    );
    assert_eq!(default.search_mode, PairedSearchMode::Default);
    assert_eq!(default.read_output, ReadOutputMode::Complete);

    let mut default_mapped_only = required();
    default_mapped_only.push(OsString::from("--mapped-only"));
    let default_mapped_only =
        parse_options_from(default_mapped_only).expect("mapped-only default output");
    assert_eq!(
        strategy_id(&default_mapped_only),
        "balanced-d5-adapter-recovery-mapq0-all-v2"
    );
    assert_eq!(default_mapped_only.read_output, ReadOutputMode::MappedOnly);

    let mut sensitive = required();
    sensitive.push(OsString::from("--sensitive"));
    let sensitive = parse_options_from(sensitive).expect("sensitive defaults");
    assert_eq!(
        strategy_id(&sensitive),
        sensitive_read_complete_strategy_id()
    );
    assert_eq!(sensitive.search_mode, PairedSearchMode::Sensitive);
    assert_eq!(sensitive.read_output, ReadOutputMode::Complete);

    let mut mapped_only = required();
    mapped_only.extend(["--sensitive", "--mapped-only"].map(OsString::from));
    let mapped_only = parse_options_from(mapped_only).expect("mapped-only sensitive output");
    assert_eq!(mapped_only.read_output, ReadOutputMode::MappedOnly);
    assert_eq!(strategy_id(&mapped_only), sensitive_mapq_zero_strategy_id());

    let mut duplicate = required();
    duplicate.extend(["--sensitive", "--sensitive"].map(OsString::from));
    assert_eq!(
        parse_options_from(duplicate)
            .expect_err("sensitive may be selected only once")
            .to_string(),
        "--sensitive may be specified only once"
    );

    let mut duplicate_mapped_only = required();
    duplicate_mapped_only.extend(["--mapped-only", "--mapped-only"].map(OsString::from));
    assert_eq!(
        parse_options_from(duplicate_mapped_only)
            .expect_err("mapped-only may be selected only once")
            .to_string(),
        "--mapped-only may be specified only once"
    );

    let mut retired_mapq_policy = required();
    retired_mapq_policy.extend(["--mapq-policy", "qualified"].map(OsString::from));
    assert_eq!(
        parse_options_from(retired_mapq_policy)
            .expect_err("the retired MAPQ-policy selector must stay unavailable")
            .to_string(),
        "unknown option --mapq-policy"
    );

    for hidden in [
        "--fast",
        "--repeat-frontier-v1",
        "--repeat-frontier-v2",
        "--sensitive-audited",
        "--max-edit-distance",
        "--soft-clip-fallback",
        "--mate-rescue",
        "--insert-prior",
        "--records",
        "--skip-records",
        "--staging-bam",
    ] {
        let mut args = required();
        args.push(OsString::from(hidden));
        assert_eq!(
            parse_options_from(args)
                .expect_err("development-only option must not enter a implementation build")
                .to_string(),
            format!("unknown option {hidden}")
        );
    }
}
