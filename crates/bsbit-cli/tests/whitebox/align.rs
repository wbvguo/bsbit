//! White-box tests for the standard alignment command orchestration.
//!
//! Kept outside implementation `src/` while remaining a child module so private
//! invariants can be tested without widening the crate API.

use super::{
    HELP, MetricsTimer, ReadLayout, ReadOutputMode, SearchMode, parse_options_from,
    sensitive_mapq_zero_strategy_id, sensitive_read_complete_strategy_id, strategy_id,
    throughput_thread_split,
};
use bsbit_align::library::LibraryProfile;
use bsbit_cpu::BackendRequest;
use bsbit_hts::AlignmentAuxiliaryMode;
use std::ffi::OsString;

#[test]
fn disabled_metrics_timer_never_starts_a_clock() {
    let timer = MetricsTimer::start(false);
    assert!(timer.0.is_none());
    assert_eq!(timer.elapsed_ns(), 0);
}

#[cfg(unix)]
#[test]
fn alignment_parser_preserves_non_utf8_filesystem_paths() {
    use std::os::unix::ffi::{OsStrExt, OsStringExt};

    let raw = |stem: &[u8]| {
        let mut bytes = stem.to_vec();
        bytes.push(0xff);
        OsString::from_vec(bytes)
    };
    let index = raw(b"index-");
    let read = raw(b"read-");
    let output = raw(b"output-");
    let parsed = parse_options_from(vec![
        OsString::from("--index"),
        index.clone(),
        OsString::from("--read1"),
        read.clone(),
        OsString::from("--output"),
        output.clone(),
    ])
    .expect("filesystem path bytes do not need to be UTF-8");
    assert_eq!(parsed.index.as_os_str().as_bytes(), index.as_bytes());
    assert_eq!(parsed.read1.as_os_str().as_bytes(), read.as_bytes());
    assert_eq!(parsed.output.as_os_str().as_bytes(), output.as_bytes());
}

#[test]
fn help_describes_supported_alignment_modes_and_inputs() {
    assert!(HELP.contains("Without --sensitive, default mode"));
    assert!(HELP.contains("--sensitive"));
    assert!(HELP.contains("--mapped-only"));
    assert!(HELP.contains("minimal|bismark"));
    assert!(HELP.contains("--non-directional"));
    assert!(HELP.contains("--read1 only"));
    assert!(HELP.contains("-1, --read1 PATH"));
    assert!(HELP.contains("-2, --read2 PATH"));
    assert!(HELP.contains("same persisted combined index and bounded d3/d5"));
    assert!(HELP.contains("does not apply"));
    assert!(HELP.contains("profile-specific MAPQ promotion tables"));
    assert!(HELP.contains("A recovered pair that cannot earn"));
    assert!(HELP.contains("coordinate-bearing, ambiguous, and MAPQ 0"));
}

#[test]
fn standard_single_input_supports_shared_controls_and_pair_only_options_fail_closed() {
    let single = [
        "--index",
        "reference.bsbit",
        "--read1",
        "reads.fastq.gz",
        "--output",
        "output.bam",
    ]
    .map(OsString::from)
    .to_vec();
    let parsed = parse_options_from(single.clone()).expect("single-end input parses");
    assert_eq!(parsed.layout(), ReadLayout::SingleEnd);
    assert!(parsed.read2.is_none());
    assert_eq!(parsed.threads, 1);
    assert_eq!(parsed.compression_threads, 1);
    assert_eq!(parsed.compression_level, Some(1));
    assert_eq!(parsed.batch_size, 1_000);
    assert_eq!(parsed.queue_batches, 2);
    assert_eq!(parsed.tie_break_seed, 0);
    assert_eq!(parsed.maximum_edit_distance, 5);
    assert_eq!(parsed.simd_backend, BackendRequest::Auto);

    let mut threaded = single.clone();
    threaded.extend(["--threads", "4"].map(OsString::from));
    assert_eq!(
        parse_options_from(threaded)
            .expect("single-end threads parse")
            .threads,
        4
    );

    let mut forced = single.clone();
    forced.extend(["--simd-backend", "scalar"].map(OsString::from));
    assert_eq!(
        parse_options_from(forced)
            .expect("forced baseline parses")
            .simd_backend,
        BackendRequest::Scalar
    );

    let mut forced_neon = single.clone();
    forced_neon.extend(["--simd-backend", "neon"].map(OsString::from));
    assert_eq!(
        parse_options_from(forced_neon)
            .expect("forced NEON backend parses")
            .simd_backend,
        BackendRequest::Neon
    );

    let mut invalid_backend = single.clone();
    invalid_backend.extend(["--simd-backend", "native"].map(OsString::from));
    assert!(
        parse_options_from(invalid_backend)
            .expect_err("unbounded native backend is rejected")
            .to_string()
            .contains("expected auto, scalar, sse2, sse4.2, avx2, avx512, or neon")
    );

    let mut sensitive = single.clone();
    sensitive.push(OsString::from("--sensitive"));
    let parsed = parse_options_from(sensitive).expect("single-end sensitive input parses");
    assert_eq!(parsed.layout(), ReadLayout::SingleEnd);
    assert_eq!(parsed.search_mode, SearchMode::Sensitive);

    let mut total_budget = single.clone();
    total_budget.extend(["--total-threads", "10"].map(OsString::from));
    let parsed = parse_options_from(total_budget).expect("single-end total budget parses");
    assert_eq!(parsed.total_thread_budget, Some(10));

    let mut seeded = single.clone();
    seeded.extend(["--tie-break-seed", "18446744073709551615"].map(OsString::from));
    assert_eq!(
        parse_options_from(seeded)
            .expect("full-width reporting tie-break seed parses")
            .tie_break_seed,
        u64::MAX
    );

    let mut shared = single.clone();
    shared.extend(
        [
            "--non-directional",
            "--output-contract",
            "bismark",
            "--mapped-only",
            "--metrics",
        ]
        .map(OsString::from),
    );
    let parsed = parse_options_from(shared).expect("single-end shared controls parse");
    assert_eq!(parsed.library_profile, LibraryProfile::NonDirectional);
    assert_eq!(parsed.output_contract, AlignmentAuxiliaryMode::Bismark);
    assert_eq!(parsed.read_output, ReadOutputMode::MappedOnly);
    assert!(parsed.emit_metrics);

    for pair_only in ["--min-template-span", "--max-template-span"] {
        let arguments = [
            single.clone(),
            vec![OsString::from(pair_only), OsString::from("2")],
        ]
        .concat();
        assert!(
            parse_options_from(arguments)
                .expect_err("paired-only option must reject single input")
                .to_string()
                .contains("requires paired input via --read2")
        );
    }
}

#[test]
fn edit_budget_is_explicit_bounded_and_shared_by_both_layouts() {
    let required = [
        "--index",
        "reference.bsbit",
        "--read1",
        "reads.fastq.gz",
        "--output",
        "output.bam",
    ]
    .map(OsString::from)
    .to_vec();
    let mut bounded = required.clone();
    bounded.extend(["--max-edit-distance", "3"].map(OsString::from));
    assert_eq!(
        parse_options_from(bounded)
            .expect("single-end edit budget parses")
            .maximum_edit_distance,
        3
    );

    let mut excessive = required.clone();
    excessive.extend(["--max-edit-distance", "6"].map(OsString::from));
    assert_eq!(
        parse_options_from(excessive)
            .expect_err("edit budget above the implementation bound is rejected")
            .to_string(),
        "--max-edit-distance must be in 0..=5"
    );

    let mut paired = required;
    paired.extend(["--read2", "reads2.fastq.gz", "--max-edit-distance", "3"].map(OsString::from));
    let paired = parse_options_from(paired).expect("paired-end edit budget parses");
    assert_eq!(paired.layout(), ReadLayout::PairedEnd);
    assert_eq!(paired.maximum_edit_distance, 3);
}

#[test]
fn adapter_and_soft_clip_policy_is_shared_validated_and_reproducible() {
    let required = [
        "--index",
        "reference.bsbit",
        "--read1",
        "reads.fastq.gz",
        "--output",
        "output.bam",
    ]
    .map(OsString::from)
    .to_vec();
    let defaults = parse_options_from(required.clone()).expect("default clipping policy parses");
    assert_eq!(
        defaults.output_policy.adapter(),
        bsbit_align::AlignmentOutputPolicy::default().adapter()
    );
    assert_eq!(defaults.output_policy.adapter_minimum_overlap(), 8);
    assert_eq!(defaults.output_policy.adapter_maximum_clip_bases(), 30);
    assert_eq!(
        defaults.output_policy.soft_clip_mode(),
        bsbit_align::SoftClipMode::Auto
    );
    assert_eq!(defaults.output_policy.maximum_soft_clip_bases(), 30);

    let mut custom = required.clone();
    custom.extend(
        [
            "--adapter",
            "acgtacgt",
            "--adapter-min-overlap",
            "4",
            "--adapter-max-clip",
            "11",
            "--soft-clip",
            "adapter",
            "--max-soft-clip",
            "12",
        ]
        .map(OsString::from),
    );
    let custom = parse_options_from(custom).expect("custom clipping policy parses");
    assert_eq!(custom.output_policy.adapter(), Some(b"ACGTACGT".as_slice()));
    assert_eq!(custom.output_policy.adapter_minimum_overlap(), 4);
    assert_eq!(custom.output_policy.adapter_maximum_clip_bases(), 11);
    assert_eq!(
        custom.output_policy.soft_clip_mode(),
        bsbit_align::SoftClipMode::Adapter
    );
    assert_eq!(custom.output_policy.maximum_soft_clip_bases(), 12);

    let mut disabled = required.clone();
    disabled.extend(
        [
            "--adapter",
            "none",
            "--soft-clip",
            "none",
            "--max-soft-clip",
            "0",
        ]
        .map(OsString::from),
    );
    let disabled = parse_options_from(disabled).expect("disabled clipping policy parses");
    assert_eq!(disabled.output_policy.adapter(), None);
    assert_eq!(
        disabled.output_policy.soft_clip_mode(),
        bsbit_align::SoftClipMode::None
    );

    let mut invalid_adapter = required.clone();
    invalid_adapter.extend(["--adapter", "ACGT-X"].map(OsString::from));
    assert!(
        parse_options_from(invalid_adapter)
            .expect_err("non-DNA adapter is rejected")
            .to_string()
            .contains("outside A/C/G/T/N")
    );

    let mut excessive_clip = required;
    excessive_clip.extend(["--max-soft-clip", "193"].map(OsString::from));
    assert!(
        parse_options_from(excessive_clip)
            .expect_err("clip bound outside the read domain is rejected")
            .to_string()
            .contains("exceeds 192")
    );
}

#[test]
fn single_input_accepts_canonical_compression_and_pipeline_controls() {
    let parsed = parse_options_from(
        [
            "--index",
            "reference.bsbit",
            "--read1",
            "reads.fastq.gz",
            "--output",
            "output.bam",
            "--compression-threads",
            "2",
            "--compression-level",
            "default",
            "--batch-size",
            "2048",
            "--queue-batches",
            "3",
        ]
        .map(OsString::from),
    )
    .expect("shared pipeline controls parse");
    assert_eq!(parsed.compression_threads, 2);
    assert_eq!(parsed.compression_level, None);
    assert_eq!(parsed.batch_size, 2_048);
    assert_eq!(parsed.queue_batches, 3);

    let expanded_worker_count = parse_options_from(
        [
            "--index",
            "reference.bsbit",
            "--read1",
            "reads.fastq.gz",
            "--output",
            "output.bam",
            "--compression-threads",
            "65",
        ]
        .map(OsString::from),
    )
    .expect("compression workers are not arbitrarily capped at 64");
    assert_eq!(expanded_worker_count.compression_threads, 65);

    let error = parse_options_from(
        [
            "--index",
            "reference.bsbit",
            "--read1",
            "reads.fastq.gz",
            "--output",
            "output.bam",
            "--compression-threads",
            "2147483648",
        ]
        .map(OsString::from),
    )
    .expect_err("compression workers must fit the native API domain");
    assert!(error.to_string().contains("signed 32-bit worker domain"));
}

#[test]
fn canonical_short_options_and_total_thread_budget_parse() {
    let parsed = parse_options_from(
        [
            "-x",
            "reference.bsbit",
            "-1",
            "r1.fq",
            "-2",
            "r2.fq",
            "-o",
            "output.bam",
            "--total-threads",
            "12",
            "--compression-level",
            "default",
        ]
        .map(OsString::from),
    )
    .expect("short options and total budget parse");
    assert_eq!(parsed.total_thread_budget, Some(12));
    assert_eq!(parsed.compression_level, None);

    let conflict = parse_options_from(
        [
            "-x",
            "reference.bsbit",
            "-1",
            "r1.fq",
            "-2",
            "r2.fq",
            "-o",
            "output.bam",
            "--total-threads",
            "12",
            "-t",
            "8",
        ]
        .map(OsString::from),
    )
    .expect_err("total and explicit mapping threads conflict");
    assert!(conflict.to_string().contains("--total-threads conflicts"));
}

#[test]
fn total_thread_split_scales_the_output_budget() {
    assert_eq!(throughput_thread_split(1, false), (1, 0));
    assert_eq!(throughput_thread_split(2, false), (1, 1));
    assert_eq!(throughput_thread_split(10, false), (8, 2));
    assert_eq!(throughput_thread_split(14, false), (11, 3));
    assert_eq!(throughput_thread_split(14, true), (10, 4));
    assert_eq!(throughput_thread_split(64, false), (51, 13));
    assert_eq!(throughput_thread_split(64, true), (48, 16));
}

#[test]
fn read_layout_accepts_canonical_short_forms_and_rejects_duplicates() {
    let short_single = [
        "--index",
        "reference.bsbit",
        "-1",
        "single.fastq.gz",
        "--output",
        "single.bam",
    ]
    .map(OsString::from);
    let parsed = parse_options_from(short_single).expect("-1 selects single-end input");
    assert_eq!(parsed.layout(), ReadLayout::SingleEnd);
    assert_eq!(parsed.read1, std::path::PathBuf::from("single.fastq.gz"));
    assert!(parsed.read2.is_none());

    let short_pair = [
        "--index",
        "reference.bsbit",
        "-1",
        "r1.fastq.gz",
        "-2",
        "r2.fastq.gz",
        "--output",
        "paired.bam",
    ]
    .map(OsString::from);
    let parsed = parse_options_from(short_pair).expect("-1 and -2 select paired input");
    assert_eq!(parsed.layout(), ReadLayout::PairedEnd);
    assert_eq!(parsed.read1, std::path::PathBuf::from("r1.fastq.gz"));
    assert_eq!(
        parsed.read2.as_deref(),
        Some(std::path::Path::new("r2.fastq.gz"))
    );

    let duplicate_form = [
        "--index",
        "reference.bsbit",
        "--read1",
        "first.fastq.gz",
        "-1",
        "second.fastq.gz",
        "--output",
        "output.bam",
    ]
    .map(OsString::from);
    assert_eq!(
        parse_options_from(duplicate_form)
            .expect_err("short and long forms identify one option")
            .to_string(),
        "--read1 may be specified only once"
    );

    let read2_only = [
        "--index",
        "reference.bsbit",
        "-2",
        "r2.fastq.gz",
        "--output",
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
            "--output",
            "output.bam",
        ]
        .map(OsString::from)
        .to_vec()
    };

    let defaults = parse_options_from(required()).expect("default output contract");
    assert_eq!(defaults.output_contract, AlignmentAuxiliaryMode::Minimal);
    assert_eq!(defaults.library_profile, LibraryProfile::Directional);
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
        LibraryProfile::NonDirectional
    );
    assert_eq!(
        non_directional.output_contract,
        AlignmentAuxiliaryMode::Minimal
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

    let mut invalid_contract = required();
    invalid_contract.extend(["--output-contract", "unknown"].map(OsString::from));
    assert_eq!(
        parse_options_from(invalid_contract)
            .expect_err("unknown output contract is rejected")
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
        "--output",
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
    assert!(mapped_only.ends_with("-mapq0-hash-tie-v1"));
    assert!(complete.ends_with("-read-complete-hash-tie-v1"));
    assert_ne!(mapped_only, complete);
    assert_eq!(
        mapped_only,
        "sensitive-bounded-integrated-mapq0-hash-tie-v1"
    );
    assert_eq!(
        complete,
        "sensitive-bounded-integrated-read-complete-hash-tie-v1"
    );
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
            "--output",
            "output.bam",
        ]
        .map(OsString::from)
        .to_vec()
    };

    let default = parse_options_from(required()).expect("default options");
    assert_eq!(
        strategy_id(&default),
        "balanced-d5-adapter-recovery-read-complete-hash-tie-v1"
    );
    assert_eq!(default.search_mode, SearchMode::Default);
    assert_eq!(default.read_output, ReadOutputMode::Complete);

    let mut default_mapped_only = required();
    default_mapped_only.push(OsString::from("--mapped-only"));
    let default_mapped_only =
        parse_options_from(default_mapped_only).expect("mapped-only default output");
    assert_eq!(
        strategy_id(&default_mapped_only),
        "balanced-d5-adapter-recovery-mapq0-hash-tie-v1"
    );
    assert_eq!(default_mapped_only.read_output, ReadOutputMode::MappedOnly);

    let mut sensitive = required();
    sensitive.push(OsString::from("--sensitive"));
    let sensitive = parse_options_from(sensitive).expect("sensitive defaults");
    assert_eq!(
        strategy_id(&sensitive),
        sensitive_read_complete_strategy_id()
    );
    assert_eq!(sensitive.search_mode, SearchMode::Sensitive);
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

    let mut unknown = required();
    unknown.push(OsString::from("--unknown-option"));
    assert_eq!(
        parse_options_from(unknown)
            .expect_err("unknown options are rejected")
            .to_string(),
        "unknown option --unknown-option"
    );
}
