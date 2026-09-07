//! Public backend-selection contract tests.

use bsbit_cpu::{Backend, BackendRequest, CpuFeatures};

#[test]
fn automatic_policy_orders_all_x86_backends_and_requires_popcnt_independently() {
    let cases = [
        (
            CpuFeatures::new_x86_64(true, true, true).with_avx512(true, true),
            Backend::Avx512BwPopcnt,
        ),
        (
            CpuFeatures::new_x86_64(true, true, true).with_avx512(true, false),
            Backend::Avx2Popcnt,
        ),
        (
            CpuFeatures::new_x86_64(true, true, true),
            Backend::Avx2Popcnt,
        ),
        (
            CpuFeatures::new_x86_64(true, false, true),
            Backend::Avx2Popcnt,
        ),
        (
            CpuFeatures::new_x86_64(false, true, true).with_sse42(true),
            Backend::Sse42Popcnt,
        ),
        (CpuFeatures::new_x86_64(false, true, true), Backend::Sse2),
        (CpuFeatures::new_x86_64(false, false, false), Backend::Sse2),
        (
            CpuFeatures::new_x86_64(false, false, false).with_sse2(false),
            Backend::Scalar,
        ),
        (CpuFeatures::new_x86_64(false, true, false), Backend::Sse2),
        (CpuFeatures::new_x86_64(false, false, true), Backend::Sse2),
        (CpuFeatures::new_x86_64(true, true, false), Backend::Sse2),
    ];

    for (features, expected) in cases {
        assert_eq!(BackendRequest::Auto.select(features), Ok(expected));
    }
}

#[test]
fn automatic_aarch64_policy_selects_neon() {
    assert_eq!(
        BackendRequest::Auto.select(CpuFeatures::new_aarch64()),
        Ok(Backend::Neon)
    );
    assert_eq!(
        BackendRequest::Neon.select(CpuFeatures::new_aarch64()),
        Ok(Backend::Neon)
    );
}

#[test]
fn automatic_policy_uses_scalar_when_no_specialized_backend_is_available() {
    assert_eq!(
        BackendRequest::Auto.select(CpuFeatures::new_aarch64_with_neon(false)),
        Ok(Backend::Scalar)
    );
    assert_eq!(
        BackendRequest::Auto.select(CpuFeatures::new_other()),
        Ok(Backend::Scalar)
    );
    for features in [
        CpuFeatures::new_x86_64(true, true, true),
        CpuFeatures::new_aarch64(),
        CpuFeatures::new_other(),
    ] {
        assert_eq!(BackendRequest::Scalar.select(features), Ok(Backend::Scalar));
    }
}

#[cfg(target_arch = "aarch64")]
#[test]
fn detected_aarch64_contract_keeps_features_independent() {
    let features = CpuFeatures::detect();
    assert_eq!(features.architecture(), bsbit_cpu::Architecture::Aarch64);
    assert!(!features.sse2());
    assert!(!features.sse41());
    assert!(!features.sse42());
    assert!(!features.popcnt());
    assert!(!features.avx2());
    assert!(!features.avx512f());
    assert!(!features.avx512bw());
    let expected = if features.neon() {
        Backend::Neon
    } else {
        Backend::Scalar
    };
    assert_eq!(BackendRequest::Auto.select(features), Ok(expected));
}

#[test]
fn forced_backend_rejects_the_wrong_architecture() {
    let arm_error = BackendRequest::Avx2Popcnt
        .select(CpuFeatures::new_aarch64())
        .expect_err("x86 backend must be unavailable on AArch64");
    assert_eq!(
        arm_error.to_string(),
        "SIMD backend avx2 requires x86_64; current architecture is aarch64"
    );

    let x86_error = BackendRequest::Neon
        .select(CpuFeatures::new_x86_64(false, false, false))
        .expect_err("NEON backend must be unavailable on x86-64");
    assert_eq!(
        x86_error.to_string(),
        "SIMD backend neon requires aarch64; current architecture is x86_64"
    );
}

#[test]
fn forced_sse2_reports_the_missing_feature() {
    assert_unavailable_cases(&[(
        BackendRequest::Sse2,
        CpuFeatures::new_x86_64(false, false, false).with_sse2(false),
        Backend::Sse2,
        "SIMD backend sse2 requires SSE2; CPU reports SSE2=0",
    )]);
}

#[test]
fn forced_later_x86_tiers_report_each_independently_missing_feature() {
    let cases = [
        (
            BackendRequest::Avx2Popcnt,
            CpuFeatures::new_x86_64(true, true, false),
            Backend::Avx2Popcnt,
            "SIMD backend avx2 requires AVX2+POPCNT; missing POPCNT",
        ),
        (
            BackendRequest::Avx2Popcnt,
            CpuFeatures::new_x86_64(false, true, true),
            Backend::Avx2Popcnt,
            "SIMD backend avx2 requires AVX2+POPCNT; missing AVX2",
        ),
        (
            BackendRequest::Sse42Popcnt,
            CpuFeatures::new_x86_64(false, true, false).with_sse42(true),
            Backend::Sse42Popcnt,
            "SIMD backend sse4.2 requires SSE4.2+POPCNT; missing POPCNT",
        ),
        (
            BackendRequest::Sse42Popcnt,
            CpuFeatures::new_x86_64(false, true, true),
            Backend::Sse42Popcnt,
            "SIMD backend sse4.2 requires SSE4.2+POPCNT; missing SSE4.2",
        ),
        (
            BackendRequest::Sse42Popcnt,
            CpuFeatures::new_x86_64(false, false, true).with_sse42(true),
            Backend::Sse42Popcnt,
            "SIMD backend sse4.2 requires SSE4.2+POPCNT; missing SSE4.1",
        ),
        (
            BackendRequest::Avx512BwPopcnt,
            CpuFeatures::new_x86_64(true, true, true).with_avx512(false, true),
            Backend::Avx512BwPopcnt,
            "SIMD backend avx512 requires AVX-512F+AVX-512BW+AVX2+POPCNT; missing AVX-512F",
        ),
        (
            BackendRequest::Avx512BwPopcnt,
            CpuFeatures::new_x86_64(true, true, true).with_avx512(true, false),
            Backend::Avx512BwPopcnt,
            "SIMD backend avx512 requires AVX-512F+AVX-512BW+AVX2+POPCNT; missing AVX-512BW",
        ),
        (
            BackendRequest::Avx512BwPopcnt,
            CpuFeatures::new_x86_64(false, true, true).with_avx512(true, true),
            Backend::Avx512BwPopcnt,
            "SIMD backend avx512 requires AVX-512F+AVX-512BW+AVX2+POPCNT; missing AVX2",
        ),
        (
            BackendRequest::Avx512BwPopcnt,
            CpuFeatures::new_x86_64(true, true, false).with_avx512(true, true),
            Backend::Avx512BwPopcnt,
            "SIMD backend avx512 requires AVX-512F+AVX-512BW+AVX2+POPCNT; missing POPCNT",
        ),
    ];

    assert_unavailable_cases(&cases);
}

fn assert_unavailable_cases(cases: &[(BackendRequest, CpuFeatures, Backend, &str)]) {
    for &(request, features, requested, message) in cases {
        let error = request
            .select(features)
            .expect_err("unsupported forced backend must be rejected");
        assert_eq!(error.requested(), requested);
        assert_eq!(error.features(), features);
        assert_eq!(error.to_string(), message);
    }
}

#[test]
fn parser_uses_stable_public_spellings() {
    assert_eq!("auto".parse(), Ok(BackendRequest::Auto));
    assert_eq!("scalar".parse(), Ok(BackendRequest::Scalar));
    assert_eq!("sse2".parse(), Ok(BackendRequest::Sse2));
    assert_eq!("sse4.2".parse(), Ok(BackendRequest::Sse42Popcnt));
    assert_eq!("avx2".parse(), Ok(BackendRequest::Avx2Popcnt));
    assert_eq!("avx512".parse(), Ok(BackendRequest::Avx512BwPopcnt));
    assert_eq!("neon".parse(), Ok(BackendRequest::Neon));
    for removed in ["sse", "sse3", "ssse3", "sse4.1", "avx", "native"] {
        assert!(removed.parse::<BackendRequest>().is_err());
    }
}

#[test]
fn instruction_set_names_include_every_independent_requirement() {
    assert_eq!(Backend::Scalar.instruction_set(), "portable");
    assert_eq!(Backend::Sse2.instruction_set(), "sse2");
    assert_eq!(Backend::Sse42Popcnt.instruction_set(), "sse4.2+popcnt");
    assert_eq!(Backend::Avx2Popcnt.instruction_set(), "avx2+popcnt");
    assert_eq!(
        Backend::Avx512BwPopcnt.instruction_set(),
        "avx512f+avx512bw+avx2+popcnt"
    );
    assert_eq!(Backend::Neon.instruction_set(), "neon");
}
