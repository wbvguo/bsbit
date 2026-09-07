//! Deterministic narrow-kernel benchmark used by the SIMD compatibility audit.

use std::hint::black_box;
use std::time::Instant;

use bsbit_align::verification::{NarrowEndpointDistances, narrow_banded_fixed_start_batch};
use bsbit_cpu::{BackendRequest, initialize};

// One batch fills exactly one AVX-512BW byte vector, two AVX2 vectors, or four
// 128-bit SSE vectors, so per-batch comparisons do not favor a partial chunk.
const CANDIDATES: usize = 64;
const QUERY_LENGTH: usize = 150;
const MAX_DISTANCE: usize = 3;

fn main() {
    let request = std::env::args()
        .nth(1)
        .map(|value| {
            value
                .parse::<BackendRequest>()
                .expect("backend must be auto, scalar, sse2, sse4.2, avx2, avx512, or neon")
        })
        .unwrap_or_default();
    // `cargo test --all-targets` also executes harness-free bench binaries in
    // the debug profile. Keep that structural smoke run short; real benchmark
    // and CI invocations pass an explicit iteration count.
    let default_iterations = if cfg!(debug_assertions) { 1 } else { 200_000 };
    let iterations = std::env::args().nth(2).map_or(default_iterations, |value| {
        value
            .parse::<u32>()
            .expect("iteration count must be an integer")
    });
    assert!(iterations > 0, "iteration count must be positive");
    let configuration = initialize(request).expect("requested backend must be supported");
    let pattern_length = QUERY_LENGTH + 2 * MAX_DISTANCE;
    let masks = [0b0001, 0b1010, 0b0100, 0b1010, 0];
    let query = (0..QUERY_LENGTH)
        .map(|position| u8::try_from((position * 3 + position / 11) % 5).unwrap())
        .collect::<Vec<_>>();
    let patterns = (0..CANDIDATES * pattern_length)
        .map(|position| u8::try_from((position * 5 + position / 7 + position / 31) % 5).unwrap())
        .collect::<Vec<_>>();
    let mut output = [NarrowEndpointDistances::EMPTY; CANDIDATES];

    for _ in 0..iterations.min(2_000) {
        narrow_banded_fixed_start_batch(
            black_box(&masks),
            black_box(&query),
            black_box(&patterns),
            MAX_DISTANCE,
            black_box(&mut output),
        )
        .unwrap();
    }

    let started = Instant::now();
    let mut checksum = 0_u64;
    for _ in 0..iterations {
        narrow_banded_fixed_start_batch(
            black_box(&masks),
            black_box(&query),
            black_box(&patterns),
            MAX_DISTANCE,
            black_box(&mut output),
        )
        .unwrap();
        checksum = checksum.wrapping_add(u64::from(output[0].in_budget_mask()));
    }
    let elapsed = started.elapsed();
    let ns_per_batch = elapsed.as_secs_f64() * 1_000_000_000.0 / f64::from(iterations);
    println!(
        "backend={}\titerations={iterations}\tns_per_batch={ns_per_batch:.3}\tchecksum={checksum}",
        configuration.backend()
    );
}
