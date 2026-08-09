//! Pinned in-process libsais FFI adapters.
//!
//! Release builds use an in-process low-memory builder. Test-only helpers use
//! direct libsais64 to validate the current format implementation.

#![deny(unsafe_op_in_unsafe_fn)]
#![cfg(feature = "index-construction")]

#[cfg(test)]
use crate::build::SuffixArrayBuildError;
#[cfg(test)]
use crate::build::SuffixArrayBuilder;
#[cfg(test)]
use crate::storage::fm::SearchBase;
use core::ffi::c_int;
#[cfg(test)]
use core::ffi::c_longlong;
use std::ffi::c_uchar;

unsafe extern "C" {
    #[cfg(test)]
    fn libsais(
        text: *const c_uchar,
        suffix_array: *mut c_int,
        length: c_int,
        extra_space: c_int,
        frequencies: *mut c_int,
    ) -> c_int;
    pub(crate) fn libsais_omp(
        text: *const c_uchar,
        suffix_array: *mut c_int,
        length: c_int,
        extra_space: c_int,
        frequencies: *mut c_int,
        threads: c_int,
    ) -> c_int;
    #[cfg(all(test, feature = "index-construction"))]
    pub(crate) fn libsais_bwt_aux_omp(
        text: *const c_uchar,
        transformed: *mut c_uchar,
        temporary: *mut c_int,
        length: c_int,
        extra_space: c_int,
        frequencies: *mut c_int,
        sample_stride: c_int,
        sampled_rows: *mut c_int,
        threads: c_int,
    ) -> c_int;
    #[cfg(test)]
    fn libsais64(
        text: *const c_uchar,
        suffix_array: *mut c_longlong,
        length: c_longlong,
        extra_space: c_longlong,
        frequencies: *mut c_longlong,
    ) -> c_longlong;
    #[cfg(test)]
    fn libsais64_omp(
        text: *const c_uchar,
        suffix_array: *mut c_longlong,
        length: c_longlong,
        extra_space: c_longlong,
        frequencies: *mut c_longlong,
        threads: c_longlong,
    ) -> c_longlong;
    #[cfg(test)]
    pub(crate) fn libsais64_bwt_aux_omp(
        text: *const c_uchar,
        transformed: *mut c_uchar,
        temporary: *mut c_longlong,
        length: c_longlong,
        extra_space: c_longlong,
        frequencies: *mut c_longlong,
        sample_stride: c_longlong,
        sampled_rows: *mut c_longlong,
        threads: c_longlong,
    ) -> c_longlong;
}

/// Pinned libsais v2.10.4 constructor.
#[cfg(test)]
#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct LibsaisSuffixArrayBuilder;

#[cfg(test)]
impl SuffixArrayBuilder for LibsaisSuffixArrayBuilder {
    fn backend_name(&self) -> &'static str {
        "libsais-2.10.4"
    }

    fn build_suffix_array(&self, text: &[SearchBase]) -> Result<Vec<u32>, SuffixArrayBuildError> {
        build_libsais_suffix_array(text, self.backend_name(), LibsaisExecution::Serial)
    }
}

/// Pinned libsais v2.10.4 OpenMP constructor with an explicit thread count.
#[cfg(test)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct LibsaisOpenMpSuffixArrayBuilder {
    threads: c_int,
}

#[cfg(test)]
impl LibsaisOpenMpSuffixArrayBuilder {
    /// Creates an OpenMP builder with a positive C-compatible thread count.
    ///
    /// # Errors
    ///
    /// Rejects zero or a count outside the signed 32-bit libsais API domain.
    pub(crate) fn new(threads: u32) -> Result<Self, SuffixArrayBuildError> {
        let threads = c_int::try_from(threads).map_err(|_| {
            backend_error(
                "libsais-2.10.4-openmp",
                "thread count exceeds the signed 32-bit API domain",
            )
        })?;
        if threads == 0 {
            return Err(backend_error(
                "libsais-2.10.4-openmp",
                "thread count must be nonzero",
            ));
        }
        Ok(Self { threads })
    }

    /// Returns the configured OpenMP thread count.
    #[must_use]
    pub(crate) fn threads(self) -> u32 {
        self.threads.cast_unsigned()
    }
}

#[cfg(test)]
impl SuffixArrayBuilder for LibsaisOpenMpSuffixArrayBuilder {
    fn backend_name(&self) -> &'static str {
        "libsais-2.10.4-openmp"
    }

    fn build_suffix_array(&self, text: &[SearchBase]) -> Result<Vec<u32>, SuffixArrayBuildError> {
        build_libsais_suffix_array(
            text,
            self.backend_name(),
            LibsaisExecution::OpenMp(self.threads),
        )
    }
}

#[cfg(test)]
#[derive(Clone, Copy)]
enum LibsaisExecution {
    Serial,
    OpenMp(c_int),
}

#[cfg(test)]
fn build_libsais_suffix_array(
    text: &[SearchBase],
    backend: &'static str,
    execution: LibsaisExecution,
) -> Result<Vec<u32>, SuffixArrayBuildError> {
    let terminal = u32::try_from(text.len()).map_err(|_| SuffixArrayBuildError::TextExceedsU32)?;
    let input = encoded_text(text)?;
    if input.is_empty() {
        return Ok(vec![0]);
    }
    let mut suffixes = reserved_vec(text.len() + 1, "libsais suffix array")?;
    suffixes.push(terminal);
    if let Ok(length) = c_int::try_from(input.len()) {
        let mut native = reserved_vec(input.len(), "libsais i32 output")?;
        native.resize(input.len(), 0_i32);
        // SAFETY: `input` and `native` are live, non-overlapping allocations
        // of exactly `length` elements. The optional frequency pointer is
        // null, and the pinned API writes only within the supplied SA.
        let status = unsafe {
            match execution {
                LibsaisExecution::Serial => libsais(
                    input.as_ptr(),
                    native.as_mut_ptr(),
                    length,
                    0,
                    core::ptr::null_mut(),
                ),
                LibsaisExecution::OpenMp(threads) => libsais_omp(
                    input.as_ptr(),
                    native.as_mut_ptr(),
                    length,
                    0,
                    core::ptr::null_mut(),
                    threads,
                ),
            }
        };
        if status != 0 {
            return Err(backend_error(
                backend,
                format!("libsais returned status {status}"),
            ));
        }
        for offset in native {
            suffixes.push(
                u32::try_from(offset)
                    .map_err(|_| backend_error(backend, "negative or oversized i32 suffix"))?,
            );
        }
    } else {
        let length =
            c_longlong::try_from(input.len()).map_err(|_| SuffixArrayBuildError::TextExceedsU32)?;
        let mut native = reserved_vec(input.len(), "libsais i64 output")?;
        native.resize(input.len(), 0_i64);
        // SAFETY: the same allocation and length conditions as the 32-bit
        // call hold; this branch uses the pinned 64-bit API for long runs.
        let status = unsafe {
            match execution {
                LibsaisExecution::Serial => libsais64(
                    input.as_ptr(),
                    native.as_mut_ptr(),
                    length,
                    0,
                    core::ptr::null_mut(),
                ),
                LibsaisExecution::OpenMp(threads) => libsais64_omp(
                    input.as_ptr(),
                    native.as_mut_ptr(),
                    length,
                    0,
                    core::ptr::null_mut(),
                    c_longlong::from(threads),
                ),
            }
        };
        if status != 0 {
            return Err(backend_error(
                backend,
                format!("libsais64 returned status {status}"),
            ));
        }
        for offset in native {
            suffixes.push(
                u32::try_from(offset)
                    .map_err(|_| backend_error(backend, "negative or oversized i64 suffix"))?,
            );
        }
    }
    Ok(suffixes)
}

#[cfg(test)]
fn encoded_text(text: &[SearchBase]) -> Result<Vec<u8>, SuffixArrayBuildError> {
    let mut encoded = reserved_vec(text.len(), "encoded suffix-array text")?;
    encoded.extend(text.iter().map(|base| base.as_ascii()));
    Ok(encoded)
}

#[cfg(test)]
fn backend_error(backend: &'static str, message: impl Into<String>) -> SuffixArrayBuildError {
    SuffixArrayBuildError::Backend {
        backend,
        message: message.into(),
    }
}

#[cfg(test)]
pub(crate) fn reserved_vec<T>(
    elements: usize,
    component: &'static str,
) -> Result<Vec<T>, SuffixArrayBuildError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(elements)
        .map_err(|_| SuffixArrayBuildError::AllocationFailed {
            component,
            elements,
        })?;
    Ok(values)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::build::{PrefixDoublingSuffixArrayBuilder, validate_suffix_array};

    fn canonical(bytes: &[u8]) -> Vec<SearchBase> {
        bytes
            .iter()
            .map(|byte| match byte {
                b'A' => SearchBase::A,
                b'C' => SearchBase::C,
                b'G' => SearchBase::G,
                b'T' => SearchBase::T,
                _ => panic!("canonical fixture"),
            })
            .collect()
    }

    #[test]
    fn libsais_matches_rust_suffix_array() {
        let rust = PrefixDoublingSuffixArrayBuilder;
        let libsais = LibsaisSuffixArrayBuilder;
        for bytes in [
            b"".as_slice(),
            b"A",
            b"ACGTACGT",
            b"AAAAAAAAAAAAAAAA",
            b"GATTACACCTG",
        ] {
            let text = canonical(bytes);
            let expected = rust.build_suffix_array(&text).expect("Rust SA");
            let observed = libsais.build_suffix_array(&text).expect("libsais SA");
            validate_suffix_array(&text, &observed).expect("libsais validates");
            assert_eq!(observed, expected);
        }
    }

    #[test]
    fn libsais_openmp_matches_serial_and_rejects_invalid_threads() {
        assert!(LibsaisOpenMpSuffixArrayBuilder::new(0).is_err());
        let text = canonical(b"ACGTACGTGATTACACCCCCAAAAATTTTGGGG");
        let expected = LibsaisSuffixArrayBuilder
            .build_suffix_array(&text)
            .expect("serial libsais SA");
        for threads in [1, 2, 4] {
            let builder = LibsaisOpenMpSuffixArrayBuilder::new(threads).expect("thread count");
            assert_eq!(builder.threads(), threads);
            let observed = builder
                .build_suffix_array(&text)
                .expect("OpenMP libsais SA");
            validate_suffix_array(&text, &observed).expect("OpenMP libsais validates");
            assert_eq!(observed, expected);
        }
    }
}
