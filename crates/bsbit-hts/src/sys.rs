//! Audited opaque ownership wrappers over the project `HTSlib` C shim.
//!
//! All project unsafe code for `HTSlib` is confined to this module. Its safe
//! wrappers own opaque handles and never expose a raw pointer or `HTSlib`
//! structure outside the module.

#![deny(unsafe_op_in_unsafe_fn)]

use core::ffi::{c_char, c_int};
use core::fmt;
use core::marker::PhantomData;
use core::mem::MaybeUninit;
use core::ptr::NonNull;
use std::ffi::{CStr, c_void};
use std::rc::Rc;

const ERROR_CAPACITY: usize = 512;

#[repr(C)]
struct CReader {
    _private: [u8; 0],
    _alignment: [c_void; 0],
}

#[repr(C)]
struct CIndexedReader {
    _private: [u8; 0],
    _alignment: [c_void; 0],
}

#[repr(C)]
struct CIndexedFastaReader {
    _private: [u8; 0],
    _alignment: [c_void; 0],
}

#[repr(C)]
struct CBgzfWriter {
    _private: [u8; 0],
    _alignment: [c_void; 0],
}

#[repr(C)]
#[derive(Clone, Copy)]
struct CBamRecordView {
    reference_id: i32,
    position: i64,
    mapping_quality: u8,
    flag: u16,
    mate_reference_id: i32,
    mate_position: i64,
    template_length: i64,
    query_name: *const c_char,
    query_name_length: usize,
    cigar: *const u32,
    cigar_count: usize,
    sequence: *const u8,
    packed_sequence_length: usize,
    sequence_length: usize,
    quality: *const u8,
    auxiliary: *const u8,
    auxiliary_length: usize,
}

impl CBamRecordView {
    const fn empty() -> Self {
        Self {
            reference_id: -1,
            position: -1,
            mapping_quality: 0,
            flag: 0,
            mate_reference_id: -1,
            mate_position: -1,
            template_length: 0,
            query_name: core::ptr::null(),
            query_name_length: 0,
            cigar: core::ptr::null(),
            cigar_count: 0,
            sequence: core::ptr::null(),
            packed_sequence_length: 0,
            sequence_length: 0,
            quality: core::ptr::null(),
            auxiliary: core::ptr::null(),
            auxiliary_length: 0,
        }
    }
}

#[repr(C)]
struct CWriter {
    _private: [u8; 0],
    _alignment: [c_void; 0],
}

unsafe extern "C" {
    #[cfg(test)]
    fn bsbit_hts_shim_abi_version() -> u32;
    #[cfg(test)]
    fn bsbit_hts_runtime_version() -> *const c_char;
    #[cfg(test)]
    fn bsbit_hts_health_check() -> c_int;
    #[cfg(test)]
    fn bsbit_hts_tabix_index_build(
        path: *const c_char,
        index_path: *const c_char,
        preset: c_int,
        threads: u32,
        out_system_errno: *mut c_int,
        error: *mut c_char,
        error_capacity: usize,
    ) -> c_int;
    fn bsbit_hts_bam_index_build(
        path: *const c_char,
        index_path: *const c_char,
        threads: u32,
        out_system_errno: *mut c_int,
        error: *mut c_char,
        error_capacity: usize,
    ) -> c_int;
    fn bsbit_hts_reader_open(
        path: *const c_char,
        out_reader: *mut *mut CReader,
        out_system_errno: *mut c_int,
        error: *mut c_char,
        error_capacity: usize,
    ) -> c_int;
    fn bsbit_hts_reader_compression(
        reader: *const CReader,
        out_compression: *mut c_int,
        out_system_errno: *mut c_int,
        error: *mut c_char,
        error_capacity: usize,
    ) -> c_int;
    fn bsbit_hts_reader_read(
        reader: *mut CReader,
        buffer: *mut u8,
        capacity: usize,
        out_count: *mut usize,
        out_system_errno: *mut c_int,
        error: *mut c_char,
        error_capacity: usize,
    ) -> c_int;
    fn bsbit_hts_reader_close(
        reader: *mut CReader,
        out_system_errno: *mut c_int,
        error: *mut c_char,
        error_capacity: usize,
    ) -> c_int;
    fn bsbit_hts_reader_destroy(reader: *mut CReader);
    fn bsbit_hts_indexed_reader_open(
        path: *const c_char,
        out_reader: *mut *mut CIndexedReader,
        out_system_errno: *mut c_int,
        error: *mut c_char,
        error_capacity: usize,
    ) -> c_int;
    fn bsbit_hts_indexed_reader_header_text(
        reader: *const CIndexedReader,
        out_text: *mut *const c_char,
        out_length: *mut usize,
        out_system_errno: *mut c_int,
        error: *mut c_char,
        error_capacity: usize,
    ) -> c_int;
    fn bsbit_hts_indexed_reader_reference_count(
        reader: *const CIndexedReader,
        out_count: *mut i32,
        out_system_errno: *mut c_int,
        error: *mut c_char,
        error_capacity: usize,
    ) -> c_int;
    fn bsbit_hts_indexed_reader_reference(
        reader: *const CIndexedReader,
        reference_id: i32,
        out_name: *mut *const c_char,
        out_name_length: *mut usize,
        out_length: *mut i64,
        out_system_errno: *mut c_int,
        error: *mut c_char,
        error_capacity: usize,
    ) -> c_int;
    fn bsbit_hts_indexed_reader_query(
        reader: *mut CIndexedReader,
        reference_id: i32,
        start: i64,
        end: i64,
        out_system_errno: *mut c_int,
        error: *mut c_char,
        error_capacity: usize,
    ) -> c_int;
    fn bsbit_hts_indexed_reader_next(
        reader: *mut CIndexedReader,
        out_record: *mut CBamRecordView,
        out_has_record: *mut c_int,
        out_system_errno: *mut c_int,
        error: *mut c_char,
        error_capacity: usize,
    ) -> c_int;
    fn bsbit_hts_indexed_reader_close(
        reader: *mut CIndexedReader,
        out_system_errno: *mut c_int,
        error: *mut c_char,
        error_capacity: usize,
    ) -> c_int;
    fn bsbit_hts_indexed_reader_destroy(reader: *mut CIndexedReader);
    fn bsbit_hts_indexed_fasta_reader_open(
        path: *const c_char,
        out_reader: *mut *mut CIndexedFastaReader,
        out_system_errno: *mut c_int,
        error: *mut c_char,
        error_capacity: usize,
    ) -> c_int;
    fn bsbit_hts_indexed_fasta_reader_reference_count(
        reader: *const CIndexedFastaReader,
        out_count: *mut i32,
        out_system_errno: *mut c_int,
        error: *mut c_char,
        error_capacity: usize,
    ) -> c_int;
    fn bsbit_hts_indexed_fasta_reader_reference(
        reader: *const CIndexedFastaReader,
        reference_id: i32,
        out_name: *mut *const c_char,
        out_name_length: *mut usize,
        out_length: *mut i64,
        out_system_errno: *mut c_int,
        error: *mut c_char,
        error_capacity: usize,
    ) -> c_int;
    fn bsbit_hts_indexed_fasta_reader_fetch(
        reader: *mut CIndexedFastaReader,
        reference_id: i32,
        start: i64,
        end: i64,
        out_sequence: *mut *const c_char,
        out_length: *mut usize,
        out_system_errno: *mut c_int,
        error: *mut c_char,
        error_capacity: usize,
    ) -> c_int;
    fn bsbit_hts_indexed_fasta_reader_close(
        reader: *mut CIndexedFastaReader,
        out_system_errno: *mut c_int,
        error: *mut c_char,
        error_capacity: usize,
    ) -> c_int;
    fn bsbit_hts_indexed_fasta_reader_destroy(reader: *mut CIndexedFastaReader);
    fn bsbit_hts_bgzf_writer_open(
        path: *const c_char,
        compression_threads: u32,
        out_writer: *mut *mut CBgzfWriter,
        out_system_errno: *mut c_int,
        error: *mut c_char,
        error_capacity: usize,
    ) -> c_int;
    fn bsbit_hts_bgzf_writer_write(
        writer: *mut CBgzfWriter,
        data: *const u8,
        length: usize,
        out_count: *mut usize,
        out_system_errno: *mut c_int,
        error: *mut c_char,
        error_capacity: usize,
    ) -> c_int;
    fn bsbit_hts_bgzf_writer_flush(
        writer: *mut CBgzfWriter,
        out_system_errno: *mut c_int,
        error: *mut c_char,
        error_capacity: usize,
    ) -> c_int;
    fn bsbit_hts_bgzf_writer_finish(
        writer: *mut CBgzfWriter,
        out_system_errno: *mut c_int,
        error: *mut c_char,
        error_capacity: usize,
    ) -> c_int;
    fn bsbit_hts_bgzf_writer_destroy(writer: *mut CBgzfWriter);
    fn bsbit_hts_writer_open_bam_threads(
        path: *const c_char,
        header_text: *const c_char,
        header_length: usize,
        compression_threads: u32,
        out_writer: *mut *mut CWriter,
        out_system_errno: *mut c_int,
        error: *mut c_char,
        error_capacity: usize,
    ) -> c_int;
    fn bsbit_hts_writer_open_bam_threads_level(
        path: *const c_char,
        header_text: *const c_char,
        header_length: usize,
        compression_threads: u32,
        compression_level: c_int,
        out_writer: *mut *mut CWriter,
        out_system_errno: *mut c_int,
        error: *mut c_char,
        error_capacity: usize,
    ) -> c_int;
    fn bsbit_hts_writer_write_record(
        writer: *mut CWriter,
        record_text: *const c_char,
        record_length: usize,
        out_system_errno: *mut c_int,
        error: *mut c_char,
        error_capacity: usize,
    ) -> c_int;
    fn bsbit_hts_writer_write_bam_fields(
        writer: *mut CWriter,
        query_name: *const c_char,
        query_name_length: usize,
        flag: u16,
        reference_id: i32,
        position: i64,
        mapping_quality: u8,
        cigar: *const u32,
        cigar_count: usize,
        mate_reference_id: i32,
        mate_position: i64,
        template_length: i64,
        sequence: *const c_char,
        sequence_length: usize,
        quality: *const u8,
        has_mapping: c_int,
        literal_nm: u32,
        has_md: c_int,
        md: *const c_char,
        md_length: usize,
        has_xg: c_int,
        xg: *const c_char,
        has_bismark: c_int,
        bismark_xm: *const c_char,
        bismark_xm_length: usize,
        bismark_xr: *const c_char,
        out_system_errno: *mut c_int,
        error: *mut c_char,
        error_capacity: usize,
    ) -> c_int;
    fn bsbit_hts_writer_finish(
        writer: *mut CWriter,
        out_system_errno: *mut c_int,
        error: *mut c_char,
        error_capacity: usize,
    ) -> c_int;
    fn bsbit_hts_writer_destroy(writer: *mut CWriter);
}

/// One stable project shim status without exposing C integer sentinels.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum NativeStatus {
    /// A required pointer, string, length, or state argument was invalid.
    InvalidArgument,
    /// Native allocation failed.
    AllocationFailed,
    /// A source or destination could not be opened.
    OpenFailed,
    /// A SAM/BAM header was rejected or could not be written.
    HeaderFailed,
    /// A SAM record was rejected.
    RecordFailed,
    /// Closing or finalizing a native stream failed.
    CloseFailed,
    /// The native handle is already terminal or closed.
    InvalidState,
    /// Decoding input failed.
    ReadFailed,
    /// Encoding output failed.
    WriteFailed,
    /// The shim returned a status outside its accepted ABI.
    Unknown(i32),
}

impl NativeStatus {
    const fn from_raw(raw: c_int) -> Option<Self> {
        match raw {
            0 => None,
            1 => Some(Self::InvalidArgument),
            2 => Some(Self::AllocationFailed),
            3 => Some(Self::OpenFailed),
            4 => Some(Self::HeaderFailed),
            5 => Some(Self::RecordFailed),
            6 => Some(Self::CloseFailed),
            7 => Some(Self::InvalidState),
            8 => Some(Self::ReadFailed),
            9 => Some(Self::WriteFailed),
            value => Some(Self::Unknown(value)),
        }
    }
}

/// A fully copied native failure; it borrows no C storage.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct NativeError {
    status: NativeStatus,
    system_errno: i32,
    message: String,
}

impl NativeError {
    /// Returns the typed shim status.
    #[must_use]
    pub const fn status(&self) -> NativeStatus {
        self.status
    }

    /// Returns the captured `errno`, or zero when the operation had none.
    #[must_use]
    pub const fn system_errno(&self) -> i32 {
        self.system_errno
    }

    /// Returns the copied ASCII native diagnostic.
    #[must_use]
    pub fn message(&self) -> &str {
        &self.message
    }

    fn invalid_state(message: &str) -> Self {
        Self {
            status: NativeStatus::InvalidState,
            system_errno: 0,
            message: String::from(message),
        }
    }
}

impl fmt::Display for NativeError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "native HTS {:?} (errno {}): {}",
            self.status, self.system_errno, self.message
        )
    }
}

impl std::error::Error for NativeError {}

/// Compression classification returned from content-aware `HTSlib` detection.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum NativeCompression {
    /// Uncompressed bytes.
    Plain,
    /// Generic RFC 1952 gzip.
    Gzip,
    /// Block gzip.
    Bgzf,
}

impl NativeCompression {
    fn from_raw(raw: c_int) -> Result<Self, NativeError> {
        match raw {
            0 => Ok(Self::Plain),
            1 => Ok(Self::Gzip),
            2 => Ok(Self::Bgzf),
            _ => Err(NativeError::invalid_state(
                "shim returned an unknown compression class",
            )),
        }
    }
}

/// Returns the project C shim ABI version.
#[must_use]
#[cfg(test)]
pub fn shim_abi_version() -> u32 {
    // SAFETY: the function has no arguments and returns a value type.
    unsafe { bsbit_hts_shim_abi_version() }
}

/// Returns a copied `HTSlib` runtime version.
///
/// # Errors
///
/// Returns an invalid-state error if the native runtime returns a null string.
#[cfg(test)]
pub fn runtime_version() -> Result<String, NativeError> {
    // SAFETY: the runtime owns a process-lifetime NUL-terminated version string.
    let pointer = unsafe { bsbit_hts_runtime_version() };
    if pointer.is_null() {
        return Err(NativeError::invalid_state(
            "HTSlib returned a null runtime version",
        ));
    }
    // SAFETY: non-null was checked and HTSlib documents this as a C string.
    let version = unsafe { CStr::from_ptr(pointer) };
    Ok(version.to_string_lossy().into_owned())
}

/// Verifies the exact runtime expected by the shim.
///
/// # Errors
///
/// Returns the typed native status if the runtime health check fails.
#[cfg(test)]
pub fn health_check() -> Result<(), NativeError> {
    let call = NativeCall::new();
    // SAFETY: the function has no arguments and returns a status value.
    let status = unsafe { bsbit_hts_health_check() };
    call.finish(status)
}

#[cfg(test)]
pub fn build_tabix_index(
    path: &CStr,
    index_path: &CStr,
    preset: c_int,
    threads: u32,
) -> Result<(), NativeError> {
    let mut call = NativeCall::new();
    // SAFETY: both paths are NUL-terminated, all output pointers are writable,
    // and the native function retains none of them.
    let status = unsafe {
        bsbit_hts_tabix_index_build(
            path.as_ptr(),
            index_path.as_ptr(),
            preset,
            threads,
            &raw mut call.system_errno,
            call.error.as_mut_ptr(),
            call.error.len(),
        )
    };
    call.finish(status)
}

/// Builds a BAI at an explicit path for one coordinate-sorted BAM.
pub fn build_bam_index(path: &CStr, index_path: &CStr, threads: u32) -> Result<(), NativeError> {
    let mut call = NativeCall::new();
    // SAFETY: both paths are NUL-terminated, all output pointers are writable,
    // and the native function retains none of them.
    let status = unsafe {
        bsbit_hts_bam_index_build(
            path.as_ptr(),
            index_path.as_ptr(),
            threads,
            &raw mut call.system_errno,
            call.error.as_mut_ptr(),
            call.error.len(),
        )
    };
    call.finish(status)
}

/// One owned, thread-confined native decode handle.
pub struct NativeReader {
    handle: NonNull<CReader>,
    closed: bool,
    _thread_confined: PhantomData<Rc<()>>,
}

impl NativeReader {
    /// Opens one content-detected plain/gzip/BGZF path.
    ///
    /// # Errors
    ///
    /// Returns a copied native error and owns no handle on failure.
    pub fn open(path: &CStr) -> Result<Self, NativeError> {
        let mut handle = core::ptr::null_mut();
        let mut call = NativeCall::new();
        // SAFETY: all output pointers are valid for the call and `path` is NUL-terminated.
        let status = unsafe {
            bsbit_hts_reader_open(
                path.as_ptr(),
                &raw mut handle,
                &raw mut call.system_errno,
                call.error.as_mut_ptr(),
                call.error.len(),
            )
        };
        call.finish(status)?;
        let handle = NonNull::new(handle)
            .ok_or_else(|| NativeError::invalid_state("reader open returned a null handle"))?;
        Ok(Self {
            handle,
            closed: false,
            _thread_confined: PhantomData,
        })
    }

    /// Returns the content-derived compression class.
    ///
    /// # Errors
    ///
    /// Returns a copied native error for a closed or invalid handle.
    pub fn compression(&self) -> Result<NativeCompression, NativeError> {
        let mut compression = -1;
        let mut call = NativeCall::new();
        // SAFETY: `handle` remains owned and live until Drop.
        let status = unsafe {
            bsbit_hts_reader_compression(
                self.handle.as_ptr(),
                &raw mut compression,
                &raw mut call.system_errno,
                call.error.as_mut_ptr(),
                call.error.len(),
            )
        };
        call.finish(status)?;
        NativeCompression::from_raw(compression)
    }

    /// Decodes at most `buffer.len()` bytes.
    ///
    /// # Errors
    ///
    /// Returns a copied terminal native error without publishing bytes from the
    /// failing native call.
    pub fn read(&mut self, buffer: &mut [u8]) -> Result<usize, NativeError> {
        let mut count = 0;
        let mut call = NativeCall::new();
        let pointer = if buffer.is_empty() {
            core::ptr::null_mut()
        } else {
            buffer.as_mut_ptr()
        };
        // SAFETY: the buffer is writable for its reported length and the handle is exclusive.
        let status = unsafe {
            bsbit_hts_reader_read(
                self.handle.as_ptr(),
                pointer,
                buffer.len(),
                &raw mut count,
                &raw mut call.system_errno,
                call.error.as_mut_ptr(),
                call.error.len(),
            )
        };
        call.finish(status)?;
        if count > buffer.len() {
            return Err(NativeError::invalid_state(
                "reader returned more bytes than the supplied buffer",
            ));
        }
        Ok(count)
    }

    /// Closes the native stream exactly once.
    ///
    /// # Errors
    ///
    /// Returns a copied close error. The handle remains closed even on error.
    pub fn close(&mut self) -> Result<(), NativeError> {
        if self.closed {
            return Err(NativeError::invalid_state("reader is already closed"));
        }
        let mut call = NativeCall::new();
        // SAFETY: the handle is live and this method marks it closed exactly once.
        let status = unsafe {
            bsbit_hts_reader_close(
                self.handle.as_ptr(),
                &raw mut call.system_errno,
                call.error.as_mut_ptr(),
                call.error.len(),
            )
        };
        self.closed = true;
        call.finish(status)
    }
}

impl Drop for NativeReader {
    fn drop(&mut self) {
        if !self.closed {
            // SAFETY: Drop owns the live handle; null diagnostics are accepted by the shim.
            unsafe {
                let _ = bsbit_hts_reader_close(
                    self.handle.as_ptr(),
                    core::ptr::null_mut(),
                    core::ptr::null_mut(),
                    0,
                );
            }
        }
        // SAFETY: this is the sole owner and destroy is valid after close or close failure.
        unsafe { bsbit_hts_reader_destroy(self.handle.as_ptr()) };
    }
}

/// One owned, thread-confined native BGZF encode handle.
pub struct NativeBgzfWriter {
    handle: NonNull<CBgzfWriter>,
    finished: bool,
    _thread_confined: PhantomData<Rc<()>>,
}

impl NativeBgzfWriter {
    /// Opens one BGZF output path with optional private compression workers.
    ///
    /// `compression_threads == 0` selects synchronous compression. Values
    /// above 64 are rejected by the audited shim.
    ///
    /// # Errors
    ///
    /// Returns a copied native error and owns no handle on failure.
    pub fn open(path: &CStr, compression_threads: u32) -> Result<Self, NativeError> {
        let mut handle = core::ptr::null_mut();
        let mut call = NativeCall::new();
        // SAFETY: path is NUL-terminated and all output pointers are writable.
        let status = unsafe {
            bsbit_hts_bgzf_writer_open(
                path.as_ptr(),
                compression_threads,
                &raw mut handle,
                &raw mut call.system_errno,
                call.error.as_mut_ptr(),
                call.error.len(),
            )
        };
        call.finish(status)?;
        let handle = NonNull::new(handle)
            .ok_or_else(|| NativeError::invalid_state("BGZF writer open returned a null handle"))?;
        Ok(Self {
            handle,
            finished: false,
            _thread_confined: PhantomData,
        })
    }

    /// Encodes at most `data.len()` bytes and returns the accepted byte count.
    ///
    /// # Errors
    ///
    /// Returns a copied terminal native write error without claiming bytes
    /// from the failing native call.
    pub fn write(&mut self, data: &[u8]) -> Result<usize, NativeError> {
        let mut count = 0;
        let mut call = NativeCall::new();
        let pointer = if data.is_empty() {
            core::ptr::null()
        } else {
            data.as_ptr()
        };
        // SAFETY: data is readable for its exact length and the handle is exclusive.
        let status = unsafe {
            bsbit_hts_bgzf_writer_write(
                self.handle.as_ptr(),
                pointer,
                data.len(),
                &raw mut count,
                &raw mut call.system_errno,
                call.error.as_mut_ptr(),
                call.error.len(),
            )
        };
        call.finish(status)?;
        if count > data.len() {
            return Err(NativeError::invalid_state(
                "BGZF writer accepted more bytes than supplied",
            ));
        }
        Ok(count)
    }

    /// Flushes complete encoded blocks without finalizing the BGZF stream.
    ///
    /// # Errors
    ///
    /// Returns a copied terminal native write error.
    pub fn flush(&mut self) -> Result<(), NativeError> {
        let mut call = NativeCall::new();
        // SAFETY: the live handle is exclusively owned.
        let status = unsafe {
            bsbit_hts_bgzf_writer_flush(
                self.handle.as_ptr(),
                &raw mut call.system_errno,
                call.error.as_mut_ptr(),
                call.error.len(),
            )
        };
        call.finish(status)
    }

    /// Writes the canonical BGZF EOF block and finalizes exactly once.
    ///
    /// # Errors
    ///
    /// Returns a copied close or poisoned-state error. The stream remains
    /// terminal even when finalization fails.
    pub fn finish(&mut self) -> Result<(), NativeError> {
        if self.finished {
            return Err(NativeError::invalid_state(
                "BGZF writer is already finished",
            ));
        }
        let mut call = NativeCall::new();
        // SAFETY: the live handle is exclusively owned and finalization occurs once.
        let status = unsafe {
            bsbit_hts_bgzf_writer_finish(
                self.handle.as_ptr(),
                &raw mut call.system_errno,
                call.error.as_mut_ptr(),
                call.error.len(),
            )
        };
        self.finished = true;
        call.finish(status)
    }
}

impl Drop for NativeBgzfWriter {
    fn drop(&mut self) {
        if !self.finished {
            // SAFETY: Drop owns the live handle; null diagnostics are accepted.
            unsafe {
                let _ = bsbit_hts_bgzf_writer_finish(
                    self.handle.as_ptr(),
                    core::ptr::null_mut(),
                    core::ptr::null_mut(),
                    0,
                );
            }
        }
        // SAFETY: this is the sole owner and destroy accepts a finalized handle.
        unsafe { bsbit_hts_bgzf_writer_destroy(self.handle.as_ptr()) };
    }
}

/// One copied BAM reference dictionary entry.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct NativeBamReference {
    pub name: Vec<u8>,
    pub length: i64,
}

/// One BAM record borrowed from the indexed reader's reusable native storage.
pub struct NativeIndexedBamRecordView<'a> {
    pub reference_id: i32,
    pub position: i64,
    pub mapping_quality: u8,
    pub flag: u16,
    pub mate_reference_id: i32,
    pub mate_position: i64,
    pub template_length: i64,
    pub query_name: &'a [u8],
    pub cigar: &'a [u32],
    pub packed_sequence: &'a [u8],
    pub sequence_length: usize,
    pub quality: &'a [u8],
    pub auxiliary: &'a [u8],
}

/// One owned, thread-confined indexed BAM decode handle.
pub struct NativeIndexedBamReader {
    handle: NonNull<CIndexedReader>,
    closed: bool,
    _thread_confined: PhantomData<Rc<()>>,
}

impl NativeIndexedBamReader {
    /// Opens a BAM together with its adjacent BAI or CSI index.
    pub fn open(path: &CStr) -> Result<Self, NativeError> {
        let mut handle = core::ptr::null_mut();
        let mut call = NativeCall::new();
        // SAFETY: all output pointers are valid for the call and `path` is NUL-terminated.
        let status = unsafe {
            bsbit_hts_indexed_reader_open(
                path.as_ptr(),
                &raw mut handle,
                &raw mut call.system_errno,
                call.error.as_mut_ptr(),
                call.error.len(),
            )
        };
        call.finish(status)?;
        let handle = NonNull::new(handle)
            .ok_or_else(|| NativeError::invalid_state("indexed BAM open returned a null handle"))?;
        Ok(Self {
            handle,
            closed: false,
            _thread_confined: PhantomData,
        })
    }

    /// Copies the SAM header text owned by the native handle.
    pub fn header_text(&self) -> Result<Vec<u8>, NativeError> {
        let mut pointer = core::ptr::null();
        let mut length = 0_usize;
        let mut call = NativeCall::new();
        // SAFETY: the live handle owns the returned storage through this call.
        let status = unsafe {
            bsbit_hts_indexed_reader_header_text(
                self.handle.as_ptr(),
                &raw mut pointer,
                &raw mut length,
                &raw mut call.system_errno,
                call.error.as_mut_ptr(),
                call.error.len(),
            )
        };
        call.finish(status)?;
        // SAFETY: the shim guarantees `pointer` is readable for `length` bytes
        // until the next mutable operation; this method copies it immediately.
        unsafe { copy_native_slice(pointer.cast::<u8>(), length, "BAM header") }
    }

    /// Copies the complete BAM reference dictionary.
    pub fn references(&self) -> Result<Vec<NativeBamReference>, NativeError> {
        let mut count = 0_i32;
        let mut call = NativeCall::new();
        // SAFETY: the output pointer is valid and the handle remains live.
        let status = unsafe {
            bsbit_hts_indexed_reader_reference_count(
                self.handle.as_ptr(),
                &raw mut count,
                &raw mut call.system_errno,
                call.error.as_mut_ptr(),
                call.error.len(),
            )
        };
        call.finish(status)?;
        let count = usize::try_from(count)
            .map_err(|_| NativeError::invalid_state("negative BAM reference count"))?;
        let mut references = Vec::with_capacity(count);
        for ordinal in 0..count {
            let reference_id = i32::try_from(ordinal)
                .map_err(|_| NativeError::invalid_state("BAM reference id exceeds i32"))?;
            let mut name = core::ptr::null();
            let mut name_length = 0_usize;
            let mut length = 0_i64;
            let mut call = NativeCall::new();
            // SAFETY: all output pointers are valid; returned storage is copied below.
            let status = unsafe {
                bsbit_hts_indexed_reader_reference(
                    self.handle.as_ptr(),
                    reference_id,
                    &raw mut name,
                    &raw mut name_length,
                    &raw mut length,
                    &raw mut call.system_errno,
                    call.error.as_mut_ptr(),
                    call.error.len(),
                )
            };
            call.finish(status)?;
            // SAFETY: the shim guarantees the name is readable until the next
            // mutable operation; it is copied before continuing.
            let name =
                unsafe { copy_native_slice(name.cast::<u8>(), name_length, "BAM reference name")? };
            references.push(NativeBamReference { name, length });
        }
        Ok(references)
    }

    /// Selects one zero-based, half-open reference interval.
    pub fn query(&mut self, reference_id: i32, start: i64, end: i64) -> Result<(), NativeError> {
        let mut call = NativeCall::new();
        // SAFETY: the handle is exclusively borrowed for iterator replacement.
        let status = unsafe {
            bsbit_hts_indexed_reader_query(
                self.handle.as_ptr(),
                reference_id,
                start,
                end,
                &raw mut call.system_errno,
                call.error.as_mut_ptr(),
                call.error.len(),
            )
        };
        call.finish(status)
    }

    /// Borrows the next record overlapping the selected interval.
    pub fn next_record(&mut self) -> Result<Option<NativeIndexedBamRecordView<'_>>, NativeError> {
        let mut view = CBamRecordView::empty();
        let mut has_record = 0;
        let mut call = NativeCall::new();
        // SAFETY: the view and diagnostic outputs are writable and the handle
        // is exclusively borrowed for native iteration.
        let status = unsafe {
            bsbit_hts_indexed_reader_next(
                self.handle.as_ptr(),
                &raw mut view,
                &raw mut has_record,
                &raw mut call.system_errno,
                call.error.as_mut_ptr(),
                call.error.len(),
            )
        };
        call.finish(status)?;
        match has_record {
            0 => Ok(None),
            1 => {
                let expected_packed =
                    view.sequence_length.checked_add(1).ok_or_else(|| {
                        NativeError::invalid_state("BAM sequence length overflow")
                    })? / 2;
                if view.packed_sequence_length != expected_packed {
                    return Err(NativeError::invalid_state(
                        "shim returned an inconsistent packed BAM sequence length",
                    ));
                }
                // SAFETY: every pointer/length pair is borrowed from the live
                // bam1_t. The returned view is tied to this exclusive reader
                // borrow, so another mutable native call cannot invalidate it.
                let query_name = unsafe {
                    borrow_native_slice(
                        view.query_name.cast::<u8>(),
                        view.query_name_length,
                        "BAM query name",
                    )?
                };
                // SAFETY: same record-lifetime guarantee as above.
                let cigar =
                    unsafe { borrow_native_slice(view.cigar, view.cigar_count, "BAM CIGAR")? };
                // SAFETY: same record-lifetime guarantee as above.
                let packed_sequence = unsafe {
                    borrow_native_slice(
                        view.sequence,
                        view.packed_sequence_length,
                        "packed BAM sequence",
                    )?
                };
                // SAFETY: same record-lifetime guarantee as above.
                let quality = unsafe {
                    borrow_native_slice(view.quality, view.sequence_length, "BAM quality")?
                };
                // SAFETY: same record-lifetime guarantee as above.
                let auxiliary = unsafe {
                    borrow_native_slice(
                        view.auxiliary,
                        view.auxiliary_length,
                        "BAM auxiliary data",
                    )?
                };
                Ok(Some(NativeIndexedBamRecordView {
                    reference_id: view.reference_id,
                    position: view.position,
                    mapping_quality: view.mapping_quality,
                    flag: view.flag,
                    mate_reference_id: view.mate_reference_id,
                    mate_position: view.mate_position,
                    template_length: view.template_length,
                    query_name,
                    cigar,
                    packed_sequence,
                    sequence_length: view.sequence_length,
                    quality,
                    auxiliary,
                }))
            }
            _ => Err(NativeError::invalid_state(
                "shim returned an invalid record-presence sentinel",
            )),
        }
    }

    /// Closes the native BAM, header, index, iterator, and record storage.
    pub fn close(&mut self) -> Result<(), NativeError> {
        if self.closed {
            return Err(NativeError::invalid_state(
                "indexed BAM reader is already closed",
            ));
        }
        let mut call = NativeCall::new();
        // SAFETY: the handle is live and this method marks it closed exactly once.
        let status = unsafe {
            bsbit_hts_indexed_reader_close(
                self.handle.as_ptr(),
                &raw mut call.system_errno,
                call.error.as_mut_ptr(),
                call.error.len(),
            )
        };
        self.closed = true;
        call.finish(status)
    }
}

impl Drop for NativeIndexedBamReader {
    fn drop(&mut self) {
        if !self.closed {
            // SAFETY: Drop owns the live handle; null diagnostics are accepted.
            unsafe {
                let _ = bsbit_hts_indexed_reader_close(
                    self.handle.as_ptr(),
                    core::ptr::null_mut(),
                    core::ptr::null_mut(),
                    0,
                );
            }
        }
        // SAFETY: this is the sole owner and destroy accepts a closed handle.
        unsafe { bsbit_hts_indexed_reader_destroy(self.handle.as_ptr()) };
    }
}

/// One copied indexed FASTA dictionary entry.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct NativeFastaReference {
    pub name: Vec<u8>,
    pub length: i64,
}

/// One owned, thread-confined indexed FASTA handle.
pub struct NativeIndexedFastaReader {
    handle: NonNull<CIndexedFastaReader>,
    closed: bool,
    _thread_confined: PhantomData<Rc<()>>,
}

impl NativeIndexedFastaReader {
    /// Opens a FASTA with its existing adjacent `.fai` and optional `.gzi` index.
    pub fn open(path: &CStr) -> Result<Self, NativeError> {
        let mut handle = core::ptr::null_mut();
        let mut call = NativeCall::new();
        // SAFETY: all output pointers are valid and `path` is NUL-terminated.
        let status = unsafe {
            bsbit_hts_indexed_fasta_reader_open(
                path.as_ptr(),
                &raw mut handle,
                &raw mut call.system_errno,
                call.error.as_mut_ptr(),
                call.error.len(),
            )
        };
        call.finish(status)?;
        let handle = NonNull::new(handle).ok_or_else(|| {
            NativeError::invalid_state("indexed FASTA open returned a null handle")
        })?;
        Ok(Self {
            handle,
            closed: false,
            _thread_confined: PhantomData,
        })
    }

    /// Copies the FASTA dictionary owned by the `HTSlib` index handle.
    pub fn references(&self) -> Result<Vec<NativeFastaReference>, NativeError> {
        let mut count = 0_i32;
        let mut call = NativeCall::new();
        // SAFETY: the output pointer is valid and the handle remains live.
        let status = unsafe {
            bsbit_hts_indexed_fasta_reader_reference_count(
                self.handle.as_ptr(),
                &raw mut count,
                &raw mut call.system_errno,
                call.error.as_mut_ptr(),
                call.error.len(),
            )
        };
        call.finish(status)?;
        let count = usize::try_from(count)
            .map_err(|_| NativeError::invalid_state("negative FASTA reference count"))?;
        let mut references = Vec::with_capacity(count);
        for ordinal in 0..count {
            let reference_id = i32::try_from(ordinal)
                .map_err(|_| NativeError::invalid_state("FASTA reference id exceeds i32"))?;
            let mut name = core::ptr::null();
            let mut name_length = 0_usize;
            let mut length = 0_i64;
            let mut call = NativeCall::new();
            // SAFETY: all output pointers are valid; returned storage is copied below.
            let status = unsafe {
                bsbit_hts_indexed_fasta_reader_reference(
                    self.handle.as_ptr(),
                    reference_id,
                    &raw mut name,
                    &raw mut name_length,
                    &raw mut length,
                    &raw mut call.system_errno,
                    call.error.as_mut_ptr(),
                    call.error.len(),
                )
            };
            call.finish(status)?;
            // SAFETY: the name is owned by the live faidx handle and copied now.
            let name = unsafe {
                copy_native_slice(name.cast::<u8>(), name_length, "FASTA reference name")?
            };
            references.push(NativeFastaReference { name, length });
        }
        Ok(references)
    }

    /// Copies one zero-based, half-open FASTA interval.
    pub fn fetch(
        &mut self,
        reference_id: i32,
        start: i64,
        end: i64,
    ) -> Result<Vec<u8>, NativeError> {
        let mut sequence = core::ptr::null();
        let mut length = 0_usize;
        let mut call = NativeCall::new();
        // SAFETY: all output pointers are valid; returned handle storage is copied below.
        let status = unsafe {
            bsbit_hts_indexed_fasta_reader_fetch(
                self.handle.as_ptr(),
                reference_id,
                start,
                end,
                &raw mut sequence,
                &raw mut length,
                &raw mut call.system_errno,
                call.error.as_mut_ptr(),
                call.error.len(),
            )
        };
        call.finish(status)?;
        // SAFETY: the sequence remains readable until the next mutable handle call.
        unsafe { copy_native_slice(sequence.cast::<u8>(), length, "FASTA interval") }
    }

    /// Closes the FASTA and index resources.
    pub fn close(&mut self) -> Result<(), NativeError> {
        if self.closed {
            return Err(NativeError::invalid_state(
                "indexed FASTA reader is already closed",
            ));
        }
        let mut call = NativeCall::new();
        // SAFETY: the handle is live and marked closed exactly once below.
        let status = unsafe {
            bsbit_hts_indexed_fasta_reader_close(
                self.handle.as_ptr(),
                &raw mut call.system_errno,
                call.error.as_mut_ptr(),
                call.error.len(),
            )
        };
        self.closed = true;
        call.finish(status)
    }
}

impl Drop for NativeIndexedFastaReader {
    fn drop(&mut self) {
        if !self.closed {
            // SAFETY: Drop owns the live handle; null diagnostics are accepted.
            unsafe {
                let _ = bsbit_hts_indexed_fasta_reader_close(
                    self.handle.as_ptr(),
                    core::ptr::null_mut(),
                    core::ptr::null_mut(),
                    0,
                );
            }
        }
        // SAFETY: this is the sole owner and destroy accepts a closed handle.
        unsafe { bsbit_hts_indexed_fasta_reader_destroy(self.handle.as_ptr()) };
    }
}

unsafe fn copy_native_slice<T: Copy>(
    pointer: *const T,
    length: usize,
    label: &str,
) -> Result<Vec<T>, NativeError> {
    if length == 0 {
        return Ok(Vec::new());
    }
    if pointer.is_null() {
        return Err(NativeError::invalid_state(&format!(
            "shim returned a null nonempty {label}"
        )));
    }
    // SAFETY: callers establish that `pointer` is readable for `length`
    // elements during this call, and the non-null case was checked above.
    let slice = unsafe { core::slice::from_raw_parts(pointer, length) };
    Ok(slice.to_vec())
}

unsafe fn borrow_native_slice<'a, T>(
    pointer: *const T,
    length: usize,
    label: &str,
) -> Result<&'a [T], NativeError> {
    if length == 0 {
        return Ok(&[]);
    }
    if pointer.is_null() {
        return Err(NativeError::invalid_state(&format!(
            "shim returned a null nonempty {label}"
        )));
    }
    // SAFETY: callers tie the returned lifetime to the exclusive borrow of
    // the native owner and guarantee readability for `length` elements.
    Ok(unsafe { core::slice::from_raw_parts(pointer, length) })
}

/// One owned, thread-confined native BAM encode handle.
pub struct NativeBamWriter {
    handle: NonNull<CWriter>,
    finished: bool,
    _thread_confined: PhantomData<Rc<()>>,
}

/// Borrowed, already validated fields for one direct BAM record.
pub struct NativeBamRecordFields<'a> {
    /// Query name without a terminating NUL.
    pub query_name: &'a [u8],
    /// Complete BAM flag word.
    pub flag: u16,
    /// Zero-based reference dictionary ordinal, or -1 when unmapped.
    pub reference_id: i32,
    /// Zero-based reference position, or -1 when unmapped.
    pub position: i64,
    /// Numeric mapping quality.
    pub mapping_quality: u8,
    /// Packed BAM CIGAR words.
    pub cigar: &'a [u32],
    /// Zero-based mate reference dictionary ordinal, or -1 when absent.
    pub mate_reference_id: i32,
    /// Zero-based mate position, or -1 when absent.
    pub mate_position: i64,
    /// Signed observed template length.
    pub template_length: i64,
    /// Query sequence as canonical ASCII bases.
    pub sequence: &'a [u8],
    /// Printable Phred+33 qualities converted in place after native packing.
    pub quality: Option<&'a [u8]>,
    /// Literal NM and canonical MD when the record is mapped.
    pub literal_nm_and_md: Option<(u32, &'a [u8])>,
    /// Whether the canonical MD value is serialized for a mapped record.
    pub emit_md: bool,
    /// Bisulfite genome conversion (XG) when requested for a mapped record.
    pub bisulfite_genome_conversion: Option<&'static [u8; 2]>,
    /// Bismark XM and XR values when the compatibility contract is active.
    pub bismark_auxiliary: Option<(&'a [u8], &'static [u8; 2])>,
}

impl NativeBamWriter {
    /// Opens a BAM path and writes one canonical SAM header.
    ///
    /// # Errors
    ///
    /// Returns a copied native error and owns no handle on failure.
    #[cfg(test)]
    pub fn open(path: &CStr, header: &[u8]) -> Result<Self, NativeError> {
        Self::open_with_threads(path, header, 0)
    }

    /// Opens a BAM path with private `HTSlib` BGZF compression workers.
    ///
    /// `compression_threads == 0` preserves the synchronous writer. Values
    /// above 64 are rejected by the audited shim.
    ///
    /// # Errors
    ///
    /// Returns a copied native error and owns no handle on failure.
    pub fn open_with_threads(
        path: &CStr,
        header: &[u8],
        compression_threads: u32,
    ) -> Result<Self, NativeError> {
        let mut handle = core::ptr::null_mut();
        let mut call = NativeCall::new();
        // SAFETY: path/header pointers remain valid for the call and outputs are writable.
        let status = unsafe {
            bsbit_hts_writer_open_bam_threads(
                path.as_ptr(),
                header.as_ptr().cast(),
                header.len(),
                compression_threads,
                &raw mut handle,
                &raw mut call.system_errno,
                call.error.as_mut_ptr(),
                call.error.len(),
            )
        };
        call.finish(status)?;
        let handle = NonNull::new(handle)
            .ok_or_else(|| NativeError::invalid_state("writer open returned a null handle"))?;
        Ok(Self {
            handle,
            finished: false,
            _thread_confined: PhantomData,
        })
    }

    /// Opens a BAM path with private BGZF workers and an explicit compression level.
    ///
    /// # Errors
    ///
    /// Returns a copied native error when the level is outside `0..=9` or the
    /// writer cannot be created.
    pub fn open_with_threads_and_compression_level(
        path: &CStr,
        header: &[u8],
        compression_threads: u32,
        compression_level: u8,
    ) -> Result<Self, NativeError> {
        let mut handle = core::ptr::null_mut();
        let mut call = NativeCall::new();
        // SAFETY: path/header pointers remain valid for the call and outputs are writable.
        let status = unsafe {
            bsbit_hts_writer_open_bam_threads_level(
                path.as_ptr(),
                header.as_ptr().cast(),
                header.len(),
                compression_threads,
                c_int::from(compression_level),
                &raw mut handle,
                &raw mut call.system_errno,
                call.error.as_mut_ptr(),
                call.error.len(),
            )
        };
        call.finish(status)?;
        let handle = NonNull::new(handle)
            .ok_or_else(|| NativeError::invalid_state("writer open returned a null handle"))?;
        Ok(Self {
            handle,
            finished: false,
            _thread_confined: PhantomData,
        })
    }

    /// Encodes one canonical LF-terminated SAM record.
    ///
    /// # Errors
    ///
    /// Returns a copied native error; native record errors are terminal.
    pub fn write_record(&mut self, record: &[u8]) -> Result<(), NativeError> {
        let mut call = NativeCall::new();
        // SAFETY: record storage is readable for its exact length and handle ownership is exclusive.
        let status = unsafe {
            bsbit_hts_writer_write_record(
                self.handle.as_ptr(),
                record.as_ptr().cast(),
                record.len(),
                &raw mut call.system_errno,
                call.error.as_mut_ptr(),
                call.error.len(),
            )
        };
        call.finish(status)
    }

    /// Writes validated fields without a SAM text round trip.
    ///
    /// # Errors
    ///
    /// Returns a copied native error; invalid fields and write errors are
    /// terminal for the native writer.
    pub fn write_bam_fields(
        &mut self,
        fields: &NativeBamRecordFields<'_>,
    ) -> Result<(), NativeError> {
        let mut call = NativeCall::new();
        let quality = fields.quality.map_or(core::ptr::null(), <[u8]>::as_ptr);
        let (has_mapping, literal_nm, has_md, md, md_length) =
            fields
                .literal_nm_and_md
                .map_or((0, 0, 0, core::ptr::null(), 0), |(literal_nm, md)| {
                    if fields.emit_md {
                        (1, literal_nm, 1, md.as_ptr().cast(), md.len())
                    } else {
                        (1, literal_nm, 0, core::ptr::null(), 0)
                    }
                });
        let (has_xg, xg) = fields
            .bisulfite_genome_conversion
            .map_or((0, core::ptr::null()), |xg| (1, xg.as_ptr().cast()));
        let (has_bismark, bismark_xm, bismark_xm_length, bismark_read_conversion) = fields
            .bismark_auxiliary
            .map_or((0, core::ptr::null(), 0, core::ptr::null()), |(xm, xr)| {
                (1, xm.as_ptr().cast(), xm.len(), xr.as_ptr().cast())
            });
        // SAFETY: all borrowed slices remain readable for the call, optional
        // pointers are null exactly when absent, and handle ownership is exclusive.
        let status = unsafe {
            bsbit_hts_writer_write_bam_fields(
                self.handle.as_ptr(),
                fields.query_name.as_ptr().cast(),
                fields.query_name.len(),
                fields.flag,
                fields.reference_id,
                fields.position,
                fields.mapping_quality,
                fields.cigar.as_ptr(),
                fields.cigar.len(),
                fields.mate_reference_id,
                fields.mate_position,
                fields.template_length,
                fields.sequence.as_ptr().cast(),
                fields.sequence.len(),
                quality,
                has_mapping,
                literal_nm,
                has_md,
                md,
                md_length,
                has_xg,
                xg,
                has_bismark,
                bismark_xm,
                bismark_xm_length,
                bismark_read_conversion,
                &raw mut call.system_errno,
                call.error.as_mut_ptr(),
                call.error.len(),
            )
        };
        call.finish(status)
    }

    /// Finalizes the BAM stream exactly once.
    ///
    /// # Errors
    ///
    /// Returns a copied close or poisoned-state error. The stream remains
    /// terminal even when finalization fails.
    pub fn finish(&mut self) -> Result<(), NativeError> {
        if self.finished {
            return Err(NativeError::invalid_state("writer is already finished"));
        }
        let mut call = NativeCall::new();
        // SAFETY: the live handle is exclusively owned and finalization is invoked once.
        let status = unsafe {
            bsbit_hts_writer_finish(
                self.handle.as_ptr(),
                &raw mut call.system_errno,
                call.error.as_mut_ptr(),
                call.error.len(),
            )
        };
        self.finished = true;
        call.finish(status)
    }
}

impl Drop for NativeBamWriter {
    fn drop(&mut self) {
        if !self.finished {
            // SAFETY: Drop owns the live handle; null diagnostics are accepted by the shim.
            unsafe {
                let _ = bsbit_hts_writer_finish(
                    self.handle.as_ptr(),
                    core::ptr::null_mut(),
                    core::ptr::null_mut(),
                    0,
                );
            }
        }
        // SAFETY: this is the sole owner and destroy is valid after finalize or failure.
        unsafe { bsbit_hts_writer_destroy(self.handle.as_ptr()) };
    }
}

struct NativeCall {
    system_errno: c_int,
    error: NativeErrorBuffer,
}

struct NativeErrorBuffer(MaybeUninit<[c_char; ERROR_CAPACITY]>);

impl NativeErrorBuffer {
    const fn new() -> Self {
        Self(MaybeUninit::uninit())
    }

    fn as_mut_ptr(&mut self) -> *mut c_char {
        self.0.as_mut_ptr().cast()
    }

    const fn len(&self) -> usize {
        core::mem::size_of_val(&self.0)
    }

    fn message_bytes(&self) -> &[u8] {
        // SAFETY: every shim entry point calls `set_result`, which initializes
        // at least the first byte and always terminates copied error text. This
        // method is reached only after one such call returned an error status.
        unsafe { CStr::from_ptr(self.0.as_ptr().cast()).to_bytes() }
    }
}

impl NativeCall {
    const fn new() -> Self {
        Self {
            system_errno: 0,
            error: NativeErrorBuffer::new(),
        }
    }

    fn finish(self, raw_status: c_int) -> Result<(), NativeError> {
        let Some(status) = NativeStatus::from_raw(raw_status) else {
            return Ok(());
        };
        let bytes = self.error.message_bytes();
        Err(NativeError {
            status,
            system_errno: self.system_errno,
            message: String::from_utf8_lossy(bytes).into_owned(),
        })
    }
}

#[cfg(test)]
mod tests {
    use std::ffi::CString;
    use std::fs;
    use std::time::{SystemTime, UNIX_EPOCH};

    use super::*;

    fn unique_path(label: &str) -> std::path::PathBuf {
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock")
            .as_nanos();
        std::env::temp_dir().join(format!("bsbit-hts-{label}-{}-{nonce}", std::process::id()))
    }

    fn c_path(path: &std::path::Path) -> CString {
        CString::new(path.to_str().expect("test path is UTF-8")).expect("no NUL")
    }

    #[test]
    fn exact_abi_runtime_and_health_are_available() {
        assert_eq!(shim_abi_version(), 3);
        assert_eq!(runtime_version().expect("runtime version"), "1.24");
        health_check().expect("exact runtime passes");
    }

    #[test]
    fn owned_reader_decodes_plain_bytes_and_closes_once() {
        let path = unique_path("reader");
        fs::write(&path, b"ACGT").expect("fixture");
        let mut reader = NativeReader::open(&c_path(&path)).expect("open");
        assert_eq!(
            reader.compression().expect("compression"),
            NativeCompression::Plain
        );
        let mut bytes = [0_u8; 8];
        assert_eq!(reader.read(&mut bytes).expect("read"), 4);
        assert_eq!(&bytes[..4], b"ACGT");
        assert_eq!(reader.read(&mut bytes).expect("EOF"), 0);
        reader.close().expect("close");
        assert_eq!(
            reader.close().expect_err("repeat close").status(),
            NativeStatus::InvalidState
        );
        fs::remove_file(path).expect("cleanup");
    }

    #[test]
    fn owned_writer_encodes_one_canonical_record() {
        const HEADER: &[u8] = b"@HD\tVN:1.6\tSO:unknown\n@SQ\tSN:chr1\tLN:100\n";
        const RECORD: &[u8] = b"read1\t0\tchr1\t2\t255\t4M\t*\t0\t0\tACGT\tIIII\tNM:i:0\tMD:Z:4\n";
        let path = unique_path("writer.bam");
        let mut writer = NativeBamWriter::open(&c_path(&path), HEADER).expect("open");
        writer.write_record(RECORD).expect("record");
        writer.finish().expect("finish");
        assert_eq!(
            writer.finish().expect_err("repeat finish").status(),
            NativeStatus::InvalidState
        );
        assert!(fs::metadata(&path).expect("BAM exists").len() > 0);
        fs::remove_file(path).expect("cleanup");
    }

    #[test]
    fn direct_fields_are_byte_identical_to_the_canonical_sam_path() {
        const HEADER: &[u8] = b"@HD\tVN:1.6\tSO:unknown\n@SQ\tSN:chr1\tLN:100\n";
        const RECORD: &[u8] = b"read1\t0\tchr1\t2\t255\t4M\t*\t0\t0\tACGT\tIIII\tNM:i:0\tMD:Z:4\n";
        let text_path = unique_path("writer-text.bam");
        let direct_path = unique_path("writer-direct.bam");

        let mut text = NativeBamWriter::open(&c_path(&text_path), HEADER).expect("text open");
        text.write_record(RECORD).expect("text record");
        text.finish().expect("text finish");

        let cigar = [4_u32 << 4];
        let fields = NativeBamRecordFields {
            query_name: b"read1",
            flag: 0,
            reference_id: 0,
            position: 1,
            mapping_quality: 255,
            cigar: &cigar,
            mate_reference_id: -1,
            mate_position: -1,
            template_length: 0,
            sequence: b"ACGT",
            quality: Some(b"IIII"),
            literal_nm_and_md: Some((0, b"4")),
            emit_md: true,
            bisulfite_genome_conversion: None,
            bismark_auxiliary: None,
        };
        let mut direct = NativeBamWriter::open(&c_path(&direct_path), HEADER).expect("direct open");
        direct.write_bam_fields(&fields).expect("direct record");
        direct.finish().expect("direct finish");

        assert_eq!(
            fs::read(&direct_path).expect("direct BAM"),
            fs::read(&text_path).expect("text BAM")
        );
        fs::remove_file(text_path).expect("text cleanup");
        fs::remove_file(direct_path).expect("direct cleanup");
    }

    #[test]
    fn owned_writer_accepts_private_bgzf_workers() {
        const HEADER: &[u8] = b"@HD\tVN:1.6\tSO:unknown\n@SQ\tSN:chr1\tLN:100\n";
        const RECORD: &[u8] = b"read1\t0\tchr1\t2\t255\t4M\t*\t0\t0\tACGT\tIIII\tNM:i:0\tMD:Z:4\n";
        let path = unique_path("threaded-writer.bam");
        let mut writer =
            NativeBamWriter::open_with_threads(&c_path(&path), HEADER, 2).expect("open");
        writer.write_record(RECORD).expect("record");
        writer.finish().expect("finish");
        assert!(fs::metadata(&path).expect("BAM exists").len() > 0);
        fs::remove_file(path).expect("cleanup");
    }

    #[test]
    fn owned_writer_accepts_level_one_and_rejects_level_ten() {
        const HEADER: &[u8] = b"@HD\tVN:1.6\tSO:unknown\n@SQ\tSN:chr1\tLN:100\n";
        const RECORD: &[u8] = b"read1\t0\tchr1\t2\t255\t4M\t*\t0\t0\tACGT\tIIII\tNM:i:0\tMD:Z:4\n";
        let accepted = unique_path("level-one.bam");
        let mut writer = NativeBamWriter::open_with_threads_and_compression_level(
            &c_path(&accepted),
            HEADER,
            2,
            1,
        )
        .expect("level one opens");
        writer.write_record(RECORD).expect("record");
        writer.finish().expect("finish");
        assert!(fs::metadata(&accepted).expect("BAM exists").len() > 0);

        let rejected = unique_path("level-ten.bam");
        let error = NativeBamWriter::open_with_threads_and_compression_level(
            &c_path(&rejected),
            HEADER,
            2,
            10,
        )
        .err()
        .expect("level ten is rejected");
        assert_eq!(error.status(), NativeStatus::InvalidArgument);
        assert!(!rejected.exists());
        fs::remove_file(accepted).expect("cleanup");
    }
}
