//! Safe ownership for immutable Linux file mappings.

use core::ffi::{c_int, c_long, c_void};
use core::mem::size_of;
use core::ptr::NonNull;
use std::fs::File;
use std::io;
use std::os::fd::{AsRawFd, RawFd};

const PROT_READ: c_int = 0x1;
const MAP_PRIVATE: c_int = 0x2;
const MADV_HUGEPAGE: c_int = 14;
const MAP_FAILED: *mut c_void = usize::MAX as *mut c_void;

unsafe extern "C" {
    fn mmap(
        address: *mut c_void,
        length: usize,
        protection: c_int,
        flags: c_int,
        descriptor: RawFd,
        offset: c_long,
    ) -> *mut c_void;

    fn madvise(address: *mut c_void, length: usize, advice: c_int) -> c_int;
    fn munmap(address: *mut c_void, length: usize) -> c_int;
}

/// An owned, immutable mapping of one complete nonempty local file.
///
/// The mapping owns no Rust references and may be shared between worker
/// threads. Accessors perform a release-build bounds assertion before the
/// private pointer read, so malformed offsets cannot cause undefined behavior.
#[derive(Debug)]
pub struct ReadOnlyMapping {
    pointer: NonNull<u8>,
    length: usize,
}

// SAFETY: the mapping is immutable for its complete lifetime and owns no Rust
// references. `munmap` runs only after ownership is dropped.
unsafe impl Send for ReadOnlyMapping {}
// SAFETY: all shared accessors validate their range before reading immutable
// bytes from the live mapping.
unsafe impl Sync for ReadOnlyMapping {}

impl ReadOnlyMapping {
    /// Maps one complete nonempty file read-only and privately.
    ///
    /// # Errors
    ///
    /// Returns `InvalidData` for an empty file, `FileTooLarge` when its length
    /// cannot fit this process, or the direct mapping error.
    pub fn map(file: &File) -> io::Result<Self> {
        let length = usize::try_from(file.metadata()?.len())
            .map_err(|_| io::Error::from(io::ErrorKind::FileTooLarge))?;
        if length == 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "cannot map an empty file",
            ));
        }
        // SAFETY: the descriptor is live for this call, `length` came from
        // that descriptor, and no writable or shared mapping is requested.
        let mapped = unsafe {
            mmap(
                core::ptr::null_mut(),
                length,
                PROT_READ,
                MAP_PRIVATE,
                file.as_raw_fd(),
                0,
            )
        };
        if mapped == MAP_FAILED {
            return Err(io::Error::last_os_error());
        }
        let Some(pointer) = NonNull::new(mapped.cast::<u8>()) else {
            // SAFETY: even a null successful address must be released before
            // this safe abstraction rejects it as unusable for Rust slices.
            let _ = unsafe { munmap(mapped, length) };
            return Err(io::Error::other("mmap returned a null address"));
        };
        // Advisory only. Unsupported kernels or filesystems retain ordinary
        // demand paging without changing mapping semantics.
        // SAFETY: this is the live mapping returned immediately above.
        let _ = unsafe { madvise(mapped, length, MADV_HUGEPAGE) };
        Ok(Self { pointer, length })
    }

    /// Returns the mapped byte length.
    #[must_use]
    pub const fn len(&self) -> usize {
        self.length
    }

    /// Returns whether the mapping is empty.
    ///
    /// Complete file mappings are never empty; this method is provided for
    /// ordinary slice-like inspection.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        false
    }

    /// Borrows the complete immutable byte mapping.
    #[must_use]
    #[inline]
    pub fn as_slice(&self) -> &[u8] {
        // SAFETY: construction accepts only a non-null pointer and nonzero
        // length from mmap, and the borrow cannot outlive this owner.
        unsafe { core::slice::from_raw_parts(self.pointer.as_ptr(), self.length) }
    }

    /// Reads one byte after checking the mapped range.
    ///
    /// # Panics
    ///
    /// Panics when `offset` is outside this mapping.
    #[must_use]
    #[inline]
    pub fn read_u8(&self, offset: usize) -> u8 {
        assert!(offset < self.length, "mapped byte offset is out of bounds");
        // SAFETY: the assertion proves this byte lies inside the live mapping.
        unsafe { self.pointer.as_ptr().add(offset).read() }
    }

    /// Reads one unaligned little-endian `u32` after checking the range.
    ///
    /// # Panics
    ///
    /// Panics unless the complete four-byte value lies inside this mapping.
    #[must_use]
    #[inline]
    pub fn read_u32(&self, offset: usize) -> u32 {
        self.assert_range(offset, size_of::<u32>());
        // SAFETY: the checked range contains four live bytes; unaligned reads
        // impose no stronger address requirement.
        u32::from_le(unsafe {
            self.pointer
                .as_ptr()
                .add(offset)
                .cast::<u32>()
                .read_unaligned()
        })
    }

    /// Reads one unaligned little-endian `u64` after checking the range.
    ///
    /// # Panics
    ///
    /// Panics unless the complete eight-byte value lies inside this mapping.
    #[must_use]
    #[inline]
    pub fn read_u64(&self, offset: usize) -> u64 {
        self.assert_range(offset, size_of::<u64>());
        // SAFETY: the checked range contains eight live bytes.
        u64::from_le(unsafe {
            self.pointer
                .as_ptr()
                .add(offset)
                .cast::<u64>()
                .read_unaligned()
        })
    }

    #[inline]
    fn assert_range(&self, offset: usize, width: usize) {
        assert!(
            offset
                .checked_add(width)
                .is_some_and(|end| end <= self.length),
            "mapped integer range is out of bounds"
        );
    }
}

impl Drop for ReadOnlyMapping {
    fn drop(&mut self) {
        // SAFETY: this is the exact pointer and nonzero length returned by the
        // successful mmap call, released exactly once by this owner.
        let _ = unsafe { munmap(self.pointer.as_ptr().cast(), self.length) };
    }
}
