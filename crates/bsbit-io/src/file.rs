//! Identity-aware filesystem primitives shared by format and index code.
//!
//! This module owns only local-file mechanics. It does not know which bytes a
//! file contains or which domain produced them.

use core::ffi::{c_char, c_int, c_long, c_void};
use std::ffi::CString;
use std::fs::{self, File, Metadata, OpenOptions};
use std::io;
use std::os::fd::{AsRawFd, BorrowedFd};
use std::os::unix::ffi::OsStrExt;
use std::os::unix::fs::{MetadataExt, OpenOptionsExt};
use std::path::{Path, PathBuf};

#[cfg(not(target_os = "linux"))]
compile_error!("bsbit-io currently supports only the audited Linux profile");

const AT_CURRENT_WORKING_DIRECTORY: c_int = -100;
const AT_SYMLINK_FOLLOW: c_int = 0x400;
const O_NOFOLLOW: c_int = 0x20_000;
const O_NONBLOCK: c_int = 0x800;
const V9FS_SUPER_MAGIC: c_long = 0x0102_1997;
const STATFS_STORAGE_WORDS: usize = 32;

unsafe extern "C" {
    fn fstatfs(descriptor: c_int, statistics: *mut c_void) -> c_int;

    fn linkat(
        old_directory: c_int,
        old_path: *const c_char,
        new_directory: c_int,
        new_path: *const c_char,
        flags: c_int,
    ) -> c_int;
}

/// Stable identity of one Unix filesystem object.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct FileIdentity {
    device: u64,
    inode: u64,
}

impl FileIdentity {
    /// Captures the identity represented by metadata.
    #[must_use]
    pub fn from_metadata(metadata: &Metadata) -> Self {
        Self {
            device: metadata.dev(),
            inode: metadata.ino(),
        }
    }

    /// Captures the identity of an open file descriptor.
    ///
    /// # Errors
    ///
    /// Returns the underlying metadata error.
    pub fn from_file(file: &File) -> io::Result<Self> {
        file.metadata()
            .map(|metadata| Self::from_metadata(&metadata))
    }

    /// Reports whether `path` still names this exact filesystem object.
    ///
    /// A missing path returns `false`; a replacement is never treated as the
    /// object originally opened by the caller.
    ///
    /// # Errors
    ///
    /// Returns metadata errors other than a missing path.
    pub fn matches_path(self, path: &Path) -> io::Result<bool> {
        match fs::symlink_metadata(path) {
            Ok(metadata) => Ok(self == Self::from_metadata(&metadata)),
            Err(source) if source.kind() == io::ErrorKind::NotFound => Ok(false),
            Err(source) => Err(source),
        }
    }
}

/// Converts a caller path to a lexical absolute path without resolving links.
///
/// # Errors
///
/// Returns the current-directory lookup error for a relative path.
pub fn absolute_path(path: &Path) -> io::Result<PathBuf> {
    if path.is_absolute() {
        Ok(path.to_path_buf())
    } else {
        std::env::current_dir().map(|directory| directory.join(path))
    }
}

/// Verifies that a future create-only target is absent, including symlinks.
///
/// # Errors
///
/// Returns `AlreadyExists` for any existing directory entry and propagates
/// other metadata errors.
pub fn validate_absent(path: &Path) -> io::Result<()> {
    match fs::symlink_metadata(path) {
        Err(source) if source.kind() == io::ErrorKind::NotFound => Ok(()),
        Ok(_) => Err(io::Error::from(io::ErrorKind::AlreadyExists)),
        Err(source) => Err(source),
    }
}

/// Verifies that an existing local path is a regular file while permitting a
/// missing path to pass through to the format-specific opener.
///
/// This is useful at adapter boundaries that must reject directories and
/// special files without replacing the richer missing-file diagnostic from a
/// downstream codec or native library.
///
/// # Errors
///
/// Returns `Unsupported` when an existing path is not a regular file and
/// propagates metadata failures other than `NotFound`.
pub fn validate_regular_file_or_absent(path: &Path) -> io::Result<()> {
    match fs::metadata(path) {
        Ok(metadata) if metadata.is_file() => Ok(()),
        Ok(_) => Err(io::Error::new(
            io::ErrorKind::Unsupported,
            "existing path is not a regular file",
        )),
        Err(source) if source.kind() == io::ErrorKind::NotFound => Ok(()),
        Err(source) => Err(source),
    }
}

/// Verifies that a future create-only file target is absent and has an
/// existing directory as its parent.
///
/// The parent lookup follows directory symlinks in the same way as an
/// eventual file creation. The target itself is inspected without following
/// links, so a dangling symlink still counts as occupied.
///
/// # Errors
///
/// Returns `AlreadyExists` for an occupied target, `NotADirectory` when the
/// parent exists but is not a directory, and propagates absolutization or
/// metadata errors otherwise.
pub fn validate_create_target(path: &Path) -> io::Result<()> {
    let target = absolute_path(path)?;
    validate_absent(&target)?;
    let parent = target.parent().ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            "create-only target has no parent directory",
        )
    })?;
    let metadata = fs::metadata(parent)?;
    if metadata.is_dir() {
        Ok(())
    } else {
        Err(io::Error::new(
            io::ErrorKind::NotADirectory,
            "create-only target parent is not a directory",
        ))
    }
}

/// Verifies that a future replaceable file target has an existing directory
/// as its parent.
///
/// A missing target and an existing regular file or symbolic link are valid.
/// Directories and other special filesystem objects are never replaced.
///
/// # Errors
///
/// Returns `Unsupported` for an existing non-file target, `NotADirectory` when
/// the parent exists but is not a directory, and propagates absolutization or
/// metadata errors otherwise.
pub fn validate_replace_target(path: &Path) -> io::Result<()> {
    let target = absolute_path(path)?;
    match fs::symlink_metadata(&target) {
        Ok(metadata) if metadata.is_file() || metadata.file_type().is_symlink() => {}
        Ok(_) => {
            return Err(io::Error::new(
                io::ErrorKind::Unsupported,
                "existing output target is not a regular file or symbolic link",
            ));
        }
        Err(source) if source.kind() == io::ErrorKind::NotFound => {}
        Err(source) => return Err(source),
    }
    let parent = target.parent().ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            "replacement target has no parent directory",
        )
    })?;
    let metadata = fs::metadata(parent)?;
    if metadata.is_dir() {
        Ok(())
    } else {
        Err(io::Error::new(
            io::ErrorKind::NotADirectory,
            "replacement target parent is not a directory",
        ))
    }
}

/// Verifies that two paths are lexically distinct and, when both exist, do not
/// identify the same filesystem object.
///
/// # Errors
///
/// Returns `InvalidInput` for the same absolute path or filesystem object, or
/// an absolutization/metadata error.
pub fn validate_distinct_paths(first: &Path, second: &Path) -> io::Result<()> {
    if absolute_path(first)? == absolute_path(second)? || paths_share_identity(first, second)? {
        Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "paths resolve to the same filesystem object",
        ))
    } else {
        Ok(())
    }
}

fn paths_share_identity(first: &Path, second: &Path) -> io::Result<bool> {
    let first = match fs::metadata(first) {
        Ok(metadata) => metadata,
        Err(source) if source.kind() == io::ErrorKind::NotFound => return Ok(false),
        Err(source) => return Err(source),
    };
    let second = match fs::metadata(second) {
        Ok(metadata) => metadata,
        Err(source) if source.kind() == io::ErrorKind::NotFound => return Ok(false),
        Err(source) => return Err(source),
    };
    Ok(FileIdentity::from_metadata(&first) == FileIdentity::from_metadata(&second))
}

/// Creates one absent local file without following or replacing an entry.
///
/// # Errors
///
/// Returns the direct `create_new` or metadata error.
pub fn create_new(path: &Path) -> io::Result<(File, FileIdentity)> {
    let file = OpenOptions::new()
        .read(true)
        .write(true)
        .create_new(true)
        .open(path)?;
    let identity = FileIdentity::from_file(&file)?;
    Ok((file, identity))
}

/// Opens a regular output file for direct writing, creating it when absent or
/// truncating it when present.
///
/// Symbolic links, directories, and special files are rejected. The Linux
/// `O_NOFOLLOW` flag closes the race in which the final path becomes a symbolic
/// link between validation and open; `O_NONBLOCK` prevents a raced FIFO from
/// blocking the process before its file type is checked.
///
/// # Errors
///
/// Returns path, parent-directory, permission, open, or file-type errors.
pub fn open_direct_output(path: &Path) -> io::Result<File> {
    open_direct_output_distinct_from(path, &[])
}

/// Opens a regular output for direct writing after proving that it is not any
/// protected input path or existing filesystem object.
///
/// The target is opened without `O_TRUNC`, compared with every protected path
/// by lexical path and file identity, and truncated only after those checks.
///
/// # Errors
///
/// Returns path, parent-directory, permission, file-type, identity, or input
/// collision errors.
pub fn open_direct_output_distinct_from(
    path: &Path,
    protected_paths: &[&Path],
) -> io::Result<File> {
    let target = absolute_path(path)?;
    for protected in protected_paths {
        if target == absolute_path(protected)? {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "output target is also a protected input: {}",
                    protected.display()
                ),
            ));
        }
    }
    let parent = target.parent().ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            "output target has no parent directory",
        )
    })?;
    if !fs::metadata(parent)?.is_dir() {
        return Err(io::Error::new(
            io::ErrorKind::NotADirectory,
            "output target parent is not a directory",
        ));
    }
    match fs::symlink_metadata(&target) {
        Ok(metadata) if metadata.is_file() => {}
        Ok(_) => {
            return Err(io::Error::new(
                io::ErrorKind::Unsupported,
                "existing output target is not a regular file",
            ));
        }
        Err(source) if source.kind() == io::ErrorKind::NotFound => {}
        Err(source) => return Err(source),
    }
    let file = OpenOptions::new()
        .write(true)
        .create(true)
        .mode(0o666)
        .custom_flags(O_NOFOLLOW | O_NONBLOCK)
        .open(&target)?;
    let metadata = file.metadata()?;
    if !metadata.is_file() {
        return Err(io::Error::new(
            io::ErrorKind::Unsupported,
            "opened output target is not a regular file",
        ));
    }
    let identity = FileIdentity::from_metadata(&metadata);
    if !identity.matches_path(&target)? {
        return Err(io::Error::other(
            "output target changed while it was opened",
        ));
    }
    for protected in protected_paths {
        let protected_metadata = match fs::metadata(protected) {
            Ok(metadata) => metadata,
            Err(source) if source.kind() == io::ErrorKind::NotFound => continue,
            Err(source) => return Err(source),
        };
        if identity == FileIdentity::from_metadata(&protected_metadata) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "output target is also a protected input: {}",
                    protected.display()
                ),
            ));
        }
    }
    file.set_len(0)?;
    Ok(file)
}

/// Reopens a live descriptor as an independent read-write file description.
///
/// Unlike [`File::try_clone`], the returned handle has its own cursor. The
/// `/proc/self/fd` lookup binds to the held filesystem object rather than a
/// caller-controlled namespace path, so replacing the original path cannot
/// redirect subsequent writes.
///
/// # Errors
///
/// Returns the direct procfs descriptor-open error.
pub fn reopen_read_write(file: &File) -> io::Result<File> {
    OpenOptions::new()
        .read(true)
        .write(true)
        .open(format!("/proc/self/fd/{}", file.as_raw_fd()))
}

/// Removes `path` only while it still names `identity`.
///
/// # Errors
///
/// Returns `NotFound`/`Other` when ownership was lost and propagates metadata
/// or removal errors. A replacement is never removed.
pub fn remove_if_identity_matches(path: &Path, identity: FileIdentity) -> io::Result<()> {
    if identity.matches_path(path)? {
        fs::remove_file(path)
    } else {
        Err(io::Error::other(
            "path no longer identifies the file owned by this operation",
        ))
    }
}

/// Atomically creates a hard-link target for one live Linux file descriptor.
///
/// The source is resolved through `/proc/self/fd` with
/// `AT_SYMLINK_FOLLOW`. Linux 9p is rejected because its descriptor-link
/// behavior cannot be proven identity-safe. The target is create-only by
/// `linkat(2)` definition.
///
/// # Errors
///
/// Returns `InvalidInput` for an embedded NUL, `Unsupported` for Linux 9p or
/// an unverified procfs view, and the direct `linkat(2)` error otherwise.
pub(crate) fn hard_link_descriptor_create_new(
    source: BorrowedFd<'_>,
    target: &Path,
) -> io::Result<()> {
    let source_path = CString::new(format!("/proc/self/fd/{}", source.as_raw_fd()))
        .map_err(|_| io::Error::from(io::ErrorKind::InvalidInput))?;
    let target = CString::new(target.as_os_str().as_bytes())
        .map_err(|_| io::Error::from(io::ErrorKind::InvalidInput))?;
    reject_unsafe_file_system(source)?;
    let descriptor_metadata = File::from(source.try_clone_to_owned()?).metadata()?;
    let procfs_path = Path::new(std::ffi::OsStr::from_bytes(source_path.as_bytes()));
    let resolved_path = fs::read_link(procfs_path)?;
    let resolved_metadata = match fs::symlink_metadata(resolved_path) {
        Ok(metadata) => metadata,
        Err(source) if source.kind() == io::ErrorKind::NotFound => {
            return Err(io::Error::new(
                io::ErrorKind::Unsupported,
                "procfs descriptor path has no live namespace entry",
            ));
        }
        Err(source) => return Err(source),
    };
    if FileIdentity::from_metadata(&descriptor_metadata)
        != FileIdentity::from_metadata(&resolved_metadata)
    {
        return Err(io::Error::new(
            io::ErrorKind::Unsupported,
            "procfs descriptor path does not identify the borrowed file",
        ));
    }
    // SAFETY: both arguments are live NUL-terminated strings for this call;
    // `linkat` retains neither pointer.
    let result = unsafe {
        linkat(
            AT_CURRENT_WORKING_DIRECTORY,
            source_path.as_ptr(),
            AT_CURRENT_WORKING_DIRECTORY,
            target.as_ptr(),
            AT_SYMLINK_FOLLOW,
        )
    };
    if result == 0 {
        Ok(())
    } else {
        Err(io::Error::last_os_error())
    }
}

fn reject_unsafe_file_system(source: BorrowedFd<'_>) -> io::Result<()> {
    let mut statistics = [0 as c_long; STATFS_STORAGE_WORDS];
    // SAFETY: the initialized buffer is aligned for Linux `c_long` and is at
    // least as large as the supported Linux `struct statfs` layouts.
    let result = unsafe { fstatfs(source.as_raw_fd(), statistics.as_mut_ptr().cast::<c_void>()) };
    if result != 0 {
        return Err(io::Error::last_os_error());
    }
    if statistics[0] == V9FS_SUPER_MAGIC {
        return Err(io::Error::new(
            io::ErrorKind::Unsupported,
            "Linux 9p does not provide audited descriptor-link semantics",
        ));
    }
    Ok(())
}
