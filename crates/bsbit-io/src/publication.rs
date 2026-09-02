//! Staging, completion, create-only or replacement publication, and rollback.

use std::fs::{self, File};
use std::io;
use std::os::fd::AsFd;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use crate::file::{
    FileIdentity, absolute_path, create_new, hard_link_descriptor_create_new,
    remove_if_identity_matches, validate_absent, validate_create_target, validate_replace_target,
};

static STAGING_SEQUENCE: AtomicU64 = AtomicU64::new(0);

/// Selects an unused private sibling staging path for a future create-only
/// operation that cannot accept an already-open descriptor.
///
/// The returned path is absolute, belongs to the target's directory, and is
/// absent at the instant it is selected. It is deliberately not reserved;
/// the receiving format or index implementation must still use exclusive
/// creation and treat a concurrent winner as an ordinary `AlreadyExists`
/// failure. Code that can write through an owned descriptor should use
/// [`StagedFile::create_sibling`] instead.
///
/// `label` is sanitized to ASCII alphanumeric, dash, and underscore bytes.
///
/// # Errors
///
/// Returns target validation, parent-directory, or staging-path inspection
/// failures, including `AlreadyExists` after 64 occupied candidates.
pub fn select_sibling_staging_path(target: &Path, label: &str) -> io::Result<PathBuf> {
    let target = absolute_path(target)?;
    if target.file_name().is_none() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "create-only target does not name a file",
        ));
    }
    validate_create_target(&target)?;
    let parent = target.parent().ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            "create-only target has no parent directory",
        )
    })?;
    let label = sanitized_staging_label(label);
    for _ in 0..64 {
        let staging = sibling_staging_candidate(parent, &label);
        match validate_absent(&staging) {
            Ok(()) => return Ok(staging),
            Err(source) if source.kind() == io::ErrorKind::AlreadyExists => {}
            Err(source) => return Err(source),
        }
    }
    Err(io::Error::new(
        io::ErrorKind::AlreadyExists,
        "could not select an unused sibling staging path",
    ))
}

/// Selects an unused private sibling staging path for a future replacement.
///
/// Unlike [`select_sibling_staging_path`], an existing regular-file or
/// symbolic-link target is accepted. The returned path is not reserved; the
/// receiving writer must still create it exclusively.
///
/// # Errors
///
/// Returns target validation, parent-directory, or staging-path inspection
/// failures, including `AlreadyExists` after 64 occupied candidates.
pub fn select_sibling_staging_path_replace(target: &Path, label: &str) -> io::Result<PathBuf> {
    let target = absolute_path(target)?;
    if target.file_name().is_none() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "replacement target does not name a file",
        ));
    }
    validate_replace_target(&target)?;
    let parent = target.parent().ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            "replacement target has no parent directory",
        )
    })?;
    let label = sanitized_staging_label(label);
    for _ in 0..64 {
        let staging = sibling_staging_candidate(parent, &label);
        match validate_absent(&staging) {
            Ok(()) => return Ok(staging),
            Err(source) if source.kind() == io::ErrorKind::AlreadyExists => {}
            Err(source) => return Err(source),
        }
    }
    Err(io::Error::new(
        io::ErrorKind::AlreadyExists,
        "could not select an unused sibling staging path",
    ))
}

fn sanitized_staging_label(label: &str) -> String {
    let label: String = label
        .chars()
        .filter(|character| character.is_ascii_alphanumeric() || matches!(character, '-' | '_'))
        .collect();
    if label.is_empty() {
        String::from("output")
    } else {
        label
    }
}

fn sibling_staging_candidate(parent: &Path, label: &str) -> PathBuf {
    let sequence = STAGING_SEQUENCE.fetch_add(1, Ordering::Relaxed);
    parent.join(format!(
        ".bsbit-{label}-{:08x}-{sequence:016x}.tmp",
        std::process::id()
    ))
}

/// Lifecycle phase associated with a generic publication failure.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum PublicationPhase {
    /// Validate target/staging path policy.
    ValidatePaths,
    /// Reserve the private create-only staging path.
    CreateStaging,
    /// Verify descriptor and namespace identity.
    ValidateStaging,
    /// Synchronize complete bytes.
    Sync,
    /// Create the final target without replacement.
    Publish,
    /// Remove an owned staging or published path.
    Cleanup,
    /// Remove a just-published target during transactional rollback.
    Rollback,
}

/// Error from the generic file-publication lifecycle.
#[derive(Debug)]
pub struct PublicationError {
    phase: PublicationPhase,
    path: PathBuf,
    source: io::Error,
    cleanup_warning: Option<io::ErrorKind>,
}

impl PublicationError {
    fn new(phase: PublicationPhase, path: &Path, source: io::Error) -> Self {
        Self {
            phase,
            path: path.to_path_buf(),
            source,
            cleanup_warning: None,
        }
    }

    fn with_cleanup_warning(mut self, cleanup_warning: Option<io::ErrorKind>) -> Self {
        self.cleanup_warning = cleanup_warning;
        self
    }

    /// Returns the state transition that failed.
    #[must_use]
    pub const fn phase(&self) -> PublicationPhase {
        self.phase
    }

    /// Returns the path involved in the failed transition.
    #[must_use]
    pub fn path(&self) -> &Path {
        &self.path
    }

    /// Returns the underlying filesystem error class.
    #[must_use]
    pub fn kind(&self) -> io::ErrorKind {
        self.source.kind()
    }

    /// Returns a secondary failure from identity-safe staging cleanup.
    ///
    /// The primary phase and error kind remain authoritative. A warning is
    /// present only when the lifecycle also failed to remove the staging name
    /// that it had created exclusively.
    #[must_use]
    pub const fn cleanup_warning(&self) -> Option<io::ErrorKind> {
        self.cleanup_warning
    }

    /// Consumes this value and returns its underlying filesystem error.
    #[must_use]
    pub fn into_io_error(self) -> io::Error {
        self.source
    }
}

impl core::fmt::Display for PublicationError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            formatter,
            "file publication failed in {:?} for {}: {}",
            self.phase,
            self.path.display(),
            self.source
        )
    }
}

impl std::error::Error for PublicationError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        Some(&self.source)
    }
}

/// One exclusively created private file before format-specific finalization.
#[derive(Debug)]
pub struct StagedFile {
    target: Option<PathBuf>,
    staging: PathBuf,
    identity: FileIdentity,
    file: Option<File>,
    owns_staging: bool,
}

impl StagedFile {
    /// Reserves an explicit create-only staging path.
    ///
    /// # Errors
    ///
    /// Returns path absolutization, creation, or metadata failures.
    pub fn create_new(staging: impl AsRef<Path>) -> Result<Self, PublicationError> {
        let staging = absolute_path(staging.as_ref()).map_err(|source| {
            PublicationError::new(PublicationPhase::ValidatePaths, staging.as_ref(), source)
        })?;
        let (file, identity) = create_new(&staging).map_err(|source| {
            PublicationError::new(PublicationPhase::CreateStaging, &staging, source)
        })?;
        Ok(Self {
            target: None,
            staging,
            identity,
            file: Some(file),
            owns_staging: true,
        })
    }

    /// Reserves a unique private sibling of an absent final target.
    ///
    /// `label` is sanitized to ASCII alphanumeric, dash, and underscore bytes
    /// before it becomes part of the hidden staging name.
    ///
    /// # Errors
    ///
    /// Returns target validation or staging creation failures.
    pub fn create_sibling(target: impl AsRef<Path>, label: &str) -> Result<Self, PublicationError> {
        Self::create_sibling_with_policy(target.as_ref(), label, false)
    }

    /// Reserves a unique private sibling of a missing or replaceable target.
    ///
    /// Existing regular files and symbolic links are accepted; directories
    /// and other special filesystem objects are rejected.
    ///
    /// # Errors
    ///
    /// Returns target validation or staging creation failures.
    pub fn create_sibling_replace(
        target: impl AsRef<Path>,
        label: &str,
    ) -> Result<Self, PublicationError> {
        Self::create_sibling_with_policy(target.as_ref(), label, true)
    }

    fn create_sibling_with_policy(
        target: &Path,
        label: &str,
        replace: bool,
    ) -> Result<Self, PublicationError> {
        let target = absolute_path(target).map_err(|source| {
            PublicationError::new(PublicationPhase::ValidatePaths, target, source)
        })?;
        let validation = if replace {
            validate_replace_target(&target)
        } else {
            validate_absent(&target)
        };
        validation.map_err(|source| {
            PublicationError::new(PublicationPhase::ValidatePaths, &target, source)
        })?;
        let parent = target.parent().ok_or_else(|| {
            PublicationError::new(
                PublicationPhase::ValidatePaths,
                &target,
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "target has no parent directory",
                ),
            )
        })?;
        let label = sanitized_staging_label(label);
        for _ in 0..64 {
            let staging = sibling_staging_candidate(parent, &label);
            match create_new(&staging) {
                Ok((file, identity)) => {
                    return Ok(Self {
                        target: Some(target),
                        staging,
                        identity,
                        file: Some(file),
                        owns_staging: true,
                    });
                }
                Err(source) if source.kind() == io::ErrorKind::AlreadyExists => {}
                Err(source) => {
                    return Err(PublicationError::new(
                        PublicationPhase::CreateStaging,
                        &staging,
                        source,
                    ));
                }
            }
        }
        Err(PublicationError::new(
            PublicationPhase::CreateStaging,
            &target,
            io::Error::new(
                io::ErrorKind::AlreadyExists,
                "could not reserve a unique sibling staging path",
            ),
        ))
    }

    /// Returns the private staging path.
    #[must_use]
    pub fn path(&self) -> &Path {
        &self.staging
    }

    /// Borrows the initially reserved descriptor.
    #[must_use]
    pub fn file(&self) -> Option<&File> {
        self.file.as_ref()
    }

    /// Transfers the initially reserved descriptor to a format-specific
    /// encoder while this value retains namespace cleanup ownership.
    ///
    /// # Errors
    ///
    /// Returns `BrokenPipe` when the descriptor was already taken.
    pub fn take_file(&mut self) -> Result<File, PublicationError> {
        self.file.take().ok_or_else(|| {
            PublicationError::new(
                PublicationPhase::ValidateStaging,
                &self.staging,
                io::Error::new(
                    io::ErrorKind::BrokenPipe,
                    "staging descriptor already taken",
                ),
            )
        })
    }

    /// Transfers a finalized descriptor into the completed state.
    ///
    /// # Errors
    ///
    /// Returns an identity error if the descriptor or namespace path no longer
    /// identifies the exclusively created staging file.
    pub fn complete(mut self, file: File) -> Result<CompletedFile, PublicationError> {
        let descriptor_identity = match FileIdentity::from_file(&file) {
            Ok(identity) => identity,
            Err(source) => {
                let error =
                    PublicationError::new(PublicationPhase::ValidateStaging, &self.staging, source);
                return Err(self.error_with_cleanup(error));
            }
        };
        if descriptor_identity != self.identity {
            let error = PublicationError::new(
                PublicationPhase::ValidateStaging,
                &self.staging,
                io::Error::other("completed descriptor does not identify reserved staging file"),
            );
            return Err(self.error_with_cleanup(error));
        }
        let staging_matches = match self.identity.matches_path(&self.staging) {
            Ok(matches) => matches,
            Err(source) => {
                let error =
                    PublicationError::new(PublicationPhase::ValidateStaging, &self.staging, source);
                return Err(self.error_with_cleanup(error));
            }
        };
        if !staging_matches {
            self.owns_staging = false;
            return Err(PublicationError::new(
                PublicationPhase::ValidateStaging,
                &self.staging,
                io::Error::other("staging identity changed"),
            ));
        }
        self.file.take();
        self.owns_staging = false;
        Ok(CompletedFile {
            target: self.target.take(),
            staging: self.staging.clone(),
            identity: self.identity,
            file,
            owns_staging: true,
        })
    }

    /// Explicitly removes the owned staging path.
    ///
    /// # Errors
    ///
    /// Returns identity or removal failures without deleting a replacement.
    pub fn abort(mut self) -> Result<(), PublicationError> {
        self.file.take();
        self.cleanup(PublicationPhase::Cleanup)
    }

    /// Removes the owned staging path using a caller-supplied filesystem
    /// operation.
    ///
    /// This is an operation seam for deterministic filesystem-failure tests;
    /// ordinary callers should use [`Self::abort`]. Identity ownership and
    /// state transitions remain enforced by this type.
    ///
    /// # Errors
    ///
    /// Returns identity or removal failures without deleting a replacement.
    #[doc(hidden)]
    pub fn abort_with<Remove>(mut self, mut remove: Remove) -> Result<(), PublicationError>
    where
        Remove: FnMut(&Path, FileIdentity) -> io::Result<()>,
    {
        self.file.take();
        self.cleanup_with(PublicationPhase::Cleanup, &mut remove)
    }

    fn cleanup(&mut self, phase: PublicationPhase) -> Result<(), PublicationError> {
        self.cleanup_with(phase, &mut remove_if_identity_matches)
    }

    fn error_with_cleanup(&mut self, error: PublicationError) -> PublicationError {
        let cleanup_warning = self
            .cleanup(PublicationPhase::Cleanup)
            .err()
            .map(|cleanup| cleanup.kind());
        error.with_cleanup_warning(cleanup_warning)
    }

    fn cleanup_with<Remove>(
        &mut self,
        phase: PublicationPhase,
        remove: &mut Remove,
    ) -> Result<(), PublicationError>
    where
        Remove: FnMut(&Path, FileIdentity) -> io::Result<()>,
    {
        if !self.owns_staging {
            return Ok(());
        }
        self.owns_staging = false;
        remove(&self.staging, self.identity)
            .map_err(|source| PublicationError::new(phase, &self.staging, source))
    }
}

impl Drop for StagedFile {
    fn drop(&mut self) {
        self.file.take();
        let _ = self.cleanup(PublicationPhase::Cleanup);
    }
}

/// One format-complete file that has not yet been published.
#[derive(Debug)]
pub struct CompletedFile {
    target: Option<PathBuf>,
    staging: PathBuf,
    identity: FileIdentity,
    file: File,
    owns_staging: bool,
}

impl CompletedFile {
    /// Returns the completed private path.
    #[must_use]
    pub fn staging_path(&self) -> &Path {
        &self.staging
    }

    /// Borrows the completed descriptor for format-specific validation.
    ///
    /// The descriptor remains owned by this lifecycle state and therefore
    /// cannot outlive identity-safe staging cleanup.
    #[must_use]
    pub const fn file(&self) -> &File {
        &self.file
    }

    /// Verifies that the private path still names the completed descriptor.
    ///
    /// # Errors
    ///
    /// Returns namespace metadata failures. A missing or replaced path is
    /// reported as `Ok(false)` and is never removed by this check.
    pub fn staging_identity_matches(&self) -> Result<bool, PublicationError> {
        self.identity.matches_path(&self.staging).map_err(|source| {
            PublicationError::new(PublicationPhase::ValidateStaging, &self.staging, source)
        })
    }

    /// Relinquishes cleanup ownership and returns the still-verified staging
    /// path to a caller that assumes responsibility for it.
    ///
    /// # Errors
    ///
    /// Returns an identity error when the staging namespace entry was removed
    /// or replaced.
    pub fn into_path(mut self) -> Result<PathBuf, PublicationError> {
        if !self
            .identity
            .matches_path(&self.staging)
            .map_err(|source| {
                PublicationError::new(PublicationPhase::ValidateStaging, &self.staging, source)
            })?
        {
            self.owns_staging = false;
            return Err(PublicationError::new(
                PublicationPhase::ValidateStaging,
                &self.staging,
                io::Error::other("staging identity changed"),
            ));
        }
        self.owns_staging = false;
        Ok(self.staging.clone())
    }

    /// Publishes at the target recorded by [`StagedFile::create_sibling`].
    ///
    /// # Errors
    ///
    /// Returns synchronization, identity, target-race, or descriptor-link
    /// failures. Publication never replaces an existing target. If an error
    /// follows a successful link, the implementation attempts rollback only
    /// while the target still identifies the completed descriptor; a foreign
    /// replacement is preserved.
    pub fn publish_create_new(self) -> Result<PublishedFile, PublicationError> {
        let target = self.target.clone().ok_or_else(|| {
            PublicationError::new(
                PublicationPhase::ValidatePaths,
                &self.staging,
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "no publication target was recorded",
                ),
            )
        })?;
        self.publish_to(target)
    }

    /// Publishes an explicitly staged file at a distinct sibling target.
    ///
    /// # Errors
    ///
    /// Returns path-policy, synchronization, identity, or link failures.
    pub fn publish_create_new_at(
        self,
        target: impl AsRef<Path>,
    ) -> Result<PublishedFile, PublicationError> {
        let target = absolute_path(target.as_ref()).map_err(|source| {
            PublicationError::new(PublicationPhase::ValidatePaths, target.as_ref(), source)
        })?;
        self.publish_to(target)
    }

    /// Atomically publishes at the target recorded by
    /// [`StagedFile::create_sibling_replace`].
    ///
    /// An existing regular file or symbolic link is replaced only after the
    /// completed staging descriptor is synchronized and revalidated. The old
    /// target remains available for [`PublishedFile::rollback`] until the
    /// returned publication is dropped.
    ///
    /// # Errors
    ///
    /// Returns path-policy, synchronization, identity, backup, or rename
    /// failures. A failure before replacement leaves the old target intact.
    pub fn publish_replace(self) -> Result<PublishedFile, PublicationError> {
        let target = self.target.clone().ok_or_else(|| {
            PublicationError::new(
                PublicationPhase::ValidatePaths,
                &self.staging,
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "no publication target was recorded",
                ),
            )
        })?;
        self.replace_to(target)
    }

    /// Atomically publishes an explicitly staged file at a missing or
    /// replaceable sibling target.
    ///
    /// # Errors
    ///
    /// Returns the same lifecycle failures as [`Self::publish_replace`].
    pub fn publish_replace_at(
        self,
        target: impl AsRef<Path>,
    ) -> Result<PublishedFile, PublicationError> {
        let target = absolute_path(target.as_ref()).map_err(|source| {
            PublicationError::new(PublicationPhase::ValidatePaths, target.as_ref(), source)
        })?;
        self.replace_to(target)
    }

    /// Publishes with caller-supplied synchronization, linking, and cleanup
    /// operations while retaining this type's path and identity policy.
    ///
    /// This is an operation seam for deterministic filesystem-failure tests;
    /// ordinary callers should use [`Self::publish_create_new_at`].
    ///
    /// # Errors
    ///
    /// Returns the same lifecycle failures as
    /// [`Self::publish_create_new_at`], including any secondary cleanup
    /// warning.
    #[doc(hidden)]
    pub fn publish_create_new_at_with<Sync, Link, Remove>(
        self,
        target: impl AsRef<Path>,
        sync: Sync,
        link: Link,
        remove: Remove,
    ) -> Result<PublishedFile, PublicationError>
    where
        Sync: FnMut(&File) -> io::Result<()>,
        Link: FnMut(&File, &Path) -> io::Result<()>,
        Remove: FnMut(&Path, FileIdentity) -> io::Result<()>,
    {
        let target = absolute_path(target.as_ref()).map_err(|source| {
            PublicationError::new(PublicationPhase::ValidatePaths, target.as_ref(), source)
        })?;
        self.publish_to_with(target, sync, link, remove)
    }

    fn publish_to(self, target: PathBuf) -> Result<PublishedFile, PublicationError> {
        let staging = self.staging.clone();
        let identity = self.identity;
        self.publish_to_with(
            target,
            File::sync_all,
            move |file, target| match hard_link_descriptor_create_new(file.as_fd(), target) {
                Ok(()) => Ok(()),
                Err(source) if source.kind() == io::ErrorKind::Unsupported => {
                    // Some mounted filesystems cannot resolve a procfs
                    // descriptor for linkat. The path fallback remains safe
                    // because both the source and newly created target are
                    // checked against the held descriptor identity.
                    if !identity.matches_path(&staging)? {
                        return Err(io::Error::other(
                            "staging identity changed before path-based publication",
                        ));
                    }
                    fs::hard_link(&staging, target)
                }
                Err(source) => Err(source),
            },
            remove_if_identity_matches,
        )
    }

    fn replace_to(mut self, target: PathBuf) -> Result<PublishedFile, PublicationError> {
        if let Err(source) = validate_sibling_publication_paths(&target, &self.staging) {
            let error = PublicationError::new(PublicationPhase::ValidatePaths, &target, source);
            return Err(self.error_with_cleanup(error, &mut remove_if_identity_matches));
        }
        if let Err(source) = validate_replace_target(&target) {
            let error = PublicationError::new(PublicationPhase::ValidatePaths, &target, source);
            return Err(self.error_with_cleanup(error, &mut remove_if_identity_matches));
        }
        if let Err(source) = self.file.sync_all() {
            let error = PublicationError::new(PublicationPhase::Sync, &self.staging, source);
            return Err(self.error_with_cleanup(error, &mut remove_if_identity_matches));
        }
        let staging_matches = match self.identity.matches_path(&self.staging) {
            Ok(matches) => matches,
            Err(source) => {
                let error =
                    PublicationError::new(PublicationPhase::ValidateStaging, &self.staging, source);
                return Err(self.error_with_cleanup(error, &mut remove_if_identity_matches));
            }
        };
        if !staging_matches {
            self.owns_staging = false;
            return Err(PublicationError::new(
                PublicationPhase::ValidateStaging,
                &self.staging,
                io::Error::other("staging identity changed"),
            ));
        }
        let backup = match ReplacementBackup::capture(&target) {
            Ok(backup) => backup,
            Err(source) => {
                let error = PublicationError::new(PublicationPhase::Publish, &target, source);
                return Err(self.error_with_cleanup(error, &mut remove_if_identity_matches));
            }
        };
        if let Err(source) = fs::rename(&self.staging, &target) {
            let cleanup_warning = backup
                .and_then(|backup| backup.remove().err())
                .map(|warning| warning.kind());
            let error = PublicationError::new(PublicationPhase::Publish, &target, source)
                .with_cleanup_warning(cleanup_warning);
            return Err(self.error_with_cleanup(error, &mut remove_if_identity_matches));
        }
        self.owns_staging = false;
        let target_matches = self
            .identity
            .matches_path(&target)
            .map_err(|source| PublicationError::new(PublicationPhase::Publish, &target, source))?;
        if !target_matches {
            return Err(PublicationError::new(
                PublicationPhase::Publish,
                &target,
                io::Error::other("replaced target does not identify the completed descriptor"),
            ));
        }
        Ok(PublishedFile {
            target,
            staging: self.staging.clone(),
            identity: self.identity,
            cleanup_warning: None,
            replacement_backup: backup,
        })
    }

    fn publish_to_with<Sync, Link, Remove>(
        mut self,
        target: PathBuf,
        mut sync: Sync,
        mut link: Link,
        mut remove: Remove,
    ) -> Result<PublishedFile, PublicationError>
    where
        Sync: FnMut(&File) -> io::Result<()>,
        Link: FnMut(&File, &Path) -> io::Result<()>,
        Remove: FnMut(&Path, FileIdentity) -> io::Result<()>,
    {
        if let Err(source) = validate_sibling_publication_paths(&target, &self.staging) {
            let error = PublicationError::new(PublicationPhase::ValidatePaths, &target, source);
            return Err(self.error_with_cleanup(error, &mut remove));
        }
        if let Err(source) = validate_absent(&target) {
            let error = PublicationError::new(PublicationPhase::ValidatePaths, &target, source);
            return Err(self.error_with_cleanup(error, &mut remove));
        }
        if let Err(source) = sync(&self.file) {
            let error = PublicationError::new(PublicationPhase::Sync, &self.staging, source);
            return Err(self.error_with_cleanup(error, &mut remove));
        }
        let staging_matches = match self.identity.matches_path(&self.staging) {
            Ok(matches) => matches,
            Err(source) => {
                let error =
                    PublicationError::new(PublicationPhase::ValidateStaging, &self.staging, source);
                return Err(self.error_with_cleanup(error, &mut remove));
            }
        };
        if !staging_matches {
            self.owns_staging = false;
            return Err(PublicationError::new(
                PublicationPhase::ValidateStaging,
                &self.staging,
                io::Error::other("staging identity changed"),
            ));
        }
        if let Err(source) = link(&self.file, &target) {
            let error = PublicationError::new(PublicationPhase::Publish, &target, source);
            return Err(self.error_with_cleanup(error, &mut remove));
        }
        let target_matches = match self.identity.matches_path(&target) {
            Ok(matches) => matches,
            Err(source) => {
                // The link operation already reported success, so try to
                // retract only the object represented by our held descriptor.
                // Never derive deletion authority from a fresh observation of
                // a namespace entry that another publisher may have replaced.
                let _ = remove(&target, self.identity);
                let error = PublicationError::new(PublicationPhase::Publish, &target, source);
                return Err(self.error_with_cleanup(error, &mut remove));
            }
        };
        if !target_matches {
            // A concurrent actor may have displaced our link and installed a
            // replacement before this validation. Conditional rollback must
            // use the identity owned before publication, never the identity of
            // whatever happens to occupy `target` now.
            let _ = remove(&target, self.identity);
            let error = PublicationError::new(
                PublicationPhase::Publish,
                &target,
                io::Error::other("published target does not identify the completed descriptor"),
            );
            return Err(self.error_with_cleanup(error, &mut remove));
        }
        let cleanup_warning = remove(&self.staging, self.identity)
            .err()
            .map(|source| source.kind());
        self.owns_staging = false;
        Ok(PublishedFile {
            target,
            staging: self.staging.clone(),
            identity: self.identity,
            cleanup_warning,
            replacement_backup: None,
        })
    }

    /// Removes the completed private file without publishing it.
    ///
    /// # Errors
    ///
    /// Returns identity or removal failures without deleting a replacement.
    pub fn remove(mut self) -> Result<(), PublicationError> {
        self.cleanup()
    }

    fn cleanup(&mut self) -> Result<(), PublicationError> {
        self.cleanup_with(PublicationPhase::Cleanup, &mut remove_if_identity_matches)
    }

    fn error_with_cleanup<Remove>(
        &mut self,
        error: PublicationError,
        remove: &mut Remove,
    ) -> PublicationError
    where
        Remove: FnMut(&Path, FileIdentity) -> io::Result<()>,
    {
        let cleanup_warning = self
            .cleanup_with(PublicationPhase::Cleanup, remove)
            .err()
            .map(|cleanup| cleanup.kind());
        error.with_cleanup_warning(cleanup_warning)
    }

    fn cleanup_with<Remove>(
        &mut self,
        phase: PublicationPhase,
        remove: &mut Remove,
    ) -> Result<(), PublicationError>
    where
        Remove: FnMut(&Path, FileIdentity) -> io::Result<()>,
    {
        if !self.owns_staging {
            return Ok(());
        }
        self.owns_staging = false;
        remove(&self.staging, self.identity)
            .map_err(|source| PublicationError::new(phase, &self.staging, source))
    }
}

/// Verifies the shared path policy for an explicitly staged publication.
///
/// The paths are compared after lexical absolutization. They must be distinct
/// and have the same parent so that publication and cleanup remain confined to
/// one directory. Target absence is deliberately checked at the create-only
/// publication transition, where a concurrent winner is reported.
///
/// # Errors
///
/// Returns absolutization or `InvalidInput` path-policy failures.
pub fn validate_sibling_publication_paths(target: &Path, staging: &Path) -> io::Result<()> {
    let target = absolute_path(target)?;
    let staging = absolute_path(staging)?;
    if target == staging || target.parent() != staging.parent() {
        Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "target and staging must be distinct siblings",
        ))
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use std::fs;
    use std::io::{self, Write};
    use std::os::fd::AsFd;
    use std::path::{Path, PathBuf};
    use std::time::{SystemTime, UNIX_EPOCH};

    use super::{CompletedFile, PublicationPhase, StagedFile};
    use crate::{hard_link_descriptor_create_new, remove_if_identity_matches};

    fn unique_directory(label: &str) -> PathBuf {
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock")
            .as_nanos();
        std::env::temp_dir().join(format!(
            "bsbit-io-publication-unit-{label}-{}-{nonce}",
            std::process::id()
        ))
    }

    fn completed(target: &Path, bytes: &[u8]) -> CompletedFile {
        let mut staged = StagedFile::create_sibling(target, "fault-contract").expect("stage");
        let mut file = staged.take_file().expect("descriptor");
        file.write_all(bytes).expect("write complete bytes");
        staged.complete(file).expect("complete")
    }

    #[test]
    fn injected_sync_link_and_cleanup_failures_have_exact_publication_semantics() {
        let directory = unique_directory("faults");
        fs::create_dir(&directory).expect("directory");

        let sync_target = directory.join("sync-target");
        let sync_staging = {
            let completed = completed(&sync_target, b"sync bytes");
            let staging = completed.staging_path().to_path_buf();
            let error = completed
                .publish_to_with(
                    sync_target.clone(),
                    |_| Err(io::Error::from(io::ErrorKind::WriteZero)),
                    |_, _| panic!("link must not follow failed sync"),
                    remove_if_identity_matches,
                )
                .expect_err("sync failure");
            assert_eq!(error.phase(), PublicationPhase::Sync);
            assert_eq!(error.kind(), io::ErrorKind::WriteZero);
            assert_eq!(error.cleanup_warning(), None);
            staging
        };
        assert!(!sync_target.exists());
        assert!(!sync_staging.exists());

        let link_target = directory.join("link-target");
        let link_staging = {
            let completed = completed(&link_target, b"link bytes");
            let staging = completed.staging_path().to_path_buf();
            let error = completed
                .publish_to_with(
                    link_target.clone(),
                    |_| Ok(()),
                    |_, _| Err(io::Error::from(io::ErrorKind::PermissionDenied)),
                    remove_if_identity_matches,
                )
                .expect_err("link failure");
            assert_eq!(error.phase(), PublicationPhase::Publish);
            assert_eq!(error.kind(), io::ErrorKind::PermissionDenied);
            assert_eq!(error.cleanup_warning(), None);
            staging
        };
        assert!(!link_target.exists());
        assert!(!link_staging.exists());

        let cleanup_target = directory.join("cleanup-target");
        let completed = completed(&cleanup_target, b"cleanup bytes");
        let cleanup_staging = completed.staging_path().to_path_buf();
        let published = match completed.publish_to_with(
            cleanup_target.clone(),
            |_| Ok(()),
            |file, target| hard_link_descriptor_create_new(file.as_fd(), target),
            |_, _| Err(io::Error::from(io::ErrorKind::PermissionDenied)),
        ) {
            Ok(published) => published,
            Err(error) if error.kind() == io::ErrorKind::Unsupported => {
                fs::remove_dir_all(directory).expect("unsupported cleanup");
                return;
            }
            Err(error) => panic!("publish before injected cleanup failure: {error}"),
        };
        assert_eq!(
            published.cleanup_warning(),
            Some(io::ErrorKind::PermissionDenied)
        );
        assert_eq!(fs::read(&cleanup_target).expect("target"), b"cleanup bytes");
        assert_eq!(
            fs::read(&cleanup_staging).expect("retained staging"),
            b"cleanup bytes"
        );
        published
            .rollback()
            .expect("rollback target and retained owned staging");
        assert!(!cleanup_target.exists());
        assert!(!cleanup_staging.exists());
        fs::remove_dir(directory).expect("directory cleanup");
    }

    #[test]
    fn descriptor_link_never_publishes_a_post_validation_replacement() {
        let directory = unique_directory("descriptor-link");
        fs::create_dir(&directory).expect("directory");
        let target = directory.join("target");
        let completed = completed(&target, b"owned bytes");
        let staging = completed.staging_path().to_path_buf();
        let displaced = directory.join("displaced");
        let replacement = b"replacement bytes";

        let published = match completed.publish_to_with(
            target.clone(),
            |_| Ok(()),
            |file, target| {
                fs::rename(&staging, &displaced)?;
                fs::write(&staging, replacement)?;
                hard_link_descriptor_create_new(file.as_fd(), target)
            },
            remove_if_identity_matches,
        ) {
            Ok(published) => published,
            Err(error) if error.kind() == io::ErrorKind::Unsupported => {
                assert!(!target.exists());
                fs::remove_file(staging).expect("replacement cleanup");
                fs::remove_file(displaced).expect("displaced cleanup");
                fs::remove_dir(directory).expect("directory cleanup");
                return;
            }
            Err(error) => panic!("descriptor publication: {error}"),
        };

        assert_eq!(published.cleanup_warning(), Some(io::ErrorKind::Other));
        assert_eq!(fs::read(&target).expect("published bytes"), b"owned bytes");
        assert_eq!(fs::read(&staging).expect("replacement"), replacement);
        fs::remove_file(target).expect("target cleanup");
        fs::remove_file(staging).expect("replacement cleanup");
        fs::remove_file(displaced).expect("displaced cleanup");
        fs::remove_dir(directory).expect("directory cleanup");
    }

    #[test]
    fn post_link_target_replacement_is_preserved_during_failed_validation() {
        let directory = unique_directory("post-link-target-replacement");
        fs::create_dir(&directory).expect("directory");
        let target = directory.join("target");
        let displaced_target = directory.join("displaced-target");
        let completed = completed(&target, b"owned bytes");
        let staging = completed.staging_path().to_path_buf();
        let replacement = b"caller-owned replacement";

        let error = completed
            .publish_to_with(
                target.clone(),
                |_| Ok(()),
                |_, target| {
                    fs::hard_link(&staging, target)?;
                    fs::rename(target, &displaced_target)?;
                    fs::write(target, replacement)
                },
                remove_if_identity_matches,
            )
            .expect_err("post-link replacement prevents confirmed publication");

        assert_eq!(error.phase(), PublicationPhase::Publish);
        assert_eq!(error.kind(), io::ErrorKind::Other);
        assert_eq!(error.cleanup_warning(), None);
        assert_eq!(
            fs::read(&target).expect("replacement survives"),
            replacement
        );
        assert_eq!(
            fs::read(&displaced_target).expect("owned link survives displacement"),
            b"owned bytes"
        );
        assert!(!staging.exists());

        fs::remove_file(target).expect("replacement cleanup");
        fs::remove_file(displaced_target).expect("owned link cleanup");
        fs::remove_dir(directory).expect("directory cleanup");
    }

    #[test]
    fn descriptor_mismatch_removes_only_the_reserved_staging_path() {
        let directory = unique_directory("descriptor-mismatch");
        fs::create_dir(&directory).expect("directory");
        let staging = directory.join("staging");
        let foreign = directory.join("foreign");
        let mut staged = StagedFile::create_new(&staging).expect("stage");
        let reserved = staged.take_file().expect("reserved descriptor");
        let foreign_descriptor = fs::File::create(&foreign).expect("foreign descriptor");

        let error = staged
            .complete(foreign_descriptor)
            .expect_err("foreign descriptor cannot complete reserved staging");

        assert_eq!(error.phase(), PublicationPhase::ValidateStaging);
        assert_eq!(error.kind(), io::ErrorKind::Other);
        assert_eq!(error.cleanup_warning(), None);
        assert!(!staging.exists());
        assert!(foreign.exists());

        drop(reserved);
        fs::remove_file(foreign).expect("foreign cleanup");
        fs::remove_dir(directory).expect("directory cleanup");
    }

    #[test]
    fn staging_metadata_error_reports_identity_safe_cleanup_warning() {
        use std::os::unix::fs::symlink;

        let directory = unique_directory("staging-metadata-error");
        let displaced_directory = directory.with_extension("displaced");
        fs::create_dir(&directory).expect("directory");
        let staging = directory.join("staging");
        let mut staged = StagedFile::create_new(&staging).expect("stage");
        let descriptor = staged.take_file().expect("reserved descriptor");
        fs::rename(&directory, &displaced_directory).expect("displace parent directory");
        symlink(&directory, &directory).expect("self-referential parent symlink");

        let error = staged
            .complete(descriptor)
            .expect_err("parent symlink loop prevents staging metadata validation");

        assert_eq!(error.phase(), PublicationPhase::ValidateStaging);
        assert_eq!(error.cleanup_warning(), Some(error.kind()));
        assert!(displaced_directory.join("staging").exists());

        fs::remove_file(&directory).expect("symlink cleanup");
        fs::remove_file(displaced_directory.join("staging")).expect("staging cleanup");
        fs::remove_dir(displaced_directory).expect("directory cleanup");
    }
}

impl Drop for CompletedFile {
    fn drop(&mut self) {
        let _ = self.cleanup();
    }
}

#[derive(Debug, Eq, PartialEq)]
struct ReplacementBackup {
    path: PathBuf,
    identity: FileIdentity,
}

impl ReplacementBackup {
    fn capture(target: &Path) -> io::Result<Option<Self>> {
        let metadata = match fs::symlink_metadata(target) {
            Ok(metadata) => metadata,
            Err(source) if source.kind() == io::ErrorKind::NotFound => return Ok(None),
            Err(source) => return Err(source),
        };
        if !metadata.is_file() && !metadata.file_type().is_symlink() {
            return Err(io::Error::new(
                io::ErrorKind::Unsupported,
                "existing output target is not a regular file or symbolic link",
            ));
        }
        let identity = FileIdentity::from_metadata(&metadata);
        let parent = target.parent().ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                "replacement target has no parent directory",
            )
        })?;
        for _ in 0..64 {
            let path = sibling_staging_candidate(parent, "backup");
            match fs::hard_link(target, &path) {
                Ok(()) => {
                    let backup = Self { path, identity };
                    if backup.identity.matches_path(&backup.path)?
                        && backup.identity.matches_path(target)?
                    {
                        return Ok(Some(backup));
                    }
                    let _ = backup.remove();
                    return Err(io::Error::other(
                        "output target changed while its rollback backup was created",
                    ));
                }
                Err(source) if source.kind() == io::ErrorKind::AlreadyExists => {}
                Err(source) => return Err(source),
            }
        }
        Err(io::Error::new(
            io::ErrorKind::AlreadyExists,
            "could not reserve an unused replacement-backup path",
        ))
    }

    fn remove(self) -> io::Result<()> {
        remove_if_identity_matches(&self.path, self.identity)
    }

    fn restore(self, target: &Path) -> io::Result<()> {
        if !self.identity.matches_path(&self.path)? {
            return Err(io::Error::other("replacement backup identity changed"));
        }
        fs::hard_link(&self.path, target)?;
        if !self.identity.matches_path(target)? {
            let _ = remove_if_identity_matches(target, self.identity);
            return Err(io::Error::other(
                "restored target does not identify the replacement backup",
            ));
        }
        remove_if_identity_matches(&self.path, self.identity)
    }
}

/// One successfully published file with rollback authority.
#[derive(Debug, Eq, PartialEq)]
pub struct PublishedFile {
    target: PathBuf,
    staging: PathBuf,
    identity: FileIdentity,
    cleanup_warning: Option<io::ErrorKind>,
    replacement_backup: Option<ReplacementBackup>,
}

impl PublishedFile {
    /// Returns the final target path.
    #[must_use]
    pub fn target_path(&self) -> &Path {
        &self.target
    }

    /// Returns the private path used before publication.
    #[must_use]
    pub fn staging_path(&self) -> &Path {
        &self.staging
    }

    /// Returns a non-fatal post-publication staging cleanup warning.
    #[must_use]
    pub const fn cleanup_warning(&self) -> Option<io::ErrorKind> {
        self.cleanup_warning
    }

    /// Removes the target only while it still names the published object and
    /// retries cleanup of an owned staging link retained after publication.
    ///
    /// # Errors
    ///
    /// Returns identity or removal failures without deleting a replacement.
    pub fn rollback(mut self) -> Result<(), PublicationError> {
        let replacement_backup = self.replacement_backup.take();
        let target_result = remove_if_identity_matches(&self.target, self.identity);
        let staging_result = self.cleanup_warning.map_or(Ok(()), |_| {
            remove_if_identity_matches(&self.staging, self.identity)
        });
        target_result.map_err(|source| {
            PublicationError::new(PublicationPhase::Rollback, &self.target, source)
        })?;
        staging_result.map_err(|source| {
            PublicationError::new(PublicationPhase::Rollback, &self.staging, source)
        })?;
        replacement_backup.map_or(Ok(()), |backup| {
            let backup_path = backup.path.clone();
            backup.restore(&self.target).map_err(|source| {
                PublicationError::new(PublicationPhase::Rollback, &backup_path, source)
            })
        })
    }
}

impl Drop for PublishedFile {
    fn drop(&mut self) {
        if let Some(backup) = self.replacement_backup.take() {
            let _ = backup.remove();
        }
    }
}
