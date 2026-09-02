//! Streaming text/BGZF encoding over the generic file-publication lifecycle.

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

use bsbit_io::{CompletedFile, PublicationError, PublicationPhase, PublishedFile, StagedFile};

use super::{BgzfWriter, HtsError, HtsErrorKind, HtsOperation, io_error, simple_error};

const OUTPUT_BUFFER_BYTES: usize = 1 << 20;

/// Compression used by one streaming text output.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum TextOutputCompression {
    /// Write ordinary uncompressed bytes.
    Plain,
    /// Write block gzip through `HTSlib`, including its canonical EOF block.
    Bgzf,
}

enum TextBackend {
    Plain(BufWriter<File>),
    Bgzf(BufWriter<BgzfWriter>),
}

impl TextBackend {
    fn finish(self) -> io::Result<File> {
        match self {
            Self::Plain(mut writer) => {
                writer.flush()?;
                writer
                    .into_inner()
                    .map_err(std::io::IntoInnerError::into_error)
            }
            Self::Bgzf(mut writer) => {
                writer.flush()?;
                let writer = writer
                    .into_inner()
                    .map_err(std::io::IntoInnerError::into_error)?;
                writer.finish()
            }
        }
    }

    fn writer(&mut self) -> &mut dyn Write {
        match self {
            Self::Plain(writer) => writer,
            Self::Bgzf(writer) => writer,
        }
    }
}

/// Streaming format encoder backed by one generic private staging file.
pub struct TextStagingWriter {
    staged: Option<StagedFile>,
    backend: Option<TextBackend>,
    terminal: bool,
}

impl TextStagingWriter {
    /// Creates a private sibling staging file beside an absent final target.
    ///
    /// `compression_threads == 0` performs synchronous BGZF compression.
    ///
    /// # Errors
    ///
    /// Returns target/staging publication errors or a BGZF-open error.
    pub fn create_sibling(
        target: impl AsRef<Path>,
        compression: TextOutputCompression,
        compression_threads: u32,
    ) -> Result<Self, HtsError> {
        Self::create_sibling_with_policy(target.as_ref(), compression, compression_threads, false)
    }

    /// Creates a private sibling staging file beside a missing or replaceable
    /// final target.
    ///
    /// # Errors
    ///
    /// Returns target/staging publication errors or a BGZF-open error.
    pub fn create_sibling_replace(
        target: impl AsRef<Path>,
        compression: TextOutputCompression,
        compression_threads: u32,
    ) -> Result<Self, HtsError> {
        Self::create_sibling_with_policy(target.as_ref(), compression, compression_threads, true)
    }

    fn create_sibling_with_policy(
        target: &Path,
        compression: TextOutputCompression,
        compression_threads: u32,
        replace: bool,
    ) -> Result<Self, HtsError> {
        let mut staged = if replace {
            StagedFile::create_sibling_replace(target, "text")
        } else {
            StagedFile::create_sibling(target, "text")
        }
        .map_err(map_publication_error)?;
        let path = staged.path().to_path_buf();
        let file = staged.take_file().map_err(map_publication_error)?;
        let backend = match compression {
            TextOutputCompression::Plain => {
                TextBackend::Plain(BufWriter::with_capacity(OUTPUT_BUFFER_BYTES, file))
            }
            TextOutputCompression::Bgzf => {
                let writer =
                    BgzfWriter::from_file(file, compression_threads).map_err(|source| {
                        io_error(HtsOperation::OpenTextOutput, &path, None, source)
                    })?;
                TextBackend::Bgzf(BufWriter::with_capacity(OUTPUT_BUFFER_BYTES, writer))
            }
        };
        Ok(Self {
            staged: Some(staged),
            backend: Some(backend),
            terminal: false,
        })
    }

    /// Returns the private staging path currently owned by this writer.
    ///
    /// # Panics
    ///
    /// Panics only if the writer's internal live-state invariant is violated.
    /// A caller cannot produce that state through the public API.
    #[must_use]
    pub fn staging_path(&self) -> &Path {
        self.staged
            .as_ref()
            .expect("live text writer retains staging state")
            .path()
    }

    /// Finalizes the encoding and transfers the file into the completed state.
    ///
    /// # Errors
    ///
    /// Returns a terminal, compression-finalization, or identity error.
    pub fn finish(mut self) -> Result<CompletedTextOutput, HtsError> {
        let path = self.staging_path().to_path_buf();
        if self.terminal {
            return Err(simple_error(
                HtsOperation::FinishTextOutput,
                &path,
                None,
                HtsErrorKind::Terminal,
            ));
        }
        let backend = self.backend.take().ok_or_else(|| {
            simple_error(
                HtsOperation::FinishTextOutput,
                &path,
                None,
                HtsErrorKind::Terminal,
            )
        })?;
        let file = backend
            .finish()
            .map_err(|source| io_error(HtsOperation::FinishTextOutput, &path, None, source))?;
        let staged = self.staged.take().ok_or_else(|| {
            simple_error(
                HtsOperation::FinishTextOutput,
                &path,
                None,
                HtsErrorKind::Terminal,
            )
        })?;
        let completed = staged.complete(file).map_err(map_publication_error)?;
        Ok(CompletedTextOutput { completed })
    }
}

impl Write for TextStagingWriter {
    fn write(&mut self, buffer: &[u8]) -> io::Result<usize> {
        if self.terminal {
            return Err(io::Error::other("text staging writer is terminal"));
        }
        let result = self
            .backend
            .as_mut()
            .ok_or_else(|| io::Error::other("text staging writer is already finished"))?
            .writer()
            .write(buffer);
        if result.is_err() {
            self.terminal = true;
        }
        result
    }

    fn flush(&mut self) -> io::Result<()> {
        if self.terminal {
            return Err(io::Error::other("text staging writer is terminal"));
        }
        let result = self
            .backend
            .as_mut()
            .ok_or_else(|| io::Error::other("text staging writer is already finished"))?
            .writer()
            .flush();
        if result.is_err() {
            self.terminal = true;
        }
        result
    }
}

impl Drop for TextStagingWriter {
    fn drop(&mut self) {
        self.backend.take();
        self.staged.take();
    }
}

/// One completely finalized text/BGZF staging file.
pub struct CompletedTextOutput {
    completed: CompletedFile,
}

impl CompletedTextOutput {
    /// Publishes the completed bytes at the originally requested absent target.
    ///
    /// # Errors
    ///
    /// Returns synchronization, identity, target-race, or link failures.
    pub fn publish_create_new(self) -> Result<TextPublication, HtsError> {
        self.completed
            .publish_create_new()
            .map(|published| TextPublication { published })
            .map_err(map_publication_error)
    }

    /// Atomically publishes the completed bytes, replacing an existing file.
    ///
    /// # Errors
    ///
    /// Returns synchronization, identity, backup, or rename failures.
    pub fn publish_replace(self) -> Result<TextPublication, HtsError> {
        self.completed
            .publish_replace()
            .map(|published| TextPublication { published })
            .map_err(map_publication_error)
    }
}

/// Details and rollback authority for one published text output.
#[derive(Debug, Eq, PartialEq)]
pub struct TextPublication {
    published: PublishedFile,
}

impl TextPublication {
    /// Returns the absolute final target path.
    #[must_use]
    pub fn target_path(&self) -> &Path {
        self.published.target_path()
    }

    /// Returns the absolute private staging path.
    #[must_use]
    pub fn staging_path(&self) -> &Path {
        self.published.staging_path()
    }

    /// Returns a non-fatal post-publication staging cleanup warning.
    #[must_use]
    pub fn cleanup_warning(&self) -> Option<HtsErrorKind> {
        self.published.cleanup_warning().map(HtsErrorKind::Io)
    }

    /// Removes the just-published target if it still has the published identity.
    ///
    /// # Errors
    ///
    /// Returns an identity or filesystem failure without removing a replacement.
    pub fn rollback(self) -> Result<(), HtsError> {
        self.published.rollback().map_err(map_publication_error)
    }
}

fn map_publication_error(error: PublicationError) -> HtsError {
    let operation = match error.phase() {
        PublicationPhase::ValidatePaths => HtsOperation::ValidatePublicationPaths,
        PublicationPhase::CreateStaging => HtsOperation::CreateStaging,
        PublicationPhase::ValidateStaging => HtsOperation::ValidateStaging,
        PublicationPhase::Sync => HtsOperation::SyncOutput,
        PublicationPhase::Publish => HtsOperation::PublishOutput,
        PublicationPhase::Cleanup => HtsOperation::Cleanup,
        PublicationPhase::Rollback => HtsOperation::RollbackOutput,
    };
    let path = error.path().to_path_buf();
    io_error(operation, &path, None, error.into_io_error())
}
