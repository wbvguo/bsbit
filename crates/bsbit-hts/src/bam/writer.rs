//! BAM staging, encoding, finalization, and publication.

use std::ffi::CString;
use std::io;
use std::os::fd::AsRawFd;
use std::path::{Path, PathBuf};

use bsbit_core::cigar::CoreCigarOp;
use bsbit_io::{CompletedFile, PublicationError, PublishedFile, StagedFile};

use crate::alignment_record::{
    AlignmentAuxiliaryMode, AlignmentCigarOp, AlignmentRecord, AlignmentRecordLimits,
    BorrowedAlignmentRecord,
};
use crate::htslib::{
    HtsError, HtsErrorKind, HtsOperation, absolute_path, encode_error, io_error, native_error,
    nul_error, simple_error,
};
use crate::sam::{SamHeader, sam_flag, sam_header_bytes, sam_record_bytes};
use crate::sys::{NativeBamRecordFields, NativeBamWriter};

fn reference_id_from_ordinal(ordinal: u64) -> Option<i32> {
    i32::try_from(ordinal).ok()
}

pub(crate) fn core_cigar_word(length: u64, operation: CoreCigarOp) -> Option<u32> {
    let operation = match operation {
        CoreCigarOp::M => AlignmentCigarOp::Match,
        CoreCigarOp::I => AlignmentCigarOp::Insertion,
        CoreCigarOp::D => AlignmentCigarOp::Deletion,
    };
    alignment_cigar_word(length, operation)
}

pub(crate) fn alignment_cigar_word(length: u64, operation: AlignmentCigarOp) -> Option<u32> {
    const MAX_BAM_CIGAR_LENGTH: u64 = (1_u64 << 28) - 1;
    if length == 0 || length > MAX_BAM_CIGAR_LENGTH {
        return None;
    }
    let operation = match operation {
        AlignmentCigarOp::Match => 0_u32,
        AlignmentCigarOp::Insertion => 1_u32,
        AlignmentCigarOp::Deletion => 2_u32,
        AlignmentCigarOp::SoftClip => 4_u32,
    };
    Some(u32::try_from(length).ok()? << 4 | operation)
}

/// A terminal-on-error owner of one exclusive native BAM staging file.
pub struct BamStagingWriter {
    path: PathBuf,
    staged: Option<StagedFile>,
    native: Option<NativeBamWriter>,
    direct_cigar: Vec<u32>,
    records_written: u64,
    terminal: bool,
}

impl BamStagingWriter {
    fn open_staged(
        staged: StagedFile,
        header_bytes: &[u8],
        compression_threads: u32,
        compression_level: Option<u8>,
    ) -> Result<Self, HtsError> {
        let path = staged.path().to_path_buf();
        let descriptor = staged.file().ok_or_else(|| {
            simple_error(
                HtsOperation::CreateStaging,
                &path,
                None,
                HtsErrorKind::Terminal,
            )
        })?;
        let descriptor_path = CString::new(format!("/proc/self/fd/{}", descriptor.as_raw_fd()))
            .map_err(|source| nul_error(&path, source))?;
        let native = match compression_level {
            Some(level) => NativeBamWriter::open_with_threads_and_compression_level(
                &descriptor_path,
                header_bytes,
                compression_threads,
                level,
            ),
            None => NativeBamWriter::open_with_threads(
                &descriptor_path,
                header_bytes,
                compression_threads,
            ),
        };
        match native {
            Ok(native) => Ok(Self {
                path,
                staged: Some(staged),
                native: Some(native),
                direct_cigar: Vec::new(),
                records_written: 0,
                terminal: false,
            }),
            Err(source) => {
                let primary = native_error(HtsOperation::OpenBam, &path, None, source);
                drop(staged);
                Err(primary)
            }
        }
    }

    /// Reserves an absent staging path and writes the canonical alignment header.
    ///
    /// Header encoding completes before path creation. The path is created with
    /// `create_new`; native open may only truncate the file this call just
    /// reserved. Failure removes only that adapter-owned path.
    ///
    /// # Errors
    ///
    /// Returns a path, header encoding, staging creation, or native header/open
    /// error. An existing path is never changed.
    pub fn create_new(
        path: impl AsRef<Path>,
        header: &SamHeader,
        limits: AlignmentRecordLimits,
    ) -> Result<Self, HtsError> {
        Self::create_new_with_threads(path, header, limits, 0)
    }

    /// Reserves a staging path and enables private `HTSlib` BGZF workers.
    ///
    /// `compression_threads == 0` preserves synchronous compression. The
    /// native shim rejects values above 64.
    ///
    /// # Errors
    ///
    /// Returns the same path, encoding, identity, and native errors as
    /// [`Self::create_new`], including failure to create compression workers.
    pub fn create_new_with_threads(
        path: impl AsRef<Path>,
        header: &SamHeader,
        limits: AlignmentRecordLimits,
        compression_threads: u32,
    ) -> Result<Self, HtsError> {
        let path = absolute_path(path.as_ref(), HtsOperation::ValidatePath)?;
        let header_bytes = sam_header_bytes(header, limits)
            .map_err(|source| encode_error(HtsOperation::EncodeHeader, &path, None, source))?;
        let staged = StagedFile::create_new(&path).map_err(map_bam_publication_error)?;
        Self::open_staged(staged, &header_bytes, compression_threads, None)
    }

    /// Reserves a staging path with private BGZF workers and an explicit
    /// compression level in `0..=9`.
    ///
    /// # Errors
    ///
    /// Returns the same path, encoding, identity, and native errors as
    /// [`Self::create_new_with_threads`].
    pub fn create_new_with_threads_and_compression_level(
        path: impl AsRef<Path>,
        header: &SamHeader,
        limits: AlignmentRecordLimits,
        compression_threads: u32,
        compression_level: u8,
    ) -> Result<Self, HtsError> {
        let path = absolute_path(path.as_ref(), HtsOperation::ValidatePath)?;
        let header_bytes = sam_header_bytes(header, limits)
            .map_err(|source| encode_error(HtsOperation::EncodeHeader, &path, None, source))?;
        let staged = StagedFile::create_new(&path).map_err(map_bam_publication_error)?;
        Self::open_staged(
            staged,
            &header_bytes,
            compression_threads,
            Some(compression_level),
        )
    }

    /// Creates a private sibling staging file beside an absent BAM target.
    ///
    /// The staging name is selected and reserved by the shared publication
    /// lifecycle. Callers provide only the final target and cannot collide on
    /// a fixed, user-visible staging name.
    ///
    /// # Errors
    ///
    /// Returns target validation, header encoding, staging reservation, or
    /// native BAM-open errors.
    pub fn create_sibling(
        target: impl AsRef<Path>,
        header: &SamHeader,
        limits: AlignmentRecordLimits,
    ) -> Result<Self, HtsError> {
        Self::create_sibling_with_threads(target, header, limits, 0)
    }

    /// Creates a private sibling staging file with BGZF worker threads.
    ///
    /// # Errors
    ///
    /// Returns the same errors as [`Self::create_sibling`], including native
    /// compression-worker setup failures.
    pub fn create_sibling_with_threads(
        target: impl AsRef<Path>,
        header: &SamHeader,
        limits: AlignmentRecordLimits,
        compression_threads: u32,
    ) -> Result<Self, HtsError> {
        let target = absolute_path(target.as_ref(), HtsOperation::ValidatePath)?;
        let header_bytes = sam_header_bytes(header, limits)
            .map_err(|source| encode_error(HtsOperation::EncodeHeader, &target, None, source))?;
        let staged =
            StagedFile::create_sibling(&target, "bam").map_err(map_bam_publication_error)?;
        Self::open_staged(staged, &header_bytes, compression_threads, None)
    }

    /// Creates a private sibling staging file with BGZF workers and an
    /// explicit compression level in `0..=9`.
    ///
    /// # Errors
    ///
    /// Returns the same errors as [`Self::create_sibling_with_threads`].
    pub fn create_sibling_with_threads_and_compression_level(
        target: impl AsRef<Path>,
        header: &SamHeader,
        limits: AlignmentRecordLimits,
        compression_threads: u32,
        compression_level: u8,
    ) -> Result<Self, HtsError> {
        let target = absolute_path(target.as_ref(), HtsOperation::ValidatePath)?;
        let header_bytes = sam_header_bytes(header, limits)
            .map_err(|source| encode_error(HtsOperation::EncodeHeader, &target, None, source))?;
        let staged =
            StagedFile::create_sibling(&target, "bam").map_err(map_bam_publication_error)?;
        Self::open_staged(
            staged,
            &header_bytes,
            compression_threads,
            Some(compression_level),
        )
    }

    /// Returns the number of complete native records written after the header.
    #[must_use]
    pub const fn records_written(&self) -> u64 {
        self.records_written
    }

    /// Returns the owned staging path.
    #[must_use]
    pub fn path(&self) -> &Path {
        &self.path
    }

    /// Encodes and writes one canonical alignment record.
    ///
    /// # Errors
    ///
    /// Returns an encoding, counter, terminal-state, or native write failure.
    /// Any error makes this writer terminal and it cannot yield a completed BAM.
    pub fn write_record(
        &mut self,
        record: &AlignmentRecord,
        limits: AlignmentRecordLimits,
    ) -> Result<(), HtsError> {
        let ordinal = self.records_written.checked_add(1).ok_or_else(|| {
            simple_error(
                HtsOperation::EncodeRecord,
                &self.path,
                None,
                HtsErrorKind::RecordCountOverflow,
            )
        })?;
        if self.terminal {
            return Err(simple_error(
                HtsOperation::WriteRecord,
                &self.path,
                Some(ordinal),
                HtsErrorKind::Terminal,
            ));
        }
        let record_bytes = match sam_record_bytes(record, limits) {
            Ok(bytes) => bytes,
            Err(source) => {
                self.terminal = true;
                return Err(encode_error(
                    HtsOperation::EncodeRecord,
                    &self.path,
                    Some(ordinal),
                    source,
                ));
            }
        };
        let Some(native) = self.native.as_mut() else {
            self.terminal = true;
            return Err(simple_error(
                HtsOperation::WriteRecord,
                &self.path,
                Some(ordinal),
                HtsErrorKind::Terminal,
            ));
        };
        if let Err(source) = native.write_record(&record_bytes) {
            self.terminal = true;
            return Err(native_error(
                HtsOperation::WriteRecord,
                &self.path,
                Some(ordinal),
                source,
            ));
        }
        self.records_written = ordinal;
        Ok(())
    }

    /// Writes one validated owned record as BAM fields.
    ///
    /// This audit boundary retains validated FLAG, coordinates, CIGAR,
    /// sequence, quality, NM, and XG semantics while avoiding canonical SAM
    /// rendering followed by native SAM parsing. SAM text-line limits do not
    /// apply because this path never constructs SAM text; native BAM limits are
    /// checked by the pinned shim.
    ///
    /// # Errors
    ///
    /// Returns a validation, representability, allocation, counter, terminal,
    /// or native write failure. Any error makes this writer terminal.
    #[doc(hidden)]
    #[allow(clippy::too_many_lines)]
    pub fn write_record_as_bam(&mut self, record: &AlignmentRecord) -> Result<(), HtsError> {
        let ordinal = self.records_written.checked_add(1).ok_or_else(|| {
            simple_error(
                HtsOperation::EncodeRecord,
                &self.path,
                None,
                HtsErrorKind::RecordCountOverflow,
            )
        })?;
        if self.terminal {
            return Err(simple_error(
                HtsOperation::WriteRecord,
                &self.path,
                Some(ordinal),
                HtsErrorKind::Terminal,
            ));
        }
        self.direct_cigar.clear();
        let auxiliary_mode = record
            .mapping()
            .map_or(AlignmentAuxiliaryMode::Minimal, |mapping| {
                mapping.auxiliary_mode()
            });
        let (
            reference_id,
            position,
            literal_nm_and_md,
            bisulfite_genome_conversion,
            bismark_auxiliary,
        ) = if let Some(mapping) = record.mapping() {
            let Some(reference_id) = reference_id_from_ordinal(mapping.reference().ordinal())
            else {
                return self.direct_encoding_failure(ordinal);
            };
            if self
                .direct_cigar
                .try_reserve(mapping.cigar().run_count())
                .is_err()
            {
                return self.direct_encoding_failure(ordinal);
            }
            for run in mapping.cigar().runs() {
                let Some(word) = core_cigar_word(run.length(), run.operation()) else {
                    return self.direct_encoding_failure(ordinal);
                };
                self.direct_cigar.push(word);
            }
            let (md, bismark_auxiliary) = match auxiliary_mode {
                AlignmentAuxiliaryMode::Minimal => (&[][..], None),
                AlignmentAuxiliaryMode::Bismark => {
                    let (Some(md), Some(xm)) = (mapping.md(), mapping.bismark_xm()) else {
                        return self.direct_encoding_failure(ordinal);
                    };
                    (md, Some((xm, mapping.bismark_xr())))
                }
            };
            (
                reference_id,
                i64::from(mapping.reference().position()) - 1,
                Some((mapping.literal_nm(), md)),
                Some(mapping.bismark_xg()),
                bismark_auxiliary,
            )
        } else {
            (-1, -1, None, None, None)
        };
        let (mate_reference_id, mate_position) = if let Some(mate) = record.mate() {
            let Some(reference_id) = reference_id_from_ordinal(mate.reference().ordinal()) else {
                return self.direct_encoding_failure(ordinal);
            };
            (reference_id, i64::from(mate.reference().position()) - 1)
        } else {
            (-1, -1)
        };
        let fields = NativeBamRecordFields {
            query_name: record.query_name(),
            flag: sam_flag(record),
            reference_id,
            position,
            mapping_quality: record.mapping_quality().sam_value(),
            cigar: &self.direct_cigar,
            mate_reference_id,
            mate_position,
            template_length: i64::from(record.template_length()),
            sequence: record.sequence(),
            quality: record.quality(),
            literal_nm_and_md,
            emit_md: matches!(auxiliary_mode, AlignmentAuxiliaryMode::Bismark),
            bisulfite_genome_conversion,
            bismark_auxiliary,
        };
        let Some(native) = self.native.as_mut() else {
            self.terminal = true;
            return Err(simple_error(
                HtsOperation::WriteRecord,
                &self.path,
                Some(ordinal),
                HtsErrorKind::Terminal,
            ));
        };
        if let Err(source) = native.write_bam_fields(&fields) {
            self.terminal = true;
            return Err(native_error(
                HtsOperation::WriteRecord,
                &self.path,
                Some(ordinal),
                source,
            ));
        }
        self.records_written = ordinal;
        Ok(())
    }

    /// Writes one compact batch-backed primary record without intermediate ownership.
    ///
    /// # Errors
    ///
    /// Returns validation, representability, terminal-state, or native write
    /// failures and makes this writer terminal on error.
    #[doc(hidden)]
    #[allow(clippy::too_many_lines)]
    pub fn write_borrowed_alignment_record(
        &mut self,
        record: &BorrowedAlignmentRecord<'_>,
    ) -> Result<(), HtsError> {
        let ordinal = self.records_written.checked_add(1).ok_or_else(|| {
            simple_error(
                HtsOperation::EncodeRecord,
                &self.path,
                None,
                HtsErrorKind::RecordCountOverflow,
            )
        })?;
        if self.terminal {
            return Err(simple_error(
                HtsOperation::WriteRecord,
                &self.path,
                Some(ordinal),
                HtsErrorKind::Terminal,
            ));
        }
        self.direct_cigar.clear();
        if self.direct_cigar.try_reserve(record.cigar().len()).is_err() {
            return self.direct_encoding_failure(ordinal);
        }
        for run in record.cigar() {
            let Some(word) = alignment_cigar_word(run.length(), run.operation()) else {
                return self.direct_encoding_failure(ordinal);
            };
            self.direct_cigar.push(word);
        }
        let reference_id = match record.reference_ordinal() {
            Some(reference_ordinal) => {
                let Some(reference_id) = reference_id_from_ordinal(reference_ordinal) else {
                    return self.direct_encoding_failure(ordinal);
                };
                reference_id
            }
            None => -1,
        };
        let mate_reference_id = match record.mate_reference_ordinal() {
            Some(reference_ordinal) => {
                let Some(reference_id) = reference_id_from_ordinal(reference_ordinal) else {
                    return self.direct_encoding_failure(ordinal);
                };
                reference_id
            }
            None => -1,
        };
        let (md, bismark_auxiliary) = match record.auxiliary_mode() {
            AlignmentAuxiliaryMode::Minimal => (None, None),
            AlignmentAuxiliaryMode::Bismark => match (record.md(), record.bismark_xm()) {
                (Some(md), Some(xm)) => (Some(md), Some((xm, record.bismark_xr()))),
                _ => return self.direct_encoding_failure(ordinal),
            },
        };
        let fields = NativeBamRecordFields {
            query_name: record.query_name(),
            flag: record.flag(),
            reference_id,
            position: if reference_id < 0 {
                -1
            } else {
                i64::from(record.position()) - 1
            },
            mapping_quality: record.mapping_quality(),
            cigar: &self.direct_cigar,
            mate_reference_id,
            mate_position: if mate_reference_id < 0 {
                -1
            } else {
                i64::from(record.mate_position()) - 1
            },
            template_length: i64::from(record.template_length()),
            sequence: record.sequence(),
            quality: record.quality(),
            literal_nm_and_md: (reference_id >= 0)
                .then(|| (record.literal_nm(), md.unwrap_or(&[]))),
            emit_md: reference_id >= 0
                && matches!(record.auxiliary_mode(), AlignmentAuxiliaryMode::Bismark),
            bisulfite_genome_conversion: (reference_id >= 0).then(|| record.bismark_xg()),
            bismark_auxiliary,
        };
        let Some(native) = self.native.as_mut() else {
            self.terminal = true;
            return Err(simple_error(
                HtsOperation::WriteRecord,
                &self.path,
                Some(ordinal),
                HtsErrorKind::Terminal,
            ));
        };
        if let Err(source) = native.write_bam_fields(&fields) {
            self.terminal = true;
            return Err(native_error(
                HtsOperation::WriteRecord,
                &self.path,
                Some(ordinal),
                source,
            ));
        }
        self.records_written = ordinal;
        Ok(())
    }

    fn direct_encoding_failure(&mut self, ordinal: u64) -> Result<(), HtsError> {
        self.terminal = true;
        Err(simple_error(
            HtsOperation::EncodeRecord,
            &self.path,
            Some(ordinal),
            HtsErrorKind::Encode,
        ))
    }

    /// Finalizes the BAM and transfers ownership to a completed staging value.
    ///
    /// # Errors
    ///
    /// Returns a terminal or native finalize failure and removes only the
    /// adapter-owned staging path.
    pub fn finish(mut self) -> Result<CompletedBam, HtsError> {
        if self.terminal {
            self.native.take();
            self.cleanup_owned();
            return Err(simple_error(
                HtsOperation::FinishBam,
                &self.path,
                None,
                HtsErrorKind::Terminal,
            ));
        }
        let Some(mut native) = self.native.take() else {
            self.cleanup_owned();
            return Err(simple_error(
                HtsOperation::FinishBam,
                &self.path,
                None,
                HtsErrorKind::Terminal,
            ));
        };
        if let Err(source) = native.finish() {
            drop(native);
            let primary = native_error(HtsOperation::FinishBam, &self.path, None, source);
            self.cleanup_owned();
            return Err(primary);
        }
        drop(native);
        let Some(mut staged) = self.staged.take() else {
            self.cleanup_owned();
            return Err(simple_error(
                HtsOperation::FinishBam,
                &self.path,
                None,
                HtsErrorKind::Terminal,
            ));
        };
        let anchor = staged.take_file().map_err(map_bam_publication_error)?;
        let completed = staged.complete(anchor).map_err(map_bam_publication_error)?;
        Ok(CompletedBam {
            completed,
            records_written: self.records_written,
        })
    }

    /// Aborts the writer and removes only its private staging path.
    ///
    /// # Errors
    ///
    /// Returns a cleanup error after the native handle has been destroyed.
    pub fn abort(mut self) -> Result<(), HtsError> {
        self.native.take();
        self.remove_owned()
    }

    fn cleanup_owned(&mut self) {
        self.staged.take();
    }

    fn remove_owned(&mut self) -> Result<(), HtsError> {
        match self.staged.take() {
            Some(staged) => staged.abort().map_err(map_bam_publication_error),
            None => Ok(()),
        }
    }
}

impl Drop for BamStagingWriter {
    fn drop(&mut self) {
        self.native.take();
        self.cleanup_owned();
    }
}

/// One completely finalized BAM staging file still owned by the adapter.
pub struct CompletedBam {
    completed: CompletedFile,
    records_written: u64,
}

impl CompletedBam {
    /// Returns the completed staging path.
    #[must_use]
    pub fn path(&self) -> &Path {
        self.completed.staging_path()
    }

    /// Returns the number of records finalized after the header.
    #[must_use]
    pub const fn records_written(&self) -> u64 {
        self.records_written
    }

    /// Transfers a still-identity-verified path to a later publication owner.
    ///
    /// # Errors
    ///
    /// Returns an identity or metadata failure if the path was removed or
    /// replaced after finalization.
    pub fn into_path(self) -> Result<PathBuf, HtsError> {
        self.completed
            .into_path()
            .map_err(map_bam_publication_error)
    }

    /// Synchronizes and atomically publishes this BAM at an absent sibling path.
    ///
    /// Publication links the retained staging descriptor rather than reopening
    /// its pathname. A concurrent target creator wins without being modified.
    /// Once the target link succeeds, staging cleanup failure is returned as a
    /// warning in [`BamPublication`] rather than converting success into an
    /// error with a visible target.
    ///
    /// # Errors
    ///
    /// Returns a path, target-existence, staging-identity, synchronization, or
    /// descriptor-link failure. An error return never creates the target.
    pub fn publish_create_new(self, target: impl AsRef<Path>) -> Result<BamPublication, HtsError> {
        let published = self
            .completed
            .publish_create_new_at(target)
            .map_err(map_bam_publication_error)?;
        Ok(BamPublication {
            published,
            records_written: self.records_written,
        })
    }

    /// Explicitly removes this completed staging file.
    ///
    /// # Errors
    ///
    /// Returns a direct filesystem cleanup failure.
    pub fn remove(self) -> Result<(), HtsError> {
        self.completed.remove().map_err(map_bam_publication_error)
    }
}

/// Successful create-only BAM publication details.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct BamPublication {
    published: PublishedFile,
    records_written: u64,
}

impl BamPublication {
    /// Returns the absolute final target path.
    #[must_use]
    pub fn target_path(&self) -> &Path {
        self.published.target_path()
    }

    /// Returns the absolute staging path used before publication.
    #[must_use]
    pub fn staging_path(&self) -> &Path {
        self.published.staging_path()
    }

    /// Returns the number of complete records in the published BAM.
    #[must_use]
    pub const fn records_written(&self) -> u64 {
        self.records_written
    }

    /// Returns a post-publication staging cleanup warning.
    #[must_use]
    pub fn cleanup_warning(&self) -> Option<HtsErrorKind> {
        self.published.cleanup_warning().map(|kind| {
            if kind == io::ErrorKind::Other {
                HtsErrorKind::StagingIdentityChanged
            } else {
                HtsErrorKind::Io(kind)
            }
        })
    }
}

fn map_bam_publication_error(error: PublicationError) -> HtsError {
    let phase = error.phase();
    let operation = match phase {
        bsbit_io::PublicationPhase::ValidatePaths => HtsOperation::ValidatePublicationPaths,
        bsbit_io::PublicationPhase::CreateStaging => HtsOperation::CreateStaging,
        bsbit_io::PublicationPhase::ValidateStaging => HtsOperation::ValidateStaging,
        bsbit_io::PublicationPhase::Sync => HtsOperation::SyncBam,
        bsbit_io::PublicationPhase::Publish => HtsOperation::PublishBam,
        bsbit_io::PublicationPhase::Cleanup => HtsOperation::Cleanup,
        bsbit_io::PublicationPhase::Rollback => HtsOperation::RollbackOutput,
    };
    let path = error.path().to_path_buf();
    let kind = error.kind();
    if matches!(
        (phase, kind),
        (
            bsbit_io::PublicationPhase::ValidateStaging | bsbit_io::PublicationPhase::Cleanup,
            io::ErrorKind::Other
        )
    ) {
        return simple_error(operation, &path, None, HtsErrorKind::StagingIdentityChanged);
    }
    if phase == bsbit_io::PublicationPhase::ValidatePaths && kind == io::ErrorKind::InvalidInput {
        return simple_error(
            operation,
            &path,
            None,
            HtsErrorKind::PublicationPathMismatch,
        );
    }
    io_error(operation, &path, None, error.into_io_error())
}
