//! Generic identity-safe local-file creation and publication.
//!
//! This crate deliberately has no knowledge of biological formats, indexes,
//! alignments, native libraries, or command orchestration.

#![deny(unsafe_code)]
#![deny(unsafe_op_in_unsafe_fn)]

#[allow(unsafe_code)]
mod file;
mod publication;

pub use file::{
    FileIdentity, absolute_path, create_new, hard_link_descriptor_create_new,
    remove_if_identity_matches, reopen_read_write, validate_absent, validate_create_target,
    validate_distinct_paths, validate_regular_file_or_absent, validate_replace_target,
};
pub use publication::{
    CompletedFile, PublicationError, PublicationPhase, PublishedFile, StagedFile,
    select_sibling_staging_path, select_sibling_staging_path_replace,
    validate_sibling_publication_paths,
};
