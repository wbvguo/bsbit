//! BAM reader and writer adapters with separate native lifecycles.

mod reader;
mod writer;

pub use reader::{
    BamAlignmentColumn, BamCigarOperation, BamCigarRun, BamRecordDecodeWorkspace,
    BamRecordFieldError, IndexedBamHeader, IndexedBamReader, IndexedBamRecord, IndexedBamReference,
    build_bam_index_create_new,
};
pub use writer::{BamPublication, BamStagingWriter, CompletedBam};
