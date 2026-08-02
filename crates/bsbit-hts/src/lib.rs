//! Sequencing-format models, codecs, and the isolated `HTSlib` adapter.
//!
//! File creation and publication are delegated to `bsbit-io`; alignment and
//! index algorithms are deliberately outside this crate.

#![deny(unsafe_code)]

mod alignment_record;
mod bam;
mod bed;
mod bed_methyl;
mod fasta;
mod fastq;
mod htslib;
mod sam;
#[allow(unsafe_code)]
mod sys;
mod text;
mod text_output;

pub use alignment_record::{
    AlignmentAuxiliaryMode, AlignmentCigarOp, AlignmentCigarRun, AlignmentPlacement, AlignmentRead,
    AlignmentRecord, AlignmentRecordAllocation, AlignmentRecordBatch, AlignmentRecordError,
    AlignmentRecordField, AlignmentRecordLimits, AlignmentRecordResource, BorrowedAlignmentRead,
    BorrowedAlignmentRecord, MappedAlignmentRecord, RecordMappingQuality, RecordMateLocation,
    RecordReference, RecordSegment, SAM_MAX_QUERY_NAME_BYTES, SAM_MAX_REFERENCE_LENGTH,
};
pub use bam::{
    BamAlignmentColumn, BamCigarOperation, BamCigarRun, BamPublication, BamRecordDecodeWorkspace,
    BamRecordFieldError, BamStagingWriter, CompletedBam, IndexedBamHeader, IndexedBamReader,
    IndexedBamRecord, IndexedBamReference, build_bam_index_create_new,
};
pub use bed::{BedError, BedInterval};
pub use bed_methyl::{BedMethylContext, BedMethylError, BedMethylRecord, BedMethylStrand};
pub use fasta::{
    DecodedFastaReader, FastaReader, FastaRecord, IndexedFastaReader, IndexedFastaReference,
};
pub use fastq::{
    BorrowedFastqRecord, DecodedFastqReader, DecodedPairedFastqReader, FastqReader, FastqRecord,
    FastqRecordBatch, PairSourceSide, PairedFastqReader, PairedFastqRecord,
};
pub use htslib::{
    BgzfWriter, Compression, DecodedReader, HtsError, HtsErrorKind, HtsOperation, NativeError,
    NativeStatus,
};
pub use sam::{
    BsbitAlignmentMode, BsbitProgramProvenance, BsbitProgramProvenanceError, SamFileError,
    SamFilePhase, SamFilePublication, SamFileWriter, SamHeader, SamHeaderReference, SamSortOrder,
    SamWriteError, SamWritePhase, sam_borrowed_record_bytes, sam_flag, sam_header_bytes,
    sam_record_bytes, write_sam_header, write_sam_record,
};
pub use text::{
    RecordField, RecordName, RecordOrdinal, TextRecordAllocation, TextRecordError,
    TextRecordErrorKind, TextRecordFormat, TextRecordLimits, TextRecordResource,
};
pub use text_output::{
    CompletedTextOutput, TextOutputCompression, TextPublication, TextStagingWriter,
};

pub(crate) use htslib::{io_error, simple_error};
