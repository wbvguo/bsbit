//! Shared FASTQ resource contract for both alignment layouts.

use bsbit_align::MAX_READ_BASES;
use bsbit_hts::{SAM_MAX_QUERY_NAME_BYTES, TextRecordLimits};

const MAX_FASTQ_DESCRIPTION_BYTES: u64 = 1_000_000;
const MAX_FASTQ_LINE_BYTES: u64 = MAX_FASTQ_DESCRIPTION_BYTES + SAM_MAX_QUERY_NAME_BYTES + 1;

pub(crate) const fn alignment_fastq_limits() -> TextRecordLimits {
    TextRecordLimits::MAX
        .with_max_line_bytes(MAX_FASTQ_LINE_BYTES)
        .with_max_name_bytes(SAM_MAX_QUERY_NAME_BYTES)
        .with_max_description_bytes(MAX_FASTQ_DESCRIPTION_BYTES)
        .with_max_bases_per_record(MAX_READ_BASES as u64)
        .with_max_quality_bytes(MAX_READ_BASES as u64)
}
