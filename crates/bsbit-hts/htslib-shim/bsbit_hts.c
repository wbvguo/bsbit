#define _GNU_SOURCE
#include "bsbit_hts.h"

#include <errno.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>

#include <htslib/bgzf.h>
#include <htslib/faidx.h>
#include <htslib/hts.h>
#include <htslib/kstring.h>
#include <htslib/sam.h>
#include <htslib/tbx.h>

#if HTS_VERSION != 102400
#error "bsbit HTS shim requires HTSlib 1.24 exactly"
#endif

struct bsbit_hts_reader {
    BGZF *file;
    int compression;
    int failed;
    int closed;
};

struct bsbit_hts_indexed_reader {
    samFile *file;
    sam_hdr_t *header;
    hts_idx_t *index;
    hts_itr_t *iterator;
    bam1_t *record;
    int failed;
    int closed;
};

struct bsbit_hts_indexed_fasta_reader {
    faidx_t *index;
    char *sequence;
    int failed;
    int closed;
};

struct bsbit_hts_bgzf_writer {
    BGZF *file;
    int failed;
    int finished;
};

struct bsbit_hts_writer {
    samFile *file;
    sam_hdr_t *header;
    bam1_t *record;
    int failed;
    int finished;
};

static int set_result(int status,
                      int system_errno,
                      const char *message,
                      int *out_system_errno,
                      char *error,
                      size_t error_capacity) {
    size_t copy_length = 0;

    if (out_system_errno != NULL) {
        *out_system_errno = system_errno;
    }
    if (error != NULL && error_capacity > 0) {
        error[0] = '\0';
        if (message != NULL) {
            copy_length = strlen(message);
            if (copy_length >= error_capacity) {
                copy_length = error_capacity - 1;
            }
            if (copy_length > 0) {
                memcpy(error, message, copy_length);
            }
            error[copy_length] = '\0';
        }
    }
    return status;
}

static void cleanup_writer(bsbit_hts_writer *writer) {
    if (writer == NULL) {
        return;
    }
    if (writer->file != NULL) {
        (void)hts_close(writer->file);
        writer->file = NULL;
    }
    bam_destroy1(writer->record);
    writer->record = NULL;
    sam_hdr_destroy(writer->header);
    writer->header = NULL;
}

static int cleanup_indexed_reader(bsbit_hts_indexed_reader *reader) {
    int close_result = 0;

    if (reader == NULL) {
        return 0;
    }
    hts_itr_destroy(reader->iterator);
    reader->iterator = NULL;
    hts_idx_destroy(reader->index);
    reader->index = NULL;
    bam_destroy1(reader->record);
    reader->record = NULL;
    sam_hdr_destroy(reader->header);
    reader->header = NULL;
    if (reader->file != NULL) {
        close_result = hts_close(reader->file);
        reader->file = NULL;
    }
    return close_result;
}

static void cleanup_indexed_fasta_reader(
    bsbit_hts_indexed_fasta_reader *reader) {
    if (reader == NULL) {
        return;
    }
    free(reader->sequence);
    reader->sequence = NULL;
    fai_destroy(reader->index);
    reader->index = NULL;
}

uint32_t bsbit_hts_shim_abi_version(void) {
    return UINT32_C(3);
}

const char *bsbit_hts_runtime_version(void) {
    return hts_version();
}

int bsbit_hts_health_check(void) {
    const char *version = hts_version();
    return version != NULL && strcmp(version, "1.24") == 0 ? BSBIT_HTS_OK
                                                            : BSBIT_HTS_INVALID_STATE;
}

int bsbit_hts_tabix_index_build(const char *path,
                                const char *index_path,
                                int preset,
                                uint32_t threads,
                                int *out_system_errno,
                                char *error,
                                size_t error_capacity) {
    const tbx_conf_t *configuration = NULL;
    int result = 0;
    int saved_errno = 0;

    if (path == NULL || path[0] == '\0' || index_path == NULL ||
        index_path[0] == '\0' || threads > UINT32_C(64)) {
        return set_result(BSBIT_HTS_INVALID_ARGUMENT, 0,
                          "invalid tabix index arguments", out_system_errno,
                          error, error_capacity);
    }
    switch (preset) {
        case BSBIT_HTS_TABIX_VCF:
            configuration = &tbx_conf_vcf;
            break;
        case BSBIT_HTS_TABIX_BED:
            configuration = &tbx_conf_bed;
            break;
        default:
            return set_result(BSBIT_HTS_INVALID_ARGUMENT, 0,
                              "unknown tabix index preset", out_system_errno,
                              error, error_capacity);
    }
    errno = 0;
    result = tbx_index_build3(path, index_path, 0, (int)threads,
                              configuration);
    saved_errno = errno;
    if (result != 0) {
        return set_result(BSBIT_HTS_WRITE_FAILED, saved_errno,
                          "HTSlib could not build the tabix index",
                          out_system_errno, error, error_capacity);
    }
    return set_result(BSBIT_HTS_OK, 0, NULL, out_system_errno, error,
                      error_capacity);
}

int bsbit_hts_bam_index_build(const char *path,
                              const char *index_path,
                              uint32_t threads,
                              int *out_system_errno,
                              char *error,
                              size_t error_capacity) {
    int result = 0;
    int saved_errno = 0;

    if (path == NULL || path[0] == '\0' || index_path == NULL ||
        index_path[0] == '\0' || threads > UINT32_C(64)) {
        return set_result(BSBIT_HTS_INVALID_ARGUMENT, 0,
                          "invalid BAM index arguments", out_system_errno,
                          error, error_capacity);
    }
    errno = 0;
    result = sam_index_build3(path, index_path, 0, (int)threads);
    saved_errno = errno;
    if (result != 0) {
        return set_result(BSBIT_HTS_WRITE_FAILED, saved_errno,
                          "HTSlib could not build the BAM index",
                          out_system_errno, error, error_capacity);
    }
    return set_result(BSBIT_HTS_OK, 0, NULL, out_system_errno, error,
                      error_capacity);
}

int bsbit_hts_reader_open(const char *path,
                          bsbit_hts_reader **out_reader,
                          int *out_system_errno,
                          char *error,
                          size_t error_capacity) {
    bsbit_hts_reader *reader = NULL;
    int saved_errno = 0;

    if (out_reader == NULL) {
        return set_result(BSBIT_HTS_INVALID_ARGUMENT, 0, "out_reader is null",
                          out_system_errno, error, error_capacity);
    }
    *out_reader = NULL;
    if (path == NULL || path[0] == '\0') {
        return set_result(BSBIT_HTS_INVALID_ARGUMENT, 0,
                          "a nonempty path is required", out_system_errno, error,
                          error_capacity);
    }
    reader = (bsbit_hts_reader *)calloc(1, sizeof(*reader));
    if (reader == NULL) {
        return set_result(BSBIT_HTS_ALLOCATION_FAILED, errno,
                          "failed to allocate reader", out_system_errno, error,
                          error_capacity);
    }
    errno = 0;
    reader->file = bgzf_open(path, "r");
    saved_errno = errno;
    if (reader->file == NULL) {
        free(reader);
        return set_result(BSBIT_HTS_OPEN_FAILED, saved_errno,
                          "HTSlib could not open input path", out_system_errno,
                          error, error_capacity);
    }
    reader->compression = bgzf_compression(reader->file);
    if (reader->compression < BSBIT_HTS_PLAIN ||
        reader->compression > BSBIT_HTS_BGZF) {
        (void)bgzf_close(reader->file);
        free(reader);
        return set_result(BSBIT_HTS_OPEN_FAILED, 0,
                          "HTSlib returned an unknown compression class",
                          out_system_errno, error, error_capacity);
    }
    *out_reader = reader;
    return set_result(BSBIT_HTS_OK, 0, NULL, out_system_errno, error,
                      error_capacity);
}

int bsbit_hts_reader_compression(const bsbit_hts_reader *reader,
                                 int *out_compression,
                                 int *out_system_errno,
                                 char *error,
                                 size_t error_capacity) {
    if (reader == NULL || out_compression == NULL) {
        return set_result(BSBIT_HTS_INVALID_ARGUMENT, 0,
                          "reader and out_compression are required",
                          out_system_errno, error, error_capacity);
    }
    if (reader->closed) {
        return set_result(BSBIT_HTS_INVALID_STATE, 0, "reader is closed",
                          out_system_errno, error, error_capacity);
    }
    *out_compression = reader->compression;
    return set_result(BSBIT_HTS_OK, 0, NULL, out_system_errno, error,
                      error_capacity);
}

int bsbit_hts_reader_read(bsbit_hts_reader *reader,
                          uint8_t *buffer,
                          size_t capacity,
                          size_t *out_count,
                          int *out_system_errno,
                          char *error,
                          size_t error_capacity) {
    ssize_t count = 0;
    int saved_errno = 0;

    if (out_count != NULL) {
        *out_count = 0;
    }
    if (reader == NULL || out_count == NULL ||
        (capacity > 0 && buffer == NULL)) {
        return set_result(BSBIT_HTS_INVALID_ARGUMENT, 0,
                          "reader, out_count, and nonempty buffer are required",
                          out_system_errno, error, error_capacity);
    }
    if (reader->closed) {
        return set_result(BSBIT_HTS_INVALID_STATE, 0, "reader is closed",
                          out_system_errno, error, error_capacity);
    }
    if (reader->failed) {
        return set_result(BSBIT_HTS_READ_FAILED, 0,
                          "reader is terminal after an earlier read failure",
                          out_system_errno, error, error_capacity);
    }
    if (capacity == 0) {
        return set_result(BSBIT_HTS_OK, 0, NULL, out_system_errno, error,
                          error_capacity);
    }
    if (capacity > (size_t)SSIZE_MAX) {
        reader->failed = 1;
        return set_result(BSBIT_HTS_INVALID_ARGUMENT, 0,
                          "read capacity exceeds ssize_t", out_system_errno,
                          error, error_capacity);
    }
    errno = 0;
    count = bgzf_read(reader->file, buffer, capacity);
    saved_errno = errno;
    if (count < 0) {
        reader->failed = 1;
        return set_result(BSBIT_HTS_READ_FAILED, saved_errno,
                          "HTSlib failed while decoding input", out_system_errno,
                          error, error_capacity);
    }
    *out_count = (size_t)count;
    return set_result(BSBIT_HTS_OK, 0, NULL, out_system_errno, error,
                      error_capacity);
}

int bsbit_hts_reader_close(bsbit_hts_reader *reader,
                           int *out_system_errno,
                           char *error,
                           size_t error_capacity) {
    int result = 0;
    int saved_errno = 0;

    if (reader == NULL) {
        return set_result(BSBIT_HTS_INVALID_ARGUMENT, 0, "reader is null",
                          out_system_errno, error, error_capacity);
    }
    if (reader->closed) {
        return set_result(BSBIT_HTS_INVALID_STATE, 0,
                          "reader is already closed", out_system_errno, error,
                          error_capacity);
    }
    reader->closed = 1;
    errno = 0;
    result = bgzf_close(reader->file);
    saved_errno = errno;
    reader->file = NULL;
    if (result < 0) {
        return set_result(BSBIT_HTS_CLOSE_FAILED, saved_errno,
                          "HTSlib failed while closing input", out_system_errno,
                          error, error_capacity);
    }
    return set_result(BSBIT_HTS_OK, 0, NULL, out_system_errno, error,
                      error_capacity);
}

void bsbit_hts_reader_destroy(bsbit_hts_reader *reader) {
    if (reader == NULL) {
        return;
    }
    if (reader->file != NULL) {
        (void)bgzf_close(reader->file);
        reader->file = NULL;
    }
    free(reader);
}

int bsbit_hts_indexed_reader_open(const char *path,
                                  bsbit_hts_indexed_reader **out_reader,
                                  int *out_system_errno,
                                  char *error,
                                  size_t error_capacity) {
    bsbit_hts_indexed_reader *reader = NULL;
    const htsFormat *format = NULL;
    int saved_errno = 0;

    if (out_reader == NULL) {
        return set_result(BSBIT_HTS_INVALID_ARGUMENT, 0, "out_reader is null",
                          out_system_errno, error, error_capacity);
    }
    *out_reader = NULL;
    if (path == NULL || path[0] == '\0') {
        return set_result(BSBIT_HTS_INVALID_ARGUMENT, 0,
                          "a nonempty path is required", out_system_errno,
                          error, error_capacity);
    }
    reader = (bsbit_hts_indexed_reader *)calloc(1, sizeof(*reader));
    if (reader == NULL) {
        return set_result(BSBIT_HTS_ALLOCATION_FAILED, errno,
                          "failed to allocate indexed BAM reader",
                          out_system_errno, error, error_capacity);
    }
    errno = 0;
    reader->file = sam_open(path, "rb");
    saved_errno = errno;
    if (reader->file == NULL) {
        free(reader);
        return set_result(BSBIT_HTS_OPEN_FAILED, saved_errno,
                          "HTSlib could not open indexed BAM input",
                          out_system_errno, error, error_capacity);
    }
    format = hts_get_format(reader->file);
    if (format == NULL || format->format != bam) {
        (void)cleanup_indexed_reader(reader);
        free(reader);
        return set_result(BSBIT_HTS_OPEN_FAILED, 0,
                          "indexed input is not BAM", out_system_errno, error,
                          error_capacity);
    }
    reader->header = sam_hdr_read(reader->file);
    if (reader->header == NULL) {
        (void)cleanup_indexed_reader(reader);
        free(reader);
        return set_result(BSBIT_HTS_HEADER_FAILED, 0,
                          "HTSlib could not read BAM header", out_system_errno,
                          error, error_capacity);
    }
    errno = 0;
    reader->index = sam_index_load(reader->file, path);
    saved_errno = errno;
    if (reader->index == NULL) {
        (void)cleanup_indexed_reader(reader);
        free(reader);
        return set_result(BSBIT_HTS_OPEN_FAILED, saved_errno,
                          "HTSlib could not load a BAM index",
                          out_system_errno, error, error_capacity);
    }
    reader->record = bam_init1();
    if (reader->record == NULL) {
        (void)cleanup_indexed_reader(reader);
        free(reader);
        return set_result(BSBIT_HTS_ALLOCATION_FAILED, errno,
                          "failed to allocate indexed BAM record",
                          out_system_errno, error, error_capacity);
    }
    *out_reader = reader;
    return set_result(BSBIT_HTS_OK, 0, NULL, out_system_errno, error,
                      error_capacity);
}

static int indexed_reader_is_live(const bsbit_hts_indexed_reader *reader,
                                  int *out_system_errno,
                                  char *error,
                                  size_t error_capacity) {
    if (reader == NULL) {
        return set_result(BSBIT_HTS_INVALID_ARGUMENT, 0, "reader is null",
                          out_system_errno, error, error_capacity);
    }
    if (reader->closed) {
        return set_result(BSBIT_HTS_INVALID_STATE, 0, "reader is closed",
                          out_system_errno, error, error_capacity);
    }
    if (reader->failed) {
        return set_result(BSBIT_HTS_READ_FAILED, 0,
                          "reader is terminal after an earlier indexed read failure",
                          out_system_errno, error, error_capacity);
    }
    return BSBIT_HTS_OK;
}

int bsbit_hts_indexed_reader_header_text(
    const bsbit_hts_indexed_reader *reader,
    const char **out_text,
    size_t *out_length,
    int *out_system_errno,
    char *error,
    size_t error_capacity) {
    int status = indexed_reader_is_live(reader, out_system_errno, error,
                                        error_capacity);
    if (out_text != NULL) {
        *out_text = NULL;
    }
    if (out_length != NULL) {
        *out_length = 0;
    }
    if (status != BSBIT_HTS_OK) {
        return status;
    }
    if (out_text == NULL || out_length == NULL) {
        return set_result(BSBIT_HTS_INVALID_ARGUMENT, 0,
                          "out_text and out_length are required",
                          out_system_errno, error, error_capacity);
    }
    *out_text = sam_hdr_str(reader->header);
    *out_length = sam_hdr_length(reader->header);
    if (*out_text == NULL && *out_length != 0) {
        return set_result(BSBIT_HTS_INVALID_STATE, 0,
                          "HTSlib returned a null nonempty header",
                          out_system_errno, error, error_capacity);
    }
    return set_result(BSBIT_HTS_OK, 0, NULL, out_system_errno, error,
                      error_capacity);
}

int bsbit_hts_indexed_reader_reference_count(
    const bsbit_hts_indexed_reader *reader,
    int32_t *out_count,
    int *out_system_errno,
    char *error,
    size_t error_capacity) {
    int count = 0;
    int status = indexed_reader_is_live(reader, out_system_errno, error,
                                        error_capacity);
    if (out_count != NULL) {
        *out_count = 0;
    }
    if (status != BSBIT_HTS_OK) {
        return status;
    }
    if (out_count == NULL) {
        return set_result(BSBIT_HTS_INVALID_ARGUMENT, 0,
                          "out_count is required", out_system_errno, error,
                          error_capacity);
    }
    count = sam_hdr_nref(reader->header);
    if (count < 0) {
        return set_result(BSBIT_HTS_HEADER_FAILED, 0,
                          "HTSlib returned a negative reference count",
                          out_system_errno, error, error_capacity);
    }
    *out_count = (int32_t)count;
    return set_result(BSBIT_HTS_OK, 0, NULL, out_system_errno, error,
                      error_capacity);
}

int bsbit_hts_indexed_reader_reference(
    const bsbit_hts_indexed_reader *reader,
    int32_t reference_id,
    const char **out_name,
    size_t *out_name_length,
    int64_t *out_length,
    int *out_system_errno,
    char *error,
    size_t error_capacity) {
    const char *name = NULL;
    hts_pos_t length = 0;
    int count = 0;
    int status = indexed_reader_is_live(reader, out_system_errno, error,
                                        error_capacity);
    if (out_name != NULL) {
        *out_name = NULL;
    }
    if (out_name_length != NULL) {
        *out_name_length = 0;
    }
    if (out_length != NULL) {
        *out_length = 0;
    }
    if (status != BSBIT_HTS_OK) {
        return status;
    }
    if (out_name == NULL || out_name_length == NULL || out_length == NULL) {
        return set_result(BSBIT_HTS_INVALID_ARGUMENT, 0,
                          "reference outputs are required", out_system_errno,
                          error, error_capacity);
    }
    count = sam_hdr_nref(reader->header);
    if (reference_id < 0 || reference_id >= count) {
        return set_result(BSBIT_HTS_INVALID_ARGUMENT, 0,
                          "reference id is outside the BAM dictionary",
                          out_system_errno, error, error_capacity);
    }
    name = sam_hdr_tid2name(reader->header, reference_id);
    length = sam_hdr_tid2len(reader->header, reference_id);
    if (name == NULL || length < 0) {
        return set_result(BSBIT_HTS_HEADER_FAILED, 0,
                          "HTSlib could not resolve a BAM reference",
                          out_system_errno, error, error_capacity);
    }
    *out_name = name;
    *out_name_length = strlen(name);
    *out_length = (int64_t)length;
    return set_result(BSBIT_HTS_OK, 0, NULL, out_system_errno, error,
                      error_capacity);
}

int bsbit_hts_indexed_reader_query(bsbit_hts_indexed_reader *reader,
                                   int32_t reference_id,
                                   int64_t start,
                                   int64_t end,
                                   int *out_system_errno,
                                   char *error,
                                   size_t error_capacity) {
    hts_pos_t reference_length = 0;
    int count = 0;
    int status = indexed_reader_is_live(reader, out_system_errno, error,
                                        error_capacity);
    if (status != BSBIT_HTS_OK) {
        return status;
    }
    count = sam_hdr_nref(reader->header);
    if (reference_id < 0 || reference_id >= count || start < 0 || end <= start) {
        return set_result(BSBIT_HTS_INVALID_ARGUMENT, 0,
                          "query requires an in-range reference and nonempty interval",
                          out_system_errno, error, error_capacity);
    }
    reference_length = sam_hdr_tid2len(reader->header, reference_id);
    if (reference_length < 0 || end > (int64_t)reference_length) {
        return set_result(BSBIT_HTS_INVALID_ARGUMENT, 0,
                          "query interval is outside the BAM reference",
                          out_system_errno, error, error_capacity);
    }
    hts_itr_destroy(reader->iterator);
    reader->iterator = sam_itr_queryi(reader->index, reference_id,
                                      (hts_pos_t)start, (hts_pos_t)end);
    if (reader->iterator == NULL) {
        return set_result(BSBIT_HTS_INVALID_STATE, 0,
                          "HTSlib could not create the BAM region iterator",
                          out_system_errno, error, error_capacity);
    }
    return set_result(BSBIT_HTS_OK, 0, NULL, out_system_errno, error,
                      error_capacity);
}

int bsbit_hts_indexed_reader_next(bsbit_hts_indexed_reader *reader,
                                  bsbit_hts_bam_record_view *out_record,
                                  int *out_has_record,
                                  int *out_system_errno,
                                  char *error,
                                  size_t error_capacity) {
    bam1_t *record = NULL;
    size_t query_name_length = 0;
    size_t sequence_length = 0;
    ptrdiff_t auxiliary_offset = 0;
    size_t auxiliary_length = 0;
    int result = 0;
    int saved_errno = 0;
    int status = indexed_reader_is_live(reader, out_system_errno, error,
                                        error_capacity);
    if (out_record != NULL) {
        memset(out_record, 0, sizeof(*out_record));
    }
    if (out_has_record != NULL) {
        *out_has_record = 0;
    }
    if (status != BSBIT_HTS_OK) {
        return status;
    }
    if (out_record == NULL || out_has_record == NULL) {
        return set_result(BSBIT_HTS_INVALID_ARGUMENT, 0,
                          "out_record and out_has_record are required",
                          out_system_errno, error, error_capacity);
    }
    if (reader->iterator == NULL) {
        return set_result(BSBIT_HTS_INVALID_STATE, 0,
                          "a region query is required before reading records",
                          out_system_errno, error, error_capacity);
    }
    errno = 0;
    result = sam_itr_next(reader->file, reader->iterator, reader->record);
    saved_errno = errno;
    if (result == -1) {
        return set_result(BSBIT_HTS_OK, 0, NULL, out_system_errno, error,
                          error_capacity);
    }
    if (result < -1) {
        reader->failed = 1;
        return set_result(BSBIT_HTS_READ_FAILED, saved_errno,
                          "HTSlib failed while reading a BAM region",
                          out_system_errno, error, error_capacity);
    }
    record = reader->record;
    if (record->core.l_qseq < 0 ||
        record->core.l_qname <= (uint16_t)record->core.l_extranul) {
        reader->failed = 1;
        return set_result(BSBIT_HTS_RECORD_FAILED, 0,
                          "HTSlib returned invalid BAM record dimensions",
                          out_system_errno, error, error_capacity);
    }
    query_name_length = (size_t)(record->core.l_qname -
                                 (uint16_t)record->core.l_extranul -
                                 UINT16_C(1));
    sequence_length = (size_t)record->core.l_qseq;
    auxiliary_offset = bam_get_aux(record) - record->data;
    if (record->l_data < 0 || auxiliary_offset < 0 ||
        (size_t)auxiliary_offset > (size_t)record->l_data) {
        reader->failed = 1;
        return set_result(BSBIT_HTS_RECORD_FAILED, 0,
                          "HTSlib returned invalid BAM auxiliary dimensions",
                          out_system_errno, error, error_capacity);
    }
    auxiliary_length = (size_t)record->l_data - (size_t)auxiliary_offset;
    out_record->reference_id = record->core.tid;
    out_record->position = (int64_t)record->core.pos;
    out_record->mapping_quality = record->core.qual;
    out_record->flag = record->core.flag;
    out_record->mate_reference_id = record->core.mtid;
    out_record->mate_position = (int64_t)record->core.mpos;
    out_record->template_length = (int64_t)record->core.isize;
    out_record->query_name = bam_get_qname(record);
    out_record->query_name_length = query_name_length;
    out_record->cigar = bam_get_cigar(record);
    out_record->cigar_count = (size_t)record->core.n_cigar;
    out_record->sequence = bam_get_seq(record);
    out_record->packed_sequence_length = (sequence_length + 1U) / 2U;
    out_record->sequence_length = sequence_length;
    out_record->quality = bam_get_qual(record);
    out_record->auxiliary = bam_get_aux(record);
    out_record->auxiliary_length = auxiliary_length;
    *out_has_record = 1;
    return set_result(BSBIT_HTS_OK, 0, NULL, out_system_errno, error,
                      error_capacity);
}

int bsbit_hts_indexed_reader_close(bsbit_hts_indexed_reader *reader,
                                   int *out_system_errno,
                                   char *error,
                                   size_t error_capacity) {
    int result = 0;
    int saved_errno = 0;

    if (reader == NULL) {
        return set_result(BSBIT_HTS_INVALID_ARGUMENT, 0, "reader is null",
                          out_system_errno, error, error_capacity);
    }
    if (reader->closed) {
        return set_result(BSBIT_HTS_INVALID_STATE, 0,
                          "reader is already closed", out_system_errno, error,
                          error_capacity);
    }
    reader->closed = 1;
    errno = 0;
    result = cleanup_indexed_reader(reader);
    saved_errno = errno;
    if (result < 0) {
        return set_result(BSBIT_HTS_CLOSE_FAILED, saved_errno,
                          "HTSlib failed while closing indexed BAM input",
                          out_system_errno, error, error_capacity);
    }
    return set_result(BSBIT_HTS_OK, 0, NULL, out_system_errno, error,
                      error_capacity);
}

void bsbit_hts_indexed_reader_destroy(bsbit_hts_indexed_reader *reader) {
    if (reader == NULL) {
        return;
    }
    (void)cleanup_indexed_reader(reader);
    free(reader);
}

static int indexed_fasta_reader_is_live(
    const bsbit_hts_indexed_fasta_reader *reader,
    int *out_system_errno,
    char *error,
    size_t error_capacity) {
    if (reader == NULL) {
        return set_result(BSBIT_HTS_INVALID_ARGUMENT, 0, "reader is null",
                          out_system_errno, error, error_capacity);
    }
    if (reader->closed || reader->index == NULL) {
        return set_result(BSBIT_HTS_INVALID_STATE, 0,
                          "indexed FASTA reader is closed", out_system_errno,
                          error, error_capacity);
    }
    if (reader->failed) {
        return set_result(BSBIT_HTS_INVALID_STATE, 0,
                          "indexed FASTA reader is terminal after a read failure",
                          out_system_errno, error, error_capacity);
    }
    return BSBIT_HTS_OK;
}

int bsbit_hts_indexed_fasta_reader_open(
    const char *path,
    bsbit_hts_indexed_fasta_reader **out_reader,
    int *out_system_errno,
    char *error,
    size_t error_capacity) {
    bsbit_hts_indexed_fasta_reader *reader = NULL;
    int saved_errno = 0;

    if (out_reader == NULL) {
        return set_result(BSBIT_HTS_INVALID_ARGUMENT, 0, "out_reader is null",
                          out_system_errno, error, error_capacity);
    }
    *out_reader = NULL;
    if (path == NULL || path[0] == '\0') {
        return set_result(BSBIT_HTS_INVALID_ARGUMENT, 0,
                          "a nonempty path is required", out_system_errno,
                          error, error_capacity);
    }
    reader = (bsbit_hts_indexed_fasta_reader *)calloc(1, sizeof(*reader));
    if (reader == NULL) {
        return set_result(BSBIT_HTS_ALLOCATION_FAILED, errno,
                          "allocate indexed FASTA reader", out_system_errno,
                          error, error_capacity);
    }
    errno = 0;
    reader->index = fai_load3(path, NULL, NULL, 0);
    saved_errno = errno;
    if (reader->index == NULL) {
        cleanup_indexed_fasta_reader(reader);
        free(reader);
        return set_result(
            BSBIT_HTS_OPEN_FAILED, saved_errno,
            "HTSlib could not load FASTA with its adjacent .fai/.gzi indexes",
            out_system_errno, error, error_capacity);
    }
    *out_reader = reader;
    return set_result(BSBIT_HTS_OK, 0, NULL, out_system_errno, error,
                      error_capacity);
}

int bsbit_hts_indexed_fasta_reader_reference_count(
    const bsbit_hts_indexed_fasta_reader *reader,
    int32_t *out_count,
    int *out_system_errno,
    char *error,
    size_t error_capacity) {
    int count = 0;
    int status = indexed_fasta_reader_is_live(
        reader, out_system_errno, error, error_capacity);
    if (status != BSBIT_HTS_OK) {
        return status;
    }
    if (out_count == NULL) {
        return set_result(BSBIT_HTS_INVALID_ARGUMENT, 0,
                          "out_count is null", out_system_errno, error,
                          error_capacity);
    }
    count = faidx_nseq(reader->index);
    if (count < 0) {
        return set_result(BSBIT_HTS_HEADER_FAILED, 0,
                          "HTSlib returned a negative FASTA reference count",
                          out_system_errno, error, error_capacity);
    }
    *out_count = (int32_t)count;
    return set_result(BSBIT_HTS_OK, 0, NULL, out_system_errno, error,
                      error_capacity);
}

int bsbit_hts_indexed_fasta_reader_reference(
    const bsbit_hts_indexed_fasta_reader *reader,
    int32_t reference_id,
    const char **out_name,
    size_t *out_name_length,
    int64_t *out_length,
    int *out_system_errno,
    char *error,
    size_t error_capacity) {
    const char *name = NULL;
    hts_pos_t length = -1;
    int count = 0;
    int status = indexed_fasta_reader_is_live(
        reader, out_system_errno, error, error_capacity);
    if (status != BSBIT_HTS_OK) {
        return status;
    }
    if (out_name == NULL || out_name_length == NULL || out_length == NULL) {
        return set_result(BSBIT_HTS_INVALID_ARGUMENT, 0,
                          "indexed FASTA reference output is null",
                          out_system_errno, error, error_capacity);
    }
    count = faidx_nseq(reader->index);
    if (reference_id < 0 || reference_id >= count) {
        return set_result(BSBIT_HTS_INVALID_ARGUMENT, 0,
                          "indexed FASTA reference id is out of range",
                          out_system_errno, error, error_capacity);
    }
    name = faidx_iseq(reader->index, reference_id);
    if (name == NULL) {
        return set_result(BSBIT_HTS_HEADER_FAILED, 0,
                          "HTSlib returned a null FASTA reference name",
                          out_system_errno, error, error_capacity);
    }
    length = faidx_seq_len64(reader->index, name);
    if (length < 0) {
        return set_result(BSBIT_HTS_HEADER_FAILED, 0,
                          "HTSlib returned an invalid FASTA reference length",
                          out_system_errno, error, error_capacity);
    }
    *out_name = name;
    *out_name_length = strlen(name);
    *out_length = (int64_t)length;
    return set_result(BSBIT_HTS_OK, 0, NULL, out_system_errno, error,
                      error_capacity);
}

int bsbit_hts_indexed_fasta_reader_fetch(
    bsbit_hts_indexed_fasta_reader *reader,
    int32_t reference_id,
    int64_t start,
    int64_t end,
    const char **out_sequence,
    size_t *out_length,
    int *out_system_errno,
    char *error,
    size_t error_capacity) {
    const char *name = NULL;
    hts_pos_t reference_length = -1;
    hts_pos_t fetched_length = -1;
    uint64_t requested_length = 0;
    int count = 0;
    int saved_errno = 0;
    int status = indexed_fasta_reader_is_live(
        reader, out_system_errno, error, error_capacity);
    if (status != BSBIT_HTS_OK) {
        return status;
    }
    if (out_sequence == NULL || out_length == NULL) {
        return set_result(BSBIT_HTS_INVALID_ARGUMENT, 0,
                          "indexed FASTA fetch output is null",
                          out_system_errno, error, error_capacity);
    }
    *out_sequence = NULL;
    *out_length = 0;
    count = faidx_nseq(reader->index);
    if (reference_id < 0 || reference_id >= count || start < 0 || end <= start) {
        return set_result(BSBIT_HTS_INVALID_ARGUMENT, 0,
                          "invalid indexed FASTA interval", out_system_errno,
                          error, error_capacity);
    }
    name = faidx_iseq(reader->index, reference_id);
    reference_length = name == NULL ? -1 : faidx_seq_len64(reader->index, name);
    if (reference_length < 0 || end > reference_length) {
        return set_result(BSBIT_HTS_INVALID_ARGUMENT, 0,
                          "indexed FASTA interval exceeds the reference",
                          out_system_errno, error, error_capacity);
    }
    requested_length = (uint64_t)(end - start);
    if (requested_length > (uint64_t)SIZE_MAX) {
        return set_result(BSBIT_HTS_INVALID_ARGUMENT, 0,
                          "indexed FASTA interval exceeds size_t",
                          out_system_errno, error, error_capacity);
    }
    free(reader->sequence);
    reader->sequence = NULL;
    errno = 0;
    reader->sequence = faidx_fetch_seq64(reader->index, name, start, end - 1,
                                         &fetched_length);
    saved_errno = errno;
    if (reader->sequence == NULL || fetched_length < 0 ||
        (uint64_t)fetched_length != requested_length) {
        free(reader->sequence);
        reader->sequence = NULL;
        reader->failed = 1;
        return set_result(BSBIT_HTS_READ_FAILED, saved_errno,
                          "HTSlib failed to fetch the complete FASTA interval",
                          out_system_errno, error, error_capacity);
    }
    *out_sequence = reader->sequence;
    *out_length = (size_t)fetched_length;
    return set_result(BSBIT_HTS_OK, 0, NULL, out_system_errno, error,
                      error_capacity);
}

int bsbit_hts_indexed_fasta_reader_close(
    bsbit_hts_indexed_fasta_reader *reader,
    int *out_system_errno,
    char *error,
    size_t error_capacity) {
    if (reader == NULL) {
        return set_result(BSBIT_HTS_INVALID_ARGUMENT, 0, "reader is null",
                          out_system_errno, error, error_capacity);
    }
    if (reader->closed) {
        return set_result(BSBIT_HTS_INVALID_STATE, 0,
                          "indexed FASTA reader is already closed",
                          out_system_errno, error, error_capacity);
    }
    reader->closed = 1;
    cleanup_indexed_fasta_reader(reader);
    return set_result(BSBIT_HTS_OK, 0, NULL, out_system_errno, error,
                      error_capacity);
}

void bsbit_hts_indexed_fasta_reader_destroy(
    bsbit_hts_indexed_fasta_reader *reader) {
    if (reader == NULL) {
        return;
    }
    cleanup_indexed_fasta_reader(reader);
    free(reader);
}

int bsbit_hts_bgzf_writer_open(const char *path,
                               uint32_t compression_threads,
                               bsbit_hts_bgzf_writer **out_writer,
                               int *out_system_errno,
                               char *error,
                               size_t error_capacity) {
    bsbit_hts_bgzf_writer *writer = NULL;
    int saved_errno = 0;

    if (out_writer == NULL) {
        return set_result(BSBIT_HTS_INVALID_ARGUMENT, 0, "out_writer is null",
                          out_system_errno, error, error_capacity);
    }
    *out_writer = NULL;
    if (path == NULL || path[0] == '\0') {
        return set_result(BSBIT_HTS_INVALID_ARGUMENT, 0,
                          "a nonempty path is required", out_system_errno,
                          error, error_capacity);
    }
    if (compression_threads > UINT32_C(64)) {
        return set_result(BSBIT_HTS_INVALID_ARGUMENT, 0,
                          "compression_threads must be in 0..=64",
                          out_system_errno, error, error_capacity);
    }
    writer = (bsbit_hts_bgzf_writer *)calloc(1, sizeof(*writer));
    if (writer == NULL) {
        return set_result(BSBIT_HTS_ALLOCATION_FAILED, errno,
                          "failed to allocate BGZF writer", out_system_errno,
                          error, error_capacity);
    }
    errno = 0;
    writer->file = bgzf_open(path, "w");
    saved_errno = errno;
    if (writer->file == NULL) {
        free(writer);
        return set_result(BSBIT_HTS_OPEN_FAILED, saved_errno,
                          "HTSlib could not open BGZF output",
                          out_system_errno, error, error_capacity);
    }
    if (compression_threads > 0 &&
        bgzf_mt(writer->file, (int)compression_threads, 0) < 0) {
        (void)bgzf_close(writer->file);
        free(writer);
        return set_result(BSBIT_HTS_OPEN_FAILED, 0,
                          "HTSlib could not create BGZF compression threads",
                          out_system_errno, error, error_capacity);
    }
    *out_writer = writer;
    return set_result(BSBIT_HTS_OK, 0, NULL, out_system_errno, error,
                      error_capacity);
}

int bsbit_hts_bgzf_writer_write(bsbit_hts_bgzf_writer *writer,
                                const uint8_t *data,
                                size_t length,
                                size_t *out_count,
                                int *out_system_errno,
                                char *error,
                                size_t error_capacity) {
    ssize_t count = 0;
    int saved_errno = 0;

    if (out_count != NULL) {
        *out_count = 0;
    }
    if (writer == NULL || out_count == NULL ||
        (length > 0 && data == NULL)) {
        return set_result(BSBIT_HTS_INVALID_ARGUMENT, 0,
                          "writer, out_count, and nonempty data are required",
                          out_system_errno, error, error_capacity);
    }
    if (writer->finished || writer->file == NULL) {
        return set_result(BSBIT_HTS_INVALID_STATE, 0,
                          "BGZF writer is already finished", out_system_errno,
                          error, error_capacity);
    }
    if (writer->failed) {
        return set_result(BSBIT_HTS_WRITE_FAILED, 0,
                          "BGZF writer is terminal after an earlier failure",
                          out_system_errno, error, error_capacity);
    }
    if (length == 0) {
        return set_result(BSBIT_HTS_OK, 0, NULL, out_system_errno, error,
                          error_capacity);
    }
    if (length > (size_t)SSIZE_MAX) {
        return set_result(BSBIT_HTS_INVALID_ARGUMENT, 0,
                          "BGZF write length exceeds ssize_t",
                          out_system_errno, error, error_capacity);
    }
    errno = 0;
    count = bgzf_write(writer->file, data, length);
    saved_errno = errno;
    if (count < 0) {
        writer->failed = 1;
        return set_result(BSBIT_HTS_WRITE_FAILED, saved_errno,
                          "HTSlib failed while encoding BGZF output",
                          out_system_errno, error, error_capacity);
    }
    *out_count = (size_t)count;
    return set_result(BSBIT_HTS_OK, 0, NULL, out_system_errno, error,
                      error_capacity);
}

int bsbit_hts_bgzf_writer_flush(bsbit_hts_bgzf_writer *writer,
                                int *out_system_errno,
                                char *error,
                                size_t error_capacity) {
    int result = 0;
    int saved_errno = 0;

    if (writer == NULL) {
        return set_result(BSBIT_HTS_INVALID_ARGUMENT, 0, "writer is null",
                          out_system_errno, error, error_capacity);
    }
    if (writer->finished || writer->file == NULL) {
        return set_result(BSBIT_HTS_INVALID_STATE, 0,
                          "BGZF writer is already finished", out_system_errno,
                          error, error_capacity);
    }
    if (writer->failed) {
        return set_result(BSBIT_HTS_WRITE_FAILED, 0,
                          "BGZF writer is terminal after an earlier failure",
                          out_system_errno, error, error_capacity);
    }
    errno = 0;
    result = bgzf_flush(writer->file);
    saved_errno = errno;
    if (result < 0) {
        writer->failed = 1;
        return set_result(BSBIT_HTS_WRITE_FAILED, saved_errno,
                          "HTSlib failed while flushing BGZF output",
                          out_system_errno, error, error_capacity);
    }
    return set_result(BSBIT_HTS_OK, 0, NULL, out_system_errno, error,
                      error_capacity);
}

int bsbit_hts_bgzf_writer_finish(bsbit_hts_bgzf_writer *writer,
                                 int *out_system_errno,
                                 char *error,
                                 size_t error_capacity) {
    int result = 0;
    int saved_errno = 0;
    int failed = 0;

    if (writer == NULL) {
        return set_result(BSBIT_HTS_INVALID_ARGUMENT, 0, "writer is null",
                          out_system_errno, error, error_capacity);
    }
    if (writer->finished || writer->file == NULL) {
        return set_result(BSBIT_HTS_INVALID_STATE, 0,
                          "BGZF writer is already finished", out_system_errno,
                          error, error_capacity);
    }
    failed = writer->failed;
    writer->finished = 1;
    errno = 0;
    result = bgzf_close(writer->file);
    saved_errno = errno;
    writer->file = NULL;
    if (result < 0) {
        return set_result(BSBIT_HTS_CLOSE_FAILED, saved_errno,
                          "HTSlib failed while closing BGZF output",
                          out_system_errno, error, error_capacity);
    }
    if (failed) {
        return set_result(BSBIT_HTS_INVALID_STATE, 0,
                          "BGZF writer cannot succeed after a write failure",
                          out_system_errno, error, error_capacity);
    }
    return set_result(BSBIT_HTS_OK, 0, NULL, out_system_errno, error,
                      error_capacity);
}

void bsbit_hts_bgzf_writer_destroy(bsbit_hts_bgzf_writer *writer) {
    if (writer == NULL) {
        return;
    }
    if (writer->file != NULL) {
        (void)bgzf_close(writer->file);
        writer->file = NULL;
    }
    free(writer);
}

int bsbit_hts_writer_open_bam_threads_level(const char *path,
                                            const char *header_text,
                                            size_t header_length,
                                            uint32_t compression_threads,
                                            int compression_level,
                                            bsbit_hts_writer **out_writer,
                                            int *out_system_errno,
                                            char *error,
                                            size_t error_capacity) {
    bsbit_hts_writer *writer = NULL;
    int saved_errno = 0;

    if (out_writer == NULL) {
        return set_result(BSBIT_HTS_INVALID_ARGUMENT, 0, "out_writer is null",
                          out_system_errno, error, error_capacity);
    }
    *out_writer = NULL;
    if (path == NULL || path[0] == '\0' || header_text == NULL ||
        header_length == 0) {
        return set_result(BSBIT_HTS_INVALID_ARGUMENT, 0,
                          "path and nonempty header are required",
                          out_system_errno, error, error_capacity);
    }
    if (compression_threads > UINT32_C(64)) {
        return set_result(BSBIT_HTS_INVALID_ARGUMENT, 0,
                          "compression_threads must be in 0..=64",
                          out_system_errno, error, error_capacity);
    }
    if (compression_level < -1 || compression_level > 9) {
        return set_result(BSBIT_HTS_INVALID_ARGUMENT, 0,
                          "compression_level must be in -1..=9",
                          out_system_errno, error, error_capacity);
    }
    if (memchr(header_text, '\0', header_length) != NULL ||
        header_text[header_length - 1] != '\n') {
        return set_result(BSBIT_HTS_HEADER_FAILED, 0,
                          "header must be NUL-free and newline-terminated",
                          out_system_errno, error, error_capacity);
    }
    writer = (bsbit_hts_writer *)calloc(1, sizeof(*writer));
    if (writer == NULL) {
        return set_result(BSBIT_HTS_ALLOCATION_FAILED, errno,
                          "failed to allocate writer", out_system_errno, error,
                          error_capacity);
    }
    writer->header = sam_hdr_parse(header_length, header_text);
    if (writer->header == NULL) {
        cleanup_writer(writer);
        free(writer);
        return set_result(BSBIT_HTS_HEADER_FAILED, 0,
                          "HTSlib rejected the SAM header", out_system_errno,
                          error, error_capacity);
    }
    writer->record = bam_init1();
    if (writer->record == NULL) {
        cleanup_writer(writer);
        free(writer);
        return set_result(BSBIT_HTS_ALLOCATION_FAILED, errno,
                          "failed to allocate HTSlib record", out_system_errno,
                          error, error_capacity);
    }
    errno = 0;
    writer->file = sam_open(path, "wb");
    saved_errno = errno;
    if (writer->file == NULL) {
        cleanup_writer(writer);
        free(writer);
        return set_result(BSBIT_HTS_OPEN_FAILED, saved_errno,
                          "HTSlib could not open output path", out_system_errno,
                          error, error_capacity);
    }
    if (compression_level >= 0 &&
        hts_set_opt(writer->file, HTS_OPT_COMPRESSION_LEVEL,
                    compression_level) < 0) {
        cleanup_writer(writer);
        free(writer);
        return set_result(BSBIT_HTS_OPEN_FAILED, 0,
                          "HTSlib rejected the BAM compression level",
                          out_system_errno, error, error_capacity);
    }
    if (compression_threads > 0 &&
        hts_set_threads(writer->file, (int)compression_threads) < 0) {
        cleanup_writer(writer);
        free(writer);
        return set_result(BSBIT_HTS_OPEN_FAILED, 0,
                          "HTSlib could not create BAM compression threads",
                          out_system_errno, error, error_capacity);
    }
    errno = 0;
    if (sam_hdr_write(writer->file, writer->header) < 0) {
        saved_errno = errno;
        cleanup_writer(writer);
        free(writer);
        return set_result(BSBIT_HTS_HEADER_FAILED, saved_errno,
                          "HTSlib could not write the BAM header",
                          out_system_errno, error, error_capacity);
    }
    *out_writer = writer;
    return set_result(BSBIT_HTS_OK, 0, NULL, out_system_errno, error,
                      error_capacity);
}

int bsbit_hts_writer_open_bam_threads(const char *path,
                                      const char *header_text,
                                      size_t header_length,
                                      uint32_t compression_threads,
                                      bsbit_hts_writer **out_writer,
                                      int *out_system_errno,
                                      char *error,
                                      size_t error_capacity) {
    return bsbit_hts_writer_open_bam_threads_level(
        path, header_text, header_length, compression_threads, -1, out_writer,
        out_system_errno, error, error_capacity);
}

int bsbit_hts_writer_open_bam(const char *path,
                              const char *header_text,
                              size_t header_length,
                              bsbit_hts_writer **out_writer,
                              int *out_system_errno,
                              char *error,
                              size_t error_capacity) {
    return bsbit_hts_writer_open_bam_threads(
        path, header_text, header_length, UINT32_C(0), out_writer,
        out_system_errno, error, error_capacity);
}

int bsbit_hts_writer_write_record(bsbit_hts_writer *writer,
                                  const char *record_text,
                                  size_t record_length,
                                  int *out_system_errno,
                                  char *error,
                                  size_t error_capacity) {
    char *copy = NULL;
    kstring_t line = KS_INITIALIZE;
    int saved_errno = 0;

    if (writer == NULL) {
        return set_result(BSBIT_HTS_INVALID_ARGUMENT, 0, "writer is null",
                          out_system_errno, error, error_capacity);
    }
    if (writer->finished || writer->file == NULL || writer->failed) {
        return set_result(BSBIT_HTS_INVALID_STATE, 0,
                          "writer is terminal or already finished",
                          out_system_errno, error, error_capacity);
    }
    if (record_text == NULL || record_length < 2) {
        writer->failed = 1;
        return set_result(BSBIT_HTS_INVALID_ARGUMENT, 0,
                          "a nonempty record is required", out_system_errno,
                          error, error_capacity);
    }
    if (record_text[record_length - 1] != '\n' ||
        memchr(record_text, '\0', record_length) != NULL ||
        memchr(record_text, '\r', record_length) != NULL ||
        memchr(record_text, '\n', record_length - 1) != NULL) {
        writer->failed = 1;
        return set_result(BSBIT_HTS_RECORD_FAILED, 0,
                          "record must contain exactly one LF-terminated SAM line",
                          out_system_errno, error, error_capacity);
    }
    copy = (char *)malloc(record_length);
    if (copy == NULL) {
        writer->failed = 1;
        return set_result(BSBIT_HTS_ALLOCATION_FAILED, errno,
                          "failed to allocate SAM record copy",
                          out_system_errno, error, error_capacity);
    }
    memcpy(copy, record_text, record_length - 1);
    copy[record_length - 1] = '\0';
    line.l = record_length - 1;
    line.m = record_length;
    line.s = copy;
    errno = 0;
    if (sam_parse1(&line, writer->header, writer->record) < 0) {
        saved_errno = errno;
        writer->failed = 1;
        ks_free(&line);
        return set_result(BSBIT_HTS_RECORD_FAILED, saved_errno,
                          "HTSlib rejected the SAM record", out_system_errno,
                          error, error_capacity);
    }
    errno = 0;
    if (sam_write1(writer->file, writer->header, writer->record) < 0) {
        saved_errno = errno;
        writer->failed = 1;
        ks_free(&line);
        return set_result(BSBIT_HTS_WRITE_FAILED, saved_errno,
                          "HTSlib could not write the BAM record",
                          out_system_errno, error, error_capacity);
    }
    ks_free(&line);
    return set_result(BSBIT_HTS_OK, 0, NULL, out_system_errno, error,
                      error_capacity);
}

int bsbit_hts_writer_write_bam_fields(bsbit_hts_writer *writer,
                                      const char *query_name,
                                      size_t query_name_length,
                                      uint16_t flag,
                                      int32_t reference_id,
                                      int64_t position,
                                      uint8_t mapping_quality,
                                      const uint32_t *cigar,
                                      size_t cigar_count,
                                      int32_t mate_reference_id,
                                      int64_t mate_position,
                                      int64_t template_length,
                                      const char *sequence,
                                      size_t sequence_length,
                                      const uint8_t *quality,
                                      int has_mapping,
                                      uint32_t literal_nm,
                                      int has_md,
                                      const char *md,
                                      size_t md_length,
                                      int has_xg,
                                      const char *xg,
                                      int has_bismark,
                                      const char *bismark_xm,
                                      size_t bismark_xm_length,
                                      const char *bismark_xr,
                                      int *out_system_errno,
                                      char *error,
                                      size_t error_capacity) {
    size_t integer_bytes = 0;
    size_t auxiliary_bytes = 0;
    size_t bismark_bytes = 0;
    size_t xg_bytes = 0;
    int saved_errno = 0;
    int reference_count = 0;

    if (writer == NULL) {
        return set_result(BSBIT_HTS_INVALID_ARGUMENT, 0, "writer is null",
                          out_system_errno, error, error_capacity);
    }
    if (writer->finished || writer->file == NULL || writer->failed) {
        return set_result(BSBIT_HTS_INVALID_STATE, 0,
                          "writer is terminal or already finished",
                          out_system_errno, error, error_capacity);
    }
    reference_count = sam_hdr_nref(writer->header);
    if (query_name == NULL || query_name_length == 0 ||
        query_name_length > UINT8_MAX - 1 ||
        memchr(query_name, '\0', query_name_length) != NULL ||
        (sequence_length > 0 && sequence == NULL) ||
        (cigar_count > 0 && cigar == NULL) ||
        (has_mapping != 0 && has_mapping != 1) ||
        (has_md != 0 && has_md != 1) ||
        (has_xg != 0 && has_xg != 1) ||
        (has_bismark != 0 && has_bismark != 1) ||
        reference_id < -1 || reference_id >= reference_count ||
        mate_reference_id < -1 || mate_reference_id >= reference_count ||
        position < -1 || mate_position < -1 ||
        (has_mapping != 0 &&
         (reference_id < 0 || position < 0 || cigar_count == 0 ||
          (has_md != 0 &&
           (md == NULL || md_length == 0 || md_length > (size_t)INT_MAX ||
            memchr(md, '\0', md_length) != NULL)) ||
          (has_md == 0 && (md != NULL || md_length != 0)) ||
          (has_xg != 0 &&
           (xg == NULL ||
            !((xg[0] == 'C' && xg[1] == 'T') ||
              (xg[0] == 'G' && xg[1] == 'A')))) ||
          (has_xg == 0 && xg != NULL) ||
          (has_bismark != 0 &&
           (has_md == 0 || has_xg == 0 || bismark_xm == NULL ||
            bismark_xm_length == 0 ||
            bismark_xm_length != sequence_length ||
            bismark_xm_length > (size_t)INT_MAX ||
            memchr(bismark_xm, '\0', bismark_xm_length) != NULL ||
            bismark_xr == NULL ||
            !((bismark_xr[0] == 'C' && bismark_xr[1] == 'T') ||
              (bismark_xr[0] == 'G' && bismark_xr[1] == 'A')))) ||
          (has_bismark == 0 &&
           (bismark_xm != NULL || bismark_xm_length != 0 ||
            bismark_xr != NULL)))) ||
        (has_mapping == 0 &&
         (reference_id != -1 || position != -1 || cigar_count != 0 ||
          has_md != 0 || md != NULL || md_length != 0 ||
          has_xg != 0 || xg != NULL ||
          has_bismark != 0 || bismark_xm != NULL ||
          bismark_xm_length != 0 || bismark_xr != NULL))) {
        writer->failed = 1;
        return set_result(BSBIT_HTS_RECORD_FAILED, 0,
                          "validated BAM fields are inconsistent",
                          out_system_errno, error, error_capacity);
    }

    if (has_mapping != 0) {
        integer_bytes = literal_nm < UINT8_MAX
                            ? (size_t)1
                            : (literal_nm < UINT16_MAX ? (size_t)2 : (size_t)4);
        const size_t md_bytes = has_md != 0 ? md_length + (size_t)4 : 0;
        xg_bytes = has_xg != 0 ? (size_t)6 : 0;
        if (has_bismark != 0) {
            if (bismark_xm_length > SIZE_MAX - (size_t)10) {
                writer->failed = 1;
                return set_result(BSBIT_HTS_RECORD_FAILED, 0,
                                  "BAM auxiliary field length overflowed",
                                  out_system_errno, error, error_capacity);
            }
            bismark_bytes = bismark_xm_length + (size_t)10;
        }
        if (md_bytes > SIZE_MAX - integer_bytes - (size_t)3 ||
            xg_bytes > SIZE_MAX - integer_bytes - (size_t)3 - md_bytes ||
            bismark_bytes >
                SIZE_MAX - integer_bytes - (size_t)3 - md_bytes - xg_bytes) {
            writer->failed = 1;
            return set_result(BSBIT_HTS_RECORD_FAILED, 0,
                              "BAM auxiliary field length overflowed",
                              out_system_errno, error, error_capacity);
        }
        auxiliary_bytes =
            md_bytes + bismark_bytes + xg_bytes + integer_bytes + (size_t)3;
    }

    errno = 0;
    if (bam_set1(writer->record, query_name_length, query_name, flag,
                 reference_id, (hts_pos_t)position, mapping_quality,
                 cigar_count, cigar, mate_reference_id,
                 (hts_pos_t)mate_position, (hts_pos_t)template_length,
                 sequence_length, sequence, (const char *)quality,
                 auxiliary_bytes) < 0) {
        saved_errno = errno;
        writer->failed = 1;
        return set_result(saved_errno == ENOMEM ? BSBIT_HTS_ALLOCATION_FAILED
                                                : BSBIT_HTS_RECORD_FAILED,
                          saved_errno, "HTSlib rejected direct BAM fields",
                          out_system_errno, error, error_capacity);
    }
    if (quality != NULL) {
        uint8_t *encoded_quality = bam_get_qual(writer->record);
        size_t quality_index = 0;
        for (quality_index = 0; quality_index < sequence_length; quality_index++) {
            if (quality[quality_index] < (uint8_t)'!' ||
                quality[quality_index] > (uint8_t)'~') {
                writer->failed = 1;
                return set_result(BSBIT_HTS_RECORD_FAILED, 0,
                                  "direct BAM quality is outside printable Phred+33",
                                  out_system_errno, error, error_capacity);
            }
            encoded_quality[quality_index] =
                quality[quality_index] - (uint8_t)'!';
        }
    }
    if (has_mapping != 0 &&
        (bam_aux_update_int(writer->record, "NM", (int64_t)literal_nm) < 0 ||
         (has_md != 0 &&
          bam_aux_update_str(writer->record, "MD", (int)md_length, md) < 0) ||
         (has_bismark != 0 &&
          (bam_aux_update_str(writer->record, "XM",
                              (int)bismark_xm_length, bismark_xm) < 0 ||
           bam_aux_update_str(writer->record, "XR", 2, bismark_xr) < 0)) ||
         (has_xg != 0 &&
          bam_aux_update_str(writer->record, "XG", 2, xg) < 0))) {
        saved_errno = errno;
        writer->failed = 1;
        return set_result(saved_errno == ENOMEM ? BSBIT_HTS_ALLOCATION_FAILED
                                                : BSBIT_HTS_RECORD_FAILED,
                          saved_errno, "HTSlib rejected direct BAM auxiliary fields",
                          out_system_errno, error, error_capacity);
    }
    errno = 0;
    if (sam_write1(writer->file, writer->header, writer->record) < 0) {
        saved_errno = errno;
        writer->failed = 1;
        return set_result(BSBIT_HTS_WRITE_FAILED, saved_errno,
                          "HTSlib could not write the direct BAM record",
                          out_system_errno, error, error_capacity);
    }
    return set_result(BSBIT_HTS_OK, 0, NULL, out_system_errno, error,
                      error_capacity);
}

int bsbit_hts_writer_finish(bsbit_hts_writer *writer,
                            int *out_system_errno,
                            char *error,
                            size_t error_capacity) {
    int result = 0;
    int saved_errno = 0;
    int failed = 0;

    if (writer == NULL) {
        return set_result(BSBIT_HTS_INVALID_ARGUMENT, 0, "writer is null",
                          out_system_errno, error, error_capacity);
    }
    if (writer->finished || writer->file == NULL) {
        return set_result(BSBIT_HTS_INVALID_STATE, 0,
                          "writer is already finished", out_system_errno, error,
                          error_capacity);
    }
    failed = writer->failed;
    writer->finished = 1;
    errno = 0;
    result = hts_close(writer->file);
    saved_errno = errno;
    writer->file = NULL;
    if (result < 0) {
        return set_result(BSBIT_HTS_CLOSE_FAILED, saved_errno,
                          "HTSlib failed while closing output", out_system_errno,
                          error, error_capacity);
    }
    if (failed) {
        return set_result(BSBIT_HTS_INVALID_STATE, 0,
                          "writer cannot succeed after a record failure",
                          out_system_errno, error, error_capacity);
    }
    return set_result(BSBIT_HTS_OK, 0, NULL, out_system_errno, error,
                      error_capacity);
}

void bsbit_hts_writer_destroy(bsbit_hts_writer *writer) {
    if (writer == NULL) {
        return;
    }
    cleanup_writer(writer);
    free(writer);
}
