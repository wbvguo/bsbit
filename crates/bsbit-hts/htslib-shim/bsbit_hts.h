#ifndef BSBIT_HTS_H
#define BSBIT_HTS_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct bsbit_hts_reader bsbit_hts_reader;
typedef struct bsbit_hts_indexed_reader bsbit_hts_indexed_reader;
typedef struct bsbit_hts_indexed_fasta_reader bsbit_hts_indexed_fasta_reader;
typedef struct bsbit_hts_bgzf_writer bsbit_hts_bgzf_writer;
typedef struct bsbit_hts_writer bsbit_hts_writer;

typedef struct bsbit_hts_bam_record_view {
    int32_t reference_id;
    int64_t position;
    uint8_t mapping_quality;
    uint16_t flag;
    int32_t mate_reference_id;
    int64_t mate_position;
    int64_t template_length;
    const char *query_name;
    size_t query_name_length;
    const uint32_t *cigar;
    size_t cigar_count;
    const uint8_t *sequence;
    size_t packed_sequence_length;
    size_t sequence_length;
    const uint8_t *quality;
    const uint8_t *auxiliary;
    size_t auxiliary_length;
} bsbit_hts_bam_record_view;

enum bsbit_hts_status {
    BSBIT_HTS_OK = 0,
    BSBIT_HTS_INVALID_ARGUMENT = 1,
    BSBIT_HTS_ALLOCATION_FAILED = 2,
    BSBIT_HTS_OPEN_FAILED = 3,
    BSBIT_HTS_HEADER_FAILED = 4,
    BSBIT_HTS_RECORD_FAILED = 5,
    BSBIT_HTS_CLOSE_FAILED = 6,
    BSBIT_HTS_INVALID_STATE = 7,
    BSBIT_HTS_READ_FAILED = 8,
    BSBIT_HTS_WRITE_FAILED = 9
};

enum bsbit_hts_compression {
    BSBIT_HTS_PLAIN = 0,
    BSBIT_HTS_GZIP = 1,
    BSBIT_HTS_BGZF = 2
};

enum bsbit_hts_tabix_preset {
    BSBIT_HTS_TABIX_VCF = 0,
    BSBIT_HTS_TABIX_BED = 1
};

uint32_t bsbit_hts_shim_abi_version(void);
const char *bsbit_hts_runtime_version(void);
int bsbit_hts_health_check(void);

int bsbit_hts_tabix_index_build(const char *path,
                                const char *index_path,
                                int preset,
                                uint32_t threads,
                                int *out_system_errno,
                                char *error,
                                size_t error_capacity);

int bsbit_hts_bam_index_build(const char *path,
                              const char *index_path,
                              uint32_t threads,
                              int *out_system_errno,
                              char *error,
                              size_t error_capacity);

int bsbit_hts_reader_open(const char *path,
                          bsbit_hts_reader **out_reader,
                          int *out_system_errno,
                          char *error,
                          size_t error_capacity);

int bsbit_hts_reader_compression(const bsbit_hts_reader *reader,
                                 int *out_compression,
                                 int *out_system_errno,
                                 char *error,
                                 size_t error_capacity);

int bsbit_hts_reader_read(bsbit_hts_reader *reader,
                          uint8_t *buffer,
                          size_t capacity,
                          size_t *out_count,
                          int *out_system_errno,
                          char *error,
                          size_t error_capacity);

int bsbit_hts_reader_close(bsbit_hts_reader *reader,
                           int *out_system_errno,
                           char *error,
                           size_t error_capacity);

void bsbit_hts_reader_destroy(bsbit_hts_reader *reader);

int bsbit_hts_indexed_reader_open(const char *path,
                                  bsbit_hts_indexed_reader **out_reader,
                                  int *out_system_errno,
                                  char *error,
                                  size_t error_capacity);

int bsbit_hts_indexed_reader_header_text(
    const bsbit_hts_indexed_reader *reader,
    const char **out_text,
    size_t *out_length,
    int *out_system_errno,
    char *error,
    size_t error_capacity);

int bsbit_hts_indexed_reader_reference_count(
    const bsbit_hts_indexed_reader *reader,
    int32_t *out_count,
    int *out_system_errno,
    char *error,
    size_t error_capacity);

int bsbit_hts_indexed_reader_reference(
    const bsbit_hts_indexed_reader *reader,
    int32_t reference_id,
    const char **out_name,
    size_t *out_name_length,
    int64_t *out_length,
    int *out_system_errno,
    char *error,
    size_t error_capacity);

int bsbit_hts_indexed_reader_query(bsbit_hts_indexed_reader *reader,
                                   int32_t reference_id,
                                   int64_t start,
                                   int64_t end,
                                   int *out_system_errno,
                                   char *error,
                                   size_t error_capacity);

int bsbit_hts_indexed_reader_next(bsbit_hts_indexed_reader *reader,
                                  bsbit_hts_bam_record_view *out_record,
                                  int *out_has_record,
                                  int *out_system_errno,
                                  char *error,
                                  size_t error_capacity);

int bsbit_hts_indexed_reader_close(bsbit_hts_indexed_reader *reader,
                                   int *out_system_errno,
                                   char *error,
                                   size_t error_capacity);

void bsbit_hts_indexed_reader_destroy(bsbit_hts_indexed_reader *reader);

int bsbit_hts_indexed_fasta_reader_open(
    const char *path,
    bsbit_hts_indexed_fasta_reader **out_reader,
    int *out_system_errno,
    char *error,
    size_t error_capacity);

int bsbit_hts_indexed_fasta_reader_reference_count(
    const bsbit_hts_indexed_fasta_reader *reader,
    int32_t *out_count,
    int *out_system_errno,
    char *error,
    size_t error_capacity);

int bsbit_hts_indexed_fasta_reader_reference(
    const bsbit_hts_indexed_fasta_reader *reader,
    int32_t reference_id,
    const char **out_name,
    size_t *out_name_length,
    int64_t *out_length,
    int *out_system_errno,
    char *error,
    size_t error_capacity);

int bsbit_hts_indexed_fasta_reader_fetch(
    bsbit_hts_indexed_fasta_reader *reader,
    int32_t reference_id,
    int64_t start,
    int64_t end,
    const char **out_sequence,
    size_t *out_length,
    int *out_system_errno,
    char *error,
    size_t error_capacity);

int bsbit_hts_indexed_fasta_reader_close(
    bsbit_hts_indexed_fasta_reader *reader,
    int *out_system_errno,
    char *error,
    size_t error_capacity);

void bsbit_hts_indexed_fasta_reader_destroy(
    bsbit_hts_indexed_fasta_reader *reader);

int bsbit_hts_bgzf_writer_open(const char *path,
                               uint32_t compression_threads,
                               bsbit_hts_bgzf_writer **out_writer,
                               int *out_system_errno,
                               char *error,
                               size_t error_capacity);

int bsbit_hts_bgzf_writer_write(bsbit_hts_bgzf_writer *writer,
                                const uint8_t *data,
                                size_t length,
                                size_t *out_count,
                                int *out_system_errno,
                                char *error,
                                size_t error_capacity);

int bsbit_hts_bgzf_writer_flush(bsbit_hts_bgzf_writer *writer,
                                int *out_system_errno,
                                char *error,
                                size_t error_capacity);

int bsbit_hts_bgzf_writer_finish(bsbit_hts_bgzf_writer *writer,
                                 int *out_system_errno,
                                 char *error,
                                 size_t error_capacity);

void bsbit_hts_bgzf_writer_destroy(bsbit_hts_bgzf_writer *writer);

int bsbit_hts_writer_open_bam(const char *path,
                              const char *header_text,
                              size_t header_length,
                              bsbit_hts_writer **out_writer,
                              int *out_system_errno,
                              char *error,
                              size_t error_capacity);

int bsbit_hts_writer_open_bam_threads(const char *path,
                                      const char *header_text,
                                      size_t header_length,
                                      uint32_t compression_threads,
                                      bsbit_hts_writer **out_writer,
                                      int *out_system_errno,
                                      char *error,
                                      size_t error_capacity);

int bsbit_hts_writer_open_bam_threads_level(const char *path,
                                            const char *header_text,
                                            size_t header_length,
                                            uint32_t compression_threads,
                                            int compression_level,
                                            bsbit_hts_writer **out_writer,
                                            int *out_system_errno,
                                            char *error,
                                            size_t error_capacity);

int bsbit_hts_writer_write_record(bsbit_hts_writer *writer,
                                  const char *record_text,
                                  size_t record_length,
                                  int *out_system_errno,
                                  char *error,
                                  size_t error_capacity);

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
                                      size_t error_capacity);

int bsbit_hts_writer_finish(bsbit_hts_writer *writer,
                            int *out_system_errno,
                            char *error,
                            size_t error_capacity);

void bsbit_hts_writer_destroy(bsbit_hts_writer *writer);

#ifdef __cplusplus
}
#endif

#endif
