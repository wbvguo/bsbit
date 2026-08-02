#include "bsbit_hts.h"

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include <htslib/bgzf.h>
#include <htslib/faidx.h>
#include <htslib/kstring.h>
#include <htslib/sam.h>
#include <zlib.h>

static const char HEADER[] =
    "@HD\tVN:1.6\tSO:unknown\n"
    "@SQ\tSN:chr1\tLN:100\n";
static const char RECORD[] =
    "read1\t0\tchr1\t2\t255\t4M\t*\t0\t0\tACGT\tIIII\tNM:i:0\tMD:Z:4\n";
static const char RECORD_WITHOUT_LF[] =
    "read1\t0\tchr1\t2\t255\t4M\t*\t0\t0\tACGT\tIIII\tNM:i:0\tMD:Z:4";
static const uint8_t PAYLOAD[] = "@r1\nACGT\n+\nIIII\n";
static const char FASTA[] = ">chr1\nACGTCGTA\n>chr2\nTTAA\n";

#define REQUIRE(condition)                                                        \
    do {                                                                          \
        if (!(condition)) {                                                        \
            fprintf(stderr, "requirement failed at %s:%d: %s\n", __FILE__,       \
                    __LINE__, #condition);                                         \
            goto done;                                                             \
        }                                                                          \
    } while (0)

static int path_join(char *destination,
                     size_t capacity,
                     const char *directory,
                     const char *name) {
    int length = snprintf(destination, capacity, "%s/%s", directory, name);
    return length >= 0 && (size_t)length < capacity ? 0 : -1;
}

static int write_plain(const char *path) {
    FILE *output = fopen(path, "wb");
    size_t length = sizeof(PAYLOAD) - 1;
    int write_failed = 0;
    int close_failed = 0;
    if (output == NULL) {
        return -1;
    }
    write_failed = fwrite(PAYLOAD, 1, length, output) != length;
    close_failed = fclose(output) != 0;
    if (write_failed || close_failed) {
        return -1;
    }
    return 0;
}

static int write_gzip(const char *path) {
    gzFile output = gzopen(path, "wb6");
    unsigned int length = (unsigned int)(sizeof(PAYLOAD) - 1);
    if (output == NULL) {
        return -1;
    }
    if (gzwrite(output, PAYLOAD, length) != (int)length) {
        (void)gzclose(output);
        return -1;
    }
    return gzclose(output) == Z_OK ? 0 : -1;
}

static int write_bgzf(const char *path) {
    BGZF *output = bgzf_open(path, "w");
    size_t length = sizeof(PAYLOAD) - 1;
    if (output == NULL) {
        return -1;
    }
    if (bgzf_write(output, PAYLOAD, length) != (ssize_t)length) {
        (void)bgzf_close(output);
        return -1;
    }
    return bgzf_close(output);
}

static int write_bgzf_with_shim(const char *path) {
    int rc = 1;
    bsbit_hts_bgzf_writer *writer = NULL;
    const size_t split = 5;
    const size_t length = sizeof(PAYLOAD) - 1;
    size_t count = 99;
    int system_errno = -1;
    char error[128];

    REQUIRE(bsbit_hts_bgzf_writer_open(path, 2, &writer, &system_errno, error,
                                       sizeof(error)) == BSBIT_HTS_OK);
    REQUIRE(writer != NULL);
    REQUIRE(bsbit_hts_bgzf_writer_write(writer, PAYLOAD, split, &count,
                                        &system_errno, error,
                                        sizeof(error)) == BSBIT_HTS_OK);
    REQUIRE(count == split);
    REQUIRE(bsbit_hts_bgzf_writer_flush(writer, &system_errno, error,
                                        sizeof(error)) == BSBIT_HTS_OK);
    REQUIRE(bsbit_hts_bgzf_writer_write(writer, PAYLOAD + split, length - split,
                                        &count, &system_errno, error,
                                        sizeof(error)) == BSBIT_HTS_OK);
    REQUIRE(count == length - split);
    REQUIRE(bsbit_hts_bgzf_writer_finish(writer, &system_errno, error,
                                         sizeof(error)) == BSBIT_HTS_OK);
    REQUIRE(bsbit_hts_bgzf_writer_finish(writer, &system_errno, error,
                                         sizeof(error)) ==
            BSBIT_HTS_INVALID_STATE);
    rc = 0;

done:
    bsbit_hts_bgzf_writer_destroy(writer);
    return rc;
}

static int write_indexed_fasta(const char *path) {
    FILE *output = fopen(path, "wb");
    size_t length = sizeof(FASTA) - 1;
    int write_failed = 0;
    int close_failed = 0;
    if (output == NULL) {
        return -1;
    }
    write_failed = fwrite(FASTA, 1, length, output) != length;
    close_failed = fclose(output) != 0;
    if (write_failed || close_failed) {
        return -1;
    }
    return fai_build(path);
}

static int read_indexed_fasta(const char *path) {
    int rc = 1;
    bsbit_hts_indexed_fasta_reader *reader = NULL;
    int32_t count = -1;
    const char *name = NULL;
    size_t name_length = 0;
    int64_t reference_length = -1;
    const char *sequence = NULL;
    size_t sequence_length = 0;
    int system_errno = -1;
    char error[128];

    REQUIRE(bsbit_hts_indexed_fasta_reader_open(
                path, &reader, &system_errno, error, sizeof(error)) ==
            BSBIT_HTS_OK);
    REQUIRE(reader != NULL);
    REQUIRE(bsbit_hts_indexed_fasta_reader_reference_count(
                reader, &count, &system_errno, error, sizeof(error)) ==
            BSBIT_HTS_OK);
    REQUIRE(count == 2);
    REQUIRE(bsbit_hts_indexed_fasta_reader_reference(
                reader, 0, &name, &name_length, &reference_length,
                &system_errno, error, sizeof(error)) == BSBIT_HTS_OK);
    REQUIRE(name_length == 4 && memcmp(name, "chr1", 4) == 0);
    REQUIRE(reference_length == 8);
    REQUIRE(bsbit_hts_indexed_fasta_reader_reference(
                reader, 1, &name, &name_length, &reference_length,
                &system_errno, error, sizeof(error)) == BSBIT_HTS_OK);
    REQUIRE(name_length == 4 && memcmp(name, "chr2", 4) == 0);
    REQUIRE(reference_length == 4);
    REQUIRE(bsbit_hts_indexed_fasta_reader_fetch(
                reader, 0, 1, 6, &sequence, &sequence_length, &system_errno,
                error, sizeof(error)) == BSBIT_HTS_OK);
    REQUIRE(sequence_length == 5 && memcmp(sequence, "CGTCG", 5) == 0);
    REQUIRE(bsbit_hts_indexed_fasta_reader_fetch(
                reader, 1, 0, 4, &sequence, &sequence_length, &system_errno,
                error, sizeof(error)) == BSBIT_HTS_OK);
    REQUIRE(sequence_length == 4 && memcmp(sequence, "TTAA", 4) == 0);
    REQUIRE(bsbit_hts_indexed_fasta_reader_fetch(
                reader, 1, 4, 4, &sequence, &sequence_length, &system_errno,
                error, sizeof(error)) == BSBIT_HTS_INVALID_ARGUMENT);
    REQUIRE(bsbit_hts_indexed_fasta_reader_close(
                reader, &system_errno, error, sizeof(error)) == BSBIT_HTS_OK);
    REQUIRE(bsbit_hts_indexed_fasta_reader_close(
                reader, &system_errno, error, sizeof(error)) ==
            BSBIT_HTS_INVALID_STATE);
    rc = 0;

done:
    bsbit_hts_indexed_fasta_reader_destroy(reader);
    return rc;
}

static int read_source(const char *path, int expected_compression) {
    int rc = 1;
    bsbit_hts_reader *reader = NULL;
    uint8_t decoded[sizeof(PAYLOAD) - 1];
    uint8_t chunk[3];
    size_t decoded_length = 0;
    size_t count = 99;
    int compression = -1;
    int system_errno = -1;
    char error[128] = "not-cleared";

    REQUIRE(bsbit_hts_reader_open(path, &reader, &system_errno, error,
                                  sizeof(error)) == BSBIT_HTS_OK);
    REQUIRE(reader != NULL);
    REQUIRE(system_errno == 0);
    REQUIRE(error[0] == '\0');
    REQUIRE(bsbit_hts_reader_compression(reader, &compression, &system_errno,
                                         error, sizeof(error)) == BSBIT_HTS_OK);
    REQUIRE(compression == expected_compression);
    REQUIRE(bsbit_hts_reader_read(reader, NULL, 0, &count, &system_errno, error,
                                  sizeof(error)) == BSBIT_HTS_OK);
    REQUIRE(count == 0);
    while (decoded_length < sizeof(decoded)) {
        REQUIRE(bsbit_hts_reader_read(reader, chunk, sizeof(chunk), &count,
                                      &system_errno, error,
                                      sizeof(error)) == BSBIT_HTS_OK);
        REQUIRE(count > 0 && count <= sizeof(chunk));
        REQUIRE(decoded_length + count <= sizeof(decoded));
        memcpy(decoded + decoded_length, chunk, count);
        decoded_length += count;
    }
    REQUIRE(memcmp(decoded, PAYLOAD, sizeof(decoded)) == 0);
    REQUIRE(bsbit_hts_reader_read(reader, chunk, sizeof(chunk), &count,
                                  &system_errno, error,
                                  sizeof(error)) == BSBIT_HTS_OK);
    REQUIRE(count == 0);
    REQUIRE(bsbit_hts_reader_read(reader, chunk, sizeof(chunk), &count,
                                  &system_errno, error,
                                  sizeof(error)) == BSBIT_HTS_OK);
    REQUIRE(count == 0);
    REQUIRE(bsbit_hts_reader_close(reader, &system_errno, error,
                                   sizeof(error)) == BSBIT_HTS_OK);
    REQUIRE(bsbit_hts_reader_close(reader, &system_errno, error,
                                   sizeof(error)) == BSBIT_HTS_INVALID_STATE);
    rc = 0;

done:
    bsbit_hts_reader_destroy(reader);
    return rc;
}

static int copy_truncated(const char *source, const char *destination) {
    int rc = -1;
    FILE *input = fopen(source, "rb");
    FILE *output = NULL;
    uint8_t buffer[256];
    size_t count = 0;
    long length = 0;

    if (input == NULL || fseek(input, 0, SEEK_END) != 0) {
        goto done;
    }
    length = ftell(input);
    if (length < 8 || fseek(input, 0, SEEK_SET) != 0) {
        goto done;
    }
    output = fopen(destination, "wb");
    if (output == NULL) {
        goto done;
    }
    while ((count = fread(buffer, 1, sizeof(buffer), input)) > 0) {
        if (fwrite(buffer, 1, count, output) != count) {
            goto done;
        }
    }
    if (ferror(input)) {
        goto done;
    }
    if (fclose(output) != 0) {
        output = NULL;
        goto done;
    }
    output = NULL;
    if (truncate(destination, (off_t)(length - 4)) != 0) {
        goto done;
    }
    rc = 0;

done:
    if (input != NULL) {
        (void)fclose(input);
    }
    if (output != NULL) {
        (void)fclose(output);
    }
    return rc;
}

static int reject_truncated_gzip(const char *path) {
    int rc = 1;
    bsbit_hts_reader *reader = NULL;
    uint8_t buffer[64];
    size_t count = 99;
    int system_errno = 0;
    int status = BSBIT_HTS_OK;
    char error[128];

    REQUIRE(bsbit_hts_reader_open(path, &reader, &system_errno, error,
                                  sizeof(error)) == BSBIT_HTS_OK);
    do {
        status = bsbit_hts_reader_read(reader, buffer, sizeof(buffer), &count,
                                       &system_errno, error, sizeof(error));
    } while (status == BSBIT_HTS_OK && count > 0);
    REQUIRE(status == BSBIT_HTS_READ_FAILED);
    REQUIRE(count == 0);
    REQUIRE(error[0] != '\0');
    REQUIRE(bsbit_hts_reader_read(reader, buffer, sizeof(buffer), &count,
                                  &system_errno, error,
                                  sizeof(error)) == BSBIT_HTS_READ_FAILED);
    REQUIRE(count == 0);
    REQUIRE(bsbit_hts_reader_close(reader, &system_errno, error,
                                   sizeof(error)) == BSBIT_HTS_CLOSE_FAILED);
    REQUIRE(error[0] != '\0');
    rc = 0;

done:
    bsbit_hts_reader_destroy(reader);
    return rc;
}

static int decode_bam(const char *path) {
    int rc = 1;
    samFile *input = NULL;
    sam_hdr_t *header = NULL;
    bam1_t *record = NULL;
    kstring_t rendered = KS_INITIALIZE;

    input = sam_open(path, "rb");
    REQUIRE(input != NULL);
    header = sam_hdr_read(input);
    REQUIRE(header != NULL);
    REQUIRE(strcmp(sam_hdr_str(header), HEADER) == 0);
    record = bam_init1();
    REQUIRE(record != NULL);
    REQUIRE(sam_read1(input, header, record) >= 0);
    REQUIRE(sam_format1(header, record, &rendered) >= 0);
    REQUIRE(rendered.s != NULL);
    REQUIRE(strcmp(rendered.s, RECORD_WITHOUT_LF) == 0);
    REQUIRE(sam_read1(input, header, record) == -1);
    REQUIRE(sam_close(input) == 0);
    input = NULL;
    rc = 0;

done:
    if (input != NULL) {
        (void)sam_close(input);
    }
    bam_destroy1(record);
    sam_hdr_destroy(header);
    ks_free(&rendered);
    return rc;
}

static int write_bam(const char *path) {
    int rc = 1;
    bsbit_hts_writer *writer = NULL;
    int system_errno = -1;
    char error[128];

    REQUIRE(bsbit_hts_writer_open_bam(path, HEADER, sizeof(HEADER) - 1, &writer,
                                      &system_errno, error,
                                      sizeof(error)) == BSBIT_HTS_OK);
    REQUIRE(writer != NULL);
    REQUIRE(bsbit_hts_writer_write_record(writer, RECORD, sizeof(RECORD) - 1,
                                          &system_errno, error,
                                          sizeof(error)) == BSBIT_HTS_OK);
    REQUIRE(bsbit_hts_writer_finish(writer, &system_errno, error,
                                    sizeof(error)) == BSBIT_HTS_OK);
    REQUIRE(bsbit_hts_writer_finish(writer, &system_errno, error,
                                    sizeof(error)) == BSBIT_HTS_INVALID_STATE);
    REQUIRE(decode_bam(path) == 0);
    rc = 0;

done:
    bsbit_hts_writer_destroy(writer);
    return rc;
}

static int writer_failure_is_terminal(const char *path) {
    int rc = 1;
    bsbit_hts_writer *writer = NULL;
    int system_errno = 0;
    char error[128];
    static const char bad_record[] = "bad\nsecond\n";

    REQUIRE(bsbit_hts_writer_open_bam(path, HEADER, sizeof(HEADER) - 1, &writer,
                                      &system_errno, error,
                                      sizeof(error)) == BSBIT_HTS_OK);
    REQUIRE(bsbit_hts_writer_write_record(writer, bad_record,
                                          sizeof(bad_record) - 1, &system_errno,
                                          error,
                                          sizeof(error)) == BSBIT_HTS_RECORD_FAILED);
    REQUIRE(bsbit_hts_writer_write_record(writer, RECORD, sizeof(RECORD) - 1,
                                          &system_errno, error,
                                          sizeof(error)) == BSBIT_HTS_INVALID_STATE);
    REQUIRE(bsbit_hts_writer_finish(writer, &system_errno, error,
                                    sizeof(error)) == BSBIT_HTS_INVALID_STATE);
    rc = 0;

done:
    bsbit_hts_writer_destroy(writer);
    return rc;
}

static int writer_invalid_argument_is_terminal(const char *path) {
    int rc = 1;
    bsbit_hts_writer *writer = NULL;
    int system_errno = 0;
    char error[128];

    REQUIRE(bsbit_hts_writer_open_bam(path, HEADER, sizeof(HEADER) - 1, &writer,
                                      &system_errno, error,
                                      sizeof(error)) == BSBIT_HTS_OK);
    REQUIRE(bsbit_hts_writer_write_record(writer, NULL, 0, &system_errno, error,
                                          sizeof(error)) ==
            BSBIT_HTS_INVALID_ARGUMENT);
    REQUIRE(bsbit_hts_writer_write_record(writer, RECORD, sizeof(RECORD) - 1,
                                          &system_errno, error,
                                          sizeof(error)) == BSBIT_HTS_INVALID_STATE);
    REQUIRE(bsbit_hts_writer_finish(writer, &system_errno, error,
                                    sizeof(error)) == BSBIT_HTS_INVALID_STATE);
    rc = 0;

done:
    bsbit_hts_writer_destroy(writer);
    return rc;
}

int main(int argc, char **argv) {
    int rc = 1;
    char directory[4096] = {0};
    char plain_path[4096] = {0};
    char gzip_path[4096] = {0};
    char bgzf_path[4096] = {0};
    char shim_bgzf_path[4096] = {0};
    char truncated_path[4096] = {0};
    char bam_path[4096] = {0};
    char failed_bam_path[4096] = {0};
    char invalid_bam_path[4096] = {0};
    char fasta_path[4096] = {0};
    char fasta_index_path[4096] = {0};
    char missing_path[4096] = {0};
    char missing_output_path[4096] = {0};
    bsbit_hts_reader *reader = (bsbit_hts_reader *)(uintptr_t)1;
    bsbit_hts_writer *writer = (bsbit_hts_writer *)(uintptr_t)1;
    bsbit_hts_indexed_fasta_reader *fasta_reader =
        (bsbit_hts_indexed_fasta_reader *)(uintptr_t)1;
    int system_errno = 0;
    char error[128];

    if (argc != 2) {
        fprintf(stderr, "usage: %s SCRATCH-PREFIX\n", argv[0]);
        return 64;
    }
    REQUIRE(snprintf(directory, sizeof(directory), "%s-%ld", argv[1],
                     (long)getpid()) > 0);
    REQUIRE(mkdir(directory, 0700) == 0);
    REQUIRE(path_join(plain_path, sizeof(plain_path), directory,
                      "plain-with-gz-suffix.gz") == 0);
    REQUIRE(path_join(gzip_path, sizeof(gzip_path), directory,
                      "gzip-with-text-suffix.txt") == 0);
    REQUIRE(path_join(bgzf_path, sizeof(bgzf_path), directory,
                      "bgzf-with-neutral-suffix.data") == 0);
    REQUIRE(path_join(shim_bgzf_path, sizeof(shim_bgzf_path), directory,
                      "shim-bgzf-with-neutral-suffix.data") == 0);
    REQUIRE(path_join(truncated_path, sizeof(truncated_path), directory,
                      "truncated.gz") == 0);
    REQUIRE(path_join(bam_path, sizeof(bam_path), directory, "accepted.bam") == 0);
    REQUIRE(path_join(failed_bam_path, sizeof(failed_bam_path), directory,
                      "failed.bam") == 0);
    REQUIRE(path_join(invalid_bam_path, sizeof(invalid_bam_path), directory,
                      "invalid-argument.bam") == 0);
    REQUIRE(path_join(fasta_path, sizeof(fasta_path), directory,
                      "reference.fa") == 0);
    REQUIRE(snprintf(fasta_index_path, sizeof(fasta_index_path), "%s.fai",
                     fasta_path) > 0);
    REQUIRE(path_join(missing_path, sizeof(missing_path), directory,
                      "missing.fastq") == 0);
    REQUIRE(path_join(missing_output_path, sizeof(missing_output_path), directory,
                      "missing/output.bam") == 0);

    REQUIRE(bsbit_hts_shim_abi_version() == 3);
    REQUIRE(strcmp(bsbit_hts_runtime_version(), "1.24") == 0);
    REQUIRE(bsbit_hts_health_check() == BSBIT_HTS_OK);

    REQUIRE(bsbit_hts_reader_open(NULL, &reader, &system_errno, error,
                                  sizeof(error)) == BSBIT_HTS_INVALID_ARGUMENT);
    REQUIRE(reader == NULL);
    REQUIRE(bsbit_hts_reader_open(missing_path, &reader, &system_errno, error,
                                  sizeof(error)) == BSBIT_HTS_OPEN_FAILED);
    REQUIRE(reader == NULL);
    REQUIRE(error[0] != '\0');

    REQUIRE(bsbit_hts_indexed_fasta_reader_open(
                NULL, &fasta_reader, &system_errno, error, sizeof(error)) ==
            BSBIT_HTS_INVALID_ARGUMENT);
    REQUIRE(fasta_reader == NULL);
    REQUIRE(write_indexed_fasta(fasta_path) == 0);
    REQUIRE(read_indexed_fasta(fasta_path) == 0);

    REQUIRE(write_plain(plain_path) == 0);
    REQUIRE(write_gzip(gzip_path) == 0);
    REQUIRE(write_bgzf(bgzf_path) == 0);
    REQUIRE(write_bgzf_with_shim(shim_bgzf_path) == 0);
    REQUIRE(read_source(plain_path, BSBIT_HTS_PLAIN) == 0);
    REQUIRE(read_source(gzip_path, BSBIT_HTS_GZIP) == 0);
    REQUIRE(read_source(bgzf_path, BSBIT_HTS_BGZF) == 0);
    REQUIRE(read_source(shim_bgzf_path, BSBIT_HTS_BGZF) == 0);
    REQUIRE(copy_truncated(gzip_path, truncated_path) == 0);
    REQUIRE(reject_truncated_gzip(truncated_path) == 0);

    REQUIRE(write_bam(bam_path) == 0);
    REQUIRE(writer_failure_is_terminal(failed_bam_path) == 0);
    REQUIRE(writer_invalid_argument_is_terminal(invalid_bam_path) == 0);
    REQUIRE(bsbit_hts_writer_open_bam(missing_path, "bad", 3, &writer,
                                      &system_errno, error,
                                      sizeof(error)) == BSBIT_HTS_HEADER_FAILED);
    REQUIRE(writer == NULL);
    REQUIRE(bsbit_hts_writer_open_bam(missing_output_path, HEADER,
                                      sizeof(HEADER) - 1, &writer, &system_errno,
                                      error,
                                      sizeof(error)) == BSBIT_HTS_OPEN_FAILED);
    REQUIRE(writer == NULL);

    rc = 0;

done:
    bsbit_hts_reader_destroy(reader == (bsbit_hts_reader *)(uintptr_t)1 ? NULL
                                                                        : reader);
    bsbit_hts_writer_destroy(writer == (bsbit_hts_writer *)(uintptr_t)1 ? NULL
                                                                        : writer);
    bsbit_hts_indexed_fasta_reader_destroy(
        fasta_reader == (bsbit_hts_indexed_fasta_reader *)(uintptr_t)1
            ? NULL
            : fasta_reader);
    (void)unlink(plain_path);
    (void)unlink(gzip_path);
    (void)unlink(bgzf_path);
    (void)unlink(shim_bgzf_path);
    (void)unlink(truncated_path);
    (void)unlink(bam_path);
    (void)unlink(failed_bam_path);
    (void)unlink(invalid_bam_path);
    (void)unlink(fasta_index_path);
    (void)unlink(fasta_path);
    (void)rmdir(directory);
    if (rc == 0) {
        puts("bsbit_htslib_shim_smoke=PASS");
    }
    return rc;
}
