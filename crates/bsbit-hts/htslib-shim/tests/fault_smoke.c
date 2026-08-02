#include "bsbit_hts.h"

#include <errno.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include <htslib/bgzf.h>
#include <htslib/hts.h>
#include <htslib/kstring.h>
#include <htslib/sam.h>

static const char HEADER[] =
    "@HD\tVN:1.6\tSO:unknown\n"
    "@SQ\tSN:chr1\tLN:100\n";
static const char RECORD[] =
    "read1\t0\tchr1\t2\t255\t4M\t*\t0\t0\tACGT\tIIII\tNM:i:0\tMD:Z:4\n";
static const uint8_t PAYLOAD[] = "@r1\nACGT\n+\nIIII\n";

enum failpoint {
    FAIL_NONE = 0,
    FAIL_CALLOC,
    FAIL_MALLOC,
    FAIL_BGZF_OPEN,
    FAIL_BGZF_COMPRESSION,
    FAIL_BGZF_READ,
    FAIL_BGZF_WRITE,
    FAIL_BGZF_FLUSH,
    FAIL_BGZF_CLOSE,
    FAIL_SAM_HDR_PARSE,
    FAIL_BAM_INIT,
    FAIL_HTS_OPEN,
    FAIL_SAM_HDR_WRITE,
    FAIL_SAM_PARSE,
    FAIL_SAM_WRITE,
    FAIL_HTS_CLOSE,
    FAILPOINT_COUNT
};

static enum failpoint active_failpoint = FAIL_NONE;
static unsigned int trigger_count = 0;
static unsigned int call_counts[FAILPOINT_COUNT];

#define REQUIRE(condition)                                                        \
    do {                                                                          \
        if (!(condition)) {                                                        \
            fprintf(stderr, "requirement failed at %s:%d: %s\n", __FILE__,       \
                    __LINE__, #condition);                                         \
            goto done;                                                             \
        }                                                                          \
    } while (0)

static void arm(enum failpoint failpoint) {
    memset(call_counts, 0, sizeof(call_counts));
    active_failpoint = failpoint;
    trigger_count = 0;
    errno = 0;
}

static int should_fail(enum failpoint failpoint) {
    call_counts[failpoint] += 1;
    if (active_failpoint != failpoint) {
        return 0;
    }
    active_failpoint = FAIL_NONE;
    trigger_count += 1;
    errno = EIO;
    return 1;
}

void *__real_calloc(size_t count, size_t size);
void *__wrap_calloc(size_t count, size_t size) {
    return should_fail(FAIL_CALLOC) ? NULL : __real_calloc(count, size);
}

void *__real_malloc(size_t size);
void *__wrap_malloc(size_t size) {
    return should_fail(FAIL_MALLOC) ? NULL : __real_malloc(size);
}

BGZF *__real_bgzf_open(const char *path, const char *mode);
BGZF *__wrap_bgzf_open(const char *path, const char *mode) {
    return should_fail(FAIL_BGZF_OPEN) ? NULL : __real_bgzf_open(path, mode);
}

int __real_bgzf_compression(BGZF *file);
int __wrap_bgzf_compression(BGZF *file) {
    return should_fail(FAIL_BGZF_COMPRESSION) ? 99
                                               : __real_bgzf_compression(file);
}

ssize_t __real_bgzf_read(BGZF *file, void *data, size_t length);
ssize_t __wrap_bgzf_read(BGZF *file, void *data, size_t length) {
    return should_fail(FAIL_BGZF_READ) ? -1
                                       : __real_bgzf_read(file, data, length);
}

ssize_t __real_bgzf_write(BGZF *file, const void *data, size_t length);
ssize_t __wrap_bgzf_write(BGZF *file, const void *data, size_t length) {
    return should_fail(FAIL_BGZF_WRITE) ? -1
                                        : __real_bgzf_write(file, data, length);
}

int __real_bgzf_flush(BGZF *file);
int __wrap_bgzf_flush(BGZF *file) {
    return should_fail(FAIL_BGZF_FLUSH) ? -1 : __real_bgzf_flush(file);
}

int __real_bgzf_close(BGZF *file);
int __wrap_bgzf_close(BGZF *file) {
    int result = __real_bgzf_close(file);
    return should_fail(FAIL_BGZF_CLOSE) ? -1 : result;
}

sam_hdr_t *__real_sam_hdr_parse(size_t length, const char *text);
sam_hdr_t *__wrap_sam_hdr_parse(size_t length, const char *text) {
    return should_fail(FAIL_SAM_HDR_PARSE) ? NULL
                                           : __real_sam_hdr_parse(length, text);
}

bam1_t *__real_bam_init1(void);
bam1_t *__wrap_bam_init1(void) {
    return should_fail(FAIL_BAM_INIT) ? NULL : __real_bam_init1();
}

htsFile *__real_hts_open(const char *path, const char *mode);
htsFile *__wrap_hts_open(const char *path, const char *mode) {
    return should_fail(FAIL_HTS_OPEN) ? NULL : __real_hts_open(path, mode);
}

int __real_sam_hdr_write(samFile *file, const sam_hdr_t *header);
int __wrap_sam_hdr_write(samFile *file, const sam_hdr_t *header) {
    return should_fail(FAIL_SAM_HDR_WRITE) ? -1
                                           : __real_sam_hdr_write(file, header);
}

int __real_sam_parse1(kstring_t *line, sam_hdr_t *header, bam1_t *record);
int __wrap_sam_parse1(kstring_t *line, sam_hdr_t *header, bam1_t *record) {
    return should_fail(FAIL_SAM_PARSE) ? -1
                                       : __real_sam_parse1(line, header, record);
}

int __real_sam_write1(samFile *file,
                      const sam_hdr_t *header,
                      const bam1_t *record);
int __wrap_sam_write1(samFile *file,
                      const sam_hdr_t *header,
                      const bam1_t *record) {
    return should_fail(FAIL_SAM_WRITE)
               ? -1
               : __real_sam_write1(file, header, record);
}

int __real_hts_close(htsFile *file);
int __wrap_hts_close(htsFile *file) {
    int result = __real_hts_close(file);
    return should_fail(FAIL_HTS_CLOSE) ? -1 : result;
}

static int write_plain(const char *path) {
    FILE *output = fopen(path, "wb");
    size_t length = sizeof(PAYLOAD) - 1;
    int failed = 0;
    if (output == NULL) {
        return -1;
    }
    failed = fwrite(PAYLOAD, 1, length, output) != length;
    if (fclose(output) != 0) {
        failed = 1;
    }
    return failed ? -1 : 0;
}

static int open_writer(const char *path, bsbit_hts_writer **writer) {
    int system_errno = 0;
    char error[128];
    return bsbit_hts_writer_open_bam(path, HEADER, sizeof(HEADER) - 1, writer,
                                     &system_errno, error, sizeof(error));
}

static int open_bgzf_writer(const char *path,
                            bsbit_hts_bgzf_writer **writer) {
    int system_errno = 0;
    char error[128];
    return bsbit_hts_bgzf_writer_open(path, 0, writer, &system_errno, error,
                                      sizeof(error));
}

int main(int argc, char **argv) {
    int rc = 1;
    char directory[4096] = {0};
    char input_path[4096] = {0};
    char output_path[4096] = {0};
    bsbit_hts_reader *reader = NULL;
    bsbit_hts_bgzf_writer *bgzf_writer = NULL;
    bsbit_hts_writer *writer = NULL;
    uint8_t buffer[16];
    size_t count = 99;
    int system_errno = 0;
    char error[128];

    if (argc != 2) {
        fprintf(stderr, "usage: %s SCRATCH-PREFIX\n", argv[0]);
        return 64;
    }
    REQUIRE(snprintf(directory, sizeof(directory), "%s-%ld", argv[1],
                     (long)getpid()) > 0);
    REQUIRE(mkdir(directory, 0700) == 0);
    REQUIRE(snprintf(input_path, sizeof(input_path), "%s/input.fastq", directory) >
            0);
    REQUIRE(snprintf(output_path, sizeof(output_path), "%s/output.bam", directory) >
            0);
    REQUIRE(write_plain(input_path) == 0);

    arm(FAIL_CALLOC);
    REQUIRE(bsbit_hts_reader_open(input_path, &reader, &system_errno, error,
                                  sizeof(error)) == BSBIT_HTS_ALLOCATION_FAILED);
    REQUIRE(reader == NULL && trigger_count == 1);

    arm(FAIL_BGZF_OPEN);
    REQUIRE(bsbit_hts_reader_open(input_path, &reader, &system_errno, error,
                                  sizeof(error)) == BSBIT_HTS_OPEN_FAILED);
    REQUIRE(reader == NULL && system_errno == EIO && trigger_count == 1);

    arm(FAIL_BGZF_COMPRESSION);
    REQUIRE(bsbit_hts_reader_open(input_path, &reader, &system_errno, error,
                                  sizeof(error)) == BSBIT_HTS_OPEN_FAILED);
    REQUIRE(reader == NULL && trigger_count == 1);
    REQUIRE(call_counts[FAIL_BGZF_CLOSE] == 1);

    arm(FAIL_NONE);
    REQUIRE(bsbit_hts_reader_open(input_path, &reader, &system_errno, error,
                                  sizeof(error)) == BSBIT_HTS_OK);
    arm(FAIL_BGZF_READ);
    REQUIRE(bsbit_hts_reader_read(reader, buffer, sizeof(buffer), &count,
                                  &system_errno, error,
                                  sizeof(error)) == BSBIT_HTS_READ_FAILED);
    REQUIRE(count == 0 && system_errno == EIO && trigger_count == 1);
    REQUIRE(bsbit_hts_reader_read(reader, buffer, sizeof(buffer), &count,
                                  &system_errno, error,
                                  sizeof(error)) == BSBIT_HTS_READ_FAILED);
    REQUIRE(call_counts[FAIL_BGZF_READ] == 1);
    REQUIRE(bsbit_hts_reader_close(reader, &system_errno, error,
                                   sizeof(error)) == BSBIT_HTS_OK);
    REQUIRE(call_counts[FAIL_BGZF_CLOSE] == 1);
    bsbit_hts_reader_destroy(reader);
    reader = NULL;

    arm(FAIL_NONE);
    REQUIRE(bsbit_hts_reader_open(input_path, &reader, &system_errno, error,
                                  sizeof(error)) == BSBIT_HTS_OK);
    arm(FAIL_BGZF_CLOSE);
    REQUIRE(bsbit_hts_reader_close(reader, &system_errno, error,
                                   sizeof(error)) == BSBIT_HTS_CLOSE_FAILED);
    REQUIRE(trigger_count == 1 && call_counts[FAIL_BGZF_CLOSE] == 1);
    bsbit_hts_reader_destroy(reader);
    reader = NULL;
    REQUIRE(call_counts[FAIL_BGZF_CLOSE] == 1);

    arm(FAIL_CALLOC);
    REQUIRE(open_bgzf_writer(output_path, &bgzf_writer) ==
            BSBIT_HTS_ALLOCATION_FAILED);
    REQUIRE(bgzf_writer == NULL && trigger_count == 1);

    arm(FAIL_BGZF_OPEN);
    REQUIRE(open_bgzf_writer(output_path, &bgzf_writer) == BSBIT_HTS_OPEN_FAILED);
    REQUIRE(bgzf_writer == NULL && system_errno == EIO && trigger_count == 1);

    arm(FAIL_NONE);
    REQUIRE(open_bgzf_writer(output_path, &bgzf_writer) == BSBIT_HTS_OK);
    arm(FAIL_BGZF_WRITE);
    REQUIRE(bsbit_hts_bgzf_writer_write(
                bgzf_writer, PAYLOAD, sizeof(PAYLOAD) - 1, &count,
                &system_errno, error, sizeof(error)) == BSBIT_HTS_WRITE_FAILED);
    REQUIRE(count == 0 && system_errno == EIO && trigger_count == 1);
    REQUIRE(bsbit_hts_bgzf_writer_write(
                bgzf_writer, PAYLOAD, sizeof(PAYLOAD) - 1, &count,
                &system_errno, error, sizeof(error)) == BSBIT_HTS_WRITE_FAILED);
    REQUIRE(call_counts[FAIL_BGZF_WRITE] == 1);
    REQUIRE(bsbit_hts_bgzf_writer_finish(bgzf_writer, &system_errno, error,
                                         sizeof(error)) ==
            BSBIT_HTS_INVALID_STATE);
    bsbit_hts_bgzf_writer_destroy(bgzf_writer);
    bgzf_writer = NULL;
    REQUIRE(call_counts[FAIL_BGZF_CLOSE] == 1);
    (void)unlink(output_path);

    arm(FAIL_NONE);
    REQUIRE(open_bgzf_writer(output_path, &bgzf_writer) == BSBIT_HTS_OK);
    arm(FAIL_BGZF_FLUSH);
    REQUIRE(bsbit_hts_bgzf_writer_flush(bgzf_writer, &system_errno, error,
                                        sizeof(error)) ==
            BSBIT_HTS_WRITE_FAILED);
    REQUIRE(system_errno == EIO && trigger_count == 1);
    REQUIRE(bsbit_hts_bgzf_writer_flush(bgzf_writer, &system_errno, error,
                                        sizeof(error)) ==
            BSBIT_HTS_WRITE_FAILED);
    REQUIRE(call_counts[FAIL_BGZF_FLUSH] == 1);
    REQUIRE(bsbit_hts_bgzf_writer_finish(bgzf_writer, &system_errno, error,
                                         sizeof(error)) ==
            BSBIT_HTS_INVALID_STATE);
    bsbit_hts_bgzf_writer_destroy(bgzf_writer);
    bgzf_writer = NULL;
    (void)unlink(output_path);

    arm(FAIL_NONE);
    REQUIRE(open_bgzf_writer(output_path, &bgzf_writer) == BSBIT_HTS_OK);
    arm(FAIL_BGZF_CLOSE);
    REQUIRE(bsbit_hts_bgzf_writer_finish(bgzf_writer, &system_errno, error,
                                         sizeof(error)) ==
            BSBIT_HTS_CLOSE_FAILED);
    REQUIRE(trigger_count == 1 && call_counts[FAIL_BGZF_CLOSE] == 1);
    bsbit_hts_bgzf_writer_destroy(bgzf_writer);
    bgzf_writer = NULL;
    REQUIRE(call_counts[FAIL_BGZF_CLOSE] == 1);
    (void)unlink(output_path);

    arm(FAIL_CALLOC);
    REQUIRE(open_writer(output_path, &writer) == BSBIT_HTS_ALLOCATION_FAILED);
    REQUIRE(writer == NULL && trigger_count == 1);

    arm(FAIL_SAM_HDR_PARSE);
    REQUIRE(open_writer(output_path, &writer) == BSBIT_HTS_HEADER_FAILED);
    REQUIRE(writer == NULL && trigger_count == 1);

    arm(FAIL_BAM_INIT);
    REQUIRE(open_writer(output_path, &writer) == BSBIT_HTS_ALLOCATION_FAILED);
    REQUIRE(writer == NULL && trigger_count == 1);

    arm(FAIL_HTS_OPEN);
    REQUIRE(open_writer(output_path, &writer) == BSBIT_HTS_OPEN_FAILED);
    REQUIRE(writer == NULL && trigger_count == 1);

    arm(FAIL_SAM_HDR_WRITE);
    REQUIRE(open_writer(output_path, &writer) == BSBIT_HTS_HEADER_FAILED);
    REQUIRE(writer == NULL && trigger_count == 1);
    REQUIRE(call_counts[FAIL_HTS_CLOSE] == 1);
    (void)unlink(output_path);

    arm(FAIL_NONE);
    REQUIRE(open_writer(output_path, &writer) == BSBIT_HTS_OK);
    arm(FAIL_MALLOC);
    REQUIRE(bsbit_hts_writer_write_record(writer, RECORD, sizeof(RECORD) - 1,
                                          &system_errno, error,
                                          sizeof(error)) ==
            BSBIT_HTS_ALLOCATION_FAILED);
    REQUIRE(trigger_count == 1);
    REQUIRE(bsbit_hts_writer_finish(writer, &system_errno, error,
                                    sizeof(error)) == BSBIT_HTS_INVALID_STATE);
    bsbit_hts_writer_destroy(writer);
    writer = NULL;
    REQUIRE(call_counts[FAIL_HTS_CLOSE] == 1);
    (void)unlink(output_path);

    arm(FAIL_NONE);
    REQUIRE(open_writer(output_path, &writer) == BSBIT_HTS_OK);
    arm(FAIL_SAM_PARSE);
    REQUIRE(bsbit_hts_writer_write_record(writer, RECORD, sizeof(RECORD) - 1,
                                          &system_errno, error,
                                          sizeof(error)) == BSBIT_HTS_RECORD_FAILED);
    REQUIRE(trigger_count == 1);
    REQUIRE(bsbit_hts_writer_finish(writer, &system_errno, error,
                                    sizeof(error)) == BSBIT_HTS_INVALID_STATE);
    bsbit_hts_writer_destroy(writer);
    writer = NULL;
    (void)unlink(output_path);

    arm(FAIL_NONE);
    REQUIRE(open_writer(output_path, &writer) == BSBIT_HTS_OK);
    arm(FAIL_SAM_WRITE);
    REQUIRE(bsbit_hts_writer_write_record(writer, RECORD, sizeof(RECORD) - 1,
                                          &system_errno, error,
                                          sizeof(error)) == BSBIT_HTS_WRITE_FAILED);
    REQUIRE(trigger_count == 1);
    REQUIRE(bsbit_hts_writer_finish(writer, &system_errno, error,
                                    sizeof(error)) == BSBIT_HTS_INVALID_STATE);
    bsbit_hts_writer_destroy(writer);
    writer = NULL;
    (void)unlink(output_path);

    arm(FAIL_NONE);
    REQUIRE(open_writer(output_path, &writer) == BSBIT_HTS_OK);
    arm(FAIL_HTS_CLOSE);
    REQUIRE(bsbit_hts_writer_finish(writer, &system_errno, error,
                                    sizeof(error)) == BSBIT_HTS_CLOSE_FAILED);
    REQUIRE(trigger_count == 1 && call_counts[FAIL_HTS_CLOSE] == 1);
    bsbit_hts_writer_destroy(writer);
    writer = NULL;
    REQUIRE(call_counts[FAIL_HTS_CLOSE] == 1);
    (void)unlink(output_path);

    rc = 0;

done:
    bsbit_hts_reader_destroy(reader);
    bsbit_hts_bgzf_writer_destroy(bgzf_writer);
    bsbit_hts_writer_destroy(writer);
    (void)unlink(input_path);
    (void)unlink(output_path);
    (void)rmdir(directory);
    if (rc == 0) {
        puts("bsbit_htslib_shim_fault_smoke=PASS");
    }
    return rc;
}
