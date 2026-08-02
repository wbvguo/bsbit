#include "bsbit_hts.h"

#include <errno.h>
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include <htslib/hts.h>

#define BSBIT_NATIVE_FUZZ_READER 1
#define BSBIT_NATIVE_FUZZ_RECORD 2
#define BSBIT_NATIVE_FUZZ_HEADER 3

#ifndef BSBIT_NATIVE_FUZZ_MODE
#error "BSBIT_NATIVE_FUZZ_MODE must select one native fuzz target"
#endif

#if BSBIT_NATIVE_FUZZ_MODE != BSBIT_NATIVE_FUZZ_READER &&                    \
    BSBIT_NATIVE_FUZZ_MODE != BSBIT_NATIVE_FUZZ_RECORD &&                    \
    BSBIT_NATIVE_FUZZ_MODE != BSBIT_NATIVE_FUZZ_HEADER
#error "unsupported BSBIT_NATIVE_FUZZ_MODE"
#endif

#define MAX_INPUT_BYTES 4096u
#define MAX_DECODED_BYTES (1024u * 1024u)
#define READ_BUFFER_BYTES 64u

#if BSBIT_NATIVE_FUZZ_MODE == BSBIT_NATIVE_FUZZ_RECORD
static const char FIXED_HEADER[] =
    "@HD\tVN:1.6\tSO:unknown\n"
    "@SQ\tSN:chr1\tLN:1000\n";
#endif

static char input_path[4096];
static char output_path[4096];
static int paths_initialized = 0;

static void fail_requirement(const char *expression, int line) {
    fprintf(stderr, "native fuzz requirement failed at line %d: %s\n", line,
            expression);
    abort();
}

#define REQUIRE(condition)                                                       \
    do {                                                                         \
        if (!(condition)) {                                                       \
            fail_requirement(#condition, __LINE__);                              \
        }                                                                         \
    } while (0)

static int diagnostic_is_empty(const char *diagnostic, size_t capacity) {
    return capacity > 0u && diagnostic[0] == '\0';
}

static int diagnostic_is_nonempty_terminated(const char *diagnostic,
                                              size_t capacity) {
    return capacity > 1u && diagnostic[0] != '\0' &&
           memchr(diagnostic, '\0', capacity) != NULL;
}

static void initialize_paths(void) {
    const char *scratch = NULL;
    int written = 0;

    if (paths_initialized) {
        return;
    }
    scratch = getenv("BSBIT_NATIVE_FUZZ_SCRATCH");
    REQUIRE(scratch != NULL && scratch[0] != '\0');
    written = snprintf(input_path, sizeof(input_path), "%s/input-%ld.bin",
                       scratch, (long)getpid());
    REQUIRE(written > 0 && (size_t)written < sizeof(input_path));
    written = snprintf(output_path, sizeof(output_path), "%s/output-%ld.bam",
                       scratch, (long)getpid());
    REQUIRE(written > 0 && (size_t)written < sizeof(output_path));
    hts_set_log_level(HTS_LOG_OFF);
    paths_initialized = 1;
}

#if BSBIT_NATIVE_FUZZ_MODE == BSBIT_NATIVE_FUZZ_READER
static void write_input(const uint8_t *data, size_t size) {
    FILE *file = fopen(input_path, "wb");
    REQUIRE(file != NULL);
    if (size > 0u) {
        REQUIRE(fwrite(data, 1u, size, file) == size);
    }
    REQUIRE(fclose(file) == 0);
}

enum reader_terminal {
    READER_EOF = 1,
    READER_FAILED = 2,
    READER_RESOURCE_LIMIT = 3
};

struct reader_outcome {
    int open_status;
    int compression;
    int terminal;
    int close_status;
    size_t decoded_bytes;
};

static struct reader_outcome run_reader_once(size_t input_size) {
    struct reader_outcome outcome = {
        BSBIT_HTS_OK, -1, 0, BSBIT_HTS_OK, 0u,
    };
    bsbit_hts_reader *reader = NULL;
    uint8_t buffer[READ_BUFFER_BYTES];
    char error[256];
    int system_errno = 0;
    size_t iteration = 0u;

    memset(error, 'X', sizeof(error));
    outcome.open_status = bsbit_hts_reader_open(
        input_path, &reader, &system_errno, error, sizeof(error));
    REQUIRE(outcome.open_status == BSBIT_HTS_OK ||
            outcome.open_status == BSBIT_HTS_OPEN_FAILED);
    if (outcome.open_status != BSBIT_HTS_OK) {
        REQUIRE(reader == NULL);
        REQUIRE(diagnostic_is_nonempty_terminated(error, sizeof(error)));
        return outcome;
    }

    REQUIRE(reader != NULL);
    REQUIRE(system_errno == 0);
    REQUIRE(diagnostic_is_empty(error, sizeof(error)));
    REQUIRE(bsbit_hts_reader_compression(reader, &outcome.compression,
                                         &system_errno, error,
                                         sizeof(error)) == BSBIT_HTS_OK);
    REQUIRE(outcome.compression >= BSBIT_HTS_PLAIN &&
            outcome.compression <= BSBIT_HTS_BGZF);

    for (iteration = 0u;
         iteration <= MAX_DECODED_BYTES / READ_BUFFER_BYTES; ++iteration) {
        int status = BSBIT_HTS_OK;
        size_t count = SIZE_MAX;

        memset(error, 'X', sizeof(error));
        system_errno = 0;
        status = bsbit_hts_reader_read(reader, buffer, sizeof(buffer), &count,
                                       &system_errno, error, sizeof(error));
        if (status == BSBIT_HTS_OK) {
            REQUIRE(system_errno == 0);
            REQUIRE(diagnostic_is_empty(error, sizeof(error)));
            REQUIRE(count <= sizeof(buffer));
            if (count == 0u) {
                outcome.terminal = READER_EOF;
                break;
            }
            if (count > MAX_DECODED_BYTES - outcome.decoded_bytes) {
                outcome.terminal = READER_RESOURCE_LIMIT;
                break;
            }
            outcome.decoded_bytes += count;
            if (outcome.decoded_bytes == MAX_DECODED_BYTES) {
                outcome.terminal = READER_RESOURCE_LIMIT;
                break;
            }
            continue;
        }

        REQUIRE(status == BSBIT_HTS_READ_FAILED);
        REQUIRE(count == 0u);
        REQUIRE(diagnostic_is_nonempty_terminated(error, sizeof(error)));
        count = SIZE_MAX;
        REQUIRE(bsbit_hts_reader_read(reader, buffer, sizeof(buffer), &count,
                                      &system_errno, error,
                                      sizeof(error)) == BSBIT_HTS_READ_FAILED);
        REQUIRE(count == 0u);
        REQUIRE(diagnostic_is_nonempty_terminated(error, sizeof(error)));
        outcome.terminal = READER_FAILED;
        break;
    }
    REQUIRE(outcome.terminal != 0);

    memset(error, 'X', sizeof(error));
    outcome.close_status = bsbit_hts_reader_close(
        reader, &system_errno, error, sizeof(error));
    REQUIRE(outcome.close_status == BSBIT_HTS_OK ||
            outcome.close_status == BSBIT_HTS_CLOSE_FAILED);
    if (outcome.close_status == BSBIT_HTS_OK) {
        REQUIRE(diagnostic_is_empty(error, sizeof(error)));
    } else {
        REQUIRE(diagnostic_is_nonempty_terminated(error, sizeof(error)));
    }
    REQUIRE(bsbit_hts_reader_close(reader, &system_errno, error,
                                   sizeof(error)) == BSBIT_HTS_INVALID_STATE);
    REQUIRE(diagnostic_is_nonempty_terminated(error, sizeof(error)));
    bsbit_hts_reader_destroy(reader);

    if (outcome.compression == BSBIT_HTS_PLAIN) {
        REQUIRE(outcome.terminal == READER_EOF);
        REQUIRE(outcome.decoded_bytes == input_size);
    }
    return outcome;
}

static void require_equal_reader(struct reader_outcome left,
                                 struct reader_outcome right) {
    REQUIRE(left.open_status == right.open_status);
    REQUIRE(left.compression == right.compression);
    REQUIRE(left.terminal == right.terminal);
    REQUIRE(left.close_status == right.close_status);
    REQUIRE(left.decoded_bytes == right.decoded_bytes);
}
#endif

#if BSBIT_NATIVE_FUZZ_MODE == BSBIT_NATIVE_FUZZ_RECORD ||                    \
    BSBIT_NATIVE_FUZZ_MODE == BSBIT_NATIVE_FUZZ_HEADER
static uintmax_t output_size(void) {
    struct stat metadata;
    REQUIRE(stat(output_path, &metadata) == 0);
    REQUIRE(metadata.st_size >= 0);
    return (uintmax_t)metadata.st_size;
}

struct writer_outcome {
    int input_status;
    int finish_status;
    uintmax_t file_bytes;
};

#if BSBIT_NATIVE_FUZZ_MODE == BSBIT_NATIVE_FUZZ_RECORD
static struct writer_outcome run_record_once(const uint8_t *data, size_t size) {
    struct writer_outcome outcome = {BSBIT_HTS_OK, BSBIT_HTS_OK, 0u};
    bsbit_hts_writer *writer = NULL;
    const char *record = size == 0u ? "" : (const char *)data;
    char error[256];
    int system_errno = 0;

    (void)unlink(output_path);
    REQUIRE(bsbit_hts_writer_open_bam(
                output_path, FIXED_HEADER, sizeof(FIXED_HEADER) - 1u, &writer,
                &system_errno, error, sizeof(error)) == BSBIT_HTS_OK);
    REQUIRE(writer != NULL);

    memset(error, 'X', sizeof(error));
    outcome.input_status = bsbit_hts_writer_write_record(
        writer, record, size, &system_errno, error, sizeof(error));
    REQUIRE(outcome.input_status == BSBIT_HTS_OK ||
            outcome.input_status == BSBIT_HTS_INVALID_ARGUMENT ||
            outcome.input_status == BSBIT_HTS_RECORD_FAILED);
    if (outcome.input_status == BSBIT_HTS_OK) {
        REQUIRE(system_errno == 0);
        REQUIRE(diagnostic_is_empty(error, sizeof(error)));
    } else {
        REQUIRE(diagnostic_is_nonempty_terminated(error, sizeof(error)));
    }

    outcome.finish_status = bsbit_hts_writer_finish(
        writer, &system_errno, error, sizeof(error));
    REQUIRE(outcome.finish_status ==
            (outcome.input_status == BSBIT_HTS_OK ? BSBIT_HTS_OK
                                                   : BSBIT_HTS_INVALID_STATE));
    REQUIRE(bsbit_hts_writer_finish(writer, &system_errno, error,
                                    sizeof(error)) == BSBIT_HTS_INVALID_STATE);
    bsbit_hts_writer_destroy(writer);
    outcome.file_bytes = output_size();
    REQUIRE(unlink(output_path) == 0);
    return outcome;
}
#endif

#if BSBIT_NATIVE_FUZZ_MODE == BSBIT_NATIVE_FUZZ_HEADER
static struct writer_outcome run_header_once(const uint8_t *data, size_t size) {
    struct writer_outcome outcome = {BSBIT_HTS_OK, BSBIT_HTS_OK, 0u};
    bsbit_hts_writer *writer = NULL;
    const char *header = size == 0u ? "" : (const char *)data;
    char error[256];
    int system_errno = 0;

    (void)unlink(output_path);
    memset(error, 'X', sizeof(error));
    outcome.input_status = bsbit_hts_writer_open_bam(
        output_path, header, size, &writer, &system_errno, error, sizeof(error));
    REQUIRE(outcome.input_status == BSBIT_HTS_OK ||
            outcome.input_status == BSBIT_HTS_INVALID_ARGUMENT ||
            outcome.input_status == BSBIT_HTS_HEADER_FAILED);
    if (outcome.input_status == BSBIT_HTS_OK) {
        REQUIRE(writer != NULL);
        REQUIRE(system_errno == 0);
        REQUIRE(diagnostic_is_empty(error, sizeof(error)));
        outcome.finish_status = bsbit_hts_writer_finish(
            writer, &system_errno, error, sizeof(error));
        REQUIRE(outcome.finish_status == BSBIT_HTS_OK);
        outcome.file_bytes = output_size();
    } else {
        REQUIRE(writer == NULL);
        REQUIRE(diagnostic_is_nonempty_terminated(error, sizeof(error)));
        outcome.finish_status = BSBIT_HTS_INVALID_STATE;
    }
    bsbit_hts_writer_destroy(writer);
    if (outcome.input_status == BSBIT_HTS_OK) {
        REQUIRE(unlink(output_path) == 0);
    } else {
        REQUIRE(unlink(output_path) == 0 || errno == ENOENT);
    }
    return outcome;
}
#endif

static void require_equal_writer(struct writer_outcome left,
                                 struct writer_outcome right) {
    REQUIRE(left.input_status == right.input_status);
    REQUIRE(left.finish_status == right.finish_status);
    REQUIRE(left.file_bytes == right.file_bytes);
}
#endif

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    if (size > MAX_INPUT_BYTES) {
        return 0;
    }
    initialize_paths();

#if BSBIT_NATIVE_FUZZ_MODE == BSBIT_NATIVE_FUZZ_READER
    write_input(data, size);
    require_equal_reader(run_reader_once(size), run_reader_once(size));
#elif BSBIT_NATIVE_FUZZ_MODE == BSBIT_NATIVE_FUZZ_RECORD
    require_equal_writer(run_record_once(data, size), run_record_once(data, size));
#elif BSBIT_NATIVE_FUZZ_MODE == BSBIT_NATIVE_FUZZ_HEADER
    require_equal_writer(run_header_once(data, size), run_header_once(data, size));
#endif
    return 0;
}
