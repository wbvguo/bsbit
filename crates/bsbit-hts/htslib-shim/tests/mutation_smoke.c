#include "bsbit_hts.h"

#include <errno.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include <htslib/hts.h>
#include <zlib.h>

#ifndef READER_CASES
#define READER_CASES 512u
#endif
#ifndef RECORD_CASES
#define RECORD_CASES 512u
#endif
#ifndef HEADER_CASES
#define HEADER_CASES 256u
#endif
#define MAX_FIXTURE_BYTES 4096u
#define MAX_DECODED_BYTES (1024u * 1024u)

static const char HEADER[] =
    "@HD\tVN:1.6\tSO:unknown\n"
    "@SQ\tSN:chr1\tLN:100\n";
static const char RECORD[] =
    "read1\t0\tchr1\t2\t255\t4M\t*\t0\t0\tACGT\tIIII\tNM:i:0\tMD:Z:4\n";
static const uint8_t PAYLOAD[] = "@r1\nACGT\n+\nIIII\n";

#define REQUIRE(condition)                                                        \
    do {                                                                          \
        if (!(condition)) {                                                        \
            fprintf(stderr, "requirement failed at %s:%d: %s\n", __FILE__,       \
                    __LINE__, #condition);                                         \
            goto done;                                                             \
        }                                                                          \
    } while (0)

static int format_path(char *destination,
                       size_t capacity,
                       const char *directory,
                       const char *name) {
    int written = snprintf(destination, capacity, "%s/%s", directory, name);
    return written > 0 && (size_t)written < capacity ? 0 : -1;
}

static int write_bytes(const char *path, const uint8_t *bytes, size_t length) {
    FILE *output = fopen(path, "wb");
    int failed = 0;
    if (output == NULL) {
        return -1;
    }
    if (length > 0 && fwrite(bytes, 1, length, output) != length) {
        failed = 1;
    }
    if (fclose(output) != 0) {
        failed = 1;
    }
    return failed ? -1 : 0;
}

static int write_gzip(const char *path) {
    gzFile output = gzopen(path, "wb6");
    unsigned int length = (unsigned int)(sizeof(PAYLOAD) - 1u);
    if (output == NULL) {
        return -1;
    }
    if (gzwrite(output, PAYLOAD, length) != (int)length) {
        (void)gzclose(output);
        return -1;
    }
    return gzclose(output) == Z_OK ? 0 : -1;
}

static int load_file(const char *path,
                     uint8_t *bytes,
                     size_t capacity,
                     size_t *out_length) {
    FILE *input = fopen(path, "rb");
    long length = 0;
    int failed = 0;
    if (input == NULL || out_length == NULL) {
        if (input != NULL) {
            (void)fclose(input);
        }
        return -1;
    }
    if (fseek(input, 0, SEEK_END) != 0) {
        failed = 1;
    }
    if (!failed) {
        length = ftell(input);
        if (length < 0 || (uintmax_t)length > (uintmax_t)capacity) {
            failed = 1;
        }
    }
    if (!failed && fseek(input, 0, SEEK_SET) != 0) {
        failed = 1;
    }
    if (!failed && length > 0 &&
        fread(bytes, 1, (size_t)length, input) != (size_t)length) {
        failed = 1;
    }
    if (fclose(input) != 0) {
        failed = 1;
    }
    if (failed) {
        return -1;
    }
    *out_length = (size_t)length;
    return 0;
}

static uint32_t next_random(uint32_t *state) {
    uint32_t value = *state;
    value ^= value << 13;
    value ^= value >> 17;
    value ^= value << 5;
    *state = value;
    return value;
}

static void mutate_bytes(uint8_t *bytes, size_t length, uint32_t seed) {
    uint32_t state = seed;
    unsigned int mutations = 1u + next_random(&state) % 3u;
    unsigned int index = 0;
    if (length == 0) {
        return;
    }
    for (index = 0; index < mutations; ++index) {
        size_t position = (size_t)next_random(&state) % length;
        unsigned int bit = next_random(&state) % 8u;
        bytes[position] ^= (uint8_t)(1u << bit);
    }
}

static int exercise_reader(const char *path) {
    int rc = -1;
    int status = BSBIT_HTS_OK;
    int close_status = BSBIT_HTS_OK;
    int system_errno = INT_MAX;
    int compression = -1;
    int terminated = 0;
    char error[128];
    uint8_t buffer[64];
    size_t count = SIZE_MAX;
    size_t decoded = 0;
    size_t iteration = 0;
    bsbit_hts_reader *reader = NULL;

    memset(error, 'X', sizeof(error));
    status = bsbit_hts_reader_open(path, &reader, &system_errno, error,
                                   sizeof(error));
    if (status != BSBIT_HTS_OK) {
        REQUIRE(status == BSBIT_HTS_OPEN_FAILED);
        REQUIRE(reader == NULL);
        REQUIRE(error[0] != '\0');
        rc = 0;
        goto done;
    }
    REQUIRE(reader != NULL);
    REQUIRE(system_errno == 0 && error[0] == '\0');
    REQUIRE(bsbit_hts_reader_compression(reader, &compression, &system_errno,
                                         error, sizeof(error)) == BSBIT_HTS_OK);
    REQUIRE(compression >= BSBIT_HTS_PLAIN && compression <= BSBIT_HTS_BGZF);

    for (iteration = 0; iteration < MAX_DECODED_BYTES / sizeof(buffer); ++iteration) {
        count = SIZE_MAX;
        system_errno = INT_MAX;
        memset(error, 'X', sizeof(error));
        status = bsbit_hts_reader_read(reader, buffer, sizeof(buffer), &count,
                                       &system_errno, error, sizeof(error));
        if (status == BSBIT_HTS_OK) {
            REQUIRE(system_errno == 0 && error[0] == '\0');
            REQUIRE(count <= sizeof(buffer));
            if (count == 0) {
                terminated = 1;
                break;
            }
            decoded += count;
            REQUIRE(decoded <= MAX_DECODED_BYTES);
            continue;
        }

        REQUIRE(status == BSBIT_HTS_READ_FAILED);
        REQUIRE(count == 0 && error[0] != '\0');
        count = SIZE_MAX;
        REQUIRE(bsbit_hts_reader_read(reader, buffer, sizeof(buffer), &count,
                                      &system_errno, error,
                                      sizeof(error)) == BSBIT_HTS_READ_FAILED);
        REQUIRE(count == 0 && error[0] != '\0');
        terminated = 1;
        break;
    }
    REQUIRE(terminated);

    close_status = bsbit_hts_reader_close(reader, &system_errno, error,
                                          sizeof(error));
    REQUIRE(close_status == BSBIT_HTS_OK ||
            close_status == BSBIT_HTS_CLOSE_FAILED);
    REQUIRE(bsbit_hts_reader_close(reader, &system_errno, error,
                                   sizeof(error)) == BSBIT_HTS_INVALID_STATE);
    bsbit_hts_reader_destroy(reader);
    reader = NULL;
    rc = 0;

done:
    bsbit_hts_reader_destroy(reader);
    return rc;
}

static int exercise_record(const char *path,
                           const uint8_t *record,
                           size_t record_length) {
    int rc = -1;
    int status = BSBIT_HTS_OK;
    int finish_status = BSBIT_HTS_OK;
    int system_errno = INT_MAX;
    char error[128];
    bsbit_hts_writer *writer = NULL;

    REQUIRE(bsbit_hts_writer_open_bam(path, HEADER, sizeof(HEADER) - 1u,
                                      &writer, &system_errno, error,
                                      sizeof(error)) == BSBIT_HTS_OK);
    REQUIRE(writer != NULL);
    system_errno = INT_MAX;
    memset(error, 'X', sizeof(error));
    status = bsbit_hts_writer_write_record(
        writer, (const char *)record, record_length, &system_errno, error,
        sizeof(error));
    REQUIRE(status == BSBIT_HTS_OK || status == BSBIT_HTS_INVALID_ARGUMENT ||
            status == BSBIT_HTS_RECORD_FAILED);
    if (status == BSBIT_HTS_OK) {
        REQUIRE(system_errno == 0 && error[0] == '\0');
    } else {
        REQUIRE(error[0] != '\0');
    }
    finish_status = bsbit_hts_writer_finish(writer, &system_errno, error,
                                            sizeof(error));
    REQUIRE(finish_status == (status == BSBIT_HTS_OK ? BSBIT_HTS_OK
                                                      : BSBIT_HTS_INVALID_STATE));
    REQUIRE(bsbit_hts_writer_finish(writer, &system_errno, error,
                                    sizeof(error)) == BSBIT_HTS_INVALID_STATE);
    bsbit_hts_writer_destroy(writer);
    writer = NULL;
    rc = 0;

done:
    bsbit_hts_writer_destroy(writer);
    (void)unlink(path);
    return rc;
}

static int exercise_header(const char *path,
                           const uint8_t *header,
                           size_t header_length) {
    int rc = -1;
    int status = BSBIT_HTS_OK;
    int system_errno = INT_MAX;
    char error[128];
    bsbit_hts_writer *writer = NULL;

    memset(error, 'X', sizeof(error));
    status = bsbit_hts_writer_open_bam(
        path, (const char *)header, header_length, &writer, &system_errno, error,
        sizeof(error));
    REQUIRE(status == BSBIT_HTS_OK || status == BSBIT_HTS_INVALID_ARGUMENT ||
            status == BSBIT_HTS_HEADER_FAILED);
    if (status == BSBIT_HTS_OK) {
        REQUIRE(writer != NULL);
        REQUIRE(system_errno == 0 && error[0] == '\0');
        REQUIRE(bsbit_hts_writer_finish(writer, &system_errno, error,
                                        sizeof(error)) == BSBIT_HTS_OK);
    } else {
        REQUIRE(writer == NULL && error[0] != '\0');
    }
    bsbit_hts_writer_destroy(writer);
    writer = NULL;
    rc = 0;

done:
    bsbit_hts_writer_destroy(writer);
    (void)unlink(path);
    return rc;
}

int main(int argc, char **argv) {
    int rc = 1;
    int written = 0;
    char directory[4096] = {0};
    char gzip_path[4096] = {0};
    char mutation_path[4096] = {0};
    char output_path[4096] = {0};
    uint8_t gzip_bytes[MAX_FIXTURE_BYTES];
    uint8_t mutated[MAX_FIXTURE_BYTES];
    size_t gzip_length = 0;
    size_t case_index = 0;
    size_t mutated_length = 0;

    if (argc != 2) {
        fprintf(stderr, "usage: %s SCRATCH-PREFIX\n", argv[0]);
        return 64;
    }
    hts_set_log_level(HTS_LOG_OFF);
    written = snprintf(directory, sizeof(directory), "%s-%ld", argv[1],
                       (long)getpid());
    REQUIRE(written > 0 && (size_t)written < sizeof(directory));
    REQUIRE(mkdir(directory, 0700) == 0);
    REQUIRE(format_path(gzip_path, sizeof(gzip_path), directory,
                        "source.gzip") == 0);
    REQUIRE(format_path(mutation_path, sizeof(mutation_path), directory,
                        "mutated-input.data") == 0);
    REQUIRE(format_path(output_path, sizeof(output_path), directory,
                        "mutated-output.bam") == 0);
    REQUIRE(write_gzip(gzip_path) == 0);
    REQUIRE(load_file(gzip_path, gzip_bytes, sizeof(gzip_bytes),
                      &gzip_length) == 0);
    REQUIRE(gzip_length > 0 && gzip_length <= sizeof(mutated));

    for (case_index = 0; case_index < (size_t)READER_CASES; ++case_index) {
        memcpy(mutated, gzip_bytes, gzip_length);
        if (case_index < gzip_length) {
            mutated_length = case_index;
        } else {
            mutated_length = gzip_length;
            mutate_bytes(mutated, mutated_length,
                         UINT32_C(0xa511e9b3) ^ (uint32_t)case_index);
        }
        REQUIRE(write_bytes(mutation_path, mutated, mutated_length) == 0);
        REQUIRE(exercise_reader(mutation_path) == 0);
    }

    for (case_index = 0; case_index < (size_t)RECORD_CASES; ++case_index) {
        mutated_length = sizeof(RECORD) - 1u;
        memcpy(mutated, RECORD, mutated_length);
        if (case_index < mutated_length) {
            mutated_length = case_index;
        } else {
            mutate_bytes(mutated, mutated_length,
                         UINT32_C(0x65d2c4f1) ^ (uint32_t)case_index);
        }
        REQUIRE(exercise_record(output_path, mutated, mutated_length) == 0);
    }

    for (case_index = 0; case_index < (size_t)HEADER_CASES; ++case_index) {
        mutated_length = sizeof(HEADER) - 1u;
        memcpy(mutated, HEADER, mutated_length);
        if (case_index < mutated_length) {
            mutated_length = case_index;
        } else {
            mutate_bytes(mutated, mutated_length,
                         UINT32_C(0xc3147a25) ^ (uint32_t)case_index);
        }
        REQUIRE(exercise_header(output_path, mutated, mutated_length) == 0);
    }

    rc = 0;

done:
    (void)unlink(gzip_path);
    (void)unlink(mutation_path);
    (void)unlink(output_path);
    (void)rmdir(directory);
    if (rc == 0) {
        printf("bsbit_htslib_shim_mutation_smoke=PASS reader_cases=%u "
               "record_cases=%u header_cases=%u\n",
               READER_CASES, RECORD_CASES, HEADER_CASES);
    }
    return rc;
}
