/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */
/* test_cli_common.h - stdout capture and assertion helpers for CLI tests */

#ifndef CFX_TEST_CLI_COMMON_H
#define CFX_TEST_CLI_COMMON_H

#include "cfx/macros.h"
#include "cfx_cmd.h"

#include <stdio.h>
#include <string.h>
#include <unistd.h>

/* ---- stdout capture */

static int g_saved_fd;
static FILE *g_capture_file;

static inline void capture_begin(void) {
    fflush(stdout);
    g_saved_fd = dup(STDOUT_FILENO);
    g_capture_file = tmpfile();
    dup2(fileno(g_capture_file), STDOUT_FILENO);
}

static inline size_t capture_end(char *buf, size_t bufsz) {
    fflush(stdout);
    dup2(g_saved_fd, STDOUT_FILENO);
    close(g_saved_fd);
    rewind(g_capture_file);
    size_t n = fread(buf, 1, bufsz - 1, g_capture_file);
    buf[n] = '\0';
    fclose(g_capture_file);
    return n;
}

/* run a cfx_<name>_run function, capture its stdout into buf.
 * returns the function's exit code. */
static inline int run_capture(int (*fn)(int, char **),
                               int argc, char **argv,
                               char *buf, size_t bufsz) {
    capture_begin();
    int rc = fn(argc, argv);
    capture_end(buf, bufsz);
    return rc;
}

static inline void assert_output_eq(const char *got, const char *expected) {
    if (strcmp(got, expected) != 0) {
        fprintf(stderr, "  expected: \"%s\"\n  got:      \"%s\"\n", expected, got);
        CFX_ASSERT(0);
    }
}

static inline void assert_output_contains(const char *got, const char *needle) {
    if (strstr(got, needle) == NULL) {
        fprintf(stderr, "  output missing \"%s\"\n  got: \"%s\"\n", needle, got);
        CFX_ASSERT(0);
    }
}

static inline void assert_output_starts_with(const char *got, const char *prefix) {
    if (strncmp(got, prefix, strlen(prefix)) != 0) {
        fprintf(stderr, "  expected prefix: \"%s\"\n  got: \"%s\"\n", prefix, got);
        CFX_ASSERT(0);
    }
}

#endif /* CFX_TEST_CLI_COMMON_H */
