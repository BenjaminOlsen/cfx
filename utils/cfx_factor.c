#include "cfx/big.h"
#include "cfx/fac.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "cfx_cmd.h"
#include "common.h"

static void usage(const char* prog) {
    fprintf(stderr, "Usage: %s [options] [n]\n", prog);
    fprintf(stderr, "  Calculate prime factorization of n.\n");
    fprintf(stderr, "  If n is omitted, reads from stdin.\n\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -x, -ix   Interpret input as hex\n");
    fprintf(stderr, "  -d, -id   Interpret input as decimal (default)\n");
    fprintf(stderr, "  -h, --help  Show this help\n\n");
    fprintf(stderr, "Examples:\n");
    fprintf(stderr, "  %s 12345678901234567\n", prog);
    fprintf(stderr, "  %s -x 0x123abc\n", prog);
    fprintf(stderr, "  cfx primegen 64 | %s\n", prog);
}

int cfx_factor_run(int argc, char* argv[]) {
    if (argc == 2 && (strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0)) {
        usage(argv[0]);
        return 0;
    }

    int argi = 1;
    int base = 0; /* auto-detect */

    while (argi < argc && argv[argi][0] == '-') {
        if (strcmp(argv[argi], "-ix") == 0 || strcmp(argv[argi], "-x") == 0) {
            base = 16;
            argi++;
        } else if (strcmp(argv[argi], "-id") == 0 || strcmp(argv[argi], "-d") == 0) {
            base = 10;
            argi++;
        } else if (strcmp(argv[argi], "-h") == 0 || strcmp(argv[argi], "--help") == 0) {
            usage(argv[0]);
            return 0;
        } else {
            fprintf(stderr, "Unknown option: %s\n", argv[argi]);
            usage(argv[0]);
            return 1;
        }
    }

    char* nstr = NULL;
    int free_nstr = 0;

    if (argi < argc) {
        nstr = argv[argi];
    } else {
        nstr = cfx_read_line_stdin();
        if (!nstr || *nstr == '\0') {
            fprintf(stderr, "Error: no input\n");
            free(nstr);
            return 1;
        }
        free_nstr = 1;
    }

    /* auto-detect hex if starts with 0x */
    if (base == 0) {
        if (strncmp(nstr, "0x", 2) == 0 || strncmp(nstr, "0X", 2) == 0) {
            base = 16;
        } else {
            base = 10;
        }
    }

    cfx_big_t n;
    cfx_big_init(&n);

    /* skip 0x prefix for hex */
    const char* numstr = nstr;
    if (base == 16 && (strncmp(numstr, "0x", 2) == 0 || strncmp(numstr, "0X", 2) == 0)) {
        numstr += 2;
    }

    if (base == 16) {
        cfx_big_from_hex(&n, numstr);
    } else {
        cfx_big_from_dec(&n, numstr);
    }

    if (free_nstr) free(nstr);

    cfx_fac_t fac;
    cfx_big_t remainder;
    cfx_big_init(&remainder);

    int rc = cfx_big_to_fac(&fac, &n, &remainder);
    cfx_fac_print(&fac);

    if (rc == 1) {
        char* rem_str = cfx_big_dec_alloc(&remainder, NULL);
        fprintf(stderr, "(incomplete: unfactored composite = %s)\n", rem_str);
        free(rem_str);
    }

    cfx_fac_free(&fac);
    cfx_big_free(&remainder);
    cfx_big_free(&n);
    return 0;
}
