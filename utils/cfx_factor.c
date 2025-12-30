#include "cfx/big.h"
#include "cfx/fac.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "cfx_utils.h"

static void usage(const char* prog) {
    fprintf(stderr, "Usage: %s [-x|-d] [-ix|-id] n\n", prog);
    fprintf(stderr, "    calculate prime factorization of n\n");
    fprintf(stderr, "   -ix : interpret input as hex\n");
    fprintf(stderr, "   -id : interpret input as decimal\n");
}


int cfx_factor_run(int argc, char* argv[]) {

    if (argc < 2 || (argc == 2 && (strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0))) {
        usage(argv[0]);
        return argc < 2 ? 1 : 0;
    }

    /* Parse flags */
    int argi = 1;
    int base = 10; /* default */

    if (strcmp(argv[argi], "-ix") == 0) {
        base = 16;
        argi++;
    } else if (strcmp(argv[argi], "-id") == 0) {
        base = 10;
        argi++;
    }

    if (argi > argc - 1) {
        usage(argv[0]);
        return 1;
    }

    char* nstr = argv[argi];
    cfx_big_t n;
    cfx_big_init(&n);

    if (base == 16) {
        cfx_big_from_hex(&n, nstr);
    } else if (base == 10) {
        cfx_big_from_dec(&n, nstr);
    }

    cfx_fac_t fac;
    cfx_big_t remainder;
    cfx_big_init(&remainder);

    int rc = cfx_big_to_fac(&fac, &n, &remainder);
    cfx_fac_print(&fac);

    if (rc == 1) {
        char* rem_str = cfx_big_to_str(&remainder, NULL);
        fprintf(stderr, "(incomplete: unfactored composite = %s)\n", rem_str);
        free(rem_str);
    }

    cfx_fac_free(&fac);
    cfx_big_free(&remainder);
    cfx_big_free(&n);
    return 0;
}
