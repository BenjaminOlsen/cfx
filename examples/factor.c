#include "cfx/big.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

static void usage(const char* prog) {
    fprintf(stderr, "Usage: %s [-x|-d] [-ix|-id] n\n", prog);
    fprintf(stderr, "    calculate prime factorization of n\n");
    fprintf(stderr, "   -ix : interpret input as hex\n");
    fprintf(stderr, "   -id : interpret input as decimal\n");
}


int main(int argc, char* argv[]) {

    if (argc < 2) {
        usage(argv[0]);
        return 1;
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

    size_t sz = 0;
    char* s =  cfx_big_to_hex(&n, &sz);
    printf("%s\n", s);

    cfx_fac_t fac;
    cfx_fac_init(&fac);
    cfx_big_to_fac(&fac, &n);
    cfx_fac_print(&fac);
    cfx_fac_free(&fac);
    cfx_big_free(&n);
    return 0;
}
