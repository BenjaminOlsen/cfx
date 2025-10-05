#include "cfx/big.h"


#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>


static void usage(const char* prog) {
    fprintf(stderr, "Usage: %s [-x|-d] n\n", prog);
    fprintf(stderr, "    miller-rabin primality test of n\n");
    fprintf(stderr, "   -x : interpret input as hex\n");
    fprintf(stderr, "   -d : interpret input as decimal\n");
}


int main(int argc, char* argv[]) {
    if (argc < 2) {
        usage(argv[0]);
        return 1;
    }

     /* Parse flags */
    int argi = 1;
    int base = 10; /* default */

    if (strcmp(argv[argi], "-x") == 0) {
        printf("reading hex\n");
        base = 16;
        ++argi;
    } else if (strcmp(argv[argi], "-d") == 0) {
        printf("reading dec\n");
        base = 10;
        ++argi;
    }

    if (argi > argc - 1) {
        printf("stopping: %d >= %d-1 hex\n", argi, argc);
        usage(argv[0]);
        return 1;
    }
    char* nstr = argv[argi];

    cfx_big_t n;
    cfx_big_init(&n);
    if (base == 16) {
        cfx_big_from_hex(&n, nstr);
    } else if (base == 10) {
        cfx_big_from_str(&n, nstr);
    }

    int isprime = cfx_big_is_prime(&n);
    printf("prime? %s\n", isprime ? "YES!" : "no");
    return 0;
}
