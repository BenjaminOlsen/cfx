#include "cfx/algo.h"

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "cfx_utils.h"

static int dec_digits(cfx_limb_t n) {
    if (n < 10ULL) return 1;
    if (n < 100ULL) return 2;
    if (n < 1000ULL) return 3;
    if (n < 10000ULL) return 4;
    if (n < 100000ULL) return 5;
    if (n < 1000000ULL) return 6;
    if (n < 10000000ULL) return 7;
    if (n < 100000000ULL) return 8;
    if (n < 1000000000ULL) return 9;
    if (n < 10000000000ULL) return 10;
    if (n < 100000000000ULL) return 11;
    if (n < 1000000000000ULL) return 12;
    if (n < 10000000000000ULL) return 13;
    if (n < 100000000000000ULL) return 14;
    if (n < 1000000000000000ULL) return 15;
    if (n < 10000000000000000ULL) return 16;
    if (n < 100000000000000000ULL) return 17;
    if (n < 1000000000000000000ULL) return 18;
    return 19;
}

static void usage(const char* prog) {
    fprintf(stderr, "Usage: %s <n> [col_width]\n", prog);
    fprintf(stderr, "  List all primes up to n using sieve of Eratosthenes\n");
    fprintf(stderr, "  col_width  optional line width for formatting output\n");
}

int cfx_primes_run(int argc, char* argv[]) {
    if (argc < 2 || (argc == 2 && (strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0))) {
        usage(argv[0]);
        return argc < 2 ? 1 : 0;
    }

    int dmax = INT32_MAX;
    if (argc == 3) dmax = (size_t)strtol(argv[2], NULL, 10);

    cfx_limb_t n = (cfx_limb_t)strtol(argv[1], NULL, 10);
    cfx_vec_t primes = cfx_sieve_primes(n);
    int dcnt = 0;
    for (size_t k = 0; k < primes.size-1; ++k) {
        int c = dec_digits(primes.data[k]) + 2; /* +2 for comma + space */
        if (dcnt + c > dmax) {
            printf("\n");
            dcnt = c;
        } else {
            dcnt += c;
        }
        printf("" CFX_PRIuLIMB ", ", primes.data[k]);
    }
    printf("" CFX_PRIuLIMB "\n", primes.data[primes.size-1]);

    printf("found %zu primes until " CFX_PRIuLIMB " (0x" CFX_PRIuLIMB ")\n", primes.size, n, n);

    cfx_vec_free(&primes);
    return 0;
}
