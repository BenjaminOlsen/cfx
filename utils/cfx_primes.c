#include "cfx/algo.h"

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "cfx_cmd.h"

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
    #if CFX_LIMB_BITS == 64
    if (n < 10000000000ULL) return 10;
    if (n < 100000000000ULL) return 11;
    if (n < 1000000000000ULL) return 12;
    if (n < 10000000000000ULL) return 13;
    if (n < 100000000000000ULL) return 14;
    if (n < 1000000000000000ULL) return 15;
    if (n < 10000000000000000ULL) return 16;
    if (n < 100000000000000000ULL) return 17;
    if (n < 1000000000000000000ULL) return 18;
    #endif
    return 19;
}

static void usage(const char* prog) {
    fprintf(stderr, "Usage: %s [options] <n>\n", prog);
    fprintf(stderr, "  List all primes up to n using sieve of Eratosthenes\n\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  --width=N    Line width for formatting output (default: no wrap)\n");
    fprintf(stderr, "  --sep=STR    Separator string (default: ', ')\n");
    fprintf(stderr, "  -h, --help   Show this help\n");
}

int cfx_primes_run(int argc, char* argv[]) {
    int dmax = INT32_MAX;
    const char* sep = ", ";
    cfx_limb_t n = 0;
    int have_n = 0;

    for (int i = 1; i < argc; i++) {
        const char* arg = argv[i];
        if (strcmp(arg, "--help") == 0 || strcmp(arg, "-h") == 0) {
            usage(argv[0]);
            return 0;
        } else if (strncmp(arg, "--width=", 8) == 0) {
            dmax = (int)strtol(arg + 8, NULL, 10);
        } else if (strncmp(arg, "--sep=", 6) == 0) {
            sep = arg + 6;
        } else {
            n = (cfx_limb_t)strtol(arg, NULL, 10);
            have_n = 1;
        }
    }

    if (!have_n) {
        usage(argv[0]);
        return 1;
    }

    size_t sep_len = strlen(sep);
    cfx_vec_t primes = cfx_sieve_primes(n);
    int dcnt = 0;
    for (size_t k = 0; k < primes.size-1; ++k) {
        int c = dec_digits(primes.data[k]) + (int)sep_len;
        if (dcnt + c > dmax) {
            printf("\n");
            dcnt = c;
        } else {
            dcnt += c;
        }
        printf(CFX_PRIuLIMB "%s", primes.data[k], sep);
    }
    printf(CFX_PRIuLIMB "\n", primes.data[primes.size-1]);

    printf("found %zu primes until " CFX_PRIuLIMB " (0x" CFX_PRIuLIMB ")\n", primes.size, n, n);

    cfx_vec_free(&primes);
    return 0;
}
