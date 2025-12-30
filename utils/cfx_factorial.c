#include "cfx/fac.h"
#include "cfx/big.h"
#include "cfx/algo.h"

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "cfx_utils.h"

static void usage(const char* prog) {
    fprintf(stderr, "Usage: %s <n> [-x]\n", prog);
    fprintf(stderr, "  Compute n! (factorial)\n");
    fprintf(stderr, "  -x  print result in hex\n");
}

int cfx_factorial_run(int argc, char* argv[]) {
    if (argc < 2 || (argc == 2 && (strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0))) {
        usage(argv[0]);
        return argc < 2 ? 1 : 0;
    }

    int cfx_print_hex = 0;

    if ((argc == 3) && (strcmp(argv[2], "-x")==0)) {
        cfx_print_hex = 1;
    }

    cfx_limb_t n = (cfx_limb_t)strtol(argv[1], NULL, 10);
    printf("" CFX_PRIuLIMB "\n", n);

    cfx_vec_t primes = cfx_sieve_primes(n);
    printf("found %zu primes to " CFX_PRIuLIMB "\n", primes.size, n);

    cfx_fac_t fac;
    cfx_fac_init(&fac);
    printf("calculating factorial... \n");
    cfx_fac_factorial(&fac, n, &primes);
    printf("done, fac size: %zu\n", fac.len);

    cfx_big_t b;
    cfx_big_init(&b);
    printf("making big int... ");
    fflush(stdout);
    cfx_big_from_fac_faster(&b, &fac);
    printf("done, limbs: %zu\n", b.n);
    size_t sz = 0;
    char* s;

    if (cfx_print_hex) {
        s = cfx_big_to_hex(&b, &sz);
    } else {
        s = cfx_big_to_str(&b, &sz);
    }

    printf("" CFX_PRIuLIMB "! = %s\ndigits: %zu\nlimbs: %zu\n", n, s, sz, b.n);
    double lg = cfx_big_log(&b, 10);
    printf("log10(" CFX_PRIuLIMB "!) = %.4f\n", n, lg);
    char sci[128];
    cfx_big_to_sci(&b, 10, 5, sci, sizeof(sci));
    printf("%s\n", sci);
    char sci2[128];
    cfx_big_to_sci(&b, 100, 3, sci2, sizeof(sci2));
    printf("%s\n", sci2);

    char sci3[128];
    cfx_big_to_sci(&b, 3, 4, sci3, sizeof(sci2));
    printf("%s\n", sci3);
    free(s);
    cfx_big_free(&b);
    cfx_fac_free(&fac);
    cfx_vec_free(&primes);

    return 0;
}
