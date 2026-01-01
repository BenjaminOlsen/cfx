#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "cfx/fac.h"
#include "cfx/big.h"
#include "cfx/algo.h"
#include "cfx_cmd.h"


static void print_binom(cfx_limb_t n, cfx_limb_t k){
    cfx_fac_t f = cfx_fac_binom(n,k);
    cfx_big_t B;
    cfx_big_init(&B);
    cfx_big_from_fac(&B, &f);
    char *s = cfx_big_dec_alloc(&B, NULL);
    size_t sz = strlen(s);
    printf("C(" CFX_PRIuLIMB ", " CFX_PRIuLIMB ") = %s\nlen: %zu\n", n, k, s, sz);
    free(s);
    cfx_big_free(&B);
    cfx_fac_free(&f);
}

static void usage(const char* prog) {
    fprintf(stderr, "Usage: %s n [k]\n", prog);
    fprintf(stderr, "  Compute binomial coefficient C(n,k)\n");
    fprintf(stderr, "  If k is omitted, prints Pascal's triangle up to row n\n");
}

int cfx_choose_run(int argc, char* argv[]) {
    if (argc < 2 || (argc == 2 && (strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0))) {
        usage(argv[0]);
        return argc < 2 ? 1 : 0;
    }

    cfx_limb_t n, k;
    if (argc == 2) {
        n = (cfx_limb_t)strtol(argv[1], NULL, 10);
        for (cfx_limb_t i = 1; i <= n; ++i) {
            for (cfx_limb_t j = 0; j <= i; ++j) {
                print_binom(i, j);
            }
            printf(".....\n");
        }
    } else if (argc == 3) {
        n = (cfx_limb_t)strtol(argv[1], NULL, 10);
        k = (cfx_limb_t)strtol(argv[2], NULL, 10);
        print_binom(n, k);
    }

    return 0;
}
