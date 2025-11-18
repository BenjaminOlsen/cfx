#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "cfx/fac.h"
#include "cfx/big.h"
#include "cfx/algo.h"


static void print_binom(cfx_limb_t n, cfx_limb_t k){
    cfx_fac_t f = cfx_fac_binom(n,k);
    cfx_big_t B;
    cfx_big_init(&B);
    cfx_big_from_fac(&B, &f);
    char *s = cfx_big_to_str(&B, NULL);
    size_t sz = strlen(s);
    printf("C("CFX_PRIuLIMB", "CFX_PRIuLIMB") = %s\nlen: %zu\n", n, k, s, sz);
    free(s);
    cfx_big_free(&B);
    cfx_fac_free(&f);
}

int main(int argc, char* argv[]) {
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
