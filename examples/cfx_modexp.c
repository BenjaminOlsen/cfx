#include <stdio.h>

#include "cfx/big.h"
#include "cfx/macros.h"

void usage(const char* prog) {
    fprintf(stderr, "Usage: %s base exp mod \n", prog);
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        usage(argv[0]);
        return 1;
    }

    cfx_big_t r, base, exp, mod;

    cfx_big_init(&base);
    cfx_big_init(&exp);
    cfx_big_init(&mod);

    cfx_big_from_dec(&base, argv[1]);
    cfx_big_from_dec(&exp, argv[2]);
    cfx_big_from_dec(&mod, argv[3]);
    cfx_big_mod_exp(&r, &base, &exp, &mod);

    PRINT_BIG("ayy", &r);
    char* s = cfx_big_to_str(&r, NULL);
    printf("%s\n", s);
    cfx_big_free(&mod);
    cfx_big_free(&base);
    cfx_big_free(&exp);
    free(s);
    return 0;
}