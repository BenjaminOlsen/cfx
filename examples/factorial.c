#include "cfx/fac.h"
#include "cfx/big.h"
#include "cfx/algo.h"

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>


int main(int argc, char* argv[]) {
    if (argc < 2) {
        puts("use: factorial <number> <optional: '-x' to print hex>\n");
        return 1;
    }
    
    int print_hex = 0;
    
    if ((argc == 3) && (strcmp(argv[2], "-x")==0)) {
        print_hex = 1;
    }
    
    uint64_t n = (uint64_t)strtol(argv[1], NULL, 10);
    printf("%llu\n", n);
    
    cfx_vec_t primes = cfx_sieve_primes(n);
    printf("found %zu primes to %llu\n", primes.size, n);
    
    cfx_fac_t fac;
    cfx_fac_init(&fac);
    printf("calculating factorial... \n");
    cfx_fac_factorial(&fac, n, &primes);
    printf("done, fac size: %zu\n", fac.len);
    
    cfx_big_t b;
    cfx_big_init(&b);
    printf("making big int... ");
    fflush(stdout);
    cfx_big_from_fac(&b, &fac);
    printf("done, limbs: %zu\n", b.n);
    size_t sz = 0;
    char* s;
    
    if (print_hex) {
        s = cfx_big_to_hex(&b, &sz);
    } else {
        s = cfx_big_to_str(&b, &sz);
    }
    
    printf("%llu! = %s\ndigits: %zu\nlimbs: %zu\n", n, s, sz, b.n);
    free(s);
    cfx_big_free(&b);
    cfx_fac_free(&fac);
    cfx_vec_free(&primes);

    return 0;
}
