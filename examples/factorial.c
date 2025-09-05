#include "cfx/fac.h"
#include "cfx/big.h"
#include "cfx/primes.h"

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char* argv[]) {
    if (argc < 2) return 1;
    uint64_t n = (uint64_t)strtol(argv[1], NULL, 10);
    printf("%llu\n", n);
    cfx_vec_t primes = cfx_sieve_primes(n);
    printf("found %zu primes to %llu\n", primes.size, n);
    cfx_fac_t fac;
    cfx_fac_init(&fac);
    cfx_fac_factorial(&fac, n, &primes);

    cfx_big_t b;
    cfx_big_init(&b);
    cfx_big_from_fac(&b, &fac);
    size_t sz = 0;
    char* s = cfx_big_to_str(&b, &sz);
    printf("%llu! = %s\nlen: %zu\n", n, s, sz);
    free(s);
    cfx_big_free(&b);
    cfx_fac_clear(&fac);

    return 0;
}
