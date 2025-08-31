#include "cfx/primes.h"
#include <stdlib.h>

/* find all primes up to n using sieve of eratosthenes */
cfx_vec_u32 cfx_sieve_primes(uint32_t n) {
    cfx_vec_u32 primes = {0};
    if (n < 2) return primes;

    /* mark all composites, multiples of i until sqrt(n) */
    /* mark[x] == 0: x prime, mark[x] == 1: x composite */
    uint8_t* mark = (uint8_t*)calloc(n + 1, 1);

    for (uint32_t i = 2; (i * i) <= n; ++i) {
        if (!mark[i]) {
            for (uint32_t j = i*i; j <= n; j += i) {
                // printf("marking %u composite\n", j);
                mark[j] = 1;
            }
        }
    }
    
    size_t p_cnt = 0;
    for (uint32_t i = 2; i < n + 1; ++i) {
        if (!mark[i]) ++p_cnt;
    }

    primes.data = (uint32_t*)malloc(p_cnt * sizeof(uint32_t));
    primes.size = p_cnt;

    size_t k = 0;
    for (uint32_t i = 2; i <= n; ++i) {
        if (!mark[i]) {
            primes.data[k] = i;
            ++k;
        }
    }

    free(mark);
    return primes;
}

/* uses legendre's formula to give th power of p that divides n! */
uint32_t cfx_legendre(uint32_t n, uint32_t p) {
    uint32_t e = 0;
    while (n) {
        n /= p;
        e += n;
    }
    return e;
}