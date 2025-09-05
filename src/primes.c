#include "cfx/primes.h"
#include "cfx/macros.h"

#include <stdlib.h>

/* find all primes up to n using sieve of eratosthenes */
cfx_vec_t cfx_sieve_primes(uint64_t n) {
    cfx_vec_t primes = {0};
    if (n < 2) return primes;

    /* mark all composites, multiples of i until sqrt(n) */
    /* mark[x] == 0: x prime, mark[x] == 1: x composite */
    uint8_t* mark = (uint8_t*)calloc(n + 1, 1);

    for (uint64_t i = 2; (i * i) <= n; ++i) {
        if (!mark[i]) {
            for (uint64_t j = i*i; j <= n; j += i) {
                // CFX_PRINT_DBG("marking %u composite\n", j);
                mark[j] = 1;
            }
        }
    }
    
    size_t p_cnt = 0;
    for (uint64_t i = 2; i < n + 1; ++i) {
        if (!mark[i]) ++p_cnt;
    }

    primes.data = (uint64_t*)malloc(p_cnt * sizeof(uint64_t));
    primes.size = p_cnt;

    size_t k = 0;
    for (uint64_t i = 2; i <= n; ++i) {
        if (!mark[i]) {
            primes.data[k] = i;
            ++k;
        }
    }

    free(mark);
    return primes;
}

/* uses legendre's formula to give th power of p that divides n! */
uint64_t cfx_legendre(uint64_t n, uint64_t p) {
    uint64_t e = 0;
    while (n) {
        n /= p;
        e += n;
    }
    return e;
}
