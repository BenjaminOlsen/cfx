#include "cfx/algo.h"

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

int dec_digits(uint64_t n) {
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

int main(int argc, char* argv[]) {
    if (argc < 2) {
        printf("use: primes <n> <col width (optional)>:\n"
            " -- find primes up to n, printed within col width columns for formatting\n"
            "e.g. 'primes 800 100': find primes to 800 printed with line width of 100");
        return 1;
    }
    
    int dmax = INT32_MAX;
    if (argc == 3) dmax = (size_t)strtol(argv[2], NULL, 10);

    uint64_t n = (uint64_t)strtol(argv[1], NULL, 10);
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
        printf("%llu, ", primes.data[k]);
    }
    printf("%llu\n", primes.data[primes.size-1]);

    printf("found %zu primes until %llu (0x%llx)\n", primes.size, n, n);

    cfx_vec_free(&primes);
    return 0;
}
