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
    
    cfx_limb_t n = (cfx_limb_t)strtol(argv[1], NULL, 10);
    printf(""CFX_PRIuLIMB"\n", n);
    
    cfx_vec_t primes = cfx_sieve_primes(n);
    printf("found %zu primes to "CFX_PRIuLIMB"\n", primes.size, n);
    
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
    
    printf(""CFX_PRIuLIMB"! = %s\ndigits: %zu\nlimbs: %zu\n", n, s, sz, b.n);
    double lg = cfx_big_log(&b, 10);
    printf("log10("CFX_PRIuLIMB"!) = %.4f\n", n, lg);
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
