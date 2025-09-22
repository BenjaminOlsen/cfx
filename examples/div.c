#include "cfx/big.h"
#include "cfx/macros.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

int main(int argc, char* argv[]) {
    
    char* u_in, *v_in;
    
    int hex_in = 0;  /* default to decimal */
    int hex_out = 0;

    if (argc < 3) { 
        puts("Use: ./div <optional -x for hex input> <optional -ox for hex output> <dividend> <divisor>");
    } else {
        for (int i = 0; i < argc; ++i) {
            if (strcmp(argv[i], "-x") == 0) hex_in = 1;
            else if (strcmp(argv[i], "-ox") == 0) hex_out = 1;
        }
        
        u_in = argv[argc-2];
        v_in = argv[argc-1];
    }

    cfx_big_t u, v;
    cfx_big_init(&u);
    cfx_big_init(&v);

    if (hex_in) {
        printf("hex in\n");
        cfx_big_from_hex(&u, u_in);
        cfx_big_from_hex(&v, v_in);
    } else {
        cfx_big_from_str(&u, u_in);
        cfx_big_from_str(&v, v_in);
    }

    cfx_big_t q,r;
    cfx_big_init(&q);
    cfx_big_init(&r);
    cfx_big_copy(&q, &u);

    cfx_big_divrem(&q, &r, &u, &v);

    char* u_str, *q_str, *r_str, *v_str;
    if (hex_out) {
        u_str = cfx_big_to_hex(&u, NULL);
        q_str = cfx_big_to_hex(&q, NULL);
        r_str = cfx_big_to_hex(&r, NULL);
        v_str = cfx_big_to_hex(&v, NULL);
    } else {
        u_str = cfx_big_to_str(&u, NULL);
        q_str = cfx_big_to_str(&q, NULL);
        r_str = cfx_big_to_str(&r, NULL);
        v_str = cfx_big_to_str(&v, NULL);
    }

    printf("u: %s\nd: %s\nu/v: %s\nr: %s\n",
        u_str, v_str, q_str, r_str);

    /* double check: */
    cfx_big_t p;
    cfx_big_init(&p);
    cfx_big_copy(&p, &v);
    cfx_big_mul(&p, &q); /* p = v*q */
    cfx_big_add(&p, &r); /* p = v*q + r */

    int eq = cfx_big_eq(&p, &u);
    printf("ok? %s\n", eq ? "yes!" : "NO.");
    //if (!eq) 
    {
        int cmp = cfx_big_cmp(&p, &u);
        cfx_big_t diff;
        cfx_big_init(&diff);
        if (cmp < 0) {
            cfx_big_copy(&diff, &u);
            cfx_big_sub(&diff, &p);
        } else {
            cfx_big_copy(&diff, &p);
            cfx_big_sub(&diff, &u);
        }
        char* diff_str = hex_out ? cfx_big_to_hex(&diff, NULL):cfx_big_to_str(&diff, NULL);
        printf("diff: %s\n", diff_str);
        free(diff_str);
        cfx_big_free(&diff);
    }

    char* p_str = hex_out ? cfx_big_to_hex(&p, NULL):cfx_big_to_str(&p, NULL);
    printf("p = v*q + r: %s\n", p_str);
    
    /* sanity check: */
    int scmp = strcmp(p_str, u_str);
    
    printf("strcmp: %d\n", scmp);

    PRINT_BIG("u", &u);
    printf("u, %s\n", u_str);
    PRINT_BIG("v", &v);
    printf("v, %s\n", hex_out ? cfx_big_to_hex(&v, NULL) : cfx_big_to_str(&v, NULL));
    PRINT_BIG("p", &p);
    printf("p, %s\n", p_str);
    PRINT_BIG("q", &q);
    PRINT_BIG("r", &r);

    free(u_str);
    free(q_str);
    free(v_str);
    free(p_str);
    
    cfx_big_free(&u);
    cfx_big_free(&v);
    cfx_big_free(&q);
    cfx_big_free(&r);
    cfx_big_free(&p);
    
    return 0;
}
