#include "cfx/big.h"
#include "cfx/macros.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "cfx_cmd.h"

static void usage(const char* prog) {
    fprintf(stderr, "Usage: %s [-x] [-ox] <dividend> <divisor>\n", prog);
    fprintf(stderr, "  -x   interpret input as hex\n");
    fprintf(stderr, "  -ox  output in hex\n");
}

int cfx_div_run(int argc, char* argv[]) {

    char* u_in, *v_in;

    int hex_in = 0;  /* default to decimal */
    int hex_out = 0;

    if (argc >= 2 && (strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0)) {
        usage(argv[0]);
        return 0;
    }

    if (argc < 3) {
        usage(argv[0]);
        return 1;
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
        cfx_big_from_dec(&u, u_in);
        cfx_big_from_dec(&v, v_in);
    }

    cfx_big_t q,r;
    cfx_big_init(&q);
    cfx_big_init(&r);
    cfx_big_copy(&q, &u);

    cfx_big_divrem(&q, &r, &u, &v);

    char* u_str, *q_str, *r_str, *v_str;
    if (hex_out) {
        u_str = cfx_big_hex_alloc(&u, NULL);
        q_str = cfx_big_hex_alloc(&q, NULL);
        r_str = cfx_big_hex_alloc(&r, NULL);
        v_str = cfx_big_hex_alloc(&v, NULL);
    } else {
        u_str = cfx_big_dec_alloc(&u, NULL);
        q_str = cfx_big_dec_alloc(&q, NULL);
        r_str = cfx_big_dec_alloc(&r, NULL);
        v_str = cfx_big_dec_alloc(&v, NULL);
    }

    printf("u: %s\nd: %s\nu/v: %s\nr: %s\n",
        u_str, v_str, q_str, r_str);

    /* double check: */
    cfx_big_t p;
    cfx_big_init(&p);
    cfx_big_mul(&p, &v, &q); /* p = v*q */
    cfx_big_add_eq(&p, &r); /* p = v*q + r */

    int eq = cfx_big_eq(&p, &u);
    printf("ok? %s\n", eq ? "yes!" : "NO.");
    /*if (!eq)  */
    {
        int cmp = cfx_big_cmp(&p, &u);
        cfx_big_t diff;
        cfx_big_init(&diff);
        if (cmp < 0) {
            cfx_big_copy(&diff, &u);
            cfx_big_sub_eq(&diff, &p);
        } else {
            cfx_big_copy(&diff, &p);
            cfx_big_sub_eq(&diff, &u);
        }
        char* diff_str = hex_out ? cfx_big_hex_alloc(&diff, NULL):cfx_big_dec_alloc(&diff, NULL);
        printf("diff: %s\n", diff_str);
        free(diff_str);
        cfx_big_free(&diff);
    }

    char* p_str = hex_out ? cfx_big_hex_alloc(&p, NULL):cfx_big_dec_alloc(&p, NULL);
    printf("p = v*q + r: %s\n", p_str);

    /* sanity check: */
    int scmp = strcmp(p_str, u_str);

    printf("strcmp: %d\n", scmp);

    PRINT_BIG("u", &u);
    printf("u, %s\n", u_str);
    PRINT_BIG("v", &v);
    printf("v, %s\n", hex_out ? cfx_big_hex_alloc(&v, NULL) : cfx_big_dec_alloc(&v, NULL));
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
