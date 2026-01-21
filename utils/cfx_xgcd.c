/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#include "cfx/big.h"
#include "cfx/sbig.h"
#include "cfx/macros.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "cfx_cmd.h"

static void usage(const char* prog) {
    fprintf(stderr, "Usage: %s [-x] [-ox] <a> <b>\n", prog);
    fprintf(stderr, "Compute extended GCD: finds g, x, y such that a*x + b*y = g\n");
    fprintf(stderr, "  -x   interpret input as hex\n");
    fprintf(stderr, "  -ox  output in hex\n");
}

int cfx_xgcd_run(int argc, char* argv[]) {
    char* a_in;
    char* b_in;

    int hex_in = 0;
    int hex_out = 0;

    if (argc >= 2 && (strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0)) {
        usage(argv[0]);
        return 0;
    }

    if (argc < 3) {
        usage(argv[0]);
        return 1;
    }

    for (int i = 0; i < argc; ++i) {
        if (strcmp(argv[i], "-x") == 0) hex_in = 1;
        else if (strcmp(argv[i], "-ox") == 0) hex_out = 1;
    }

    a_in = argv[argc - 2];
    b_in = argv[argc - 1];

    cfx_big_t a, b, g;
    cfx_sbig_t x, y;
    cfx_big_init(&a);
    cfx_big_init(&b);
    cfx_big_init(&g);
    cfx_sbig_init(&x);
    cfx_sbig_init(&y);

    if (hex_in) {
        cfx_big_from_hex(&a, a_in);
        cfx_big_from_hex(&b, b_in);
    } else {
        cfx_big_from_dec(&a, a_in);
        cfx_big_from_dec(&b, b_in);
    }

    cfx_big_xgcd(&g, &x, &y, &a, &b);

    char* a_str = hex_out ? cfx_big_hex_alloc(&a, NULL) : cfx_big_dec_alloc(&a, NULL);
    char* b_str = hex_out ? cfx_big_hex_alloc(&b, NULL) : cfx_big_dec_alloc(&b, NULL);
    char* g_str = hex_out ? cfx_big_hex_alloc(&g, NULL) : cfx_big_dec_alloc(&g, NULL);
    char* x_str = cfx_sbig_dec_alloc(&x, NULL);
    char* y_str = cfx_sbig_dec_alloc(&y, NULL);

    printf("a: %s\n", a_str);
    printf("b: %s\n", b_str);
    printf("gcd(a, b): %s\n", g_str);
    printf("x: %s\n", x_str);
    printf("y: %s\n", y_str);
    printf("a*x + b*y = gcd\n");

    /* Verify: compute a*x + b*y and check it equals g */
    cfx_sbig_t sa, sb, ax, by, sum;
    cfx_sbig_init(&sa);
    cfx_sbig_init(&sb);
    cfx_sbig_init(&ax);
    cfx_sbig_init(&by);
    cfx_sbig_init(&sum);

    cfx_sbig_assign_big(&sa, &a, 1);
    cfx_sbig_assign_big(&sb, &b, 1);
    cfx_sbig_mul(&ax, &sa, &x);
    cfx_sbig_mul(&by, &sb, &y);
    cfx_sbig_add(&sum, &ax, &by);

    char* sum_str = cfx_sbig_dec_alloc(&sum, NULL);
    int ok = (strcmp(sum_str, g_str) == 0);
    printf("verification (a*x + b*y): %s - %s\n", sum_str, ok ? "OK" : "FAIL");

    free(a_str);
    free(b_str);
    free(g_str);
    free(x_str);
    free(y_str);
    free(sum_str);

    cfx_big_free(&a);
    cfx_big_free(&b);
    cfx_big_free(&g);
    cfx_sbig_free(&x);
    cfx_sbig_free(&y);
    cfx_sbig_free(&sa);
    cfx_sbig_free(&sb);
    cfx_sbig_free(&ax);
    cfx_sbig_free(&by);
    cfx_sbig_free(&sum);

    return ok ? 0 : 1;
}
