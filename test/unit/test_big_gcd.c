/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */
/* test_big_gcd.c - Tests for big integer GCD */

#include "test_common.h"

static void test_big_gcd_basic(void) {
    cfx_big_t a, b, g;
    cfx_big_init(&a);
    cfx_big_init(&b);
    cfx_big_init(&g);

    /* gcd(12, 8) = 4 */
    cfx_big_from_limb(&a, 12);
    cfx_big_from_limb(&b, 8);
    cfx_big_gcd(&g, &a, &b);
    CFX_ASSERT(g.n == 1 && g.limb[0] == 4);

    /* gcd(17, 13) = 1 (coprime) */
    cfx_big_from_limb(&a, 17);
    cfx_big_from_limb(&b, 13);
    cfx_big_gcd(&g, &a, &b);
    CFX_ASSERT(g.n == 1 && g.limb[0] == 1);

    /* gcd(0, 5) = 5 */
    cfx_big_from_limb(&a, 0);
    cfx_big_from_limb(&b, 5);
    cfx_big_gcd(&g, &a, &b);
    CFX_ASSERT(g.n == 1 && g.limb[0] == 5);

    /* gcd(100, 0) = 100 */
    cfx_big_from_limb(&a, 100);
    cfx_big_from_limb(&b, 0);
    cfx_big_gcd(&g, &a, &b);
    CFX_ASSERT(g.n == 1 && g.limb[0] == 100);

    /* gcd(48, 18) = 6 */
    cfx_big_from_limb(&a, 48);
    cfx_big_from_limb(&b, 18);
    cfx_big_gcd(&g, &a, &b);
    CFX_ASSERT(g.n == 1 && g.limb[0] == 6);

    cfx_big_free(&a);
    cfx_big_free(&b);
    cfx_big_free(&g);
}

static void test_big_gcd_large(void) {
    cfx_big_t a, b, g;
    cfx_big_init(&a);
    cfx_big_init(&b);
    cfx_big_init(&g);

    /* a = 123456789 * 1000, b = 123456789 * 777 => gcd = 123456789 */
    cfx_big_from_str(&a, "123456789000");
    cfx_big_from_str(&b, "95925925053");
    cfx_big_gcd(&g, &a, &b);

    cfx_big_t expected;
    cfx_big_init(&expected);
    cfx_big_from_limb(&expected, 123456789);
    CFX_ASSERT(cfx_big_cmp(&g, &expected) == 0);

    cfx_big_free(&a);
    cfx_big_free(&b);
    cfx_big_free(&g);
    cfx_big_free(&expected);
}

int main(void) {
    CFX_TEST(test_big_gcd_basic);
    CFX_TEST(test_big_gcd_large);
    puts("OK");
    return 0;
}
