/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */
/* test_big_init.c - Tests for big integer initialization and memory operations */

#include "test_common.h"

static void test_cfx_big_init(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    CFX_ASSERT(b.limb == NULL);
    CFX_ASSERT(b.n == 0);
    CFX_ASSERT(b.cap == 0);
    PRINT_TEST(1);
}

static void test_cfx_big_reserve(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    size_t rcap1 = 55;
    cfx_big_reserve(&b, rcap1);
    CFX_ASSERT(b.cap >= rcap1);
    CFX_ASSERT(b.n == 0);
    CFX_ASSERT(b.limb != NULL);
    size_t rcap2 = rcap1 / 2;
    cfx_big_reserve(&b, rcap2);  /* shouldn't reserve less space. */
    CFX_ASSERT(b.cap >= rcap1);
    cfx_big_free(&b);
    PRINT_TEST(1);
}

static void test_cfx_big_assign(void) {
    cfx_big_t a, b;
    cfx_big_init(&a);
    cfx_big_init(&b);
    cfx_limb_t al[] = {0x12371237, 0x0, 0x2, CFX_LIMB_MAX-1};
    cfx_big_from_limbs(&a, al, sizeof(al)/sizeof(al[0]));
    cfx_big_assign(&b, &a);
    int c = cfx_big_cmp(&a, &b);
    CFX_ASSERT(c == 0);
    cfx_big_free(&a);
    cfx_big_free(&b);
    PRINT_TEST(1);
}

static void test_copy_swap(void) {
    cfx_big_t a, b;
    cfx_big_init(&a);
    cfx_big_init(&b);
    cfx_limb_t la[] = {1, 2, 3, 4};
    cfx_limb_t lb[] = {CFX_LIMB_MAX-1, CFX_LIMB_MAX-2, CFX_LIMB_MAX-3, CFX_LIMB_MAX-4};

    cfx_big_from_limbs(&a, la, 4);
    cfx_big_from_limbs(&b, lb, 4);
    CFX_ASSERT(a.n == b.n);

    for (size_t i = 0; i < sizeof(la)/sizeof(cfx_limb_t); ++i) {
        CFX_ASSERT(a.limb[i] == la[i]);
        CFX_ASSERT(b.limb[i] == lb[i]);
    }
    cfx_big_t aa, bb;
    cfx_big_init(&aa);
    cfx_big_init(&bb);
    cfx_big_copy(&aa, &a);
    cfx_big_copy(&bb, &b);

    CFX_ASSERT(cfx_big_eq(&aa, &a));
    CFX_ASSERT(cfx_big_eq(&bb, &b));

    cfx_big_swap(&a, &b);

    CFX_ASSERT(cfx_big_eq(&aa, &b));
    CFX_ASSERT(cfx_big_eq(&bb, &a));

    cfx_big_free(&a);
    cfx_big_free(&b);
    cfx_big_free(&aa);
    cfx_big_free(&bb);
    PRINT_TEST(1);
}

static void test_copy_self(void) {
    cfx_big_t a;
    cfx_big_init(&a);
    cfx_big_from_limb(&a, 42);
    cfx_big_copy(&a, &a);
    CFX_ASSERT(cfx_big_eq_u64(&a, 42));
    cfx_big_free(&a);
}

static void test_assign_self(void) {
    cfx_big_t a;
    cfx_big_init(&a);
    cfx_big_from_limb(&a, 99);
    cfx_big_assign(&a, &a);
    CFX_ASSERT(cfx_big_eq_u64(&a, 99));
    cfx_big_free(&a);
}

static void test_swap_self(void) {
    cfx_big_t a;
    cfx_big_init(&a);
    cfx_big_from_limb(&a, 123);
    cfx_big_swap(&a, &a);
    CFX_ASSERT(cfx_big_eq_u64(&a, 123));
    cfx_big_free(&a);
}

static void test_move_self(void) {
    cfx_big_t a;
    cfx_big_init(&a);
    cfx_big_from_limb(&a, 456);
    cfx_big_move(&a, &a);
    CFX_ASSERT(cfx_big_eq_u64(&a, 456));
    cfx_big_free(&a);
}

static void test_assign_sm_zero(void) {
    cfx_big_t a;
    cfx_big_init(&a);
    cfx_big_from_limb(&a, 100);
    cfx_big_assign_sm(&a, 0);
    CFX_ASSERT(cfx_big_is_zero(&a));
    cfx_big_free(&a);
}

static void test_copy_zero(void) {
    cfx_big_t a, b;
    cfx_big_init(&a);
    cfx_big_init(&b);
    cfx_big_from_limb(&b, 100);
    cfx_big_copy(&b, &a);
    CFX_ASSERT(cfx_big_is_zero(&b));
    cfx_big_free(&b);
    cfx_big_free(&a);
}

static void test_is_zero_variations(void) {
    cfx_big_t a;
    cfx_big_init(&a);
    CFX_ASSERT(cfx_big_is_zero(&a));  /* n == 0 */
    cfx_big_reserve(&a, 1);
    a.n = 1;
    a.limb[0] = 0;
    CFX_ASSERT(cfx_big_is_zero(&a));  /* n == 1 && limb[0] == 0 */
    a.limb[0] = 1;
    CFX_ASSERT(!cfx_big_is_zero(&a));
    cfx_big_free(&a);
}

static void test_copy_needs_realloc(void) {
    cfx_big_t a, b;
    cfx_big_init(&a);
    cfx_big_init(&b);
    cfx_big_from_hex(&b, "123456789ABCDEF0123456789ABCDEF0123456789ABCDEF0");
    int rc = cfx_big_copy(&a, &b);
    CFX_ASSERT(rc == 0);
    CFX_ASSERT(cfx_big_cmp(&a, &b) == 0);
    cfx_big_free(&b);
    cfx_big_free(&a);
}

static void test_copy_fast_path(void) {
    cfx_big_t a, b;
    cfx_big_init(&a);
    cfx_big_init(&b);
    cfx_big_reserve(&a, 100);
    cfx_big_from_limb(&b, 42);
    int rc = cfx_big_copy(&a, &b);
    CFX_ASSERT(rc == 0);
    CFX_ASSERT(cfx_big_eq_u64(&a, 42));
    cfx_big_free(&b);
    cfx_big_free(&a);
}

static void test_from_limb_clears_old(void) {
    cfx_big_t a;
    cfx_big_init(&a);
    cfx_big_from_hex(&a, "123456789ABCDEF0123456789ABCDEF0");
    CFX_ASSERT(a.n > 1);
    cfx_big_from_limb(&a, 42);
    CFX_ASSERT(a.n == 1);
    CFX_ASSERT(a.limb[0] == 42);
    cfx_big_free(&a);
}

static void test_cswap(void) {
    cfx_big_t a, b, a_orig, b_orig;
    cfx_big_init(&a);
    cfx_big_init(&b);
    cfx_big_init(&a_orig);
    cfx_big_init(&b_orig);

    /* cswap with condition=1 should swap */
    cfx_big_from_dec(&a, "123456789012345678901234567890");
    cfx_big_from_dec(&b, "987654321098765432109876543210");
    cfx_big_copy(&a_orig, &a);
    cfx_big_copy(&b_orig, &b);

    cfx_big_cswap(&a, &b, 1);
    CFX_ASSERT(cfx_big_eq(&a, &b_orig));
    CFX_ASSERT(cfx_big_eq(&b, &a_orig));

    /* condition=0 should NOT swap */
    cfx_big_copy(&a, &a_orig);
    cfx_big_copy(&b, &b_orig);

    cfx_big_cswap(&a, &b, 0);
    CFX_ASSERT(cfx_big_eq(&a, &a_orig));
    CFX_ASSERT(cfx_big_eq(&b, &b_orig));

    /* different sized operands */
    cfx_big_from_dec(&a, "42");  /* small */
    cfx_big_from_dec(&b, "99999999999999999999999999999999999999");  /* large */
    cfx_big_copy(&a_orig, &a);
    cfx_big_copy(&b_orig, &b);

    cfx_big_cswap(&a, &b, 1);
    CFX_ASSERT(cfx_big_eq(&a, &b_orig));
    CFX_ASSERT(cfx_big_eq(&b, &a_orig));

    /* same pointer (no-op) */
    cfx_big_from_dec(&a, "12345");
    cfx_big_copy(&a_orig, &a);
    cfx_big_cswap(&a, &a, 1);
    CFX_ASSERT(cfx_big_eq(&a, &a_orig));

    cfx_big_from_dec(&a, "0");
    cfx_big_from_dec(&b, "999");
    cfx_big_copy(&a_orig, &a);
    cfx_big_copy(&b_orig, &b);

    cfx_big_cswap(&a, &b, 1);
    CFX_ASSERT(cfx_big_eq(&a, &b_orig));
    CFX_ASSERT(cfx_big_eq(&b, &a_orig));

    cfx_big_free(&a);
    cfx_big_free(&b);
    cfx_big_free(&a_orig);
    cfx_big_free(&b_orig);
    PRINT_TEST(1);
}

int main(void) {
    CFX_TEST(test_cfx_big_init);
    CFX_TEST(test_cfx_big_reserve);
    CFX_TEST(test_cfx_big_assign);
    CFX_TEST(test_copy_swap);
    CFX_TEST(test_copy_self);
    CFX_TEST(test_assign_self);
    CFX_TEST(test_swap_self);
    CFX_TEST(test_move_self);
    CFX_TEST(test_assign_sm_zero);
    CFX_TEST(test_copy_zero);
    CFX_TEST(test_is_zero_variations);
    CFX_TEST(test_copy_needs_realloc);
    CFX_TEST(test_copy_fast_path);
    CFX_TEST(test_from_limb_clears_old);
    CFX_TEST(test_cswap);
    puts("OK");
    return 0;
}
