/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */
/* test_big_arith.c - Tests for basic big integer arithmetic */

#include "test_common.h"

static void test_mul_by_zero(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    cfx_big_from_limb(&b, 123);
    cfx_big_mul_sm_eq(&b, 2838);
    cfx_big_mul_sm_eq(&b, 1928);
    cfx_big_mul_sm_eq(&b, 9);
    cfx_big_mul_sm_eq(&b, 123765);
    cfx_big_mul_sm_eq(&b, 0);
    size_t sz = 0;
    char *s = cfx_big_dec_alloc(&b, &sz);
    CFX_ASSERT(sz == 1);
    CFX_ASSERT(strcmp(s, "0") == 0);
    free(s);
    cfx_big_free(&b);
    PRINT_TEST(1);
}

static void test_add_sm(void) {
    cfx_big_t b;
    cfx_big_init(&b);

    cfx_big_from_limb(&b, 123);
    cfx_big_add_sm_eq(&b, 321);
    CFX_ASSERT(b.limb[0] == 444);

    cfx_big_from_limb(&b, CFX_LIMB_MAX);
    CFX_ASSERT(b.limb[0] == CFX_LIMB_MAX);
    cfx_big_add_sm_eq(&b, 1);
    CFX_ASSERT(b.limb[0] == 0);
    CFX_ASSERT(b.limb[1] == 1);

    cfx_big_free(&b);
    PRINT_TEST(1);
}

static void test_sub_sm(void) {
    cfx_big_t b;
    cfx_big_init(&b);

    cfx_big_from_limb(&b, 1);
    cfx_big_sub_sm_eq(&b, 1);
    CFX_ASSERT(b.limb[0] == 0);
    cfx_big_sub_sm_eq(&b, 1);
    CFX_ASSERT(b.limb[0] == 0);

    cfx_big_from_limb(&b, CFX_LIMB_MAX);
    CFX_ASSERT(b.limb[0] == CFX_LIMB_MAX);
    cfx_big_add_sm_eq(&b, 1);
    CFX_ASSERT(b.limb[0] == 0);
    CFX_ASSERT(b.limb[1] == 1);

    cfx_big_sub_sm_eq(&b, 1);
    CFX_ASSERT(b.limb[0] == CFX_LIMB_MAX);
    CFX_ASSERT(b.n == 1);

    const cfx_limb_t N = 12;
    for (cfx_limb_t n = 0; n < N; ++n) {
        cfx_big_sub_sm_eq(&b, 1);
    }
    CFX_ASSERT(b.limb[0] == CFX_LIMB_MAX - N);

    cfx_limb_t q = 100;
    cfx_limb_t orig = b.limb[0];
    for (cfx_limb_t n = 0; n < 2; ++n) {
        cfx_acc_t s = (cfx_acc_t)b.limb[0] + q;
        cfx_big_add_sm_eq(&b, q);
        CFX_ASSERT(b.limb[0] == (cfx_limb_t)s);
        CFX_ASSERT(b.n > 1);
        CFX_ASSERT(b.limb[1] == (cfx_limb_t)(s >> CFX_LIMB_BITS));
        cfx_big_sub_sm_eq(&b, q);
        CFX_ASSERT(b.limb[0] == orig);
    }

    cfx_big_free(&b);
    PRINT_TEST(1);
}

static void test_sub(void) {
    cfx_big_t a, b;
    cfx_big_init(&a);
    cfx_big_init(&b);

    cfx_limb_t limbs[] = {0x1, 0x2, 0x333, 0x4444, 0x55555, 0x666666, 0x7777777, 0x88888888, 0xFF};
    cfx_big_from_limbs(&a, limbs, sizeof(limbs)/sizeof(limbs[0]));
    cfx_big_assign(&b, &a);

    cfx_big_t t;
    cfx_big_init(&t);
    cfx_big_copy(&t, &a);

    CFX_ASSERT(cfx_big_eq(&a, &b));
    CFX_ASSERT(cfx_big_eq(&a, &t));

    cfx_big_sub_eq(&t, &a);
    CFX_ASSERT(cfx_big_is_zero(&t));

    cfx_big_copy(&t, &a);
    cfx_big_add_sm_eq(&t, 0xFABADA);
    cfx_big_sub_eq(&t, &a);
    CFX_ASSERT(cfx_big_eq_u64(&t, 0xFABADA));

    cfx_big_copy(&t, &a);
    cfx_big_mul_sm_eq(&t, 2);
    cfx_big_sub_eq(&t, &a);
    CFX_ASSERT(cfx_big_eq(&t, &a));

    cfx_big_t two;
    cfx_big_init(&two);
    cfx_big_from_limb(&two, 2);

    cfx_big_copy(&t, &a);
    cfx_big_mul_eq(&t, &two);
    cfx_big_sub_eq(&t, &a);
    CFX_ASSERT(cfx_big_eq(&t, &a));

    cfx_big_free(&a);
    cfx_big_free(&b);
    cfx_big_free(&t);
    cfx_big_free(&two);
    PRINT_TEST(1);
}

static void test_zero_right(void) {
    cfx_big_t b, m;
    cfx_big_init(&b);
    cfx_big_from_limb(&b, 123);
    cfx_big_init(&m);
    cfx_big_from_limb(&m, 0);

    cfx_big_mul_eq(&b, &m);
    CFX_ASSERT(cfx_big_is_zero(&b));

    cfx_big_free(&b);
    cfx_big_free(&m);
    PRINT_TEST(1);
}

static void test_zero_left(void) {
    cfx_big_t b, m;
    cfx_big_init(&b);
    cfx_big_from_limb(&b, 0);
    cfx_big_init(&m);
    cfx_big_from_limb(&m, 123);

    cfx_big_mul_eq(&b, &m);
    CFX_ASSERT(cfx_big_is_zero(&b));

    cfx_big_from_limb(&b, 123);
    cfx_big_from_limb(&m, 0);
    cfx_big_mul_eq(&b, &m);
    CFX_ASSERT(cfx_big_is_zero(&b));

    cfx_big_free(&b);
    cfx_big_free(&m);
    PRINT_TEST(1);
}

static void test_cmp_sm_branches(void) {
    cfx_big_t a;
    cfx_big_init(&a);

    CFX_ASSERT(cfx_big_cmp_sm(&a, 0) == 0);   /* a==0, b==0 */
    CFX_ASSERT(cfx_big_cmp_sm(&a, 5) == -1);  /* a==0, b>0 */

    cfx_big_from_hex(&a, "10000000000000000");  /* > 64 bits */
    CFX_ASSERT(cfx_big_cmp_sm(&a, CFX_LIMB_MAX) == 1);  /* a > 1 limb */

    cfx_big_from_limb(&a, 5);
    CFX_ASSERT(cfx_big_cmp_sm(&a, 10) == -1);  /* a < b */
    CFX_ASSERT(cfx_big_cmp_sm(&a, 3) == 1);    /* a > b */
    CFX_ASSERT(cfx_big_cmp_sm(&a, 5) == 0);    /* a == b */

    cfx_big_free(&a);
}

static void test_cmp_u64(void) {
    cfx_big_t a;
    cfx_big_init(&a);

    /* zero cases */
    CFX_ASSERT(cfx_big_cmp_u64(&a, 0) == 0);
    CFX_ASSERT(cfx_big_cmp_u64(&a, 1) == -1);
    CFX_ASSERT(cfx_big_cmp_u64(&a, UINT64_MAX) == -1);

    /* single limb */
    cfx_big_from_u64(&a, 100);
    CFX_ASSERT(cfx_big_cmp_u64(&a, 100) == 0);
    CFX_ASSERT(cfx_big_cmp_u64(&a, 99) == 1);
    CFX_ASSERT(cfx_big_cmp_u64(&a, 101) == -1);

    /* large u64 value */
    cfx_big_from_u64(&a, 0xFFFFFFFFFFFFFFFFULL);
    CFX_ASSERT(cfx_big_cmp_u64(&a, UINT64_MAX) == 0);
    CFX_ASSERT(cfx_big_cmp_u64(&a, UINT64_MAX - 1) == 1);

    /* bigger than any u64 */
    cfx_big_from_hex(&a, "10000000000000000");  /* 2^64 */
    CFX_ASSERT(cfx_big_cmp_u64(&a, UINT64_MAX) == 1);
    CFX_ASSERT(cfx_big_cmp_u64(&a, 0) == 1);

    cfx_big_free(&a);
    PRINT_TEST(1);
}

static void test_sub_u64(void) {
    cfx_big_t a;
    cfx_big_init(&a);

    /* basic subtraction */
    cfx_big_from_u64(&a, 1000);
    cfx_big_sub_u64_eq(&a, 100);
    CFX_ASSERT(cfx_big_cmp_u64(&a, 900) == 0);

    /* subtract to zero */
    cfx_big_from_u64(&a, 42);
    cfx_big_sub_u64_eq(&a, 42);
    CFX_ASSERT(cfx_big_is_zero(&a));

    /* subtract zero does nothing */
    cfx_big_from_u64(&a, 123);
    cfx_big_sub_u64_eq(&a, 0);
    CFX_ASSERT(cfx_big_cmp_u64(&a, 123) == 0);

    /* borrow across limbs: 2^64 - 1 + 6 = 2^64 + 5, then subtract 10 => 2^64 - 5 */
    cfx_big_from_u64(&a, UINT64_MAX);
    cfx_big_add_sm_eq(&a, 1);  /* now 2^64 */
    cfx_big_sub_u64_eq(&a, 1);
    CFX_ASSERT(cfx_big_cmp_u64(&a, UINT64_MAX) == 0);

    cfx_big_free(&a);
    PRINT_TEST(1);
}

static void test_eq_u64_branches(void) {
    cfx_big_t a;
    cfx_big_init(&a);

    CFX_ASSERT(cfx_big_eq_u64(&a, 0));   /* both zero */
    CFX_ASSERT(!cfx_big_eq_u64(&a, 5));  /* a zero, n nonzero */

    cfx_big_from_limb(&a, 5);
    CFX_ASSERT(!cfx_big_eq_u64(&a, 0));  /* a nonzero, n zero */
    CFX_ASSERT(cfx_big_eq_u64(&a, 5));   /* equal */
    CFX_ASSERT(!cfx_big_eq_u64(&a, 6));  /* not equal */

    cfx_big_free(&a);
}

static void test_is_even_zero(void) {
    cfx_big_t a;
    cfx_big_init(&a);
    CFX_ASSERT(cfx_big_is_even(&a));  /* zero is even */
    cfx_big_free(&a);
}

static void test_add_non_inplace(void) {
    cfx_big_t a, b, out;
    cfx_big_init(&a);
    cfx_big_init(&b);
    cfx_big_init(&out);

    /* basic addition */
    cfx_big_from_limb(&a, 100);
    cfx_big_from_limb(&b, 200);
    cfx_big_add(&out, &a, &b);
    CFX_ASSERT(cfx_big_eq_u64(&out, 300));
    CFX_ASSERT(cfx_big_eq_u64(&a, 100));  /* a unchanged */
    CFX_ASSERT(cfx_big_eq_u64(&b, 200));  /* b unchanged */

    /* addition with carry */
    cfx_big_from_limb(&a, CFX_LIMB_MAX);
    cfx_big_from_limb(&b, 1);
    cfx_big_add(&out, &a, &b);
    CFX_ASSERT(out.n == 2);
    CFX_ASSERT(out.limb[0] == 0);
    CFX_ASSERT(out.limb[1] == 1);

    /* aliasing: out == a */
    cfx_big_from_limb(&a, 50);
    cfx_big_from_limb(&b, 25);
    cfx_big_add(&a, &a, &b);
    CFX_ASSERT(cfx_big_eq_u64(&a, 75));

    /* aliasing: out == b */
    cfx_big_from_limb(&a, 50);
    cfx_big_from_limb(&b, 25);
    cfx_big_add(&b, &a, &b);
    CFX_ASSERT(cfx_big_eq_u64(&b, 75));

    /* large numbers */
    cfx_big_from_hex(&a, "ffffffffffffffffffffffffffffffff");
    cfx_big_from_hex(&b, "1");
    cfx_big_add(&out, &a, &b);
    char *s = cfx_big_hex_alloc(&out, NULL);
    CFX_ASSERT(strcmp(s, "100000000000000000000000000000000") == 0);
    free(s);

    cfx_big_free(&a);
    cfx_big_free(&b);
    cfx_big_free(&out);
    PRINT_TEST(1);
}

static void test_sub_non_inplace(void) {
    cfx_big_t a, b, out;
    cfx_big_init(&a);
    cfx_big_init(&b);
    cfx_big_init(&out);

    /* basic subtraction */
    cfx_big_from_limb(&a, 300);
    cfx_big_from_limb(&b, 200);
    cfx_big_sub(&out, &a, &b);
    CFX_ASSERT(cfx_big_eq_u64(&out, 100));
    CFX_ASSERT(cfx_big_eq_u64(&a, 300));  /* a unchanged */
    CFX_ASSERT(cfx_big_eq_u64(&b, 200));  /* b unchanged */

    /* subtraction with borrow */
    cfx_big_t big;
    cfx_big_init(&big);
    cfx_big_from_limb(&big, 1);
    cfx_big_shl_bits_eq(&big, CFX_LIMB_BITS);  /* big = 2^64 or 2^32 */
    cfx_big_from_limb(&b, 1);
    cfx_big_sub(&out, &big, &b);
    CFX_ASSERT(out.n == 1);
    CFX_ASSERT(out.limb[0] == CFX_LIMB_MAX);
    cfx_big_free(&big);

    /* a - a = 0 */
    cfx_big_from_limb(&a, 12345);
    cfx_big_sub(&out, &a, &a);
    CFX_ASSERT(cfx_big_is_zero(&out));

    /* aliasing: out == a */
    cfx_big_from_limb(&a, 100);
    cfx_big_from_limb(&b, 30);
    cfx_big_sub(&a, &a, &b);
    CFX_ASSERT(cfx_big_eq_u64(&a, 70));

    /* aliasing: out == b */
    cfx_big_from_limb(&a, 100);
    cfx_big_from_limb(&b, 30);
    cfx_big_sub(&b, &a, &b);
    CFX_ASSERT(cfx_big_eq_u64(&b, 70));

    /* large numbers */
    cfx_big_from_hex(&a, "100000000000000000000000000000000");
    cfx_big_from_hex(&b, "1");
    cfx_big_sub(&out, &a, &b);
    char *s = cfx_big_hex_alloc(&out, NULL);
    CFX_ASSERT(strcmp(s, "ffffffffffffffffffffffffffffffff") == 0);
    free(s);

    cfx_big_free(&a);
    cfx_big_free(&b);
    cfx_big_free(&out);
    PRINT_TEST(1);
}

int main(void) {
    CFX_TEST(test_mul_by_zero);
    CFX_TEST(test_add_sm);
    CFX_TEST(test_sub_sm);
    CFX_TEST(test_sub);
    CFX_TEST(test_zero_right);
    CFX_TEST(test_zero_left);
    CFX_TEST(test_cmp_sm_branches);
    CFX_TEST(test_cmp_u64);
    CFX_TEST(test_sub_u64);
    CFX_TEST(test_eq_u64_branches);
    CFX_TEST(test_is_even_zero);
    CFX_TEST(test_add_non_inplace);
    CFX_TEST(test_sub_non_inplace);
    puts("OK");
    return 0;
}
