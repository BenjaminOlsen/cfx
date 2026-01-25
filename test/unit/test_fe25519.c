#include "cfx/fe25519.h"
#include "cfx/macros.h"

#include <stdio.h>
#include <string.h>

/* helper to print a field element for debugging */
static void fe_print(const char *name, const fe25519_t *f) {
    printf("%s: [%016llx, %016llx, %016llx, %016llx, %016llx]\n",
        name,
        (unsigned long long)f->v[0],
        (unsigned long long)f->v[1],
        (unsigned long long)f->v[2],
        (unsigned long long)f->v[3],
        (unsigned long long)f->v[4]);
}

static void fe_print_bytes(const char *name, const uint8_t s[32]) {
    printf("%s: ", name);
    for (int i = 0; i < 32; i++) printf("%02x", s[i]);
    printf("\n");
}

static int fe_limbs_eq(const fe25519_t *a, const fe25519_t *b) {
    return a->v[0] == b->v[0] &&
           a->v[1] == b->v[1] &&
           a->v[2] == b->v[2] &&
           a->v[3] == b->v[3] &&
           a->v[4] == b->v[4];
}

/* ============================================================ */
/* Basic operations */
/* ============================================================ */

static void test_zero(void) {
    fe25519_t z;
    cfx_fe25519_0(&z);
    CFX_ASSERT(z.v[0] == 0);
    CFX_ASSERT(z.v[1] == 0);
    CFX_ASSERT(z.v[2] == 0);
    CFX_ASSERT(z.v[3] == 0);
    CFX_ASSERT(z.v[4] == 0);
    CFX_ASSERT(cfx_fe25519_iszero(&z) == 1);
}

static void test_one(void) {
    fe25519_t one;
    cfx_fe25519_1(&one);
    CFX_ASSERT(one.v[0] == 1);
    CFX_ASSERT(one.v[1] == 0);
    CFX_ASSERT(one.v[2] == 0);
    CFX_ASSERT(one.v[3] == 0);
    CFX_ASSERT(one.v[4] == 0);
    CFX_ASSERT(cfx_fe25519_iszero(&one) == 0);
}

static void test_copy(void) {
    fe25519_t a = {{ 0x123, 0x456, 0x789, 0xabc, 0xdef }};
    fe25519_t b;
    cfx_fe25519_copy(&b, &a);
    CFX_ASSERT(fe_limbs_eq(&a, &b));
}

static void test_constants(void) {
    CFX_ASSERT(cfx_fe25519_iszero(&cfx_fe25519_zero) == 1);
    CFX_ASSERT(cfx_fe25519_iszero(&cfx_fe25519_one) == 0);
    CFX_ASSERT(cfx_fe25519_eq(&cfx_fe25519_zero, &cfx_fe25519_zero) == 1);
    CFX_ASSERT(cfx_fe25519_eq(&cfx_fe25519_one, &cfx_fe25519_one) == 1);
    CFX_ASSERT(cfx_fe25519_eq(&cfx_fe25519_zero, &cfx_fe25519_one) == 0);
}

/* ============================================================ */
/* Addition and subtraction */
/* ============================================================ */

static void test_add_basic(void) {
    fe25519_t a, b, c;
    cfx_fe25519_1(&a);
    cfx_fe25519_1(&b);
    cfx_fe25519_add(&c, &a, &b);

    CFX_ASSERT(c.v[0] == 2);
    CFX_ASSERT(c.v[1] == 0);
}

static void test_add_zero(void) {
    fe25519_t a = {{ 0x1234567890abc, 0x7edcba9876543, 0, 0, 0 }};
    fe25519_t z;
    cfx_fe25519_0(&z);
    fe25519_t c;
    cfx_fe25519_add(&c, &a, &z);
    CFX_ASSERT(fe_limbs_eq(&a, &c));
}

static void test_sub_basic(void) {
    fe25519_t two = {{ 2, 0, 0, 0, 0 }};
    fe25519_t one;
    cfx_fe25519_1(&one);
    fe25519_t c;
    cfx_fe25519_sub(&c, &two, &one);
    cfx_fe25519_carry(&c);
    CFX_ASSERT(cfx_fe25519_eq(&c, &one));
}

static void test_sub_self(void) {
    fe25519_t a = {{ 0x123456789, 0xabcdef, 0x111, 0x222, 0x333 }};
    fe25519_t c;
    cfx_fe25519_sub(&c, &a, &a);
    CFX_ASSERT(cfx_fe25519_iszero(&c));
}

static void test_neg_basic(void) {
    fe25519_t one, neg_one, sum;
    cfx_fe25519_1(&one);
    cfx_fe25519_neg(&neg_one, &one);
    cfx_fe25519_add(&sum, &one, &neg_one);
    CFX_ASSERT(cfx_fe25519_iszero(&sum));
}

static void test_neg_zero(void) {
    fe25519_t z, neg_z;
    cfx_fe25519_0(&z);
    cfx_fe25519_neg(&neg_z, &z);
    CFX_ASSERT(cfx_fe25519_iszero(&neg_z));
}

/* ============================================================ */
/* Multiplication */
/* ============================================================ */

static void test_mul_by_zero(void) {
    fe25519_t a = {{ 0x123, 0x456, 0x789, 0xabc, 0xdef }};
    fe25519_t z, c;
    cfx_fe25519_0(&z);
    cfx_fe25519_mul(&c, &a, &z);
    CFX_ASSERT(cfx_fe25519_iszero(&c));
}

static void test_mul_by_one(void) {
    fe25519_t a = {{ 0x123456789, 0xabcdef, 0x111, 0x222, 0x333 }};
    fe25519_t one, c;
    cfx_fe25519_1(&one);
    cfx_fe25519_mul(&c, &a, &one);
    CFX_ASSERT(cfx_fe25519_eq(&a, &c));
}

static void test_mul_small(void) {
    fe25519_t two = {{ 2, 0, 0, 0, 0 }};
    fe25519_t three = {{ 3, 0, 0, 0, 0 }};
    fe25519_t six = {{ 6, 0, 0, 0, 0 }};
    fe25519_t c;
    cfx_fe25519_mul(&c, &two, &three);
    CFX_ASSERT(cfx_fe25519_eq(&c, &six));
}

static void test_mul_commutative(void) {
    /* limbs must be < 2^51 = 0x8000000000000 */
    fe25519_t a = {{ 0x1234567890abc, 0x7edcba9876543, 0x111, 0x222, 0x333 }};
    fe25519_t b = {{ 0xaabbccdd, 0x11223344, 0x55667788, 0x4aabbcc, 0xddeeff }};
    fe25519_t ab, ba;
    cfx_fe25519_mul(&ab, &a, &b);
    cfx_fe25519_mul(&ba, &b, &a);
    CFX_ASSERT(cfx_fe25519_eq(&ab, &ba));
}

static void test_mul_associative(void) {
    fe25519_t a = {{ 0x1234, 0x5678, 0x9abc, 0xdef0, 0x1111 }};
    fe25519_t b = {{ 0xaaaa, 0xbbbb, 0xcccc, 0xdddd, 0xeeee }};
    fe25519_t c = {{ 0x2222, 0x3333, 0x4444, 0x5555, 0x6666 }};
    fe25519_t ab, ab_c, bc, a_bc;

    cfx_fe25519_mul(&ab, &a, &b);
    cfx_fe25519_mul(&ab_c, &ab, &c);

    cfx_fe25519_mul(&bc, &b, &c);
    cfx_fe25519_mul(&a_bc, &a, &bc);

    CFX_ASSERT(cfx_fe25519_eq(&ab_c, &a_bc));
}

static void test_mul_distributive(void) {
    fe25519_t a = {{ 0x1234, 0x5678, 0x9abc, 0, 0 }};
    fe25519_t b = {{ 0xaaaa, 0xbbbb, 0, 0, 0 }};
    fe25519_t c = {{ 0x2222, 0x3333, 0, 0, 0 }};
    fe25519_t b_plus_c, a_times_bc, ab, ac, ab_plus_ac;

    cfx_fe25519_add(&b_plus_c, &b, &c);
    cfx_fe25519_mul(&a_times_bc, &a, &b_plus_c);

    cfx_fe25519_mul(&ab, &a, &b);
    cfx_fe25519_mul(&ac, &a, &c);
    cfx_fe25519_add(&ab_plus_ac, &ab, &ac);

    CFX_ASSERT(cfx_fe25519_eq(&a_times_bc, &ab_plus_ac));
}

/* ============================================================ */
/* Squaring */
/* ============================================================ */

static void test_sqr_zero(void) {
    fe25519_t z, z2;
    cfx_fe25519_0(&z);
    cfx_fe25519_sqr(&z2, &z);
    CFX_ASSERT(cfx_fe25519_iszero(&z2));
}

static void test_sqr_one(void) {
    fe25519_t one, one2;
    cfx_fe25519_1(&one);
    cfx_fe25519_sqr(&one2, &one);
    CFX_ASSERT(cfx_fe25519_eq(&one2, &one));
}

static void test_sqr_vs_mul(void) {
    fe25519_t a = {{ 0x123456789, 0xabcdef123, 0x456789abc, 0xdef012345, 0x6789abcde }};
    fe25519_t sqr_a, mul_aa;
    cfx_fe25519_sqr(&sqr_a, &a);
    cfx_fe25519_mul(&mul_aa, &a, &a);
    CFX_ASSERT(cfx_fe25519_eq(&sqr_a, &mul_aa));
}

static void test_sqr_n(void) {
    fe25519_t a = {{ 0x12345, 0x67890, 0xabcde, 0xf0123, 0x45678 }};
    fe25519_t sqr_n, sqr_loop;

    cfx_fe25519_sqr_n(&sqr_n, &a, 5);

    cfx_fe25519_copy(&sqr_loop, &a);
    for (int i = 0; i < 5; i++) {
        cfx_fe25519_sqr(&sqr_loop, &sqr_loop);
    }

    CFX_ASSERT(cfx_fe25519_eq(&sqr_n, &sqr_loop));
}

/* ============================================================ */
/* mul121666 (used in X25519) */
/* ============================================================ */

static void test_mul121666_zero(void) {
    fe25519_t z, c;
    cfx_fe25519_0(&z);
    cfx_fe25519_mul121666(&c, &z);
    CFX_ASSERT(cfx_fe25519_iszero(&c));
}

static void test_mul121666_one(void) {
    fe25519_t one, c, expected;
    cfx_fe25519_1(&one);
    cfx_fe25519_mul121666(&c, &one);
    expected = (fe25519_t){{ 121666, 0, 0, 0, 0 }};
    CFX_ASSERT(cfx_fe25519_eq(&c, &expected));
}

static void test_mul121666_vs_mul(void) {
    fe25519_t a = {{ 0x123456789, 0xabcdef, 0x456, 0x789, 0xabc }};
    fe25519_t k = {{ 121666, 0, 0, 0, 0 }};
    fe25519_t via_mul121666, via_mul;

    cfx_fe25519_mul121666(&via_mul121666, &a);
    cfx_fe25519_mul(&via_mul, &a, &k);

    CFX_ASSERT(cfx_fe25519_eq(&via_mul121666, &via_mul));
}

/* ============================================================ */
/* Inversion */
/* ============================================================ */

static void test_inv_one(void) {
    fe25519_t one, inv_one;
    cfx_fe25519_1(&one);
    cfx_fe25519_inv(&inv_one, &one);
    CFX_ASSERT(cfx_fe25519_eq(&inv_one, &one));
}

static void test_inv_mul_identity(void) {
    fe25519_t a = {{ 0x1234567890abc, 0x7edcba9876543, 0x111222333, 0x444555666, 0x777888 }};
    fe25519_t inv_a, product, one;

    cfx_fe25519_inv(&inv_a, &a);
    cfx_fe25519_mul(&product, &a, &inv_a);
    cfx_fe25519_1(&one);

    CFX_ASSERT(cfx_fe25519_eq(&product, &one));
}

static void test_inv_inv(void) {
    fe25519_t a = {{ 0xabcdef, 0x123456, 0x789abc, 0xdef012, 0x345678 }};
    fe25519_t inv_a, inv_inv_a;

    cfx_fe25519_inv(&inv_a, &a);
    cfx_fe25519_inv(&inv_inv_a, &inv_a);

    CFX_ASSERT(cfx_fe25519_eq(&a, &inv_inv_a));
}

/* ============================================================ */
/* Serialization */
/* ============================================================ */

static void test_bytes_roundtrip(void) {
    uint8_t original[32] = {
        0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08,
        0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f, 0x10,
        0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17, 0x18,
        0x19, 0x1a, 0x1b, 0x1c, 0x1d, 0x1e, 0x1f, 0x20
    };

    fe25519_t f;
    uint8_t output[32];

    cfx_fe25519_frombytes(&f, original);
    cfx_fe25519_tobytes(output, &f);

    CFX_ASSERT(memcmp(original, output, 32) == 0);
}

static void test_bytes_zero(void) {
    uint8_t zero_bytes[32] = {0};
    fe25519_t f;
    uint8_t output[32];

    cfx_fe25519_frombytes(&f, zero_bytes);
    CFX_ASSERT(cfx_fe25519_iszero(&f));

    cfx_fe25519_tobytes(output, &f);
    CFX_ASSERT(memcmp(zero_bytes, output, 32) == 0);
}

static void test_bytes_one(void) {
    uint8_t one_bytes[32] = { 1, 0 };
    fe25519_t f, one;
    uint8_t output[32];

    cfx_fe25519_frombytes(&f, one_bytes);
    cfx_fe25519_1(&one);
    CFX_ASSERT(cfx_fe25519_eq(&f, &one));

    cfx_fe25519_tobytes(output, &one);
    CFX_ASSERT(memcmp(one_bytes, output, 32) == 0);
}

static void test_bytes_clears_top_bit(void) {
    uint8_t with_top_bit[32] = {0};
    with_top_bit[31] = 0x80;  /* set bit 255 */

    fe25519_t f;
    uint8_t output[32];

    cfx_fe25519_frombytes(&f, with_top_bit);
    cfx_fe25519_tobytes(output, &f);

    /* bit 255 should be cleared */
    CFX_ASSERT((output[31] & 0x80) == 0);
}

/* ============================================================ */
/* Constant-time operations */
/* ============================================================ */

static void test_cswap_zero(void) {
    fe25519_t a = {{ 1, 2, 3, 4, 5 }};
    fe25519_t b = {{ 10, 20, 30, 40, 50 }};
    fe25519_t a_orig, b_orig;
    cfx_fe25519_copy(&a_orig, &a);
    cfx_fe25519_copy(&b_orig, &b);

    cfx_fe25519_cswap(&a, &b, 0);  /* no swap */

    CFX_ASSERT(fe_limbs_eq(&a, &a_orig));
    CFX_ASSERT(fe_limbs_eq(&b, &b_orig));
}

static void test_cswap_one(void) {
    fe25519_t a = {{ 1, 2, 3, 4, 5 }};
    fe25519_t b = {{ 10, 20, 30, 40, 50 }};
    fe25519_t a_orig, b_orig;
    cfx_fe25519_copy(&a_orig, &a);
    cfx_fe25519_copy(&b_orig, &b);

    cfx_fe25519_cswap(&a, &b, 1);  /* swap */

    CFX_ASSERT(fe_limbs_eq(&a, &b_orig));
    CFX_ASSERT(fe_limbs_eq(&b, &a_orig));
}

static void test_cmov_zero(void) {
    fe25519_t h = {{ 1, 2, 3, 4, 5 }};
    fe25519_t f = {{ 10, 20, 30, 40, 50 }};
    fe25519_t h_orig;
    cfx_fe25519_copy(&h_orig, &h);

    cfx_fe25519_cmov(&h, &f, 0);  /* no move */

    CFX_ASSERT(fe_limbs_eq(&h, &h_orig));
}

static void test_cmov_one(void) {
    fe25519_t h = {{ 1, 2, 3, 4, 5 }};
    fe25519_t f = {{ 10, 20, 30, 40, 50 }};

    cfx_fe25519_cmov(&h, &f, 1);  /* move */

    CFX_ASSERT(fe_limbs_eq(&h, &f));
}

/* ============================================================ */
/* Comparison and predicates */
/* ============================================================ */

static void test_eq_same(void) {
    fe25519_t a = {{ 0x123, 0x456, 0x789, 0xabc, 0xdef }};
    CFX_ASSERT(cfx_fe25519_eq(&a, &a) == 1);
}

static void test_eq_different(void) {
    fe25519_t a = {{ 0x123, 0x456, 0x789, 0xabc, 0xdef }};
    fe25519_t b = {{ 0x123, 0x456, 0x789, 0xabc, 0xdee }};  /* last limb different */
    CFX_ASSERT(cfx_fe25519_eq(&a, &b) == 0);
}

static void test_eq_reduced_vs_unreduced(void) {
    /* p (the prime) should equal 0 */
    const uint64_t mask51 = (1ULL << 51) - 1;
    fe25519_t p = {{ mask51 - 18, mask51, mask51, mask51, mask51 }};
    fe25519_t zero;
    cfx_fe25519_0(&zero);
    /* p ≡ 0 (mod p) */
    CFX_ASSERT(cfx_fe25519_eq(&p, &zero) == 1);

    /* also: 1 + p should equal 1 */
    fe25519_t one;
    cfx_fe25519_1(&one);
    fe25519_t one_plus_p;
    cfx_fe25519_add(&one_plus_p, &one, &p);
    CFX_ASSERT(cfx_fe25519_eq(&one_plus_p, &one) == 1);
}

static void test_isnegative(void) {
    fe25519_t zero;
    cfx_fe25519_0(&zero);
    CFX_ASSERT(cfx_fe25519_isnegative(&zero) == 0);

    fe25519_t one;
    cfx_fe25519_1(&one);
    CFX_ASSERT(cfx_fe25519_isnegative(&one) == 1);  /* 1 is "negative" (LSB = 1) */

    fe25519_t two = {{ 2, 0, 0, 0, 0 }};
    CFX_ASSERT(cfx_fe25519_isnegative(&two) == 0);  /* 2 is "positive" (LSB = 0) */
}

/* ============================================================ */
/* Reduction */
/* ============================================================ */

static void test_reduce_small(void) {
    fe25519_t a = {{ 100, 200, 300, 400, 500 }};
    fe25519_t b;
    cfx_fe25519_copy(&b, &a);
    cfx_fe25519_reduce(&b);
    CFX_ASSERT(cfx_fe25519_eq(&a, &b));
}

static void test_reduce_large_limb(void) {
    /* limb 0 has value requiring carry */
    fe25519_t a = {{ (1ULL << 52), 0, 0, 0, 0 }};
    cfx_fe25519_reduce(&a);
    /* after reduction, limbs should be < 2^51 */
    CFX_ASSERT(a.v[0] < (1ULL << 51));
}

static void test_carry_propagation(void) {
    /* cascading carries */
    fe25519_t a = {{
                       (1ULL << 52) - 1,
                       (1ULL << 52) - 1,
                       (1ULL << 52) - 1,
                       (1ULL << 52) - 1,
                       (1ULL << 52) - 1
                   }};
    cfx_fe25519_carry(&a);
    cfx_fe25519_carry(&a);

    /* all limbs should now be small */
    for (int i = 0; i < 5; i++) {
        CFX_ASSERT(a.v[i] < (1ULL << 52));
    }
}

/* ============================================================ */
/* RFC 7748 test vectors */
/* ============================================================ */

/*
 * From RFC 7748 Section 6.1:
 * basepoint u = 9
 * For X25519, the u-coordinate of the base point is 9.
 */
static void test_basepoint(void) {
    uint8_t basepoint_bytes[32] = { 9 };
    fe25519_t bp;
    uint8_t output[32];

    cfx_fe25519_frombytes(&bp, basepoint_bytes);
    cfx_fe25519_tobytes(output, &bp);

    CFX_ASSERT(memcmp(basepoint_bytes, output, 32) == 0);
}

/*
 * Test that (a + b) - b = a for random-ish values
 */
static void test_add_sub_inverse(void) {
    fe25519_t a = {{ 0x1234567890abc, 0x7edcba9876543, 0xaabbccdd, 0x11223344, 0x55667788 }};
    fe25519_t b = {{ 0x1111111111111, 0x2222222222222, 0x333333333, 0x444444444, 0x555555555 }};
    fe25519_t sum, diff;

    cfx_fe25519_add(&sum, &a, &b);
    cfx_fe25519_sub(&diff, &sum, &b);

    CFX_ASSERT(cfx_fe25519_eq(&a, &diff));
}

/*
 * Test that a * b / b = a (using inversion)
 */
static void test_mul_inv_identity(void) {
    fe25519_t a = {{ 0x12345, 0x67890, 0xabcde, 0xf0123, 0x45678 }};
    fe25519_t b = {{ 0x11111, 0x22222, 0x33333, 0x44444, 0x55555 }};
    fe25519_t product, inv_b, result;

    cfx_fe25519_mul(&product, &a, &b);
    cfx_fe25519_inv(&inv_b, &b);
    cfx_fe25519_mul(&result, &product, &inv_b);

    CFX_ASSERT(cfx_fe25519_eq(&a, &result));
}

/* ============================================================ */
/* Edge cases */
/* ============================================================ */

static void test_max_limb_values(void) {
    /* limbs at maximum 51-bit value */
    const uint64_t mask51 = (1ULL << 51) - 1;
    fe25519_t a = {{ mask51, mask51, mask51, mask51, mask51 }};
    fe25519_t b = {{ mask51, mask51, mask51, mask51, mask51 }};
    fe25519_t c;

    cfx_fe25519_mul(&c, &a, &b);
    /* should not crash, result should be valid */
    cfx_fe25519_reduce(&c);

    for (int i = 0; i < 5; i++) {
        CFX_ASSERT(c.v[i] < (1ULL << 51));
    }
}

static void test_sqr_large(void) {
    const uint64_t mask51 = (1ULL << 51) - 1;
    fe25519_t a = {{ mask51, mask51, mask51, mask51, mask51 }};
    fe25519_t c;

    cfx_fe25519_sqr(&c, &a);
    cfx_fe25519_reduce(&c);

    for (int i = 0; i < 5; i++) {
        CFX_ASSERT(c.v[i] < (1ULL << 51));
    }
}

/* Test that p ≡ 0 (mod p) */
static void test_p_is_zero(void) {
    /* p = 2^255 - 19 in radix 2^51:
     * limb 0: 2^51 - 19
     * limbs 1-4: 2^51 - 1
     */
    const uint64_t mask51 = (1ULL << 51) - 1;
    fe25519_t p = {{ mask51 - 18, mask51, mask51, mask51, mask51 }};

    CFX_ASSERT(cfx_fe25519_iszero(&p) == 1);
}

/* Test 2p ≡ 0 (mod p) */
static void test_2p_is_zero(void) {
    const uint64_t mask51 = (1ULL << 51) - 1;
    fe25519_t p = {{ mask51 - 18, mask51, mask51, mask51, mask51 }};
    fe25519_t two_p;

    cfx_fe25519_add(&two_p, &p, &p);
    CFX_ASSERT(cfx_fe25519_iszero(&two_p) == 1);
}

/* ============================================================ */
/* pow22523 (used for square roots) */
/* ============================================================ */

static void test_pow22523_basic(void) {
    fe25519_t a = {{ 0x12345, 0x67890, 0xabcde, 0, 0 }};
    fe25519_t result;

    /* pow22523 computes a^((p-5)/8) */
    cfx_fe25519_pow22523(&result, &a);

    /* just verify it doesn't crash and produces valid output */
    cfx_fe25519_reduce(&result);
    for (int i = 0; i < 5; i++) {
        CFX_ASSERT(result.v[i] < (1ULL << 51));
    }
}

/* ============================================================ */
/* Aliasing tests (output == input) */
/* ============================================================ */

static void test_add_aliased(void) {
    fe25519_t a = {{ 0x1234, 0x5678, 0x9abc, 0xdef0, 0x1111 }};
    fe25519_t a_copy;
    cfx_fe25519_copy(&a_copy, &a);
    fe25519_t expected;
    cfx_fe25519_add(&expected, &a, &a);

    cfx_fe25519_add(&a, &a, &a);  /* a = a + a */
    CFX_ASSERT(cfx_fe25519_eq(&a, &expected));
}

static void test_sub_aliased(void) {
    fe25519_t a = {{ 0x1234, 0x5678, 0x9abc, 0xdef0, 0x1111 }};
    fe25519_t b = {{ 0x111, 0x222, 0x333, 0x444, 0x555 }};
    fe25519_t expected;
    cfx_fe25519_sub(&expected, &a, &b);

    cfx_fe25519_sub(&a, &a, &b);  /* a = a - b */
    CFX_ASSERT(cfx_fe25519_eq(&a, &expected));
}

static void test_mul_aliased(void) {
    fe25519_t a = {{ 0x1234, 0x5678, 0x9abc, 0xdef0, 0x1111 }};
    fe25519_t b = {{ 0xaaaa, 0xbbbb, 0xcccc, 0xdddd, 0xeeee }};
    fe25519_t expected;
    cfx_fe25519_mul(&expected, &a, &b);

    cfx_fe25519_mul(&a, &a, &b);  /* a = a * b */
    CFX_ASSERT(cfx_fe25519_eq(&a, &expected));
}

static void test_mul_self_aliased(void) {
    fe25519_t a = {{ 0x1234, 0x5678, 0x9abc, 0xdef0, 0x1111 }};
    fe25519_t expected;
    cfx_fe25519_mul(&expected, &a, &a);

    cfx_fe25519_mul(&a, &a, &a);  /* a = a * a */
    CFX_ASSERT(cfx_fe25519_eq(&a, &expected));
}

static void test_sqr_aliased(void) {
    fe25519_t a = {{ 0x1234, 0x5678, 0x9abc, 0xdef0, 0x1111 }};
    fe25519_t expected;
    cfx_fe25519_sqr(&expected, &a);

    cfx_fe25519_sqr(&a, &a);  /* a = a^2 */
    CFX_ASSERT(cfx_fe25519_eq(&a, &expected));
}

static void test_neg_aliased(void) {
    fe25519_t a = {{ 0x1234, 0x5678, 0x9abc, 0xdef0, 0x1111 }};
    fe25519_t expected;
    cfx_fe25519_neg(&expected, &a);

    cfx_fe25519_neg(&a, &a);  /* a = -a */
    CFX_ASSERT(cfx_fe25519_eq(&a, &expected));
}

static void test_inv_aliased(void) {
    fe25519_t a = {{ 0x1234, 0x5678, 0x9abc, 0xdef0, 0x1111 }};
    fe25519_t expected;
    cfx_fe25519_inv(&expected, &a);

    cfx_fe25519_inv(&a, &a);  /* a = 1/a */
    CFX_ASSERT(cfx_fe25519_eq(&a, &expected));
}

/* ============================================================ */
/* Fermat's little theorem tests */
/* ============================================================ */

static void test_fermat_little_theorem(void) {
    /* a^p = a (mod p), so a^(p-1) = 1 for a != 0 */
    /* Since inv uses a^(p-2), we verify a * a^(p-2) = 1 */
    fe25519_t a = {{ 0x7777777, 0x8888888, 0x9999999, 0xaaaaaaa, 0xbbbbbbb }};
    fe25519_t inv_a, product, one;

    cfx_fe25519_inv(&inv_a, &a);
    cfx_fe25519_mul(&product, &a, &inv_a);
    cfx_fe25519_1(&one);

    CFX_ASSERT(cfx_fe25519_eq(&product, &one));
}

static void test_fermat_various_values(void) {
    /* test several values - all limbs must be < 2^51 */
    fe25519_t values[] = {
        {{ 2, 0, 0, 0, 0 }},
        {{ 0x123456789abc, 0, 0, 0, 0 }},
        {{ 1, 1, 1, 1, 1 }},
        {{ 0x7edcba9876543, 0x123456789abcd, 0xaabbccdd, 0x11223344, 0x55667788 }},
    };
    fe25519_t one;
    cfx_fe25519_1(&one);

    for (size_t i = 0; i < sizeof(values)/sizeof(values[0]); i++) {
        fe25519_t inv, product;
        cfx_fe25519_inv(&inv, &values[i]);
        cfx_fe25519_mul(&product, &values[i], &inv);
        CFX_ASSERT(cfx_fe25519_eq(&product, &one));
    }
}

/* ============================================================ */
/* Near-prime edge cases */
/* ============================================================ */

static void test_p_minus_1(void) {
    /* p - 1 = -1 (mod p) */
    const uint64_t mask51 = (1ULL << 51) - 1;
    fe25519_t p_minus_1 = {{ mask51 - 19, mask51, mask51, mask51, mask51 }};
    fe25519_t one, neg_one;

    cfx_fe25519_1(&one);
    cfx_fe25519_neg(&neg_one, &one);

    CFX_ASSERT(cfx_fe25519_eq(&p_minus_1, &neg_one));
}

static void test_p_plus_1(void) {
    /* p + 1 = 1 (mod p) */
    const uint64_t mask51 = (1ULL << 51) - 1;
    fe25519_t p = {{ mask51 - 18, mask51, mask51, mask51, mask51 }};
    fe25519_t one, p_plus_1;

    cfx_fe25519_1(&one);
    cfx_fe25519_add(&p_plus_1, &p, &one);

    CFX_ASSERT(cfx_fe25519_eq(&p_plus_1, &one));
}

static void test_p_minus_19(void) {
    /* p - 19 + 19 = p = 0 (mod p) */
    const uint64_t mask51 = (1ULL << 51) - 1;
    fe25519_t p_minus_19 = {{ mask51 - 37, mask51, mask51, mask51, mask51 }};
    fe25519_t nineteen = {{ 19, 0, 0, 0, 0 }};
    fe25519_t sum;

    cfx_fe25519_add(&sum, &p_minus_19, &nineteen);
    CFX_ASSERT(cfx_fe25519_iszero(&sum));
}

/* ============================================================ */
/* More serialization tests */
/* ============================================================ */

static void test_bytes_max_value(void) {
    /* max value before reduction: all 0xFF except top bit */
    uint8_t max_bytes[32];
    memset(max_bytes, 0xFF, 32);
    max_bytes[31] = 0x7F;  /* clear top bit */

    fe25519_t f;
    uint8_t output[32];

    cfx_fe25519_frombytes(&f, max_bytes);
    cfx_fe25519_tobytes(output, &f);

    /* should roundtrip, though reduced */
    fe25519_t f2;
    cfx_fe25519_frombytes(&f2, output);
    CFX_ASSERT(cfx_fe25519_eq(&f, &f2));
}

static void test_bytes_alternating(void) {
    uint8_t alt_bytes[32];
    for (int i = 0; i < 32; i++) {
        alt_bytes[i] = (i & 1) ? 0xAA : 0x55;
    }
    alt_bytes[31] &= 0x7F;

    fe25519_t f;
    uint8_t output[32];

    cfx_fe25519_frombytes(&f, alt_bytes);
    cfx_fe25519_tobytes(output, &f);

    CFX_ASSERT(memcmp(alt_bytes, output, 32) == 0);
}

static void test_bytes_single_bits(void) {
    /* test each single bit position */
    for (int bit = 0; bit < 255; bit++) {
        uint8_t bytes[32] = {0};
        bytes[bit / 8] = 1 << (bit % 8);

        fe25519_t f;
        uint8_t output[32];

        cfx_fe25519_frombytes(&f, bytes);
        cfx_fe25519_tobytes(output, &f);

        CFX_ASSERT(memcmp(bytes, output, 32) == 0);
    }
}

/* ============================================================ */
/* Double negation and other identities */
/* ============================================================ */

static void test_double_neg(void) {
    fe25519_t a = {{ 0x123456, 0x789abc, 0xdef012, 0x345678, 0x9abcde }};
    fe25519_t neg_a, neg_neg_a;

    cfx_fe25519_neg(&neg_a, &a);
    cfx_fe25519_neg(&neg_neg_a, &neg_a);

    CFX_ASSERT(cfx_fe25519_eq(&a, &neg_neg_a));
}

static void test_add_neg_is_sub(void) {
    fe25519_t a = {{ 0x123456, 0x789abc, 0xdef012, 0, 0 }};
    fe25519_t b = {{ 0x111111, 0x222222, 0x333333, 0, 0 }};
    fe25519_t neg_b, a_plus_negb, a_minus_b;

    cfx_fe25519_neg(&neg_b, &b);
    cfx_fe25519_add(&a_plus_negb, &a, &neg_b);
    cfx_fe25519_sub(&a_minus_b, &a, &b);

    CFX_ASSERT(cfx_fe25519_eq(&a_plus_negb, &a_minus_b));
}

static void test_sub_neg_is_add(void) {
    fe25519_t a = {{ 0x123456, 0x789abc, 0xdef012, 0, 0 }};
    fe25519_t b = {{ 0x111111, 0x222222, 0x333333, 0, 0 }};
    fe25519_t neg_b, a_minus_negb, a_plus_b;

    cfx_fe25519_neg(&neg_b, &b);
    cfx_fe25519_sub(&a_minus_negb, &a, &neg_b);
    cfx_fe25519_add(&a_plus_b, &a, &b);

    CFX_ASSERT(cfx_fe25519_eq(&a_minus_negb, &a_plus_b));
}

/* ============================================================ */
/* Power identities */
/* ============================================================ */

static void test_sqr_is_mul_self(void) {
    /* verify sqr(a) == mul(a, a) for several values */
    fe25519_t values[] = {
        {{ 0, 0, 0, 0, 0 }},
        {{ 1, 0, 0, 0, 0 }},
        {{ 2, 0, 0, 0, 0 }},
        {{ 0x123456789, 0xabcdef012, 0x3456789ab, 0xcdef01234, 0x56789abcd }},
        {{ 0x7ffffffffffff, 0x7ffffffffffff, 0x7ffffffffffff, 0x7ffffffffffff, 0x7ffffffffffff }},
    };

    for (size_t i = 0; i < sizeof(values)/sizeof(values[0]); i++) {
        fe25519_t sqr_result, mul_result;
        cfx_fe25519_sqr(&sqr_result, &values[i]);
        cfx_fe25519_mul(&mul_result, &values[i], &values[i]);
        CFX_ASSERT(cfx_fe25519_eq(&sqr_result, &mul_result));
    }
}

static void test_sqr4_is_sqr2_sqr2(void) {
    fe25519_t a = {{ 0x12345, 0x67890, 0xabcde, 0xf0123, 0x45678 }};
    fe25519_t sqr4, sqr2;

    cfx_fe25519_sqr_n(&sqr4, &a, 4);

    cfx_fe25519_sqr_n(&sqr2, &a, 2);
    cfx_fe25519_sqr_n(&sqr2, &sqr2, 2);

    CFX_ASSERT(cfx_fe25519_eq(&sqr4, &sqr2));
}

/* ============================================================ */
/* cswap/cmov stress */
/* ============================================================ */

static void test_cswap_multiple(void) {
    fe25519_t a = {{ 1, 2, 3, 4, 5 }};
    fe25519_t b = {{ 10, 20, 30, 40, 50 }};
    fe25519_t a_orig, b_orig;
    cfx_fe25519_copy(&a_orig, &a);
    cfx_fe25519_copy(&b_orig, &b);

    /* swap twice = no change */
    cfx_fe25519_cswap(&a, &b, 1);
    cfx_fe25519_cswap(&a, &b, 1);

    CFX_ASSERT(fe_limbs_eq(&a, &a_orig));
    CFX_ASSERT(fe_limbs_eq(&b, &b_orig));
}

static void test_cmov_chain(void) {
    fe25519_t h = {{ 1, 2, 3, 4, 5 }};
    fe25519_t f1 = {{ 10, 20, 30, 40, 50 }};
    fe25519_t f2 = {{ 100, 200, 300, 400, 500 }};

    cfx_fe25519_cmov(&h, &f1, 0);  /* no change */
    cfx_fe25519_cmov(&h, &f2, 1);  /* h = f2 */

    CFX_ASSERT(fe_limbs_eq(&h, &f2));
}

/* ============================================================ */
/* mul121666 extensive tests */
/* ============================================================ */

static void test_mul121666_large(void) {
    fe25519_t a = {{ 0x7ffffffffffff, 0x7ffffffffffff, 0x7ffffffffffff, 0x7ffffffffffff, 0x7ffffffffffff }};
    fe25519_t k = {{ 121666, 0, 0, 0, 0 }};
    fe25519_t via_mul121666, via_mul;

    cfx_fe25519_mul121666(&via_mul121666, &a);
    cfx_fe25519_mul(&via_mul, &a, &k);

    CFX_ASSERT(cfx_fe25519_eq(&via_mul121666, &via_mul));
}

static void test_mul121666_chain(void) {
    fe25519_t a = {{ 1, 0, 0, 0, 0 }};

    for (int i = 0; i < 10; i++) {
        cfx_fe25519_mul121666(&a, &a);
    }

    /* verify result matches 121666^10 mod p */
    fe25519_t k = {{ 121666, 0, 0, 0, 0 }};
    fe25519_t expected;
    cfx_fe25519_1(&expected);
    for (int i = 0; i < 10; i++) {
        cfx_fe25519_mul(&expected, &expected, &k);
    }

    CFX_ASSERT(cfx_fe25519_eq(&a, &expected));
}

/* ============================================================ */
/* Stress tests */
/* ============================================================ */

static void test_many_operations(void) {
    fe25519_t a, b, c, d;
    cfx_fe25519_1(&a);
    cfx_fe25519_copy(&b, &a);

    /* do many operations */
    for (int i = 0; i < 100; i++) {
        cfx_fe25519_add(&c, &a, &b);
        cfx_fe25519_mul(&d, &a, &b);
        cfx_fe25519_sqr(&a, &c);
        cfx_fe25519_add(&b, &d, &a);
    }

    /* just verify we got some non-zero result */
    CFX_ASSERT(cfx_fe25519_iszero(&a) == 0 || cfx_fe25519_iszero(&b) == 0);
}

static void test_inv_stress(void) {
    fe25519_t a = {{ 0x12345678, 0x9abcdef0, 0x12345678, 0x9abcdef0, 0x12345678 }};
    fe25519_t inv, product, one;

    for (int i = 0; i < 10; i++) {
        cfx_fe25519_inv(&inv, &a);
        cfx_fe25519_mul(&product, &a, &inv);
        cfx_fe25519_1(&one);
        CFX_ASSERT(cfx_fe25519_eq(&product, &one));

        /* next iteration: use inv as new a */
        cfx_fe25519_add(&a, &a, &inv);
        cfx_fe25519_sqr(&a, &a);
    }
}

int main(void) {
    /* basic */
    CFX_TEST(test_zero);
    CFX_TEST(test_one);
    CFX_TEST(test_copy);
    CFX_TEST(test_constants);

    /* add/sub */
    CFX_TEST(test_add_basic);
    CFX_TEST(test_add_zero);
    CFX_TEST(test_sub_basic);
    CFX_TEST(test_sub_self);
    CFX_TEST(test_neg_basic);
    CFX_TEST(test_neg_zero);

    /* mul */
    CFX_TEST(test_mul_by_zero);
    CFX_TEST(test_mul_by_one);
    CFX_TEST(test_mul_small);
    CFX_TEST(test_mul_commutative);
    CFX_TEST(test_mul_associative);
    CFX_TEST(test_mul_distributive);

    /* sqr */
    CFX_TEST(test_sqr_zero);
    CFX_TEST(test_sqr_one);
    CFX_TEST(test_sqr_vs_mul);
    CFX_TEST(test_sqr_n);

    /* mul121666 */
    CFX_TEST(test_mul121666_zero);
    CFX_TEST(test_mul121666_one);
    CFX_TEST(test_mul121666_vs_mul);

    /* inv */
    CFX_TEST(test_inv_one);
    CFX_TEST(test_inv_mul_identity);
    CFX_TEST(test_inv_inv);

    /* serialization */
    CFX_TEST(test_bytes_roundtrip);
    CFX_TEST(test_bytes_zero);
    CFX_TEST(test_bytes_one);
    CFX_TEST(test_bytes_clears_top_bit);

    /* constant-time ops */
    CFX_TEST(test_cswap_zero);
    CFX_TEST(test_cswap_one);
    CFX_TEST(test_cmov_zero);
    CFX_TEST(test_cmov_one);

    /* comparison */
    CFX_TEST(test_eq_same);
    CFX_TEST(test_eq_different);
    CFX_TEST(test_eq_reduced_vs_unreduced);
    CFX_TEST(test_isnegative);

    /* reduction */
    CFX_TEST(test_reduce_small);
    CFX_TEST(test_reduce_large_limb);
    CFX_TEST(test_carry_propagation);

    /* RFC 7748 */
    CFX_TEST(test_basepoint);
    CFX_TEST(test_add_sub_inverse);
    CFX_TEST(test_mul_inv_identity);

    /* edge cases */
    CFX_TEST(test_max_limb_values);
    CFX_TEST(test_sqr_large);
    CFX_TEST(test_p_is_zero);
    CFX_TEST(test_2p_is_zero);

    /* pow22523 */
    CFX_TEST(test_pow22523_basic);

    /* aliasing */
    CFX_TEST(test_add_aliased);
    CFX_TEST(test_sub_aliased);
    CFX_TEST(test_mul_aliased);
    CFX_TEST(test_mul_self_aliased);
    CFX_TEST(test_sqr_aliased);
    CFX_TEST(test_neg_aliased);
    CFX_TEST(test_inv_aliased);

    /* fermat */
    CFX_TEST(test_fermat_little_theorem);
    CFX_TEST(test_fermat_various_values);

    /* near-prime */
    CFX_TEST(test_p_minus_1);
    CFX_TEST(test_p_plus_1);
    CFX_TEST(test_p_minus_19);

    /* more serialization */
    CFX_TEST(test_bytes_max_value);
    CFX_TEST(test_bytes_alternating);
    CFX_TEST(test_bytes_single_bits);

    /* identities */
    CFX_TEST(test_double_neg);
    CFX_TEST(test_add_neg_is_sub);
    CFX_TEST(test_sub_neg_is_add);

    /* power identities */
    CFX_TEST(test_sqr_is_mul_self);
    CFX_TEST(test_sqr4_is_sqr2_sqr2);

    /* cswap/cmov stress */
    CFX_TEST(test_cswap_multiple);
    CFX_TEST(test_cmov_chain);

    /* mul121666 extensive */
    CFX_TEST(test_mul121666_large);
    CFX_TEST(test_mul121666_chain);

    /* stress */
    CFX_TEST(test_many_operations);
    CFX_TEST(test_inv_stress);

    return 0;
}
