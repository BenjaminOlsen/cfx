/* test_ge25519.c - Tests for Ed25519 group element operations */

#include "cfx/ge25519.h"
#include "cfx/macros.h"
#include <stdio.h>
#include <string.h>

/* test identity is on curve and has correct properties */
static void test_identity(void) {
    ge25519_t id;
    cfx_ge25519_0(&id);
    CFX_ASSERT(cfx_ge25519_is_identity(&id));
    CFX_ASSERT(cfx_ge25519_is_on_curve(&id));
}

/* test basepoint is on curve */
static void test_basepoint_on_curve(void) {
    CFX_ASSERT(cfx_ge25519_is_on_curve(&cfx_ge25519_base));
}

/* test basepoint is not identity */
static void test_basepoint_not_identity(void) {
    CFX_ASSERT(!cfx_ge25519_is_identity(&cfx_ge25519_base));
}

/* test pack/unpack roundtrip */
static void test_pack_unpack_basepoint(void) {
    uint8_t packed[32];
    ge25519_t unpacked;

    cfx_ge25519_pack(packed, &cfx_ge25519_base);
    CFX_ASSERT(cfx_ge25519_unpack(&unpacked, packed) == 0);
    CFX_ASSERT(cfx_ge25519_is_on_curve(&unpacked));

    uint8_t repacked[32];
    cfx_ge25519_pack(repacked, &unpacked);
    CFX_ASSERT(memcmp(packed, repacked, 32) == 0);
}

/* test 1*B = B */
static void test_scalarmult_one(void) {
    uint8_t scalar[32] = {1};
    ge25519_t result;

    cfx_ge25519_scalarmult_base(&result, scalar);

    uint8_t base_packed[32], result_packed[32];
    cfx_ge25519_pack(base_packed, &cfx_ge25519_base);
    cfx_ge25519_pack(result_packed, &result);
    CFX_ASSERT(memcmp(base_packed, result_packed, 32) == 0);
}

/* test 0*B = identity */
static void test_scalarmult_zero(void) {
    uint8_t scalar[32] = {0};
    ge25519_t result;

    cfx_ge25519_scalarmult_base(&result, scalar);
    CFX_ASSERT(cfx_ge25519_is_identity(&result));
}

/* test 2*B via scalarmult vs double */
static void test_scalarmult_two_vs_double(void) {
    uint8_t scalar[32] = {2};
    ge25519_t result_scalar, result_double;

    cfx_ge25519_scalarmult_base(&result_scalar, scalar);
    cfx_ge25519_double(&result_double, &cfx_ge25519_base);

    uint8_t s_packed[32], d_packed[32];
    cfx_ge25519_pack(s_packed, &result_scalar);
    cfx_ge25519_pack(d_packed, &result_double);
    CFX_ASSERT(memcmp(s_packed, d_packed, 32) == 0);
}

/* test B + B = 2*B via add_cached */
static void test_add_equals_double(void) {
    ge25519_t sum, doubled;
    ge25519_cached_t b_cached;

    cfx_ge25519_to_cached(&b_cached, &cfx_ge25519_base);
    cfx_ge25519_add_cached(&sum, &cfx_ge25519_base, &b_cached);
    cfx_ge25519_double(&doubled, &cfx_ge25519_base);

    uint8_t sum_packed[32], doubled_packed[32];
    cfx_ge25519_pack(sum_packed, &sum);
    cfx_ge25519_pack(doubled_packed, &doubled);
    CFX_ASSERT(memcmp(sum_packed, doubled_packed, 32) == 0);
}

/* test P - P = identity */
static void test_sub_self_is_identity(void) {
    ge25519_t result;
    ge25519_cached_t b_cached;

    cfx_ge25519_to_cached(&b_cached, &cfx_ge25519_base);
    cfx_ge25519_sub_cached(&result, &cfx_ge25519_base, &b_cached);
    CFX_ASSERT(cfx_ge25519_is_identity(&result));
}

/* test negation: P + (-P) = identity */
static void test_negation(void) {
    ge25519_t neg_b, sum;
    ge25519_cached_t neg_cached;

    cfx_ge25519_neg(&neg_b, &cfx_ge25519_base);
    cfx_ge25519_to_cached(&neg_cached, &neg_b);
    cfx_ge25519_add_cached(&sum, &cfx_ge25519_base, &neg_cached);
    CFX_ASSERT(cfx_ge25519_is_identity(&sum));
}

/* test associativity: (B + B) + B = B + (B + B) */
static void test_add_associative(void) {
    ge25519_t two_b, three_b_left, three_b_right;
    ge25519_cached_t b_cached, two_b_cached;

    cfx_ge25519_to_cached(&b_cached, &cfx_ge25519_base);
    cfx_ge25519_add_cached(&two_b, &cfx_ge25519_base, &b_cached);
    cfx_ge25519_to_cached(&two_b_cached, &two_b);

    cfx_ge25519_add_cached(&three_b_left, &two_b, &b_cached);
    cfx_ge25519_add_cached(&three_b_right, &cfx_ge25519_base, &two_b_cached);

    uint8_t left[32], right[32];
    cfx_ge25519_pack(left, &three_b_left);
    cfx_ge25519_pack(right, &three_b_right);
    CFX_ASSERT(memcmp(left, right, 32) == 0);
}

/* test commutativity: B + 2B = 2B + B */
static void test_add_commutative(void) {
    ge25519_t two_b, sum1, sum2;
    ge25519_cached_t b_cached, two_b_cached;

    cfx_ge25519_double(&two_b, &cfx_ge25519_base);
    cfx_ge25519_to_cached(&b_cached, &cfx_ge25519_base);
    cfx_ge25519_to_cached(&two_b_cached, &two_b);

    cfx_ge25519_add_cached(&sum1, &cfx_ge25519_base, &two_b_cached);
    cfx_ge25519_add_cached(&sum2, &two_b, &b_cached);

    uint8_t s1[32], s2[32];
    cfx_ge25519_pack(s1, &sum1);
    cfx_ge25519_pack(s2, &sum2);
    CFX_ASSERT(memcmp(s1, s2, 32) == 0);
}

/* test 3*B via scalarmult vs explicit additions */
static void test_scalarmult_three(void) {
    uint8_t scalar[32] = {3};
    ge25519_t result_scalar;
    cfx_ge25519_scalarmult_base(&result_scalar, scalar);

    ge25519_t two_b, three_b;
    ge25519_cached_t b_cached;
    cfx_ge25519_to_cached(&b_cached, &cfx_ge25519_base);
    cfx_ge25519_add_cached(&two_b, &cfx_ge25519_base, &b_cached);
    cfx_ge25519_add_cached(&three_b, &two_b, &b_cached);

    uint8_t s[32], e[32];
    cfx_ge25519_pack(s, &result_scalar);
    cfx_ge25519_pack(e, &three_b);
    CFX_ASSERT(memcmp(s, e, 32) == 0);
}

/* test known scalar: L (group order) * B = identity */
static void test_group_order(void) {
    uint8_t L[32] = {
        0xed, 0xd3, 0xf5, 0x5c, 0x1a, 0x63, 0x12, 0x58,
        0xd6, 0x9c, 0xf7, 0xa2, 0xde, 0xf9, 0xde, 0x14,
        0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x10
    };

    ge25519_t result;
    cfx_ge25519_scalarmult_base(&result, L);
    CFX_ASSERT(cfx_ge25519_is_identity(&result));
}

/* test known vector from RFC 8032 */
static void test_rfc8032_basepoint_y(void) {
    static const uint8_t expected[32] = {
        0x58, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66,
        0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66,
        0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66,
        0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66
    };

    uint8_t packed[32];
    cfx_ge25519_pack(packed, &cfx_ge25519_base);
    CFX_ASSERT(memcmp(packed, expected, 32) == 0);
}

/* test doubling produces point on curve */
static void test_double_on_curve(void) {
    ge25519_t doubled;
    cfx_ge25519_double(&doubled, &cfx_ge25519_base);
    CFX_ASSERT(cfx_ge25519_is_on_curve(&doubled));
}

/* test repeated doubling: 4B, 8B, 16B all on curve */
static void test_repeated_doubling_on_curve(void) {
    ge25519_t p;
    cfx_ge25519_copy(&p, &cfx_ge25519_base);

    for (int i = 0; i < 4; i++) {
        cfx_ge25519_double(&p, &p);
        CFX_ASSERT(cfx_ge25519_is_on_curve(&p));
    }
}

/* test 4*B via scalarmult vs double(double(B)) */
static void test_scalarmult_four(void) {
    uint8_t scalar[32] = {4};
    ge25519_t result_scalar, result_double;

    cfx_ge25519_scalarmult_base(&result_scalar, scalar);
    cfx_ge25519_double(&result_double, &cfx_ge25519_base);
    cfx_ge25519_double(&result_double, &result_double);

    uint8_t s[32], d[32];
    cfx_ge25519_pack(s, &result_scalar);
    cfx_ge25519_pack(d, &result_double);
    CFX_ASSERT(memcmp(s, d, 32) == 0);
}

/* test scalarmult with larger value */
static void test_scalarmult_255(void) {
    uint8_t scalar[32] = {0xff};
    ge25519_t result;

    cfx_ge25519_scalarmult_base(&result, scalar);
    CFX_ASSERT(cfx_ge25519_is_on_curve(&result));
    CFX_ASSERT(!cfx_ge25519_is_identity(&result));
}

/* test copy operation */
static void test_copy(void) {
    ge25519_t copy;
    cfx_ge25519_copy(&copy, &cfx_ge25519_base);

    uint8_t orig[32], cpy[32];
    cfx_ge25519_pack(orig, &cfx_ge25519_base);
    cfx_ge25519_pack(cpy, &copy);
    CFX_ASSERT(memcmp(orig, cpy, 32) == 0);
}

/* test cneg (conditional negate) */
static void test_cneg(void) {
    ge25519_t p, p_orig;
    cfx_ge25519_copy(&p, &cfx_ge25519_base);
    cfx_ge25519_copy(&p_orig, &cfx_ge25519_base);

    cfx_ge25519_cneg(&p, 0);
    uint8_t before[32], after[32];
    cfx_ge25519_pack(before, &p_orig);
    cfx_ge25519_pack(after, &p);
    CFX_ASSERT(memcmp(before, after, 32) == 0);

    cfx_ge25519_copy(&p, &cfx_ge25519_base);
    cfx_ge25519_cneg(&p, 1);
    ge25519_t neg;
    cfx_ge25519_neg(&neg, &cfx_ge25519_base);
    uint8_t cneg_result[32], neg_result[32];
    cfx_ge25519_pack(cneg_result, &p);
    cfx_ge25519_pack(neg_result, &neg);
    CFX_ASSERT(memcmp(cneg_result, neg_result, 32) == 0);
}

/* test cmov (conditional move) */
static void test_cmov(void) {
    ge25519_t id, p;
    cfx_ge25519_0(&id);
    cfx_ge25519_copy(&p, &id);

    cfx_ge25519_cmov(&p, &cfx_ge25519_base, 0);
    CFX_ASSERT(cfx_ge25519_is_identity(&p));

    cfx_ge25519_cmov(&p, &cfx_ge25519_base, 1);
    uint8_t p_packed[32], b_packed[32];
    cfx_ge25519_pack(p_packed, &p);
    cfx_ge25519_pack(b_packed, &cfx_ge25519_base);
    CFX_ASSERT(memcmp(p_packed, b_packed, 32) == 0);
}

/* test unpack with invalid point (not on curve) */
static void test_unpack_invalid(void) {
    uint8_t invalid[32] = {2};
    ge25519_t p;
    CFX_ASSERT(cfx_ge25519_unpack(&p, invalid) != 0);
}

/* test unpack identity encoding: (0, 1) */
static void test_unpack_identity(void) {
    uint8_t id_bytes[32] = {1};
    ge25519_t p;
    CFX_ASSERT(cfx_ge25519_unpack(&p, id_bytes) == 0);
    CFX_ASSERT(cfx_ge25519_is_identity(&p));
}

/* test scalar multiplication distributivity: (a+b)*B = a*B + b*B */
static void test_scalarmult_distributive(void) {
    uint8_t a[32] = {5};
    uint8_t b[32] = {7};
    uint8_t a_plus_b[32] = {12};

    ge25519_t aB, bB, sum, direct;
    ge25519_cached_t bB_cached;

    cfx_ge25519_scalarmult_base(&aB, a);
    cfx_ge25519_scalarmult_base(&bB, b);
    cfx_ge25519_to_cached(&bB_cached, &bB);
    cfx_ge25519_add_cached(&sum, &aB, &bB_cached);
    cfx_ge25519_scalarmult_base(&direct, a_plus_b);

    uint8_t sum_packed[32], direct_packed[32];
    cfx_ge25519_pack(sum_packed, &sum);
    cfx_ge25519_pack(direct_packed, &direct);
    CFX_ASSERT(memcmp(sum_packed, direct_packed, 32) == 0);
}

/* test larger scalar */
static void test_scalarmult_0x12345(void) {
    uint8_t scalar[32] = {0x45, 0x23, 0x01};
    ge25519_t result;

    cfx_ge25519_scalarmult_base(&result, scalar);
    CFX_ASSERT(cfx_ge25519_is_on_curve(&result));
    CFX_ASSERT(!cfx_ge25519_is_identity(&result));
}

/* test that doubled identity is still identity */
static void test_double_identity(void) {
    ge25519_t id, doubled;
    cfx_ge25519_0(&id);
    cfx_ge25519_double(&doubled, &id);
    CFX_ASSERT(cfx_ge25519_is_identity(&doubled));
}

/* test identity + B = B */
static void test_identity_add(void) {
    ge25519_t id, sum;
    ge25519_cached_t b_cached;

    cfx_ge25519_0(&id);
    cfx_ge25519_to_cached(&b_cached, &cfx_ge25519_base);
    cfx_ge25519_add_cached(&sum, &id, &b_cached);

    uint8_t sum_packed[32], b_packed[32];
    cfx_ge25519_pack(sum_packed, &sum);
    cfx_ge25519_pack(b_packed, &cfx_ge25519_base);
    CFX_ASSERT(memcmp(sum_packed, b_packed, 32) == 0);
}

/* test B + identity = B */
static void test_add_identity_right(void) {
    ge25519_t id, sum;
    ge25519_cached_t id_cached;

    cfx_ge25519_0(&id);
    cfx_ge25519_to_cached(&id_cached, &id);
    cfx_ge25519_add_cached(&sum, &cfx_ge25519_base, &id_cached);

    uint8_t sum_packed[32], b_packed[32];
    cfx_ge25519_pack(sum_packed, &sum);
    cfx_ge25519_pack(b_packed, &cfx_ge25519_base);
    CFX_ASSERT(memcmp(sum_packed, b_packed, 32) == 0);
}

/* test 5*B via scalarmult vs additions */
static void test_scalarmult_five(void) {
    uint8_t scalar[32] = {5};
    ge25519_t result_scalar;
    cfx_ge25519_scalarmult_base(&result_scalar, scalar);

    ge25519_t acc;
    ge25519_cached_t b_cached;
    cfx_ge25519_copy(&acc, &cfx_ge25519_base);
    cfx_ge25519_to_cached(&b_cached, &cfx_ge25519_base);
    for (int i = 0; i < 4; i++) {
        cfx_ge25519_add_cached(&acc, &acc, &b_cached);
    }

    uint8_t s[32], e[32];
    cfx_ge25519_pack(s, &result_scalar);
    cfx_ge25519_pack(e, &acc);
    CFX_ASSERT(memcmp(s, e, 32) == 0);
}

/* test 8*B (cofactor) */
static void test_scalarmult_eight(void) {
    uint8_t scalar[32] = {8};
    ge25519_t result;
    cfx_ge25519_scalarmult_base(&result, scalar);
    CFX_ASSERT(cfx_ge25519_is_on_curve(&result));
    CFX_ASSERT(!cfx_ge25519_is_identity(&result));
}

/* test 16*B */
static void test_scalarmult_sixteen(void) {
    uint8_t scalar[32] = {16};
    ge25519_t result_scalar, result_double;

    cfx_ge25519_scalarmult_base(&result_scalar, scalar);
    cfx_ge25519_copy(&result_double, &cfx_ge25519_base);
    for (int i = 0; i < 4; i++) {
        cfx_ge25519_double(&result_double, &result_double);
    }

    uint8_t s[32], d[32];
    cfx_ge25519_pack(s, &result_scalar);
    cfx_ge25519_pack(d, &result_double);
    CFX_ASSERT(memcmp(s, d, 32) == 0);
}

/* test 256*B */
static void test_scalarmult_256(void) {
    uint8_t scalar[32] = {0, 1};
    ge25519_t result;
    cfx_ge25519_scalarmult_base(&result, scalar);
    CFX_ASSERT(cfx_ge25519_is_on_curve(&result));
    CFX_ASSERT(!cfx_ge25519_is_identity(&result));
}

/* test large scalar 0xffffffff */
static void test_scalarmult_large(void) {
    uint8_t scalar[32] = {0xff, 0xff, 0xff, 0xff};
    ge25519_t result;
    cfx_ge25519_scalarmult_base(&result, scalar);
    CFX_ASSERT(cfx_ge25519_is_on_curve(&result));
    CFX_ASSERT(!cfx_ge25519_is_identity(&result));
}

/* test (L-1)*B is not identity but (L-1)*B + B = identity */
static void test_L_minus_1(void) {
    uint8_t L_minus_1[32] = {
        0xec, 0xd3, 0xf5, 0x5c, 0x1a, 0x63, 0x12, 0x58,
        0xd6, 0x9c, 0xf7, 0xa2, 0xde, 0xf9, 0xde, 0x14,
        0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x10
    };

    ge25519_t result;
    cfx_ge25519_scalarmult_base(&result, L_minus_1);
    CFX_ASSERT(!cfx_ge25519_is_identity(&result));
    CFX_ASSERT(cfx_ge25519_is_on_curve(&result));

    ge25519_cached_t b_cached;
    cfx_ge25519_to_cached(&b_cached, &cfx_ge25519_base);
    ge25519_t sum;
    cfx_ge25519_add_cached(&sum, &result, &b_cached);
    CFX_ASSERT(cfx_ge25519_is_identity(&sum));
}

/* test double negation: -(-P) = P */
static void test_double_neg(void) {
    ge25519_t neg1, neg2;
    cfx_ge25519_neg(&neg1, &cfx_ge25519_base);
    cfx_ge25519_neg(&neg2, &neg1);

    uint8_t p1[32], p2[32];
    cfx_ge25519_pack(p1, &cfx_ge25519_base);
    cfx_ge25519_pack(p2, &neg2);
    CFX_ASSERT(memcmp(p1, p2, 32) == 0);
}

/* test negation of identity is identity */
static void test_neg_identity(void) {
    ge25519_t id, neg_id;
    cfx_ge25519_0(&id);
    cfx_ge25519_neg(&neg_id, &id);
    CFX_ASSERT(cfx_ge25519_is_identity(&neg_id));
}

/* test subtraction: B - B = identity via sub_cached */
static void test_sub_basic(void) {
    ge25519_t result;
    ge25519_cached_t b_cached;
    cfx_ge25519_to_cached(&b_cached, &cfx_ge25519_base);
    cfx_ge25519_sub_cached(&result, &cfx_ge25519_base, &b_cached);
    CFX_ASSERT(cfx_ge25519_is_identity(&result));
}

/* test 2B - B = B */
static void test_sub_2B_minus_B(void) {
    ge25519_t two_b, result;
    ge25519_cached_t b_cached;

    cfx_ge25519_double(&two_b, &cfx_ge25519_base);
    cfx_ge25519_to_cached(&b_cached, &cfx_ge25519_base);
    cfx_ge25519_sub_cached(&result, &two_b, &b_cached);

    uint8_t r[32], b[32];
    cfx_ge25519_pack(r, &result);
    cfx_ge25519_pack(b, &cfx_ge25519_base);
    CFX_ASSERT(memcmp(r, b, 32) == 0);
}

/* test 3B - 2B = B */
static void test_sub_3B_minus_2B(void) {
    uint8_t s2[32] = {2}, s3[32] = {3};
    ge25519_t two_b, three_b, result;
    ge25519_cached_t two_b_cached;

    cfx_ge25519_scalarmult_base(&two_b, s2);
    cfx_ge25519_scalarmult_base(&three_b, s3);
    cfx_ge25519_to_cached(&two_b_cached, &two_b);
    cfx_ge25519_sub_cached(&result, &three_b, &two_b_cached);

    uint8_t r[32], b[32];
    cfx_ge25519_pack(r, &result);
    cfx_ge25519_pack(b, &cfx_ge25519_base);
    CFX_ASSERT(memcmp(r, b, 32) == 0);
}

/* test pack/unpack roundtrip for 2B */
static void test_pack_unpack_2B(void) {
    ge25519_t two_b, unpacked;
    uint8_t packed[32];

    cfx_ge25519_double(&two_b, &cfx_ge25519_base);
    cfx_ge25519_pack(packed, &two_b);
    CFX_ASSERT(cfx_ge25519_unpack(&unpacked, packed) == 0);

    uint8_t repacked[32];
    cfx_ge25519_pack(repacked, &unpacked);
    CFX_ASSERT(memcmp(packed, repacked, 32) == 0);
}

/* test pack/unpack roundtrip for various scalars */
static void test_pack_unpack_various(void) {
    uint8_t scalars[5][32] = {
        {7}, {42}, {100}, {0xff}, {0x12, 0x34}
    };

    for (int i = 0; i < 5; i++) {
        ge25519_t p, unpacked;
        uint8_t packed[32];

        cfx_ge25519_scalarmult_base(&p, scalars[i]);
        cfx_ge25519_pack(packed, &p);
        CFX_ASSERT(cfx_ge25519_unpack(&unpacked, packed) == 0);

        uint8_t repacked[32];
        cfx_ge25519_pack(repacked, &unpacked);
        CFX_ASSERT(memcmp(packed, repacked, 32) == 0);
    }
}

/* test unpack rejects random bytes that aren't on curve */
static void test_unpack_random_invalid(void) {
    uint8_t invalid[32] = {
        0x12, 0x34, 0x56, 0x78, 0x9a, 0xbc, 0xde, 0xf0,
        0x11, 0x22, 0x33, 0x44, 0x55, 0x66, 0x77, 0x88,
        0xaa, 0xbb, 0xcc, 0xdd, 0xee, 0xff, 0x00, 0x11,
        0x22, 0x33, 0x44, 0x55, 0x66, 0x77, 0x88, 0x00
    };

    ge25519_t p;
    CFX_ASSERT(cfx_ge25519_unpack(&p, invalid) != 0);
}

/* test multiple additions produce points on curve */
static void test_many_additions_on_curve(void) {
    ge25519_t acc;
    ge25519_cached_t b_cached;

    cfx_ge25519_copy(&acc, &cfx_ge25519_base);
    cfx_ge25519_to_cached(&b_cached, &cfx_ge25519_base);

    for (int i = 0; i < 20; i++) {
        cfx_ge25519_add_cached(&acc, &acc, &b_cached);
        CFX_ASSERT(cfx_ge25519_is_on_curve(&acc));
    }
}

/* test many doublings produce points on curve */
static void test_many_doublings_on_curve(void) {
    ge25519_t p;
    cfx_ge25519_copy(&p, &cfx_ge25519_base);

    for (int i = 0; i < 20; i++) {
        cfx_ge25519_double(&p, &p);
        CFX_ASSERT(cfx_ge25519_is_on_curve(&p));
    }
}

/* test scalarmult with random-ish scalars all produce points on curve */
static void test_scalarmult_on_curve_stress(void) {
    uint8_t scalars[10][32] = {
        {0x17}, {0x23, 0x45}, {0xab, 0xcd, 0xef},
        {0x11, 0x22, 0x33, 0x44}, {0xff, 0xee, 0xdd, 0xcc, 0xbb},
        {0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08},
        {0x99, 0x88, 0x77, 0x66, 0x55, 0x44, 0x33, 0x22, 0x11},
        {0xde, 0xad, 0xbe, 0xef, 0xca, 0xfe, 0xba, 0xbe},
        {0x12, 0x34, 0x56, 0x78, 0x9a, 0xbc, 0xde, 0xf0, 0x12, 0x34},
        {0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x0f}
    };

    for (int i = 0; i < 10; i++) {
        ge25519_t result;
        cfx_ge25519_scalarmult_base(&result, scalars[i]);
        CFX_ASSERT(cfx_ge25519_is_on_curve(&result));
    }
}

/* test (a*b)*B = a*(b*B) for various a, b */
static void test_scalarmult_associative(void) {
    uint8_t s3[32] = {3}, s5[32] = {5}, s15[32] = {15};
    ge25519_t five_b, three_times_five_b, fifteen_b;

    cfx_ge25519_scalarmult_base(&five_b, s5);
    cfx_ge25519_scalarmult(&three_times_five_b, s3, &five_b);
    cfx_ge25519_scalarmult_base(&fifteen_b, s15);

    uint8_t p1[32], p2[32];
    cfx_ge25519_pack(p1, &three_times_five_b);
    cfx_ge25519_pack(p2, &fifteen_b);
    CFX_ASSERT(memcmp(p1, p2, 32) == 0);
}

/* test scalarmult with non-basepoint */
static void test_scalarmult_non_base(void) {
    uint8_t s2[32] = {2}, s3[32] = {3}, s6[32] = {6};
    ge25519_t two_b, result1, result2;

    cfx_ge25519_scalarmult_base(&two_b, s2);
    cfx_ge25519_scalarmult(&result1, s3, &two_b);  /* 3 * (2B) */
    cfx_ge25519_scalarmult_base(&result2, s6);     /* 6 * B */

    uint8_t p1[32], p2[32];
    cfx_ge25519_pack(p1, &result1);
    cfx_ge25519_pack(p2, &result2);
    CFX_ASSERT(memcmp(p1, p2, 32) == 0);
}

/* test cneg multiple times */
static void test_cneg_multiple(void) {
    ge25519_t p;
    cfx_ge25519_copy(&p, &cfx_ge25519_base);

    cfx_ge25519_cneg(&p, 1);
    cfx_ge25519_cneg(&p, 1);

    uint8_t p1[32], p2[32];
    cfx_ge25519_pack(p1, &p);
    cfx_ge25519_pack(p2, &cfx_ge25519_base);
    CFX_ASSERT(memcmp(p1, p2, 32) == 0);
}

/* test cmov chain */
static void test_cmov_chain(void) {
    ge25519_t p, two_b;
    uint8_t s2[32] = {2};

    cfx_ge25519_0(&p);
    cfx_ge25519_scalarmult_base(&two_b, s2);

    cfx_ge25519_cmov(&p, &cfx_ge25519_base, 1);  /* p = B */
    cfx_ge25519_cmov(&p, &two_b, 0);         /* p stays B */

    uint8_t p1[32], p2[32];
    cfx_ge25519_pack(p1, &p);
    cfx_ge25519_pack(p2, &cfx_ge25519_base);
    CFX_ASSERT(memcmp(p1, p2, 32) == 0);

    cfx_ge25519_cmov(&p, &two_b, 1);         /* p = 2B */
    cfx_ge25519_pack(p1, &p);
    cfx_ge25519_pack(p2, &two_b);
    CFX_ASSERT(memcmp(p1, p2, 32) == 0);
}

/* test to_cached preserves point info */
static void test_to_cached_preserves(void) {
    ge25519_t p, result;
    ge25519_cached_t cached;

    cfx_ge25519_0(&p);
    cfx_ge25519_to_cached(&cached, &cfx_ge25519_base);
    cfx_ge25519_add_cached(&result, &p, &cached);

    uint8_t p1[32], p2[32];
    cfx_ge25519_pack(p1, &result);
    cfx_ge25519_pack(p2, &cfx_ge25519_base);
    CFX_ASSERT(memcmp(p1, p2, 32) == 0);
}

/* test known 2B encoding */
static void test_known_2B(void) {
    ge25519_t two_b;
    cfx_ge25519_double(&two_b, &cfx_ge25519_base);

    uint8_t packed[32];
    cfx_ge25519_pack(packed, &two_b);

    /* 2B should be a specific known value (from reference impl) */
    static const uint8_t expected_2B[32] = {
        0xc9, 0xa3, 0xf8, 0x6a, 0xae, 0x46, 0x5f, 0x0e,
        0x56, 0x51, 0x38, 0x64, 0x51, 0x0f, 0x39, 0x97,
        0x56, 0x1f, 0xa2, 0xc9, 0xe8, 0x5e, 0xa2, 0x1d,
        0xc2, 0x29, 0x23, 0x09, 0xf3, 0xcd, 0x60, 0x22
    };
    CFX_ASSERT(memcmp(packed, expected_2B, 32) == 0);
}

/* test 8B (cofactor) is on curve, not identity, and distinct from B */
static void test_eight_B_properties(void) {
    uint8_t s8[32] = {8};
    ge25519_t eight_b;
    cfx_ge25519_scalarmult_base(&eight_b, s8);

    CFX_ASSERT(cfx_ge25519_is_on_curve(&eight_b));
    CFX_ASSERT(!cfx_ge25519_is_identity(&eight_b));

    /* 8B should be different from B */
    uint8_t b_packed[32], eight_packed[32];
    cfx_ge25519_pack(b_packed, &cfx_ge25519_base);
    cfx_ge25519_pack(eight_packed, &eight_b);
    CFX_ASSERT(memcmp(b_packed, eight_packed, 32) != 0);

    /* 8B should pack/unpack correctly */
    ge25519_t unpacked;
    CFX_ASSERT(cfx_ge25519_unpack(&unpacked, eight_packed) == 0);

    uint8_t repacked[32];
    cfx_ge25519_pack(repacked, &unpacked);
    CFX_ASSERT(memcmp(eight_packed, repacked, 32) == 0);
}

/* test scalarmult with scalar 2^128 */
static void test_scalarmult_2_128(void) {
    uint8_t scalar[32] = {0};
    scalar[16] = 1;  /* 2^128 */

    ge25519_t result;
    cfx_ge25519_scalarmult_base(&result, scalar);

    CFX_ASSERT(cfx_ge25519_is_on_curve(&result));
    CFX_ASSERT(!cfx_ge25519_is_identity(&result));
}

/* test scalarmult with scalar 2^200 */
static void test_scalarmult_2_200(void) {
    uint8_t scalar[32] = {0};
    scalar[25] = 1;  /* 2^200 */

    ge25519_t result;
    cfx_ge25519_scalarmult_base(&result, scalar);

    CFX_ASSERT(cfx_ge25519_is_on_curve(&result));
    CFX_ASSERT(!cfx_ge25519_is_identity(&result));
}

/* test all points from addition chain are distinct */
static void test_addition_chain_distinct(void) {
    uint8_t points[10][32];
    ge25519_t acc;
    ge25519_cached_t b_cached;

    cfx_ge25519_copy(&acc, &cfx_ge25519_base);
    cfx_ge25519_pack(points[0], &acc);
    cfx_ge25519_to_cached(&b_cached, &cfx_ge25519_base);

    for (int i = 1; i < 10; i++) {
        cfx_ge25519_add_cached(&acc, &acc, &b_cached);
        cfx_ge25519_pack(points[i], &acc);
    }

    /* check all pairs are distinct */
    for (int i = 0; i < 10; i++) {
        for (int j = i + 1; j < 10; j++) {
            CFX_ASSERT(memcmp(points[i], points[j], 32) != 0);
        }
    }
}

/* test negated point has different x but same y */
static void test_neg_changes_x_only(void) {
    ge25519_t neg;
    cfx_ge25519_neg(&neg, &cfx_ge25519_base);

    /* after packing, y should be same, sign bit different */
    uint8_t b_packed[32], neg_packed[32];
    cfx_ge25519_pack(b_packed, &cfx_ge25519_base);
    cfx_ge25519_pack(neg_packed, &neg);

    /* y coordinates (lower 255 bits) should match */
    b_packed[31] &= 0x7f;
    neg_packed[31] &= 0x7f;
    CFX_ASSERT(memcmp(b_packed, neg_packed, 32) == 0);
}

int main(void) {
    CFX_TEST(test_identity);
    CFX_TEST(test_basepoint_on_curve);
    CFX_TEST(test_basepoint_not_identity);
    CFX_TEST(test_rfc8032_basepoint_y);

    CFX_TEST(test_pack_unpack_basepoint);
    CFX_TEST(test_pack_unpack_2B);
    CFX_TEST(test_pack_unpack_various);
    CFX_TEST(test_unpack_identity);
    CFX_TEST(test_unpack_invalid);
    CFX_TEST(test_unpack_random_invalid);

    CFX_TEST(test_scalarmult_zero);
    CFX_TEST(test_scalarmult_one);
    CFX_TEST(test_scalarmult_two_vs_double);
    CFX_TEST(test_scalarmult_three);
    CFX_TEST(test_scalarmult_four);
    CFX_TEST(test_scalarmult_five);
    CFX_TEST(test_scalarmult_eight);
    CFX_TEST(test_scalarmult_sixteen);
    CFX_TEST(test_scalarmult_255);
    CFX_TEST(test_scalarmult_256);
    CFX_TEST(test_scalarmult_0x12345);
    CFX_TEST(test_scalarmult_large);
    CFX_TEST(test_scalarmult_2_128);
    CFX_TEST(test_scalarmult_2_200);
    CFX_TEST(test_scalarmult_distributive);
    CFX_TEST(test_scalarmult_associative);
    CFX_TEST(test_scalarmult_non_base);
    CFX_TEST(test_scalarmult_on_curve_stress);
    CFX_TEST(test_group_order);
    CFX_TEST(test_L_minus_1);

    CFX_TEST(test_add_equals_double);
    CFX_TEST(test_sub_self_is_identity);
    CFX_TEST(test_sub_basic);
    CFX_TEST(test_sub_2B_minus_B);
    CFX_TEST(test_sub_3B_minus_2B);
    CFX_TEST(test_negation);
    CFX_TEST(test_double_neg);
    CFX_TEST(test_neg_identity);
    CFX_TEST(test_neg_changes_x_only);
    CFX_TEST(test_add_associative);
    CFX_TEST(test_add_commutative);
    CFX_TEST(test_double_on_curve);
    CFX_TEST(test_repeated_doubling_on_curve);
    CFX_TEST(test_double_identity);
    CFX_TEST(test_identity_add);
    CFX_TEST(test_add_identity_right);
    CFX_TEST(test_many_additions_on_curve);
    CFX_TEST(test_many_doublings_on_curve);
    CFX_TEST(test_addition_chain_distinct);

    CFX_TEST(test_known_2B);
    CFX_TEST(test_eight_B_properties);

    CFX_TEST(test_copy);
    CFX_TEST(test_cneg);
    CFX_TEST(test_cneg_multiple);
    CFX_TEST(test_cmov);
    CFX_TEST(test_cmov_chain);
    CFX_TEST(test_to_cached_preserves);

    return 0;
}
