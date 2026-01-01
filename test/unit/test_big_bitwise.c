/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */
/* test_big_bitwise.c - Tests for big integer bitwise operations */

#include "test_common.h"

/* ========== Binary String Tests ========== */

static void test_binstring(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    cfx_big_from_bin(&b, "111");
    CFX_ASSERT(b.n == 1);
    CFX_ASSERT(b.limb[0] == 7);
    cfx_big_free(&b);
}

static void test_binstring2(void) {
    const char* s =
        "0b110101010101010101010101010101010101010101010101010101010101010101010101"
        "01010101010101010101010101010101010101010101010101010101010101010101010101"
        "010101010101010101010101001010101010101010101010101001100111010101010101010"
        "101010101010101010101010101010101010101010101010101010101010101010101010110"
        "101010101010101010101010101010101010101010101010111101010010100000010000000"
        "010100000100101001001111001010101010100000101011111111111111111111111111111"
        "111111111110";
    cfx_big_t b;
    cfx_big_init(&b);
    cfx_big_from_bin(&b,s);
    char* hex = cfx_big_hex_alloc(&b, NULL);
    const char* expected = "35555555555555555555555555555555555555555552aaaaaa6755555"
        "55555555555555555aaaaaaaaaaaabd4a040282527955415fffffffffe";
    int cmp = strcmp(hex, expected);
    CFX_ASSERT(cmp == 0);
    free(hex);
    cfx_big_free(&b);
}

static void test_binstring3(void) {
    const char* expected =
        "110101010101010101010101010101010101010101010101010101010101010101010101"
        "01010101010101010101010101010101010101010101010101010101010101010101010101"
        "010101010101010101010101001010101010101010101010101001100111010101010101010"
        "101010101010101010101010101010101010101010101010101010101010101010101010110"
        "101010101010101010101010101010101010101010101010111101010010100000010000000"
        "010100000100101001001111001010101010100000101011111111111111111111111111111"
        "111111111110";

    const char* s = "35555555555555555555555555555555555555555552aaaaaa6755555"
        "55555555555555555aaaaaaaaaaaabd4a040282527955415fffffffffe";

    cfx_big_t b;
    cfx_big_init(&b);
    cfx_big_from_hex(&b,s);
    char* bin = cfx_big_bin_alloc(&b, NULL);

    int cmp = strcmp(bin, expected);
    CFX_ASSERT(cmp == 0);
    free(bin);
    cfx_big_free(&b);
}

static void test_bitsset(void) {
    cfx_big_t p;
    cfx_big_init(&p);
    cfx_big_from_u64(&p, 1);

    for (size_t k = 1; k <= 512; ++k) {
        cfx_big_shl_bits_eq(&p, 1);
        CFX_ASSERT(cfx_big_bit_is_set(&p, k));
    }
    cfx_big_free(&p);
}

static void test_bit_set(void) {
    cfx_big_t x;
    cfx_big_init(&x);
    cfx_big_assign_zero(&x);

    cfx_big_bit_set(&x, 0);
    CFX_ASSERT(cfx_big_eq_u64(&x, 1));

    cfx_big_bit_set(&x, 3);
    CFX_ASSERT(cfx_big_eq_u64(&x, 9));

    cfx_big_assign_zero(&x);
    cfx_big_bit_set(&x, 100);
    CFX_ASSERT(cfx_big_bit_is_set(&x, 100));
    CFX_ASSERT(!cfx_big_bit_is_set(&x, 99));
    CFX_ASSERT(!cfx_big_bit_is_set(&x, 101));

    cfx_big_bit_set(&x, 100);  /* idempotent */
    CFX_ASSERT(cfx_big_popcount(&x) == 1);

    cfx_big_free(&x);
}

static void test_bit_clear(void) {
    cfx_big_t x;
    cfx_big_init(&x);
    cfx_big_from_u64(&x, 0xFF);

    cfx_big_bit_clear(&x, 0);
    CFX_ASSERT(cfx_big_eq_u64(&x, 0xFE));

    cfx_big_bit_clear(&x, 7);
    CFX_ASSERT(cfx_big_eq_u64(&x, 0x7E));

    cfx_big_bit_clear(&x, 0);  /* already clear */
    CFX_ASSERT(cfx_big_eq_u64(&x, 0x7E));

    cfx_big_bit_clear(&x, 1000);  /* beyond size, no-op */
    CFX_ASSERT(cfx_big_eq_u64(&x, 0x7E));

    cfx_big_free(&x);
}

static void test_bit_flip(void) {
    cfx_big_t x;
    cfx_big_init(&x);
    cfx_big_from_u64(&x, 0b1010);

    cfx_big_bit_flip(&x, 0);
    CFX_ASSERT(cfx_big_eq_u64(&x, 0b1011));

    cfx_big_bit_flip(&x, 0);
    CFX_ASSERT(cfx_big_eq_u64(&x, 0b1010));

    cfx_big_bit_flip(&x, 3);
    CFX_ASSERT(cfx_big_eq_u64(&x, 0b0010));

    cfx_big_assign_zero(&x);
    cfx_big_bit_flip(&x, 200);
    CFX_ASSERT(cfx_big_bit_is_set(&x, 200));
    cfx_big_bit_flip(&x, 200);
    CFX_ASSERT(cfx_big_is_zero(&x));

    cfx_big_free(&x);
}

static void test_popcount(void) {
    cfx_big_t x;
    cfx_big_init(&x);

    cfx_big_assign_zero(&x);
    CFX_ASSERT(cfx_big_popcount(&x) == 0);

    cfx_big_from_u64(&x, 1);
    CFX_ASSERT(cfx_big_popcount(&x) == 1);

    cfx_big_from_u64(&x, 0xFF);
    CFX_ASSERT(cfx_big_popcount(&x) == 8);

    cfx_big_from_u64(&x, 0xFFFFFFFFFFFFFFFFULL);
    CFX_ASSERT(cfx_big_popcount(&x) == 64);

    cfx_big_assign_zero(&x);
    cfx_big_bit_set(&x, 0);
    cfx_big_bit_set(&x, 64);
    cfx_big_bit_set(&x, 128);
    cfx_big_bit_set(&x, 200);
    CFX_ASSERT(cfx_big_popcount(&x) == 4);

    cfx_big_free(&x);
}

static void test_bit_ops_cross_limb(void) {
    cfx_big_t x;
    cfx_big_init(&x);

    size_t bits[] = {CFX_LIMB_BITS - 1, CFX_LIMB_BITS, CFX_LIMB_BITS + 1,
                     2 * CFX_LIMB_BITS - 1, 2 * CFX_LIMB_BITS};
    size_t n = sizeof(bits) / sizeof(bits[0]);

    for (size_t i = 0; i < n; ++i) {
        cfx_big_assign_zero(&x);
        cfx_big_bit_set(&x, bits[i]);
        CFX_ASSERT(cfx_big_bit_is_set(&x, bits[i]));
        CFX_ASSERT(cfx_big_popcount(&x) == 1);

        cfx_big_bit_clear(&x, bits[i]);
        CFX_ASSERT(!cfx_big_bit_is_set(&x, bits[i]));
        CFX_ASSERT(cfx_big_is_zero(&x));
    }

    cfx_big_free(&x);
}

/* ========== XOR Tests ========== */

static void test_xor_with_zero_is_identity(void) {
    cfx_big_t a = make_u64(0x12345678abcdef00ULL);

    cfx_big_t z;
    cfx_big_init(&z);
    cfx_big_assign_zero(&z);

    cfx_big_t out;
    cfx_big_init(&out);

    cfx_big_xor(&out, &a, &z);
    CFX_ASSERT(cfx_big_eq(&out, &a));

    cfx_big_xor(&out, &z, &a);
    CFX_ASSERT(cfx_big_eq(&out, &a));

    cfx_big_free(&out);
    cfx_big_free(&z);
    cfx_big_free(&a);
}

static void test_xor_self_is_zero(void) {
    cfx_big_t a = make_u64(0xdeadbeefcafebabeULL);
    cfx_big_t b;
    cfx_big_init(&b);
    cfx_big_assign(&b, &a);

    cfx_big_t out;
    cfx_big_init(&out);

    cfx_big_xor(&out, &a, &b);

    cfx_big_t z;
    cfx_big_init(&z);
    cfx_big_assign_zero(&z);

    CFX_ASSERT(cfx_big_eq(&out, &z));

    cfx_big_free(&z);
    cfx_big_free(&out);
    cfx_big_free(&b);
    cfx_big_free(&a);
}

static void test_xor_commutative(void) {
    cfx_big_t a = make_u64(0x0f0e0d0c0b0a0908ULL);
    cfx_big_t b = make_u64(0x0102030405060708ULL);

    cfx_big_t out1, out2;
    cfx_big_init(&out1);
    cfx_big_init(&out2);

    cfx_big_xor(&out1, &a, &b);
    cfx_big_xor(&out2, &b, &a);

    CFX_ASSERT(cfx_big_eq(&out1, &out2));

    cfx_big_free(&out2);
    cfx_big_free(&out1);
    cfx_big_free(&b);
    cfx_big_free(&a);
}

static void test_xor_different_lengths_keep_high_limb_from_larger(void) {
    cfx_big_t a = make_u64(0x0102030405060708ULL);

    cfx_big_t b;
    cfx_big_init(&b);
    cfx_big_from_hex(&b, "f0e0d0c0b0a0908");

    cfx_big_t x;
    cfx_big_init(&x);
    cfx_big_assign_one(&x);
    cfx_big_shl_bits_eq(&x, 80);

    CFX_ASSERT(!cfx_big_bit_is_set(&x, 79));
    CFX_ASSERT(!cfx_big_bit_is_set(&x, 81));
    cfx_big_or(&b, &b, &x);

    cfx_big_t out;
    cfx_big_init(&out);

    cfx_big_xor(&out, &a, &b);

    CFX_ASSERT(cfx_big_endswith_u64(&out, 0x0e0c0e080e0c0e00ULL));
    CFX_ASSERT(cfx_big_bit_is_set(&out, 80));

    cfx_big_free(&x);
    cfx_big_free(&out);
    cfx_big_free(&b);
    cfx_big_free(&a);
}

static void test_xor_out_aliases_a(void) {
    cfx_big_t a = make_u64(0x1111111111111111ULL);
    cfx_big_t b = make_u64(0x0102030405060708ULL);

    cfx_big_xor(&a, &a, &b);

    cfx_big_t expected = make_u64(0x1013121514171619ULL);
    CFX_ASSERT(cfx_big_eq(&a, &expected));

    cfx_big_free(&expected);
    cfx_big_free(&b);
    cfx_big_free(&a);
}

static void test_xor_out_aliases_b(void) {
    cfx_big_t a = make_u64(0x1111111111111111ULL);
    cfx_big_t b = make_u64(0x0102030405060708ULL);

    cfx_big_xor(&b, &a, &b);

    cfx_big_t expected = make_u64(0x1013121514171619ULL);
    CFX_ASSERT(cfx_big_eq(&b, &expected));

    cfx_big_free(&expected);
    cfx_big_free(&b);
    cfx_big_free(&a);
}

/* ========== AND Tests ========== */

static void test_and_basic_u64(void){
    cfx_big_t a, b, out;
    cfx_big_init(&a);
    cfx_big_init(&b);
    cfx_big_init(&out);

    cfx_big_from_u64(&a, 0xF0F0ULL);
    cfx_big_from_u64(&b, 0x0FF0ULL);

    cfx_big_and(&out, &a, &b);
    /* 0xF0F0 & 0x0FF0 = 0x00F0 */
    assert_big_eq_hex(&out, "f0");

    cfx_big_free(&out);
    cfx_big_free(&b);
    cfx_big_free(&a);
}

static void test_and_commutative(void) {
    cfx_big_t a, b, out1, out2;
    cfx_big_init(&a);
    cfx_big_init(&b);
    cfx_big_init(&out1);
    cfx_big_init(&out2);

    cfx_big_from_u64(&a, 0x123456789abcdef0ULL);
    cfx_big_from_u64(&b, 0x0f0f0f0f0f0f0f0fULL);

    cfx_big_and(&out1, &a, &b);
    cfx_big_and(&out2, &b, &a);

    CFX_ASSERT(cfx_big_eq(&out1, &out2));

    cfx_big_free(&out2);
    cfx_big_free(&out1);
    cfx_big_free(&b);
    cfx_big_free(&a);
}

static void test_and_different_lengths_expected_zero(void) {
    cfx_big_t a, b, out;
    cfx_big_init(&a);
    cfx_big_init(&b);
    cfx_big_init(&out);

    cfx_big_from_u64(&a, 0xFFFFFFFFFFFFFFFFULL);
    big_from_hex(&b, "100000000000000000000000000000000"); /* 1<<128 */

    cfx_big_and(&out, &a, &b);
    CFX_ASSERT(cfx_big_is_zero(&out));

    cfx_big_free(&out);
    cfx_big_free(&b);
    cfx_big_free(&a);
}

static void test_and_alias_out_is_larger_operand(void) {
    cfx_big_t small, large, expected;
    cfx_big_init(&small);
    cfx_big_init(&large);
    cfx_big_init(&expected);

    big_from_hex(&small, "00ff00ff");
    big_from_hex(&large, "123400ff00aa");
    big_from_hex(&expected, "00ff00aa");

    cfx_big_and(&large, &small, &large);

    CFX_ASSERT(cfx_big_eq(&large, &expected));

    cfx_big_free(&expected);
    cfx_big_free(&large);
    cfx_big_free(&small);
}

static void test_and_alias_out_is_smaller_operand(void) {
    cfx_big_t small, large, expected;
    cfx_big_init(&small);
    cfx_big_init(&large);
    cfx_big_init(&expected);

    big_from_hex(&small, "00ff00ff");
    big_from_hex(&large, "123400ff00aa");
    big_from_hex(&expected, "00ff00aa");

    cfx_big_and(&small, &small, &large);

    CFX_ASSERT(cfx_big_eq(&small, &expected));

    cfx_big_free(&expected);
    cfx_big_free(&large);
    cfx_big_free(&small);
}

/* ========== OR Tests ========== */

static void test_or_basic_u64(void) {
    cfx_big_t a, b, out;
    cfx_big_init(&a);
    cfx_big_init(&b);
    cfx_big_init(&out);

    cfx_big_from_u64(&a, 0xF0F0ULL);
    cfx_big_from_u64(&b, 0x0FF0ULL);

    cfx_big_or(&out, &a, &b);
    /* 0xF0F0 | 0x0FF0 = 0xFFF0 */
    assert_big_eq_hex(&out, "fff0");

    cfx_big_free(&out);
    cfx_big_free(&b);
    cfx_big_free(&a);
}

static void test_or_commutative(void) {
    cfx_big_t a, b, out1, out2;
    cfx_big_init(&a);
    cfx_big_init(&b);
    cfx_big_init(&out1);
    cfx_big_init(&out2);

    big_from_hex(&a, "123456789abcdef0");
    big_from_hex(&b, "0f0f0f0f0f0f0f0f");

    cfx_big_or(&out1, &a, &b);
    cfx_big_or(&out2, &b, &a);

    CFX_ASSERT(cfx_big_eq(&out1, &out2));

    cfx_big_free(&out2);
    cfx_big_free(&out1);
    cfx_big_free(&b);
    cfx_big_free(&a);
}

static void test_or_different_lengths_keeps_high_limbs(void) {
    cfx_big_t small, large, out, expected;
    cfx_big_init(&small);
    cfx_big_init(&large);
    cfx_big_init(&out);
    cfx_big_init(&expected);

    big_from_hex(&small, "00ff00ff");
    big_from_hex(&large, "123400ff00aa");
    big_from_hex(&expected, "123400ff00ff");

    cfx_big_or(&out, &small, &large);
    CFX_ASSERT(cfx_big_eq(&out, &expected));

    cfx_big_free(&expected);
    cfx_big_free(&out);
    cfx_big_free(&large);
    cfx_big_free(&small);
}

static void test_or_alias_out_is_smaller_operand(void) {
    cfx_big_t small, large, expected;
    cfx_big_init(&small);
    cfx_big_init(&large);
    cfx_big_init(&expected);

    big_from_hex(&small, "00ff00ff");
    big_from_hex(&large, "123400ff00aa");
    big_from_hex(&expected, "123400ff00ff");

    cfx_big_or(&small, &small, &large);
    CFX_ASSERT(cfx_big_eq(&small, &expected));

    cfx_big_free(&expected);
    cfx_big_free(&large);
    cfx_big_free(&small);
}

static void test_or_alias_out_is_larger_operand(void) {
    cfx_big_t small, large, expected;
    cfx_big_init(&small);
    cfx_big_init(&large);
    cfx_big_init(&expected);

    big_from_hex(&small, "00ff00ff");
    big_from_hex(&large, "123400ff00aa");
    big_from_hex(&expected, "123400ff00ff");

    cfx_big_or(&large, &small, &large);
    CFX_ASSERT(cfx_big_eq(&large, &expected));

    cfx_big_free(&expected);
    cfx_big_free(&large);
    cfx_big_free(&small);
}

/* ========== Rotate Tests ========== */

static void test_rot_zero_and_one(void) {
    cfx_big_t x, y;
    cfx_big_init(&x);
    cfx_big_init(&y);

    cfx_big_from_limb(&x, 0);
    cfx_big_rotl(&y, &x, 5);
    ASSERT_BIG_HEX(&y, "0");

    cfx_big_from_limb(&x, 1);
    cfx_big_rotl(&y, &x, 0);
    ASSERT_BIG_HEX(&y, "1");

    cfx_big_rotr(&y, &x, 0);
    ASSERT_BIG_HEX(&y, "1");

    cfx_big_free(&x);
    cfx_big_free(&y);
}

static void test_rot_8bit(void) {
    cfx_big_t x, y;
    cfx_big_init(&x);
    cfx_big_init(&y);

    /* x = 0b10000001 */
    cfx_big_from_bin(&x, "10000001");
    ASSERT_BIG_HEX(&x, "81");

    cfx_big_rotl(&y, &x, 1);
    ASSERT_BIG_HEX(&y, "3"); /* 00000011 */

    cfx_big_rotr(&y, &x, 1);
    ASSERT_BIG_HEX(&y, "c0"); /* 11000000 */

    cfx_big_free(&x);
    cfx_big_free(&y);
}

static void test_rot_by_width(void) {
    cfx_big_t x, y;
    cfx_big_init(&x);
    cfx_big_init(&y);

    cfx_big_from_hex(&x, "deadbeef");
    unsigned w = cfx_big_bitlen(&x);

    cfx_big_rotl(&y, &x, w);
    CFX_ASSERT(cfx_big_cmp(&x, &y) == 0);

    cfx_big_rotr(&y, &x, w);
    CFX_ASSERT(cfx_big_cmp(&x, &y) == 0);

    cfx_big_free(&x);
    cfx_big_free(&y);
}

static void test_rot_duality(void) {
    cfx_big_t x, a, b;
    cfx_big_init(&x);
    cfx_big_init(&a);
    cfx_big_init(&b);

    cfx_big_from_hex(&x, "123456789ABCDEF");
    unsigned w = cfx_big_bitlen(&x);
    unsigned r = 13;

    cfx_big_rotl(&a, &x, r);
    cfx_big_rotr(&b, &x, w - r);

    CFX_ASSERT(cfx_big_cmp(&a, &b) == 0);

    cfx_big_free(&x);
    cfx_big_free(&a);
    cfx_big_free(&b);
}

static void test_rot_composition(void) {
    cfx_big_t x, a, b;
    cfx_big_init(&x);
    cfx_big_init(&a);
    cfx_big_init(&b);

    cfx_big_from_hex(&x, "CAFEBABEDEADC0DE");
    unsigned w = cfx_big_bitlen(&x);

    cfx_big_rotl_w(&a, &x, 7, w);
    cfx_big_rotl_w(&a, &a, 19, w);

    cfx_big_rotl_w(&b, &x, (7 + 19) % w, w);

    CFX_ASSERT(cfx_big_cmp(&a, &b) == 0);

    cfx_big_free(&x);
    cfx_big_free(&a);
    cfx_big_free(&b);
}

static void test_rot_aliasing(void) {
    cfx_big_t x, y;
    cfx_big_init(&x);
    cfx_big_init(&y);

    cfx_big_from_hex(&x, "F0F0F0F0");
    cfx_big_copy(&y, &x);

    cfx_big_rotl(&x, &x, 4);
    cfx_big_rotl(&y, &y, 4);

    CFX_ASSERT(cfx_big_cmp(&x, &y) == 0);

    cfx_big_free(&x);
    cfx_big_free(&y);
}

static void test_rot_cross_limb(void) {
    cfx_big_t x, y;
    cfx_big_init(&x);
    cfx_big_init(&y);

    /* bit exactly at limb boundary */
    cfx_big_from_limb(&x, 1);
    cfx_big_shl_bits(&x, &x, CFX_LIMB_BITS);
    cfx_big_mask_bits(&x, CFX_LIMB_BITS + 1);

    cfx_big_rotl(&y, &x, 1);

    /* should wrap the top bit back to LSB */
    ASSERT_BIG_HEX(&y, "1");

    cfx_big_free(&x);
    cfx_big_free(&y);
}

static void test_rotl_w_basic(void) {
    cfx_big_t x, y;
    cfx_big_init(&x);
    cfx_big_init(&y);

    cfx_big_from_hex(&x, "ab");
    cfx_big_rotl_w(&y, &x, 4, 8);
    ASSERT_BIG_HEX(&y, "ba");

    cfx_big_free(&x);
    cfx_big_free(&y);
}

static void test_rotr_w_basic(void) {
    cfx_big_t x, y;
    cfx_big_init(&x);
    cfx_big_init(&y);

    cfx_big_from_hex(&x, "ab");
    cfx_big_rotr_w(&y, &x, 4, 8);
    ASSERT_BIG_HEX(&y, "ba");

    cfx_big_free(&x);
    cfx_big_free(&y);
}

static void test_rotl_w_larger_width(void) {
    cfx_big_t x, y, z;
    cfx_big_init(&x);
    cfx_big_init(&y);
    cfx_big_init(&z);

    cfx_big_from_hex(&x, "0f");
    cfx_big_rotl_w(&y, &x, 4, 8);
    ASSERT_BIG_HEX(&y, "f0");

    cfx_big_rotr_w(&z, &y, 4, 8);
    ASSERT_BIG_HEX(&z, "f");

    cfx_big_free(&x);
    cfx_big_free(&y);
    cfx_big_free(&z);
}

static void test_rot_w_zero_rotation(void) {
    cfx_big_t x, y;
    cfx_big_init(&x);
    cfx_big_init(&y);

    cfx_big_from_hex(&x, "DEADBEEF");
    cfx_big_rotl_w(&y, &x, 0, 32);
    CFX_ASSERT(cfx_big_cmp(&x, &y) == 0);

    cfx_big_rotr_w(&y, &x, 0, 32);
    CFX_ASSERT(cfx_big_cmp(&x, &y) == 0);

    cfx_big_free(&x);
    cfx_big_free(&y);
}

static void test_rot_w_full_rotation(void) {
    cfx_big_t x, y;
    cfx_big_init(&x);
    cfx_big_init(&y);

    cfx_big_from_hex(&x, "DEADBEEF");
    unsigned w = 32;

    cfx_big_rotl_w(&y, &x, w, w);
    CFX_ASSERT(cfx_big_cmp(&x, &y) == 0);

    cfx_big_rotr_w(&y, &x, w, w);
    CFX_ASSERT(cfx_big_cmp(&x, &y) == 0);

    cfx_big_rotl_w(&y, &x, 2 * w, w);
    CFX_ASSERT(cfx_big_cmp(&x, &y) == 0);

    cfx_big_free(&x);
    cfx_big_free(&y);
}

static void test_rot_w_aliasing(void) {
    cfx_big_t x, y;
    cfx_big_init(&x);
    cfx_big_init(&y);

    cfx_big_from_hex(&x, "CAFEBABE");
    cfx_big_copy(&y, &x);

    cfx_big_rotl_w(&x, &x, 12, 32);
    cfx_big_rotl_w(&y, &y, 12, 32);
    CFX_ASSERT(cfx_big_cmp(&x, &y) == 0);

    cfx_big_rotr_w(&x, &x, 12, 32);
    cfx_big_from_hex(&y, "CAFEBABE");
    CFX_ASSERT(cfx_big_cmp(&x, &y) == 0);

    cfx_big_free(&x);
    cfx_big_free(&y);
}

static void test_rot_w_duality(void) {
    cfx_big_t x, a, b;
    cfx_big_init(&x);
    cfx_big_init(&a);
    cfx_big_init(&b);

    cfx_big_from_hex(&x, "123456789ABCDEF0");
    unsigned w = 64;
    unsigned r = 17;

    cfx_big_rotl_w(&a, &x, r, w);
    cfx_big_rotr_w(&b, &x, w - r, w);
    CFX_ASSERT(cfx_big_cmp(&a, &b) == 0);

    cfx_big_free(&x);
    cfx_big_free(&a);
    cfx_big_free(&b);
}

static void test_rot_w_power_of_two_width(void) {
    cfx_big_t x, y;
    cfx_big_init(&x);
    cfx_big_init(&y);

    cfx_big_from_hex(&x, "FF");
    cfx_big_rotl_w(&y, &x, 8, 16);
    ASSERT_BIG_HEX(&y, "FF00");

    cfx_big_from_hex(&x, "FF");
    cfx_big_rotl_w(&y, &x, 24, 32);
    ASSERT_BIG_HEX(&y, "FF000000");

    cfx_big_free(&x);
    cfx_big_free(&y);
}

static void test_and_eq(void) {
    cfx_big_t a, b;
    cfx_big_init(&a);
    cfx_big_init(&b);
    cfx_big_from_limb(&a, 0xFF);
    cfx_big_from_limb(&b, 0x0F);
    cfx_big_and_eq(&a, &b);
    CFX_ASSERT(cfx_big_eq_u64(&a, 0x0F));
    cfx_big_free(&b);
    cfx_big_free(&a);
}

static void test_or_eq(void) {
    cfx_big_t a, b;
    cfx_big_init(&a);
    cfx_big_init(&b);
    cfx_big_from_limb(&a, 0xF0);
    cfx_big_from_limb(&b, 0x0F);
    cfx_big_or_eq(&a, &b);
    CFX_ASSERT(cfx_big_eq_u64(&a, 0xFF));
    cfx_big_free(&b);
    cfx_big_free(&a);
}

static void test_xor_eq(void) {
    cfx_big_t a, b;
    cfx_big_init(&a);
    cfx_big_init(&b);
    cfx_big_from_limb(&a, 0xFF);
    cfx_big_from_limb(&b, 0x0F);
    cfx_big_xor_eq(&a, &b);
    CFX_ASSERT(cfx_big_eq_u64(&a, 0xF0));
    cfx_big_free(&b);
    cfx_big_free(&a);
}

static void test_rotl_w_zero_width(void) {
    cfx_big_t a, out;
    cfx_big_init(&a);
    cfx_big_init(&out);
    cfx_big_from_limb(&a, 0xFF);
    cfx_big_rotl_w(&out, &a, 3, 0);  /* w == 0 */
    CFX_ASSERT(cfx_big_is_zero(&out));
    cfx_big_free(&out);
    cfx_big_free(&a);
}

static void test_rotr_w_zero_width(void) {
    cfx_big_t a, out;
    cfx_big_init(&a);
    cfx_big_init(&out);
    cfx_big_from_limb(&a, 0xFF);
    cfx_big_rotr_w(&out, &a, 3, 0);  /* w == 0 */
    CFX_ASSERT(cfx_big_is_zero(&out));
    cfx_big_free(&out);
    cfx_big_free(&a);
}

static void test_mask_bits_zero(void) {
    cfx_big_t a;
    cfx_big_init(&a);
    cfx_big_from_limb(&a, 0xFFFF);
    cfx_big_mask_bits(&a, 0);
    CFX_ASSERT(cfx_big_is_zero(&a));
    cfx_big_free(&a);
}

int main(void) {
    /* Binary string tests */
    CFX_TEST(test_binstring);
    CFX_TEST(test_binstring2);
    CFX_TEST(test_binstring3);
    CFX_TEST(test_bitsset);

    /* Bit manipulation */
    CFX_TEST(test_bit_set);
    CFX_TEST(test_bit_clear);
    CFX_TEST(test_bit_flip);
    CFX_TEST(test_popcount);
    CFX_TEST(test_bit_ops_cross_limb);

    /* XOR tests */
    CFX_TEST(test_xor_with_zero_is_identity);
    CFX_TEST(test_xor_self_is_zero);
    CFX_TEST(test_xor_commutative);
    CFX_TEST(test_xor_different_lengths_keep_high_limb_from_larger);
    CFX_TEST(test_xor_out_aliases_a);
    CFX_TEST(test_xor_out_aliases_b);

    /* AND tests */
    CFX_TEST(test_and_basic_u64);
    CFX_TEST(test_and_commutative);
    CFX_TEST(test_and_different_lengths_expected_zero);
    CFX_TEST(test_and_alias_out_is_larger_operand);
    CFX_TEST(test_and_alias_out_is_smaller_operand);

    /* OR tests */
    CFX_TEST(test_or_basic_u64);
    CFX_TEST(test_or_commutative);
    CFX_TEST(test_or_different_lengths_keeps_high_limbs);
    CFX_TEST(test_or_alias_out_is_smaller_operand);
    CFX_TEST(test_or_alias_out_is_larger_operand);

    /* Rotate tests */
    CFX_TEST(test_rot_zero_and_one);
    CFX_TEST(test_rot_8bit);
    CFX_TEST(test_rot_by_width);
    CFX_TEST(test_rot_duality);
    CFX_TEST(test_rot_composition);
    CFX_TEST(test_rot_aliasing);
    CFX_TEST(test_rot_cross_limb);
    CFX_TEST(test_rotl_w_basic);
    CFX_TEST(test_rotr_w_basic);
    CFX_TEST(test_rotl_w_larger_width);
    CFX_TEST(test_rot_w_zero_rotation);
    CFX_TEST(test_rot_w_full_rotation);
    CFX_TEST(test_rot_w_aliasing);
    CFX_TEST(test_rot_w_duality);
    CFX_TEST(test_rot_w_power_of_two_width);

    /* In-place bitwise ops */
    CFX_TEST(test_and_eq);
    CFX_TEST(test_or_eq);
    CFX_TEST(test_xor_eq);

    /* Rotate edge cases */
    CFX_TEST(test_rotl_w_zero_width);
    CFX_TEST(test_rotr_w_zero_width);

    /* Mask bits */
    CFX_TEST(test_mask_bits_zero);

    puts("OK");
    return 0;
}
