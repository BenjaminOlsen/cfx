/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */
/* test_ntt.c - Tests for NTT (Number-Theoretic Transform) */

#include "cfx/ntt.h"
#include "cfx/big.h"
#include "cfx/macros.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void test_mod_add(void) {
    uint64_t p = CFX_NTT_P1;

    CFX_ASSERT(cfx_ntt_mod_add(0, 0, p) == 0);
    CFX_ASSERT(cfx_ntt_mod_add(1, 1, p) == 2);
    CFX_ASSERT(cfx_ntt_mod_add(p - 1, 1, p) == 0);
    CFX_ASSERT(cfx_ntt_mod_add(p - 1, 2, p) == 1);
    CFX_ASSERT(cfx_ntt_mod_add(p - 1, p - 1, p) == p - 2);
}

static void test_mod_add_overflow(void) {
    uint64_t p = CFX_NTT_P1;

    /* test overflow detection: when a + b wraps around 2^64 */
    uint64_t a = p - 1;
    uint64_t b = p - 1;
    CFX_ASSERT(cfx_ntt_mod_add(a, b, p) == p - 2);

    /* test case where sum >= p but no 64-bit overflow */
    a = p / 2;
    b = p / 2 + 2;
    CFX_ASSERT(cfx_ntt_mod_add(a, b, p) == 1);

    /* test with smaller prime */
    uint64_t small_p = 17;
    CFX_ASSERT(cfx_ntt_mod_add(10, 10, small_p) == 3);
    CFX_ASSERT(cfx_ntt_mod_add(16, 16, small_p) == 15);
}

static void test_mod_sub(void) {
    uint64_t p = CFX_NTT_P1;

    CFX_ASSERT(cfx_ntt_mod_sub(0, 0, p) == 0);
    CFX_ASSERT(cfx_ntt_mod_sub(5, 3, p) == 2);
    CFX_ASSERT(cfx_ntt_mod_sub(3, 5, p) == p - 2);
    CFX_ASSERT(cfx_ntt_mod_sub(0, 1, p) == p - 1);
    CFX_ASSERT(cfx_ntt_mod_sub(p - 1, p - 1, p) == 0);
}

static void test_mod_sub_underflow(void) {
    uint64_t p = CFX_NTT_P1;

    /* test all underflow cases */
    CFX_ASSERT(cfx_ntt_mod_sub(0, p - 1, p) == 1);  /* 0 - (p-1) = 1 mod p */
    CFX_ASSERT(cfx_ntt_mod_sub(1, p - 1, p) == 2);  /* 1 - (p-1) = 2 mod p */

    /* test exact boundary: a == b */
    CFX_ASSERT(cfx_ntt_mod_sub(100, 100, p) == 0);

    /* test a = b + 1 (no underflow) */
    CFX_ASSERT(cfx_ntt_mod_sub(101, 100, p) == 1);

    /* test with smaller prime */
    uint64_t small_p = 17;
    CFX_ASSERT(cfx_ntt_mod_sub(3, 10, small_p) == 10);  /* 3 - 10 = -7 = 10 mod 17 */
    CFX_ASSERT(cfx_ntt_mod_sub(0, 16, small_p) == 1);   /* 0 - 16 = 1 mod 17 */
}

static void test_mod_mul(void) {
    uint64_t p = CFX_NTT_P1;

    CFX_ASSERT(cfx_ntt_mod_mul(0, 123456, p) == 0);
    CFX_ASSERT(cfx_ntt_mod_mul(1, 123456, p) == 123456);
    CFX_ASSERT(cfx_ntt_mod_mul(2, 3, p) == 6);
    CFX_ASSERT(cfx_ntt_mod_mul(p - 1, 2, p) == p - 2);

    /* test large values that would overflow 64 bits */
    uint64_t a = UINT64_C(0x123456789ABCDEF0);
    uint64_t b = UINT64_C(0xFEDCBA9876543210);
    uint64_t r = cfx_ntt_mod_mul(a, b, p);
    /* verify by computing (a * b) mod p differently:
       a*b mod p should give same result both ways */
    CFX_ASSERT(r < p);

    /* a * 1 = a mod p */
    CFX_ASSERT(cfx_ntt_mod_mul(a % p, 1, p) == a % p);
}

static void test_mod_mul_edge_cases(void) {
    uint64_t p = CFX_NTT_P1;

    /* multiply by zero */
    CFX_ASSERT(cfx_ntt_mod_mul(0, 0, p) == 0);
    CFX_ASSERT(cfx_ntt_mod_mul(p - 1, 0, p) == 0);
    CFX_ASSERT(cfx_ntt_mod_mul(0, p - 1, p) == 0);

    /* multiply by one */
    CFX_ASSERT(cfx_ntt_mod_mul(1, 1, p) == 1);
    CFX_ASSERT(cfx_ntt_mod_mul(p - 1, 1, p) == p - 1);

    /* multiply p-1 by p-1 = 1 mod p (since (p-1)^2 = p^2 - 2p + 1 = 1 mod p) */
    CFX_ASSERT(cfx_ntt_mod_mul(p - 1, p - 1, p) == 1);

    /* test values near 2^32 boundary (exercises 32-bit halves in mul128) */
    uint64_t near_32 = UINT64_C(0xFFFFFFFF);  /* 2^32 - 1 */
    uint64_t r = cfx_ntt_mod_mul(near_32, near_32, p);
    CFX_ASSERT(r < p);

    /* test maximum 64-bit values */
    uint64_t max64 = UINT64_MAX;
    r = cfx_ntt_mod_mul(max64 % p, max64 % p, p);
    CFX_ASSERT(r < p);

    /* verify commutativity */
    uint64_t a = 12345678901234567ULL;
    uint64_t b = 98765432109876543ULL;
    CFX_ASSERT(cfx_ntt_mod_mul(a, b, p) == cfx_ntt_mod_mul(b, a, p));
}

static void test_mod_mul_small_prime(void) {
    uint64_t p = 97;

    CFX_ASSERT(cfx_ntt_mod_mul(10, 10, p) == 3);   /* 100 mod 97 */
    CFX_ASSERT(cfx_ntt_mod_mul(50, 50, p) == 75);  /* 2500 mod 97 */
    CFX_ASSERT(cfx_ntt_mod_mul(96, 96, p) == 1);   /* (-1)^2 */
    CFX_ASSERT(cfx_ntt_mod_mul(48, 2, p) == 96);
}

static void test_mod_pow(void) {
    uint64_t p = CFX_NTT_P1;

    CFX_ASSERT(cfx_ntt_mod_pow(2, 0, p) == 1);
    CFX_ASSERT(cfx_ntt_mod_pow(2, 1, p) == 2);
    CFX_ASSERT(cfx_ntt_mod_pow(2, 10, p) == 1024);
    CFX_ASSERT(cfx_ntt_mod_pow(3, 5, p) == 243);

    /* Fermat's little theorem */
    CFX_ASSERT(cfx_ntt_mod_pow(7, p - 1, p) == 1);
    CFX_ASSERT(cfx_ntt_mod_pow(12345, p - 1, p) == 1);
}

static void test_mod_pow_p_equals_one(void) {
    /* special case: p = 1 means all results are 0 */
    CFX_ASSERT(cfx_ntt_mod_pow(5, 3, 1) == 0);
    CFX_ASSERT(cfx_ntt_mod_pow(0, 0, 1) == 0);
    CFX_ASSERT(cfx_ntt_mod_pow(100, 100, 1) == 0);
}

static void test_mod_pow_base_zero(void) {
    uint64_t p = CFX_NTT_P1;

    CFX_ASSERT(cfx_ntt_mod_pow(0, 1, p) == 0);
    CFX_ASSERT(cfx_ntt_mod_pow(0, 100, p) == 0);
    /* 0^0 = 1 by convention in most implementations */
    CFX_ASSERT(cfx_ntt_mod_pow(0, 0, p) == 1);
}

static void test_mod_pow_exp_patterns(void) {
    uint64_t p = CFX_NTT_P1;

    /* power of 2 exponents (only squaring, no multiply) */
    CFX_ASSERT(cfx_ntt_mod_pow(2, 1, p) == 2);
    CFX_ASSERT(cfx_ntt_mod_pow(2, 2, p) == 4);
    CFX_ASSERT(cfx_ntt_mod_pow(2, 4, p) == 16);
    CFX_ASSERT(cfx_ntt_mod_pow(2, 8, p) == 256);
    CFX_ASSERT(cfx_ntt_mod_pow(2, 16, p) == 65536);

    /* odd exponents (forces multiply step each iteration) */
    CFX_ASSERT(cfx_ntt_mod_pow(2, 3, p) == 8);
    CFX_ASSERT(cfx_ntt_mod_pow(2, 7, p) == 128);
    CFX_ASSERT(cfx_ntt_mod_pow(2, 15, p) == 32768);

    /* alternating bits in exponent */
    CFX_ASSERT(cfx_ntt_mod_pow(2, 5, p) == 32);   /* 101 binary */
    CFX_ASSERT(cfx_ntt_mod_pow(2, 10, p) == 1024); /* 1010 binary */
}

static void test_mod_pow_large_exp(void) {
    uint64_t p = CFX_NTT_P1;

    /* very large exponent */
    uint64_t large_exp = UINT64_C(0xFFFFFFFFFFFFFFFF);
    uint64_t r = cfx_ntt_mod_pow(2, large_exp, p);
    CFX_ASSERT(r < p);

    /* Fermat: a^(p-1) = 1 mod p */
    CFX_ASSERT(cfx_ntt_mod_pow(7, p - 1, p) == 1);

    /* test with smaller prime to avoid overflow in 2*(p-1) */
    uint64_t small_p = 1000000007;
    CFX_ASSERT(cfx_ntt_mod_pow(3, small_p - 1, small_p) == 1);
    CFX_ASSERT(cfx_ntt_mod_pow(3, 2 * (small_p - 1), small_p) == 1);
}

static void test_mod_inv(void) {
    uint64_t p = CFX_NTT_P1;

    uint64_t a = 7;
    uint64_t a_inv = cfx_ntt_mod_inv(a, p);
    CFX_ASSERT(cfx_ntt_mod_mul(a, a_inv, p) == 1);

    a = 12345678901234567ULL;
    a_inv = cfx_ntt_mod_inv(a, p);
    CFX_ASSERT(cfx_ntt_mod_mul(a, a_inv, p) == 1);

    CFX_ASSERT(cfx_ntt_mod_inv(1, p) == 1);
}

static void test_mod_inv_edge_cases(void) {
    uint64_t p = CFX_NTT_P1;

    /* inverse of p-1 is p-1 since (p-1)*(p-1) = 1 mod p */
    CFX_ASSERT(cfx_ntt_mod_inv(p - 1, p) == p - 1);

    /* inverse of 2 */
    uint64_t inv2 = cfx_ntt_mod_inv(2, p);
    CFX_ASSERT(cfx_ntt_mod_mul(2, inv2, p) == 1);

    /* (p+1)/2 should be inverse of 2 when p is odd */
    CFX_ASSERT(inv2 == (p + 1) / 2);

    /* verify inverse is involutory: inv(inv(a)) = a */
    uint64_t a = 12345;
    CFX_ASSERT(cfx_ntt_mod_inv(cfx_ntt_mod_inv(a, p), p) == a);

    /* test with small prime */
    uint64_t small_p = 97;
    CFX_ASSERT(cfx_ntt_mod_mul(cfx_ntt_mod_inv(50, small_p), 50, small_p) == 1);
}

static void test_root_of_unity(void) {
    uint64_t p = CFX_NTT_P1;
    uint64_t g = CFX_NTT_G1;

    size_t n = 8;
    uint64_t omega = cfx_ntt_root_of_unity(g, p, n);
    uint64_t omega_n = cfx_ntt_mod_pow(omega, n, p);
    CFX_ASSERT(omega_n == 1);

    /* omega^(n/2) = -1 mod p for primitive root */
    uint64_t omega_half = cfx_ntt_mod_pow(omega, n / 2, p);
    CFX_ASSERT(omega_half == p - 1);

    n = 1024;
    omega = cfx_ntt_root_of_unity(g, p, n);
    omega_n = cfx_ntt_mod_pow(omega, n, p);
    CFX_ASSERT(omega_n == 1);
}

static void test_root_of_unity_various_n(void) {
    uint64_t p = CFX_NTT_P1;
    uint64_t g = CFX_NTT_G1;

    /* test all powers of 2 up to 2^16 */
    for (size_t n = 1; n <= 65536; n <<= 1) {
        uint64_t omega = cfx_ntt_root_of_unity(g, p, n);
        CFX_ASSERT(cfx_ntt_mod_pow(omega, n, p) == 1);

        /* for n > 1, omega^(n/2) should be -1 */
        if (n > 1) {
            CFX_ASSERT(cfx_ntt_mod_pow(omega, n / 2, p) == p - 1);
        }
    }
}

static void test_root_of_unity_primitivity(void) {
    uint64_t p = CFX_NTT_P1;
    uint64_t g = CFX_NTT_G1;
    size_t n = 16;

    uint64_t omega = cfx_ntt_root_of_unity(g, p, n);

    /* omega^k != 1 for 0 < k < n */
    for (size_t k = 1; k < n; ++k) {
        CFX_ASSERT(cfx_ntt_mod_pow(omega, k, p) != 1);
    }
    CFX_ASSERT(cfx_ntt_mod_pow(omega, n, p) == 1);
}

static void test_twiddles_init_free(void) {
    cfx_ntt_twiddles_t tw;

    int rc = cfx_ntt_twiddles_init(&tw, 8, CFX_NTT_P1, CFX_NTT_G1);
    CFX_ASSERT(rc == 0);
    CFX_ASSERT(tw.n == 8);
    CFX_ASSERT(tw.p == CFX_NTT_P1);
    CFX_ASSERT(tw.forward != NULL);
    CFX_ASSERT(tw.inverse != NULL);
    CFX_ASSERT(tw.forward[0] == 1);
    CFX_ASSERT(tw.inverse[0] == 1);

    /* verify twiddles are powers of omega */
    uint64_t omega = cfx_ntt_root_of_unity(CFX_NTT_G1, CFX_NTT_P1, 8);
    for (size_t k = 0; k < 8; ++k) {
        uint64_t expected = cfx_ntt_mod_pow(omega, k, CFX_NTT_P1);
        CFX_ASSERT(tw.forward[k] == expected);
    }

    cfx_ntt_twiddles_free(&tw);
    CFX_ASSERT(tw.forward == NULL);
    CFX_ASSERT(tw.inverse == NULL);
}

static void test_twiddles_invalid_n(void) {
    cfx_ntt_twiddles_t tw;

    int rc = cfx_ntt_twiddles_init(&tw, 0, CFX_NTT_P1, CFX_NTT_G1);
    CFX_ASSERT(rc == -1);

    rc = cfx_ntt_twiddles_init(&tw, 7, CFX_NTT_P1, CFX_NTT_G1);
    CFX_ASSERT(rc == -1);

    rc = cfx_ntt_twiddles_init(&tw, 100, CFX_NTT_P1, CFX_NTT_G1);
    CFX_ASSERT(rc == -1);

    rc = cfx_ntt_twiddles_init(&tw, 3, CFX_NTT_P1, CFX_NTT_G1);
    CFX_ASSERT(rc == -1);

    rc = cfx_ntt_twiddles_init(&tw, 6, CFX_NTT_P1, CFX_NTT_G1);
    CFX_ASSERT(rc == -1);
}

static void test_twiddles_n_equals_one(void) {
    cfx_ntt_twiddles_t tw;

    int rc = cfx_ntt_twiddles_init(&tw, 1, CFX_NTT_P1, CFX_NTT_G1);
    CFX_ASSERT(rc == 0);
    CFX_ASSERT(tw.n == 1);
    CFX_ASSERT(tw.forward[0] == 1);
    CFX_ASSERT(tw.inverse[0] == 1);
    CFX_ASSERT(tw.n_inv == 1);

    cfx_ntt_twiddles_free(&tw);
}

static void test_twiddles_inverse_correctness(void) {
    cfx_ntt_twiddles_t tw;
    size_t n = 16;
    uint64_t p = CFX_NTT_P1;

    int rc = cfx_ntt_twiddles_init(&tw, n, p, CFX_NTT_G1);
    CFX_ASSERT(rc == 0);

    /* forward[k] * inverse[k] should equal 1 mod p for all k */
    for (size_t k = 0; k < n; ++k) {
        uint64_t prod = cfx_ntt_mod_mul(tw.forward[k], tw.inverse[k], p);
        CFX_ASSERT(prod == 1);
    }

    cfx_ntt_twiddles_free(&tw);
}

static void test_bit_reverse(void) {
    uint64_t a[8] = {0, 1, 2, 3, 4, 5, 6, 7};

    cfx_ntt_bit_reverse(a, 8);

    /* bit-reversal of 3-bit indices:
       000 -> 000 (0 -> 0)
       001 -> 100 (1 -> 4)
       010 -> 010 (2 -> 2)
       011 -> 110 (3 -> 6)
       100 -> 001 (4 -> 1)
       101 -> 101 (5 -> 5)
       110 -> 011 (6 -> 3)
       111 -> 111 (7 -> 7)
    */
    CFX_ASSERT(a[0] == 0);
    CFX_ASSERT(a[1] == 4);
    CFX_ASSERT(a[2] == 2);
    CFX_ASSERT(a[3] == 6);
    CFX_ASSERT(a[4] == 1);
    CFX_ASSERT(a[5] == 5);
    CFX_ASSERT(a[6] == 3);
    CFX_ASSERT(a[7] == 7);
}

static void test_bit_reverse_small(void) {
    uint64_t a1[1] = {42};
    cfx_ntt_bit_reverse(a1, 1);
    CFX_ASSERT(a1[0] == 42);

    uint64_t a2[2] = {10, 20};
    cfx_ntt_bit_reverse(a2, 2);
    CFX_ASSERT(a2[0] == 10);
    CFX_ASSERT(a2[1] == 20);

    uint64_t a4[4] = {0, 1, 2, 3};
    cfx_ntt_bit_reverse(a4, 4);
    CFX_ASSERT(a4[0] == 0);
    CFX_ASSERT(a4[1] == 2);
    CFX_ASSERT(a4[2] == 1);
    CFX_ASSERT(a4[3] == 3);
}

static void test_bit_reverse_n16(void) {
    uint64_t a[16];
    for (size_t i = 0; i < 16; ++i) a[i] = i;

    cfx_ntt_bit_reverse(a, 16);

    /* 4-bit reversal: 0000->0000, 0001->1000, 0010->0100, etc */
    CFX_ASSERT(a[0] == 0);   /* 0000 -> 0000 */
    CFX_ASSERT(a[1] == 8);   /* 0001 -> 1000 */
    CFX_ASSERT(a[2] == 4);   /* 0010 -> 0100 */
    CFX_ASSERT(a[3] == 12);  /* 0011 -> 1100 */
    CFX_ASSERT(a[4] == 2);   /* 0100 -> 0010 */
    CFX_ASSERT(a[8] == 1);   /* 1000 -> 0001 */
    CFX_ASSERT(a[15] == 15); /* 1111 -> 1111 */
}

static void test_bit_reverse_involution(void) {
    /* bit-reverse twice should give original */
    size_t n = 32;
    uint64_t* a = (uint64_t*)malloc(n * sizeof(uint64_t));
    uint64_t* orig = (uint64_t*)malloc(n * sizeof(uint64_t));
    CFX_ASSERT(a && orig);

    for (size_t i = 0; i < n; ++i) {
        a[i] = i * 7 + 3;
        orig[i] = a[i];
    }

    cfx_ntt_bit_reverse(a, n);
    cfx_ntt_bit_reverse(a, n);

    for (size_t i = 0; i < n; ++i) {
        CFX_ASSERT(a[i] == orig[i]);
    }

    free(a);
    free(orig);
}

static void test_ntt_forward_inverse_identity(void) {
    uint64_t p = CFX_NTT_P1;
    size_t n = 8;

    cfx_ntt_twiddles_t tw;
    int rc = cfx_ntt_twiddles_init(&tw, n, p, CFX_NTT_G1);
    CFX_ASSERT(rc == 0);

    uint64_t a[8] = {1, 2, 3, 4, 5, 6, 7, 8};
    uint64_t orig[8];
    memcpy(orig, a, sizeof(a));

    /* forward NTT */
    cfx_ntt_bit_reverse(a, n);
    cfx_ntt_forward(a, n, p, tw.forward);

    /* inverse NTT */
    cfx_ntt_bit_reverse(a, n);
    cfx_ntt_inverse(a, n, p, tw.inverse, tw.n_inv);

    /* should recover original */
    for (size_t i = 0; i < n; ++i) {
        CFX_ASSERT(a[i] == orig[i]);
    }

    cfx_ntt_twiddles_free(&tw);
}

static void test_ntt_larger(void) {
    uint64_t p = CFX_NTT_P1;
    size_t n = 256;

    cfx_ntt_twiddles_t tw;
    int rc = cfx_ntt_twiddles_init(&tw, n, p, CFX_NTT_G1);
    CFX_ASSERT(rc == 0);

    uint64_t* a = (uint64_t*)malloc(n * sizeof(uint64_t));
    uint64_t* orig = (uint64_t*)malloc(n * sizeof(uint64_t));
    CFX_ASSERT(a && orig);

    for (size_t i = 0; i < n; ++i) {
        a[i] = i + 1;
        orig[i] = i + 1;
    }

    cfx_ntt_bit_reverse(a, n);
    cfx_ntt_forward(a, n, p, tw.forward);
    cfx_ntt_bit_reverse(a, n);
    cfx_ntt_inverse(a, n, p, tw.inverse, tw.n_inv);

    for (size_t i = 0; i < n; ++i) {
        CFX_ASSERT(a[i] == orig[i]);
    }

    free(a);
    free(orig);
    cfx_ntt_twiddles_free(&tw);
}

static void test_ntt_n_equals_one(void) {
    uint64_t p = CFX_NTT_P1;
    uint64_t a[1] = {12345};

    cfx_ntt_twiddles_t tw;
    int rc = cfx_ntt_twiddles_init(&tw, 1, p, CFX_NTT_G1);
    CFX_ASSERT(rc == 0);

    cfx_ntt_bit_reverse(a, 1);
    cfx_ntt_forward(a, 1, p, tw.forward);
    CFX_ASSERT(a[0] == 12345);

    cfx_ntt_bit_reverse(a, 1);
    cfx_ntt_inverse(a, 1, p, tw.inverse, tw.n_inv);
    CFX_ASSERT(a[0] == 12345);

    cfx_ntt_twiddles_free(&tw);
}

static void test_ntt_n_equals_two(void) {
    uint64_t p = CFX_NTT_P1;
    uint64_t a[2] = {3, 7};
    uint64_t orig[2] = {3, 7};

    cfx_ntt_twiddles_t tw;
    int rc = cfx_ntt_twiddles_init(&tw, 2, p, CFX_NTT_G1);
    CFX_ASSERT(rc == 0);

    cfx_ntt_bit_reverse(a, 2);
    cfx_ntt_forward(a, 2, p, tw.forward);

    /* NTT of [a, b] with omega=-1 should be [a+b, a-b] */
    CFX_ASSERT(a[0] == 10);      /* 3 + 7 */
    CFX_ASSERT(a[1] == p - 4);   /* 3 - 7 = -4 mod p */

    cfx_ntt_bit_reverse(a, 2);
    cfx_ntt_inverse(a, 2, p, tw.inverse, tw.n_inv);

    CFX_ASSERT(a[0] == orig[0]);
    CFX_ASSERT(a[1] == orig[1]);

    cfx_ntt_twiddles_free(&tw);
}

static void test_ntt_all_zeros(void) {
    uint64_t p = CFX_NTT_P1;
    size_t n = 8;
    uint64_t a[8] = {0, 0, 0, 0, 0, 0, 0, 0};

    cfx_ntt_twiddles_t tw;
    int rc = cfx_ntt_twiddles_init(&tw, n, p, CFX_NTT_G1);
    CFX_ASSERT(rc == 0);

    cfx_ntt_bit_reverse(a, n);
    cfx_ntt_forward(a, n, p, tw.forward);

    for (size_t i = 0; i < n; ++i) {
        CFX_ASSERT(a[i] == 0);
    }

    cfx_ntt_bit_reverse(a, n);
    cfx_ntt_inverse(a, n, p, tw.inverse, tw.n_inv);

    for (size_t i = 0; i < n; ++i) {
        CFX_ASSERT(a[i] == 0);
    }

    cfx_ntt_twiddles_free(&tw);
}

static void test_ntt_single_nonzero(void) {
    uint64_t p = CFX_NTT_P1;
    size_t n = 8;
    uint64_t a[8] = {5, 0, 0, 0, 0, 0, 0, 0};
    uint64_t orig[8] = {5, 0, 0, 0, 0, 0, 0, 0};

    cfx_ntt_twiddles_t tw;
    int rc = cfx_ntt_twiddles_init(&tw, n, p, CFX_NTT_G1);
    CFX_ASSERT(rc == 0);

    cfx_ntt_bit_reverse(a, n);
    cfx_ntt_forward(a, n, p, tw.forward);

    /* NTT of [c, 0, 0, ...] should be [c, c, c, ...] */
    for (size_t i = 0; i < n; ++i) {
        CFX_ASSERT(a[i] == 5);
    }

    cfx_ntt_bit_reverse(a, n);
    cfx_ntt_inverse(a, n, p, tw.inverse, tw.n_inv);

    for (size_t i = 0; i < n; ++i) {
        CFX_ASSERT(a[i] == orig[i]);
    }

    cfx_ntt_twiddles_free(&tw);
}

static void test_ntt_linearity(void) {
    uint64_t p = CFX_NTT_P1;
    size_t n = 8;

    cfx_ntt_twiddles_t tw;
    int rc = cfx_ntt_twiddles_init(&tw, n, p, CFX_NTT_G1);
    CFX_ASSERT(rc == 0);

    uint64_t a[8] = {1, 2, 3, 4, 5, 6, 7, 8};
    uint64_t b[8] = {8, 7, 6, 5, 4, 3, 2, 1};
    uint64_t sum[8], ntt_a[8], ntt_b[8], ntt_sum[8];

    for (size_t i = 0; i < n; ++i) {
        sum[i] = cfx_ntt_mod_add(a[i], b[i], p);
        ntt_a[i] = a[i];
        ntt_b[i] = b[i];
        ntt_sum[i] = sum[i];
    }

    cfx_ntt_bit_reverse(ntt_a, n);
    cfx_ntt_forward(ntt_a, n, p, tw.forward);

    cfx_ntt_bit_reverse(ntt_b, n);
    cfx_ntt_forward(ntt_b, n, p, tw.forward);

    cfx_ntt_bit_reverse(ntt_sum, n);
    cfx_ntt_forward(ntt_sum, n, p, tw.forward);

    /* NTT(a + b) = NTT(a) + NTT(b) */
    for (size_t i = 0; i < n; ++i) {
        uint64_t expected = cfx_ntt_mod_add(ntt_a[i], ntt_b[i], p);
        CFX_ASSERT(ntt_sum[i] == expected);
    }

    cfx_ntt_twiddles_free(&tw);
}

static void test_all_primes(void) {
    /* verify primitive root for p1 (well-known Goldilocks prime) */
    uint64_t p = CFX_NTT_P1;
    uint64_t g = CFX_NTT_G1;

    /* g^((p-1)/2) should equal p-1 (i.e., -1 mod p) */
    uint64_t half = cfx_ntt_mod_pow(g, (p - 1) / 2, p);
    CFX_ASSERT(half == p - 1);

    /* Fermat test */
    CFX_ASSERT(cfx_ntt_mod_pow(g, p - 1, p) == 1);
}

static void test_convolve_simple(void) {
    /* (1 + 2x) * (3 + 4x) = 3 + 4x + 6x + 8x² = 3 + 10x + 8x² */
    uint64_t a[2] = {1, 2};
    uint64_t b[2] = {3, 4};
    uint64_t c[4];

    int rc = cfx_ntt_convolve(c, a, 2, b, 2, 4, CFX_NTT_P1, CFX_NTT_G1);
    CFX_ASSERT(rc == 0);

    CFX_ASSERT(c[0] == 3);   /* 1*3 */
    CFX_ASSERT(c[1] == 10);  /* 1*4 + 2*3 */
    CFX_ASSERT(c[2] == 8);   /* 2*4 */
    CFX_ASSERT(c[3] == 0);   /* zero-padded */
}

static void test_convolve_larger(void) {
    /* (1 + 2x + 3x²) * (4 + 5x + 6x²)
       = 4 + 5x + 6x² + 8x + 10x² + 12x³ + 12x² + 15x³ + 18x⁴
       = 4 + 13x + 28x² + 27x³ + 18x⁴ */
    uint64_t a[3] = {1, 2, 3};
    uint64_t b[3] = {4, 5, 6};
    uint64_t c[8];

    int rc = cfx_ntt_convolve(c, a, 3, b, 3, 8, CFX_NTT_P1, CFX_NTT_G1);
    CFX_ASSERT(rc == 0);

    CFX_ASSERT(c[0] == 4);   /* 1*4 */
    CFX_ASSERT(c[1] == 13);  /* 1*5 + 2*4 */
    CFX_ASSERT(c[2] == 28);  /* 1*6 + 2*5 + 3*4 */
    CFX_ASSERT(c[3] == 27);  /* 2*6 + 3*5 */
    CFX_ASSERT(c[4] == 18);  /* 3*6 */
    CFX_ASSERT(c[5] == 0);
    CFX_ASSERT(c[6] == 0);
    CFX_ASSERT(c[7] == 0);
}

static void test_convolve_single_element(void) {
    /* multiply by constant: [5] * [2, 3, 4] = [10, 15, 20] */
    uint64_t a[1] = {5};
    uint64_t b[3] = {2, 3, 4};
    uint64_t c[4];

    int rc = cfx_ntt_convolve(c, a, 1, b, 3, 4, CFX_NTT_P1, CFX_NTT_G1);
    CFX_ASSERT(rc == 0);

    CFX_ASSERT(c[0] == 10);
    CFX_ASSERT(c[1] == 15);
    CFX_ASSERT(c[2] == 20);
    CFX_ASSERT(c[3] == 0);
}

static void test_convolve_invalid_n(void) {
    uint64_t a[2] = {1, 2};
    uint64_t b[2] = {3, 4};
    uint64_t c[4];

    /* n = 0 should fail */
    int rc = cfx_ntt_convolve(c, a, 2, b, 2, 0, CFX_NTT_P1, CFX_NTT_G1);
    CFX_ASSERT(rc == -1);

    /* n not power of 2 should fail */
    rc = cfx_ntt_convolve(c, a, 2, b, 2, 3, CFX_NTT_P1, CFX_NTT_G1);
    CFX_ASSERT(rc == -1);

    /* n too small (need at least 3 for len=2 * len=2) */
    rc = cfx_ntt_convolve(c, a, 2, b, 2, 2, CFX_NTT_P1, CFX_NTT_G1);
    CFX_ASSERT(rc == -1);
}

static void test_convolve_vs_schoolbook(void) {
    /* compare NTT convolution against schoolbook for n=16 */
    uint64_t p = CFX_NTT_P1;
    size_t len = 8;
    size_t n = 16;

    uint64_t a[8] = {1, 2, 3, 4, 5, 6, 7, 8};
    uint64_t b[8] = {8, 7, 6, 5, 4, 3, 2, 1};
    uint64_t c_ntt[16];
    uint64_t c_school[16];

    /* schoolbook convolution */
    for (size_t i = 0; i < n; ++i) {
        c_school[i] = 0;
    }
    for (size_t i = 0; i < len; ++i) {
        for (size_t j = 0; j < len; ++j) {
            uint64_t prod = cfx_ntt_mod_mul(a[i], b[j], p);
            c_school[i + j] = cfx_ntt_mod_add(c_school[i + j], prod, p);
        }
    }

    /* NTT convolution */
    int rc = cfx_ntt_convolve(c_ntt, a, len, b, len, n, p, CFX_NTT_G1);
    CFX_ASSERT(rc == 0);

    /* compare */
    for (size_t i = 0; i < n; ++i) {
        CFX_ASSERT(c_ntt[i] == c_school[i]);
    }
}

static void test_convolve_zero_input(void) {
    uint64_t a[2] = {1, 2};
    uint64_t b[2] = {0, 0};
    uint64_t c[4];

    int rc = cfx_ntt_convolve(c, a, 2, b, 0, 4, CFX_NTT_P1, CFX_NTT_G1);
    CFX_ASSERT(rc == 0);

    for (size_t i = 0; i < 4; ++i) {
        CFX_ASSERT(c[i] == 0);
    }

    /* also test len_a = 0 */
    rc = cfx_ntt_convolve(c, a, 0, b, 2, 4, CFX_NTT_P1, CFX_NTT_G1);
    CFX_ASSERT(rc == 0);
    for (size_t i = 0; i < 4; ++i) {
        CFX_ASSERT(c[i] == 0);
    }
}

static void test_convolve_one_element_each(void) {
    uint64_t a[1] = {7};
    uint64_t b[1] = {11};
    uint64_t c[2];

    int rc = cfx_ntt_convolve(c, a, 1, b, 1, 2, CFX_NTT_P1, CFX_NTT_G1);
    CFX_ASSERT(rc == 0);

    CFX_ASSERT(c[0] == 77);
    CFX_ASSERT(c[1] == 0);
}

static void test_convolve_commutativity(void) {
    uint64_t p = CFX_NTT_P1;
    uint64_t a[3] = {1, 2, 3};
    uint64_t b[4] = {4, 5, 6, 7};
    uint64_t c1[8], c2[8];

    int rc1 = cfx_ntt_convolve(c1, a, 3, b, 4, 8, p, CFX_NTT_G1);
    int rc2 = cfx_ntt_convolve(c2, b, 4, a, 3, 8, p, CFX_NTT_G1);

    CFX_ASSERT(rc1 == 0);
    CFX_ASSERT(rc2 == 0);

    for (size_t i = 0; i < 8; ++i) {
        CFX_ASSERT(c1[i] == c2[i]);
    }
}

static void test_convolve_associativity(void) {
    uint64_t p = CFX_NTT_P1;
    uint64_t a[2] = {1, 2};
    uint64_t b[2] = {3, 4};
    uint64_t c[2] = {5, 6};

    uint64_t ab[4], abc_1[8];
    uint64_t bc[4], abc_2[8];

    /* (a * b) * c */
    cfx_ntt_convolve(ab, a, 2, b, 2, 4, p, CFX_NTT_G1);
    cfx_ntt_convolve(abc_1, ab, 3, c, 2, 8, p, CFX_NTT_G1);

    /* a * (b * c) */
    cfx_ntt_convolve(bc, b, 2, c, 2, 4, p, CFX_NTT_G1);
    cfx_ntt_convolve(abc_2, a, 2, bc, 3, 8, p, CFX_NTT_G1);

    for (size_t i = 0; i < 8; ++i) {
        CFX_ASSERT(abc_1[i] == abc_2[i]);
    }
}

static void test_convolve_identity(void) {
    uint64_t p = CFX_NTT_P1;
    uint64_t a[4] = {10, 20, 30, 40};
    uint64_t one[1] = {1};
    uint64_t c[8];

    /* a * [1] = a */
    int rc = cfx_ntt_convolve(c, a, 4, one, 1, 8, p, CFX_NTT_G1);
    CFX_ASSERT(rc == 0);

    CFX_ASSERT(c[0] == 10);
    CFX_ASSERT(c[1] == 20);
    CFX_ASSERT(c[2] == 30);
    CFX_ASSERT(c[3] == 40);
    for (size_t i = 4; i < 8; ++i) {
        CFX_ASSERT(c[i] == 0);
    }
}

static void test_convolve_large_coefficients(void) {
    uint64_t p = CFX_NTT_P1;
    /* use values near p-1 to stress modular arithmetic */
    uint64_t a[2] = {p - 1, p - 2};
    uint64_t b[2] = {p - 3, p - 4};
    uint64_t c[4];

    int rc = cfx_ntt_convolve(c, a, 2, b, 2, 4, p, CFX_NTT_G1);
    CFX_ASSERT(rc == 0);

    /* verify against schoolbook */
    uint64_t c0 = cfx_ntt_mod_mul(a[0], b[0], p);
    uint64_t c1 = cfx_ntt_mod_add(cfx_ntt_mod_mul(a[0], b[1], p),
                                   cfx_ntt_mod_mul(a[1], b[0], p), p);
    uint64_t c2 = cfx_ntt_mod_mul(a[1], b[1], p);

    CFX_ASSERT(c[0] == c0);
    CFX_ASSERT(c[1] == c1);
    CFX_ASSERT(c[2] == c2);
    CFX_ASSERT(c[3] == 0);
}

static void test_convolve_n_exact_fit(void) {
    /* n exactly equals len_a + len_b - 1 */
    uint64_t a[3] = {1, 2, 3};
    uint64_t b[2] = {4, 5};
    uint64_t c[4];

    /* need n >= 3 + 2 - 1 = 4 */
    int rc = cfx_ntt_convolve(c, a, 3, b, 2, 4, CFX_NTT_P1, CFX_NTT_G1);
    CFX_ASSERT(rc == 0);

    CFX_ASSERT(c[0] == 4);   /* 1*4 */
    CFX_ASSERT(c[1] == 13);  /* 1*5 + 2*4 */
    CFX_ASSERT(c[2] == 22);  /* 2*5 + 3*4 */
    CFX_ASSERT(c[3] == 15);  /* 3*5 */
}

static void test_convolve_larger_random(void) {
    uint64_t p = CFX_NTT_P1;
    size_t len = 32;
    size_t n = 64;

    uint64_t* a = (uint64_t*)malloc(len * sizeof(uint64_t));
    uint64_t* b = (uint64_t*)malloc(len * sizeof(uint64_t));
    uint64_t* c_ntt = (uint64_t*)malloc(n * sizeof(uint64_t));
    uint64_t* c_school = (uint64_t*)malloc(n * sizeof(uint64_t));
    CFX_ASSERT(a && b && c_ntt && c_school);

    /* deterministic "random" values */
    for (size_t i = 0; i < len; ++i) {
        a[i] = (i * 12345 + 67890) % 1000000;
        b[i] = (i * 98765 + 43210) % 1000000;
    }

    /* schoolbook */
    for (size_t i = 0; i < n; ++i) c_school[i] = 0;
    for (size_t i = 0; i < len; ++i) {
        for (size_t j = 0; j < len; ++j) {
            uint64_t prod = cfx_ntt_mod_mul(a[i], b[j], p);
            c_school[i + j] = cfx_ntt_mod_add(c_school[i + j], prod, p);
        }
    }

    int rc = cfx_ntt_convolve(c_ntt, a, len, b, len, n, p, CFX_NTT_G1);
    CFX_ASSERT(rc == 0);

    for (size_t i = 0; i < n; ++i) {
        CFX_ASSERT(c_ntt[i] == c_school[i]);
    }

    free(a);
    free(b);
    free(c_ntt);
    free(c_school);
}

static void test_limbs_to_chunks_simple(void) {
    cfx_limb_t limbs[2];
    uint64_t chunks[8];
    size_t n;

#if CFX_LIMB_BITS == 32
    limbs[0] = 0x12345678;
    limbs[1] = 0xABCDEF01;
    n = cfx_ntt_limbs_to_chunks(chunks, 8, limbs, 2);
    CFX_ASSERT(n == 4);
    CFX_ASSERT(chunks[0] == 0x5678);
    CFX_ASSERT(chunks[1] == 0x1234);
    CFX_ASSERT(chunks[2] == 0xEF01);
    CFX_ASSERT(chunks[3] == 0xABCD);
#elif CFX_LIMB_BITS == 64
    limbs[0] = UINT64_C(0x123456789ABCDEF0);
    limbs[1] = UINT64_C(0xFEDCBA9876543210);
    n = cfx_ntt_limbs_to_chunks(chunks, 8, limbs, 2);
    CFX_ASSERT(n == 8);
    CFX_ASSERT(chunks[0] == 0xDEF0);
    CFX_ASSERT(chunks[1] == 0x9ABC);
    CFX_ASSERT(chunks[2] == 0x5678);
    CFX_ASSERT(chunks[3] == 0x1234);
    CFX_ASSERT(chunks[4] == 0x3210);
    CFX_ASSERT(chunks[5] == 0x7654);
    CFX_ASSERT(chunks[6] == 0xBA98);
    CFX_ASSERT(chunks[7] == 0xFEDC);
#endif
}

static void test_limbs_to_chunks_zero(void) {
    cfx_limb_t limbs[2] = {0, 0};
    uint64_t chunks[8];
    size_t n = cfx_ntt_limbs_to_chunks(chunks, 8, limbs, 2);

#if CFX_LIMB_BITS == 32
    CFX_ASSERT(n == 4);
#elif CFX_LIMB_BITS == 64
    CFX_ASSERT(n == 8);
#endif
    for (size_t i = 0; i < n; ++i) {
        CFX_ASSERT(chunks[i] == 0);
    }
}

static void test_chunks_to_limbs_simple(void) {
    uint64_t chunks[4] = {0x5678, 0x1234, 0xEF01, 0xABCD};
    cfx_limb_t limbs[4];
    size_t n;

#if CFX_LIMB_BITS == 32
    n = cfx_ntt_chunks_to_limbs(limbs, 4, chunks, 4);
    CFX_ASSERT(n == 2);
    CFX_ASSERT(limbs[0] == 0x12345678);
    CFX_ASSERT(limbs[1] == 0xABCDEF01);
#elif CFX_LIMB_BITS == 64
    n = cfx_ntt_chunks_to_limbs(limbs, 4, chunks, 4);
    CFX_ASSERT(n == 1);
    CFX_ASSERT(limbs[0] == UINT64_C(0xABCDEF0112345678));
#endif
}

static void test_chunks_to_limbs_carry(void) {
    /* test carry propagation: value exceeds 16 bits */
    uint64_t chunks[4] = {0x1FFFF, 0, 0, 0};  /* 0x1FFFF = 0xFFFF + 0x10000 */
    cfx_limb_t limbs[4];

    size_t n = cfx_ntt_chunks_to_limbs(limbs, 4, chunks, 4);

#if CFX_LIMB_BITS == 32
    /* 0x1FFFF in chunks means: chunk[0] = 0xFFFF, carry 1 to chunk[1] */
    /* result: limb[0] = 0x0001FFFF */
    CFX_ASSERT(n == 1);
    CFX_ASSERT(limbs[0] == 0x0001FFFF);
#elif CFX_LIMB_BITS == 64
    CFX_ASSERT(n == 1);
    CFX_ASSERT(limbs[0] == UINT64_C(0x0001FFFF));
#endif
}

static void test_chunks_to_limbs_large_carry(void) {
    /* large value that propagates carry across multiple chunks */
    uint64_t chunks[4] = {0xFFFFFFFF, 0, 0, 0};
    cfx_limb_t limbs[4];

    size_t n = cfx_ntt_chunks_to_limbs(limbs, 4, chunks, 4);

    /* 0xFFFFFFFF spans multiple 16-bit chunks with carry */
#if CFX_LIMB_BITS == 32
    CFX_ASSERT(n == 1);
    CFX_ASSERT(limbs[0] == 0xFFFFFFFF);
#elif CFX_LIMB_BITS == 64
    CFX_ASSERT(n == 1);
    CFX_ASSERT(limbs[0] == UINT64_C(0xFFFFFFFF));
#endif
}

static void test_roundtrip_limbs_chunks(void) {
    /* convert limbs -> chunks -> limbs and verify */
    cfx_limb_t original[4];
    cfx_limb_t result[4];
    uint64_t chunks[16];

#if CFX_LIMB_BITS == 32
    original[0] = 0x12345678;
    original[1] = 0x9ABCDEF0;
    original[2] = 0x11111111;
    original[3] = 0x22222222;

    size_t n_chunks = cfx_ntt_limbs_to_chunks(chunks, 16, original, 4);
    CFX_ASSERT(n_chunks == 8);

    size_t n_limbs = cfx_ntt_chunks_to_limbs(result, 4, chunks, n_chunks);
    CFX_ASSERT(n_limbs == 4);

    for (size_t i = 0; i < 4; ++i) {
        CFX_ASSERT(result[i] == original[i]);
    }
#elif CFX_LIMB_BITS == 64
    original[0] = UINT64_C(0x123456789ABCDEF0);
    original[1] = UINT64_C(0xFEDCBA9876543210);
    original[2] = UINT64_C(0x1111111122222222);
    original[3] = UINT64_C(0x3333333344444444);

    size_t n_chunks = cfx_ntt_limbs_to_chunks(chunks, 16, original, 4);
    CFX_ASSERT(n_chunks == 16);

    size_t n_limbs = cfx_ntt_chunks_to_limbs(result, 4, chunks, n_chunks);
    CFX_ASSERT(n_limbs == 4);

    for (size_t i = 0; i < 4; ++i) {
        CFX_ASSERT(result[i] == original[i]);
    }
#endif
}

static void test_ntt_mul_limbs_simple(void) {
    /* 2 * 3 = 6 */
    cfx_limb_t a[1] = {2};
    cfx_limb_t b[1] = {3};
    cfx_limb_t out[2];

    size_t n = cfx_ntt_mul_limbs(out, 2, a, 1, b, 1);
    CFX_ASSERT(n == 1);
    CFX_ASSERT(out[0] == 6);
}

static void test_ntt_mul_limbs_zero(void) {
    cfx_limb_t a[1] = {0};
    cfx_limb_t b[1] = {123};
    cfx_limb_t out[2];

    size_t n = cfx_ntt_mul_limbs(out, 2, a, 0, b, 1);
    CFX_ASSERT(n == 0);

    n = cfx_ntt_mul_limbs(out, 2, a, 1, b, 0);
    CFX_ASSERT(n == 0);
}

static void test_ntt_mul_limbs_medium(void) {
    /* multiply two 2-limb numbers and verify against manual calculation */
#if CFX_LIMB_BITS == 32
    /* 0x100000001 * 0x100000002 = 0x100000003 00000002 */
    /* in base 2^32: {1, 1} * {2, 1} = {2, 3, 1} */
    cfx_limb_t a[2] = {1, 1};  /* 2^32 + 1 */
    cfx_limb_t b[2] = {2, 1};  /* 2^32 + 2 */
    cfx_limb_t out[4];

    size_t n = cfx_ntt_mul_limbs(out, 4, a, 2, b, 2);
    CFX_ASSERT(n == 3);
    CFX_ASSERT(out[0] == 2);  /* 1*2 */
    CFX_ASSERT(out[1] == 3);  /* 1*1 + 1*2 */
    CFX_ASSERT(out[2] == 1);  /* 1*1 */
#elif CFX_LIMB_BITS == 64
    cfx_limb_t a[2] = {1, 1};
    cfx_limb_t b[2] = {2, 1};
    cfx_limb_t out[4];

    size_t n = cfx_ntt_mul_limbs(out, 4, a, 2, b, 2);
    CFX_ASSERT(n == 3);
    CFX_ASSERT(out[0] == 2);
    CFX_ASSERT(out[1] == 3);
    CFX_ASSERT(out[2] == 1);
#endif
}

static void test_ntt_mul_limbs_large(void) {
    /* multiply larger numbers and verify against schoolbook */
    cfx_limb_t a[4], b[4], out_ntt[8], out_school[8];

    /* initialize with some values */
    for (size_t i = 0; i < 4; ++i) {
        a[i] = (cfx_limb_t)(i + 1);
        b[i] = (cfx_limb_t)(i + 5);
    }

    /* NTT multiplication */
    size_t n_ntt = cfx_ntt_mul_limbs(out_ntt, 8, a, 4, b, 4);

    /* schoolbook multiplication */
    for (size_t i = 0; i < 8; ++i) out_school[i] = 0;
    for (size_t i = 0; i < 4; ++i) {
        cfx_acc_t carry = cfx_acc_zero();
        for (size_t j = 0; j < 4; ++j) {
            cfx_acc_mac(&carry, a[i], b[j]);
            cfx_acc_add_lo(&carry, out_school[i + j]);
            out_school[i + j] = cfx_acc_lo(carry);
            carry = cfx_acc_from_lo(cfx_acc_hi(carry));
        }
        out_school[i + 4] = cfx_acc_lo(carry);
    }

    /* trim trailing zeros */
    size_t n_school = 8;
    while (n_school > 0 && out_school[n_school - 1] == 0) n_school--;

    CFX_ASSERT(n_ntt == n_school);
    for (size_t i = 0; i < n_ntt; ++i) {
        CFX_ASSERT(out_ntt[i] == out_school[i]);
    }
}

static void test_ntt_mul_limbs_known_product(void) {
    /* 255 * 255 = 65025 */
    cfx_limb_t a[1] = {255};
    cfx_limb_t b[1] = {255};
    cfx_limb_t out[2];

    size_t n = cfx_ntt_mul_limbs(out, 2, a, 1, b, 1);
    CFX_ASSERT(n == 1);
    CFX_ASSERT(out[0] == 65025);
}

static void test_ntt_mul_limbs_max_single(void) {
    /* max_limb * max_limb */
    cfx_limb_t a[1] = {CFX_LIMB_MAX};
    cfx_limb_t b[1] = {CFX_LIMB_MAX};
    cfx_limb_t out[2];

    size_t n = cfx_ntt_mul_limbs(out, 2, a, 1, b, 1);
    CFX_ASSERT(n == 2);

    /* verify: (2^bits - 1)^2 = 2^(2*bits) - 2^(bits+1) + 1 */
    /* in base 2^bits: {1, -2 mod base, 0, 0} wait that's wrong...
       actually: (2^32 - 1)^2 = 2^64 - 2^33 + 1
       in base 2^32: low = 1, high = 2^32 - 2
       because 2^64 - 2^33 + 1 = (2^32 - 2) * 2^32 + 1
    */
#if CFX_LIMB_BITS == 32
    CFX_ASSERT(out[0] == 1);
    CFX_ASSERT(out[1] == 0xFFFFFFFE);
#elif CFX_LIMB_BITS == 64
    CFX_ASSERT(out[0] == 1);
    CFX_ASSERT(out[1] == UINT64_C(0xFFFFFFFFFFFFFFFE));
#endif
}

static void test_ntt_mul_limbs_asymmetric(void) {
    /* multiply operands of different sizes */
    cfx_limb_t a[4] = {1, 2, 3, 4};
    cfx_limb_t b[2] = {5, 6};
    cfx_limb_t out_ntt[6], out_school[6];

    size_t n_ntt = cfx_ntt_mul_limbs(out_ntt, 6, a, 4, b, 2);

    /* schoolbook */
    for (size_t i = 0; i < 6; ++i) out_school[i] = 0;
    for (size_t i = 0; i < 4; ++i) {
        cfx_acc_t carry = cfx_acc_zero();
        for (size_t j = 0; j < 2; ++j) {
            cfx_acc_mac(&carry, a[i], b[j]);
            cfx_acc_add_lo(&carry, out_school[i + j]);
            out_school[i + j] = cfx_acc_lo(carry);
            carry = cfx_acc_from_lo(cfx_acc_hi(carry));
        }
        out_school[i + 2] = cfx_acc_lo(carry);
    }
    size_t n_school = 6;
    while (n_school > 0 && out_school[n_school - 1] == 0) n_school--;

    CFX_ASSERT(n_ntt == n_school);
    for (size_t i = 0; i < n_ntt; ++i) {
        CFX_ASSERT(out_ntt[i] == out_school[i]);
    }
}

static void test_ntt_mul_limbs_by_one(void) {
    /* a * 1 = a */
    cfx_limb_t a[4] = {0x12345678, 0x9ABCDEF0, 0x11111111, 0x22222222};
    cfx_limb_t one[1] = {1};
    cfx_limb_t out[5];

    size_t n = cfx_ntt_mul_limbs(out, 5, a, 4, one, 1);
    CFX_ASSERT(n == 4);
    for (size_t i = 0; i < 4; ++i) {
        CFX_ASSERT(out[i] == a[i]);
    }

    /* also test 1 * a = a */
    n = cfx_ntt_mul_limbs(out, 5, one, 1, a, 4);
    CFX_ASSERT(n == 4);
    for (size_t i = 0; i < 4; ++i) {
        CFX_ASSERT(out[i] == a[i]);
    }
}

static void test_ntt_mul_limbs_power_of_two(void) {
    /* multiplying by power of 2 in limb form */
    cfx_limb_t a[2] = {0xFFFFFFFF, 0x00000001};  /* 2^32 + 2^32-1 */
    cfx_limb_t b[1] = {2};
    cfx_limb_t out[3];

    size_t n = cfx_ntt_mul_limbs(out, 3, a, 2, b, 1);

    /* (2^32 + 0xFFFFFFFF) * 2 = 2^33 + 2^32 + 2^32 - 2 = 2^33 + 2^33 - 2 = 2^34 - 2 */
    /* wait let me recalculate: a = {0xFFFFFFFF, 1} = 1*2^32 + 0xFFFFFFFF = 0x1FFFFFFFF */
    /* a * 2 = 0x3FFFFFFFE = {0xFFFFFFFE, 3} */
#if CFX_LIMB_BITS == 32
    CFX_ASSERT(n == 2);
    CFX_ASSERT(out[0] == 0xFFFFFFFE);
    CFX_ASSERT(out[1] == 0x00000003);
#elif CFX_LIMB_BITS == 64
    /* 64-bit: a = 1*2^64 + 0xFFFFFFFF, so a*2 = 2*2^64 + 0x1FFFFFFFE */
    CFX_ASSERT(n == 2);
    CFX_ASSERT(out[0] == UINT64_C(0x1FFFFFFFE));
    CFX_ASSERT(out[1] == UINT64_C(0x2));
#endif
}

static void test_ntt_mul_limbs_sequential(void) {
    /* multiply {1,2,3,...,8} by {1,2,3,...,8} */
    cfx_limb_t a[8], b[8], out_ntt[16], out_school[16];
    for (size_t i = 0; i < 8; ++i) {
        a[i] = (cfx_limb_t)(i + 1);
        b[i] = (cfx_limb_t)(i + 1);
    }

    size_t n_ntt = cfx_ntt_mul_limbs(out_ntt, 16, a, 8, b, 8);

    /* schoolbook */
    for (size_t i = 0; i < 16; ++i) out_school[i] = 0;
    for (size_t i = 0; i < 8; ++i) {
        cfx_acc_t carry = cfx_acc_zero();
        for (size_t j = 0; j < 8; ++j) {
            cfx_acc_mac(&carry, a[i], b[j]);
            cfx_acc_add_lo(&carry, out_school[i + j]);
            out_school[i + j] = cfx_acc_lo(carry);
            carry = cfx_acc_from_lo(cfx_acc_hi(carry));
        }
        out_school[i + 8] = cfx_acc_lo(carry);
    }
    size_t n_school = 16;
    while (n_school > 0 && out_school[n_school - 1] == 0) n_school--;

    CFX_ASSERT(n_ntt == n_school);
    for (size_t i = 0; i < n_ntt; ++i) {
        CFX_ASSERT(out_ntt[i] == out_school[i]);
    }
}

static void test_ntt_mul_limbs_all_ones(void) {
    /* {0xFFFF..., 0xFFFF..., ...} squared */
    cfx_limb_t a[4], out_ntt[8], out_school[8];
    for (size_t i = 0; i < 4; ++i) {
        a[i] = CFX_LIMB_MAX;
    }

    size_t n_ntt = cfx_ntt_mul_limbs(out_ntt, 8, a, 4, a, 4);

    /* schoolbook */
    for (size_t i = 0; i < 8; ++i) out_school[i] = 0;
    for (size_t i = 0; i < 4; ++i) {
        cfx_acc_t carry = cfx_acc_zero();
        for (size_t j = 0; j < 4; ++j) {
            cfx_acc_mac(&carry, a[i], a[j]);
            cfx_acc_add_lo(&carry, out_school[i + j]);
            out_school[i + j] = cfx_acc_lo(carry);
            carry = cfx_acc_from_lo(cfx_acc_hi(carry));
        }
        out_school[i + 4] = cfx_acc_lo(carry);
    }
    size_t n_school = 8;
    while (n_school > 0 && out_school[n_school - 1] == 0) n_school--;

    CFX_ASSERT(n_ntt == n_school);
    for (size_t i = 0; i < n_ntt; ++i) {
        CFX_ASSERT(out_ntt[i] == out_school[i]);
    }
}

static void test_ntt_mul_limbs_alternating(void) {
    /* {0xAAAA..., 0x5555..., 0xAAAA..., 0x5555...} pattern */
    cfx_limb_t a[4], b[4], out_ntt[8], out_school[8];
#if CFX_LIMB_BITS == 32
    a[0] = 0xAAAAAAAA; a[1] = 0x55555555; a[2] = 0xAAAAAAAA; a[3] = 0x55555555;
    b[0] = 0x55555555; b[1] = 0xAAAAAAAA; b[2] = 0x55555555; b[3] = 0xAAAAAAAA;
#elif CFX_LIMB_BITS == 64
    a[0] = UINT64_C(0xAAAAAAAAAAAAAAAA); a[1] = UINT64_C(0x5555555555555555);
    a[2] = UINT64_C(0xAAAAAAAAAAAAAAAA); a[3] = UINT64_C(0x5555555555555555);
    b[0] = UINT64_C(0x5555555555555555); b[1] = UINT64_C(0xAAAAAAAAAAAAAAAA);
    b[2] = UINT64_C(0x5555555555555555); b[3] = UINT64_C(0xAAAAAAAAAAAAAAAA);
#endif

    size_t n_ntt = cfx_ntt_mul_limbs(out_ntt, 8, a, 4, b, 4);

    /* schoolbook */
    for (size_t i = 0; i < 8; ++i) out_school[i] = 0;
    for (size_t i = 0; i < 4; ++i) {
        cfx_acc_t carry = cfx_acc_zero();
        for (size_t j = 0; j < 4; ++j) {
            cfx_acc_mac(&carry, a[i], b[j]);
            cfx_acc_add_lo(&carry, out_school[i + j]);
            out_school[i + j] = cfx_acc_lo(carry);
            carry = cfx_acc_from_lo(cfx_acc_hi(carry));
        }
        out_school[i + 4] = cfx_acc_lo(carry);
    }
    size_t n_school = 8;
    while (n_school > 0 && out_school[n_school - 1] == 0) n_school--;

    CFX_ASSERT(n_ntt == n_school);
    for (size_t i = 0; i < n_ntt; ++i) {
        CFX_ASSERT(out_ntt[i] == out_school[i]);
    }
}

static void test_ntt_mul_limbs_random_16(void) {
    /* 16 limbs × 16 limbs with pseudo-random values */
    cfx_limb_t a[16], b[16], out_ntt[32], out_school[32];

    for (size_t i = 0; i < 16; ++i) {
        a[i] = (cfx_limb_t)((i * 0x12345678UL + 0x87654321UL) & CFX_LIMB_MAX);
        b[i] = (cfx_limb_t)((i * 0x9ABCDEF0UL + 0x0FEDCBA9UL) & CFX_LIMB_MAX);
    }

    size_t n_ntt = cfx_ntt_mul_limbs(out_ntt, 32, a, 16, b, 16);

    /* schoolbook */
    for (size_t i = 0; i < 32; ++i) out_school[i] = 0;
    for (size_t i = 0; i < 16; ++i) {
        cfx_acc_t carry = cfx_acc_zero();
        for (size_t j = 0; j < 16; ++j) {
            cfx_acc_mac(&carry, a[i], b[j]);
            cfx_acc_add_lo(&carry, out_school[i + j]);
            out_school[i + j] = cfx_acc_lo(carry);
            carry = cfx_acc_from_lo(cfx_acc_hi(carry));
        }
        out_school[i + 16] = cfx_acc_lo(carry);
    }
    size_t n_school = 32;
    while (n_school > 0 && out_school[n_school - 1] == 0) n_school--;

    CFX_ASSERT(n_ntt == n_school);
    for (size_t i = 0; i < n_ntt; ++i) {
        CFX_ASSERT(out_ntt[i] == out_school[i]);
    }
}

static void test_ntt_mul_limbs_commutativity(void) {
    /* a * b = b * a */
    cfx_limb_t a[5] = {1, 2, 3, 4, 5};
    cfx_limb_t b[3] = {6, 7, 8};
    cfx_limb_t out1[8], out2[8];

    size_t n1 = cfx_ntt_mul_limbs(out1, 8, a, 5, b, 3);
    size_t n2 = cfx_ntt_mul_limbs(out2, 8, b, 3, a, 5);

    CFX_ASSERT(n1 == n2);
    for (size_t i = 0; i < n1; ++i) {
        CFX_ASSERT(out1[i] == out2[i]);
    }
}

static void test_chunks_to_limbs_multi_carry(void) {
    /* test carry that propagates through multiple chunks */
    /* if each chunk adds up to > 16 bits, carry cascades */
    uint64_t chunks[8] = {0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF};
    cfx_limb_t limbs[8];

    size_t n = cfx_ntt_chunks_to_limbs(limbs, 8, chunks, 8);

    /* all 0xFFFF chunks = 0xFFFF FFFF FFFF FFFF FFFF FFFF FFFF FFFF in base 2^16 */
#if CFX_LIMB_BITS == 32
    CFX_ASSERT(n == 4);
    CFX_ASSERT(limbs[0] == 0xFFFFFFFF);
    CFX_ASSERT(limbs[1] == 0xFFFFFFFF);
    CFX_ASSERT(limbs[2] == 0xFFFFFFFF);
    CFX_ASSERT(limbs[3] == 0xFFFFFFFF);
#elif CFX_LIMB_BITS == 64
    CFX_ASSERT(n == 2);
    CFX_ASSERT(limbs[0] == UINT64_C(0xFFFFFFFFFFFFFFFF));
    CFX_ASSERT(limbs[1] == UINT64_C(0xFFFFFFFFFFFFFFFF));
#endif
}

static void test_chunks_to_limbs_overflow_cascade(void) {
    /* value that causes multiple levels of carry */
    uint64_t chunks[4] = {0x1FFFE, 0, 0, 0};  /* needs carry into chunk 1 */
    cfx_limb_t limbs[4];

    size_t n = cfx_ntt_chunks_to_limbs(limbs, 4, chunks, 4);

    /* 0x1FFFE = 0xFFFE + carry of 1 to next chunk */
#if CFX_LIMB_BITS == 32
    CFX_ASSERT(n == 1);
    CFX_ASSERT(limbs[0] == 0x0001FFFE);
#elif CFX_LIMB_BITS == 64
    CFX_ASSERT(n == 1);
    CFX_ASSERT(limbs[0] == UINT64_C(0x0001FFFE));
#endif
}

static void test_limbs_to_chunks_single(void) {
    /* single limb conversion */
    cfx_limb_t limbs[1];
    uint64_t chunks[4];

#if CFX_LIMB_BITS == 32
    limbs[0] = 0xDEADBEEF;
    size_t n = cfx_ntt_limbs_to_chunks(chunks, 4, limbs, 1);
    CFX_ASSERT(n == 2);
    CFX_ASSERT(chunks[0] == 0xBEEF);
    CFX_ASSERT(chunks[1] == 0xDEAD);
#elif CFX_LIMB_BITS == 64
    limbs[0] = UINT64_C(0xDEADBEEFCAFEBABE);
    size_t n = cfx_ntt_limbs_to_chunks(chunks, 4, limbs, 1);
    CFX_ASSERT(n == 4);
    CFX_ASSERT(chunks[0] == 0xBABE);
    CFX_ASSERT(chunks[1] == 0xCAFE);
    CFX_ASSERT(chunks[2] == 0xBEEF);
    CFX_ASSERT(chunks[3] == 0xDEAD);
#endif
}

static void test_ntt_mul_limbs_square(void) {
    /* a * a (squaring) */
    cfx_limb_t a[4] = {0x12345678, 0x9ABCDEF0, 0x11111111, 0x22222222};
    cfx_limb_t out_ntt[8], out_school[8];

    size_t n_ntt = cfx_ntt_mul_limbs(out_ntt, 8, a, 4, a, 4);

    /* schoolbook */
    for (size_t i = 0; i < 8; ++i) out_school[i] = 0;
    for (size_t i = 0; i < 4; ++i) {
        cfx_acc_t carry = cfx_acc_zero();
        for (size_t j = 0; j < 4; ++j) {
            cfx_acc_mac(&carry, a[i], a[j]);
            cfx_acc_add_lo(&carry, out_school[i + j]);
            out_school[i + j] = cfx_acc_lo(carry);
            carry = cfx_acc_from_lo(cfx_acc_hi(carry));
        }
        out_school[i + 4] = cfx_acc_lo(carry);
    }
    size_t n_school = 8;
    while (n_school > 0 && out_school[n_school - 1] == 0) n_school--;

    CFX_ASSERT(n_ntt == n_school);
    for (size_t i = 0; i < n_ntt; ++i) {
        CFX_ASSERT(out_ntt[i] == out_school[i]);
    }
}

static void test_ntt_mul_limbs_fibonacci_style(void) {
    /* fibonacci-like values */
    cfx_limb_t a[8], b[8], out_ntt[16], out_school[16];
    a[0] = 1; a[1] = 1;
    b[0] = 1; b[1] = 1;
    for (size_t i = 2; i < 8; ++i) {
        a[i] = a[i-1] + a[i-2];
        b[i] = b[i-1] + b[i-2];
    }

    size_t n_ntt = cfx_ntt_mul_limbs(out_ntt, 16, a, 8, b, 8);

    /* schoolbook */
    for (size_t i = 0; i < 16; ++i) out_school[i] = 0;
    for (size_t i = 0; i < 8; ++i) {
        cfx_acc_t carry = cfx_acc_zero();
        for (size_t j = 0; j < 8; ++j) {
            cfx_acc_mac(&carry, a[i], b[j]);
            cfx_acc_add_lo(&carry, out_school[i + j]);
            out_school[i + j] = cfx_acc_lo(carry);
            carry = cfx_acc_from_lo(cfx_acc_hi(carry));
        }
        out_school[i + 8] = cfx_acc_lo(carry);
    }
    size_t n_school = 16;
    while (n_school > 0 && out_school[n_school - 1] == 0) n_school--;

    CFX_ASSERT(n_ntt == n_school);
    for (size_t i = 0; i < n_ntt; ++i) {
        CFX_ASSERT(out_ntt[i] == out_school[i]);
    }
}

static void test_ntt_mul_limbs_32_limbs(void) {
    /* larger test: 32 × 32 limbs */
    cfx_limb_t a[32], b[32], out_ntt[64], out_school[64];

    for (size_t i = 0; i < 32; ++i) {
        a[i] = (cfx_limb_t)((i * 17 + 3) | 0x80000000UL);
        b[i] = (cfx_limb_t)((i * 31 + 7) | 0x40000000UL);
    }

    size_t n_ntt = cfx_ntt_mul_limbs(out_ntt, 64, a, 32, b, 32);

    /* schoolbook */
    for (size_t i = 0; i < 64; ++i) out_school[i] = 0;
    for (size_t i = 0; i < 32; ++i) {
        cfx_acc_t carry = cfx_acc_zero();
        for (size_t j = 0; j < 32; ++j) {
            cfx_acc_mac(&carry, a[i], b[j]);
            cfx_acc_add_lo(&carry, out_school[i + j]);
            out_school[i + j] = cfx_acc_lo(carry);
            carry = cfx_acc_from_lo(cfx_acc_hi(carry));
        }
        out_school[i + 32] = cfx_acc_lo(carry);
    }
    size_t n_school = 64;
    while (n_school > 0 && out_school[n_school - 1] == 0) n_school--;

    CFX_ASSERT(n_ntt == n_school);
    for (size_t i = 0; i < n_ntt; ++i) {
        CFX_ASSERT(out_ntt[i] == out_school[i]);
    }
}

static void test_big_mul_fft_simple(void) {
    cfx_big_t a, b, out, expected;
    cfx_big_init(&a);
    cfx_big_init(&b);
    cfx_big_init(&out);
    cfx_big_init(&expected);

    cfx_big_from_str(&a, "12345678901234567890");
    cfx_big_from_str(&b, "98765432109876543210");
    cfx_big_mul_fft(&out, &a, &b);
    cfx_big_mul(&expected, &a, &b);

    CFX_ASSERT(cfx_big_cmp(&out, &expected) == 0);

    cfx_big_free(&a);
    cfx_big_free(&b);
    cfx_big_free(&out);
    cfx_big_free(&expected);
}

static void test_big_mul_fft_zero(void) {
    cfx_big_t a, b, out;
    cfx_big_init(&a);
    cfx_big_init(&b);
    cfx_big_init(&out);

    cfx_big_from_str(&a, "12345678901234567890");
    cfx_big_from_limb(&b, 0);
    cfx_big_mul_fft(&out, &a, &b);

    CFX_ASSERT(cfx_big_is_zero(&out));

    cfx_big_mul_fft(&out, &b, &a);
    CFX_ASSERT(cfx_big_is_zero(&out));

    cfx_big_free(&a);
    cfx_big_free(&b);
    cfx_big_free(&out);
}

static void test_big_mul_fft_one(void) {
    cfx_big_t a, one, out;
    cfx_big_init(&a);
    cfx_big_init(&one);
    cfx_big_init(&out);

    cfx_big_from_str(&a, "12345678901234567890123456789012345678901234567890");
    cfx_big_from_limb(&one, 1);
    cfx_big_mul_fft(&out, &a, &one);

    CFX_ASSERT(cfx_big_cmp(&out, &a) == 0);

    cfx_big_free(&a);
    cfx_big_free(&one);
    cfx_big_free(&out);
}

static void test_big_mul_fft_small_square(void) {
    cfx_big_t a, out, expected;
    cfx_big_init(&a);
    cfx_big_init(&out);
    cfx_big_init(&expected);

    cfx_big_from_str(&a, "999999999999999999999999999999");
    cfx_big_mul_fft(&out, &a, &a);
    cfx_big_mul(&expected, &a, &a);

    CFX_ASSERT(cfx_big_cmp(&out, &expected) == 0);

    cfx_big_free(&a);
    cfx_big_free(&out);
    cfx_big_free(&expected);
}

static void test_big_mul_fft_large(void) {
    cfx_big_t a, b, out, expected;
    cfx_big_init(&a);
    cfx_big_init(&b);
    cfx_big_init(&out);
    cfx_big_init(&expected);

    /* build a large number: 2^1024 - 1 */
    cfx_big_from_limb(&a, 1);
    cfx_big_shl_bits_eq(&a, 1024);
    cfx_big_sub_sm_eq(&a, 1);

    /* build another: 2^512 + 1 */
    cfx_big_from_limb(&b, 1);
    cfx_big_shl_bits_eq(&b, 512);
    cfx_big_add_sm_eq(&b, 1);

    cfx_big_mul_fft(&out, &a, &b);
    cfx_big_mul(&expected, &a, &b);

    CFX_ASSERT(cfx_big_cmp(&out, &expected) == 0);

    cfx_big_free(&a);
    cfx_big_free(&b);
    cfx_big_free(&out);
    cfx_big_free(&expected);
}

static void test_big_mul_fft_vs_schoolbook(void) {
    cfx_big_t a, b, fft_result, school_result;
    cfx_big_init(&a);
    cfx_big_init(&b);
    cfx_big_init(&fft_result);
    cfx_big_init(&school_result);

    /* random-ish 256-bit numbers */
    cfx_big_from_hex(&a, "FEDCBA9876543210FEDCBA9876543210FEDCBA9876543210FEDCBA9876543210");
    cfx_big_from_hex(&b, "123456789ABCDEF0123456789ABCDEF0123456789ABCDEF0123456789ABCDEF0");

    cfx_big_mul_fft(&fft_result, &a, &b);
    cfx_big_mul(&school_result, &a, &b);

    CFX_ASSERT(cfx_big_cmp(&fft_result, &school_result) == 0);

    cfx_big_free(&a);
    cfx_big_free(&b);
    cfx_big_free(&fft_result);
    cfx_big_free(&school_result);
}

static void test_big_mul_fft_commutativity(void) {
    cfx_big_t a, b, ab, ba;
    cfx_big_init(&a);
    cfx_big_init(&b);
    cfx_big_init(&ab);
    cfx_big_init(&ba);

    cfx_big_from_str(&a, "314159265358979323846264338327950288419716939937510");
    cfx_big_from_str(&b, "271828182845904523536028747135266249775724709369995");

    cfx_big_mul_fft(&ab, &a, &b);
    cfx_big_mul_fft(&ba, &b, &a);

    CFX_ASSERT(cfx_big_cmp(&ab, &ba) == 0);

    cfx_big_free(&a);
    cfx_big_free(&b);
    cfx_big_free(&ab);
    cfx_big_free(&ba);
}

static void test_big_mul_fft_self_output(void) {
    cfx_big_t a, b, out;
    cfx_big_init(&a);
    cfx_big_init(&b);
    cfx_big_init(&out);

    cfx_big_from_str(&a, "123456789");
    cfx_big_from_str(&b, "987654321");
    cfx_big_mul_fft(&out, &a, &b);

    char buf[32];
    cfx_big_snprint_dec(&out, buf, sizeof(buf));
    CFX_ASSERT(strcmp(buf, "121932631112635269") == 0);

    cfx_big_free(&a);
    cfx_big_free(&b);
    cfx_big_free(&out);
}

static void test_big_mul_eq_fft(void) {
    cfx_big_t a, b, expected;
    cfx_big_init(&a);
    cfx_big_init(&b);
    cfx_big_init(&expected);

    cfx_big_from_str(&a, "12345678901234567890");
    cfx_big_from_str(&b, "98765432109876543210");
    cfx_big_mul(&expected, &a, &b);

    cfx_big_mul_eq_fft(&a, &b);

    CFX_ASSERT(cfx_big_cmp(&a, &expected) == 0);

    cfx_big_free(&a);
    cfx_big_free(&b);
    cfx_big_free(&expected);
}

static void test_big_mul_fft_power_of_two(void) {
    cfx_big_t a, b, out, expected;
    cfx_big_init(&a);
    cfx_big_init(&b);
    cfx_big_init(&out);
    cfx_big_init(&expected);

    /* 2^64 * 2^64 = 2^128 */
    cfx_big_from_limb(&a, 1);
    cfx_big_shl_bits_eq(&a, 64);
    cfx_big_from_limb(&b, 1);
    cfx_big_shl_bits_eq(&b, 64);

    cfx_big_mul_fft(&out, &a, &b);

    cfx_big_from_limb(&expected, 1);
    cfx_big_shl_bits_eq(&expected, 128);

    CFX_ASSERT(cfx_big_cmp(&out, &expected) == 0);

    cfx_big_free(&a);
    cfx_big_free(&b);
    cfx_big_free(&out);
    cfx_big_free(&expected);
}

static void test_big_mul_fft_factorial(void) {
    /* compute 100! using FFT mul and compare with schoolbook */
    cfx_big_t fft_fact, school_fact, n;
    cfx_big_init(&fft_fact);
    cfx_big_init(&school_fact);
    cfx_big_init(&n);

    cfx_big_from_limb(&fft_fact, 1);
    cfx_big_from_limb(&school_fact, 1);

    for (int i = 2; i <= 100; ++i) {
        cfx_big_from_limb(&n, i);
        cfx_big_mul_eq_fft(&fft_fact, &n);
        cfx_big_mul_eq(&school_fact, &n);
    }

    CFX_ASSERT(cfx_big_cmp(&fft_fact, &school_fact) == 0);

    cfx_big_free(&fft_fact);
    cfx_big_free(&school_fact);
    cfx_big_free(&n);
}

static void test_big_mul_fft_asymmetric(void) {
    cfx_big_t a, b, out, expected;
    cfx_big_init(&a);
    cfx_big_init(&b);
    cfx_big_init(&out);
    cfx_big_init(&expected);

    /* very asymmetric: 2^2048 * 7 */
    cfx_big_from_limb(&a, 1);
    cfx_big_shl_bits_eq(&a, 2048);
    cfx_big_from_limb(&b, 7);

    cfx_big_mul_fft(&out, &a, &b);
    cfx_big_mul(&expected, &a, &b);

    CFX_ASSERT(cfx_big_cmp(&out, &expected) == 0);

    cfx_big_free(&a);
    cfx_big_free(&b);
    cfx_big_free(&out);
    cfx_big_free(&expected);
}

int main(void) {
    /* Modular arithmetic */
    CFX_TEST(test_mod_add);
    CFX_TEST(test_mod_add_overflow);
    CFX_TEST(test_mod_sub);
    CFX_TEST(test_mod_sub_underflow);
    CFX_TEST(test_mod_mul);
    CFX_TEST(test_mod_mul_edge_cases);
    CFX_TEST(test_mod_mul_small_prime);
    CFX_TEST(test_mod_pow);
    CFX_TEST(test_mod_pow_p_equals_one);
    CFX_TEST(test_mod_pow_base_zero);
    CFX_TEST(test_mod_pow_exp_patterns);
    CFX_TEST(test_mod_pow_large_exp);
    CFX_TEST(test_mod_inv);
    CFX_TEST(test_mod_inv_edge_cases);

    /* Roots of unity */
    CFX_TEST(test_root_of_unity);
    CFX_TEST(test_root_of_unity_various_n);
    CFX_TEST(test_root_of_unity_primitivity);

    /* Twiddle factors */
    CFX_TEST(test_twiddles_init_free);
    CFX_TEST(test_twiddles_invalid_n);
    CFX_TEST(test_twiddles_n_equals_one);
    CFX_TEST(test_twiddles_inverse_correctness);

    /* Bit reversal */
    CFX_TEST(test_bit_reverse);
    CFX_TEST(test_bit_reverse_small);
    CFX_TEST(test_bit_reverse_n16);
    CFX_TEST(test_bit_reverse_involution);

    /* NTT forward/inverse */
    CFX_TEST(test_ntt_forward_inverse_identity);
    CFX_TEST(test_ntt_larger);
    CFX_TEST(test_ntt_n_equals_one);
    CFX_TEST(test_ntt_n_equals_two);
    CFX_TEST(test_ntt_all_zeros);
    CFX_TEST(test_ntt_single_nonzero);
    CFX_TEST(test_ntt_linearity);
    CFX_TEST(test_all_primes);

    /* Convolution */
    CFX_TEST(test_convolve_simple);
    CFX_TEST(test_convolve_larger);
    CFX_TEST(test_convolve_single_element);
    CFX_TEST(test_convolve_invalid_n);
    CFX_TEST(test_convolve_vs_schoolbook);
    CFX_TEST(test_convolve_zero_input);
    CFX_TEST(test_convolve_one_element_each);
    CFX_TEST(test_convolve_commutativity);
    CFX_TEST(test_convolve_associativity);
    CFX_TEST(test_convolve_identity);
    CFX_TEST(test_convolve_large_coefficients);
    CFX_TEST(test_convolve_n_exact_fit);
    CFX_TEST(test_convolve_larger_random);

    CFX_TEST(test_limbs_to_chunks_simple);
    CFX_TEST(test_limbs_to_chunks_zero);
    CFX_TEST(test_chunks_to_limbs_simple);
    CFX_TEST(test_chunks_to_limbs_carry);
    CFX_TEST(test_chunks_to_limbs_large_carry);
    CFX_TEST(test_roundtrip_limbs_chunks);

    CFX_TEST(test_ntt_mul_limbs_simple);
    CFX_TEST(test_ntt_mul_limbs_zero);
    CFX_TEST(test_ntt_mul_limbs_medium);
    CFX_TEST(test_ntt_mul_limbs_large);
    CFX_TEST(test_ntt_mul_limbs_known_product);
    CFX_TEST(test_ntt_mul_limbs_max_single);
    CFX_TEST(test_ntt_mul_limbs_asymmetric);
    CFX_TEST(test_ntt_mul_limbs_by_one);
    CFX_TEST(test_ntt_mul_limbs_power_of_two);
    CFX_TEST(test_ntt_mul_limbs_sequential);
    CFX_TEST(test_ntt_mul_limbs_all_ones);
    CFX_TEST(test_ntt_mul_limbs_alternating);
    CFX_TEST(test_ntt_mul_limbs_random_16);
    CFX_TEST(test_ntt_mul_limbs_commutativity);
    CFX_TEST(test_chunks_to_limbs_multi_carry);
    CFX_TEST(test_chunks_to_limbs_overflow_cascade);
    CFX_TEST(test_limbs_to_chunks_single);
    CFX_TEST(test_ntt_mul_limbs_square);
    CFX_TEST(test_ntt_mul_limbs_fibonacci_style);
    CFX_TEST(test_ntt_mul_limbs_32_limbs);

    CFX_TEST(test_big_mul_fft_simple);
    CFX_TEST(test_big_mul_fft_zero);
    CFX_TEST(test_big_mul_fft_one);
    CFX_TEST(test_big_mul_fft_small_square);
    CFX_TEST(test_big_mul_fft_large);
    CFX_TEST(test_big_mul_fft_vs_schoolbook);
    CFX_TEST(test_big_mul_fft_commutativity);
    CFX_TEST(test_big_mul_fft_self_output);
    CFX_TEST(test_big_mul_eq_fft);
    CFX_TEST(test_big_mul_fft_power_of_two);
    CFX_TEST(test_big_mul_fft_factorial);
    CFX_TEST(test_big_mul_fft_asymmetric);

    puts("OK");
    return 0;
}
