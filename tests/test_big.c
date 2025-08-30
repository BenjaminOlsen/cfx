/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#include "cfx/big.h"
#include "cfx/macros.h"
#include "cfx/error.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>


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
    PRINT_TEST(1);
}

static void test_cfx_big_assign(void) {
    cfx_big_t a, b;
    cfx_big_init(&a);
    cfx_big_init(&b);
    cfx_limb_t al[] = {0x12371237, 0x0, 0x2, 0xFFFFFFFFFFFFFFFF};
    cfx_big_from_limbs(&a, al, sizeof(al)/sizeof(al[0]));
    cfx_big_assign(&b, &a);
    int c = cfx_big_cmp(&a, &b);
    CFX_ASSERT(c == 0);
}

static void test_copy_swap(void) {
    cfx_big_t a, b;
    cfx_big_init(&a);
    cfx_big_init(&b);
    cfx_limb_t la[] = {1, 2, 3, 4};
    cfx_limb_t lb[] = {UINT64_MAX-1, UINT64_MAX-2, UINT64_MAX-3, UINT64_MAX-4};

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
    PRINT_TEST(1);
}

static void big_init_from_limbs_base_1e9(cfx_big_t* b, const cfx_limb_t *limbs, size_t n) {
    cfx_big_init(b);
    if (n == 0) return;
    for (size_t i = n; i--;) {
        cfx_big_mul_sm(b, 1000000000ull);
        cfx_big_add_sm(b, limbs[i]);
    }
    PRINT_TEST(1);
}

static void test_mul_by_zero(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    cfx_big_from_u64(&b, 123);
    cfx_big_mul_sm(&b, 2838);
    cfx_big_mul_sm(&b, 1928);
    cfx_big_mul_sm(&b, 9);
    cfx_big_mul_sm(&b, 123765);
    cfx_big_mul_sm(&b, 0);
    size_t sz = 0;
    char* s = cfx_big_to_str(&b, &sz);
    CFX_PRINT_DBG("s: %s, sz: %zu\n", s, sz);
    CFX_ASSERT(sz == 1);
    CFX_ASSERT(strcmp(s, "0") == 0);
    free(s);
    PRINT_TEST(1);
}

static void test_add_sm(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    CFX_BIG_PRINTF(&b, "init: b.n:%ld; ", b.n);
    
    cfx_big_from_u64(&b, 123);
    CFX_BIG_PRINTF(&b, "after setting val: ");

    cfx_big_add_sm(&b, 321);
    CFX_BIG_PRINTF(&b, "after add: ");
    CFX_ASSERT(b.limb[0] == 444);
    cfx_big_from_u64(&b, UINT64_MAX);
    CFX_BIG_PRINTF(&b, "after set:");
    
    CFX_ASSERT(b.limb[0] == UINT64_MAX);
    cfx_big_add_sm(&b, 1);
    CFX_BIG_PRINTF(&b, "after carry: ");
    CFX_ASSERT(b.limb[0] == 0);
    CFX_ASSERT(b.limb[1] == 1);
    PRINT_TEST(1);
}


static void test_sub_sm(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    CFX_BIG_PRINTF(&b, "init: b.n:%ld; ", b.n);

    cfx_big_from_u64(&b, 1);
    cfx_big_sub_sm(&b, 1);
    CFX_ASSERT(b.limb[0] == 0);
    cfx_big_sub_sm(&b, 1);
    CFX_ASSERT(b.limb[0] == 0);
    
    cfx_big_from_u64(&b, UINT64_MAX);
    CFX_BIG_PRINTF(&b, "set to UINT64_MAX: ");
    CFX_ASSERT(b.limb[0] == UINT64_MAX);
    cfx_big_add_sm(&b, 1);
    CFX_BIG_PRINTF(&b, "add 1: ");
    CFX_ASSERT(b.limb[0] == 0);
    CFX_ASSERT(b.limb[1] == 1);

    cfx_big_sub_sm(&b, 1);
    CFX_BIG_PRINTF(&b, "sub 1: ");
    CFX_ASSERT(b.limb[0] == UINT64_MAX);
    CFX_ASSERT(b.n == 1);

    const cfx_limb_t N = 12;
    for (cfx_limb_t n = 0; n < N; ++n) {
        cfx_big_sub_sm(&b, 1);
        CFX_BIG_PRINTF(&b, "sub 1: ");
    }
    CFX_ASSERT(b.limb[0] == UINT64_MAX - N);

    cfx_limb_t q = 100;
    cfx_limb_t orig = b.limb[0];
    for (cfx_limb_t n = 0; n < 2; ++n) {
        cfx_acc_t s = (cfx_acc_t)b.limb[0] + q;
        cfx_big_add_sm(&b, q);
        CFX_BIG_PRINTF(&b, "add "CFX_PRIuLIMB": ", q);
        CFX_ASSERT(b.limb[0] == (cfx_limb_t)s);
        CFX_ASSERT(b.n > 1);
        CFX_ASSERT(b.limb[1] == (cfx_limb_t)(s >> CFX_LIMB_BITS));
        cfx_big_sub_sm(&b, q);
        CFX_BIG_PRINTF(&b, "sub "CFX_PRIuLIMB": ", q);
        CFX_ASSERT(b.limb[0] == orig);
    }
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

    cfx_big_sub(&t, &a);
    CFX_ASSERT(cfx_big_is_zero(&t));

    cfx_big_copy(&t, &a);
    cfx_big_add_sm(&t, 0xFABADA);
    cfx_big_sub(&t, &a);
    CFX_ASSERT(cfx_big_eq_u64(&t, 0xFABADA));

    cfx_big_copy(&t, &a);
    cfx_big_mul_sm(&t, 2);
    cfx_big_sub(&t, &a);
    CFX_ASSERT(cfx_big_eq(&t, &a));

    cfx_big_t two;
    cfx_big_init(&two);
    cfx_big_from_u64(&two, 2);

    cfx_big_copy(&t, &a);
    cfx_big_mul(&t, &two);
    cfx_big_sub(&t, &a);
    CFX_ASSERT(cfx_big_eq(&t, &a));
    PRINT_TEST(1);
}

/* helper to run one limb test */
static void check(const char *label, const cfx_limb_t *limbs, size_t n, const char *expect) {
    cfx_big_t b;
    big_init_from_limbs_base_1e9(&b, limbs, n);
    size_t len = 0;
    char *s = cfx_big_to_str(&b, &len);
    // CFX_BIG_PRINTF(&b, "str is:\n");
    // CFX_PRINT_DBG("str should be:\n%s\n", expect);
    CFX_ASSERT(strcmp(s, expect) == 0);
    CFX_ASSERT(len == strlen(expect));
    CFX_ASSERT(s[len] == '\0');
    CFX_PRINT_DBG("[ok] %s -> %s\n", label, s);
    free(s);
    cfx_big_free(&b);
}


// Single small limb
static void test_limb1(void) {
    cfx_limb_t L[] = { 123456789u };
    check("single limb", L, 1, "123456789");
    PRINT_TEST(1);
}

// -->) Single 0 limb but n>0 (should normally be normalized away; if you allow it, behavior is “000…000”?)
// Prefer: represent zero as n==0. You can skip this if you enforce normalization.

// Two limbs with inner zero-padding needed:
// value = limb[1]*1e9 + limb[0] = 42*1e9 + 123456789
static void test_limb2(void) {
    cfx_limb_t L[] = { 123456789u, 4200u };
    check("two limbs pad", L, 2, "4200123456789");
    PRINT_TEST(1);
}

// Inner limb exact zero-padding boundary: limb[0] has fewer than 9 digits
static void test_limb3(void) {
    cfx_limb_t L[] = { 1u, 1u };
    check("two limbs tiny low", L, 2, "1000000001");
    PRINT_TEST(1);
}

// Max limb values
static void test_limb4(void) {
    cfx_limb_t L[] = { 999999999u, 999999999u };
    check("two limbs max", L, 2, "999999999999999999");
    PRINT_TEST(1);
}

// Four limbs mixed
// value = 1*1e27 + 7*1e18 + 42*1e9 + 5 → "1 000000007 000000042 000000005"
static void test_limb5(void) {
    cfx_limb_t L[] = { 5u, 42u, 7u, 1u };
    check("four limbs pad", L, 4, "1000000007000000042000000005");
    PRINT_TEST(1);
}

static void test_limb6(void) {
    cfx_limb_t L[] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    check("9 limbs pad", L, 9,
        "900000000"
        "800000000"
        "700000000"
        "600000000"
        "500000000"
        "400000000"
        "300000000"
        "200000000"
        "1"
    );
    PRINT_TEST(1);
}

// Large ndigits sanity: build via mul to exercise carry
static void test_limb7(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    cfx_big_from_u64(&b, 1);
    for (int i = 0; i < 10; ++i) cfx_big_mul_sm(&b, 1000000000u - 1u); // (1e9-1)^10
    char *s = cfx_big_to_str(&b, NULL);
    // spot checks: starts with '9' and length >= 9
    CFX_PRINT_DBG("%s\n", s);
    CFX_ASSERT(s[0] == '9');
    CFX_ASSERT(strlen(s) >= 9);
    char* expect = "999999990000000044999999880000000209999999"
        "748000000209999999880000000044999999990000000001";
    CFX_ASSERT(strcmp(s, expect) == 0);
    free(s);
    cfx_big_free(&b);
    PRINT_TEST(1);
}

static void test_str1(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    const char *sin =   "99911231231238761239876981273469128374169283476129"
                        "38471629384761250389172603459812630498672312387123"
                        "87123981723918273912891238719248719238719248169723"
                        "00091203901290909090911100091231283761000101023882";
    cfx_big_from_str(&b, sin);
    char *sout = cfx_big_to_str(&b, NULL);
    int ok = (strcmp(sin, sout) == 0); 
    CFX_PRINT_DBG("test str1: in:\n%s \n%s\nout.. %s\n", sin, sout,
        ok ? "ok":"NOT ok");
    CFX_ASSERT(ok);
    PRINT_TEST(1);
}

static void test_str2(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    const char *sin = "9218";
    cfx_big_from_str(&b, sin);
    char *sout = cfx_big_to_str(&b, NULL);
    int ok = (strcmp(sin, sout) == 0); 
    CFX_PRINT_DBG("test str: in:\n%s \n%s\nout.. %s\n", sin, sout,
        ok ? "ok":"NOT ok");
    CFX_ASSERT(ok);
    PRINT_TEST(1);
}


// --- hex tests ---
static void test_hex_zero_empty_n(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    // n == 0 should yield "0"
    size_t len = 12345;
    char* s = cfx_big_to_hex(&b, &len);
    CFX_ASSERT(s);
    CFX_ASSERT(strcmp(s, "0") == 0);
    CFX_ASSERT(len == 1);
    CFX_ASSERT(s[len] == '\0');
    printf("[%s] - %s\n", __func__, s);
    free(s);
    PRINT_TEST(1);
}

static void test_hex_zero_explicit_limb_zero(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    cfx_limb_t limbs[] = {0};
    cfx_big_from_limbs(&b, limbs, 1);
    size_t len = 0;
    char* s = cfx_big_to_hex(&b, &len);
    CFX_ASSERT(s);
    CFX_ASSERT(strcmp(s, "0") == 0);
    CFX_ASSERT(len == 1);
    CFX_ASSERT(s[len] == '\0');
    cfx_big_free(&b);
    printf("[%s] - %s\n", __func__, s);
    free(s);
    PRINT_TEST(1);
}

static void test_hex_single_limb_basic(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    // 0x1 -> "1"
    cfx_limb_t limbs1[] = {0x1ull};
    cfx_big_from_limbs(&b, limbs1, 1);
    size_t len1 = 0;
    char* s1 = cfx_big_to_hex(&b, &len1);
    CFX_ASSERT(strcmp(s1, "1") == 0);
    CFX_ASSERT(len1 == strlen("1"));
    CFX_ASSERT(s1[len1] == '\0');
    printf("[%s] - %s\n", __func__, s1);
    free(s1);

    // 0xabcdef -> "abcdef" (lowercase)
    b.limb[0] = 0xabcdefull;
    size_t len2 = 0;
    char* s2 = cfx_big_to_hex(&b, &len2);
    CFX_ASSERT(strcmp(s2, "abcdef") == 0);
    CFX_ASSERT(len2 == strlen("abcdef"));
    CFX_ASSERT(s2[len2] == '\0');
    cfx_big_free(&b);
    printf("[%s] - %s\n", __func__, s2);
    free(s2);
    PRINT_TEST(1);
}

static void test_hex_single_limb_hex_digit_count(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    // 2^60 = 0x1000000000000000 -> 1 followed by 15 zeros (16 digits total)
    cfx_limb_t limbs[] = {0x1000000000000000ull};
    cfx_big_from_limbs(&b, limbs, 1);
    size_t len = 0;
    char* s = cfx_big_to_hex(&b, &len);
    CFX_ASSERT(strcmp(s, "1000000000000000") == 0);
    CFX_ASSERT(len == 16);
    CFX_ASSERT(s[len] == '\0');
    cfx_big_free(&b);
    printf("[%s] - %s\n", __func__, s);
    free(s);
    PRINT_TEST(1);
}

static void test_hex_two_limbs_padding(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    // high = 0x1, low = 0x1 -> "1" + "0000000000000001"
    cfx_limb_t limbs[] = {
        0x0000000000000001ull, // low
        0x0000000000000001ull  // high
    };
    cfx_big_from_limbs(&b, limbs, 2);
    size_t len = 0;
    char* s = cfx_big_to_hex(&b, &len);
    CFX_ASSERT(strcmp(s, "10000000000000001") == 0);
    CFX_ASSERT(len == 17);                  // 1 + 16
    CFX_ASSERT(s[len] == '\0');
    cfx_big_free(&b);
    printf("[%s] - %s\n", __func__, s);
    free(s);
    PRINT_TEST(1);
}

static void test_hex_two_limbs_mixed_digits(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    // high = 0xABC, low = 0x0011223344556677
    // expect: "abc" + "0011223344556677"
    cfx_limb_t limbs[] = {
        0x0011223344556677ull, // low
        0x0000000000000ABCull  // high
    };
    cfx_big_from_limbs(&b, limbs, 2);
    char* s = cfx_big_to_hex(&b, NULL); // also test sz_out == NULL
    CFX_ASSERT(strcmp(s, "abc0011223344556677") == 0);
    CFX_ASSERT(s[strlen(s)] == '\0');
    cfx_big_free(&b);
    printf("[%s] - %s\n", __func__, s);
    free(s);
    PRINT_TEST(1);
}

static void test_hex_leading_zero_limb_skipped(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    // limbs: [low=X, mid=Y, high=0] -> should behave like just [low=X, mid=Y]
    cfx_limb_t limbs3[] = {
        0xDEADBEEFCAFEBABEuLL, // low
        0x0000000000000123uLL, // mid
        0x0000000000000000uLL  // high (zero)
    };
    cfx_big_from_limbs(&b, limbs3, 3);
    size_t len = 0;
    char* s = cfx_big_to_hex(&b, &len);
    // expected: "123" + "%016" of low
    const char* expect_low = "deadbeefcafebabe"; // lowercase
    // const char* expect = "1230000000000000" "deadbeefcafebabe";
    // But careful: mid=0x123 -> "123", low padded to 16
    char expect_buf[3 + 16 + 1];
    snprintf(expect_buf, sizeof(expect_buf), "%s%s", "123", expect_low);
    CFX_ASSERT(strcmp(s, expect_buf) == 0);
    CFX_ASSERT(len == strlen(expect_buf));
    CFX_ASSERT(s[len] == '\0');
    cfx_big_free(&b);
    printf("[%s] - %s\n", __func__, s);
    free(s);
    PRINT_TEST(1);
}

static void test_hex_no_leading_zeros_on_msl(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    // high = 0x00000000000000ab -> "ab", low zero-padded
    cfx_limb_t limbs[] = {
        0ULL,
        0xABULL
    };
    cfx_big_from_limbs(&b, limbs, 2);
    size_t len = 0;
    char* s = cfx_big_to_hex(&b, &len);
    CFX_ASSERT(strcmp(s, "ab0000000000000000") == 0);
    CFX_ASSERT(len == 2 + 16);
    CFX_ASSERT(s[len] == '\0');
    cfx_big_free(&b);
    free(s);
    PRINT_TEST(1);
}
static void test_cache(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    CFX_ASSERT(b.cache == NULL);
    cfx_big_enable_cache(&b);
    CFX_ASSERT(b.cache != NULL);

    cfx_big_from_u64(&b, 1);
    // CFX_ASSERT(b.cache->primes.data == NULL);
    // CFX_ASSERT(b.cache->state == CFX_FAC_FULL); todo
    PRINT_TEST(1);
}

static void test_zero_right(void) {
    cfx_big_t b, m;
    cfx_big_init(&b);
    cfx_big_from_u64(&b, 123);
    cfx_big_init(&m);
    cfx_big_from_u64(&m, 0);

    cfx_big_mul(&b, &m);
    CFX_ASSERT(cfx_big_is_zero(&b));

    cfx_big_free(&b);
    cfx_big_free(&m);
    PRINT_TEST(1);
}

static void test_zero_left(void) {
    cfx_big_t b, m;
    cfx_big_init(&b);
    cfx_big_from_u64(&b, 0);
    cfx_big_init(&m);
    cfx_big_from_u64(&m, 987);

    cfx_big_mul(&b, &m);
    
    CFX_ASSERT(cfx_big_is_zero(&b));

    cfx_big_free(&b);
    cfx_big_free(&m);
    PRINT_TEST(1);
}

static void big_expect_limbs(const char* s, const cfx_big_t* b, const cfx_limb_t* limbs, size_t n) {
    if (b->n != n) {
        printf("[%s]: size mismatch! n: %zu, b->n: %zu!\n", s, n, b->n);
    }
    CFX_ASSERT(b->n == n);
    int ok = 1;
    for (size_t i = 0; i < n; ++i) {
        if (b->limb[i] != limbs[i]) {
            printf("[%s]: limb mismatch at idx %zu!\n", s, i);
            ok = 0;
            break;
        }
    }
    CFX_ASSERT(ok);
    PRINT_TEST(1);
}

static void test_mul1(void) {
    cfx_big_t b, m;
    cfx_big_init(&b);
    cfx_big_from_u64(&b, 5);
    cfx_big_init(&m);
    cfx_big_from_u64(&m, 7);

    cfx_big_mul(&b, &m);
    cfx_limb_t expect[] = {35};
    CFX_BIG_PRINT(b);
    big_expect_limbs(__func__, &b, expect, 1);

    cfx_big_free(&b);
    cfx_big_free(&m);
    PRINT_TEST(1);
}

static void test_mul_adduiv(void) {
    cfx_big_t a, m;
    cfx_big_init(&a);
    cfx_big_from_u64(&a, 0);
    cfx_big_init(&m);
    cfx_big_from_u64(&m, 1);
    cfx_limb_t K = 0x1B2CDE;
    cfx_big_mul_sm(&m, K);

    for (cfx_limb_t i = 0; i < K; ++i) {
        cfx_big_add_sm(&a, 1);
    }
    CFX_ASSERT(cfx_big_eq(&a, &m));
    PRINT_TEST(1);
}

static void test_carry_two_limbs_times_2(void) {
    cfx_big_t b, m;
    cfx_big_init(&b);
    cfx_limb_t limbs_b[] = {0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFFFFFFull};
    cfx_big_from_limbs(&b, limbs_b, 2);
    cfx_big_init(&m);
    cfx_big_from_u64(&m, 2);

    cfx_big_mul(&b, &m);
    cfx_limb_t expect[] = {0xFFFFFFFFFFFFFFFEull, 0xFFFFFFFFFFFFFFFFull, 1ull};
    big_expect_limbs(__func__, &b, expect, 3);

    cfx_big_free(&b);
    cfx_big_free(&m);
    PRINT_TEST(1);
}


static void test_mul_by_base_2_64_shift(void) {
    // Multiply by 2^64 (limbs = [0,1]) should shift by one limb
    cfx_big_t b, m;
    cfx_big_init(&b);
    cfx_big_init(&m);
    cfx_limb_t limbs_b[] = {0x0123456789ABCDEFull, 0x0FEDCBA987654321ull, 0x0000000000000001ull};
    size_t sz0 = sizeof(limbs_b)/sizeof(limbs_b[0]);
    cfx_big_from_limbs(&b, limbs_b, sz0);
    cfx_limb_t limbs_m[] = {0ull, 1ull};
    cfx_big_from_limbs(&m, limbs_m, 2);

    const size_t N = 1111;

    for (size_t n = 1; n < N; ++n) {

        cfx_big_mul(&b, &m);
        cfx_limb_t* expect = (cfx_limb_t*)malloc((n + sz0)*sizeof(cfx_limb_t));
        for (size_t k = 0; k < n; ++k) {
            expect[k] = 0;
        }
        expect[n+sz0-3] = 0x0123456789ABCDEFull,
        expect[n+sz0-2] = 0x0FEDCBA987654321ull;
        expect[n+sz0-1] = 0x0000000000000001ull;

        big_expect_limbs(__func__, &b, expect, sz0 + n);
        // printf("[test_mul_by_base_2_64_shift]: tested shift of %zu OK\n", n);
        free(expect);
        expect = NULL;
    }

    cfx_big_free(&b);
    cfx_big_free(&m);
    PRINT_TEST(1);
}

static void test_self_multiply_square(void) {
    // (2^64 - 1)^2 = [0xFFFFFFFFFFFFFFFE, 1] in base 2^64
    cfx_big_t b;
    cfx_big_init(&b);
    cfx_limb_t limbs[] = {0xFFFFFFFFFFFFFFFFull};
    cfx_big_from_limbs(&b, limbs, 1);
    big_expect_limbs(__func__, &b, limbs, 1);
    cfx_big_mul(&b, &b); // self-mul path
    cfx_limb_t expect[] = {1ull, 0xFFFFFFFFFFFFFFFEull};
    big_expect_limbs(__func__, &b, expect, 2);
    cfx_big_free(&b);
    PRINT_TEST(1);
}

static void test_self_multiply_big(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    size_t N = 10;
    cfx_limb_t* limbs = (cfx_limb_t*)malloc(N * sizeof(cfx_limb_t));
    for (size_t i = 0; i < N; ++i) {
        limbs[i] = 2;
    }
    cfx_big_from_limbs(&b, limbs, N);
    printf("before: "); CFX_BIG_PRINT(b);
    cfx_big_mul(&b, &b);
    printf("after: "); CFX_BIG_PRINT(b);
    cfx_big_free(&b);
    free(limbs);
    limbs = NULL;
    PRINT_TEST(1);
}


void test_mul(cfx_big_t* b, const cfx_big_t* m) {
    if(cfx_big_is_zero(m)) {
        printf("multiplying b by zero!\n");
        cfx_big_free(b);
        cfx_big_from_u64(b, 0);
        return;
    }
    // if (cfx_big_eq(b, m)) {
    //     cfx_big_sq(b);
    //     return;
    // }
    PRINT_TEST(1);
    
}

// static void write_string_wrapped(const char* s, const char* fn, size_t w) {
//     FILE* f = fopen(fn, "w");
//     if (!f) {
//         perror("fopen");
//         return;
//     }
//     size_t len = strlen(s);
//     for (size_t i = 0; i < len; i += w) {
//         fprintf(f, "\"");
//         for (size_t j = i; j < i + w && j < len; ++j) {
//             fputc(s[j], f);
//         }
//         fprintf(f, "\"\n");
//     }
//     fclose(f);
// }

static void test_known_squares(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    cfx_big_from_str(&b, 
        "12554203470773361528352143580257209"
        "759168353591939024551938");
    cfx_big_sq(&b);
    char* s = cfx_big_to_str(&b, NULL);
    char* expect = 
        "15760802478557791686620405668794173"
        "58808686358020659979337980952437065"
        "51032029871143396883518477623727858"
        "361659555844";
    CFX_ASSERT(strcmp(s, expect) == 0);
    
    /* ----------------------------------- */
    cfx_big_from_str(&b,
        "14536774485912137811774516281687980"
        "27136112775646765168338161504493023"
        "20618275753867907499968765052767290"
        "9553511701551770336148695547906");
    cfx_big_sq(&b);
    free(s);
    s = cfx_big_to_str(&b, NULL);
    
    expect = 
        "211317812454266098563847037101386611"
        "572936979272149330366075109813640724"
        "474345860994392765176988025566104455"
        "837520931631539188377309778866973982"
        "401217516000816262383671909413419933"
        "179663682199874639617901251522391335"
        "258229693927038827878426679857709370"
        "5798799065540984836";
    CFX_ASSERT(strcmp(s, expect) == 0);

    /* ----------------------------------- */
    cfx_big_from_str(&b,
        "788040123927889584288300542721290477751811"
        "061317446299559616920266928982601233965513"
        "84510299169195903266945438318594");
    cfx_big_sq(&b);
    free(s);
    s = cfx_big_to_str(&b, NULL);
    
    expect = 
        "621007236920283574126921534924820203000914"
        "006088690925637430299255754323320524243949"
        "820727721413359346225814236510862646127210"
        "854935998980345500359712802179223868607115"
        "902506409479615101839926583004283044656249"
        "4079870273849846136836";
    CFX_ASSERT(strcmp(s, expect) == 0);

    /* ----------------------------------- */
    cfx_big_from_str(&b,
        "4946608029462090681478206578991795742708644"
        "2564742658586426229545514803499564697000372"
        "3095350971345437292114654548843072761868784"
        "674125049315509629339381027496416991005114370");
    cfx_big_sq(&b); // 1
    free(s);
    s = cfx_big_to_str(&b, NULL);
    
    expect = 
        "24468930997138827791465884302231872460739875"
        "68207399045598121707346936565872814004700782"
        "89425671088010093095401304027496760620656589"
        "55442387983480255066626620918677531232865400"
        "98771942569641370987258548756784517935582567"
        "13069390644656846533223333372022161629347849"
        "16188372769675748466345955863374240885700087"
        "1949664098243244687781332547496780496900";
    CFX_ASSERT(strcmp(s, expect) == 0);

    /* ----------------------------------- */
    cfx_big_sq(&b); // 2
    free(s);
    s = cfx_big_to_str(&b, NULL);
    
    expect = 
        "598728584142741349308708530097477655730195121"
        "8921561456478992373941874807025827989117005224"
        "4544159576837861157818585384355732464079298824"
        "8840307403854743209208295687980560227400969404"
        "2099334329107992954561531802255073931559398056"
        "4137316892818829069898815609946436650172187960"
        "4903063742870707747793924033318816999047826082"
        "6483886794388751538092273101337790920556525050"
        "8202721113875166122496565273190156790054685080"
        "2591709257787294752819636451016941276946114924"
        "0557203469761158225922311367092258610357178141"
        "4512816029237243505436725602625497437379829908"
        "1001420696314530701776191559963115206965790843"
        "4600574609622913370900822575527296202868593414"
        "9525207968313075153130098691732396070700210909"
        "610000";
    CFX_ASSERT(strcmp(s, expect) == 0);

    cfx_big_sq(&b); // 4
    free(s);
    s = cfx_big_to_str(&b, NULL);
    
    expect = 
        "3584759174695717079200799669904391027371015034"
        "5591263067167166620789779901562032648556807601"
        "2950877610090545985260292606958082794183173841"
        "0531271568492233387532134452461086273878119171"
        "5847641737220633768273029588942704704129884200"
        "9958682084547837574688649786501681990292232657"
        "7415652116159576158769762214147356000386767514"
        "9511670598205070085582557818076333629441697422"
        "9875224640887417011148169060064563559950573932"
        "5226205227402913380194309571812791018912125575"
        "1791447739605363127819124149182487530334381706"
        "8896663768263012311430826689280351346703702625"
        "8078023280801979632626881628873401907467070917"
        "5815446662382885870340476544847226508933219513"
        "1908572998226578819696551060038441941135034737"
        "0068773357966050317170281108465947743312221154"
        "8985206867725840086507000464147209225780870815"
        "2470140579958636771833023777601016674168494135"
        "5848715099748171254374285813092185256574867981"
        "4922088922467495708073886349364048519341765660"
        "2078730664244522169366402136000455068825485572"
        "3466936983629748500464140600538157767988061984"
        "0563990463715200572971418532092893645704169342"
        "4426635419183083905555687209236852205715503369"
        "2336623040698494499300396883453818464450499457"
        "8627966191480124089611256089345895606330526380"
        "0963369509588825439132963159814190607878092811"
        "6978547138214917414204956687046025729039394600"
        "8260306269245606621871407784590216117503519011"
        "6801526199652563056546897824927378333686359035"
        "2100000000";
    printf("\n\n%s\n\n", s);
    CFX_ASSERT(strcmp(s, expect) == 0);
    // CFX_BIG_PRINT_LIMBS(b);

    
    // cfx_big_sq(&b); // 8
    cfx_big_mul(&b, &b);
    free(s);
    s = cfx_big_to_str(&b, NULL);
    //// sanity check:
    cfx_big_t B;
    cfx_big_init(&B);
    cfx_big_from_str(&B, s);
    CFX_ASSERT(cfx_big_eq(&B, &b));
    char* sanity = cfx_big_to_str(&B, NULL);
    CFX_ASSERT(strcmp(s, sanity) == 0);
    expect = "12850498340565118640831124663943566637163477270761965804775357209550765230785703949035729472531170991846950765062231843931603991577837122702202578619413163343802199724298959199688594791894654764450428039375530903550022460719636135025802643105174385371702288121193336984280156229260813405468551910887207681430170841269088143925314707610300410631962624300469472858167288519605269872627946462373021700902442150967139315116172545474219540016600664971616358393885751109046842533050183826430384246617880667515303122705936214635228967259835005388450384900687056571344933537245392429971062926387006758181202634581358151550436165392630946855812938662698021351517729227827971635271649805374176176940281753018739171849301949232131565817532850243174625134450254366859930465467865935871669065946539235916390634180117949435696071838514048795111101267436279103470592457465410839430058374274997666026789906916220271637971862818408213637491993426120588003613013514121426918778425609885983132067309940493824037246336471087139347566932289527147579507444510833470558987620105627392491106342931918229603332167519231845701934759523730028587564778522689440283636622636742434431091107839566802176801151089094414637433966684982281939928036568083598873318392485777034090015016265841172493029035363177944812390618909768667573763207436642271055957361540865636837910531318607009603088734768184385262471743568922118461465472820244894124358177759866070832621634448723670657830089717393358547747703735801523850299292208462216521430657525373415183422463571048104152842920726565698681208502981595141456253796350916896570590072884174607249064969441514018160468583314653443352671437479825779995459466312201730735736110098653787331865116531446616641987026900818215864623614884767160893264417488979936319828793408610211926232514931040704633236494245315339349224851289392387332474248568479216841581124860328079900988652092392614618200827173162666704740632419460637566175022280383967839231370086128119069239566679344622785639260673870367122861560404871658960189532546990946727895958575281330314215816959430714079369740156022375856316371677709024505475091778006456248457155666092229065835357913477405119329372494342721686545800319902824829662799906780886313825440342423658177641941252419690985063319979924638260305722020541317833802446921285816583681614997935512759860955018042317788791616381896041300997716291956236555872471276333548522003017903198166576927518763140985360706323427487549508484530837381768891367647567673106372095645829146242426938921111540219157093806251149323044120727645913143734374450245354473392643178616466237265916832878586939031607485456329408461025836997224109513396266781143113993671870210102616033541932236431301110708951036119790977069093536028698906465137616493507440201974410000000000000000";
    cfx_big_from_str(&B, expect);
    free(sanity);
    sanity = cfx_big_to_str(&B, NULL);
    CFX_ASSERT(strcmp(sanity, expect) == 0);

    for (size_t i = 0; i < b.n; ++i) {
        printf("calculated: b.limb[%zu]: "CFX_PRIuLIMB"; correct: B.limb[%zu]: "CFX_PRIuLIMB": %s: diff %d\n",
        i, b.limb[i], i, B.limb[i], b.limb[i] == B.limb[i] ? "ok" : "--- NOT OK",
        (int)(b.limb[i] & 0xFFFF) - (int)(B.limb[i] & 0xFFFF));
    }
    printf("\n\n");
    CFX_ASSERT(cfx_big_eq(&B, &b));
    CFX_ASSERT(strcmp(s, expect) == 0);
    int cnt = 0;

    cfx_big_mul(&b, &b); // 16
    printf("mul %d len: %zu \n", ++cnt, b.n);

    cfx_big_mul(&b, &b); // 32
    printf("mul %d len: %zu \n", ++cnt, b.n);
    cfx_big_mul(&b, &b); // 64
    printf("mul %d len: %zu \n", ++cnt, b.n);
    cfx_big_mul(&b, &b); // 128
    printf("mul %d len: %zu \n", ++cnt, b.n);
    cfx_big_mul(&b, &b); // 256
    printf("mul %d len: %zu \n", ++cnt, b.n);
    cfx_big_mul(&b, &b); // 512
    printf("mul %d len: %zu \n", ++cnt, b.n);
    // cfx_big_mul(&b, &b); // 1024
    // printf("mul %d len: %zu \n", ++cnt, b.n);
    // cfx_big_mul(&b, &b); // 2048
    // printf("mul %d len: %zu \n", ++cnt, b.n);
    // cfx_big_mul(&b, &b); // 4096
    // printf("mul %d len: %zu \n", ++cnt, b.n);
    // cfx_big_mul(&b, &b); // 8192
    // printf("mul %d len: %zu \n", ++cnt, b.n);

    // for (size_t i = 0; i < b.n; ++i) {
    //     printf("calculated: b.limb[%zu]: "CFX_PRIuLIMB" (0x"CFX_PRI0xLIMB")\n", i, b.limb[i], b.limb[i]);
    // }
    
    char* huge = cfx_big_to_hex(&b, NULL);
    // printf("%s\n", huge);
    printf("digits: %zu\n", strlen(huge));
    free(huge);
    free(s);
    cfx_big_free(&b);
    PRINT_TEST(1);
}

static void test_known_squares_2(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    cfx_big_from_str(&b,
        "4946608029462090681478206578991795742708644"
        "2564742658586426229545514803499564697000372"
        "3095350971345437292114654548843072761868784"
        "674125049315509629339381027496416991005114370");
    cfx_big_mul_csa(&b, &b); // 1
    char* s = cfx_big_to_str(&b, NULL);
    
    char* expect = 
        "24468930997138827791465884302231872460739875"
        "68207399045598121707346936565872814004700782"
        "89425671088010093095401304027496760620656589"
        "55442387983480255066626620918677531232865400"
        "98771942569641370987258548756784517935582567"
        "13069390644656846533223333372022161629347849"
        "16188372769675748466345955863374240885700087"
        "1949664098243244687781332547496780496900";
    printf("calculated\n%s\nexpected\n%s\n", s, expect);
    CFX_ASSERT(strcmp(s, expect) == 0);

    /* ----------------------------------- */
    cfx_big_mul_csa(&b, &b); // 2
    free(s);
    s = cfx_big_to_str(&b, NULL);
    
    expect = 
        "598728584142741349308708530097477655730195121"
        "8921561456478992373941874807025827989117005224"
        "4544159576837861157818585384355732464079298824"
        "8840307403854743209208295687980560227400969404"
        "2099334329107992954561531802255073931559398056"
        "4137316892818829069898815609946436650172187960"
        "4903063742870707747793924033318816999047826082"
        "6483886794388751538092273101337790920556525050"
        "8202721113875166122496565273190156790054685080"
        "2591709257787294752819636451016941276946114924"
        "0557203469761158225922311367092258610357178141"
        "4512816029237243505436725602625497437379829908"
        "1001420696314530701776191559963115206965790843"
        "4600574609622913370900822575527296202868593414"
        "9525207968313075153130098691732396070700210909"
        "610000";
    CFX_ASSERT(strcmp(s, expect) == 0);

    cfx_big_mul_csa(&b, &b); // 4
    free(s);
    s = cfx_big_to_str(&b, NULL);
    // write_string_wrapped(s, "", 80);
    expect = 
        "3584759174695717079200799669904391027371015034"
        "5591263067167166620789779901562032648556807601"
        "2950877610090545985260292606958082794183173841"
        "0531271568492233387532134452461086273878119171"
        "5847641737220633768273029588942704704129884200"
        "9958682084547837574688649786501681990292232657"
        "7415652116159576158769762214147356000386767514"
        "9511670598205070085582557818076333629441697422"
        "9875224640887417011148169060064563559950573932"
        "5226205227402913380194309571812791018912125575"
        "1791447739605363127819124149182487530334381706"
        "8896663768263012311430826689280351346703702625"
        "8078023280801979632626881628873401907467070917"
        "5815446662382885870340476544847226508933219513"
        "1908572998226578819696551060038441941135034737"
        "0068773357966050317170281108465947743312221154"
        "8985206867725840086507000464147209225780870815"
        "2470140579958636771833023777601016674168494135"
        "5848715099748171254374285813092185256574867981"
        "4922088922467495708073886349364048519341765660"
        "2078730664244522169366402136000455068825485572"
        "3466936983629748500464140600538157767988061984"
        "0563990463715200572971418532092893645704169342"
        "4426635419183083905555687209236852205715503369"
        "2336623040698494499300396883453818464450499457"
        "8627966191480124089611256089345895606330526380"
        "0963369509588825439132963159814190607878092811"
        "6978547138214917414204956687046025729039394600"
        "8260306269245606621871407784590216117503519011"
        "6801526199652563056546897824927378333686359035"
        "2100000000";
    printf("\n\n%s\n\n", s);
    CFX_ASSERT(strcmp(s, expect) == 0);
    // CFX_BIG_PRINT_LIMBS(b);

    
    cfx_big_mul(&b, &b); // 8
    free(s);
    s = cfx_big_to_str(&b, NULL);

    /* ---------- sanity check: ---------------------------------- */
    #if 0
    cfx_big_t B;
    cfx_big_init(&B);
    cfx_big_from_str(&B, s);
    CFX_ASSERT(cfx_big_eq(&B, &b));
    char* sanity = cfx_big_to_str(&B, NULL);
    CFX_ASSERT(strcmp(s, sanity) == 0);
    expect = "12850498340565118640831124663943566637163477270761965804775357209550765230785703949035729472531170991846950765062231843931603991577837122702202578619413163343802199724298959199688594791894654764450428039375530903550022460719636135025802643105174385371702288121193336984280156229260813405468551910887207681430170841269088143925314707610300410631962624300469472858167288519605269872627946462373021700902442150967139315116172545474219540016600664971616358393885751109046842533050183826430384246617880667515303122705936214635228967259835005388450384900687056571344933537245392429971062926387006758181202634581358151550436165392630946855812938662698021351517729227827971635271649805374176176940281753018739171849301949232131565817532850243174625134450254366859930465467865935871669065946539235916390634180117949435696071838514048795111101267436279103470592457465410839430058374274997666026789906916220271637971862818408213637491993426120588003613013514121426918778425609885983132067309940493824037246336471087139347566932289527147579507444510833470558987620105627392491106342931918229603332167519231845701934759523730028587564778522689440283636622636742434431091107839566802176801151089094414637433966684982281939928036568083598873318392485777034090015016265841172493029035363177944812390618909768667573763207436642271055957361540865636837910531318607009603088734768184385262471743568922118461465472820244894124358177759866070832621634448723670657830089717393358547747703735801523850299292208462216521430657525373415183422463571048104152842920726565698681208502981595141456253796350916896570590072884174607249064969441514018160468583314653443352671437479825779995459466312201730735736110098653787331865116531446616641987026900818215864623614884767160893264417488979936319828793408610211926232514931040704633236494245315339349224851289392387332474248568479216841581124860328079900988652092392614618200827173162666704740632419460637566175022280383967839231370086128119069239566679344622785639260673870367122861560404871658960189532546990946727895958575281330314215816959430714079369740156022375856316371677709024505475091778006456248457155666092229065835357913477405119329372494342721686545800319902824829662799906780886313825440342423658177641941252419690985063319979924638260305722020541317833802446921285816583681614997935512759860955018042317788791616381896041300997716291956236555872471276333548522003017903198166576927518763140985360706323427487549508484530837381768891367647567673106372095645829146242426938921111540219157093806251149323044120727645913143734374450245354473392643178616466237265916832878586939031607485456329408461025836997224109513396266781143113993671870210102616033541932236431301110708951036119790977069093536028698906465137616493507440201974410000000000000000";
    cfx_big_from_str(&B, expect);
    free(sanity);
    sanity = cfx_big_to_str(&B, NULL);
    CFX_ASSERT(strcmp(sanity, expect) == 0);

    for (size_t i = 0; i < b.n; ++i) {
        printf("calculated: b.limb[%zu]: "CFX_PRIuLIMB"; correct: B.limb[%zu]: "CFX_PRIuLIMB": %s: diff %d\n",
        i, b.limb[i], i, B.limb[i], b.limb[i] == B.limb[i] ? "ok" : "--- NOT OK",
        (int)(b.limb[i] & 0xFFFF) - (int)(B.limb[i] & 0xFFFF));
    }
    printf("\n\n");
    CFX_ASSERT(cfx_big_eq(&B, &b));
    CFX_ASSERT(strcmp(s, expect) == 0);
    #endif
    /* ----------------------------------------------------------- */

    int cnt = 0;
    cfx_big_mul_csa(&b, &b); // 16
    printf("mul csa %d len: %zu \n", ++cnt, b.n);
    cfx_big_mul_csa(&b, &b); // 32
    printf("mul csa %d len: %zu \n", ++cnt, b.n);
    cfx_big_mul_csa(&b, &b); // 64
    printf("mul csa %d len: %zu \n", ++cnt, b.n);
    cfx_big_mul_csa(&b, &b); // 128
    printf("mul csa %d len: %zu \n", ++cnt, b.n);
    cfx_big_mul_csa(&b, &b); // 256
    printf("mul csa %d len: %zu \n", ++cnt, b.n);
    cfx_big_mul_csa(&b, &b); // 512
    printf("mul csa %d len: %zu \n", ++cnt, b.n);
    // cfx_big_mul_csa(&b, &b); // 1024
    // printf("mul csa %d len: %zu \n", ++cnt, b.n);
    // cfx_big_mul_csa(&b, &b); // 2048
    // printf("mul csa %d len: %zu \n", ++cnt, b.n);
    // cfx_big_mul_csa(&b, &b); // 4096
    // printf("mul csa %d len: %zu \n", ++cnt, b.n);

    // for (size_t i = 0; i < b.n; ++i) {
    //     printf("calculated: b.limb[%zu]: "CFX_PRIuLIMB" (0x"CFX_PRI0xLIMB")\n", i, b.limb[i], b.limb[i]);
    // }
    
    char* huge = cfx_big_to_hex(&b, NULL);
    // printf("%s\n", huge);
    printf("digits: %zu\n", strlen(huge));
    free(huge);
    free(s);
    cfx_big_free(&b);
    PRINT_TEST(1);
}
/* ------------------------------------------------------------------------ */


static void expect_dec_eq(const cfx_big_t* x, const char* dec) {
    cfx_big_t tmp; 
    cfx_big_init(&tmp);
    cfx_big_from_str(&tmp, dec);
    CFX_ASSERT(cfx_big_eq(x, &tmp));
    cfx_big_free(&tmp);
}

static void big_dec1(cfx_big_t* x) {
    cfx_limb_t borrow = 1;
    for (size_t i = 0; i < x->n && borrow; ++i) {
        cfx_limb_t old = x->limb[i];
        x->limb[i] = old - borrow;
        borrow = (x->limb[i] > old);
    }
    // rely on library ops to keep canonical form later
}

// Property check: n == q*d + r and r < d
static void assert_n_eq_qd_plus_r(const cfx_big_t* n, const cfx_big_t* q,
                                  const cfx_big_t* d, const cfx_big_t* r) {
    cfx_big_t check;
    cfx_big_init(&check);
    cfx_big_copy(&check, q);
    cfx_big_t tmp;
    cfx_big_init(&tmp);

    cfx_big_mul(&check, d);
    cfx_big_copy(&tmp, r);
    cfx_big_add(&check, &tmp);

    CFX_ASSERT(cfx_big_eq(&check, n));
    CFX_ASSERT(cfx_big_cmp(r, d) == -1);

    cfx_big_free(&check);
    cfx_big_free(&tmp);
}

// ---------- tests ----------

void test_big_div_divide_by_zero(void) {
    cfx_big_t n, d, q, r;
    cfx_big_init(&n);
    cfx_big_init(&d);
    cfx_big_init(&q);
    cfx_big_init(&r);
    cfx_big_from_u64(&n, 123);
    cfx_big_from_u64(&d, 0);

    int rc = cfx_big_divrem(&q, &r, &n, &d);
    CFX_ASSERT(rc == -1);

    cfx_big_free(&n);
    cfx_big_free(&d);
    cfx_big_free(&q);
    cfx_big_free(&r);
}

void test_big_div_zero_dividend(void) {
    cfx_big_t n, d, q, r;
    cfx_big_init(&n);
    cfx_big_init(&d);
    cfx_big_init(&q);
    cfx_big_init(&r);
    cfx_big_from_u64(&n, 0);
    cfx_big_from_u64(&d, 42);

    int rc = cfx_big_divrem(&q, &r, &n, &d);
    CFX_ASSERT(rc == 0);
    expect_dec_eq(&q, "0");
    expect_dec_eq(&r, "0");

    cfx_big_free(&n);
    cfx_big_free(&d);
    cfx_big_free(&q);
    cfx_big_free(&r);
}

void test_big_div_n_less_than_d(void) {
    cfx_big_t n, d, q, r;
    cfx_big_init(&n);
    cfx_big_init(&d);
    cfx_big_init(&q);
    cfx_big_init(&r);
    cfx_big_from_u64(&n, 123456);
    cfx_big_from_u64(&d, 123456789);

    int rc = cfx_big_divrem(&q, &r, &n, &d);
    CFX_ASSERT(rc == 0);
    expect_dec_eq(&q, "0");
    CFX_ASSERT(cfx_big_eq(&r, &n));

    cfx_big_free(&n);
    cfx_big_free(&d);
    cfx_big_free(&q);
    cfx_big_free(&r);
}

void test_big_div_equal_numbers(void) {
    cfx_big_t n, d, q, r;
    cfx_big_init(&n);
    cfx_big_init(&d);
    cfx_big_init(&q);
    cfx_big_init(&r);
    cfx_big_from_str(&n, "1234567890123456789012345678901234567890");
    cfx_big_copy(&d, &n);

    int rc = cfx_big_divrem(&q, &r, &n, &d);
    CFX_ASSERT(rc == 0);
    expect_dec_eq(&q, "1");
    expect_dec_eq(&r, "0");

    cfx_big_free(&n);
    cfx_big_free(&d);
    cfx_big_free(&q);
    cfx_big_free(&r);
}

void test_big_div_single_limb_divisor_property(void) {
    cfx_big_t n, d, q, r;
    cfx_big_init(&n);
    cfx_big_init(&d);
    cfx_big_init(&q);
    cfx_big_init(&r);
    cfx_big_from_str(&n, "340282366920938463463374607431768211455"); // 2^128 - 1
    cfx_big_from_u64(&d, 123456789ULL);

    int rc = cfx_big_divrem(&q, &r, &n, &d);
    CFX_ASSERT(rc == 0);
    assert_n_eq_qd_plus_r(&n, &q, &d, &r);

    cfx_big_free(&n);
    cfx_big_free(&d);
    cfx_big_free(&q);
    cfx_big_free(&r);
}

void test_big_div_multi_limb_divisor_exact_and_remainder(void) {
    // n = a*b + r, then n / b -> q=a, rem=r
    cfx_big_t a, b, r, n, q, rem;
    cfx_big_init(&a);
    cfx_big_init(&b);
    cfx_big_init(&n);
    cfx_big_init(&q);
    cfx_big_init(&r);
    cfx_big_init(&rem);
    cfx_big_from_str(&a, "123456789012345678901234567890123456789");
    cfx_big_from_str(&b, "987654321098765432109876543210987654321");
    cfx_big_from_str(&r, "12345678901234567890");

    cfx_big_copy(&n, &a);
    cfx_big_mul(&n, &b);
    cfx_big_add(&n, &r);

    int rc = cfx_big_divrem(&q, &rem, &n, &b);
    CFX_ASSERT(rc == 0);
    CFX_ASSERT(cfx_big_eq(&q, &a));
    CFX_ASSERT(cfx_big_eq(&rem, &r));
    assert_n_eq_qd_plus_r(&n, &q, &b, &rem);

    cfx_big_free(&a); 
    cfx_big_free(&b);
    cfx_big_free(&r);
    cfx_big_free(&n);
    cfx_big_free(&q);
    cfx_big_free(&rem);
}

void test_big_div_in_place_eq_with_remainder(void) {
    // Build n = a*b + 42, then n := n / b, rem = 42
    cfx_big_t a, b, n, rem, forty_two;
    cfx_big_init(&a);
    cfx_big_init(&b);
    cfx_big_init(&n);
    cfx_big_init(&rem);
    cfx_big_init(&forty_two);
    cfx_big_from_str(&a, "1122334455667788990011223344556677889900");
    cfx_big_from_str(&b, "18446744073709551616"); // 2^64
    cfx_big_copy(&n, &a);
    cfx_big_mul(&n, &b);
    cfx_big_from_u64(&forty_two, 42);
    cfx_big_add(&n, &forty_two);

    int rc = cfx_big_div_eq(&n, &b, &rem);
    CFX_ASSERT(rc == 0);
    CFX_ASSERT(cfx_big_eq(&n, &a));
    expect_dec_eq(&rem, "42");

    cfx_big_free(&a);
    cfx_big_free(&b);
    cfx_big_free(&n);
    cfx_big_free(&rem);
    cfx_big_free(&forty_two);
}

void test_big_div_quotient_only_and_remainder_only(void) {
    cfx_big_t a, b, n, q, r, b_minus_1;
    cfx_big_init(&a);
    cfx_big_init(&b);
    cfx_big_init(&n);
    cfx_big_init(&q);
    cfx_big_init(&r);
    cfx_big_init(&b_minus_1);
    cfx_big_from_str(&a, "3141592653589793238462643383279502884197");
    cfx_big_from_str(&b, "2718281828459045235360287471352662497757");

    // n = a*b + (b-1)
    cfx_big_copy(&n, &a);
    cfx_big_mul(&n, &b);

    cfx_big_copy(&b_minus_1, &b);
    big_dec1(&b_minus_1);

    cfx_big_add(&n, &b_minus_1);

    // quotient only
    CFX_ASSERT(cfx_big_div_out(&q, &n, &b) == 0);
    CFX_ASSERT(cfx_big_eq(&q, &a));

    // remainder only
    CFX_ASSERT(cfx_big_mod(&r, &n, &b) == 0);
    CFX_ASSERT(cfx_big_eq(&r, &b_minus_1));

    cfx_big_free(&a);
    cfx_big_free(&b);
    cfx_big_free(&n);
    cfx_big_free(&q);
    cfx_big_free(&r);
    cfx_big_free(&b_minus_1);
}

void test_big_div_alias_remainder_eq_src(void) {
    CFX_ASSERT(1);
    /* TODO! */
    #if 0
    // Verify cfx_big_div_eq supports r == b (if your impl promises this).
    cfx_big_t b, d;
    cfx_big_init(&b);
    cfx_big_init(&d);
    cfx_big_from_str(&b, "123456789012345678901234567890");
    cfx_big_from_str(&d, "987654321");

    cfx_big_t orig;
    cfx_big_init(&orig);
    
    cfx_big_copy(&orig, &b);

    int rc = cfx_big_div_eq(&b, &d, &b);   // r aliases src
    CFX_ASSERT(rc == 0);

    // Check property with a fresh recompute: orig == q*d + r, where q is in 'b' (after div)
    cfx_big_t q_copy;
    cfx_big_init(&q_copy);
    cfx_big_copy(&q_copy, &b);  // b now holds q
    assert_n_eq_qd_plus_r(&orig, &q_copy, &d, &b);

    cfx_big_free(&orig);
    cfx_big_free(&q_copy);
    cfx_big_free(&b);
    cfx_big_free(&d);
    #endif
}

static void test_shifts(void) {
    cfx_big_t a, b;
    cfx_big_init(&a);
    cfx_big_init(&b);

    cfx_limb_t lsb = 123, msb = 283912;
    cfx_limb_t limbs[] = {lsb, msb};
    cfx_big_from_limbs(&a, limbs, sizeof(limbs)/sizeof(limbs[0]));
    cfx_big_copy(&b, &a);

    cfx_big_shl_bits(&b, &b, 1);
    CFX_ASSERT(b.limb[1] == 2 * msb);
    CFX_ASSERT(b.limb[0] == 2 * lsb);
    
    cfx_big_shl_bits(&b, &b, 1);
    CFX_ASSERT(b.limb[1] == 4 * msb);
    CFX_ASSERT(b.limb[0] == 4 * lsb);
    
    cfx_big_shr_bits(&b, &b, 1);
    CFX_ASSERT(b.limb[1] == 2 * msb);
    CFX_ASSERT(b.limb[0] == 2 * lsb);
    
    cfx_big_shr_bits(&b, &b, 1);
    CFX_ASSERT(cfx_big_cmp(&a, &b) == 0);


}

/* ---- helpers ----------------------------------------------------------- */

static void assert_hex_eq(const char* tag, const cfx_big_t* x, const char* hex_exp) {
    char* got = cfx_big_to_hex(x, NULL);
    int ok = (strcmp(got, hex_exp) == 0);
    if (!ok) {
        fprintf(stderr, "[%s] expected 0x%s, got 0x%s\n", tag, hex_exp, got);
        // fflush(stderr);
    }
    CFX_ASSERT_PRINT(ok);
    free(got);
}

static void check_shl_case(const char* msg, const char* hex_in, unsigned s, const char* hex_exp) {
    /* out-of-place */
    cfx_big_t a, out;
    cfx_big_init(&a); 
    cfx_big_init(&out);
    cfx_big_from_hex(&a, hex_in);
    cfx_big_shl_bits(&out, &a, s);
    char obuf[128];
    snprintf(obuf, sizeof obuf, "%s - %s", msg, "shl oopl");
    assert_hex_eq(obuf, &out, hex_exp);
    cfx_big_free(&out);

    /* in-place (aliasing) */
    cfx_big_shl_bits(&a, &a, s);
    char ibuf[128];
    snprintf(ibuf, sizeof ibuf, "%s - %s", msg, "shl inpl");
    assert_hex_eq(ibuf, &a, hex_exp);
    cfx_big_free(&a);
}

static void check_shr_case(const char* msg, const char* hex_in, unsigned s, const char* hex_exp) {
    /* out-of-place */
    cfx_big_t a, out;
    cfx_big_init(&a);
    cfx_big_init(&out);
    cfx_big_from_hex(&a, hex_in);
    cfx_big_shr_bits(&out, &a, s);
    char obuf[128];
    snprintf(obuf, sizeof obuf, "%s - %s", msg, "shr oopl");
    assert_hex_eq(obuf, &out, hex_exp);
    cfx_big_free(&out);

    /* in-place (aliasing) */
    cfx_big_shr_bits(&a, &a, s);
    char ibuf[128];
    snprintf(ibuf, sizeof ibuf, "%s - %s", msg, "shr inpl");
    assert_hex_eq(ibuf, &a, hex_exp);
    cfx_big_free(&a);
}

#define CHECK_SHL_CASE(hex_in, s, hex_out) check_shl_case(__func__, hex_in, s, hex_out)
#define CHECK_SHR_CASE(hex_in, s, hex_out) check_shr_case(__func__, hex_in, s, hex_out)

/* ---- tests ------------------------------------------------------------- */

static void test_shl_basic_identity(void) {
    /* shift by 0: identity */
    CHECK_SHL_CASE("0", 0, "0");
    CHECK_SHL_CASE("1", 0, "1");
    CHECK_SHL_CASE("deadbeef", 0, "deadbeef");
}

static void test_shl_create_top_limb(void) {
    /* new top limb created when MSB carries out */
    CHECK_SHL_CASE("8000000000000000", 1, "10000000000000000");    // 2^63 << 1 -> 2^64
    CHECK_SHL_CASE("1", 64, "10000000000000000");                   // append 16 hex zeros
    CHECK_SHL_CASE("1", 68, "100000000000000000");                 // 16 zeros + 1 hex nibble
}

static void test_shl_cross_limb_1bit(void) {
    /* two-limb pattern: hi=1, lo=0x8000...0001, shift left by 1
       expected: hi' = 3, lo' = 2  => hex "3" || "000...0002" */
    CHECK_SHL_CASE("100000000000000008000000000000001", 1, "200000000000000010000000000000002");
}

static void test_shr_basic_identity(void) {
    /* shift by 0: identity */
    CHECK_SHR_CASE("0", 0, "0");
    CHECK_SHR_CASE("1", 0, "1");
    CHECK_SHR_CASE("deadbeef", 0, "deadbeef");
}

static void test_shr_drop_whole_limb(void) {
    /* exact 64-bit drops */
    CHECK_SHR_CASE("10000000000000000", 64, "1");         // 2^64 >> CFX_LIMB_BITS
    CHECK_SHR_CASE("100000000000000000", 68, "1");       // (2^68) >> 68
}

static void test_shr_to_zero(void) {
    /* shifting past total bit-length -> zero */
    CHECK_SHR_CASE("1", 65, "0");
    CHECK_SHR_CASE("ffffffffffffffff", 128, "0");
}

static void test_shr_cross_limb_carry_6bits(void) {
    /* This reproduces the cross-limb carry you just debugged:
       Input limbs: hi=0x3e, lo=0x0c40  -> hex "3e0000000000000c40"
       Right shift by 6:
         r1 = 0x3e >> 6 = 0
         carry = (0x3e & 0x3f) << 58 = 0xf800000000000000
         r0 = (0x0c40 >> 6) | carry = 0x31 | 0xf800000000000000 = 0xf800000000000031
       Trim top zero limb => "f800000000000031"
    */
    CHECK_SHR_CASE("3e0000000000000c40", 6, "f800000000000031");
}

static void test_shr_mixed_cases(void) {
    /* assorted mixes with odd s and multiple limbs */
    CHECK_SHR_CASE("100000000000000000000", 4, "10000000000000000000");  // divide by 16
    CHECK_SHR_CASE("abcdef0123456789", 4, "abcdef012345678");          // >>4 removes low nibble
    CHECK_SHR_CASE("abcdef0123456789", 60, "a");                         // >>68 = >>64 then >>4
}

/* ------------------------------------------------------------------ */
static void expect_limb_pattern(const cfx_big_t* a, size_t n,
                                cfx_limb_t l2, cfx_limb_t l1, cfx_limb_t l0) {
    CFX_ASSERT(a->n == n);
    if (n >= 1) CFX_ASSERT(a->limb[0] == l0);
    if (n >= 2) CFX_ASSERT(a->limb[1] == l1);
    if (n >= 3) CFX_ASSERT(a->limb[2] == l2);
}

void test_exp_edge_cases(void) {
    cfx_big_t n, p, out;
    cfx_big_init(&n); cfx_big_init(&p); cfx_big_init(&out);

    // 0^0 := 1
    cfx_big_from_u64(&n, 0);
    cfx_big_from_u64(&p, 0);
    cfx_big_exp(&out, &n, &p);
    CFX_ASSERT(cfx_big_eq_u64(&out, 1));

    // 0^k = 0 (k>0)
    cfx_big_from_u64(&n, 0);
    cfx_big_from_u64(&p, 5);
    cfx_big_exp(&out, &n, &p);
    CFX_ASSERT(cfx_big_eq_u64(&out, 0));

    // 1^p = 1
    cfx_big_from_u64(&n, 1);
    cfx_big_from_u64(&p, 1234567);
    cfx_big_exp(&out, &n, &p);
    CFX_ASSERT(cfx_big_eq_u64(&out, 1));

    // n^0 = 1
    cfx_big_from_u64(&n, 42);
    cfx_big_from_u64(&p, 0);
    cfx_big_exp(&out, &n, &p);
    CFX_ASSERT(cfx_big_eq_u64(&out, 1));

    cfx_big_free(&out); cfx_big_free(&p); cfx_big_free(&n);
}

void test_exp_small_values(void) {
    cfx_big_t n, p, out;
    cfx_big_init(&n); cfx_big_init(&p); cfx_big_init(&out);

    // n^1 = n
    cfx_big_from_u64(&n, 123456789);
    cfx_big_from_u64(&p, 1);
    cfx_big_exp(&out, &n, &p);
    CFX_ASSERT(out.n == 1 && out.limb[0] == 123456789ULL);

    // 2^10 = 1024
    cfx_big_from_u64(&n, 2);
    cfx_big_from_u64(&p, 10);
    cfx_big_exp(&out, &n, &p);
    CFX_ASSERT(cfx_big_eq_u64(&out, 1024ULL));

    // 3^5 = 243
    cfx_big_from_u64(&n, 3);
    cfx_big_from_u64(&p, 5);
    cfx_big_exp(&out, &n, &p);
    CFX_ASSERT(cfx_big_eq_u64(&out, 243ULL));

    cfx_big_free(&out); cfx_big_free(&p); cfx_big_free(&n);
}

void test_exp_powers_of_two_boundaries(void) {
    cfx_big_t n, p, out;
    cfx_big_init(&n); cfx_big_init(&p); cfx_big_init(&out);

    // 2^64 = 1<<64 -> limbs [0]=0, [1]=1
    cfx_big_from_u64(&n, 2);
    cfx_big_from_u64(&p, 64);
    cfx_big_exp(&out, &n, &p);
    CFX_ASSERT(out.n == 2);
    expect_limb_pattern(&out, 2, 0, 1, 0);

    // 2^128 = 1<<128 -> limbs [0]=0, [1]=0, [2]=1
    cfx_big_from_u64(&p, 128);
    cfx_big_exp(&out, &n, &p);
    CFX_ASSERT(out.n == 3);
    expect_limb_pattern(&out, 3, 1, 0, 0);

    // 2^127 -> highest bit in limb[1], limb[0]=0
    cfx_big_from_u64(&p, 127);
    cfx_big_exp(&out, &n, &p);
    CFX_ASSERT(out.n == 2);
    CFX_ASSERT(out.limb[0] == 0);
    CFX_ASSERT(out.limb[1] == (1ULL << 63));

    cfx_big_free(&out); cfx_big_free(&p); cfx_big_free(&n);
}

void test_exp_aliasing(void) {
    cfx_big_t n, p;
    cfx_big_init(&n); cfx_big_init(&p);

    // out aliases base (out == n)
    cfx_big_from_u64(&n, 7);
    cfx_big_from_u64(&p, 6);      // 7^6 = 117,649
    cfx_big_exp(&n, &n, &p);
    CFX_ASSERT(n.n == 2 || n.n == 1);
    // 117,649 fits in 64 bits
    CFX_ASSERT(n.n == 1 && n.limb[0] == 117649ULL);

    // out aliases exponent (out == p)
    cfx_big_from_u64(&n, 2);
    cfx_big_from_u64(&p, 80);     // 2^80 -> limbs [0]=0, [1]=0x1000000000000000, [2]=0x0000000000000001
    cfx_big_exp(&p, &n, &p);
    CFX_ASSERT(p.n == 2 || p.n == 3);
    // exact check:
    CFX_ASSERT(p.n == 2 || p.n == 3); // (library may trim leading zero limbs)
    if (p.n == 2) {
        // 2^80 = 1<<80 => limb[0]=0, limb[1]=1<<16
        CFX_ASSERT(p.limb[0] == 0);
        CFX_ASSERT(p.limb[1] == (1ULL << 16));
    } else {
        // If your representation keeps a third limb for carry, adjust accordingly.
        CFX_ASSERT(0 && "Unexpected limb count for 2^80");
    }

    cfx_big_free(&p); cfx_big_free(&n);
}

void test_exp_compare_with_naive_mul(void) {
    cfx_big_t n, p, out1, out2;
    cfx_big_init(&n); cfx_big_init(&p);
    cfx_big_init(&out1); cfx_big_init(&out2);

    // Random-ish small cases to avoid huge numbers; compare against repeated multiply.
    cfx_limb_t bases[] = {2,3,5,10,17,1234567};
    cfx_limb_t exps[]  = {2,3,4,5,8,16};

    for (size_t i=0;i<sizeof(bases)/sizeof(bases[0]);++i) {
        for (size_t j=0;j<sizeof(exps)/sizeof(exps[0]);++j) {
            cfx_big_from_u64(&n, bases[i]);
            cfx_big_from_u64(&p, exps[j]);

            // exp under test
            cfx_big_exp(&out1, &n, &p);

            // naive: res=1; repeat e times: res*=n
            cfx_big_from_u64(&out2, 1);
            cfx_limb_t e = exps[j];
            while (e--) {
                cfx_big_mul_auto(&out2, &n);
            }

            CFX_ASSERT(cfx_big_cmp(&out1, &out2) == 0);
        }
    }

    cfx_big_free(&out2); cfx_big_free(&out1);
    cfx_big_free(&p); cfx_big_free(&n);
}

/* ------------------------------------------------------------------ */

int main(void) {
    CFX_TEST(test_cfx_big_assign);
    CFX_TEST(test_copy_swap);
    CFX_TEST(test_cfx_big_init);
    CFX_TEST(test_cfx_big_reserve);
    CFX_TEST(test_mul_by_zero);
    CFX_TEST(test_limb1);
    CFX_TEST(test_limb2);
    CFX_TEST(test_limb3);
    CFX_TEST(test_limb4);
    CFX_TEST(test_limb5);
    CFX_TEST(test_limb6);
    CFX_TEST(test_limb7);
    CFX_TEST(test_str1);
    CFX_TEST(test_str2);
    CFX_TEST(test_add_sm);
    CFX_TEST(test_sub);
    CFX_TEST(test_sub_sm);
    CFX_TEST(test_cache);
    CFX_TEST(test_zero_right);
    CFX_TEST(test_zero_left);
    CFX_TEST(test_mul1);
    CFX_TEST(test_carry_two_limbs_times_2);
    CFX_TEST(test_mul_by_base_2_64_shift);
    CFX_TEST(test_self_multiply_square);
    CFX_TEST(test_self_multiply_big);
    CFX_TEST(test_known_squares);
    CFX_TEST(test_known_squares_2);
    CFX_TEST(test_mul_adduiv);
    CFX_TEST(test_big_div_divide_by_zero);
    CFX_TEST(test_big_div_zero_dividend);
    CFX_TEST(test_big_div_n_less_than_d);
    CFX_TEST(test_big_div_equal_numbers);
    CFX_TEST(test_big_div_single_limb_divisor_property);

    CFX_TEST(test_big_div_multi_limb_divisor_exact_and_remainder);
    CFX_TEST(test_big_div_in_place_eq_with_remainder);
    CFX_TEST(test_big_div_quotient_only_and_remainder_only);
    CFX_TEST(test_big_div_alias_remainder_eq_src);
    
    CFX_TEST(test_hex_leading_zero_limb_skipped);
    CFX_TEST(test_hex_no_leading_zeros_on_msl);
    CFX_TEST(test_hex_single_limb_basic);
    CFX_TEST(test_hex_single_limb_hex_digit_count);
    CFX_TEST(test_hex_two_limbs_mixed_digits);
    CFX_TEST(test_hex_two_limbs_padding);
    CFX_TEST(test_hex_zero_empty_n);
    CFX_TEST(test_hex_zero_explicit_limb_zero);
    CFX_TEST(test_shifts);
    CFX_TEST(test_shl_basic_identity);
    CFX_TEST(test_shl_create_top_limb);
    CFX_TEST(test_shl_cross_limb_1bit);
    CFX_TEST(test_shr_basic_identity);
    CFX_TEST(test_shr_drop_whole_limb);
    CFX_TEST(test_shr_to_zero);
    CFX_TEST(test_shr_cross_limb_carry_6bits);
    CFX_TEST(test_shr_mixed_cases);

    CFX_TEST(test_exp_edge_cases);
    CFX_TEST(test_exp_small_values);
    CFX_TEST(test_exp_powers_of_two_boundaries);
    CFX_TEST(test_exp_aliasing);
    CFX_TEST(test_exp_compare_with_naive_mul);
    puts("OK");
    return 0;
}
