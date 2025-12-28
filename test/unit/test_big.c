/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#include "cfx/big.h"
#include "cfx/macros.h"
#include "cfx/primes.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdio.h>

#ifndef _MSC_VER
#include <strings.h>
#endif


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
    cfx_limb_t al[] = {0x12371237, 0x0, 0x2, CFX_LIMB_MAX-1};
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

static void big_init_from_limbs_base_1e9(cfx_big_t* b, const cfx_limb_t *limbs, size_t n) {
    cfx_big_init(b);
    if (n == 0) return;
    for (size_t i = n; i--;) {
        cfx_big_mul_sm(b, UINT64_C(1000000000));
        cfx_big_add_sm(b, limbs[i]);
    }
    PRINT_TEST(1);
}

static void test_mul_by_zero(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    cfx_big_from_limb(&b, 123);
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

    cfx_big_from_limb(&b, 123);
    CFX_BIG_PRINTF(&b, "after setting val: ");

    cfx_big_add_sm(&b, 321);
    CFX_BIG_PRINTF(&b, "after add: ");
    CFX_ASSERT(b.limb[0] == 444);
    cfx_big_from_limb(&b, CFX_LIMB_MAX);
    CFX_BIG_PRINTF(&b, "after set:");

    CFX_ASSERT(b.limb[0] == CFX_LIMB_MAX);
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

    cfx_big_from_limb(&b, 1);
    cfx_big_sub_sm(&b, 1);
    CFX_ASSERT(b.limb[0] == 0);
    cfx_big_sub_sm(&b, 1);
    CFX_ASSERT(b.limb[0] == 0);

    cfx_big_from_limb(&b, CFX_LIMB_MAX);
    CFX_BIG_PRINTF(&b, "set to CFX_LIMB_MAX: ");
    CFX_ASSERT(b.limb[0] == CFX_LIMB_MAX);
    cfx_big_add_sm(&b, 1);
    CFX_BIG_PRINTF(&b, "add 1: ");
    CFX_ASSERT(b.limb[0] == 0);
    CFX_ASSERT(b.limb[1] == 1);

    cfx_big_sub_sm(&b, 1);
    CFX_BIG_PRINTF(&b, "sub 1: ");
    CFX_ASSERT(b.limb[0] == CFX_LIMB_MAX);
    CFX_ASSERT(b.n == 1);

    const cfx_limb_t N = 12;
    for (cfx_limb_t n = 0; n < N; ++n) {
        cfx_big_sub_sm(&b, 1);
        CFX_BIG_PRINTF(&b, "sub 1: ");
    }
    CFX_ASSERT(b.limb[0] == CFX_LIMB_MAX - N);

    cfx_limb_t q = 100;
    cfx_limb_t orig = b.limb[0];
    for (cfx_limb_t n = 0; n < 2; ++n) {
        cfx_acc_t s = (cfx_acc_t)b.limb[0] + q;
        cfx_big_add_sm(&b, q);
        CFX_BIG_PRINTF(&b, "add " CFX_PRIuLIMB ": ", q);
        CFX_ASSERT(b.limb[0] == (cfx_limb_t)s);
        CFX_ASSERT(b.n > 1);
        CFX_ASSERT(b.limb[1] == (cfx_limb_t)(s >> CFX_LIMB_BITS));
        cfx_big_sub_sm(&b, q);
        CFX_BIG_PRINTF(&b, "sub " CFX_PRIuLIMB ": ", q);
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
    cfx_big_from_limb(&two, 2);

    cfx_big_copy(&t, &a);
    cfx_big_mul_eq(&t, &two);
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
    /* CFX_BIG_PRINTF(&b, "str is:\n"); */
    /* CFX_PRINT_DBG("str should be:\n%s\n", expect); */
    CFX_ASSERT(strcmp(s, expect) == 0);
    CFX_ASSERT(len == strlen(expect));
    CFX_ASSERT(s[len] == '\0');
    CFX_PRINT_DBG("[ok] %s -> %s\n", label, s);
    free(s);
    cfx_big_free(&b);
}

static void test_limb1(void) {
    cfx_limb_t L[] = { 123456789u };
    check("single limb", L, 1, "123456789");
    PRINT_TEST(1);
}

/* Two limbs with inner zero-padding needed: */
/* value = limb[1]*1e9 + limb[0] = 42*1e9 + 123456789 */
static void test_limb2(void) {
    cfx_limb_t L[] = { 123456789u, 4200u };
    check("two limbs pad", L, 2, "4200123456789");
    PRINT_TEST(1);
}

/* Inner limb exact zero-padding boundary: limb[0] has fewer than 9 digits */
static void test_limb3(void) {
    cfx_limb_t L[] = { 1u, 1u };
    check("two limbs tiny low", L, 2, "1000000001");
    PRINT_TEST(1);
}

/* Max limb values */
static void test_limb4(void) {
    cfx_limb_t L[] = { 999999999u, 999999999u };
    check("two limbs max", L, 2, "999999999999999999");
    PRINT_TEST(1);
}

/* Four limbs mixed */
/* value = 1*1e27 + 7*1e18 + 42*1e9 + 5 -> "1 000000007 000000042 000000005" */
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

/* Large ndigits sanity: build via mul to exercise carry */
static void test_limb7(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    cfx_big_from_limb(&b, 1);
    for (int i = 0; i < 10; ++i) cfx_big_mul_sm(&b, 1000000000u - 1u); /* (1e9-1)^10 */
    char *s = cfx_big_to_str(&b, NULL);
    /* spot checks: starts with '9' and length >= 9 */
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
    cfx_big_from_dec(&b, sin);
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
    cfx_big_from_dec(&b, sin);
    char *sout = cfx_big_to_str(&b, NULL);
    int ok = (strcmp(sin, sout) == 0);
    CFX_PRINT_DBG("test str: in:\n%s \n%s\nout.. %s\n", sin, sout,
        ok ? "ok":"NOT ok");
    CFX_ASSERT(ok);
    PRINT_TEST(1);
}

static void test_from_str_matches1(void) {
    cfx_big_t b1, b2;
    cfx_big_init(&b1);
    cfx_big_init(&b2);

    const char* s1 = "0xFF";
    const char* s2 = "255";

    cfx_big_from_str(&b1, s1);
    cfx_big_from_str(&b2, s2);

    CFX_ASSERT(cfx_big_eq(&b1, &b2));
}

static void test_from_str_matches2(void) {
    cfx_big_t b1, b2, b3;
    cfx_big_init(&b1);
    cfx_big_init(&b2);
    cfx_big_init(&b3);

    const char* s1 = "0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF";
    const char* s2 = "1393796574908163946345982392040522594123775";
    const char* s3 = "0b11111111111111111111111111111111111111111"
                    "11111111111111111111111111111111111111111111"
                    "1111111111111111111111111111111111111111111111111111111";


    cfx_big_from_str(&b1, s1);
    cfx_big_from_str(&b2, s2);
    cfx_big_from_str(&b3, s3);

    CFX_ASSERT(cfx_big_eq(&b1, &b2));
    CFX_ASSERT(cfx_big_eq(&b1, &b3));
}

static void test_from_str_zero_forms(void) {
    cfx_big_t z1, z2, z3, z4;
    cfx_big_init(&z1);
    cfx_big_init(&z2);
    cfx_big_init(&z3);
    cfx_big_init(&z4);

    CFX_ASSERT(cfx_big_from_str(&z1, "0") == 0);
    CFX_ASSERT(cfx_big_from_str(&z2, "0000000") == 0);
    CFX_ASSERT(cfx_big_from_str(&z3, "0x0") == 0);
    CFX_ASSERT(cfx_big_from_str(&z4, "0b0") == 0);

    CFX_ASSERT(cfx_big_eq(&z1, &z2));
    CFX_ASSERT(cfx_big_eq(&z1, &z3));
    CFX_ASSERT(cfx_big_eq(&z1, &z4));
}

static void test_from_str_whitespace_ok(void) {
    cfx_big_t a, b, c;
    cfx_big_init(&a);
    cfx_big_init(&b);
    cfx_big_init(&c);

    CFX_ASSERT(cfx_big_from_str(&a, "   255") == 0);
    CFX_ASSERT(cfx_big_from_str(&b, "255   \n\t") == 0);
    CFX_ASSERT(cfx_big_from_str(&c, " \t  0xFF \r\n") == 0);

    /* a == b == 0xFF */
    CFX_ASSERT(cfx_big_eq(&a, &b));
    CFX_ASSERT(cfx_big_eq(&a, &c));
}

static void test_from_str_hex_prefix_case(void) {
    cfx_big_t a, b, c;
    cfx_big_init(&a);
    cfx_big_init(&b);
    cfx_big_init(&c);

    CFX_ASSERT(cfx_big_from_str(&a, "0xdeadBEEF") == 0);
    CFX_ASSERT(cfx_big_from_str(&b, "0XDEADBEEF") == 0);
    CFX_ASSERT(cfx_big_from_str(&c, "3735928559") == 0); /* 0xDEADBEEF */

    CFX_ASSERT(cfx_big_eq(&a, &b));
    CFX_ASSERT(cfx_big_eq(&a, &c));
}

static void test_from_str_bin_prefix_case(void) {
    cfx_big_t a, b;
    cfx_big_init(&a);
    cfx_big_init(&b);

    CFX_ASSERT(cfx_big_from_str(&a, "0b101010") == 0);
    CFX_ASSERT(cfx_big_from_str(&b, "0B101010") == 0);

    CFX_ASSERT(cfx_big_eq(&a, &b));
}


/* Strictness: public wrapper should reject junk after the number */
static void test_from_str_rejects_trailing_junk(void) {
    cfx_big_t a;
    cfx_big_init(&a);

    CFX_ASSERT(cfx_big_from_str(&a, "255x") < 0);
    CFX_ASSERT(cfx_big_from_str(&a, "0xFFgg") < 0);
    CFX_ASSERT(cfx_big_from_str(&a, "0b10102") < 0);
}

/* Reject empty / whitespace-only */
static void test_from_str_rejects_empty(void) {
    cfx_big_t a;
    cfx_big_init(&a);

    CFX_ASSERT(cfx_big_from_str(&a, "") < 0);
    CFX_ASSERT(cfx_big_from_str(&a, "   \t\r\n") < 0);
}

/* Reject prefix with no digits */
static void test_from_str_rejects_prefix_only(void) {
    cfx_big_t a;
    cfx_big_init(&a);

    CFX_ASSERT(cfx_big_from_str(&a, "0x") < 0);
    CFX_ASSERT(cfx_big_from_str(&a, "0b") < 0);
    CFX_ASSERT(cfx_big_from_str(&a, "0X   ") < 0);
}

/* without 'legacy octal': bare leading 0 should be decimal */
static void test_from_str_no_legacy_octal(void) {
    cfx_big_t a, b;
    cfx_big_init(&a);
    cfx_big_init(&b);

    CFX_ASSERT(cfx_big_from_str(&a, "010") == 0);
    CFX_ASSERT(cfx_big_from_str(&b, "10") == 0);
    CFX_ASSERT(cfx_big_eq(&a, &b));
}

static void test_from_str_limb_boundary_all_ones(void) {
    cfx_big_t a, b;
    cfx_big_init(&a);
    cfx_big_init(&b);

    /* 2^64 - 1 */
    CFX_ASSERT(cfx_big_from_str(&a, "0xFFFFFFFFFFFFFFFF") == 0);
    CFX_ASSERT(cfx_big_from_str(&b, "18446744073709551615") == 0);

    CFX_ASSERT(cfx_big_eq(&a, &b));
}

static void test_from_str_limb_boundary_power_of_two(void) {
    cfx_big_t a, b;
    cfx_big_init(&a);
    cfx_big_init(&b);

    /* 2^64 */
    CFX_ASSERT(cfx_big_from_str(&a, "0x10000000000000000") == 0);
    CFX_ASSERT(cfx_big_from_str(&b, "18446744073709551616") == 0);

    CFX_ASSERT(cfx_big_eq(&a, &b));
}

static void test_scan_num_hex_basic(void) {
    cfx_big_t b, ref;
    cfx_big_init(&b);
    cfx_big_init(&ref);

    const char *s = "0xFFF+1";
    size_t n = 12345;

    CFX_ASSERT(cfx_big_scan_num_n(&b, (const uint8_t*)s, strlen(s), &n) == 0);
    CFX_ASSERT(n == 5); /* "0xFFF" */

    CFX_ASSERT(cfx_big_from_str(&ref, "4095") == 0);
    CFX_ASSERT(cfx_big_eq(&b, &ref));
}

static void test_scan_num_bin_basic(void) {
    cfx_big_t b, ref;
    cfx_big_init(&b);
    cfx_big_init(&ref);

    const char *s = "0b0000000001+7";
    size_t n = 0;

    int rc = cfx_big_scan_num_n(&b, (const uint8_t*)s, strlen(s), &n);
    printf("rc=%d n=%zu\n", rc, n);
    CFX_ASSERT(cfx_big_scan_num_n(&b, (const uint8_t*)s, strlen(s), &n) == 0);
    CFX_ASSERT(n == 12); /* "0b" + 10 digits */

    CFX_ASSERT(cfx_big_from_str(&ref, "1") == 0);
    CFX_ASSERT(cfx_big_eq(&b, &ref));
}

static void test_scan_num_b64_basic(void) {
    cfx_big_t b, ref;
    cfx_big_init(&b);
    cfx_big_init(&ref);

    /* "b64:AQ==" decodes to bytes {0x01} => value 1 (big-endian magnitude) */
    const char *s = "b64:AQ==+7";
    size_t n = 0;

    CFX_ASSERT(cfx_big_scan_num_n(&b, (const uint8_t*)s, strlen(s), &n) == 0);
    CFX_ASSERT(n == 8); /* "b64:" (4) + "AQ==" (4) */

    CFX_ASSERT(cfx_big_from_str(&ref, "1") == 0);
    CFX_ASSERT(cfx_big_eq(&b, &ref));
}

static void test_scan_num_b64_whitespace_inside(void) {
    cfx_big_t b, ref;
    cfx_big_init(&b);
    cfx_big_init(&ref);

    /* still decodes to 1, decoder ignores whitespace */
    const char *s = "b64: A Q = =   +1";
    size_t n = 12345;

    CFX_ASSERT(cfx_big_scan_num_n(&b, (const uint8_t*)s, strlen(s), &n) == 0);

    /* scanned portion should include the whitespace after "b64:" and inside */
    /* "b64:" = 4, then " A Q = =   " = 11  => total 15 */
    printf("n: %zu\n",n);
    CFX_ASSERT(n == 15);
    CFX_ASSERT(cfx_big_from_str(&ref, "1") == 0);
    CFX_ASSERT(cfx_big_eq(&b, &ref));
}

static void test_scan_num_b64_multibyte_value(void) {
    cfx_big_t b, ref;
    cfx_big_init(&b);
    cfx_big_init(&ref);

    /* bytes {0x01,0x00} => value 256, base64 is "AQA=" */
    const char* s = "b64:AQA=+1";
    size_t n = 0;

    CFX_ASSERT(cfx_big_scan_num_n(&b, (const uint8_t*)s, strlen(s), &n) == 0);
    CFX_ASSERT(n == 8); /* "b64:" + "AQA=" */
    CFX_ASSERT(cfx_big_from_str(&ref, "256") == 0);
    CFX_ASSERT(cfx_big_eq(&b, &ref));
}

static void test_scan_num_b64_stops_before_junk(void) {
    cfx_big_t b, ref;
    cfx_big_init(&b);
    cfx_big_init(&ref);

    /*
      The scanner should consume only the base64-ish span and stop at '!'
      so scan_num_n succeeds and leaves n pointing to end of b64 payload.
    */
    const char *s = "b64:AQ==!+7";
    size_t n = 0;

    CFX_ASSERT(cfx_big_scan_num_n(&b, (const uint8_t*)s, strlen(s), &n) == 0);
    CFX_ASSERT(n == 8); /* "b64:" + "AQ==" */

    CFX_ASSERT(cfx_big_from_str(&ref, "1") == 0);
    CFX_ASSERT(cfx_big_eq(&b, &ref));
}

static void test_from_str_b64_allows_trailing_spaces_only(void) {
    cfx_big_t b, ref;
    cfx_big_init(&b);
    cfx_big_init(&ref);

    /* should succeed: only trailing whitespace after the base64 */
    CFX_ASSERT(cfx_big_from_str(&b, "b64:AQ==   \r\n\t") == 0);

    CFX_ASSERT(cfx_big_from_str(&ref, "1") == 0);
    CFX_ASSERT(cfx_big_eq(&b, &ref));

    /* should fail: junk after base64 (non-space) */
    CFX_ASSERT(cfx_big_from_str(&b, "b64:AQ==XYZ") != 0);
}

static void test_from_str_b64_rejects_invalid(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    CFX_ASSERT(cfx_big_from_str(&b, "b64:Zm9v*") != 0);
    CFX_ASSERT(cfx_big_from_str(&b, "b64:Zg=") != 0);   /* wrong length */
    CFX_ASSERT(cfx_big_from_str(&b, "b64:Zg=A") != 0);  /* non-canonical tail bits */
    CFX_ASSERT(cfx_big_from_str(&b, "b64:   \n\t") != 0);
}


static void test_scan_num_oct_basic(void) {
    cfx_big_t b, ref;
    cfx_big_init(&b);
    cfx_big_init(&ref);

    const char *s = "0o377)";
    size_t n = 0;

    // CFX_ASSERT(cfx_big_scan_num_n(&b, (const uint8_t*)s, strlen(s), &n) == 0);
    // CFX_ASSERT(n == 5); /* "0o377" */

    // CFX_ASSERT(cfx_big_from_str(&ref, "255") == 0);
    // CFX_ASSERT(cfx_big_eq(&b, &ref));
}

static void test_scan_num_dec_basic(void) {
    cfx_big_t b, ref;
    cfx_big_init(&b);
    cfx_big_init(&ref);

    const char *s = "12345*9";
    size_t n = 0;

    CFX_ASSERT(cfx_big_scan_num_n(&b, (const uint8_t*)s, strlen(s), &n) == 0);
    CFX_ASSERT(n == 5); /* "12345" */

    CFX_ASSERT(cfx_big_from_str(&ref, "12345") == 0);
    CFX_ASSERT(cfx_big_eq(&b, &ref));
}

static void test_scan_num_stops_at_invalid_digit(void) {
    cfx_big_t b, ref;
    cfx_big_init(&b);
    cfx_big_init(&ref);

    const char *s = "0x12zz";
    size_t n = 0;

    CFX_ASSERT(cfx_big_scan_num_n(&b, (const uint8_t*)s, strlen(s), &n) == 0);
    CFX_ASSERT(n == 4); /* "0x12" */

    CFX_ASSERT(cfx_big_from_str(&ref, "18") == 0);
    CFX_ASSERT(cfx_big_eq(&b, &ref));
}

static void test_scan_num_prefix_only_fails(void) {
    cfx_big_t b;
    cfx_big_init(&b);

    size_t n = 777;

    CFX_ASSERT(cfx_big_scan_num_n(&b, (const uint8_t*)"0x", 2, &n) < 0);
    CFX_ASSERT(n == 0);

    n = 777;
    CFX_ASSERT(cfx_big_scan_num_n(&b, (const uint8_t*)"0b", 2, &n) < 0);
    CFX_ASSERT(n == 0);

    n = 777;
    // CFX_ASSERT(cfx_big_scan_num_n(&b, (const uint8_t*)"0o", 2, &n) < 0);
    // CFX_ASSERT(n == 0);
}

static void test_scan_num_non_number_fails(void) {
    cfx_big_t b;
    cfx_big_init(&b);

    size_t n = 777;
    CFX_ASSERT(cfx_big_scan_num_n(&b, (const uint8_t*)"+123", 4, &n) < 0);
    CFX_ASSERT(n == 0);

    n = 777;
    CFX_ASSERT(cfx_big_scan_num_n(&b, (const uint8_t*)"(", 1, &n) < 0);
    CFX_ASSERT(n == 0);

    n = 777;
    CFX_ASSERT(cfx_big_scan_num_n(&b, (const uint8_t*)" 123", 4, &n) < 0); /* scanner does not skip ws */
    CFX_ASSERT(n == 0);
}

static void test_scan_num_prefix_case_insensitive(void) {
    cfx_big_t a, b, c;
    cfx_big_init(&a);
    cfx_big_init(&b);
    cfx_big_init(&c);

    size_t na = 0, nb = 0, nc = 0;

    CFX_ASSERT(cfx_big_scan_num_n(&a, (const uint8_t*)"0XFF ", 5, &na) == 0);
    CFX_ASSERT(na == 4);

    CFX_ASSERT(cfx_big_scan_num_n(&b, (const uint8_t*)"0B1010", 6, &nb) == 0);
    CFX_ASSERT(nb == 6);

    // CFX_ASSERT(cfx_big_scan_num_n(&c, (const uint8_t*)"0O77", 4, &nc) == 0);
    // CFX_ASSERT(nc == 4);
}

static void test_scan_num_respects_in_len(void) {
    cfx_big_t b, ref;
    cfx_big_init(&b);
    cfx_big_init(&ref);

    const char *s = "0x1234+999";
    size_t n = 0;

    /* Only give it "0x12" */
    CFX_ASSERT(cfx_big_scan_num_n(&b, (const uint8_t*)s, 4, &n) == 0);
    CFX_ASSERT(n == 4);

    CFX_ASSERT(cfx_big_from_str(&ref, "18") == 0); /* 0x12 */
    CFX_ASSERT(cfx_big_eq(&b, &ref));
}



/* --- hex tests --- */
static void test_hex_zero_empty_n(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    /* n == 0 should yield "0" */
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
    /* 0x1 -> "1" */
    cfx_limb_t limbs1[] = {UINT64_C(0x1)};
    cfx_big_from_limbs(&b, limbs1, 1);
    size_t len1 = 0;
    char* s1 = cfx_big_to_hex(&b, &len1);
    CFX_ASSERT(strcmp(s1, "1") == 0);
    CFX_ASSERT(len1 == strlen("1"));
    CFX_ASSERT(s1[len1] == '\0');
    printf("[%s] - %s\n", __func__, s1);
    free(s1);

    /* 0xabcdef -> "abcdef" (lowercase) */
    b.limb[0] = UINT64_C(0xabcdef);
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
    /* 2^60 = 0x1000000000000000 -> 1 followed by 15 zeros (16 digits total) */
    uint64_t v = UINT64_C(0x1000000000000000);
    cfx_big_from_u64(&b, v);
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
    /* high = 0x1, low = 0x1 -> "1" + "0000000000000001" */
    cfx_limb_t limbs[2];
    limbs[0] = (cfx_limb_t)1; /* low */
    limbs[1] = (cfx_limb_t)1; /* high */
    cfx_big_from_limbs(&b, limbs, 2);
    size_t len = 0;
    char* s = cfx_big_to_hex(&b, &len);
    #if (CFX_LIMB_BITS == 64)
    CFX_ASSERT(strcmp(s, "10000000000000001") == 0);
    CFX_ASSERT(len == 17);                  /* 1 + 16 = CFX_LIMB_DIGITS_DEC */
    #elif (CFX_LIMB_BITS == 32)
    CFX_ASSERT(strcmp(s, "100000001") == 0);
    CFX_ASSERT(len == 9);                   /* 1 + 8 hex digits */
    #endif
    CFX_ASSERT(s[len] == '\0');
    cfx_big_free(&b);
    printf("[%s] - %s\n", __func__, s);
    free(s);
    PRINT_TEST(1);
}

static void test_hex_two_limbs_mixed_digits(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    /* high = 0xABC, low = 0x0011223344556677 */
    /* expect: "abc" + "0011223344556677" */
    cfx_limb_t limbs[2];
    char* s = NULL;
    size_t szout = 0;
#if (CFX_LIMB_BITS==64)
    limbs[0] = UINT64_C(0x0011223344556677);    /* low */
    limbs[1] = UINT64_C(0x0000000000000ABC);    /* high */
    cfx_big_from_limbs(&b, limbs, 2);
    s = cfx_big_to_hex(&b, &szout);
    CFX_ASSERT(strcmp(s, "abc0011223344556677") == 0);
#else
    limbs[0] = UINT32_C(0x00112233);    /* low */
    limbs[1] = UINT32_C(0x00000ABC);    /* high */
    cfx_big_from_limbs(&b, limbs, 2);
    s = cfx_big_to_hex(&b, &szout);
    CFX_ASSERT(strcmp(s, "abc00112233") == 0);
#endif
    CFX_ASSERT(s != NULL);
    CFX_ASSERT(strlen(s) == szout);
    CFX_ASSERT(s[szout] == '\0');
    printf("[%s] - %s\n", __func__, s);
    cfx_big_free(&b);
    free(s);
    PRINT_TEST(1);
}

static void test_hex_leading_zero_limb_skipped(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    /* limbs: [low=X, mid=Y, high=0] -> should behave like just [low=X, mid=Y] */
    cfx_limb_t limbs3[3];
#if CFX_LIMB_BITS == 64
    limbs3[0] = UINT64_C(0xDEADBEEFCAFEBABE); /* low */
    limbs3[1] = UINT64_C(0x0000000000000123); /* mid */
    limbs3[2] = UINT64_C(0x0000000000000000); /* high (zero) */
#elif CFX_LIMB_BITS == 32
    limbs3[0] = UINT64_C(0xDEADBEEF); /* low */
    limbs3[1] = UINT64_C(0x00000123); /* mid */
    limbs3[2] = UINT64_C(0x00000000); /* high */
#endif
    cfx_big_from_limbs(&b, limbs3, 3);
    size_t len = 0;
    char* s = cfx_big_to_hex(&b, &len);
    /* expected: "123" + "%016" of low */
#if CFX_LIMB_BITS == 64
    const char* expect_low = "deadbeefcafebabe"; /* lowercase */
    /* const char* expect = "1230000000000000" "deadbeefcafebabe"; */
    /* But careful: mid=0x123 -> "123", low padded to 16 */
    char expect_buf[3 + 16 + 1];
#elif CFX_LIMB_BITS == 32
    const char* expect_low = "deadbeef";
    char expect_buf[3 + 8 + 1];
#endif
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
    /* high = 0x00000000000000ab -> "ab", low zero-padded */
    cfx_limb_t limbs[] = {
        (cfx_limb_t)0,
        (cfx_limb_t)0xAB
    };
    cfx_big_from_limbs(&b, limbs, 2);
    size_t len = 0;
    char* s = cfx_big_to_hex(&b, &len);
#if CFX_LIMB_BITS == 64
    CFX_ASSERT(strcmp(s, "ab0000000000000000") == 0);
    CFX_ASSERT(len == 2 + 16);
#else
    CFX_ASSERT(strcmp(s, "ab00000000") == 0);
    CFX_ASSERT(len == 2 + 8);
#endif
    CFX_ASSERT(s[len] == '\0');
    cfx_big_free(&b);
    free(s);
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
    cfx_big_from_limb(&m, 987);

    cfx_big_mul_eq(&b, &m);

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
    cfx_big_from_limb(&b, 5);
    cfx_big_init(&m);
    cfx_big_from_limb(&m, 7);

    cfx_big_mul_eq(&b, &m);
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
    cfx_big_from_limb(&a, 0);
    cfx_big_init(&m);
    cfx_big_from_limb(&m, 1);
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
#if CFX_LIMB_BITS == 64
    cfx_limb_t limbs_b[] = {UINT64_C(0xFFFFFFFFFFFFFFFF), UINT64_C(0xFFFFFFFFFFFFFFFF)};
#elif CFX_LIMB_BITS == 32
    cfx_limb_t limbs_b[] = {UINT32_C(0xFFFFFFFF), UINT32_C(0xFFFFFFFF)};
#endif
    cfx_big_from_limbs(&b, limbs_b, 2);
    cfx_big_init(&m);
    cfx_big_from_limb(&m, 2);

    cfx_big_mul_eq(&b, &m);
#if CFX_LIMB_BITS == 64
    cfx_limb_t expect[] = {UINT64_C(0xFFFFFFFFFFFFFFFE), UINT64_C(0xFFFFFFFFFFFFFFFF), UINT64_C(1)};
#elif CFX_LIMB_BITS == 32
    cfx_limb_t expect[] = {UINT32_C(0xFFFFFFFE), UINT32_C(0xFFFFFFFF), UINT32_C(1)};
#endif
    big_expect_limbs(__func__, &b, expect, 3);

    cfx_big_free(&b);
    cfx_big_free(&m);
    PRINT_TEST(1);
}


static void test_mul_by_base_2_64_shift(void) {
    /* Multiply by 2^64 (limbs = [0,1]) should shift by one limb */
    cfx_big_t b, m;
    cfx_big_init(&b);
    cfx_big_init(&m);
#if CFX_LIMB_BITS == 64
    cfx_limb_t limbs_b[] = {UINT64_C(0x0123456789ABCDEF), UINT64_C(0x0FEDCBA987654321), UINT64_C(0x0000000000000001)};
#elif CFX_LIMB_BITS == 32
    cfx_limb_t limbs_b[] = {UINT32_C(0x01234567), UINT32_C(0x0FEDCBA9), UINT32_C(0x00000001)};
#endif
    size_t sz0 = sizeof(limbs_b)/sizeof(limbs_b[0]);
    cfx_big_from_limbs(&b, limbs_b, sz0);
#if CFX_LIMB_BITS == 64
    cfx_limb_t limbs_m[] = {UINT64_C(0), UINT64_C(1)};
#elif CFX_LIMB_BITS == 32
    cfx_limb_t limbs_m[] = {UINT32_C(0), UINT32_C(1)};
#endif
    cfx_big_from_limbs(&m, limbs_m, 2);

    const size_t N = 1111;

    for (size_t n = 1; n < N; ++n) {

        cfx_big_mul_eq(&b, &m);
        cfx_limb_t* expect = (cfx_limb_t*)malloc((n + sz0)*sizeof(cfx_limb_t));
        for (size_t k = 0; k < n; ++k) {
            expect[k] = 0;
        }
#if CFX_LIMB_BITS == 64
        expect[n+sz0-3] = UINT64_C(0x0123456789ABCDEF),
        expect[n+sz0-2] = UINT64_C(0x0FEDCBA987654321);
        expect[n+sz0-1] = UINT64_C(0x0000000000000001);
#elif CFX_LIMB_BITS == 32
        expect[n+sz0-3] = UINT32_C(0x01234567),
        expect[n+sz0-2] = UINT32_C(0x0FEDCBA9);
        expect[n+sz0-1] = UINT32_C(0x00000001);
#endif

        big_expect_limbs(__func__, &b, expect, sz0 + n);
        /* printf("[test_mul_by_base_2_64_shift]: tested shift of %zu OK\n", n); */
        free(expect);
        expect = NULL;
    }

    cfx_big_free(&b);
    cfx_big_free(&m);
    PRINT_TEST(1);
}

static void test_self_multiply_square(void) {
    /* (2^64 - 1)^2 = [0xFFFFFFFFFFFFFFFE, 1] in base 2^64 */
    cfx_big_t b;
    cfx_big_init(&b);
#if CFX_LIMB_BITS == 64
    cfx_limb_t limbs[] = {UINT64_C(0xFFFFFFFFFFFFFFFF)};
#elif CFX_LIMB_BITS == 32
    cfx_limb_t limbs[] = {UINT32_C(0xFFFFFFFF)};
#endif

    cfx_big_from_limbs(&b, limbs, 1);
    big_expect_limbs(__func__, &b, limbs, 1);
    cfx_big_mul_eq(&b, &b); /* self-mul path */
#if CFX_LIMB_BITS == 64
    cfx_limb_t expect[] = {UINT64_C(1), UINT64_C(0xFFFFFFFFFFFFFFFE)};
#elif CFX_LIMB_BITS == 32
    cfx_limb_t expect[] = {UINT32_C(1), UINT32_C(0xFFFFFFFE)};
#endif

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
    cfx_big_mul_eq(&b, &b);
    printf("after: "); CFX_BIG_PRINT(b);
    cfx_big_free(&b);
    free(limbs);
    limbs = NULL;
    PRINT_TEST(1);
}


static void test_known_squares(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    cfx_big_from_dec(&b,
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
    cfx_big_from_dec(&b,
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
    cfx_big_from_dec(&b,
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
    cfx_big_from_dec(&b,
        "4946608029462090681478206578991795742708644"
        "2564742658586426229545514803499564697000372"
        "3095350971345437292114654548843072761868784"
        "674125049315509629339381027496416991005114370");
    cfx_big_sq(&b); /* 1 */
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
    cfx_big_sq(&b); /* 2 */
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

    cfx_big_sq(&b); /* 4 */
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
    /* CFX_BIG_PRINT_LIMBS(b); */


    /* cfx_big_sq(&b); // 8 */
    cfx_big_mul_eq(&b, &b);
    free(s);
    s = cfx_big_to_str(&b, NULL);
    /*// sanity check: */
    cfx_big_t B;
    cfx_big_init(&B);
    cfx_big_from_dec(&B, s);
    CFX_ASSERT(cfx_big_eq(&B, &b));
    char* sanity = cfx_big_to_str(&B, NULL);
    CFX_ASSERT(strcmp(s, sanity) == 0);
    expect = "12850498340565118640831124663943566637163477270761965804775357209550765230785703949035729472531170991846950765062231843931603991577837122702202578619413163343802199724298959199688594791894654764450428039375530903550022460719636135025802643105174385371702288121193336984280156229260813405468551910887207681430170841269088143925314707610300410631962624300469472858167288519605269872627946462373021700902442150967139315116172545474219540016600664971616358393885751109046842533050183826430384246617880667515303122705936214635228967259835005388450384900687056571344933537245392429971062926387006758181202634581358151550436165392630946855812938662698021351517729227827971635271649805374176176940281753018739171849301949232131565817532850243174625134450254366859930465467865935871669065946539235916390634180117949435696071838514048795111101267436279103470592457465410839430058374274997666026789906916220271637971862818408213637491993426120588003613013514121426918778425609885983132067309940493824037246336471087139347566932289527147579507444510833470558987620105627392491106342931918229603332167519231845701934759523730028587564778522689440283636622636742434431091107839566802176801151089094414637433966684982281939928036568083598873318392485777034090015016265841172493029035363177944812390618909768667573763207436642271055957361540865636837910531318607009603088734768184385262471743568922118461465472820244894124358177759866070832621634448723670657830089717393358547747703735801523850299292208462216521430657525373415183422463571048104152842920726565698681208502981595141456253796350916896570590072884174607249064969441514018160468583314653443352671437479825779995459466312201730735736110098653787331865116531446616641987026900818215864623614884767160893264417488979936319828793408610211926232514931040704633236494245315339349224851289392387332474248568479216841581124860328079900988652092392614618200827173162666704740632419460637566175022280383967839231370086128119069239566679344622785639260673870367122861560404871658960189532546990946727895958575281330314215816959430714079369740156022375856316371677709024505475091778006456248457155666092229065835357913477405119329372494342721686545800319902824829662799906780886313825440342423658177641941252419690985063319979924638260305722020541317833802446921285816583681614997935512759860955018042317788791616381896041300997716291956236555872471276333548522003017903198166576927518763140985360706323427487549508484530837381768891367647567673106372095645829146242426938921111540219157093806251149323044120727645913143734374450245354473392643178616466237265916832878586939031607485456329408461025836997224109513396266781143113993671870210102616033541932236431301110708951036119790977069093536028698906465137616493507440201974410000000000000000";
    cfx_big_from_dec(&B, expect);
    free(sanity);
    sanity = cfx_big_to_str(&B, NULL);
    CFX_ASSERT(strcmp(sanity, expect) == 0);

    for (size_t i = 0; i < b.n; ++i) {
        printf("calculated: b.limb[%zu]: " CFX_PRIuLIMB "; correct: B.limb[%zu]: " CFX_PRIuLIMB ": %s: diff %d\n",
        i, b.limb[i], i, B.limb[i], b.limb[i] == B.limb[i] ? "ok" : "--- NOT OK",
        (int)(b.limb[i] & 0xFFFF) - (int)(B.limb[i] & 0xFFFF));
    }
    printf("\n\n");
    CFX_ASSERT(cfx_big_eq(&B, &b));
    CFX_ASSERT(strcmp(s, expect) == 0);
    int cnt = 0;

    cfx_big_mul_eq(&b, &b); /* 16 */
    printf("mul %d len: %zu \n", ++cnt, b.n);

    cfx_big_mul_eq(&b, &b); /* 32 */
    printf("mul %d len: %zu \n", ++cnt, b.n);
    cfx_big_mul_eq(&b, &b); /* 64 */
    printf("mul %d len: %zu \n", ++cnt, b.n);
    cfx_big_mul_eq(&b, &b); /* 128 */
    printf("mul %d len: %zu \n", ++cnt, b.n);
    cfx_big_mul_eq(&b, &b); /* 256 */
    printf("mul %d len: %zu \n", ++cnt, b.n);
    cfx_big_mul_eq(&b, &b); /* 512 */
    printf("mul %d len: %zu \n", ++cnt, b.n);
    /* cfx_big_mul_eq(&b, &b); // 1024 */
    /* printf("mul %d len: %zu \n", ++cnt, b.n); */
    /* cfx_big_mul_eq(&b, &b); // 2048 */
    /* printf("mul %d len: %zu \n", ++cnt, b.n); */
    /* cfx_big_mul_eq(&b, &b); // 4096 */
    /* printf("mul %d len: %zu \n", ++cnt, b.n); */
    /* cfx_big_mul_eq(&b, &b); // 8192 */
    /* printf("mul %d len: %zu \n", ++cnt, b.n); */

    /* for (size_t i = 0; i < b.n; ++i) { */
    /*     printf("calculated: b.limb[%zu]: " CFX_PRIuLIMB " (0x"CFX_PRI0xLIMB")\n", i, b.limb[i], b.limb[i]); */
    /* } */

    char* huge = cfx_big_to_hex(&b, NULL);
    /* printf("%s\n", huge); */
    printf("digits: %zu\n", strlen(huge));
    free(huge);
    free(s);
    cfx_big_free(&b);
    PRINT_TEST(1);
}

static void test_known_squares_2(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    cfx_big_from_dec(&b,
        "4946608029462090681478206578991795742708644"
        "2564742658586426229545514803499564697000372"
        "3095350971345437292114654548843072761868784"
        "674125049315509629339381027496416991005114370");
    cfx_big_mul_eq_csa(&b, &b); /* 1 */
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
    cfx_big_mul_eq_csa(&b, &b); /* 2 */
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

    cfx_big_mul_eq_csa(&b, &b); /* 4 */
    free(s);
    s = cfx_big_to_str(&b, NULL);
    /* write_string_wrapped(s, "", 80); */
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
    /* CFX_BIG_PRINT_LIMBS(b); */


    cfx_big_mul_eq(&b, &b); /* 8 */
    free(s);
    s = cfx_big_to_str(&b, NULL);

    /* ---------- sanity check: ---------------------------------- */
    #if 0
    cfx_big_t B;
    cfx_big_init(&B);
    cfx_big_from_dec(&B, s);
    CFX_ASSERT(cfx_big_eq(&B, &b));
    char* sanity = cfx_big_to_str(&B, NULL);
    CFX_ASSERT(strcmp(s, sanity) == 0);
    expect = "12850498340565118640831124663943566637163477270761965804775357209550765230785703949035729472531170991846950765062231843931603991577837122702202578619413163343802199724298959199688594791894654764450428039375530903550022460719636135025802643105174385371702288121193336984280156229260813405468551910887207681430170841269088143925314707610300410631962624300469472858167288519605269872627946462373021700902442150967139315116172545474219540016600664971616358393885751109046842533050183826430384246617880667515303122705936214635228967259835005388450384900687056571344933537245392429971062926387006758181202634581358151550436165392630946855812938662698021351517729227827971635271649805374176176940281753018739171849301949232131565817532850243174625134450254366859930465467865935871669065946539235916390634180117949435696071838514048795111101267436279103470592457465410839430058374274997666026789906916220271637971862818408213637491993426120588003613013514121426918778425609885983132067309940493824037246336471087139347566932289527147579507444510833470558987620105627392491106342931918229603332167519231845701934759523730028587564778522689440283636622636742434431091107839566802176801151089094414637433966684982281939928036568083598873318392485777034090015016265841172493029035363177944812390618909768667573763207436642271055957361540865636837910531318607009603088734768184385262471743568922118461465472820244894124358177759866070832621634448723670657830089717393358547747703735801523850299292208462216521430657525373415183422463571048104152842920726565698681208502981595141456253796350916896570590072884174607249064969441514018160468583314653443352671437479825779995459466312201730735736110098653787331865116531446616641987026900818215864623614884767160893264417488979936319828793408610211926232514931040704633236494245315339349224851289392387332474248568479216841581124860328079900988652092392614618200827173162666704740632419460637566175022280383967839231370086128119069239566679344622785639260673870367122861560404871658960189532546990946727895958575281330314215816959430714079369740156022375856316371677709024505475091778006456248457155666092229065835357913477405119329372494342721686545800319902824829662799906780886313825440342423658177641941252419690985063319979924638260305722020541317833802446921285816583681614997935512759860955018042317788791616381896041300997716291956236555872471276333548522003017903198166576927518763140985360706323427487549508484530837381768891367647567673106372095645829146242426938921111540219157093806251149323044120727645913143734374450245354473392643178616466237265916832878586939031607485456329408461025836997224109513396266781143113993671870210102616033541932236431301110708951036119790977069093536028698906465137616493507440201974410000000000000000";
    cfx_big_from_dec(&B, expect);
    free(sanity);
    sanity = cfx_big_to_str(&B, NULL);
    CFX_ASSERT(strcmp(sanity, expect) == 0);

    for (size_t i = 0; i < b.n; ++i) {
        printf("calculated: b.limb[%zu]: " CFX_PRIuLIMB "; correct: B.limb[%zu]: " CFX_PRIuLIMB ": %s: diff %d\n",
        i, b.limb[i], i, B.limb[i], b.limb[i] == B.limb[i] ? "ok" : "--- NOT OK",
        (int)(b.limb[i] & 0xFFFF) - (int)(B.limb[i] & 0xFFFF));
    }
    printf("\n\n");
    CFX_ASSERT(cfx_big_eq(&B, &b));
    CFX_ASSERT(strcmp(s, expect) == 0);
    #endif
    /* ----------------------------------------------------------- */

    int cnt = 0;
    cfx_big_mul_eq_csa(&b, &b); /* 16 */
    printf("mul csa %d len: %zu \n", ++cnt, b.n);
    cfx_big_mul_eq_csa(&b, &b); /* 32 */
    printf("mul csa %d len: %zu \n", ++cnt, b.n);
    cfx_big_mul_eq_csa(&b, &b); /* 64 */
    printf("mul csa %d len: %zu \n", ++cnt, b.n);
    cfx_big_mul_eq_csa(&b, &b); /* 128 */
    printf("mul csa %d len: %zu \n", ++cnt, b.n);
    cfx_big_mul_eq_csa(&b, &b); /* 256 */
    printf("mul csa %d len: %zu \n", ++cnt, b.n);
    cfx_big_mul_eq_csa(&b, &b); /* 512 */
    printf("mul csa %d len: %zu \n", ++cnt, b.n);
    /* cfx_big_mul_eq_csa(&b, &b); // 1024 */
    /* printf("mul csa %d len: %zu \n", ++cnt, b.n); */
    /* cfx_big_mul_eq_csa(&b, &b); // 2048 */
    /* printf("mul csa %d len: %zu \n", ++cnt, b.n); */
    /* cfx_big_mul_eq_csa(&b, &b); // 4096 */
    /* printf("mul csa %d len: %zu \n", ++cnt, b.n); */

    /* for (size_t i = 0; i < b.n; ++i) { */
    /*     printf("calculated: b.limb[%zu]: " CFX_PRIuLIMB " (0x"CFX_PRI0xLIMB")\n", i, b.limb[i], b.limb[i]); */
    /* } */

    char* huge = cfx_big_to_hex(&b, NULL);
    /* printf("%s\n", huge); */
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
    cfx_big_from_dec(&tmp, dec);
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
    /* rely on library ops to keep canonical form later */
}

/* Property check: n == q*d + r and r < d */
static void assert_n_eq_qd_plus_r(const cfx_big_t* n, const cfx_big_t* q,
                                  const cfx_big_t* d, const cfx_big_t* r) {
    cfx_big_t check;
    cfx_big_init(&check);
    cfx_big_copy(&check, q);
    cfx_big_t tmp;
    cfx_big_init(&tmp);

    cfx_big_mul_eq(&check, d);
    cfx_big_copy(&tmp, r);
    cfx_big_add(&check, &tmp);

    CFX_ASSERT(cfx_big_eq(&check, n));
    CFX_ASSERT(cfx_big_cmp(r, d) == -1);

    cfx_big_free(&check);
    cfx_big_free(&tmp);
}

/* ---------- tests ---------- */

static void test_big_div_divide_by_zero(void) {
    cfx_big_t n, d, q, r;
    cfx_big_init(&n);
    cfx_big_init(&d);
    cfx_big_init(&q);
    cfx_big_init(&r);
    cfx_big_from_limb(&n, 123);
    cfx_big_from_limb(&d, 0);

    int rc = cfx_big_divrem(&q, &r, &n, &d);
    CFX_ASSERT(rc == -1);

    cfx_big_free(&n);
    cfx_big_free(&d);
    cfx_big_free(&q);
    cfx_big_free(&r);
}

static void test_big_div_zero_dividend(void) {
    cfx_big_t n, d, q, r;
    cfx_big_init(&n);
    cfx_big_init(&d);
    cfx_big_init(&q);
    cfx_big_init(&r);
    cfx_big_from_limb(&n, 0);
    cfx_big_from_limb(&d, 42);

    int rc = cfx_big_divrem(&q, &r, &n, &d);
    CFX_ASSERT(rc == 0);
    expect_dec_eq(&q, "0");
    expect_dec_eq(&r, "0");

    cfx_big_free(&n);
    cfx_big_free(&d);
    cfx_big_free(&q);
    cfx_big_free(&r);
}

static void test_big_div_n_less_than_d(void) {
    cfx_big_t n, d, q, r;
    cfx_big_init(&n);
    cfx_big_init(&d);
    cfx_big_init(&q);
    cfx_big_init(&r);
    cfx_big_from_limb(&n, 123456);
    cfx_big_from_limb(&d, 123456789);

    int rc = cfx_big_divrem(&q, &r, &n, &d);
    CFX_ASSERT(rc == 0);
    expect_dec_eq(&q, "0");
    CFX_ASSERT(cfx_big_eq(&r, &n));

    cfx_big_free(&n);
    cfx_big_free(&d);
    cfx_big_free(&q);
    cfx_big_free(&r);
}

static void test_big_div_equal_numbers(void) {
    cfx_big_t n, d, q, r;
    cfx_big_init(&n);
    cfx_big_init(&d);
    cfx_big_init(&q);
    cfx_big_init(&r);
    cfx_big_from_dec(&n, "1234567890123456789012345678901234567890");
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

static void test_big_div_single_limb_divisor_property(void) {
    cfx_big_t n, d, q, r;
    cfx_big_init(&n);
    cfx_big_init(&d);
    cfx_big_init(&q);
    cfx_big_init(&r);
    cfx_big_from_dec(&n, "340282366920938463463374607431768211455"); /* 2^128 - 1 */
    cfx_big_from_limb(&d, UINT64_C(123456789));

    int rc = cfx_big_divrem(&q, &r, &n, &d);
    CFX_ASSERT(rc == 0);
    assert_n_eq_qd_plus_r(&n, &q, &d, &r);

    cfx_big_free(&n);
    cfx_big_free(&d);
    cfx_big_free(&q);
    cfx_big_free(&r);
}

static void test_big_div_multi_limb_divisor_exact_and_remainder(void) {
    /* n = a*b + r, then n / b -> q=a, rem=r */
    cfx_big_t a, b, r, n, q, rem;
    cfx_big_init(&a);
    cfx_big_init(&b);
    cfx_big_init(&n);
    cfx_big_init(&q);
    cfx_big_init(&r);
    cfx_big_init(&rem);
    cfx_big_from_dec(&a, "123456789012345678901234567890123456789");
    cfx_big_from_dec(&b, "987654321098765432109876543210987654321");
    cfx_big_from_dec(&r, "12345678901234567890");

    cfx_big_copy(&n, &a);
    cfx_big_mul_eq(&n, &b);
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

static void test_big_div_in_place_eq_with_remainder(void) {
    /* Build n = a*b + 42, then n := n / b, rem = 42 */
    cfx_big_t a, b, n, rem, forty_two;
    cfx_big_init(&a);
    cfx_big_init(&b);
    cfx_big_init(&n);
    cfx_big_init(&rem);
    cfx_big_init(&forty_two);
    cfx_big_from_dec(&a, "1122334455667788990011223344556677889900");
    cfx_big_from_dec(&b, "18446744073709551616"); /* 2^64 */
    cfx_big_copy(&n, &a);
    cfx_big_mul_eq(&n, &b);
    cfx_big_from_limb(&forty_two, 42);
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

static void test_big_div_quotient_only_and_remainder_only(void) {
    cfx_big_t a, b, n, q, r, b_minus_1;
    cfx_big_init(&a);
    cfx_big_init(&b);
    cfx_big_init(&n);
    cfx_big_init(&q);
    cfx_big_init(&r);
    cfx_big_init(&b_minus_1);
    cfx_big_from_dec(&a, "3141592653589793238462643383279502884197");
    cfx_big_from_dec(&b, "2718281828459045235360287471352662497757");

    /* n = a*b + (b-1) */
    cfx_big_copy(&n, &a);
    cfx_big_mul_eq(&n, &b);

    cfx_big_copy(&b_minus_1, &b);
    big_dec1(&b_minus_1);

    cfx_big_add(&n, &b_minus_1);

    /* quotient only */
    CFX_ASSERT(cfx_big_div_out(&q, &n, &b) == 0);
    CFX_ASSERT(cfx_big_eq(&q, &a));

    /* remainder only */
    CFX_ASSERT(cfx_big_mod(&r, &n, &b) == 0);
    CFX_ASSERT(cfx_big_eq(&r, &b_minus_1));

    cfx_big_free(&a);
    cfx_big_free(&b);
    cfx_big_free(&n);
    cfx_big_free(&q);
    cfx_big_free(&r);
    cfx_big_free(&b_minus_1);
}

static void test_big_div_alias_remainder_eq_src(void) {
    CFX_ASSERT(1);
    /* TODO! */
    #if 0
    /* Verify cfx_big_div_eq supports r == b */
    cfx_big_t b, d;
    cfx_big_init(&b);
    cfx_big_init(&d);
    cfx_big_from_dec(&b, "123456789012345678901234567890");
    cfx_big_from_dec(&d, "987654321");

    cfx_big_t orig;
    cfx_big_init(&orig);

    cfx_big_copy(&orig, &b);

    int rc = cfx_big_div_eq(&b, &d, &b);   /* r aliases src */
    CFX_ASSERT(rc == 0);

    /* Check property with a fresh recompute: orig == q*d + r, where q is in 'b' (after div) */
    cfx_big_t q_copy;
    cfx_big_init(&q_copy);
    cfx_big_copy(&q_copy, &b);  /* b now holds q */
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

    cfx_limb_t l0 = b.limb[0];
    cfx_limb_t l1 = b.limb[1];

    char* s = cfx_big_to_bin(&b, NULL);
    printf("---------------> B4RE %s\n", s);
    cfx_big_shl_bits(&b, &b, CFX_LIMB_BITS+1);
    free(s);
    s = cfx_big_to_bin(&b, NULL);
    printf("---------------> AFTR %s\n", s);
    free(s);

    CFX_ASSERT(b.limb[2] == l1*2);
    CFX_ASSERT(b.limb[1] == l0*2);
    CFX_ASSERT(b.limb[0] == 0);

    l0 = b.limb[0];
    l1 = b.limb[1];

    cfx_big_shr_bits(&b, &b, 1);
    CFX_ASSERT(b.limb[1] == l1/2);
    CFX_ASSERT(b.limb[0] == l0/2);
}

/* ---- helpers ----------------------------------------------------------- */

static void assert_hex_eq(const char* tag, const cfx_big_t* x, const char* hex_exp) {
    char* got = cfx_big_to_hex(x, NULL);
    int ok = (strcmp(got, hex_exp) == 0);
    if (!ok) {
        fprintf(stderr, "[%s] expected 0x%s, got 0x%s\n", tag, hex_exp, got);
        /* fflush(stderr); */
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
#if CFX_LIMB_BITS == 64
    CHECK_SHL_CASE("8000000000000000", 1, "10000000000000000");    /* 2^63 << 1 -> 2^64 */
    CHECK_SHL_CASE("1", 64, "10000000000000000");                   /* append 16 hex zeros */
    CHECK_SHL_CASE("1", 68, "100000000000000000");                 /* 16 zeros + 1 hex nibble */
#else
    CHECK_SHL_CASE("80000000", 1, "100000000");                    /* 2^31 << 1 -> 2^32 */
    CHECK_SHL_CASE("1", 32, "100000000");                          /* append 8 hex zeros */
    CHECK_SHL_CASE("1", 36, "1000000000");                         /* 8 zeros + 1 hex nibble */
#endif
}

static void test_shl_cross_limb_1bit(void) {
    /* two-limb pattern: hi=1, lo=0x8000...0001, shift left by 1
       expected: hi' = 3, lo' = 2  => hex "3" || "000...0002" */
#if CFX_LIMB_BITS == 64
    CHECK_SHL_CASE("100000000000000008000000000000001", 1, "200000000000000010000000000000002");
#else
    CHECK_SHL_CASE("180000001", 1, "300000002");  /* 32-bit: hi=1, lo=0x80000001 */
#endif
}

static void test_shr_basic_identity(void) {
    /* shift by 0: identity */
    CHECK_SHR_CASE("0", 0, "0");
    CHECK_SHR_CASE("1", 0, "1");
    CHECK_SHR_CASE("deadbeef", 0, "deadbeef");
}

static void test_shr_drop_whole_limb(void) {
    /* exact limb-width drops */
#if CFX_LIMB_BITS == 64
    CHECK_SHR_CASE("10000000000000000", 64, "1");         /* 2^64 >> 64 */
    CHECK_SHR_CASE("100000000000000000", 68, "1");       /* (2^68) >> 68 */
#else
    CHECK_SHR_CASE("100000000", 32, "1");                 /* 2^32 >> 32 */
    CHECK_SHR_CASE("1000000000", 36, "1");               /* (2^36) >> 36 */
#endif
}

static void test_shr_to_zero(void) {
    /* shifting past total bit-length -> zero */
#if CFX_LIMB_BITS == 64
    CHECK_SHR_CASE("1", 65, "0");
    CHECK_SHR_CASE("ffffffffffffffff", 128, "0");
#else
    CHECK_SHR_CASE("1", 33, "0");
    CHECK_SHR_CASE("ffffffff", 64, "0");
#endif
}

static void test_shr_cross_limb_carry_6bits(void) {
#if CFX_LIMB_BITS == 64
    /*
       Input limbs: hi=0x3e, lo=0x0c40  -> hex "3e0000000000000c40"
       Right shift by 6:
         r1 = 0x3e >> 6 = 0
         carry = (0x3e & 0x3f) << 58 = 0xf800000000000000
         r0 = (0x0c40 >> 6) | carry = 0x31 | 0xf800000000000000 = 0xf800000000000031
       Trim top zero limb => "f800000000000031"
    */
    CHECK_SHR_CASE("3e0000000000000c40", 6, "f800000000000031");
#else
    /* 32-bit version:
       Input limbs: hi=0x3e, lo=0x0c40  -> hex "3e00000c40"
       Right shift by 6:
         carry = (0x3e & 0x3f) << 26 = 0xf8000000
         r0 = (0x0c40 >> 6) | carry = 0x31 | 0xf8000000 = 0xf8000031
       Trim top zero limb => "f8000031"
    */
    CHECK_SHR_CASE("3e00000c40", 6, "f8000031");
#endif
}

static void test_shr_mixed_cases(void) {
    /* assorted mixes with odd s and multiple limbs */
#if CFX_LIMB_BITS == 64
    CHECK_SHR_CASE("100000000000000000000", 4, "10000000000000000000");  /* divide by 16 */
    CHECK_SHR_CASE("abcdef0123456789", 4, "abcdef012345678");          /* >>4 removes low nibble */
    CHECK_SHR_CASE("abcdef0123456789", 60, "a");                       /* >>60 keeps top nibble */
#else
    CHECK_SHR_CASE("10000000000", 4, "1000000000");                    /* divide by 16 */
    CHECK_SHR_CASE("abcdef01", 4, "abcdef0");                          /* >>4 removes low nibble */
    CHECK_SHR_CASE("abcdef01", 28, "a");                               /* >>28 keeps top nibble */
#endif
}

/* ------------------------------------------------------------------ */
static void expect_limb_pattern(const cfx_big_t* a, size_t n,
                                cfx_limb_t l2, cfx_limb_t l1, cfx_limb_t l0) {
    CFX_ASSERT(a->n == n);
    if (n >= 1) CFX_ASSERT(a->limb[0] == l0);
    if (n >= 2) CFX_ASSERT(a->limb[1] == l1);
    if (n >= 3) CFX_ASSERT(a->limb[2] == l2);
}

static void test_exp_edge_cases(void) {
    cfx_big_t n, p, out;
    cfx_big_init(&n);
    cfx_big_init(&p);
    cfx_big_init(&out);

    /* 0^0 := 1 */
    cfx_big_from_limb(&n, 0);
    cfx_big_from_limb(&p, 0);
    cfx_big_exp(&out, &n, &p);
    CFX_ASSERT(cfx_big_eq_u64(&out, 1));

    /* 0^k = 0 (k>0) */
    cfx_big_from_limb(&n, 0);
    cfx_big_from_limb(&p, 5);
    cfx_big_exp(&out, &n, &p);
    CFX_ASSERT(cfx_big_eq_u64(&out, 0));

    /* 1^p = 1 */
    cfx_big_from_limb(&n, 1);
    cfx_big_from_limb(&p, 1234567);
    cfx_big_exp(&out, &n, &p);
    CFX_ASSERT(cfx_big_eq_u64(&out, 1));

    /* n^0 = 1 */
    cfx_big_from_limb(&n, 42);
    cfx_big_from_limb(&p, 0);
    cfx_big_exp(&out, &n, &p);
    CFX_ASSERT(cfx_big_eq_u64(&out, 1));

    cfx_big_free(&out);
    cfx_big_free(&p);
    cfx_big_free(&n);
}

static void test_exp_small_values(void) {
    cfx_big_t n, p, out;
    cfx_big_init(&n);
    cfx_big_init(&p);
    cfx_big_init(&out);

    /* n^1 = n */
    cfx_big_from_limb(&n, 123456789);
    cfx_big_from_limb(&p, 1);
    cfx_big_exp(&out, &n, &p);
    CFX_ASSERT(out.n == 1 && out.limb[0] == 123456789);

    /* 2^10 = 1024 */
    cfx_big_from_limb(&n, 2);
    cfx_big_from_limb(&p, 10);
    cfx_big_exp(&out, &n, &p);
    CFX_ASSERT(cfx_big_eq_u64(&out, UINT64_C(1024)));

    /* 3^5 = 243 */
    cfx_big_from_limb(&n, 3);
    cfx_big_from_limb(&p, 5);
    cfx_big_exp(&out, &n, &p);
    CFX_ASSERT(cfx_big_eq_u64(&out, UINT64_C(243)));

    cfx_big_free(&out);
    cfx_big_free(&p);
    cfx_big_free(&n);
}

static void test_exp_powers_of_two_boundaries(void) {
    cfx_big_t n, p, out;
    cfx_big_init(&n);
    cfx_big_init(&p);
    cfx_big_init(&out);

#if CFX_LIMB_BITS == 64
    /* 2^64 = 1<<64 -> limbs [0]=0, [1]=1 */
    cfx_big_from_limb(&n, 2);
    cfx_big_from_limb(&p, 64);
    cfx_big_exp(&out, &n, &p);
    CFX_ASSERT(out.n == 2);
    expect_limb_pattern(&out, 2, 0, 1, 0);

    /* 2^128 = 1<<128 -> limbs [0]=0, [1]=0, [2]=1 */
    cfx_big_from_limb(&p, 128);
    cfx_big_exp(&out, &n, &p);
    CFX_ASSERT(out.n == 3);
    expect_limb_pattern(&out, 3, 1, 0, 0);

    /* 2^127 -> highest bit in limb[1], limb[0]=0 */
    cfx_big_from_limb(&p, 127);
    cfx_big_exp(&out, &n, &p);
    CFX_ASSERT(out.n == 2);
    CFX_ASSERT(out.limb[0] == 0);
    CFX_ASSERT(out.limb[1] == (UINT64_C(1) << 63));
#else
    /* 32-bit limbs: 2^32 = 1<<32 -> limbs [0]=0, [1]=1 */
    cfx_big_from_limb(&n, 2);
    cfx_big_from_limb(&p, 32);
    cfx_big_exp(&out, &n, &p);
    CFX_ASSERT(out.n == 2);
    expect_limb_pattern(&out, 2, 0, 1, 0);

    /* 2^64 = 1<<64 -> limbs [0]=0, [1]=0, [2]=1 */
    cfx_big_from_limb(&p, 64);
    cfx_big_exp(&out, &n, &p);
    CFX_ASSERT(out.n == 3);
    expect_limb_pattern(&out, 3, 1, 0, 0);

    /* 2^63 -> highest bit in limb[1], limb[0]=0 */
    cfx_big_from_limb(&p, 63);
    cfx_big_exp(&out, &n, &p);
    CFX_ASSERT(out.n == 2);
    CFX_ASSERT(out.limb[0] == 0);
    CFX_ASSERT(out.limb[1] == (UINT32_C(1) << 31));
#endif

    cfx_big_free(&out);
    cfx_big_free(&p);
    cfx_big_free(&n);
}

static void test_exp_aliasing(void) {
    cfx_big_t n, p;
    cfx_big_init(&n);
    cfx_big_init(&p);

    /* out aliases base (out == n) */
    cfx_big_from_limb(&n, 7);
    cfx_big_from_limb(&p, 6);      /* 7^6 = 117,649 */
    cfx_big_exp(&n, &n, &p);
    CFX_ASSERT(n.n == 1);
    /* 117,649 fits in both 32 and 64 bits */
    CFX_ASSERT(n.limb[0] == 117649);

#if CFX_LIMB_BITS == 64
    /* out aliases exponent (out == p) */
    cfx_big_from_limb(&n, 2);
    cfx_big_from_limb(&p, 80);     /* 2^80 */
    cfx_big_exp(&p, &n, &p);
    CFX_ASSERT(p.n == 2);
    /* 2^80 = 1<<80 => limb[0]=0, limb[1]=1<<16 */
    CFX_ASSERT(p.limb[0] == 0);
    CFX_ASSERT(p.limb[1] == (UINT64_C(1) << 16));
#else
    /* out aliases exponent (out == p) */
    cfx_big_from_limb(&n, 2);
    cfx_big_from_limb(&p, 48);     /* 2^48 */
    cfx_big_exp(&p, &n, &p);
    CFX_ASSERT(p.n == 2);
    /* 2^48 = 1<<48 => limb[0]=0, limb[1]=1<<16 */
    CFX_ASSERT(p.limb[0] == 0);
    CFX_ASSERT(p.limb[1] == (UINT32_C(1) << 16));
#endif

    cfx_big_free(&p);
    cfx_big_free(&n);
}

static void test_exp_compare_with_naive_mul(void) {
    cfx_big_t n, p, out1, out2;
    cfx_big_init(&n);
    cfx_big_init(&p);
    cfx_big_init(&out1);
    cfx_big_init(&out2);

    /* Random-ish small cases to avoid huge numbers; compare against repeated multiply. */
    cfx_limb_t bases[] = {2,3,5,10,17,1234567};
    cfx_limb_t exps[]  = {2,3,4,5,8,16};

    for (size_t i=0;i<sizeof(bases)/sizeof(bases[0]);++i) {
        for (size_t j=0;j<sizeof(exps)/sizeof(exps[0]);++j) {
            cfx_big_from_limb(&n, bases[i]);
            cfx_big_from_limb(&p, exps[j]);

            /* exp under test */
            cfx_big_exp(&out1, &n, &p);

            /* naive: res=1; repeat e times: res*=n */
            cfx_big_from_limb(&out2, 1);
            cfx_limb_t e = exps[j];
            while (e--) {
                cfx_big_mul_auto(&out2, &n);
            }

            CFX_ASSERT(cfx_big_cmp(&out1, &out2) == 0);
        }
    }

    cfx_big_free(&out2);
    cfx_big_free(&out1);
    cfx_big_free(&p);
    cfx_big_free(&n);
}

static void test_big_prime(void) {

    cfx_big_t b1, b2;
    cfx_big_init(&b1);
    cfx_big_init(&b2);

    unsigned cnt = 0;
    for (size_t i = 0; i < cfx_primes_len; i += 100) {
        ++cnt;
        cfx_limb_t val = (cfx_limb_t)cfx_primes[i];
        cfx_big_assign_sm(&b1, val);
        cfx_big_copy(&b2, &b1);
        cfx_big_add_sm(&b2, 1);
        char *s1 = cfx_big_to_str(&b1, NULL);
        char *s2 = cfx_big_to_str(&b2, NULL);
        if (!(i % 1000)) printf("[%u] %zu/%zu: %s - %s\n", cnt, i, cfx_primes_len, s1, s2);
        CFX_ASSERT(cfx_big_is_prime(&b1) == 1);
        if(b2.limb[0] != 3) CFX_ASSERT(cfx_big_is_prime(&b2) == 0);
        free(s1);
        free(s2);
    }

    cfx_big_from_dec(&b1, "1361129467683753853853498429727072845819");
    CFX_ASSERT(cfx_big_is_prime(&b1) == 1);

    /* Test 64-bit prime - largest prime < 2^64 */
    cfx_big_from_dec(&b1, "18446744073709551557");
    CFX_ASSERT(cfx_big_is_prime(&b1) == 1);

    /* 2^6756 - 3 */
    /*cfx_big_from_dec(&b1, "573654897752006159794481955898822158113920214270241209"
        "25635666276034335300856382095142846471417071668191140945094630698369875"
        "66631537436548447029478575758968992386708485640006992070480694219891934"
        "08983229865318123845679680448636474335504011750403715814625235853632108"
        "32510813384190821984771420350348401583245371154506981497337732247034659"
        "89791190497126242354416704708867761202850374319663630302534393183662102"
        "37774086306510627630105936191759170028672518096806607335355382988531084"
        "13585963855215682166412354665269420831009438931643015527054907321186480"
        "08467072224644169864407042736592703496932183838875536359164774443235772"
        "20698725918060185768692624254570270881633872872577471583780874734588286"
        "35061347591509383020246191029377755560779387137813325363778069629459980"
        "43914343832706878574687751720888568598104695522719230422380863383165507"
        "46306087371141763442604320381820184284035975250133853301399889868685078"
        "29967193069489899022218321815605311154856225966265518399582461191918084"
        "80740601908939367991006883297396812279708779391310859124668288600721777"
        "53226956265377619207398563399726073993896250088517589450437790246150190"
        "446927180095447103482425499667719313732291601456430261578372332019656677"
        "896453080055164822871840365336342438411290927543791099138036306374841591"
        "059558557877958856730302588409342153469885066048185407561394475290028330"
        "518130105586737950471943107949174683295982756994725213594680712496901475"
        "046756877495702743919758450936453095476252542722036986436973989150339214"
        "956017778894532128724317121313957388944591254440659626552888537996232368"
        "338089503544469686772932640408229670640831441116682886948777984032837767"
        "468030785628920398897109037757154108718950934160722362714341834056522445"
        "357041197542260971374241940260226120519489090741127736193132419113752015"
        "592809030821548196696397350276542321369146084922630068542903623401081786"
        "098279211918322047706495259370234867540414543464989826624620891566280407"
        "817130752574504656477565810657708387717230455941547154332894546757710272"
        "791794594236817128537076504733880398185421850279933");*/

    /* 2^5630 - 3 */
    /*cfx_big_from_hex(&b1,
        "3fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff"
        "fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff"
        "fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff"
        "fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff"
        "fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff"
        "fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff"
        "fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff"
        "fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff"
        "fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff"
        "fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff"
        "fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff"
        "fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff"
        "fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff"
        "fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff"
        "fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff"
        "fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff"
        "fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff"
        "fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff"
        "fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff"
        "fffffffffffffffffffffd");
        */
    CFX_ASSERT(cfx_big_is_prime(&b1) == 1);

    cfx_big_free(&b2);
    cfx_big_free(&b1);
}

static void test_binstring(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    cfx_big_from_bin(&b, "111");
    CFX_ASSERT(b.n == 1);
    CFX_ASSERT(b.limb[0] == 7);
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
    char* hex = cfx_big_to_hex(&b, NULL);
    const char* expected = "35555555555555555555555555555555555555555552aaaaaa6755555"
        "55555555555555555aaaaaaaaaaaabd4a040282527955415fffffffffe";
    int cmp = strcmp(hex, expected);
    printf("%s\n%s\n-- %d\n", expected, hex, cmp);
    CFX_ASSERT(cmp == 0);
    free(hex);
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
    char* bin = cfx_big_to_bin(&b, NULL);

    int cmp = strcmp(bin, expected);
    printf("%s\n%s\n-- %d\n", expected, bin, cmp);
    CFX_ASSERT(cmp == 0);
    free(bin);
}

static void test_bitsset(void) {
    cfx_big_t p;
    cfx_big_init(&p);
    cfx_big_from_u64(&p, 1);

    for (size_t k = 1; k <= 512; ++k) {
        cfx_big_shl_bits_eq(&p, 1);
        char* s = cfx_big_to_str(&p, NULL);
        printf("1 << %zu: %s\n", k, s);
        fflush(stdout);
        CFX_ASSERT(cfx_big_bit_is_set(&p, k));
    }
}


static cfx_big_t make_u64(uint64_t x) {
    cfx_big_t r;
    cfx_big_init(&r);
    cfx_big_from_u64(&r, x);
    return r;
}

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
    cfx_big_assign(&b, &a); // ensure distinct object with same value

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
    // Build b as a number requiring >64 bits.
    cfx_big_t a = make_u64(0x0102030405060708ULL);

    cfx_big_t b;
    cfx_big_init(&b);
    // Example: b = (1<<80) + 0x0f0e0d0c0b0a0908
    // - reserve limbs
    // - set limbs manually

    cfx_big_from_hex(&b, "f0e0d0c0b0a0908");

    cfx_big_t x;
    cfx_big_init(&x);
    cfx_big_assign_one(&x);
    cfx_big_shl_bits_eq(&x, 80);
    for (int k = 0; k < x.n; ++k)
        printf("%d %d\n", k, cfx_big_bit_is_set(&x, k));

    CFX_ASSERT(!cfx_big_bit_is_set(&x, 79));
    CFX_ASSERT(!cfx_big_bit_is_set(&x, 81));
    cfx_big_or(&b, &b, &x);

    cfx_big_t out;
    cfx_big_init(&out);

    cfx_big_xor(&out, &a, &b);

    cfx_big_t expected_low = make_u64( 0x0f0e0d0c0b0a0908 ^ 0x0102030405060708);

    CFX_ASSERT(cfx_big_endswith_u64(&out, 0x0e0c0e080e0c0e00ULL));
    PRINT_BIG("AAAAA", &out);
    CFX_ASSERT(cfx_big_bit_is_set(&out, 80));  // << fails

    cfx_big_free(&expected_low);
    cfx_big_free(&out);
    cfx_big_free(&b);
    cfx_big_free(&a);
}

static void test_xor_out_aliases_a(void) {
    cfx_big_t a = make_u64(0x1111111111111111ULL);
    cfx_big_t b = make_u64(0x0102030405060708ULL);

    cfx_big_xor(&a, &a, &b); // out == a

    cfx_big_t expected = make_u64(0x1013121514171619ULL);
    CFX_ASSERT(cfx_big_eq(&a, &expected));

    cfx_big_free(&expected);
    cfx_big_free(&b);
    cfx_big_free(&a);
}

static void test_xor_out_aliases_b(void) {
    cfx_big_t a = make_u64(0x1111111111111111ULL);
    cfx_big_t b = make_u64(0x0102030405060708ULL);

    cfx_big_xor(&b, &a, &b); // out == b

    cfx_big_t expected = make_u64(0x1013121514171619ULL);
    CFX_ASSERT(cfx_big_eq(&b, &expected));

    cfx_big_free(&expected);
    cfx_big_free(&b);
    cfx_big_free(&a);
}



static void assert_big_eq_hex(const cfx_big_t* x, const char* expected_hex) {
    size_t sz = 0;
    char* got = cfx_big_to_hex(x, &sz);
    CFX_ASSERT(got != NULL);

    /* Allow either "0" or "00...0" depending on your hex formatting.
       Easiest: compare after stripping leading zeros (except keep one). */
    const char* e = expected_hex;
    while (e[0] == '0' && e[1] != '\0') e++;

    char* g = got;
    while (g[0] == '0' && g[1] != '\0') g++;

    CFX_ASSERT(strcmp(g, e) == 0);

    free(got);
}

static void big_from_hex(cfx_big_t* out, const char* hex) {
    int rc = cfx_big_from_hex(out, hex);
    CFX_ASSERT(rc == 0);
}

static void big_from_u64(cfx_big_t* out, uint64_t v) {
    int rc = cfx_big_from_u64(out, v);
    CFX_ASSERT(rc == 0);
}

static void test_and_basic_u64(void){
    cfx_big_t a, b, out;
    cfx_big_init(&a);
    cfx_big_init(&b);
    cfx_big_init(&out);

    big_from_u64(&a, 0xF0F0ULL);
    big_from_u64(&b, 0x0FF0ULL);

    cfx_big_and(&out, &a, &b);
    /* 0xF0F0 & 0x0FF0 = 0x00F0 */
    assert_big_eq_hex(&out, "f0");

    cfx_big_free(&out);
    cfx_big_free(&b);
    cfx_big_free(&a);
}

static void test_and_commutative(void)
{
    cfx_big_t a, b, out1, out2;
    cfx_big_init(&a);
    cfx_big_init(&b);
    cfx_big_init(&out1);
    cfx_big_init(&out2);

    big_from_u64(&a, 0x123456789abcdef0ULL);
    big_from_u64(&b, 0x0f0f0f0f0f0f0f0fULL);

    cfx_big_and(&out1, &a, &b);
    cfx_big_and(&out2, &b, &a);

    CFX_ASSERT(cfx_big_eq(&out1, &out2));

    cfx_big_free(&out2);
    cfx_big_free(&out1);
    cfx_big_free(&b);
    cfx_big_free(&a);
}

static void test_and_different_lengths_expected_zero(void) {
    /* a = 0xffffffffffffffff
       b = 1<<128  (hex: 1 followed by 32 zeros)
       a & b should be 0 */
    cfx_big_t a, b, out;
    cfx_big_init(&a);
    cfx_big_init(&b);
    cfx_big_init(&out);

    big_from_u64(&a, 0xFFFFFFFFFFFFFFFFULL);
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

    /* small:  0x00ff00ff */
    big_from_hex(&small, "00ff00ff");

    /* large:  0x123400ff00aa (note: longer) */
    big_from_hex(&large, "123400ff00aa");

    /* expected = small & low(same-limb-count-of-small) of large
       0x00ff00ff
     & 0x00ff00aa
       ----------
       0x00ff00aa
    */
    big_from_hex(&expected, "00ff00aa");

    /* out aliases the larger operand */
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

    /* out aliases the smaller operand */
    cfx_big_and(&small, &small, &large);

    CFX_ASSERT(cfx_big_eq(&small, &expected));

    cfx_big_free(&expected);
    cfx_big_free(&large);
    cfx_big_free(&small);
}

static void test_or_basic_u64(void) {
    cfx_big_t a, b, out;
    cfx_big_init(&a);
    cfx_big_init(&b);
    cfx_big_init(&out);

    big_from_u64(&a, 0xF0F0ULL);
    big_from_u64(&b, 0x0FF0ULL);

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

    /* expected: 0x1234 | (low part ORed)
       low: 0x00ff00aa | 0x00ff00ff = 0x00ff00ff
       => 0x123400ff00ff
    */
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

static void test_or_alias_out_is_larger_operand(void)
{
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

#ifdef _MSC_VER
#define strcasecmp _stricmp
#endif

#define ASSERT_BIG_HEX(b, hex) \
    do { \
        char* s = cfx_big_to_hex(b, NULL); \
        CFX_ASSERT(strcasecmp(s, hex) == 0); \
        free(s); \
    } while (0)

void test_rot_zero_and_one(void) {
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

void test_rot_8bit(void) {
    cfx_big_t x, y;
    cfx_big_init(&x);
    cfx_big_init(&y);

    /* x = 0b10000001 */
    cfx_big_from_bin(&x, "10000001");
    ASSERT_BIG_HEX(&x, "81");
    // cfx_big_mask_bits(&x, 8);
    char* s = cfx_big_to_bin(&x, NULL);
    printf("first: %s\n", s);
    free(s);
    s = cfx_big_to_bin(&x, NULL);
    printf("before %s\n", s);
    cfx_big_rotl(&y, &x, 1);
    free(s);
    s = cfx_big_to_bin(&y, NULL);
    printf("after %s\n", s);
    ASSERT_BIG_HEX(&y, "3"); /* 00000011 */


    cfx_big_rotr(&y, &x, 1);
    s = cfx_big_to_bin(&y, NULL);
    printf("finally %s\n", s);
    free(s);
    ASSERT_BIG_HEX(&y, "c0"); /* 11000000 */

    cfx_big_free(&x);
    cfx_big_free(&y);
}

void test_rot_by_width(void) {
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

void test_rot_duality(void) {
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

void test_rot_composition(void) {
    cfx_big_t x, a, b;
    cfx_big_init(&x);
    cfx_big_init(&a);
    cfx_big_init(&b);

    cfx_big_from_hex(&x, "CAFEBABEDEADC0DE");
    unsigned w = cfx_big_bitlen(&x);

    /* Use explicit-width rotations to preserve bitlen across compositions */
    cfx_big_rotl_w(&a, &x, 7, w);
    cfx_big_rotl_w(&a, &a, 19, w);

    cfx_big_rotl_w(&b, &x, (7 + 19) % w, w);

    CFX_ASSERT(cfx_big_cmp(&a, &b) == 0);

    cfx_big_free(&x);
    cfx_big_free(&a);
    cfx_big_free(&b);
}

void test_rot_aliasing(void) {
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

void test_rot_cross_limb(void) {
    cfx_big_t x, y;
    cfx_big_init(&x);
    cfx_big_init(&y);

    /* bit exactly at limb boundary */
    cfx_big_from_limb(&x, 1);

    char* s = cfx_big_to_bin(&x, NULL);
    printf("first: %s\n", s);
    free(s);

    cfx_big_shl_bits(&x, &x, CFX_LIMB_BITS);

    s = cfx_big_to_bin(&x, NULL);
    printf("then: %s\n", s);
    free(s);

    cfx_big_mask_bits(&x, CFX_LIMB_BITS + 1);

    s = cfx_big_to_bin(&x, NULL);
    printf("masked: %s\n", s);
    free(s);

    cfx_big_rotl(&y, &x, 1);

    s = cfx_big_to_bin(&y, NULL);
    printf("rotl: %s\n", s);
    free(s);

    /* should wrap the top bit back to LSB */
    ASSERT_BIG_HEX(&y, "1");

    cfx_big_free(&x);
    cfx_big_free(&y);
}

void test_rotl_w_basic(void) {
    cfx_big_t x, y;
    cfx_big_init(&x);
    cfx_big_init(&y);

    cfx_big_from_hex(&x, "ab");
    cfx_big_rotl_w(&y, &x, 4, 8);
    ASSERT_BIG_HEX(&y, "ba");

    cfx_big_free(&x);
    cfx_big_free(&y);
}

void test_rotr_w_basic(void) {
    cfx_big_t x, y;
    cfx_big_init(&x);
    cfx_big_init(&y);

    cfx_big_from_hex(&x, "ab");
    cfx_big_rotr_w(&y, &x, 4, 8);
    ASSERT_BIG_HEX(&y, "ba");

    cfx_big_free(&x);
    cfx_big_free(&y);
}

void test_rotl_w_larger_width(void) {
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

void test_rot_w_zero_rotation(void) {
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

void test_rot_w_full_rotation(void) {
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

void test_rot_w_aliasing(void) {
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

void test_rot_w_duality(void) {
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

void test_rot_w_power_of_two_width(void) {
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

/*  big endian test helper */
static void hex_to_bytes(uint8_t* out, size_t out_len, const char *hex) {
    /* hex length must be 2*out_len, no "0x", no spaces */
    for (size_t i = 0; i < out_len; i++) {
        unsigned v = 0;
        char c1 = hex[2*i], c2 = hex[2*i + 1];
        CFX_ASSERT(c1 && c2);

        if      (c1 >= '0' && c1 <= '9') v = (unsigned)(c1 - '0') << 4;
        else if (c1 >= 'a' && c1 <= 'f') v = (unsigned)(c1 - 'a' + 10) << 4;
        else if (c1 >= 'A' && c1 <= 'F') v = (unsigned)(c1 - 'A' + 10) << 4;
        else CFX_ASSERT(0);

        if      (c2 >= '0' && c2 <= '9') v |= (unsigned)(c2 - '0');
        else if (c2 >= 'a' && c2 <= 'f') v |= (unsigned)(c2 - 'a' + 10);
        else if (c2 >= 'A' && c2 <= 'F') v |= (unsigned)(c2 - 'A' + 10);
        else CFX_ASSERT(0);

        out[i] = (uint8_t)v;
    }
}

/* -------------------------------------------------------------------------- */

static void test_big_to_bytes_be_zero(void) {
    cfx_big_t b;
    cfx_big_init(&b);

    cfx_big_from_limb(&b, 0);

    /* size query */
    size_t need = 123;
    int rc = cfx_big_to_bytes_be(NULL, &need, &b);
    CFX_ASSERT(rc == 0);
    CFX_ASSERT(need == 0);

    /* write into buffer */
    uint8_t out[8];
    memset(out, 0xAA, sizeof(out));
    size_t out_len = sizeof(out);
    rc = cfx_big_to_bytes_be(out, &out_len, &b);
    CFX_ASSERT(rc == 0);
    CFX_ASSERT(out_len == 0);

    cfx_big_free(&b);
}


static void test_big_to_bytes_be(void) {

    struct {
        const char*hex_in;
        const char *hex_be;
    } cases[] = {
        { "01", "01" },
        { "ff", "FF" },
        { "0100", "0100" },
        { "123456", "123456" },
        { "1234567890abcdef", "1234567890ABCDEF" },
        { "0102030405060708090a0b0c0d0e0f10", "0102030405060708090A0B0C0D0E0F10" },
        { "8000000000000000", "8000000000000000" },         /* top bit set */
        { "ffffffffffffffffffffffffffffffff", "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF" },
        { "010000000000000000", "010000000000000000" },
        { "010203040506070809090909090909", "010203040506070809090909090909"}
    };

    cfx_big_t b;
    cfx_big_init(&b);

    for (size_t i = 0; i < sizeof(cases)/sizeof(cases[0]); i++) {
        int ok = cfx_big_from_hex(&b, cases[i].hex_in);
        CFX_ASSERT(ok == 0);

        size_t need = 0;
        int rc = cfx_big_to_bytes_be(NULL, &need, &b);
        CFX_ASSERT(rc == 0);

        if (strcmp(cases[i].hex_in, "0") == 0) {
            CFX_ASSERT(need == 0);
            continue;
        }

        size_t exp_len = strlen(cases[i].hex_be)/2;
        CFX_ASSERT(need == exp_len);

        uint8_t exp[64];
        hex_to_bytes(exp, exp_len, cases[i].hex_be);

        uint8_t out[64];
        size_t out_len = sizeof(out);
        rc = cfx_big_to_bytes_be(out, &out_len, &b);
        CFX_ASSERT(rc == 0);
        CFX_ASSERT(out_len == exp_len);
        CFX_ASSERT(memcmp(out, exp, exp_len) == 0);
    }

    cfx_big_free(&b);

}

static void test_big_to_bytes_be_buffer_too_small(void) {
    cfx_big_t b;
    cfx_big_init(&b);

    cfx_big_from_limb(&b, 0x010203u);

    size_t need = 0;
    int rc = cfx_big_to_bytes_be(NULL, &need, &b);
    CFX_ASSERT(rc == 0);
    CFX_ASSERT(need == 3);

    uint8_t out[3];
    size_t out_len = 2; /* too small on purpose */
    rc = cfx_big_to_bytes_be(out, &out_len, &b);
    CFX_ASSERT(rc == -1);

    cfx_big_free(&b);
}

/* ------------------------------------------------------------------ */
/* Tests for cfx_big_gcd                                              */
/* ------------------------------------------------------------------ */

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

    /* gcd of two large numbers with known gcd */
    /* a = 123456789 * 1000, b = 123456789 * 777 => gcd = 123456789 */
    cfx_big_from_str(&a, "123456789000");
    cfx_big_from_str(&b, "95925925053");  /* 123456789 * 777 */
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

/* ------------------------------------------------------------------ */
/* Tests for cfx_big_pollard_rho                                      */
/* ------------------------------------------------------------------ */

static void test_big_pollard_rho_small_composite(void) {
    cfx_big_t n, factor;
    cfx_big_init(&n);
    cfx_big_init(&factor);

    /* n = 15 = 3 * 5 */
    cfx_big_from_limb(&n, 15);
    cfx_big_pollard_rho(&factor, &n);
    /* Should find either 3 or 5 */
    CFX_ASSERT(factor.n == 1);
    CFX_ASSERT(factor.limb[0] == 3 || factor.limb[0] == 5);

    /* n = 21 = 3 * 7 */
    cfx_big_from_limb(&n, 21);
    cfx_big_pollard_rho(&factor, &n);
    CFX_ASSERT(factor.n == 1);
    CFX_ASSERT(factor.limb[0] == 3 || factor.limb[0] == 7);

    cfx_big_free(&n);
    cfx_big_free(&factor);
}

static void test_big_pollard_rho_semiprime(void) {
    cfx_big_t n, factor, quotient, check;
    cfx_big_init(&n);
    cfx_big_init(&factor);
    cfx_big_init(&quotient);
    cfx_big_init(&check);

    /* n = 1000003 * 1000033 = 1000036000099 (semiprime) */
    cfx_big_from_str(&n, "1000036000099");
    cfx_big_pollard_rho(&factor, &n);

    /* factor should be non-trivial (not 1 and not n) */
    CFX_ASSERT(!cfx_big_is_one(&factor));
    CFX_ASSERT(cfx_big_cmp(&factor, &n) != 0);

    /* Verify: n should be divisible by factor */
    cfx_big_copy(&quotient, &n);
    cfx_big_div_eq(&quotient, &factor, NULL);
    cfx_big_mul(&check, &quotient, &factor);
    CFX_ASSERT(cfx_big_cmp(&check, &n) == 0);

    cfx_big_free(&n);
    cfx_big_free(&factor);
    cfx_big_free(&quotient);
    cfx_big_free(&check);
}

static void test_big_pollard_rho_prime(void) {
    cfx_big_t n, factor;
    cfx_big_init(&n);
    cfx_big_init(&factor);

    /* n = 104729 (prime) - Pollard-Rho should return n itself */
    cfx_big_from_limb(&n, 104729);
    cfx_big_pollard_rho(&factor, &n);
    /* For a prime, it should return the prime itself (can't factor) */
    CFX_ASSERT(cfx_big_cmp(&factor, &n) == 0);

    cfx_big_free(&n);
    cfx_big_free(&factor);
}

static void test_big_pollard_rho_even(void) {
    cfx_big_t n, factor;
    cfx_big_init(&n);
    cfx_big_init(&factor);

    /* Even number should return 2 */
    cfx_big_from_limb(&n, 1234567890);
    cfx_big_pollard_rho(&factor, &n);
    CFX_ASSERT(factor.n == 1 && factor.limb[0] == 2);

    cfx_big_free(&n);
    cfx_big_free(&factor);
}

static void test_big_pollard_rho_large_semiprime(void) {
    cfx_big_t n, factor, quotient, check;
    cfx_big_init(&n);
    cfx_big_init(&factor);
    cfx_big_init(&quotient);
    cfx_big_init(&check);

    /* n = 10007 * 10009 = 100160063 */
    cfx_big_from_str(&n, "100160063");
    cfx_big_pollard_rho(&factor, &n);

    /* factor should be non-trivial */
    CFX_ASSERT(!cfx_big_is_one(&factor));
    CFX_ASSERT(cfx_big_cmp(&factor, &n) != 0);

    /* Verify divisibility */
    cfx_big_copy(&quotient, &n);
    cfx_big_div_eq(&quotient, &factor, NULL);
    cfx_big_mul(&check, &quotient, &factor);
    CFX_ASSERT(cfx_big_cmp(&check, &n) == 0);

    cfx_big_free(&n);
    cfx_big_free(&factor);
    cfx_big_free(&quotient);
    cfx_big_free(&check);
}


/* ------------------------------------------------------------------ */
int main(void) {
    CFX_TEST(test_cfx_big_assign);
    CFX_TEST(test_copy_swap);
    CFX_TEST(test_cswap);
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
    CFX_TEST(test_from_str_matches1);
    CFX_TEST(test_from_str_matches2);
    CFX_TEST(test_from_str_zero_forms);
    CFX_TEST(test_from_str_whitespace_ok);
    CFX_TEST(test_from_str_hex_prefix_case);
    CFX_TEST(test_from_str_bin_prefix_case);
    CFX_TEST(test_from_str_rejects_trailing_junk);
    CFX_TEST(test_from_str_rejects_empty);
    CFX_TEST(test_from_str_rejects_prefix_only);
    CFX_TEST(test_from_str_no_legacy_octal);
    CFX_TEST(test_from_str_limb_boundary_all_ones);
    CFX_TEST(test_from_str_limb_boundary_power_of_two);
    CFX_TEST(test_scan_num_hex_basic);
    CFX_TEST(test_scan_num_bin_basic);
    CFX_TEST(test_scan_num_oct_basic);
    CFX_TEST(test_scan_num_dec_basic);
    CFX_TEST(test_scan_num_stops_at_invalid_digit);
    CFX_TEST(test_scan_num_prefix_only_fails);
    CFX_TEST(test_scan_num_non_number_fails);
    CFX_TEST(test_scan_num_prefix_case_insensitive);
    CFX_TEST(test_scan_num_respects_in_len);
    CFX_TEST(test_scan_num_b64_basic);
    CFX_TEST(test_scan_num_b64_whitespace_inside);
    CFX_TEST(test_scan_num_b64_multibyte_value);
    CFX_TEST(test_scan_num_b64_stops_before_junk);
    CFX_TEST(test_from_str_b64_allows_trailing_spaces_only);
    CFX_TEST(test_from_str_b64_rejects_invalid);
    CFX_TEST(test_add_sm);
    CFX_TEST(test_sub);
    CFX_TEST(test_sub_sm);
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
    CFX_TEST(test_big_prime);
    CFX_TEST(test_binstring);
    CFX_TEST(test_binstring2);
    CFX_TEST(test_binstring3);
    CFX_TEST(test_bitsset);
    CFX_TEST(test_xor_with_zero_is_identity);
    CFX_TEST(test_xor_self_is_zero);
    CFX_TEST(test_xor_commutative);
    CFX_TEST(test_xor_different_lengths_keep_high_limb_from_larger);
    CFX_TEST(test_xor_out_aliases_a);
    CFX_TEST(test_xor_out_aliases_b);
    CFX_TEST(test_and_basic_u64);
    CFX_TEST(test_and_commutative);
    CFX_TEST(test_and_different_lengths_expected_zero);
    CFX_TEST(test_and_alias_out_is_larger_operand);
    CFX_TEST(test_and_alias_out_is_smaller_operand);
    CFX_TEST(test_or_basic_u64);
    CFX_TEST(test_or_commutative);
    CFX_TEST(test_or_different_lengths_keeps_high_limbs);
    CFX_TEST(test_or_alias_out_is_smaller_operand);
    CFX_TEST(test_or_alias_out_is_larger_operand);
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
    CFX_TEST(test_big_to_bytes_be_zero);
    CFX_TEST(test_big_to_bytes_be);
    CFX_TEST(test_big_to_bytes_be_buffer_too_small);
    CFX_TEST(test_big_gcd_basic);
    CFX_TEST(test_big_gcd_large);
    CFX_TEST(test_big_pollard_rho_small_composite);
    CFX_TEST(test_big_pollard_rho_semiprime);
    CFX_TEST(test_big_pollard_rho_prime);
    CFX_TEST(test_big_pollard_rho_even);
    CFX_TEST(test_big_pollard_rho_large_semiprime);
    puts("OK");
    return 0;
}
