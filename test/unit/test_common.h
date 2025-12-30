/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */
/* test_common.h - Shared test helpers for cfx unit tests */

#ifndef CFX_TEST_COMMON_H
#define CFX_TEST_COMMON_H

#include "cfx/big.h"
#include "cfx/macros.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdio.h>

#ifndef _MSC_VER
#include <strings.h>
#else
#define strcasecmp _stricmp
#endif

/* ---- Assertion helpers -------------------------------------------------- */

#define ASSERT_BIG_HEX(b, hex) \
    do { \
        char* s = cfx_big_to_hex(b, NULL); \
        CFX_ASSERT(strcasecmp(s, hex) == 0); \
        free(s); \
    } while (0)

/* ---- Big integer creation helpers --------------------------------------- */

static inline cfx_big_t make_u64(uint64_t x) {
    cfx_big_t r;
    cfx_big_init(&r);
    cfx_big_from_u64(&r, x);
    return r;
}

static inline void big_from_hex(cfx_big_t* out, const char* hex) {
    int rc = cfx_big_from_hex(out, hex);
    CFX_ASSERT(rc == 0);
}

static inline void big_from_u64(cfx_big_t* out, uint64_t v) {
    int rc = cfx_big_from_u64(out, v);
    CFX_ASSERT(rc == 0);
}

static inline void big_init_from_limbs_base_1e9(cfx_big_t* b, const cfx_limb_t *limbs, size_t n) {
    cfx_big_init(b);
    if (n == 0) return;
    for (size_t i = n; i--;) {
        cfx_big_mul_sm(b, UINT64_C(1000000000));
        cfx_big_add_sm(b, limbs[i]);
    }
}

/* ---- Comparison/assertion helpers --------------------------------------- */

static inline void assert_hex_eq(const char* tag, const cfx_big_t* x, const char* hex_exp) {
    char* got = cfx_big_to_hex(x, NULL);
    int ok = (strcmp(got, hex_exp) == 0);
    if (!ok) {
        fprintf(stderr, "[%s] expected 0x%s, got 0x%s\n", tag, hex_exp, got);
    }
    CFX_ASSERT_PRINT(ok);
    free(got);
}

static inline void assert_big_eq_hex(const cfx_big_t* x, const char* expected_hex) {
    size_t sz = 0;
    char* got = cfx_big_to_hex(x, &sz);
    CFX_ASSERT(got != NULL);

    /* Allow either "0" or "00...0" depending on hex formatting.
       Compare after stripping leading zeros (except keep one). */
    const char* e = expected_hex;
    while (e[0] == '0' && e[1] != '\0') e++;

    char* g = got;
    while (g[0] == '0' && g[1] != '\0') g++;

    CFX_ASSERT(strcmp(g, e) == 0);
    free(got);
}

static inline void expect_dec_eq(const cfx_big_t* x, const char* dec) {
    cfx_big_t tmp;
    cfx_big_init(&tmp);
    cfx_big_from_dec(&tmp, dec);
    CFX_ASSERT(cfx_big_eq(x, &tmp));
    cfx_big_free(&tmp);
}

static inline void big_expect_limbs(const char* s, const cfx_big_t* b, const cfx_limb_t* limbs, size_t n) {
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
}

static inline void expect_limb_pattern(const cfx_big_t* a, size_t n,
                                cfx_limb_t l2, cfx_limb_t l1, cfx_limb_t l0) {
    CFX_ASSERT(a->n == n);
    if (n >= 1) CFX_ASSERT(a->limb[0] == l0);
    if (n >= 2) CFX_ASSERT(a->limb[1] == l1);
    if (n >= 3) CFX_ASSERT(a->limb[2] == l2);
}

/* Division verification: n == q*d + r and r < d */
static inline void assert_n_eq_qd_plus_r(const cfx_big_t* n, const cfx_big_t* q,
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

/* ---- String conversion test helper -------------------------------------- */

static inline void check_str_conversion(const char *label, const cfx_limb_t *limbs, size_t n, const char *expect) {
    cfx_big_t b;
    big_init_from_limbs_base_1e9(&b, limbs, n);
    size_t len = 0;
    char *s = cfx_big_to_str(&b, &len);
    CFX_ASSERT(strcmp(s, expect) == 0);
    CFX_ASSERT(len == strlen(expect));
    CFX_ASSERT(s[len] == '\0');
    CFX_PRINT_DBG("[ok] %s -> %s\n", label, s);
    free(s);
    cfx_big_free(&b);
}

/* ---- Shift test helpers ------------------------------------------------- */

static inline void check_shl_case(const char* msg, const char* hex_in, unsigned s, const char* hex_exp) {
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

static inline void check_shr_case(const char* msg, const char* hex_in, unsigned s, const char* hex_exp) {
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

/* ---- Byte conversion helpers -------------------------------------------- */

static inline void hex_to_bytes(uint8_t* out, size_t out_len, const char *hex) {
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

/* ---- Debug printing helpers --------------------------------------------- */

static inline void print_limbs(const char* name, const cfx_big_t* b) {
    printf("%s: n=%zu, limbs = {", name, b->n);
    for (size_t i = 0; i < b->n; i++) {
        printf("0x" CFX_PRIxLIMB "%s", b->limb[i], i < b->n - 1 ? ", " : "");
    }
    printf("}\n");
}

/* Decrement by 1 (helper for some tests) */
static inline void big_dec1(cfx_big_t* x) {
    cfx_limb_t borrow = 1;
    for (size_t i = 0; i < x->n && borrow; ++i) {
        cfx_limb_t old = x->limb[i];
        x->limb[i] = old - borrow;
        borrow = (x->limb[i] > old);
    }
}

#endif /* CFX_TEST_COMMON_H */
