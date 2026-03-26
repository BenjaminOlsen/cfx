#include "cfx/rng_lcg.h"
#include "cfx/macros.h"

#include <string.h>
#include <stdio.h>

/* LCG (Numerical Recipes) */

static void test_lcg_deterministic(void) {
    uint32_t s1 = 42, s2 = 42;
    for (int i = 0; i < 100; i++) {
        CFX_ASSERT(cfx_lcg(&s1) == cfx_lcg(&s2));
    }
}

static void test_lcg_seed_gen32(void) {
    cfx_lcg_seed(123);
    uint32_t a = cfx_lcg_gen32();
    uint32_t b = cfx_lcg_gen32();
    CFX_ASSERT(a != b);

    cfx_lcg_seed(123);
    CFX_ASSERT(cfx_lcg_gen32() == a);
    CFX_ASSERT(cfx_lcg_gen32() == b);
}

static void test_lcg_known(void) {
    /* s = s * 1664525 + 1013904223 (mod 2^32) */
    uint32_t s = 0;
    CFX_ASSERT(cfx_lcg(&s) == 1013904223u);
    CFX_ASSERT(s == 1013904223u);
}

static void test_lcg_bytes(void) {
    uint8_t buf[32];
    cfx_lcg_bytes(42, buf, sizeof(buf));

    /* deterministic */
    uint8_t buf2[32];
    cfx_lcg_bytes(42, buf2, sizeof(buf2));
    CFX_ASSERT(memcmp(buf, buf2, sizeof(buf)) == 0);

    /* matches manual: each byte = (x >> 24) where x advances one LCG step */
    uint32_t x = 42;
    for (int i = 0; i < 32; i++) {
        x = 1664525u * x + 1013904223u;
        CFX_ASSERT(buf[i] == (uint8_t)(x >> 24));
    }
}

static void test_lcg_bytes_zero_len(void) {
    uint8_t buf[4] = {0xAA, 0xAA, 0xAA, 0xAA};
    cfx_lcg_bytes(1, buf, 0);
    for (int i = 0; i < 4; i++)
        CFX_ASSERT(buf[i] == 0xAA);
}

static void test_lcg_different_seeds(void) {
    uint32_t s1 = 0, s2 = 1;
    CFX_ASSERT(cfx_lcg(&s1) != cfx_lcg(&s2));
}

/* drand48 */

static void test_drand48_deterministic(void) {
    uint64_t s1 = 42, s2 = 42;
    for (int i = 0; i < 100; i++) {
        CFX_ASSERT(cfx_drand48(&s1) == cfx_drand48(&s2));
    }
}

static void test_drand48_seed_gen32(void) {
    cfx_drand48_seed(99);
    uint32_t a = cfx_drand48_gen32();
    uint32_t b = cfx_drand48_gen32();

    cfx_drand48_seed(99);
    CFX_ASSERT(cfx_drand48_gen32() == a);
    CFX_ASSERT(cfx_drand48_gen32() == b);
}

static void test_drand48_known(void) {
    /* s = (s * 25214903917 + 11) & 0xFFFFFFFFFFFF; return s >> 16 */
    uint64_t s = 0;
    uint32_t r = cfx_drand48(&s);
    CFX_ASSERT(r == 0);      /* 11 >> 16 = 0 */
    CFX_ASSERT(s == 11);
}

static void test_drand48_state_48bit(void) {
    uint64_t s = UINT64_C(0xFFFFFFFFFFFF);
    for (int i = 0; i < 1000; i++) {
        cfx_drand48(&s);
        CFX_ASSERT(s <= UINT64_C(0xFFFFFFFFFFFF));
    }
}

static void test_drand48_bytes(void) {
    cfx_drand48_seed(77);
    uint8_t buf1[64];
    cfx_drand48_bytes(buf1, sizeof(buf1));

    cfx_drand48_seed(77);
    uint8_t buf2[64];
    cfx_drand48_bytes(buf2, sizeof(buf2));

    CFX_ASSERT(memcmp(buf1, buf2, sizeof(buf1)) == 0);
}

/* MINSTD (Park-Miller) */

static void test_minstd_deterministic(void) {
    uint64_t s1 = 42, s2 = 42;
    for (int i = 0; i < 100; i++) {
        CFX_ASSERT(cfx_minstd(&s1) == cfx_minstd(&s2));
    }
}

static void test_minstd_seed_gen32(void) {
    cfx_minstd_seed(7);
    uint32_t a = cfx_minstd_gen32();
    uint32_t b = cfx_minstd_gen32();

    cfx_minstd_seed(7);
    CFX_ASSERT(cfx_minstd_gen32() == a);
    CFX_ASSERT(cfx_minstd_gen32() == b);
}

static void test_minstd_known(void) {
    /* Park-Miller: s = s * 16807 mod 2147483647
       Known sequence from seed=1: 16807, 282475249, 1622650073 */
    uint64_t s = 1;
    CFX_ASSERT(cfx_minstd(&s) == 16807u);
    CFX_ASSERT(cfx_minstd(&s) == 282475249u);
    CFX_ASSERT(cfx_minstd(&s) == 1622650073u);
}

static void test_minstd_seed_zero(void) {
    /* seed=0 should map to state=1 to avoid absorbing state */
    cfx_minstd_seed(0);
    uint32_t a = cfx_minstd_gen32();
    cfx_minstd_seed(1);
    uint32_t b = cfx_minstd_gen32();
    CFX_ASSERT(a == b);
}

static void test_minstd_bytes(void) {
    cfx_minstd_seed(33);
    uint8_t buf1[64];
    cfx_minstd_bytes(buf1, sizeof(buf1));

    cfx_minstd_seed(33);
    uint8_t buf2[64];
    cfx_minstd_bytes(buf2, sizeof(buf2));

    CFX_ASSERT(memcmp(buf1, buf2, sizeof(buf1)) == 0);
}

static void test_minstd_range(void) {
    /* output must be in [1, 2^31 - 2] for MINSTD */
    uint64_t s = 12345;
    for (int i = 0; i < 10000; i++) {
        uint32_t v = cfx_minstd(&s);
        CFX_ASSERT(v >= 1 && v <= 2147483646u);
    }
}

/* main */

int main(void) {
    CFX_TEST(test_lcg_deterministic);
    CFX_TEST(test_lcg_seed_gen32);
    CFX_TEST(test_lcg_known);
    CFX_TEST(test_lcg_bytes);
    CFX_TEST(test_lcg_bytes_zero_len);
    CFX_TEST(test_lcg_different_seeds);

    CFX_TEST(test_drand48_deterministic);
    CFX_TEST(test_drand48_seed_gen32);
    CFX_TEST(test_drand48_known);
    CFX_TEST(test_drand48_state_48bit);
    CFX_TEST(test_drand48_bytes);

    CFX_TEST(test_minstd_deterministic);
    CFX_TEST(test_minstd_seed_gen32);
    CFX_TEST(test_minstd_known);
    CFX_TEST(test_minstd_seed_zero);
    CFX_TEST(test_minstd_bytes);
    CFX_TEST(test_minstd_range);

    puts("OK");
    return 0;
}
