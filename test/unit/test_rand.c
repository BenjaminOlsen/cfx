/* test_rand.c */

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "cfx/macros.h"
#include "cfx/rand.h"

#define N 16

static void test_rng_reproducible(const cfx_rand_desc_t *rng) {
    const uint32_t SEED = 0x12345678u;
    uint32_t seq1[N], seq2[N];

    rng->seed(SEED);
    for (size_t i = 0; i < N; ++i) seq1[i] = rng->rng32();

    rng->seed(SEED);
    for (size_t i = 0; i < N; ++i) seq2[i] = rng->rng32();

    CFX_ASSERT(memcmp(seq1, seq2, sizeof seq1) == 0);
}

/* catches "seed is ignored" bugs */
static void test_rng_different_seeds_differ(const cfx_rand_desc_t* rng) {
    uint32_t seq1[N], seq2[N];

    rng->seed(1u);
    for (size_t i = 0; i < N; ++i) seq1[i] = rng->rng32();

    rng->seed(2u);
    for (size_t i = 0; i < N; ++i) seq2[i] = rng->rng32();

    size_t equal = 0;
    for (size_t i = 0; i < N; ++i) {
        if (seq1[i] == seq2[i]) ++equal;
    }
    CFX_ASSERT(equal < N);  /* allow a few collisions, but not all */
}

/* checks output isn't degenerate (all zeros, all ones, or constant) */
static void test_rng_non_degeneracy(const cfx_rand_desc_t* rng) {
    rng->seed(42u);

    uint32_t or_acc = 0;
    uint32_t and_acc = 0xFFFFFFFFu;
    uint32_t prev = rng->rng32();
    size_t same_count = 0;

    for (size_t i = 0; i < 100; ++i) {
        uint32_t v = rng->rng32();
        or_acc |= v;
        and_acc &= v;
        if (v == prev) ++same_count;
        prev = v;
    }

    CFX_ASSERT(or_acc != 0);            /* not all zeros */
    CFX_ASSERT(and_acc != 0xFFFFFFFFu); /* not all ones */
    CFX_ASSERT(same_count < 50);        /* not constant output */
}

/* test _bytes() api with various lengths */
static void test_rng_bytes_api(const cfx_rand_desc_t* rng) {
    uint8_t buf[128];
    size_t lengths[] = {0, 1, 3, 4, 7, 8, 15, 16, 31, 32, 100};
    size_t num_lengths = sizeof(lengths) / sizeof(lengths[0]);

    for (size_t i = 0; i < num_lengths; ++i) {
        size_t len = lengths[i];
        memset(buf, 0xAA, sizeof buf);
        rng->seed(123u);
        rng->bytes(buf, len);

        /* bytes beyond len should be untouched */
        if (len < sizeof buf) {
            CFX_ASSERT(buf[len] == 0xAA);
        }

        /* for non-zero lengths, check we got some output */
        if (len >= 4) {
            uint32_t sum = 0;
            for (size_t j = 0; j < len; ++j) sum += buf[j];
            CFX_ASSERT(sum != 0);  /* very unlikely to be all zeros */
        }
    }
}

static void test_rng_zero_seed(const cfx_rand_desc_t* rng) {
    rng->seed(0u);
    uint32_t v1 = rng->rng32();
    uint32_t v2 = rng->rng32();
    /* just checking it doesn't crash and produces something */
    (void)v1;
    CFX_ASSERT(v1 != v2 || v1 != rng->rng32());  /* not stuck */
}

static void run_all_table_rng_tests(void) {
    for (size_t i = 0; i < g_rand_gen_cnt; ++i) {
        const cfx_rand_desc_t* rng = &g_rand_gens[i];
        printf("  %s\n", rng->name);

        test_rng_reproducible(rng);
        test_rng_different_seeds_differ(rng);
        test_rng_non_degeneracy(rng);
        test_rng_bytes_api(rng);
        test_rng_zero_seed(rng);
    }
}

static void test_explicit_state_independence(void) {
    /* two separate states should produce independent sequences */
    uint64_t s1 = 0x0123456789ABCDEFull;
    uint64_t s2 = 0x0123456789ABCDEFull;
    uint64_t s3 = 0xFEDCBA9876543210ull;

    uint64_t a = cfx_splitmix64(&s1);
    uint64_t b = cfx_splitmix64(&s2);
    CFX_ASSERT(a == b);  /* same seed -> same output */

    uint64_t c = cfx_splitmix64(&s3);
    CFX_ASSERT(a != c);  /* different seed -> different output */

    /* advancing s1 shouldn't affect s2 */
    cfx_splitmix64(&s1);
    cfx_splitmix64(&s1);
    uint64_t d = cfx_splitmix64(&s2);
    uint64_t e = cfx_splitmix64(&s2);
    CFX_ASSERT(d != a);  /* s2 advanced */
    CFX_ASSERT(d != e);

    /* test a few other explicit-state rngs */
    uint32_t x1 = 12345, x2 = 12345;
    CFX_ASSERT(cfx_xorshift32(&x1) == cfx_xorshift32(&x2));

    uint32_t p1 = 42, p2 = 42;
    CFX_ASSERT(cfx_splitmix32(&p1) == cfx_splitmix32(&p2));
}

static void test_rand_limb(void) {
    cfx_srand(999u);

    cfx_limb_t or_acc = 0;
    cfx_limb_t and_acc = CFX_LIMB_MAX;

    for (int i = 0; i < 100; ++i) {
        cfx_limb_t v = cfx_rand_limb();
        or_acc |= v;
        and_acc &= v;
    }

    CFX_ASSERT(or_acc != 0);
    CFX_ASSERT(and_acc != CFX_LIMB_MAX);
}

static void test_cfx_rand_api(void) {
    cfx_srand(12345u);

    /* cfx_rand returns 0..0x7fffffff */
    for (int i = 0; i < 100; ++i) {
        int r = cfx_rand();
        CFX_ASSERT(r >= 0);
    }

    /* cfx_urand returns full 32-bit range */
    uint32_t or_acc = 0;
    for (int i = 0; i < 100; ++i) {
        or_acc |= cfx_urand();
    }
    CFX_ASSERT(or_acc > 0x7FFFFFFFu);  /* should hit high bit sometimes */

    /* cfx_rand_bytes */
    uint8_t buf[32];
    memset(buf, 0, sizeof buf);
    cfx_rand_bytes(buf, sizeof buf);
    uint32_t sum = 0;
    for (size_t i = 0; i < sizeof buf; ++i) sum += buf[i];
    CFX_ASSERT(sum != 0);
}

static void test_chacha20_rng_ctx(void) {
    cfx_rng_ctx_t ctx1, ctx2;

    cfx_chacha20_rng_init(&ctx1, 42);
    cfx_chacha20_rng_init(&ctx2, 42);

    uint8_t buf1[64], buf2[64];
    cfx_chacha20_rng(&ctx1, buf1, sizeof buf1);
    cfx_chacha20_rng(&ctx2, buf2, sizeof buf2);

    CFX_ASSERT(memcmp(buf1, buf2, sizeof buf1) == 0);

    /* different seed -> different output */
    cfx_chacha20_rng_init(&ctx2, 43);
    cfx_chacha20_rng(&ctx2, buf2, sizeof buf2);
    CFX_ASSERT(memcmp(buf1, buf2, sizeof buf1) != 0);
}

/* test cfx_srand_os produces non-zero output */
static void test_srand_os_not_zero(void) {
    uint8_t buf[64];
    memset(buf, 0, sizeof buf);

    cfx_srand_os();
    cfx_rand_bytes(buf, sizeof buf);

    uint32_t sum = 0;
    for (size_t i = 0; i < sizeof buf; i++) sum += buf[i];
    CFX_ASSERT(sum != 0);
}

/* test cfx_srand_os works without calling cfx_srand first */
static void test_srand_os_standalone(void) {
    /* this is a fresh test - don't assume any prior seeding */
    uint8_t buf1[32], buf2[32];

    cfx_srand_os();
    cfx_rand_bytes(buf1, sizeof buf1);

    cfx_srand_os();
    cfx_rand_bytes(buf2, sizeof buf2);

    /* both should be non-zero */
    uint32_t sum1 = 0, sum2 = 0;
    for (size_t i = 0; i < sizeof buf1; i++) {
        sum1 += buf1[i];
        sum2 += buf2[i];
    }
    CFX_ASSERT(sum1 != 0);
    CFX_ASSERT(sum2 != 0);

    /* should be different - (overwhelmingly likely with OS entropy) */
    CFX_ASSERT(memcmp(buf1, buf2, sizeof buf1) != 0);
}

static void test_rand_bytes_os(void) {
    uint8_t buf[32];
    memset(buf, 0, sizeof buf);

    cfx_rand_bytes_os(buf, sizeof buf);

    uint32_t sum = 0;
    for (size_t i = 0; i < sizeof buf; i++) sum += buf[i];
    CFX_ASSERT(sum != 0);
}

/* test various buffer sizes with cfx_srand_os (catches SIMD alignment bugs) */
static void test_srand_os_various_sizes(void) {
    size_t sizes[] = {1, 7, 8, 15, 16, 31, 32, 63, 64, 65, 127, 128, 256, 512, 1024};
    size_t num_sizes = sizeof(sizes) / sizeof(sizes[0]);

    for (size_t i = 0; i < num_sizes; i++) {
        size_t len = sizes[i];
        uint8_t* buf = (uint8_t*)malloc(len);
        CFX_ASSERT(buf != NULL);
        memset(buf, 0, len);

        cfx_srand_os();
        cfx_rand_bytes(buf, len);

        /* check not all zeros (for any reasonable size) */
        uint32_t sum = 0;
        for (size_t j = 0; j < len; j++) sum += buf[j];
        CFX_ASSERT(sum != 0);

        free(buf);
    }
}

int main(void) {
    CFX_TEST(run_all_table_rng_tests);
    CFX_TEST(test_explicit_state_independence);
    CFX_TEST(test_rand_limb);
    CFX_TEST(test_cfx_rand_api);
    CFX_TEST(test_chacha20_rng_ctx);
    CFX_TEST(test_srand_os_not_zero);
    CFX_TEST(test_srand_os_standalone);
    CFX_TEST(test_rand_bytes_os);
    CFX_TEST(test_srand_os_various_sizes);
    puts("OK");
    return 0;
}
