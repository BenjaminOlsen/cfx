/* test_rand.c */

#include <assert.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>

#include "cfx/macros.h"
#include "cfx/rand.h"

#define N 16

static void test_rng_reproducible(const cfx_rand_desc_t *rng) {
    const uint32_t SEED = 0x12345678u;

    uint32_t seq1[N];
    uint32_t seq2[N];

    rng->seed(SEED);
    for (size_t i = 0; i < N; ++i) {
        seq1[i] = rng->gen32();
    }

    rng->seed(SEED);
    for (size_t i = 0; i < N; ++i) {
        seq2[i] = rng->gen32();
    }

    assert(memcmp(seq1, seq2, sizeof seq1) == 0);
}

/* Very weak, but catches “seed is ignored” type bugs */
static void test_rng_different_seeds_differ(const cfx_rand_desc_t* rng) {

    uint32_t seq1[N];
    uint32_t seq2[N];

    rng->seed(1u);
    for (size_t i = 0; i < N; ++i) {
        seq1[i] = rng->gen32();
    }

    rng->seed(2u);
    for (size_t i = 0; i < N; ++i) {
        seq2[i] = rng->gen32();
    }

    size_t equal = 0;
    for (size_t i = 0; i < N; ++i) {
        if (seq1[i] == seq2[i]) ++equal;
    }

    /* Allow a couple of accidental collisions, but not all of them */
    assert(equal < N);
}

static void run_all_table_rng_tests(void) {
    const size_t count = sizeof g_rand_gen_cnt;

    for (size_t i = 0; i < count; ++i) {
        const cfx_rand_desc_t* rng = &g_rand_gens[i];
        printf("Testing RNG: %s\n", rng->name);

        test_rng_reproducible(rng);
        test_rng_different_seeds_differ(rng);
    }
}

static void test_rng_state_advance(void) {
    uint64_t s = 0x0123456789ABCDEFull;
    uint64_t s2 = s;

    uint64_t a = cfx_splitmix64(&s);
    uint64_t b = cfx_splitmix64(&s);
    uint64_t c = cfx_splitmix64(&s2);

    /* Advancing one step from the same seed twice should match */
    assert(a == c);
    /* And state must advance, so second sample differs with overwhelming prob. */
    assert(a != b);
}

int main(void) {
    CFX_TEST(run_all_table_rng_tests);
    CFX_TEST(test_rng_state_advance);
    return 0;
}


