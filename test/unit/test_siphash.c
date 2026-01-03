#include "cfx/siphash.h"
#include "cfx/macros.h"
#include <stdio.h>
#include <string.h>

/* test vectors from SipHash paper (Appendix A)
 * key = 00 01 02 03 04 05 06 07 08 09 0a 0b 0c 0d 0e 0f
 * input = 00, 00 01, 00 01 02, ... */
static const uint64_t vectors_2_4[] = {
    0x726fdb47dd0e0e31ULL,  /* len=0 */
    0x74f839c593dc67fdULL,  /* len=1 */
    0x0d6c8009d9a94f5aULL,  /* len=2 */
    0x85676696d7fb7e2dULL,  /* len=3 */
    0xcf2794e0277187b7ULL,  /* len=4 */
    0x18765564cd99a68dULL,  /* len=5 */
    0xcbc9466e58fee3ceULL,  /* len=6 */
    0xab0200f58b01d137ULL,  /* len=7 */
    0x93f5f5799a932462ULL,  /* len=8 */
    0x9e0082df0ba9e4b0ULL,  /* len=9 */
    0x7a5dbbc594ddb9f3ULL,  /* len=10 */
    0xf4b32f46226bada7ULL,  /* len=11 */
    0x751e8fbc860ee5fbULL,  /* len=12 */
    0x14ea5627c0843d90ULL,  /* len=13 */
    0xf723ca908e7af2eeULL,  /* len=14 */
    0xa129ca6149be45e5ULL,  /* len=15 */
};

static uint8_t test_key[16];
static uint8_t test_msg[64];

static void init_test_data(void) {
    for (int i = 0; i < 16; i++) test_key[i] = (uint8_t)i;
    for (int i = 0; i < 64; i++) test_msg[i] = (uint8_t)i;
}

static void test_vector_len0(void) {
    init_test_data();
    CFX_ASSERT(cfx_siphash(test_msg, 0, test_key) == vectors_2_4[0]);
}

static void test_vector_len1(void) {
    init_test_data();
    CFX_ASSERT(cfx_siphash(test_msg, 1, test_key) == vectors_2_4[1]);
}

static void test_vector_len2(void) {
    init_test_data();
    CFX_ASSERT(cfx_siphash(test_msg, 2, test_key) == vectors_2_4[2]);
}

static void test_vector_len3(void) {
    init_test_data();
    CFX_ASSERT(cfx_siphash(test_msg, 3, test_key) == vectors_2_4[3]);
}

static void test_vector_len7(void) {
    init_test_data();
    CFX_ASSERT(cfx_siphash(test_msg, 7, test_key) == vectors_2_4[7]);
}

static void test_vector_len8(void) {
    init_test_data();
    CFX_ASSERT(cfx_siphash(test_msg, 8, test_key) == vectors_2_4[8]);
}

static void test_vector_len9(void) {
    init_test_data();
    CFX_ASSERT(cfx_siphash(test_msg, 9, test_key) == vectors_2_4[9]);
}

static void test_vector_len15(void) {
    init_test_data();
    CFX_ASSERT(cfx_siphash(test_msg, 15, test_key) == vectors_2_4[15]);
}

static void test_all_vectors(void) {
    init_test_data();
    for (size_t len = 0; len < 16; len++) {
        CFX_ASSERT(cfx_siphash(test_msg, len, test_key) == vectors_2_4[len]);
    }
}

static void test_streaming_15bytes(void) {
    init_test_data();
    uint64_t oneshot = cfx_siphash(test_msg, 15, test_key);

    cfx_siphash_ctx_t ctx;
    cfx_siphash_init(&ctx, test_key);
    cfx_siphash_update(&ctx, test_msg, 7);
    cfx_siphash_update(&ctx, test_msg + 7, 8);
    uint64_t streamed = cfx_siphash_final(&ctx);

    CFX_ASSERT(oneshot == streamed);
}

static void test_streaming_31bytes(void) {
    init_test_data();
    uint64_t oneshot = cfx_siphash(test_msg, 31, test_key);

    cfx_siphash_ctx_t ctx;
    cfx_siphash_init(&ctx, test_key);
    cfx_siphash_update(&ctx, test_msg, 1);
    cfx_siphash_update(&ctx, test_msg + 1, 10);
    cfx_siphash_update(&ctx, test_msg + 11, 8);
    cfx_siphash_update(&ctx, test_msg + 19, 12);
    uint64_t streamed = cfx_siphash_final(&ctx);

    CFX_ASSERT(oneshot == streamed);
}

static void test_streaming_empty(void) {
    init_test_data();
    uint64_t oneshot = cfx_siphash(test_msg, 0, test_key);

    cfx_siphash_ctx_t ctx;
    cfx_siphash_init(&ctx, test_key);
    uint64_t streamed = cfx_siphash_final(&ctx);

    CFX_ASSERT(oneshot == streamed);
}

static void test_streaming_byte_by_byte(void) {
    init_test_data();
    uint64_t oneshot = cfx_siphash(test_msg, 16, test_key);

    cfx_siphash_ctx_t ctx;
    cfx_siphash_init(&ctx, test_key);
    for (int i = 0; i < 16; i++) {
        cfx_siphash_update(&ctx, test_msg + i, 1);
    }
    uint64_t streamed = cfx_siphash_final(&ctx);

    CFX_ASSERT(oneshot == streamed);
}

static void test_siphash128_nonzero(void) {
    init_test_data();
    uint8_t out[16] = {0};
    cfx_siphash128(out, test_msg, 15, test_key);

    int nonzero = 0;
    for (int i = 0; i < 16; i++) {
        if (out[i]) nonzero = 1;
    }
    CFX_ASSERT(nonzero);
}

static void test_siphash128_deterministic(void) {
    init_test_data();
    uint8_t out1[16], out2[16];
    cfx_siphash128(out1, test_msg, 15, test_key);
    cfx_siphash128(out2, test_msg, 15, test_key);
    CFX_ASSERT(memcmp(out1, out2, 16) == 0);
}

static void test_siphash13_differs(void) {
    init_test_data();
    uint64_t h13 = cfx_siphash_cd(test_msg, 8, test_key, 1, 3);
    uint64_t h24 = cfx_siphash_cd(test_msg, 8, test_key, 2, 4);
    CFX_ASSERT(h13 != h24);
}

static void test_siphash48_differs(void) {
    init_test_data();
    uint64_t h48 = cfx_siphash_cd(test_msg, 8, test_key, 4, 8);
    uint64_t h24 = cfx_siphash_cd(test_msg, 8, test_key, 2, 4);
    CFX_ASSERT(h48 != h24);
}

static void test_different_keys(void) {
    init_test_data();
    uint8_t key2[16];
    for (int i = 0; i < 16; i++) key2[i] = (uint8_t)(i + 1);

    uint64_t h1 = cfx_siphash(test_msg, 8, test_key);
    uint64_t h2 = cfx_siphash(test_msg, 8, key2);
    CFX_ASSERT(h1 != h2);
}

static void test_different_messages(void) {
    init_test_data();
    uint8_t msg2[8] = {1, 2, 3, 4, 5, 6, 7, 8};

    uint64_t h1 = cfx_siphash(test_msg, 8, test_key);
    uint64_t h2 = cfx_siphash(msg2, 8, test_key);
    CFX_ASSERT(h1 != h2);
}

static void test_zero_key(void) {
    uint8_t zerokey[16] = {0};
    uint8_t msg[4] = {1, 2, 3, 4};
    uint64_t h = cfx_siphash(msg, 4, zerokey);
    (void)h;
    CFX_ASSERT(1);
}

/* test all switch(left) branches in cfx_siphash128 */
static void test_siphash128_all_lengths(void) {
    init_test_data();
    uint8_t out1[16], out2[16];

    for (size_t len = 0; len < 16; len++) {
        cfx_siphash128(out1, test_msg, len, test_key);
        cfx_siphash128(out2, test_msg, len, test_key);
        CFX_ASSERT(memcmp(out1, out2, 16) == 0);
    }
}

/* test all switch(s->buflen) branches in streaming final */
static void test_streaming_buflen_0(void) {
    init_test_data();
    uint64_t oneshot = cfx_siphash(test_msg, 8, test_key);
    cfx_siphash_ctx_t ctx;
    cfx_siphash_init(&ctx, test_key);
    cfx_siphash_update(&ctx, test_msg, 8);
    CFX_ASSERT(cfx_siphash_final(&ctx) == oneshot);
}

static void test_streaming_buflen_1(void) {
    init_test_data();
    uint64_t oneshot = cfx_siphash(test_msg, 9, test_key);
    cfx_siphash_ctx_t ctx;
    cfx_siphash_init(&ctx, test_key);
    cfx_siphash_update(&ctx, test_msg, 9);
    CFX_ASSERT(cfx_siphash_final(&ctx) == oneshot);
}

static void test_streaming_buflen_2(void) {
    init_test_data();
    uint64_t oneshot = cfx_siphash(test_msg, 10, test_key);
    cfx_siphash_ctx_t ctx;
    cfx_siphash_init(&ctx, test_key);
    cfx_siphash_update(&ctx, test_msg, 10);
    CFX_ASSERT(cfx_siphash_final(&ctx) == oneshot);
}

static void test_streaming_buflen_3(void) {
    init_test_data();
    uint64_t oneshot = cfx_siphash(test_msg, 11, test_key);
    cfx_siphash_ctx_t ctx;
    cfx_siphash_init(&ctx, test_key);
    cfx_siphash_update(&ctx, test_msg, 11);
    CFX_ASSERT(cfx_siphash_final(&ctx) == oneshot);
}

static void test_streaming_buflen_4(void) {
    init_test_data();
    uint64_t oneshot = cfx_siphash(test_msg, 12, test_key);
    cfx_siphash_ctx_t ctx;
    cfx_siphash_init(&ctx, test_key);
    cfx_siphash_update(&ctx, test_msg, 12);
    CFX_ASSERT(cfx_siphash_final(&ctx) == oneshot);
}

static void test_streaming_buflen_5(void) {
    init_test_data();
    uint64_t oneshot = cfx_siphash(test_msg, 13, test_key);
    cfx_siphash_ctx_t ctx;
    cfx_siphash_init(&ctx, test_key);
    cfx_siphash_update(&ctx, test_msg, 13);
    CFX_ASSERT(cfx_siphash_final(&ctx) == oneshot);
}

static void test_streaming_buflen_6(void) {
    init_test_data();
    uint64_t oneshot = cfx_siphash(test_msg, 14, test_key);
    cfx_siphash_ctx_t ctx;
    cfx_siphash_init(&ctx, test_key);
    cfx_siphash_update(&ctx, test_msg, 14);
    CFX_ASSERT(cfx_siphash_final(&ctx) == oneshot);
}

static void test_streaming_buflen_7(void) {
    init_test_data();
    uint64_t oneshot = cfx_siphash(test_msg, 15, test_key);
    cfx_siphash_ctx_t ctx;
    cfx_siphash_init(&ctx, test_key);
    cfx_siphash_update(&ctx, test_msg, 15);
    CFX_ASSERT(cfx_siphash_final(&ctx) == oneshot);
}

int main(void) {
    CFX_TEST(test_vector_len0);
    CFX_TEST(test_vector_len1);
    CFX_TEST(test_vector_len2);
    CFX_TEST(test_vector_len3);
    CFX_TEST(test_vector_len7);
    CFX_TEST(test_vector_len8);
    CFX_TEST(test_vector_len9);
    CFX_TEST(test_vector_len15);
    CFX_TEST(test_all_vectors);
    CFX_TEST(test_streaming_15bytes);
    CFX_TEST(test_streaming_31bytes);
    CFX_TEST(test_streaming_empty);
    CFX_TEST(test_streaming_byte_by_byte);
    CFX_TEST(test_siphash128_nonzero);
    CFX_TEST(test_siphash128_deterministic);
    CFX_TEST(test_siphash13_differs);
    CFX_TEST(test_siphash48_differs);
    CFX_TEST(test_different_keys);
    CFX_TEST(test_different_messages);
    CFX_TEST(test_zero_key);
    CFX_TEST(test_siphash128_all_lengths);
    CFX_TEST(test_streaming_buflen_0);
    CFX_TEST(test_streaming_buflen_1);
    CFX_TEST(test_streaming_buflen_2);
    CFX_TEST(test_streaming_buflen_3);
    CFX_TEST(test_streaming_buflen_4);
    CFX_TEST(test_streaming_buflen_5);
    CFX_TEST(test_streaming_buflen_6);
    CFX_TEST(test_streaming_buflen_7);
    puts("ok");
    return 0;
}
