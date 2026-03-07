/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#include "cfx/aead_chacha20_poly1305.h"
#include "cfx/macros.h"
#include "cfx/rand.h"

#include <string.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>


#define PR printf

static void print_hex(const char *label, const uint8_t *buf, size_t len) {
    printf("%s (len=%zu):\n", label, len);
    for (size_t i = 0; i < len; ++i) {
        printf("%02x", buf[i]);
        if ((i + 1) % 16 == 0 || i + 1 == len)
            printf("\n");
        else
            printf(" ");
    }
}

/* -- test vectors from RFC 8439 2.8.2: Example and Test Vector for AEAD_CHACHA20_POLY1305 -- */
static const uint8_t KEY[32] = {
    0x80,0x81,0x82,0x83,0x84,0x85,0x86,0x87,
    0x88,0x89,0x8a,0x8b,0x8c,0x8d,0x8e,0x8f,
    0x90,0x91,0x92,0x93,0x94,0x95,0x96,0x97,
    0x98,0x99,0x9a,0x9b,0x9c,0x9d,0x9e,0x9f
};

static const uint8_t NONCE[12] = {
    0x07,0x00,0x00,0x00,
    0x40,0x41,0x42,0x43,0x44,0x45,0x46,0x47
};

static const uint8_t AAD[12] = {
    0x50,0x51,0x52,0x53,
    0xc0,0xc1,0xc2,0xc3,
    0xc4,0xc5,0xc6,0xc7
};

static const uint8_t PT[114] = {
    0x4c,0x61,0x64,0x69,0x65,0x73,0x20,0x61,0x6e,0x64,0x20,0x47,0x65,0x6e,0x74,0x6c,
    0x65,0x6d,0x65,0x6e,0x20,0x6f,0x66,0x20,0x74,0x68,0x65,0x20,0x63,0x6c,0x61,0x73,
    0x73,0x20,0x6f,0x66,0x20,0x27,0x39,0x39,0x3a,0x20,0x49,0x66,0x20,0x49,0x20,0x63,
    0x6f,0x75,0x6c,0x64,0x20,0x6f,0x66,0x66,0x65,0x72,0x20,0x79,0x6f,0x75,0x20,0x6f,
    0x6e,0x6c,0x79,0x20,0x6f,0x6e,0x65,0x20,0x74,0x69,0x70,0x20,0x66,0x6f,0x72,0x20,
    0x74,0x68,0x65,0x20,0x66,0x75,0x74,0x75,0x72,0x65,0x2c,0x20,0x73,0x75,0x6e,0x73,
    0x63,0x72,0x65,0x65,0x6e,0x20,0x77,0x6f,0x75,0x6c,0x64,0x20,0x62,0x65,0x20,0x69,
    0x74,0x2e
};

static const uint8_t CT_EXPECTED[114] = {
    0xd3,0x1a,0x8d,0x34,0x64,0x8e,0x60,0xdb,0x7b,0x86,0xaf,0xbc,0x53,0xef,0x7e,0xc2,
    0xa4,0xad,0xed,0x51,0x29,0x6e,0x08,0xfe,0xa9,0xe2,0xb5,0xa7,0x36,0xee,0x62,0xd6,
    0x3d,0xbe,0xa4,0x5e,0x8c,0xa9,0x67,0x12,0x82,0xfa,0xfb,0x69,0xda,0x92,0x72,0x8b,
    0x1a,0x71,0xde,0x0a,0x9e,0x06,0x0b,0x29,0x05,0xd6,0xa5,0xb6,0x7e,0xcd,0x3b,0x36,
    0x92,0xdd,0xbd,0x7f,0x2d,0x77,0x8b,0x8c,0x98,0x03,0xae,0xe3,0x28,0x09,0x1b,0x58,
    0xfa,0xb3,0x24,0xe4,0xfa,0xd6,0x75,0x94,0x55,0x85,0x80,0x8b,0x48,0x31,0xd7,0xbc,
    0x3f,0xf4,0xde,0xf0,0x8e,0x4b,0x7a,0x9d,0xe5,0x76,0xd2,0x65,0x86,0xce,0xc6,0x4b,
    0x61,0x16
};

static const uint8_t TAG_EXPECTED[16] = {
    0x1a,0xe1,0x0b,0x59,0x4f,0x09,0xe2,0x6a,
    0x7e,0x90,0x2e,0xcb,0xd0,0x60,0x06,0x91
};

static void test_rfc8439_encrypt(void) {
    uint8_t ct[sizeof PT];
    uint8_t tag[16];

    printf("== test_rfc8439_encrypt ==\n");
    print_hex("key",   KEY,   sizeof KEY);
    print_hex("nonce", NONCE, sizeof NONCE);
    print_hex("aad",   AAD,   sizeof AAD);
    print_hex("pt",    PT,    sizeof PT);

    int rc = cfx_chacha20_poly1305_encrypt(
        ct, tag,
        PT, sizeof PT,
        AAD, sizeof AAD,
        KEY, NONCE
        );
    CFX_ASSERT(rc == 0);

    print_hex("ct (got)", ct, sizeof ct);
    print_hex("ct (exp)", CT_EXPECTED, sizeof CT_EXPECTED);
    print_hex("tag (got)", tag, 16);
    print_hex("tag (exp)", TAG_EXPECTED, 16);

    CFX_ASSERT(memcmp(ct,  CT_EXPECTED, sizeof ct)  == 0);
    CFX_ASSERT(memcmp(tag, TAG_EXPECTED, 16)       == 0);
}

static void test_rfc8439_decrypt(void) {
    uint8_t pt_out[sizeof PT];

    printf("== test_rfc8439_decrypt ==\n");
    print_hex("ct (in)",  CT_EXPECTED, sizeof CT_EXPECTED);
    print_hex("tag (in)", TAG_EXPECTED, 16);

    int rc = cfx_chacha20_poly1305_decrypt(
        pt_out,
        CT_EXPECTED, sizeof CT_EXPECTED,
        AAD, sizeof AAD,
        KEY, NONCE,
        TAG_EXPECTED
        );
    CFX_ASSERT(rc == 0);

    print_hex("pt_out", pt_out, sizeof pt_out);
    print_hex("pt_exp", PT,   sizeof PT);

    CFX_ASSERT(memcmp(pt_out, PT, sizeof PT) == 0);
}

static void test_rfc8439_bad_tag(void) {
    uint8_t pt_out[sizeof PT];
    uint8_t bad_tag[16];

    memcpy(bad_tag, TAG_EXPECTED, 16);
    bad_tag[0] ^= 0x01;

    printf("== test_rfc8439_bad_tag ==\n");
    print_hex("ct",      CT_EXPECTED, sizeof CT_EXPECTED);
    print_hex("tag bad", bad_tag,       16);

    int rc = cfx_chacha20_poly1305_decrypt(
        pt_out,
        CT_EXPECTED, sizeof CT_EXPECTED,
        AAD, sizeof AAD,
        KEY, NONCE,
        bad_tag
        );
    CFX_ASSERT(rc != 0);
}


static void fuzz_fill(uint8_t *buf, size_t len) {
    for (size_t i = 0; i < len; ++i) {
        buf[i] = (uint8_t)cfx_rand();
    }
}

/* draft-irtf-cfrg-xchacha A.3.1 test vector */
static const uint8_t XCHACHA_KEY[32] = {
    0x80,0x81,0x82,0x83,0x84,0x85,0x86,0x87,
    0x88,0x89,0x8a,0x8b,0x8c,0x8d,0x8e,0x8f,
    0x90,0x91,0x92,0x93,0x94,0x95,0x96,0x97,
    0x98,0x99,0x9a,0x9b,0x9c,0x9d,0x9e,0x9f
};

static const uint8_t XCHACHA_NONCE[24] = {
    0x40,0x41,0x42,0x43,0x44,0x45,0x46,0x47,
    0x48,0x49,0x4a,0x4b,0x4c,0x4d,0x4e,0x4f,
    0x50,0x51,0x52,0x53,0x54,0x55,0x56,0x57
};

static const uint8_t XCHACHA_AAD[12] = {
    0x50,0x51,0x52,0x53,0xc0,0xc1,0xc2,0xc3,
    0xc4,0xc5,0xc6,0xc7
};

static const uint8_t XCHACHA_PT[114] = {
    0x4c,0x61,0x64,0x69,0x65,0x73,0x20,0x61,0x6e,0x64,0x20,0x47,0x65,0x6e,0x74,0x6c,
    0x65,0x6d,0x65,0x6e,0x20,0x6f,0x66,0x20,0x74,0x68,0x65,0x20,0x63,0x6c,0x61,0x73,
    0x73,0x20,0x6f,0x66,0x20,0x27,0x39,0x39,0x3a,0x20,0x49,0x66,0x20,0x49,0x20,0x63,
    0x6f,0x75,0x6c,0x64,0x20,0x6f,0x66,0x66,0x65,0x72,0x20,0x79,0x6f,0x75,0x20,0x6f,
    0x6e,0x6c,0x79,0x20,0x6f,0x6e,0x65,0x20,0x74,0x69,0x70,0x20,0x66,0x6f,0x72,0x20,
    0x74,0x68,0x65,0x20,0x66,0x75,0x74,0x75,0x72,0x65,0x2c,0x20,0x73,0x75,0x6e,0x73,
    0x63,0x72,0x65,0x65,0x6e,0x20,0x77,0x6f,0x75,0x6c,0x64,0x20,0x62,0x65,0x20,0x69,
    0x74,0x2e
};

static const uint8_t XCHACHA_CT_EXPECTED[114] = {
    0xbd,0x6d,0x17,0x9d,0x3e,0x83,0xd4,0x3b,0x95,0x76,0x57,0x94,0x93,0xc0,0xe9,0x39,
    0x57,0x2a,0x17,0x00,0x25,0x2b,0xfa,0xcc,0xbe,0xd2,0x90,0x2c,0x21,0x39,0x6c,0xbb,
    0x73,0x1c,0x7f,0x1b,0x0b,0x4a,0xa6,0x44,0x0b,0xf3,0xa8,0x2f,0x4e,0xda,0x7e,0x39,
    0xae,0x64,0xc6,0x70,0x8c,0x54,0xc2,0x16,0xcb,0x96,0xb7,0x2e,0x12,0x13,0xb4,0x52,
    0x2f,0x8c,0x9b,0xa4,0x0d,0xb5,0xd9,0x45,0xb1,0x1b,0x69,0xb9,0x82,0xc1,0xbb,0x9e,
    0x3f,0x3f,0xac,0x2b,0xc3,0x69,0x48,0x8f,0x76,0xb2,0x38,0x35,0x65,0xd3,0xff,0xf9,
    0x21,0xf9,0x66,0x4c,0x97,0x63,0x7d,0xa9,0x76,0x88,0x12,0xf6,0x15,0xc6,0x8b,0x13,
    0xb5,0x2e
};

static const uint8_t XCHACHA_TAG_EXPECTED[16] = {
    0xc0,0x87,0x59,0x24,0xc1,0xc7,0x98,0x79,
    0x47,0xde,0xaf,0xd8,0x78,0x0a,0xcf,0x49
};

static void test_xchacha_encrypt(void) {
    uint8_t ct[sizeof XCHACHA_PT];
    uint8_t tag[16];

    printf("== test_xchacha_encrypt ==\n");

    int rc = cfx_xchacha20_poly1305_encrypt(
        ct, tag,
        XCHACHA_PT, sizeof XCHACHA_PT,
        XCHACHA_AAD, sizeof XCHACHA_AAD,
        XCHACHA_KEY, XCHACHA_NONCE
        );
    CFX_ASSERT(rc == 0);

    print_hex("ct (got)", ct, sizeof ct);
    print_hex("ct (exp)", XCHACHA_CT_EXPECTED, sizeof XCHACHA_CT_EXPECTED);
    print_hex("tag (got)", tag, 16);
    print_hex("tag (exp)", XCHACHA_TAG_EXPECTED, 16);

    CFX_ASSERT(memcmp(ct, XCHACHA_CT_EXPECTED, sizeof ct) == 0);
    CFX_ASSERT(memcmp(tag, XCHACHA_TAG_EXPECTED, 16) == 0);
}

static void test_xchacha_decrypt(void) {
    uint8_t pt_out[sizeof XCHACHA_PT];

    printf("== test_xchacha_decrypt ==\n");

    int rc = cfx_xchacha20_poly1305_decrypt(
        pt_out,
        XCHACHA_CT_EXPECTED, sizeof XCHACHA_CT_EXPECTED,
        XCHACHA_AAD, sizeof XCHACHA_AAD,
        XCHACHA_KEY, XCHACHA_NONCE,
        XCHACHA_TAG_EXPECTED
        );
    CFX_ASSERT(rc == 0);
    CFX_ASSERT(memcmp(pt_out, XCHACHA_PT, sizeof XCHACHA_PT) == 0);
}

static void test_xchacha_bad_tag(void) {
    uint8_t pt_out[sizeof XCHACHA_PT];
    uint8_t bad_tag[16];

    memcpy(bad_tag, XCHACHA_TAG_EXPECTED, 16);
    bad_tag[0] ^= 0x01;

    printf("== test_xchacha_bad_tag ==\n");

    int rc = cfx_xchacha20_poly1305_decrypt(
        pt_out,
        XCHACHA_CT_EXPECTED, sizeof XCHACHA_CT_EXPECTED,
        XCHACHA_AAD, sizeof XCHACHA_AAD,
        XCHACHA_KEY, XCHACHA_NONCE,
        bad_tag
        );
    CFX_ASSERT(rc != 0);
}

static void test_xchacha_fuzz(void) {
    #define XFUZZ_MAX_PT 256
    #define XFUZZ_MAX_AAD 64
    #define XFUZZ_ITERS 3000

    uint8_t key[32];
    uint8_t nonce[24];
    uint8_t aad[XFUZZ_MAX_AAD];
    uint8_t pt[XFUZZ_MAX_PT];
    uint8_t ct[XFUZZ_MAX_PT];
    uint8_t pt_out[XFUZZ_MAX_PT];
    uint8_t tag[16];
    uint8_t bad_tag[16];

    printf("== test_xchacha_fuzz ==\n");

    cfx_srand(0x12345678);

    for (int iter = 0; iter < XFUZZ_ITERS; ++iter) {
        for (size_t i = 0; i < sizeof key; ++i) key[i] = (uint8_t)cfx_rand();
        for (size_t i = 0; i < sizeof nonce; ++i) nonce[i] = (uint8_t)cfx_rand();

        size_t pt_len = cfx_rand() % (XFUZZ_MAX_PT + 1);
        size_t aad_len = cfx_rand() % (XFUZZ_MAX_AAD + 1);

        for (size_t i = 0; i < pt_len; ++i) pt[i] = (uint8_t)cfx_rand();
        for (size_t i = 0; i < aad_len; ++i) aad[i] = (uint8_t)cfx_rand();

        const uint8_t *aad_ptr = aad_len ? aad : NULL;

        int rc = cfx_xchacha20_poly1305_encrypt(ct, tag, pt, pt_len, aad_ptr, aad_len, key, nonce);
        CFX_ASSERT(rc == 0);

        rc = cfx_xchacha20_poly1305_decrypt(pt_out, ct, pt_len, aad_ptr, aad_len, key, nonce, tag);
        CFX_ASSERT(rc == 0);
        CFX_ASSERT(memcmp(pt_out, pt, pt_len) == 0);

        memcpy(bad_tag, tag, 16);
        bad_tag[cfx_rand() % 16] ^= (uint8_t)(1u << (cfx_rand() & 7u));

        rc = cfx_xchacha20_poly1305_decrypt(pt_out, ct, pt_len, aad_ptr, aad_len, key, nonce, bad_tag);
        CFX_ASSERT(rc != 0);

        if ((iter % (XFUZZ_ITERS / 10)) == 0) {
            printf("   xchacha fuzz iter %d OK\n", iter);
        }
    }
    printf("-- xchacha fuzz OK\n");
}

static void test_aead_fuzz_basic(void) {
    #define MAX_PT 256
    #define MAX_AAD 64
    #define ITERS 5000

    uint8_t key[32];
    uint8_t nonce[12];
    uint8_t aad[MAX_AAD];
    uint8_t pt[MAX_PT];
    uint8_t ct[MAX_PT];
    uint8_t pt_out[MAX_PT];
    uint8_t tag[16];
    uint8_t bad_tag[16];

    printf("== test_aead_fuzz_basic ==\n");

    cfx_srand(0x828736);
    uint32_t ok_cnt = 0;

    for (int iter = 0; iter < ITERS; ++iter) {
        fuzz_fill(key,   sizeof key);
        fuzz_fill(nonce, sizeof nonce);

        size_t pt_len  = cfx_rand() % (MAX_PT + 1);
        size_t aad_len = cfx_rand() % (MAX_AAD + 1);

        fuzz_fill(pt,  pt_len);
        fuzz_fill(aad, aad_len);

        const uint8_t *aad_ptr = aad;
        if (aad_len == 0 && (cfx_rand() & 1u)) {
            aad_ptr = NULL;
        }

        int rc = cfx_chacha20_poly1305_encrypt(
            ct, tag,
            pt, pt_len,
            aad_ptr, aad_len,
            key, nonce
            );
        CFX_ASSERT(rc == 0);

        rc = cfx_chacha20_poly1305_decrypt(
            pt_out,
            ct, pt_len,
            aad_ptr, aad_len,
            key, nonce,
            tag
            );
        CFX_ASSERT(rc == 0);
        CFX_ASSERT(memcmp(pt_out, pt, pt_len) == 0);

        memcpy(bad_tag, tag, sizeof bad_tag);
        size_t pos   =  cfx_rand() % 16;
        uint8_t bit  = (uint8_t)(1u << (cfx_rand() & 7u));
        bad_tag[pos] ^= bit;

        rc = cfx_chacha20_poly1305_decrypt(
            pt_out,
            ct, pt_len,
            aad_ptr, aad_len,
            key, nonce,
            bad_tag
            );
        CFX_ASSERT(rc != 0);

        ++ok_cnt;
        if ((iter % (ITERS / 20)) == 0) {
            printf("   fuzz iter %d OK (pt_len=%zu, aad_len=%zu)\n", iter, pt_len, aad_len);
        }
    }
    printf("-- fuzz OK cnt: %u\n", ok_cnt);
}

/* --- streaming AEAD tests --- */

static void test_stream_single_chunk(void) {
    /* single small chunk, roundtrip */
    uint8_t key[32], nonce[24];
    uint8_t pt[] = "hello streaming AEAD";
    uint8_t ct[sizeof pt], pt_out[sizeof pt];
    uint8_t tag[16];

    printf("== test_stream_single_chunk ==\n");

    cfx_srand(0xAABBCCDD);
    fuzz_fill(key, 32);
    fuzz_fill(nonce, 24);

    int rc = cfx_stream_xchacha20_poly1305_encrypt_chunk(
        ct, tag, pt, sizeof pt, 0, 1, key, nonce);
    CFX_ASSERT(rc == 0);

    rc = cfx_stream_xchacha20_poly1305_decrypt_chunk(
        pt_out, ct, sizeof pt, tag, 0, 1, key, nonce);
    CFX_ASSERT(rc == 0);
    CFX_ASSERT(memcmp(pt_out, pt, sizeof pt) == 0);
}

#ifndef CFX_BAREMETAL  /* 192KB+ RAM needed — skip on bare-metal */
static void test_stream_multi_chunk(void) {
    /* 3 chunks: 64KB + 64KB + 30 bytes */
    #define MC_CHUNK CFX_STREAM_CHUNK_SIZE
    #define MC_TAIL  30
    #define MC_TOTAL (MC_CHUNK * 2 + MC_TAIL)

    uint8_t key[32], nonce[24];
    uint8_t pt[MC_TOTAL];
    uint8_t ct0[MC_CHUNK], ct1[MC_CHUNK], ct2[MC_TAIL];
    uint8_t tag0[16], tag1[16], tag2[16];
    uint8_t dec[MC_TOTAL];

    printf("== test_stream_multi_chunk ==\n");

    cfx_srand(0x11223344);
    fuzz_fill(key, 32);
    fuzz_fill(nonce, 24);
    fuzz_fill(pt, MC_TOTAL);

    /* encrypt 3 chunks */
    int rc = cfx_stream_xchacha20_poly1305_encrypt_chunk(
        ct0, tag0, pt, MC_CHUNK, 0, 0, key, nonce);
    CFX_ASSERT(rc == 0);

    rc = cfx_stream_xchacha20_poly1305_encrypt_chunk(
        ct1, tag1, pt + MC_CHUNK, MC_CHUNK, 1, 0, key, nonce);
    CFX_ASSERT(rc == 0);

    rc = cfx_stream_xchacha20_poly1305_encrypt_chunk(
        ct2, tag2, pt + MC_CHUNK * 2, MC_TAIL, 2, 1, key, nonce);
    CFX_ASSERT(rc == 0);

    /* decrypt in order */
    rc = cfx_stream_xchacha20_poly1305_decrypt_chunk(
        dec, ct0, MC_CHUNK, tag0, 0, 0, key, nonce);
    CFX_ASSERT(rc == 0);

    rc = cfx_stream_xchacha20_poly1305_decrypt_chunk(
        dec + MC_CHUNK, ct1, MC_CHUNK, tag1, 1, 0, key, nonce);
    CFX_ASSERT(rc == 0);

    rc = cfx_stream_xchacha20_poly1305_decrypt_chunk(
        dec + MC_CHUNK * 2, ct2, MC_TAIL, tag2, 2, 1, key, nonce);
    CFX_ASSERT(rc == 0);

    CFX_ASSERT(memcmp(dec, pt, MC_TOTAL) == 0);
}
#endif /* CFX_BAREMETAL */

static void test_stream_empty(void) {
    /* empty final chunk */
    uint8_t key[32], nonce[24];
    uint8_t ct[1], pt_out[1], tag[16];

    printf("== test_stream_empty ==\n");

    cfx_srand(0x55667788);
    fuzz_fill(key, 32);
    fuzz_fill(nonce, 24);

    int rc = cfx_stream_xchacha20_poly1305_encrypt_chunk(
        ct, tag, NULL, 0, 0, 1, key, nonce);
    CFX_ASSERT(rc == 0);

    rc = cfx_stream_xchacha20_poly1305_decrypt_chunk(
        pt_out, ct, 0, tag, 0, 1, key, nonce);
    CFX_ASSERT(rc == 0);
}

static void test_stream_wrong_order(void) {
    /* swapping chunk 0 and chunk 1 must fail */
    uint8_t key[32], nonce[24];
    uint8_t pt0[32], pt1[16];
    uint8_t ct0[32], ct1[16], tag0[16], tag1[16];
    uint8_t dec[32];

    printf("== test_stream_wrong_order ==\n");

    cfx_srand(0xDEADBEEF);
    fuzz_fill(key, 32);
    fuzz_fill(nonce, 24);
    fuzz_fill(pt0, 32);
    fuzz_fill(pt1, 16);

    cfx_stream_xchacha20_poly1305_encrypt_chunk(ct0, tag0, pt0, 32, 0, 0, key, nonce);
    cfx_stream_xchacha20_poly1305_encrypt_chunk(ct1, tag1, pt1, 16, 1, 1, key, nonce);

    /* try decrypting chunk 1 data at position 0 — must fail */
    int rc = cfx_stream_xchacha20_poly1305_decrypt_chunk(
        dec, ct1, 16, tag1, 0, 0, key, nonce);
    CFX_ASSERT(rc != 0);

    /* try decrypting chunk 0 data at position 1 — must fail */
    rc = cfx_stream_xchacha20_poly1305_decrypt_chunk(
        dec, ct0, 32, tag0, 1, 1, key, nonce);
    CFX_ASSERT(rc != 0);
}

static void test_stream_truncation(void) {
    /* decrypting a non-final chunk as final must fail */
    uint8_t key[32], nonce[24];
    uint8_t pt[64];
    uint8_t ct[64], tag[16];
    uint8_t dec[64];

    printf("== test_stream_truncation ==\n");

    cfx_srand(0xCAFEBABE);
    fuzz_fill(key, 32);
    fuzz_fill(nonce, 24);
    fuzz_fill(pt, 64);

    /* encrypt as non-final (chunk 0 of a multi-chunk stream) */
    cfx_stream_xchacha20_poly1305_encrypt_chunk(ct, tag, pt, 64, 0, 0, key, nonce);

    /* try to decrypt as final — must fail (truncation attack) */
    int rc = cfx_stream_xchacha20_poly1305_decrypt_chunk(
        dec, ct, 64, tag, 0, 1, key, nonce);
    CFX_ASSERT(rc != 0);
}

static void test_stream_tamper(void) {
    /* flipping a bit in ciphertext must fail */
    uint8_t key[32], nonce[24];
    uint8_t pt[100];
    uint8_t ct[100], tag[16];
    uint8_t dec[100];

    printf("== test_stream_tamper ==\n");

    cfx_srand(0xF00DBABE);
    fuzz_fill(key, 32);
    fuzz_fill(nonce, 24);
    fuzz_fill(pt, 100);

    cfx_stream_xchacha20_poly1305_encrypt_chunk(ct, tag, pt, 100, 0, 1, key, nonce);

    /* tamper with ciphertext */
    ct[50] ^= 0x01;
    int rc = cfx_stream_xchacha20_poly1305_decrypt_chunk(
        dec, ct, 100, tag, 0, 1, key, nonce);
    CFX_ASSERT(rc != 0);

    /* restore ct, tamper with tag */
    ct[50] ^= 0x01;
    tag[8] ^= 0x01;
    rc = cfx_stream_xchacha20_poly1305_decrypt_chunk(
        dec, ct, 100, tag, 0, 1, key, nonce);
    CFX_ASSERT(rc != 0);
}

#ifndef CFX_BAREMETAL  /* 3 × 64KB malloc — skip on bare-metal */
static void test_stream_exact_chunk(void) {
    /* single chunk that is exactly CFX_STREAM_CHUNK_SIZE (64KB) */
    uint8_t key[32], nonce[24];
    uint8_t *pt  = (uint8_t *)malloc(CFX_STREAM_CHUNK_SIZE);
    uint8_t *ct  = (uint8_t *)malloc(CFX_STREAM_CHUNK_SIZE);
    uint8_t *dec = (uint8_t *)malloc(CFX_STREAM_CHUNK_SIZE);
    uint8_t tag[16];

    printf("== test_stream_exact_chunk ==\n");

    CFX_ASSERT(pt && ct && dec);

    cfx_srand(0xA1B2C3D4);
    fuzz_fill(key, 32);
    fuzz_fill(nonce, 24);
    fuzz_fill(pt, CFX_STREAM_CHUNK_SIZE);

    int rc = cfx_stream_xchacha20_poly1305_encrypt_chunk(
        ct, tag, pt, CFX_STREAM_CHUNK_SIZE, 0, 1, key, nonce);
    CFX_ASSERT(rc == 0);

    rc = cfx_stream_xchacha20_poly1305_decrypt_chunk(
        dec, ct, CFX_STREAM_CHUNK_SIZE, tag, 0, 1, key, nonce);
    CFX_ASSERT(rc == 0);
    CFX_ASSERT(memcmp(dec, pt, CFX_STREAM_CHUNK_SIZE) == 0);

    free(pt); free(ct); free(dec);
}
#endif /* CFX_BAREMETAL */

#ifndef CFX_BAREMETAL  /* 4 × 64KB malloc — skip on bare-metal */
static void test_stream_exact_two_chunks(void) {
    /* two full 64KB chunks, no tail — both at exact chunk boundary */
    #define E2C CFX_STREAM_CHUNK_SIZE
    uint8_t key[32], nonce[24];
    uint8_t *pt  = (uint8_t *)malloc(E2C * 2);
    uint8_t *ct0 = (uint8_t *)malloc(E2C);
    uint8_t *ct1 = (uint8_t *)malloc(E2C);
    uint8_t *dec = (uint8_t *)malloc(E2C * 2);
    uint8_t tag0[16], tag1[16];

    printf("== test_stream_exact_two_chunks ==\n");

    CFX_ASSERT(pt && ct0 && ct1 && dec);

    cfx_srand(0x1A2B3C4D);
    fuzz_fill(key, 32);
    fuzz_fill(nonce, 24);
    fuzz_fill(pt, E2C * 2);

    int rc = cfx_stream_xchacha20_poly1305_encrypt_chunk(
        ct0, tag0, pt, E2C, 0, 0, key, nonce);
    CFX_ASSERT(rc == 0);

    rc = cfx_stream_xchacha20_poly1305_encrypt_chunk(
        ct1, tag1, pt + E2C, E2C, 1, 1, key, nonce);
    CFX_ASSERT(rc == 0);

    rc = cfx_stream_xchacha20_poly1305_decrypt_chunk(
        dec, ct0, E2C, tag0, 0, 0, key, nonce);
    CFX_ASSERT(rc == 0);

    rc = cfx_stream_xchacha20_poly1305_decrypt_chunk(
        dec + E2C, ct1, E2C, tag1, 1, 1, key, nonce);
    CFX_ASSERT(rc == 0);

    CFX_ASSERT(memcmp(dec, pt, E2C * 2) == 0);

    free(pt); free(ct0); free(ct1); free(dec);
    #undef E2C
}
#endif /* CFX_BAREMETAL */

static void test_stream_one_byte(void) {
    /* single chunk with exactly 1 byte of plaintext */
    uint8_t key[32], nonce[24];
    uint8_t pt = 0x42, ct, dec, tag[16];

    printf("== test_stream_one_byte ==\n");

    cfx_srand(0xBEEF0001);
    fuzz_fill(key, 32);
    fuzz_fill(nonce, 24);

    int rc = cfx_stream_xchacha20_poly1305_encrypt_chunk(
        &ct, tag, &pt, 1, 0, 1, key, nonce);
    CFX_ASSERT(rc == 0);

    rc = cfx_stream_xchacha20_poly1305_decrypt_chunk(
        &dec, &ct, 1, tag, 0, 1, key, nonce);
    CFX_ASSERT(rc == 0);
    CFX_ASSERT(dec == pt);
}

#ifndef CFX_BAREMETAL  /* 2 × 64KB malloc — skip on bare-metal */
static void test_stream_boundary_plus_one(void) {
    /* 64KB + 1 byte: splits into full chunk + 1-byte final chunk */
    #define BP1 CFX_STREAM_CHUNK_SIZE
    uint8_t key[32], nonce[24];
    uint8_t *pt  = (uint8_t *)malloc(BP1 + 1);
    uint8_t *ct0 = (uint8_t *)malloc(BP1);
    uint8_t ct1;
    uint8_t *dec = (uint8_t *)malloc(BP1 + 1);
    uint8_t tag0[16], tag1[16];

    printf("== test_stream_boundary_plus_one ==\n");

    CFX_ASSERT(pt && ct0 && dec);

    cfx_srand(0xBEEF0002);
    fuzz_fill(key, 32);
    fuzz_fill(nonce, 24);
    fuzz_fill(pt, BP1 + 1);

    int rc = cfx_stream_xchacha20_poly1305_encrypt_chunk(
        ct0, tag0, pt, BP1, 0, 0, key, nonce);
    CFX_ASSERT(rc == 0);

    rc = cfx_stream_xchacha20_poly1305_encrypt_chunk(
        &ct1, tag1, pt + BP1, 1, 1, 1, key, nonce);
    CFX_ASSERT(rc == 0);

    rc = cfx_stream_xchacha20_poly1305_decrypt_chunk(
        dec, ct0, BP1, tag0, 0, 0, key, nonce);
    CFX_ASSERT(rc == 0);

    rc = cfx_stream_xchacha20_poly1305_decrypt_chunk(
        dec + BP1, &ct1, 1, tag1, 1, 1, key, nonce);
    CFX_ASSERT(rc == 0);

    CFX_ASSERT(memcmp(dec, pt, BP1 + 1) == 0);

    free(pt); free(ct0); free(dec);
    #undef BP1
}
#endif /* CFX_BAREMETAL */

static void test_stream_wrong_key(void) {
    /* decrypting with a different key must fail */
    uint8_t key[32], bad_key[32], nonce[24];
    uint8_t pt[48], ct[48], dec[48], tag[16];

    printf("== test_stream_wrong_key ==\n");

    cfx_srand(0xBEEF0003);
    fuzz_fill(key, 32);
    fuzz_fill(nonce, 24);
    fuzz_fill(pt, 48);
    memcpy(bad_key, key, 32);
    bad_key[0] ^= 0x01;

    cfx_stream_xchacha20_poly1305_encrypt_chunk(
        ct, tag, pt, 48, 0, 1, key, nonce);

    int rc = cfx_stream_xchacha20_poly1305_decrypt_chunk(
        dec, ct, 48, tag, 0, 1, bad_key, nonce);
    CFX_ASSERT(rc != 0);
}

static void test_stream_extension(void) {
    /* after a final chunk, an attacker-appended chunk must not authenticate */
    uint8_t key[32], nonce[24];
    uint8_t pt0[32], pt1[16];
    uint8_t ct0[32], ct1[16], tag0[16], tag1[16];
    uint8_t dec[32];

    printf("== test_stream_extension ==\n");

    cfx_srand(0xBEEF0004);
    fuzz_fill(key, 32);
    fuzz_fill(nonce, 24);
    fuzz_fill(pt0, 32);
    fuzz_fill(pt1, 16);

    /* encrypt a single-chunk stream (chunk 0, final) */
    cfx_stream_xchacha20_poly1305_encrypt_chunk(
        ct0, tag0, pt0, 32, 0, 1, key, nonce);

    /* encrypt chunk 1 as non-final (attacker tries to extend) */
    cfx_stream_xchacha20_poly1305_encrypt_chunk(
        ct1, tag1, pt1, 16, 1, 0, key, nonce);

    /* chunk 0 as final still works */
    int rc = cfx_stream_xchacha20_poly1305_decrypt_chunk(
        dec, ct0, 32, tag0, 0, 1, key, nonce);
    CFX_ASSERT(rc == 0);

    /* but trying to decrypt chunk 1 as non-final (counter=1, is_final=0)
     * succeeds cryptographically — the real protection is that the receiver
     * already saw is_final=1 on chunk 0 and stops.  Verify that if the
     * attacker re-labels chunk 0 as non-final, it fails: */
    rc = cfx_stream_xchacha20_poly1305_decrypt_chunk(
        dec, ct0, 32, tag0, 0, 0, key, nonce);
    CFX_ASSERT(rc != 0);
}

int main(void) {
    CFX_TEST(test_rfc8439_encrypt);
    CFX_TEST(test_rfc8439_decrypt);
    CFX_TEST(test_rfc8439_bad_tag);
    CFX_TEST(test_aead_fuzz_basic);
    CFX_TEST(test_xchacha_encrypt);
    CFX_TEST(test_xchacha_decrypt);
    CFX_TEST(test_xchacha_bad_tag);
    CFX_TEST(test_xchacha_fuzz);
    CFX_TEST(test_stream_single_chunk);
#ifndef CFX_BAREMETAL
    CFX_TEST(test_stream_multi_chunk);
#endif
    CFX_TEST(test_stream_empty);
    CFX_TEST(test_stream_wrong_order);
    CFX_TEST(test_stream_truncation);
    CFX_TEST(test_stream_tamper);
#ifndef CFX_BAREMETAL
    CFX_TEST(test_stream_exact_chunk);
    CFX_TEST(test_stream_exact_two_chunks);
#endif
    CFX_TEST(test_stream_one_byte);
#ifndef CFX_BAREMETAL
    CFX_TEST(test_stream_boundary_plus_one);
#endif
    CFX_TEST(test_stream_wrong_key);
    CFX_TEST(test_stream_extension);
    puts("ok");
    return 0;
}
