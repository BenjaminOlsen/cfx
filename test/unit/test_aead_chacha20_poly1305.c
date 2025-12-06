/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#include "cfx/aead_chacha20_poly1305.h" 
#include "cfx/macros.h"
#include "cfx/rand.h"

#include <string.h>
#include <stdint.h>
#include <stdio.h>


#define PR printf

static void print_hex(const char* label, const uint8_t* buf, size_t len) {
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


static void fuzz_fill(uint8_t* buf, size_t len) {
    for (size_t i = 0; i < len; ++i) {
        buf[i] = (uint8_t)cfx_rand();
    }
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

        const uint8_t* aad_ptr = aad;
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

int main(void) {
    CFX_TEST(test_rfc8439_encrypt);
    CFX_TEST(test_rfc8439_decrypt);
    CFX_TEST(test_rfc8439_bad_tag);
    CFX_TEST(test_aead_fuzz_basic);
    return 0;
}
