/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#include "cfx/hmac.h"
#include "cfx/macros.h"

#include <stdint.h>
#include <stdio.h>
#include <string.h>

static int hex_nibble(char c) {
    if (c >= '0' && c <= '9') return c - '0';
    if (c >= 'a' && c <= 'f') return 10 + (c - 'a');
    if (c >= 'A' && c <= 'F') return 10 + (c - 'A');
    return -1;
}

static void hex_to_bytes(const char *hex, uint8_t *out, size_t out_len) {
    CFX_ASSERT(strlen(hex) == out_len * 2);
    for (size_t i = 0; i < out_len; i++) {
        int hi = hex_nibble(hex[2 * i]);
        int lo = hex_nibble(hex[2 * i + 1]);
        CFX_ASSERT(hi >= 0 && lo >= 0);
        out[i] = (uint8_t)((hi << 4) | lo);
    }
}

static void check_hmac256(const char *name,
                           const uint8_t *key, size_t key_len,
                           const uint8_t *data, size_t data_len,
                           const char *expected_hex)
{
    uint8_t got[32], expect[32];
    hex_to_bytes(expected_hex, expect, 32);

    cfx_hmac_sha256(got, key, key_len, data, data_len);
    if (memcmp(got, expect, 32) != 0) {
        fprintf(stderr, "[FAIL] HMAC-SHA256 %s\n", name);
        abort();
    }

    /* also verify streaming gives the same result */
    cfx_hmac_sha256_ctx ctx;
    uint8_t got2[32];
    cfx_hmac_sha256_init(&ctx, key, key_len);
    /* feed one byte at a time */
    for (size_t i = 0; i < data_len; i++) {
        cfx_hmac_sha256_update(&ctx, data + i, 1);
    }
    cfx_hmac_sha256_final(&ctx, got2);
    CFX_ASSERT(memcmp(got, got2, 32) == 0);

    printf("[OK] HMAC-SHA256 %s\n", name);
}

static void check_hmac512(const char *name,
                           const uint8_t *key, size_t key_len,
                           const uint8_t *data, size_t data_len,
                           const char *expected_hex)
{
    uint8_t got[64], expect[64];
    hex_to_bytes(expected_hex, expect, 64);

    cfx_hmac_sha512(got, key, key_len, data, data_len);
    if (memcmp(got, expect, 64) != 0) {
        fprintf(stderr, "[FAIL] HMAC-SHA512 %s\n", name);
        abort();
    }

    /* streaming */
    cfx_hmac_sha512_ctx ctx;
    uint8_t got2[64];
    cfx_hmac_sha512_init(&ctx, key, key_len);
    for (size_t i = 0; i < data_len; i++) {
        cfx_hmac_sha512_update(&ctx, data + i, 1);
    }
    cfx_hmac_sha512_final(&ctx, got2);
    CFX_ASSERT(memcmp(got, got2, 64) == 0);

    printf("[OK] HMAC-SHA512 %s\n", name);
}

/* RFC 4231 test vectors */

static void test_rfc4231_case1(void) {
    /* Key = 20 bytes of 0x0b */
    uint8_t key[20];
    memset(key, 0x0b, 20);
    const uint8_t *data = (const uint8_t *)"Hi There";

    check_hmac256("RFC4231 #1", key, 20, data, 8,
        "b0344c61d8db38535ca8afceaf0bf12b"
        "881dc200c9833da726e9376c2e32cff7");

    check_hmac512("RFC4231 #1", key, 20, data, 8,
        "87aa7cdea5ef619d4ff0b4241a1d6cb0"
        "2379f4e2ce4ec2787ad0b30545e17cde"
        "daa833b7d6b8a702038b274eaea3f4e4"
        "be9d914eeb61f1702e696c203a126854");
}

static void test_rfc4231_case2(void) {
    /* Key = "Jefe" */
    const uint8_t *key = (const uint8_t *)"Jefe";
    const uint8_t *data = (const uint8_t *)"what do ya want for nothing?";

    check_hmac256("RFC4231 #2", key, 4, data, 28,
        "5bdcc146bf60754e6a042426089575c7"
        "5a003f089d2739839dec58b964ec3843");

    check_hmac512("RFC4231 #2", key, 4, data, 28,
        "164b7a7bfcf819e2e395fbe73b56e0a3"
        "87bd64222e831fd610270cd7ea250554"
        "9758bf75c05a994a6d034f65f8f0e6fd"
        "caeab1a34d4a6b4b636e070a38bce737");
}

static void test_rfc4231_case3(void) {
    /* Key = 20 bytes of 0xaa, Data = 50 bytes of 0xdd */
    uint8_t key[20], data[50];
    memset(key, 0xaa, 20);
    memset(data, 0xdd, 50);

    check_hmac256("RFC4231 #3", key, 20, data, 50,
        "773ea91e36800e46854db8ebd09181a7"
        "2959098b3ef8c122d9635514ced565fe");

    check_hmac512("RFC4231 #3", key, 20, data, 50,
        "fa73b0089d56a284efb0f0756c890be9"
        "b1b5dbdd8ee81a3655f83e33b2279d39"
        "bf3e848279a722c806b485a47e67c807"
        "b946a337bee8942674278859e13292fb");
}

static void test_rfc4231_case4(void) {
    /* Key = 0x0102...19 (25 bytes), Data = 50 bytes of 0xcd */
    uint8_t key[25], data[50];
    for (int i = 0; i < 25; i++) key[i] = (uint8_t)(i + 1);
    memset(data, 0xcd, 50);

    check_hmac256("RFC4231 #4", key, 25, data, 50,
        "82558a389a443c0ea4cc819899f2083a"
        "85f0faa3e578f8077a2e3ff46729665b");

    check_hmac512("RFC4231 #4", key, 25, data, 50,
        "b0ba465637458c6990e5a8c5f61d4af7"
        "e576d97ff94b872de76f8050361ee3db"
        "a91ca5c11aa25eb4d679275cc5788063"
        "a5f19741120c4f2de2adebeb10a298dd");
}

static void test_rfc4231_case6(void) {
    /* Key = 131 bytes of 0xaa (key > block size — tests hashing the key) */
    uint8_t key[131];
    memset(key, 0xaa, 131);
    const uint8_t *data = (const uint8_t *)
        "Test Using Larger Than Block-Size Key - Hash Key First";

    check_hmac256("RFC4231 #6", key, 131, data, 54,
        "60e431591ee0b67f0d8a26aacbf5b77f"
        "8e0bc6213728c5140546040f0ee37f54");

    check_hmac512("RFC4231 #6", key, 131, data, 54,
        "80b24263c7c1a3ebb71493c1dd7be8b4"
        "9b46d1f41b4aeec1121b013783f8f352"
        "6b56d037e05f2598bd0fd2215d6a1e52"
        "95e64f73f63f0aec8b915a985d786598");
}

static void test_rfc4231_case7(void) {
    /* Key = 131 bytes of 0xaa, longer data */
    uint8_t key[131];
    memset(key, 0xaa, 131);
    const uint8_t *data = (const uint8_t *)
        "This is a test using a larger than block-size key and a "
        "larger than block-size data. The key needs to be hashed "
        "before being used by the HMAC algorithm.";

    check_hmac256("RFC4231 #7", key, 131, data, 152,
        "9b09ffa71b942fcb27635fbcd5b0e944"
        "bfdc63644f0713938a7f51535c3a35e2");

    check_hmac512("RFC4231 #7", key, 131, data, 152,
        "e37b6a775dc87dbaa4dfa9f96e5e3ffd"
        "debd71f8867289865df5a32d20cdc944"
        "b6022cac3c4982b10d5eeb55c3e4de15"
        "134676fb6de0446065c97440fa8c6a58");
}

int main(void) {
    CFX_TEST(test_rfc4231_case1);
    CFX_TEST(test_rfc4231_case2);
    CFX_TEST(test_rfc4231_case3);
    CFX_TEST(test_rfc4231_case4);
    CFX_TEST(test_rfc4231_case6);
    CFX_TEST(test_rfc4231_case7);
    printf("All HMAC tests passed.\n");
    return 0;
}
