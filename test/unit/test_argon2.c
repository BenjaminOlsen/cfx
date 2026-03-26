/*
 * Argon2 tests — KAT (known answer tests) and parameter validation.
 *
 * Note: RFC 9106 test vectors use secret + associated data, which
 * this implementation omits (hardcoded to zero).  KAT values below
 * are verified against the reference C implementation (phc-winner-argon2)
 * with empty secret and AD.
 */

#include "cfx/argon2.h"
#include "cfx/macros.h"
#include <stdio.h>
#include <string.h>

/* helpers */

static int hex_to_bytes(const char *hex, uint8_t *out, size_t outlen) {
    for (size_t i = 0; i < outlen; i++) {
        unsigned int byte;
        if (sscanf(hex + 2 * i, "%2x", &byte) != 1) return -1;
        out[i] = (uint8_t)byte;
    }
    return 0;
}

/* KAT: argon2id "password" / "somesalt" / m=65536 t=3 p=4 */

static void test_argon2id_basic(void) {
    const char *password = "password";
    const char *salt_str = "somesalt";
    const uint8_t expected[] = {
        0x4f, 0xe1, 0xf9, 0xd5, 0x46, 0x2f, 0x63, 0xd4,
        0x3e, 0xb7, 0x98, 0xdc, 0x5e, 0xa0, 0x17, 0x14,
        0x51, 0xc4, 0xc3, 0x12, 0x97, 0xc7, 0x29, 0xcd,
        0x6d, 0x6a, 0x41, 0x4a, 0xcb, 0xa9, 0x14, 0x8a
    };

    uint8_t hash[32];
    int rc = cfx_argon2id(
        hash, sizeof(hash),
        (const uint8_t *)password, strlen(password),
        (const uint8_t *)salt_str, strlen(salt_str),
        65536, 3, 4);

    CFX_ASSERT(rc == 0);
    CFX_ASSERT(memcmp(hash, expected, 32) == 0);
}

/* KAT: argon2d  0x01*32 / 0x02*16 / m=32 t=3 p=4 */

static void test_argon2d_kat(void) {
    uint8_t password[32], salt[16];
    memset(password, 0x01, sizeof(password));
    memset(salt, 0x02, sizeof(salt));

    uint8_t expected[32];
    hex_to_bytes("9e34c31a47866ce0c30a90c69dd21022"
                 "d5329a3b75f9c513722dd2541fe93a1a",
                 expected, 32);

    uint8_t hash[32];
    int rc = cfx_argon2d(
        hash, sizeof(hash),
        password, sizeof(password),
        salt, sizeof(salt),
        32, 3, 4);

    CFX_ASSERT(rc == 0);
    CFX_ASSERT(memcmp(hash, expected, 32) == 0);
}

/* KAT: argon2i  0x01*32 / 0x02*16 / m=32 t=3 p=4 */

static void test_argon2i_kat(void) {
    uint8_t password[32], salt[16];
    memset(password, 0x01, sizeof(password));
    memset(salt, 0x02, sizeof(salt));

    uint8_t expected[32];
    hex_to_bytes("a9a7510e6db4d588ba3414cd0e094d48"
                 "0d683f97b9ccb612a544fe8ef65ba8e0",
                 expected, 32);

    uint8_t hash[32];
    int rc = cfx_argon2i(
        hash, sizeof(hash),
        password, sizeof(password),
        salt, sizeof(salt),
        32, 3, 4);

    CFX_ASSERT(rc == 0);
    CFX_ASSERT(memcmp(hash, expected, 32) == 0);
}

/* KAT: argon2id 0x01*32 / 0x02*16 / m=32 t=3 p=4 */

static void test_argon2id_kat(void) {
    uint8_t password[32], salt[16];
    memset(password, 0x01, sizeof(password));
    memset(salt, 0x02, sizeof(salt));

    uint8_t expected[32];
    hex_to_bytes("03aab965c12001c9d7d0d2de33192c04"
                 "94b684bb148196d73c1df1acaf6d0c2e",
                 expected, 32);

    uint8_t hash[32];
    int rc = cfx_argon2id(
        hash, sizeof(hash),
        password, sizeof(password),
        salt, sizeof(salt),
        32, 3, 4);

    CFX_ASSERT(rc == 0);
    CFX_ASSERT(memcmp(hash, expected, 32) == 0);
}

/* encode + verify round-trip */

static void test_argon2_encode_verify(void) {
    const char *password = "correcthorsebatterystaple";
    uint8_t salt[16] = {0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08,
                        0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f, 0x10};
    const uint8_t expected[] = {
        0x42, 0x50, 0xc6, 0x03, 0x54, 0x7a, 0x56, 0x21,
        0xbf, 0x40, 0xe7, 0x32, 0x15, 0xc4, 0x29, 0x41,
        0x7d, 0x64, 0x96, 0x5f, 0x82, 0x23, 0x99, 0x6c,
        0x72, 0xca, 0x4f, 0x95, 0x81, 0x28, 0x5a, 0x88
    };

    uint8_t hash[32];
    int rc = cfx_argon2id(
        hash, sizeof(hash),
        (const uint8_t *)password, strlen(password),
        salt, sizeof(salt),
        1024, 2, 1);
    CFX_ASSERT(rc == 0);
    CFX_ASSERT(memcmp(hash, expected, 32) == 0);

    char encoded[256];
    int len = cfx_argon2_encode(encoded, sizeof(encoded),
        hash, sizeof(hash),
        salt, sizeof(salt),
        1024, 2, 1,
        CFX_ARGON2ID);
    CFX_ASSERT(len > 0);

    /* verify correct password */
    rc = cfx_argon2_verify(encoded, (const uint8_t *)password, strlen(password));
    CFX_ASSERT(rc == 0);

    /* verify wrong password is rejected */
    const char *wrong = "wrongpassword";
    rc = cfx_argon2_verify(encoded, (const uint8_t *)wrong, strlen(wrong));
    CFX_ASSERT(rc != 0);
}

/* determinism: same inputs → same hash */

static void test_argon2id_deterministic(void) {
    const char *pwd = "determinism";
    uint8_t salt[16] = {0xde, 0xad, 0xbe, 0xef, 0xca, 0xfe, 0xba, 0xbe,
                        0x01, 0x23, 0x45, 0x67, 0x89, 0xab, 0xcd, 0xef};
    uint8_t h1[32], h2[32];

    int rc1 = cfx_argon2id(h1, 32, (const uint8_t *)pwd, strlen(pwd),
                           salt, 16, 32, 1, 1);
    int rc2 = cfx_argon2id(h2, 32, (const uint8_t *)pwd, strlen(pwd),
                           salt, 16, 32, 1, 1);

    CFX_ASSERT(rc1 == 0);
    CFX_ASSERT(rc2 == 0);
    CFX_ASSERT(memcmp(h1, h2, 32) == 0);
}

/* different passwords → different hashes */

static void test_argon2id_different_passwords(void) {
    uint8_t salt[16] = {0};
    uint8_t h1[32], h2[32];

    const char *pwd1 = "alpha";
    const char *pwd2 = "bravo";

    cfx_argon2id(h1, 32, (const uint8_t *)pwd1, strlen(pwd1), salt, 16, 32, 1, 1);
    cfx_argon2id(h2, 32, (const uint8_t *)pwd2, strlen(pwd2), salt, 16, 32, 1, 1);

    CFX_ASSERT(memcmp(h1, h2, 32) != 0);
}

/* different salts → different hashes */

static void test_argon2id_different_salts(void) {
    const char *pwd = "same_password";
    uint8_t salt1[16] = {0}, salt2[16] = {0};
    salt2[0] = 0x01;  /* one bit difference */
    uint8_t h1[32], h2[32];

    cfx_argon2id(h1, 32, (const uint8_t *)pwd, strlen(pwd), salt1, 16, 32, 1, 1);
    cfx_argon2id(h2, 32, (const uint8_t *)pwd, strlen(pwd), salt2, 16, 32, 1, 1);

    CFX_ASSERT(memcmp(h1, h2, 32) != 0);
}

/* empty password is valid */

static void test_argon2id_empty_password(void) {
    uint8_t salt[16] = {0x42};
    uint8_t hash[32];

    int rc = cfx_argon2id(hash, 32, (const uint8_t *)"", 0, salt, 16, 32, 1, 1);
    CFX_ASSERT(rc == 0);
}

/* minimum parameters (m=8, t=1, p=1) */

static void test_argon2id_min_params(void) {
    const char *pwd = "test";
    uint8_t salt[16] = {0};
    uint8_t hash[32];

    int rc = cfx_argon2id(hash, 32, (const uint8_t *)pwd, strlen(pwd),
                          salt, 16, 8, 1, 1);
    CFX_ASSERT(rc == 0);
}

/* parameter validation */

static void test_argon2_short_salt_rejected(void) {
    uint8_t hash[32];
    uint8_t pwd[] = "password";
    uint8_t salt[16] = {0};

    int rc = cfx_argon2id(hash, 32, pwd, 8, salt, 4, 32, 1, 1);
    CFX_ASSERT(rc == -1);
}

static void test_argon2_short_output_rejected(void) {
    uint8_t hash[32];
    uint8_t pwd[] = "password";
    uint8_t salt[16] = {0};

    int rc = cfx_argon2id(hash, 2, pwd, 8, salt, 16, 32, 1, 1);
    CFX_ASSERT(rc == -1);
}

static void test_argon2_zero_iterations_rejected(void) {
    uint8_t hash[32];
    uint8_t pwd[] = "password";
    uint8_t salt[16] = {0};

    int rc = cfx_argon2id(hash, 32, pwd, 8, salt, 16, 32, 0, 1);
    CFX_ASSERT(rc == -1);
}

int main(void) {
    CFX_TEST(test_argon2id_basic);
    CFX_TEST(test_argon2d_kat);
    CFX_TEST(test_argon2i_kat);
    CFX_TEST(test_argon2id_kat);
    CFX_TEST(test_argon2_encode_verify);

    CFX_TEST(test_argon2id_deterministic);
    CFX_TEST(test_argon2id_different_passwords);
    CFX_TEST(test_argon2id_different_salts);
    CFX_TEST(test_argon2id_empty_password);
    CFX_TEST(test_argon2id_min_params);

    CFX_TEST(test_argon2_short_salt_rejected);
    CFX_TEST(test_argon2_short_output_rejected);
    CFX_TEST(test_argon2_zero_iterations_rejected);

    puts("ok");
    return 0;
}
