/*
 * Argon2 tests — KAT (known answer tests) and parameter validation.
 *
 * Note: RFC 9106 test vectors use secret + associated data, which
 * this implementation omits (hardcoded to zero).  The KAT values
 * below are generated from the cfx implementation with those fields
 * absent, so they differ from the RFC.  They still catch regressions.
 */

#include "cfx/argon2.h"
#include "cfx/macros.h"
#include <stdio.h>
#include <string.h>

/* -- helpers ------------------------------------------------------- */

static int hex_to_bytes(const char *hex, uint8_t *out, size_t outlen) {
    for (size_t i = 0; i < outlen; i++) {
        unsigned int byte;
        if (sscanf(hex + 2 * i, "%2x", &byte) != 1) return -1;
        out[i] = (uint8_t)byte;
    }
    return 0;
}

/* -- KAT: argon2id "password" / "somesalt" / m=65536 t=3 p=4 ---- */

static void test_argon2id_basic(void) {
    const char *password = "password";
    const char *salt_str = "somesalt";
    const uint8_t expected[] = {
        0xf4, 0x02, 0x70, 0xaa, 0x34, 0xfd, 0x8e, 0x8c,
        0x73, 0x4c, 0xef, 0xa8, 0x5f, 0xa3, 0xbf, 0x52,
        0xab, 0x7d, 0x4c, 0x3b, 0xad, 0x0f, 0x49, 0x9d,
        0x69, 0x39, 0x2b, 0xba, 0xfd, 0xca, 0xd6, 0x20
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

/* -- KAT: argon2d  0x01*32 / 0x02*16 / m=32 t=3 p=4 ------------ */

static void test_argon2d_kat(void) {
    uint8_t password[32], salt[16];
    memset(password, 0x01, sizeof(password));
    memset(salt, 0x02, sizeof(salt));

    uint8_t expected[32];
    hex_to_bytes("802327c66e61474206ff1e52319d9dd5"
                 "e3b2013661cccd63adf5acb3098242af",
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

/* -- KAT: argon2i  0x01*32 / 0x02*16 / m=32 t=3 p=4 ------------ */

static void test_argon2i_kat(void) {
    uint8_t password[32], salt[16];
    memset(password, 0x01, sizeof(password));
    memset(salt, 0x02, sizeof(salt));

    uint8_t expected[32];
    hex_to_bytes("f459c0a2321f4c468cda14594f614cda"
                 "c3018c987dfd0000d4196573d4c40462",
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

/* -- KAT: argon2id 0x01*32 / 0x02*16 / m=32 t=3 p=4 ------------ */

static void test_argon2id_kat(void) {
    uint8_t password[32], salt[16];
    memset(password, 0x01, sizeof(password));
    memset(salt, 0x02, sizeof(salt));

    uint8_t expected[32];
    hex_to_bytes("5876ffa323150be782e707cec93aece6"
                 "807dfb2b223f30495e98e5c30370738d",
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

/* -- encode + verify round-trip ------------------------------------ */

static void test_argon2_encode_verify(void) {
    const char *password = "correcthorsebatterystaple";
    uint8_t salt[16] = {0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08,
                        0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f, 0x10};
    const uint8_t expected[] = {
        0x52, 0xbb, 0xfa, 0xf3, 0xce, 0x18, 0xcc, 0x0f,
        0x37, 0xc9, 0x15, 0xe8, 0x7a, 0x8a, 0xad, 0x20,
        0xe1, 0xff, 0xb0, 0x3c, 0x15, 0x4b, 0x11, 0x7a,
        0x4e, 0x1d, 0xd7, 0x8d, 0x77, 0xa7, 0xa2, 0xc0
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

/* -- determinism: same inputs → same hash -------------------------- */

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

/* -- different passwords → different hashes ------------------------ */

static void test_argon2id_different_passwords(void) {
    uint8_t salt[16] = {0};
    uint8_t h1[32], h2[32];

    const char *pwd1 = "alpha";
    const char *pwd2 = "bravo";

    cfx_argon2id(h1, 32, (const uint8_t *)pwd1, strlen(pwd1), salt, 16, 32, 1, 1);
    cfx_argon2id(h2, 32, (const uint8_t *)pwd2, strlen(pwd2), salt, 16, 32, 1, 1);

    CFX_ASSERT(memcmp(h1, h2, 32) != 0);
}

/* -- different salts → different hashes ---------------------------- */

static void test_argon2id_different_salts(void) {
    const char *pwd = "same_password";
    uint8_t salt1[16] = {0}, salt2[16] = {0};
    salt2[0] = 0x01;  /* one bit difference */
    uint8_t h1[32], h2[32];

    cfx_argon2id(h1, 32, (const uint8_t *)pwd, strlen(pwd), salt1, 16, 32, 1, 1);
    cfx_argon2id(h2, 32, (const uint8_t *)pwd, strlen(pwd), salt2, 16, 32, 1, 1);

    CFX_ASSERT(memcmp(h1, h2, 32) != 0);
}

/* -- empty password is valid --------------------------------------- */

static void test_argon2id_empty_password(void) {
    uint8_t salt[16] = {0x42};
    uint8_t hash[32];

    int rc = cfx_argon2id(hash, 32, (const uint8_t *)"", 0, salt, 16, 32, 1, 1);
    CFX_ASSERT(rc == 0);
}

/* -- minimum parameters (m=8, t=1, p=1) --------------------------- */

static void test_argon2id_min_params(void) {
    const char *pwd = "test";
    uint8_t salt[16] = {0};
    uint8_t hash[32];

    int rc = cfx_argon2id(hash, 32, (const uint8_t *)pwd, strlen(pwd),
                          salt, 16, 8, 1, 1);
    CFX_ASSERT(rc == 0);
}

/* -- parameter validation ------------------------------------------ */

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
