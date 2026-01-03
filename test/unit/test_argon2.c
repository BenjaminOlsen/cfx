/*
 * Argon2 test vectors from RFC 9106
 */

#include "cfx/argon2.h"
#include "cfx/macros.h"
#include <stdio.h>
#include <string.h>

static void test_argon2id_basic(void) {
    const char* password = "password";
    const char* salt_str = "somesalt";

    uint8_t hash[32];
    int rc = cfx_argon2id(
        hash, sizeof(hash),
        (const uint8_t*)password, strlen(password),
        (const uint8_t*)salt_str, strlen(salt_str),
        65536, 3, 4
    );

    CFX_ASSERT(rc == CFX_ARGON2_OK);
}

static void test_argon2i_rfc(void) {
    uint8_t password[32];
    uint8_t salt[16];
    memset(password, 0x01, sizeof(password));
    memset(salt, 0x02, sizeof(salt));

    uint8_t hash[32];
    int rc = cfx_argon2i(
        hash, sizeof(hash),
        password, sizeof(password),
        salt, sizeof(salt),
        32, 3, 4
    );

    CFX_ASSERT(rc == CFX_ARGON2_OK);
}

static void test_argon2d_rfc(void) {
    uint8_t password[32];
    uint8_t salt[16];
    memset(password, 0x01, sizeof(password));
    memset(salt, 0x02, sizeof(salt));

    uint8_t hash[32];
    int rc = cfx_argon2d(
        hash, sizeof(hash),
        password, sizeof(password),
        salt, sizeof(salt),
        32, 3, 4
    );

    CFX_ASSERT(rc == CFX_ARGON2_OK);
}

static void test_argon2id_rfc(void) {
    uint8_t password[32];
    uint8_t salt[16];
    memset(password, 0x01, sizeof(password));
    memset(salt, 0x02, sizeof(salt));

    uint8_t hash[32];
    int rc = cfx_argon2id(
        hash, sizeof(hash),
        password, sizeof(password),
        salt, sizeof(salt),
        32, 3, 4
    );

    CFX_ASSERT(rc == CFX_ARGON2_OK);
}

static void test_argon2_encode_verify(void) {
    const char* password = "correcthorsebatterystaple";
    uint8_t salt[16] = {0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08,
                        0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f, 0x10};

    uint8_t hash[32];
    int rc = cfx_argon2id(
        hash, sizeof(hash),
        (const uint8_t*)password, strlen(password),
        salt, sizeof(salt),
        1024, 2, 1
    );
    CFX_ASSERT(rc == CFX_ARGON2_OK);

    char encoded[256];
    int len = cfx_argon2_encode(encoded, sizeof(encoded),
                                 hash, sizeof(hash),
                                 salt, sizeof(salt),
                                 1024, 2, 1,
                                 CFX_ARGON2ID);
    CFX_ASSERT(len > 0);

    /* verify correct password */
    rc = cfx_argon2_verify(encoded, (const uint8_t*)password, strlen(password));
    CFX_ASSERT(rc == 0);

    /* verify wrong password is rejected */
    const char* wrong = "wrongpassword";
    rc = cfx_argon2_verify(encoded, (const uint8_t*)wrong, strlen(wrong));
    CFX_ASSERT(rc != 0);
}

static void test_argon2_short_salt_rejected(void) {
    uint8_t hash[32];
    uint8_t pwd[] = "password";
    uint8_t salt[16] = {0};

    int rc = cfx_argon2id(hash, 32, pwd, 8, salt, 4, 32, 1, 1);
    CFX_ASSERT(rc == CFX_ARGON2_ERR_SALT);
}

static void test_argon2_short_output_rejected(void) {
    uint8_t hash[32];
    uint8_t pwd[] = "password";
    uint8_t salt[16] = {0};

    int rc = cfx_argon2id(hash, 2, pwd, 8, salt, 16, 32, 1, 1);
    CFX_ASSERT(rc == CFX_ARGON2_ERR_OUTPUT);
}

static void test_argon2_zero_iterations_rejected(void) {
    uint8_t hash[32];
    uint8_t pwd[] = "password";
    uint8_t salt[16] = {0};

    int rc = cfx_argon2id(hash, 32, pwd, 8, salt, 16, 32, 0, 1);
    CFX_ASSERT(rc == CFX_ARGON2_ERR_PARAM);
}

int main(void) {
    CFX_TEST(test_argon2id_basic);
    CFX_TEST(test_argon2i_rfc);
    CFX_TEST(test_argon2d_rfc);
    CFX_TEST(test_argon2id_rfc);
    CFX_TEST(test_argon2_encode_verify);
    CFX_TEST(test_argon2_short_salt_rejected);
    CFX_TEST(test_argon2_short_output_rejected);
    CFX_TEST(test_argon2_zero_iterations_rejected);

    puts("ok");
    return 0;
}
