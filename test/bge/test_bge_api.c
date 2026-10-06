#include "bge.h"
#include "cfx/base64.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int test_armored_whitespace(const uint8_t *encrypted, size_t encrypted_len,
                                  const uint8_t *expected, size_t expected_len,
                                  const uint8_t *password, size_t password_len) {
    static const char *prefixes[] = {
        "", " \t\r\n\f\v", "\n\n", "   "
    };
    static const char header[] = "-----BEGIN BGE MESSAGE-----\r\n";
    static const char footer[] = "\r\n-----END BGE MESSAGE-----";
    static const char suffix[] = " \t\r\n\f\v";
    size_t encoded_len = cfx_base64_enc_len(encrypted_len);

    for (size_t i = 0; i < sizeof(prefixes) / sizeof(prefixes[0]); ++i) {
        size_t prefix_len = strlen(prefixes[i]);
        size_t total = prefix_len + sizeof(header) - 1 + encoded_len +
                       sizeof(footer) - 1 + sizeof(suffix) - 1;
        char *armored = malloc(total);
        if (!armored) return 1;
        char *body = armored + prefix_len + sizeof(header) - 1;
        memcpy(armored, prefixes[i], prefix_len);
        memcpy(armored + prefix_len, header, sizeof(header) - 1);
        size_t written = encoded_len;
        if (cfx_base64_encode(body, &written, encrypted, encrypted_len) != 0) {
            free(armored);
            return 1;
        }
        memcpy(body + encoded_len, footer, sizeof(footer) - 1);
        memcpy(body + encoded_len + sizeof(footer) - 1, suffix, sizeof(suffix) - 1);
        uint8_t *decoded = NULL;
        size_t decoded_len = 0;
        int rc = cfx_bge_decrypt((const uint8_t *)armored, total,
                               password, password_len, &decoded, &decoded_len);
        free(armored);
        int failed = rc != 0 || decoded_len != expected_len ||
                     memcmp(decoded, expected, expected_len) != 0;
        cfx_bge_free(decoded, decoded_len);
        if (failed) {
            fprintf(stderr, "armored whitespace test %zu failed\n", i);
            return 1;
        }
    }

    static const uint8_t whitespace[] = " \t\r\n\f\v";
    uint8_t *decoded = NULL;
    size_t decoded_len = 0;
    int rc = cfx_bge_decrypt(whitespace, sizeof(whitespace) - 1,
                           password, password_len, &decoded, &decoded_len);
    cfx_bge_free(decoded, decoded_len);
    return rc == -2 ? 0 : 1;
}

static int test_stream_roundtrip(size_t len) {
    static const uint8_t password[] = "test-password";
    FILE *input = NULL;
    FILE *encrypted = NULL;
    FILE *decrypted = NULL;
    uint8_t *expected = NULL;
    int result = 1;

#define CHECK(expr) do {                                              \
    if (!(expr)) {                                                    \
        fprintf(stderr, "stream test (%zu bytes), line %d: %s\n",     \
                len, __LINE__, #expr);                               \
        goto cleanup;                                                \
    }                                                                \
} while (0)

    input = tmpfile();
    encrypted = tmpfile();
    decrypted = tmpfile();
    CHECK(input && encrypted && decrypted);

    expected = malloc(len ? len : 1);
    CHECK(expected != NULL);

    /* adjacent chunks have different contents. */
    uint32_t state = 1;
    for (size_t i = 0; i < len; ++i) {
        state = state * UINT32_C(1664525) + UINT32_C(1013904223);
        expected[i] = (uint8_t)(state >> 24);
    }

    CHECK(fwrite(expected, 1, len, input) == len);
    CHECK(fseek(input, 0, SEEK_SET) == 0);

    CHECK(cfx_bge_encrypt_stream(
        input, encrypted, password, sizeof(password) - 1) == 0);
    CHECK(fseek(encrypted, 0, SEEK_SET) == 0);

    CHECK(cfx_bge_decrypt_stream(
        encrypted, decrypted, password, sizeof(password) - 1) == 0);
    CHECK(fseek(decrypted, 0, SEEK_SET) == 0);

    for (size_t i = 0; i < len; ++i) {
        int byte = fgetc(decrypted);
        CHECK(byte != EOF);
        CHECK(byte == expected[i]);
    }

    /* Detect extra output as well as missing/wrong bytes. */
    CHECK(fgetc(decrypted) == EOF);
    CHECK(feof(decrypted) && !ferror(decrypted));
    result = 0;

cleanup:
    free(expected);
    if (decrypted) fclose(decrypted);
    if (encrypted) fclose(encrypted);
    if (input) fclose(input);
#undef CHECK
    return result;
}

int main(void) {
    static const uint8_t message[] = {'b', 'g', 'e', 0, 'a', 'p', 'i'};
    static const uint8_t password[] = "test-password";
    uint8_t *encrypted = NULL;
    size_t encrypted_len = 0;
    uint8_t *decrypted = NULL;
    size_t decrypted_len = 0;

    if (cfx_bge_encrypt(message, sizeof(message), password, sizeof(password) - 1,
                        &encrypted, &encrypted_len) != 0) {
        return 1;
    }
    if (!encrypted || encrypted_len <= sizeof(message)) {
        return 2;
    }

    if (cfx_bge_decrypt(encrypted, encrypted_len, password, sizeof(password) - 1,
                        &decrypted, &decrypted_len) != 0) {
        cfx_bge_free(encrypted, encrypted_len);
        return 3;
    }
    if (decrypted_len != sizeof(message) ||
        memcmp(decrypted, message, sizeof(message)) != 0) {
        cfx_bge_free(decrypted, decrypted_len);
        cfx_bge_free(encrypted, encrypted_len);
        return 4;
    }

    int armor_failed = test_armored_whitespace(
        encrypted, encrypted_len, message, sizeof(message),
        password, sizeof(password) - 1);
    cfx_bge_free(decrypted, decrypted_len);
    cfx_bge_free(encrypted, encrypted_len);
    if (armor_failed) return 6;

    static const size_t lengths[] = {
        0, 1, 65535, 65536, 65537, 131072, 131073
    };
    for (size_t i = 0; i < sizeof(lengths) / sizeof(lengths[0]); ++i) {
        if (test_stream_roundtrip(lengths[i]) != 0)
            return 5;
    }
    return 0;
}
