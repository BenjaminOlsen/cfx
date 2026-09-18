#include "bge.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

    cfx_bge_free(decrypted, decrypted_len);
    cfx_bge_free(encrypted, encrypted_len);

    static const size_t lengths[] = {
        0, 1, 65535, 65536, 65537, 131072, 131073
    };
    for (size_t i = 0; i < sizeof(lengths) / sizeof(lengths[0]); ++i) {
        if (test_stream_roundtrip(lengths[i]) != 0)
            return 5;
    }
    return 0;
}
