#include "bge.h"

#include <stdint.h>
#include <string.h>

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
    return 0;
}
