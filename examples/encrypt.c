#include "cfx/chacha20.h"
#include "utils.h"

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

static void usage(const char* name) {
    fprintf(stderr,
        "Usage: %s -k <key|0x...> [-c counter] [-n 0x...] <text>\n"
        "  Encrypts <text> with ChaCha20 (RFC8439 layout: 32B key, 32-bit counter, 96-bit nonce).\n"
        "  -k  key as ASCII (padded/truncated to 32 bytes) OR 0x64_hex_chars\n"
        "  -c  initial 32-bit counter (default 0)\n"
        "  -n  nonce as 0x24_hex_chars (96-bit, default all zeros)\n"
        "  -v  verbose printing\n"
        "  text: plaintext string to encrypt\n"
        "  -------\n"
        "  example: %s -k 0xe999fb0f95c3b2e2a703ea2d55565a8a2a8a725291c4ad4d614c20c31a14708e\n"
        "              -n 0x0ccd47f7ddf772db86352163 -c 2 HELLOOOO\n"
        "  output: 8e 9a 97 a5 1e 21 fd d8\n",
        name, name);
}

int main(int argc, char** argv) {
    if (argc < 2) {
        usage(argv[0]);
        return EXIT_FAILURE;
    }

    const char* key_in   = NULL;
    const char* nonce_in = NULL;
    const char* pt       = NULL;
    uint32_t counter     = 0;
    int verbose = 0;

    #define CHECK_ARG(i) if (i >= argc) { usage(argv[0]); return EXIT_FAILURE; }

    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-k") == 0) {
            ++i;
            CHECK_ARG(i);
            key_in = argv[i];
        } else if (strcmp(argv[i], "-c") == 0) {
            ++i;
            CHECK_ARG(i);
            errno = 0;
            unsigned long v = strtoul(argv[i], NULL, 0);
            if (errno) { perror("counter"); return EXIT_FAILURE; }
            counter = (uint32_t)v;
        } else if (strcmp(argv[i], "-n") == 0) {
            ++i;
            CHECK_ARG(i);
            nonce_in = argv[i];
        } else if (strcmp(argv[i], "-v") == 0) {
            verbose = 1;
        } else if (argv[i][0] == '-' ) {
            /* unknown flag */
            usage(argv[0]); return EXIT_FAILURE;
        } else {
            pt = argv[i];
        }
    }

    if (!pt) {
        fprintf(stderr, "error: missing plaintext\n");
        usage(argv[0]);
        return EXIT_FAILURE;
    }
    const char default_key[] = "1234567890";
    if (!key_in) {
        if (verbose) fprintf(stderr, "using default key\n");
        key_in = default_key;
    }

    uint8_t key[32] = {0};
    if (strncmp(key_in, "0x", 2) == 0) {
        if (parse_hex(key_in + 2, key, sizeof(key)) != 0) {
            fprintf(stderr, "error: -k expects hex: with exactly 64 hex chars\n");
            return EXIT_FAILURE;
        }
    } else {
        size_t kl = strlen(key_in);
        if (kl > sizeof(key)) kl = sizeof(key);
        memcpy(key, key_in, kl);
        if (strlen(key_in) > sizeof(key)) {
            fprintf(stderr, "warn: ASCII key longer than 32 bytes; truncated\n");
        }
    }

    if (verbose) {
        fprintf(stderr, "key: ");
        for (size_t i = 0; i < sizeof key; ++i) {
            fprintf(stderr, "%02x ", key[i]);
        }
        fprintf(stderr, "\n");
    }

    uint8_t nonce[12] = {0};
    if (nonce_in) {
        if (strncmp(nonce_in, "0x", 2) != 0 ||
            parse_hex(nonce_in + 2, nonce, sizeof(nonce)) != 0) {
            fprintf(stderr, "error: -n expects hex: with exactly 24 hex chars, using default\n");
        }
    }

    size_t len = strlen(pt);
    uint8_t* ct = (uint8_t*)malloc(len ? len : 1);
    if (!ct) { perror("malloc"); return EXIT_FAILURE; }

    cfx_chacha20_encrypt(key, counter, nonce, (const uint8_t*)pt, len, ct);

    if (verbose) printf("ciphertext (%zu bytes):\n", len);
    for (size_t i = 0; i < len; ++i) printf("%02x%s", ct[i], (i+1==len) ? "" : " ");
    printf("\n");

    free(ct);
    return 0;
}
