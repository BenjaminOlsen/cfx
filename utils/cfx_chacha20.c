#include "cfx/chacha20.h"
#include "cfx/base64.h"
#include "misc.h"
#include "cfx_utils.h"

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

enum output_format {
    FMT_HEX,
    FMT_BASE64
};

static void usage(const char* name) {
    fprintf(stderr,
        "Usage: %s [options] <text>\n"
        "  Encrypt text with ChaCha20 (RFC8439 layout).\n\n"
        "Options:\n"
        "  -k <key>    Key as ASCII (padded/truncated to 32 bytes) OR 0x64_hex_chars\n"
        "  -c <ctr>    Initial 32-bit counter (default 0)\n"
        "  -n <nonce>  Nonce as 0x24_hex_chars (96-bit, default zeros)\n"
        "  -x          Output as hex (default)\n"
        "  -b64        Output as base64\n"
        "  -v          Verbose output\n"
        "  -h, --help  Show this help\n\n"
        "Example:\n"
        "  %s -k 0xe999fb0f95c3b2e2a703ea2d55565a8a2a8a725291c4ad4d614c20c31a14708e \\\n"
        "     -n 0x0ccd47f7ddf772db86352163 -c 2 HELLOOOO\n",
        name, name);
}

int cfx_chacha20_run(int argc, char** argv) {
    if (argc < 2 || (argc == 2 && (strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0))) {
        usage(argv[0]);
        return argc < 2 ? EXIT_FAILURE : 0;
    }

    const char* key_in   = NULL;
    const char* nonce_in = NULL;
    const char* pt       = NULL;
    uint32_t counter     = 0;
    int verbose = 0;
    enum output_format fmt = FMT_HEX;

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
        } else if (strcmp(argv[i], "-x") == 0) {
            fmt = FMT_HEX;
        } else if (strcmp(argv[i], "-b64") == 0) {
            fmt = FMT_BASE64;
        } else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            usage(argv[0]);
            return 0;
        } else if (argv[i][0] == '-') {
            fprintf(stderr, "Unknown option: %s\n", argv[i]);
            usage(argv[0]);
            return EXIT_FAILURE;
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
        if (cfx_parse_hex(key_in + 2, key, sizeof(key)) != 0) {
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
            cfx_parse_hex(nonce_in + 2, nonce, sizeof(nonce)) != 0) {
            fprintf(stderr, "error: -n expects hex: with exactly 24 hex chars, using default\n");
        }
    }

    size_t len = strlen(pt);
    uint8_t* ct = (uint8_t*)malloc(len ? len : 1);
    if (!ct) { perror("malloc"); return EXIT_FAILURE; }

    cfx_chacha20_encrypt(key, counter, nonce, (const uint8_t*)pt, len, ct);

    if (verbose) printf("ciphertext (%zu bytes):\n", len);

    if (fmt == FMT_BASE64) {
        size_t b64_len = cfx_base64_enc_len(len);
        char* b64 = (char*)malloc(b64_len + 1);
        if (b64) {
            cfx_base64_encode(b64, &b64_len, ct, len);
            b64[b64_len] = '\0';
            printf("%s\n", b64);
            free(b64);
        }
    } else {
        for (size_t i = 0; i < len; ++i) {
            printf("%02x", ct[i]);
        }
        printf("\n");
    }

    free(ct);
    return 0;
}
