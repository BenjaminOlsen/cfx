#include "cfx/rand.h"
#include "cfx/poly1305.h"
#include "cfx/base64.h"
#include "misc.h"
#include "cfx_utils.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

enum output_format {
    FMT_HEX,
    FMT_BASE64
};

static void usage(const char *prog) {
    fprintf(stderr,
        "Usage: %s [options] <message>\n"
        "  Compute Poly1305 MAC of a message.\n\n"
        "Options:\n"
        "  -k <key>    Key as ASCII (padded/truncated to 32 bytes) OR 0x64_hex_chars\n"
        "  -s <seed>   Seed for random key generation (ignored if -k specified)\n"
        "  -x          Output as hex (default)\n"
        "  -b64        Output as base64\n"
        "  -h, --help  Show this help\n",
        prog);
}

int cfx_mac_run(int argc, char** argv) {
    const char* prog = argv[0];
    if (argc < 2 || (argc == 2 && (strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0))) {
        usage(prog);
        return argc < 2 ? 1 : 0;
    }

    const char* key_in  = NULL;
    const char* seed_in = NULL;
    const char* msg_in  = NULL;
    enum output_format fmt = FMT_HEX;

    #define CHECK_ARG(i) if (i >= argc) { usage(argv[0]); return EXIT_FAILURE; }
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-k") == 0) {
            ++i;
            CHECK_ARG(i);
            key_in = argv[i];
        } else if (strcmp(argv[i], "-s") == 0) {
            ++i;
            CHECK_ARG(i);
            seed_in = argv[i];
        } else if (strcmp(argv[i], "-x") == 0) {
            fmt = FMT_HEX;
        } else if (strcmp(argv[i], "-b64") == 0) {
            fmt = FMT_BASE64;
        } else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            usage(prog);
            return 0;
        } else if (argv[i][0] == '-') {
            fprintf(stderr, "Unknown option: %s\n", argv[i]);
            usage(prog);
            return EXIT_FAILURE;
        } else {
            msg_in = argv[i];
        }
    }

    if (!msg_in) {
        fprintf(stderr, "error: missing message\n");
        usage(prog);
        return EXIT_FAILURE;
    }

    const uint8_t* msg = (const uint8_t*)msg_in;
    size_t len = strlen(msg_in);
    uint8_t tag[16];
    uint8_t key[32];
    memset(key, 0, sizeof key);

    if (!key_in) {
        if (seed_in) {
            unsigned seed = strtoull(seed_in, NULL, 0);
            printf("got seed: %u\n", seed);
            cfx_srand(seed);
            cfx_rand_bytes((void*)key, sizeof(key));
        } else {
            cfx_rand_bytes_os((void*)key, sizeof(key));
        }
    } else if (strncmp(key_in, "0x", 2) == 0) {
        if (cfx_parse_hex(key_in + 2, key, sizeof(key)) != 0) {
            fprintf(stderr, "error: -k expects hex: with exactly 32 bytes\n");
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

    cfx_poly1305_mac(key, msg, len, tag);

    /* Print key */
    printf("key: ");
    for (size_t i = 0; i < sizeof(key); ++i) {
        printf("%02x", key[i]);
    }
    printf("\n\nMAC: ");

    /* Print MAC in selected format */
    if (fmt == FMT_BASE64) {
        size_t b64_len = cfx_base64_enc_len(sizeof(tag));
        char b64[32];
        cfx_base64_encode(b64, &b64_len, tag, sizeof(tag));
        b64[b64_len] = '\0';
        printf("%s\n", b64);
    } else {
        for (size_t i = 0; i < sizeof(tag); ++i) {
            printf("%02x", tag[i]);
        }
        printf("\n");
    }

    return 0;
}
