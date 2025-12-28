#include "cfx/rand.h"
#include "cfx/poly1305.h"
#include "utils.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>


static void usage(const char *prog) {
    fprintf(stderr, "Usage: %s -k <key|0x...> -s <seed|0x...> <message>\n"
            "   outputs poly1305 signature of the message.\n"
            "   -k  key as ASCII (padded/truncated to 32 bytes) OR 0x64_hex_chars\n"
            "   -s seed for random key generation (unsigned int value, optionally hex). ignored if key specified\n",
            prog);
    exit(1);
}

const char* key_in  = NULL;
const char* seed_in = NULL;
const char* msg_in  = NULL;

int main(int argc, char** argv) {
    const char* prog = argv[0];
    if (argc < 2) {
        usage(prog);
        return EXIT_FAILURE;
    }

    #define CHECK_ARG(i) if (i >= argc) { usage(argv[0]); return EXIT_FAILURE; }
    for (int i = 0; i < argc; ++i) {
        if (strcmp(argv[i], "-k") == 0) {
            ++i;
            CHECK_ARG(i);
            key_in = argv[i];
        } else if (strcmp(argv[i], "-s") == 0) {
            ++i;
            CHECK_ARG(i);
            seed_in = argv[i];
        } else {
            msg_in = argv[i];
        }

    }


    const uint8_t* msg = (const uint8_t*)msg_in;
    size_t len = strlen(msg_in);
    uint8_t tag[16];
    uint8_t key[32];
    memset(key, 0, sizeof key);


    if (!key_in) {
        if (seed_in) {
            unsigned seed = strtoull(seed_in, NULL, 0);  /* base 0: deduces base from a */;
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
    printf("key: ");
    size_t i;
    for (i = 0; i < sizeof(key); ++i) {
        printf("%02x", key[i]);
    }
    printf("\n\n");

    printf("MAC:\n");
    for (i = 0; i < sizeof(tag); ++i) {
        printf("%02x ", tag[i]);
    }
    printf("\n");

    return 0;
}
