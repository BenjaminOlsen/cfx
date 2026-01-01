#include "cfx/rand.h"
#include "cfx/base64.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include "cfx_cmd.h"
#include "misc.h"

static void usage(const char *prog) {
    fprintf(stderr,
        "Usage: %s [options] <bytes>\n"
        "  Generate random bytes for use as keys.\n\n"
        "Options:\n"
        "  -s <seed>   Seed for RNG (default: random)\n"
        "  -x          Output as hex (default)\n"
        "  -b64        Output as base64\n"
        "  -q          Quiet mode (key only, no seed output)\n"
        "  -h, --help  Show this help\n\n"
        "Examples:\n"
        "  %s 32                      Generate 32-byte key\n"
        "  %s 32 -q | cfx mac -k -    Pipe key to mac\n",
        prog, prog, prog);
}

int cfx_keygen_run(int argc, char **argv) {
    const char* prog = argv[0];

    if (argc < 2 || (argc == 2 && (strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0))) {
        usage(prog);
        return argc < 2 ? 1 : 0;
    }

    long nbytes = -1;
    const char* seed_in = NULL;
    enum cfx_str_format fmt = CFX_STR_FMT_HEX;
    int quiet = 0;

    #define CHECK_ARG(i) if (i >= argc) { usage(argv[0]); return EXIT_FAILURE; }

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-s") == 0) {
            ++i;
            CHECK_ARG(i);
            seed_in = argv[i];
        } else if (strcmp(argv[i], "-x") == 0) {
            fmt = CFX_STR_FMT_HEX;
        } else if (strcmp(argv[i], "-b64") == 0) {
            fmt = CFX_STR_FMT_BASE64;
        } else if (strcmp(argv[i], "-q") == 0) {
            quiet = 1;
        } else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            usage(prog);
            return 0;
        } else if (argv[i][0] == '-') {
            fprintf(stderr, "Unknown option: %s\n", argv[i]);
            usage(prog);
            return EXIT_FAILURE;
        } else {
            char *end = NULL;
            nbytes = strtol(argv[i], &end, 10);
            if (*end != '\0' || nbytes <= 0) {
                fprintf(stderr, "Invalid bytes: %s\n", argv[i]);
                usage(prog);
                return EXIT_FAILURE;
            }
        }
    }

    if (nbytes <= 0) {
        fprintf(stderr, "error: <bytes> is required\n");
        usage(prog);
        return EXIT_FAILURE;
    }

    uint8_t *buf = malloc((size_t)nbytes);
    if (!buf) {
        fprintf(stderr, "Allocation failed\n");
        return 2;
    }
    unsigned seed = 0;

    if (seed_in) {
        seed = strtoull(seed_in, NULL, 0);
    } else {
        cfx_srand_os();
        seed = cfx_rand_os();
    }

    if (!quiet) printf("using seed 0x%x\n", seed);
    cfx_srand(seed);
    cfx_rand_bytes(buf, (size_t)nbytes);

    if (fmt == CFX_STR_FMT_BASE64) {
        size_t b64_len = cfx_base64_enc_len((size_t)nbytes);
        char* b64 = (char*)malloc(b64_len + 1);
        if (b64) {
            cfx_base64_encode(b64, &b64_len, buf, (size_t)nbytes);
            b64[b64_len] = '\0';
            printf("%s\n", b64);
            free(b64);
        } else {
            fprintf(stderr, "Allocation failed\n");
            free(buf);
            return 2;
        }
    } else {
        for (long i = 0; i < nbytes; i++) {
            printf("%02x", buf[i]);
        }
        printf("\n");
    }

    free(buf);
    return 0;
}
