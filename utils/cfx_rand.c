#include "cfx/rand.h"
#include "cfx/base64.h"
#include "cfx/compat.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "cfx_cmd.h"
#include "misc.h"


static void usage(const char* prog) {
    fprintf(stderr,
        "Usage: %s [options] <num_bytes>\n"
        "  Generate random bytes.\n\n"
        "Options:\n"
        "  --seed=N       Use integer N to derive key/nonce (default: random)\n"
        "  --rng=name     Select RNG (default: %s)\n"
        "  -x             Output as hex (default)\n"
        "  -b64           Output as base64\n"
        "  -v, --verbose  Verbose output\n"
        "  -h, --help     Show this help\n\n"
        "Available RNGs:\n",
        prog,
        g_rand_gens[0].name
    );
    for (size_t i = 0; i < g_rand_gen_cnt; ++i) {
        fprintf(stderr, "  %s\n", g_rand_gens[i].name);
    }
}

static void print_bytes(const uint8_t* bytes, size_t n, enum cfx_str_format fmt) {
    if (fmt == CFX_STR_FMT_BASE64) {
        size_t b64_len = cfx_base64_enc_len(n);
        char* b64 = (char*)malloc(b64_len + 1);
        if (b64) {
            cfx_base64_encode(b64, &b64_len, bytes, n);
            b64[b64_len] = '\0';
            printf("%s\n", b64);
            free(b64);
        }
    } else {
        for (size_t i = 0; i < n; ++i) {
            printf("%02x", bytes[i]);
        }
        printf("\n");
    }
}

int cfx_rand_run(int argc, char** argv) {
    uint64_t ns = cfx_time_ns();
    uint32_t seed = (uint32_t)ns;

    size_t n = 4;
    int verbose = 0;
    enum cfx_str_format fmt = CFX_STR_FMT_HEX;
    const char* prog = argv[0];

    const cfx_rand_desc_t* rand_gen = &g_rand_gens[0]; /* default */
    for (int i = 1; i < argc; i++) {
        const char* arg = argv[i];

        if (strncmp(arg, "--seed=", 7) == 0) {
            char* end = NULL;
            seed = strtoull(arg + 7, &end, 0);
            if (end == arg + 7) {
                fprintf(stderr, "Invalid seed: %s\n\n", arg + 7);
                usage(prog);
                return EXIT_FAILURE;
            }
        } else if (strncmp(arg, "--rng=", 6) == 0) {
            const char *name = arg + 6;
            int found = 0;
            for (size_t j = 0; j < g_rand_gen_cnt; ++j) {
                if (strcmp(name, g_rand_gens[j].name) == 0) {
                    rand_gen = &g_rand_gens[j];
                    found = 1;
                    break;
                }
            }
            if (!found) {
                fprintf(stderr, "Unknown RNG: %s\n\n", name);
                usage(prog);
                return EXIT_FAILURE;
            }
        } else if ((strcmp(arg, "--help") == 0) || (strcmp(arg, "-h") == 0)) {
            usage(prog);
            return 0;
        } else if ((strcmp(arg, "-v") == 0) || (strcmp(arg, "--verbose") == 0)) {
            verbose = 1;
        } else if (strcmp(arg, "-x") == 0) {
            fmt = CFX_STR_FMT_HEX;
        } else if (strcmp(arg, "-b64") == 0) {
            fmt = CFX_STR_FMT_BASE64;
        } else {
            n = strtoull(arg, NULL, 0);
            if (errno) {
                perror("num bytes");
                return EXIT_FAILURE;
            }
        }
    }

    if (verbose) printf("selected gen: %s: seed = 0x%08x\n.........\n", rand_gen->name, seed);

    rand_gen->seed(seed);
    uint8_t* bytes = (uint8_t*)malloc(n);
    if (bytes) {
        size_t i = 0;
        while (i < n) {
            uint32_t v = rand_gen->rng32();
            for (size_t vn = 0; vn < sizeof(v) && (i < n); vn++) {
                bytes[i] = v & 0xFF;
                v >>= 8;
                ++i;
            }
        }
    } else {
        return EXIT_FAILURE;
    }

    print_bytes(bytes, n, fmt);
    free(bytes);
    return EXIT_SUCCESS;
}
