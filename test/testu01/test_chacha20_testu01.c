
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/* TestU01 includes */
#include "unif01.h"
#include "bbattery.h"

#include "cfx/chacha20.h"
#include "cfx/rand.h"

static uint8_t g_key[32];
static uint8_t g_nonce[12];
static uint32_t g_counter = 0;

static void init_from_seed(uint64_t seed) {
    /* Very simple reproducible mixer from 64-bit seed
       into key[32] and nonce[12]. This is *not* crypto,
       just a way to vary the test generator deterministically.

       ref: LCG - https://en.wikipedia.org/wiki/Linear_congruential_generator
    */
    uint32_t x = (uint32_t)seed;
    uint32_t i;

    for (i = 0; i < sizeof(g_key); i++) {
        x = 1664525u * x + 1013904223u;
        g_key[i] = (uint8_t)(x >> 24);
    }

    for (i = 0; i < sizeof(g_nonce); i++) {
        x = 1664525u * x + 1013904223u;
        g_nonce[i] = (uint8_t)(x >> 24);
    }

    g_counter = 0;
}

/* ------------------------------------------------------------------ */
/* 32-bit generator interface for TestU01                             */
static uint32_t cfx_chacha32(void) {
    static uint8_t buf[64];
    static size_t idx = 64;  /* force initial refill */

    if (idx >= 64) {
        /* Fill next block */
        cfx_chacha20_block_rfc8439(g_key, g_counter++, g_nonce, buf);
        idx = 0;
    }

    uint32_t v =
        ((uint32_t)buf[idx + 0]      ) |
        ((uint32_t)buf[idx + 1] <<  8) |
        ((uint32_t)buf[idx + 2] << 16) |
        ((uint32_t)buf[idx + 3] << 24);
    idx += 4;
    return v;
}

typedef enum {
    MODE_SMALLCRUSH,
    MODE_CRUSH,
    MODE_BIGCRUSH
} test_mode_t;

static void usage(const char *prog) {
    fprintf(stderr,
        "Usage: %s [--smallcrush|--crush|--bigcrush] [--seed=<integer>]\n"
        "\n"
        "  --smallcrush   Run the SmallCrush battery (default)\n"
        "  --crush        Run the Crush battery\n"
        "  --bigcrush     Run the BigCrush battery (can take a long time)\n"
        "  --seed=N       Use integer N to derive key/nonce (default: 0xC0FFEE)\n"
        "  --help, -h     Print this menu and exit\n",
        prog
    );
}

int main(int argc, char **argv) {
    test_mode_t mode = MODE_SMALLCRUSH;
    uint64_t seed = 0x0000000000C0FFEEull;
    const char* prog = argv[0];

    for (int i = 1; i < argc; i++) {
        const char* arg = argv[i];
        
        if (strcmp(arg, "--smallcrush") == 0) {
            mode = MODE_SMALLCRUSH;
        } else if (strcmp(arg, "--crush") == 0) {
            mode = MODE_CRUSH;
        } else if (strcmp(arg, "--bigcrush") == 0) {
            mode = MODE_BIGCRUSH;
        } else if (strncmp(arg, "--seed=", 7) == 0) {
            char* end = NULL;
            seed = strtoull(arg + 7, &end, 0);  /* base 0: deduces base from a */
            if (end == arg + 8) {
                fprintf(stderr, "Invalid seed: %s\n\n", arg + 8);
                usage(prog);
                return EXIT_FAILURE;
            }
        } else if ((strcmp(arg, "--help") == 0) || (strcmp(arg, "-h") == 0)) {
            usage(prog);
            return EXIT_FAILURE;
        } else {
            fprintf(stderr, "Unknown option: %s\n\n", arg);
            usage(prog);
            return EXIT_FAILURE;
        }
    }

    init_from_seed(seed);

    printf("cfx_chacha_testu01: seed = 0x%016llX\n",
           (unsigned long long)seed);

    /* Wrap our 32-bit generator for TestU01 */
    unif01_Gen *gen = unif01_CreateExternGenBits("cfx_chacha20", cfx_chacha32);

    switch (mode) {
        case MODE_SMALLCRUSH:
            printf("Running SmallCrush...\n");
            bbattery_SmallCrush(gen);
            break;
        case MODE_CRUSH:
            printf("Running Crush...\n");
            bbattery_Crush(gen);
            break;
        case MODE_BIGCRUSH:
            printf("Running BigCrush...\n");
            bbattery_BigCrush(gen);
            break;
    }

    unif01_DeleteExternGenBits(gen);
    return EXIT_SUCCESS;
}
