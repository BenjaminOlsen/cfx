
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "cfx/rand.h"
#include "cfx/memory.h"

/* TestU01 includes */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmacro-redefined"
#include "unif01.h"
#include "bbattery.h"
#pragma GCC diagnostic pop

#if CFX_HAVE_OPENSSL

#include <openssl/rand.h>

static inline void openssl_bytes(void *buf, size_t len) {
    RAND_bytes((unsigned char*)(buf), (int)len);
}

static inline uint32_t openssl_gen32(void) {
    unsigned char buf[4];
    RAND_bytes(buf, sizeof buf);
    uint32_t v = cfx_load32_le((void*)buf);
    return v;
}

static inline void openssl_seed(uint32_t seed) {
    RAND_seed((void*)&seed, sizeof seed);
}

static const cfx_rand_desc_t g_openssl_desc = {
    "openssl_rand",
    /* gen32 = */ openssl_gen32,
    /* seed  = */ openssl_seed,
    /* bytes = */ openssl_bytes
};

#endif

static void usage(const char *prog) {
    fprintf(stderr,
        "Usage: %s [--smallcrush|--crush|--bigcrush] [--seed=<integer>] [--rng=<name>]\n"
        "\n"
        "  --smallcrush   Run the SmallCrush battery (default)\n"
        "  --crush        Run the Crush battery\n"
        "  --bigcrush     Run the BigCrush battery (can take a long time)\n"
        "  --seed=N       Use integer N to derive key/nonce (default: 0xC0FFEE)\n"
        "  --rng=name     Select RNG (default: %s)\n"
        "                 Available RNGs:\n",
        prog,
        g_rand_gens[0].name
    );
    for (size_t i = 0; i < g_rand_gen_cnt; ++i) {
        fprintf(stderr, "                   %s\n", g_rand_gens[i].name);
    }
    fprintf(stderr, "                   %s\n", g_openssl_desc.name);

}


typedef enum {
    MODE_SMALLCRUSH,
    MODE_CRUSH,
    MODE_BIGCRUSH
} test_mode_t;


int main(int argc, char **argv) {
    test_mode_t mode = MODE_SMALLCRUSH;
    uint32_t seed = 0x00C0FFEEu;
    const char* prog = argv[0];

    const cfx_rand_desc_t *selected_gen = &g_rand_gens[0]; /* default */
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
        } else if (strncmp(arg, "--rng=", 6) == 0) {
            const char *name = arg + 6;
            int found = 0;
            for (size_t j = 0; j < g_rand_gen_cnt; ++j) {
                if (strcmp(name, g_rand_gens[j].name) == 0) {
                    selected_gen = &g_rand_gens[j];
                    found = 1;
                    break;
                } else if (strcmp(name, g_openssl_desc.name) == 0) {
                    selected_gen = &g_openssl_desc;
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
            return EXIT_FAILURE;
        } else {
            fprintf(stderr, "Unknown option: %s\n\n", arg);
            usage(prog);
            return EXIT_FAILURE;
        }
    }

    printf("selected gen: %s: seed = 0x%016llX\n", selected_gen->name,
           (unsigned long long)seed);

    selected_gen->seed(seed);
    /* Wrap our 32-bit generator for TestU01 */
    unif01_Gen* gen = unif01_CreateExternGenBits((char *)selected_gen->name, selected_gen->rng32);

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
