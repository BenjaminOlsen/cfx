#include "../test/stats/rand_gen.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static void usage(const char* prog) {
    fprintf(stderr,
        "Usage: %s --seed=<seed|0x...> --rng=<generator name>\n"
        "  --seed=N       Use integer N to derive key/nonce (default: 0xC0FFEE)\n"
        "  --n=<integer>  Produce n bytes (default 4)\n"
        "  --rng=name     Select RNG (default: %s)\n"
        "                 Available RNGs:\n",
        prog,
        g_rand_gens[0].name
    );
    for (size_t i = 0; i < g_rand_gen_cnt; ++i) {
        fprintf(stderr, "                   %s\n", g_rand_gens[i].name);
    }
}


int main(int argc, char** argv) {
    struct timespec res;
    clock_gettime(CLOCK_REALTIME, &res);
    long ns = res.tv_nsec;
    uint32_t seed = (uint32_t)ns;
    size_t n = 4;
    const char* prog = argv[0];

    const rand_desc_t* rand_gen = &g_rand_gens[0]; /* default */
    for (int i = 1; i < argc; i++) {
        const char* arg = argv[i];

       if (strncmp(arg, "--seed=", 7) == 0) {
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
            return EXIT_FAILURE;
        } else if ((strncmp(arg, "--n=", 4) == 0)) {
            n = strtoull(arg+4, NULL, 0);
        } else {
            fprintf(stderr, "Unknown option: %s\n\n", arg);
            usage(prog);
            return EXIT_FAILURE;
        }
    }

    printf("selected gen: %s: seed = 0x%08x\n\n.........\n", rand_gen->name, seed);

    rand_gen->seed(seed);
    uint8_t* bytes = (uint8_t*)malloc(n);
    if (bytes) {
        size_t i = 0;
        while (i < n) {
            uint32_t v = rand_gen->gen32();
            for (size_t vn = 0; vn < sizeof(v) && (i < n); vn++) {
                bytes[i] = v & 0xFF;
                v >>= 8;
                ++i;
            }
        }
    } else {
        return EXIT_FAILURE;
    }

    for (size_t i = 0; i < n; ++i) {
        printf("%02x", bytes[i]);
    }
    printf("\n");
    free(bytes);
    return EXIT_SUCCESS;
}
