
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/* TestU01 includes */
#include "unif01.h"
#include "bbattery.h"

#include "cfx/chacha20.h"
#include "cfx/poly1305.h"
#include "cfx/rand.h"

typedef uint32_t (*test_fn)(void);
typedef void (*seed_fn)(unsigned);


/* .................................................................. */
/* chacha20 */

static uint8_t g_chacha20_key[32];
static uint8_t g_chacha20_nonce[12];
static uint32_t g_chacha20_counter = 0;

static void lcg(uint32_t seed, uint8_t *data, size_t len) {
    /* Very simple reproducible mixer from 64-bit seed
       This is *not* crypto, just a way to vary the test generator deterministically.

       ref: LCG - https://en.wikipedia.org/wiki/Linear_congruential_generator
    */
    uint32_t x = seed;
    uint32_t i;

    for (i = 0; i < len; i++) {
        x = 1664525u * x + 1013904223u;
        data[i] = (uint8_t)(x >> 24);
    }
}
static void chacha20_seed(uint32_t seed) {
    
    lcg(seed, g_chacha20_key, sizeof g_chacha20_key);
    lcg(seed, g_chacha20_nonce, sizeof g_chacha20_nonce);

    g_chacha20_counter = 0;
}

static uint32_t cfx_chacha20_gen(void) {
    static uint8_t buf[64];
    static size_t idx = 64;  /* force initial refill */

    if (idx >= 64) {
        /* Fill next block */
        cfx_chacha20_block_rfc8439(g_chacha20_key, g_chacha20_counter++, g_chacha20_nonce, buf);
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

/* .................................................................. */
/* poly1305 - this is a toy example not suitable for crypto use at all! */
uint8_t g_poly1305_key[32];
static uint64_t g_poly1305_counter = 0;

static void poly1305_seed(unsigned seed) {
    lcg(seed, g_poly1305_key, sizeof g_poly1305_key);
    g_poly1305_counter = 0;
}

static uint32_t cfx_poly1305_gen(void) {
    uint8_t tag[16];
    uint8_t msg[8];

    /* Use counter as message */
    uint64_t ctr = g_poly1305_counter++;
    for (size_t i = 0; i < 8; ++i) {
        msg[i] = (uint8_t)(ctr & 0xFFu);
        ctr >>= 8;
    }
    cfx_poly1305_mac(g_poly1305_key, msg, sizeof msg, tag);

    uint32_t r =
        ((uint32_t)tag[0]      ) |
        ((uint32_t)tag[1] <<  8) |
        ((uint32_t)tag[2] << 16) |
        ((uint32_t)tag[3] << 24);
    
    return r;
}

/* .................................................................. */
/* cfx_rand */
/* seed with cfx_srand (cfx/rand.h) */
static uint32_t cfx_rand_gen(void) {
    uint32_t r;
    cfx_randombytes((void*)&r, sizeof(r));
    return r;
}

/* .................................................................. */
/* c std library rand */
/* seed with srand (stdlib.h) */
static uint32_t rand_gen(void) {
    uint32_t r = (uint32_t)(rand() & 0x7FFFFFFF);
    r ^= 1;
    r ^= (uint32_t)(rand() & 0x7FFFFFFF);
    // debug:
    static int callcnt = 0;
    printf("[%d] rand: 0x%x\n", ++callcnt, r);
    return r;
}

typedef enum {
    MODE_SMALLCRUSH,
    MODE_CRUSH,
    MODE_BIGCRUSH
} test_mode_t;

typedef struct {
    const char *name;     /* name passed to TestU01 */
    const char *cli_opt;  /* command-line selector, e.g. "--chacha20" */
    test_fn     fn;       /* pointer to the RNG function */
    seed_fn     sfn;
} gen_desc_t;

/* For now we only have one, but this scales nicely */
static const gen_desc_t g_generators[] = {
    {"cfx_chacha20", "--chacha20", cfx_chacha20_gen,  chacha20_seed},
    {"cfx_poly1305", "--poly1305", cfx_poly1305_gen, poly1305_seed},
    {"cfx_rand", "--cfx_rand", cfx_rand_gen, cfx_srand},
    {"cfx_splitmix32", "--cfx_splitmix32", cfx_splitmix32, cfx_splitmix_seed},
    {"cfx_pcg32", "--cfx_pcg32", cfx_pcg32, cfx_pcg_seed},
    {"rand", "--rand", rand_gen, srand}

    /* todo later:
       { "cfx_xoshiro256ss", "--xoshiro256ss", cfx_xoshiro256ss_32 },
       ...
    */
};

static const size_t g_generators_count = sizeof(g_generators) / sizeof(g_generators[0]);

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
        g_generators[0].name
    );
    for (size_t i = 0; i < g_generators_count; ++i) {
        fprintf(stderr, "                   %s\n", g_generators[i].name);
    }
}


int main(int argc, char **argv) {
    test_mode_t mode = MODE_SMALLCRUSH;
    uint32_t seed = 0x00C0FFEEu;
    const char* prog = argv[0];

    const gen_desc_t *selected_gen = &g_generators[0]; /* default */
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
            for (size_t j = 0; j < g_generators_count; ++j) {
                if (strcmp(name, g_generators[j].name) == 0) {
                    selected_gen = &g_generators[j];
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

    selected_gen->sfn(seed);
    /* Wrap our 32-bit generator for TestU01 */
    unif01_Gen* gen = unif01_CreateExternGenBits((char *)selected_gen->name, selected_gen->fn);

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
