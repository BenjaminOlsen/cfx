#include "cfx/rand.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>


static void usage(const char *prog) {
    fprintf(stderr, "Usage: %s [-s seed (optional)] -n <bytes>\n", prog);
    exit(1);
}

int main(int argc, char **argv) {
    const char* prog = argv[0];
    long nbytes = -1;
    const char* seed_in = NULL;

    #define CHECK_ARG(i) if (i >= argc) { usage(argv[0]); return EXIT_FAILURE; }
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-n") == 0) {
            ++i;
            CHECK_ARG(i);
            char *end = NULL;
            nbytes = strtol(argv[i], &end, 10);
            if (*end != '\0' || nbytes <= 0) {
                fprintf(stderr, "Invalid value for -n: %s\n", argv[i]);
                usage(prog);
            }
        } else if (strcmp(argv[i], "-s") == 0) {
            ++i;
            CHECK_ARG(i);
            seed_in = argv[i];
        } else {
            usage(prog);
        }
    }

    if (nbytes <= 0) {
        usage(prog);
    }

    uint8_t *buf = malloc((size_t)nbytes);
    if (!buf) {
        fprintf(stderr, "Allocation failed\n");
        return 2;
    }

    if (seed_in) {
        unsigned seed = strtoull(seed_in, NULL, 0);  /* base 0: deduces base from a */;
        printf("got seed: %u\n", seed);
        cfx_srand(seed);
        cfx_randombytes(buf, (size_t)nbytes);
    } else {
        printf("using OS RNG\n");
        cfx_randombytes_os(buf, (size_t)nbytes);
    }

    for (long i = 0; i < nbytes; i++) {
        printf("%02x", buf[i]);
    }
    printf("\n");

    free(buf);
    return 0;
}
