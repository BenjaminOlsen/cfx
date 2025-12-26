#include "cfx/rand.h"
#include "cfx/base64.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>


static void usage(const char *prog) {
    fprintf(stderr, 
        "Usage: %s  [-x output hex] [-s seed (optional)] -n <bytes>\n"
        "   create a random string of n bytes, encoded in base64 (optionally hex)\n"
        , prog);
    exit(1);
}

int main(int argc, char **argv) {
    const char* prog = argv[0];
    long nbytes = -1;
    const char* seed_in = NULL;
    int hex = 0;

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
        } else if (strcmp(argv[i], "-x") == 0) {
            hex = 1;
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
    unsigned seed = 0;

    if (seed_in) {
        seed = strtoull(seed_in, NULL, 0);  /* base 0: deduces base from a */;
    } else {
        cfx_srand_os();
        seed = cfx_rand_os();
    }

    printf("using seed 0x%x\n", seed);
    cfx_srand(seed);
    cfx_rand_bytes(buf, (size_t)nbytes);

    if (hex) {
        for (long i = 0; i < nbytes; i++) {
            printf("%02x", buf[i]);
        }
        printf("\n");
    } else {
        size_t nchars = 0;
        cfx_base64_encode(NULL, &nchars, buf, nbytes);
        char* s = (char*)malloc(nchars);
        if (s) {
            cfx_base64_encode(s, &nchars, buf, nbytes);
            printf("%s\n", s);
        } else {
            printf("problem allocating buffer\n");
        }
        free(s);
    }

    free(buf);
    return 0;
}
