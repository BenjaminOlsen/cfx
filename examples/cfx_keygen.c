#include "cfx/rand.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>


static void usage(const char *prog) {
    fprintf(stderr, "Usage: %s -n <bytes>\n", prog);
    exit(1);
}

int main(int argc, char **argv) {
    const char* prog = argv[0];
    long nbytes = -1;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-n") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "Missing argument for -n\n");
                usage(prog);
            }
            char *end = NULL;
            nbytes = strtol(argv[++i], &end, 10);
            if (*end != '\0' || nbytes <= 0) {
                fprintf(stderr, "Invalid value for -n: %s\n", argv[i]);
                usage(prog);
            }
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

    cfx_randombytes(buf, (size_t)nbytes);

    for (long i = 0; i < nbytes; i++) {
        printf("%02x", buf[i]);
    }
    printf("\n");

    free(buf);
    return 0;
}
