#include "cfx/rand.h"
#include "cfx/crypto.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>


static void usage(const char *prog) {
    fprintf(stderr, "Usage: %s <message>\n"
            "   outputs poly1305 signature of the message.\n",
            prog);
    exit(1);
}

int main(int argc, char** argv) {
    const char* prog = argv[0];
    if (argc < 2) {
        usage(prog);
        return EXIT_FAILURE;
    }

    const uint8_t* msg = (uint8_t*)argv[1];
    size_t len = strlen(argv[1]);
    uint8_t tag[16];
    uint8_t key[32];
    cfx_randombytes((void*)key, sizeof(key));
    cfx_poly1305_auth(key, msg, len, tag);

    for (size_t i = 0; i < sizeof(tag); ++i) {
        printf("%02x ", tag[i]);
    }
    printf("\n");
    return 0;
}
