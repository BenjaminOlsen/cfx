/* cfx_keygen.c - cryptographic key generation utility */

#include "cfx/rand.h"
#include "cfx/base64.h"
#include "cfx/memory.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#ifdef _WIN32
#include <io.h>
#include <fcntl.h>
#endif

#include "cfx_cmd.h"
#include "misc.h"

static void usage(const char* prog) {
    printf("Usage: %s [options] <bytes>\n", prog);
    printf("  Generate cryptographically secure random bytes.\n\n");
    printf("Options:\n");
    printf("  -x          Output as hex (default)\n");
    printf("  -b64        Output as base64\n");
    printf("  -r          Output raw bytes (binary)\n");
    printf("  -h, --help  Show this help\n\n");
    printf("Examples:\n");
    printf("  %s 32           Generate 32-byte key (hex)\n", prog);
    printf("  %s 16 -b64      Generate 16-byte key (base64)\n", prog);
    printf("  %s 32 -r > key  Write raw key to file\n", prog);
}

int cfx_keygen_run(int argc, char** argv) {
    if (argc < 2 || (argc == 2 && (strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0))) {
        usage(argv[0]);
        return argc < 2 ? 1 : 0;
    }

    long nbytes = -1;
    enum cfx_str_format fmt = CFX_STR_FMT_HEX;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-x") == 0) {
            fmt = CFX_STR_FMT_HEX;
        } else if (strcmp(argv[i], "-b64") == 0) {
            fmt = CFX_STR_FMT_BASE64;
        } else if (strcmp(argv[i], "-r") == 0) {
            fmt = CFX_STR_FMT_BINARY;
        } else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            usage(argv[0]);
            return 0;
        } else if (argv[i][0] == '-') {
            fprintf(stderr, "Unknown option: %s\n", argv[i]);
            usage(argv[0]);
            return 1;
        } else {
            char* end = NULL;
            nbytes = strtol(argv[i], &end, 10);
            if (*end != '\0' || nbytes <= 0 || nbytes > 1024 * 1024) {
                fprintf(stderr, "Invalid byte count: %s (must be 1-%d)\n", argv[i], 1024 * 1024);
                return 1;
            }
        }
    }

    if (nbytes <= 0) {
        fprintf(stderr, "Error: <bytes> is required\n");
        usage(argv[0]);
        return 1;
    }

    uint8_t* buf = malloc((size_t)nbytes);
    if (!buf) {
        fprintf(stderr, "Allocation failed\n");
        return 1;
    }

    /* seed ChaCha20 CSPRNG from OS entropy, then generate bytes */
    cfx_srand_os();
    cfx_rand_bytes(buf, (size_t)nbytes);

    if (fmt == CFX_STR_FMT_BINARY) {
#ifdef _WIN32
        _setmode(_fileno(stdout), _O_BINARY);
#endif
        fwrite(buf, 1, (size_t)nbytes, stdout);
    } else if (fmt == CFX_STR_FMT_BASE64) {
        size_t b64_len = cfx_base64_enc_len((size_t)nbytes);
        char* b64 = malloc(b64_len + 1);
        if (!b64) {
            fprintf(stderr, "Allocation failed\n");
            free(buf);
            return 1;
        }
        cfx_base64_encode(b64, &b64_len, buf, (size_t)nbytes);
        b64[b64_len] = '\0';
        printf("%s\n", b64);
        free(b64);
    } else {
        for (long i = 0; i < nbytes; i++) {
            printf("%02x", buf[i]);
        }
        printf("\n");
    }

    /* zero sensitive data */
    memset(buf, 0, (size_t)nbytes);
    free(buf);
    return 0;
}
