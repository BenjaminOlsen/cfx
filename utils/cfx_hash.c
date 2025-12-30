/* cfx_hash.c - SHA256 hashing utility */

#include "cfx/sha256.h"
#include "cfx/base64.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
#include <io.h>
#include <fcntl.h>
#endif

#include "cfx_utils.h"

enum output_format {
    FMT_HEX,
    FMT_BASE64
};

static void usage(const char* prog) {
    fprintf(stderr, "Usage: %s [options] [file...]\n", prog);
    fprintf(stderr, "  Compute SHA256 hash of files or stdin.\n\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -s <string>   Hash the given string instead of a file\n");
    fprintf(stderr, "  -x            Output as hex (default)\n");
    fprintf(stderr, "  -b64          Output as base64\n");
    fprintf(stderr, "  -h, --help    Show this help message\n\n");
    fprintf(stderr, "Examples:\n");
    fprintf(stderr, "  %s file.txt           Hash a file\n", prog);
    fprintf(stderr, "  %s -s \"hello world\"   Hash a string\n", prog);
    fprintf(stderr, "  %s -b64 file.txt      Hash with base64 output\n", prog);
    fprintf(stderr, "  echo data | %s        Hash stdin\n", prog);
}

static void print_hash(const uint8_t hash[32], const char* name, enum output_format fmt) {
    if (fmt == FMT_BASE64) {
        size_t b64_len = cfx_base64_enc_len(32);
        char b64[64]; /* SHA256 = 32 bytes -> 44 base64 chars + null */
        cfx_base64_encode(b64, &b64_len, hash, 32);
        b64[b64_len] = '\0';
        printf("%s", b64);
    } else {
        for (int i = 0; i < 32; i++) {
            printf("%02x", hash[i]);
        }
    }
    if (name) {
        printf("  %s", name);
    }
    printf("\n");
}

static int hash_string(const char* str, enum output_format fmt) {
    cfx_sha256_ctx ctx;
    uint8_t hash[32];

    cfx_sha256_init(&ctx);
    cfx_sha256_update(&ctx, (const uint8_t*)str, strlen(str));
    cfx_sha256_final(&ctx, hash);

    print_hash(hash, NULL, fmt);
    return 0;
}

static int hash_file(FILE* f, const char* filename, enum output_format fmt) {
    cfx_sha256_ctx ctx;
    uint8_t hash[32];
    uint8_t buf[8192];
    size_t n;

    cfx_sha256_init(&ctx);
    while ((n = fread(buf, 1, sizeof(buf), f)) > 0) {
        cfx_sha256_update(&ctx, buf, n);
    }

    if (ferror(f)) {
        fprintf(stderr, "Error reading %s\n", filename ? filename : "stdin");
        return 1;
    }

    cfx_sha256_final(&ctx, hash);
    print_hash(hash, filename, fmt);
    return 0;
}

int cfx_hash_run(int argc, char** argv) {
    if (argc == 2 && (strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0)) {
        usage(argv[0]);
        return 0;
    }

    int argi = 1;
    int had_input = 0;
    int rc = 0;
    enum output_format fmt = FMT_HEX;

    while (argi < argc) {
        if (strcmp(argv[argi], "-s") == 0) {
            if (argi + 1 >= argc) {
                fprintf(stderr, "Error: -s requires a string argument\n");
                usage(argv[0]);
                return 1;
            }
            rc |= hash_string(argv[argi + 1], fmt);
            had_input = 1;
            argi += 2;
        } else if (strcmp(argv[argi], "-x") == 0) {
            fmt = FMT_HEX;
            argi++;
        } else if (strcmp(argv[argi], "-b64") == 0) {
            fmt = FMT_BASE64;
            argi++;
        } else if (argv[argi][0] == '-' && argv[argi][1] != '\0') {
            fprintf(stderr, "Unknown option: %s\n", argv[argi]);
            usage(argv[0]);
            return 1;
        } else {
            /* It's a filename */
            FILE* f = fopen(argv[argi], "rb");
            if (!f) {
                perror(argv[argi]);
                rc = 1;
            } else {
                rc |= hash_file(f, argv[argi], fmt);
                fclose(f);
            }
            had_input = 1;
            argi++;
        }
    }

    /* No input specified - read from stdin */
    if (!had_input) {
#ifdef _WIN32
        _setmode(_fileno(stdin), _O_BINARY);
#endif
        rc = hash_file(stdin, "-", fmt);
    }

    return rc;
}
