/* cfx_hmac.c - HMAC-SHA256/HMAC-SHA512 utility */

#include "cfx/hmac.h"
#include "cfx/base64.h"
#include "common.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
#include <io.h>
#include <fcntl.h>
#endif

#include "cfx_cmd.h"

typedef enum {
    HMAC_SHA256,
    HMAC_SHA512
} hmac_algo_t;

static void usage(const char *prog) {
    printf("Usage: %s [options] -k <hex_key> [file...]\n", prog);
    printf("  Compute HMAC of files or stdin.\n\n");
    printf("Algorithms:\n");
    printf("  (default)       HMAC-SHA256 (32 bytes)\n");
    printf("  --sha512        HMAC-SHA512 (64 bytes)\n\n");
    printf("Options:\n");
    printf("  -k <hex>        Key (hex-encoded, required)\n");
    printf("  -s <string>     HMAC the given string instead of a file\n");
    printf("  -f              Include filename in output\n");
    printf("  -x              Output as hex (default)\n");
    printf("  -b64            Output as base64\n");
    printf("  -h, --help      Show this help message\n\n");
    printf("Examples:\n");
    printf("  %s -k 0b0b0b0b0b0b0b0b0b0b0b0b0b0b0b0b0b0b0b0b -s \"Hi There\"\n", prog);
    printf("  %s --sha512 -k deadbeef file.txt\n", prog);
    printf("  echo -n \"data\" | %s -k 4a656665\n", prog);
}

static void print_mac(const uint8_t *mac, size_t len, const char *name,
                       enum cfx_str_format fmt) {
    if (fmt == CFX_STR_FMT_BASE64) {
        size_t b64_len = cfx_base64_enc_len(len);
        char b64[128];
        cfx_base64_encode(b64, &b64_len, mac, len);
        b64[b64_len] = '\0';
        printf("%s", b64);
    } else {
        for (size_t i = 0; i < len; i++) {
            printf("%02x", mac[i]);
        }
    }
    if (name) {
        printf("  %s", name);
    }
    printf("\n");
}

static int hmac_string(hmac_algo_t algo, const char *str,
                        const uint8_t *key, size_t keylen,
                        enum cfx_str_format fmt) {
    size_t len = strlen(str);
    if (algo == HMAC_SHA512) {
        uint8_t mac[64];
        cfx_hmac_sha512(mac, key, keylen, (const uint8_t *)str, len);
        print_mac(mac, 64, NULL, fmt);
    } else {
        uint8_t mac[32];
        cfx_hmac_sha256(mac, key, keylen, (const uint8_t *)str, len);
        print_mac(mac, 32, NULL, fmt);
    }
    return 0;
}

static int hmac_file(hmac_algo_t algo, FILE *f, const char *filename,
                      const uint8_t *key, size_t keylen,
                      enum cfx_str_format fmt) {
    uint8_t buf[8192];
    size_t n;

    if (algo == HMAC_SHA512) {
        cfx_hmac_sha512_ctx ctx;
        cfx_hmac_sha512_init(&ctx, key, keylen);
        while ((n = fread(buf, 1, sizeof(buf), f)) > 0) {
            cfx_hmac_sha512_update(&ctx, buf, n);
        }
        if (ferror(f)) {
            fprintf(stderr, "Error reading %s\n", filename ? filename : "stdin");
            return 1;
        }
        uint8_t mac[64];
        cfx_hmac_sha512_final(&ctx, mac);
        print_mac(mac, 64, filename, fmt);
    } else {
        cfx_hmac_sha256_ctx ctx;
        cfx_hmac_sha256_init(&ctx, key, keylen);
        while ((n = fread(buf, 1, sizeof(buf), f)) > 0) {
            cfx_hmac_sha256_update(&ctx, buf, n);
        }
        if (ferror(f)) {
            fprintf(stderr, "Error reading %s\n", filename ? filename : "stdin");
            return 1;
        }
        uint8_t mac[32];
        cfx_hmac_sha256_final(&ctx, mac);
        print_mac(mac, 32, filename, fmt);
    }
    return 0;
}

int cfx_hmac_run(int argc, char **argv) {
    if (argc == 2 && (strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0)) {
        usage(argv[0]);
        return 0;
    }

    int rc = 0;
    int show_filename = 0;
    enum cfx_str_format fmt = CFX_STR_FMT_HEX;
    hmac_algo_t algo = HMAC_SHA256;
    uint8_t key[128] = {0};
    size_t keylen = 0;

    const char *strings[64];
    const char *files[64];
    int nstrings = 0;
    int nfiles = 0;

    for (int argi = 1; argi < argc; argi++) {
        if (strcmp(argv[argi], "-s") == 0) {
            if (argi + 1 >= argc) {
                fprintf(stderr, "Error: -s requires a string argument\n");
                return 1;
            }
            if (nstrings < 64) strings[nstrings++] = argv[argi + 1];
            argi++;
        } else if (strcmp(argv[argi], "-k") == 0) {
            if (argi + 1 >= argc) {
                fprintf(stderr, "Error: -k requires a hex key argument\n");
                return 1;
            }
            int parsed = cfx_parse_hex_auto(argv[argi + 1], key, sizeof(key));
            if (parsed < 0) {
                fprintf(stderr, "Error: invalid hex key\n");
                return 1;
            }
            keylen = (size_t)parsed;
            argi++;
        } else if (strcmp(argv[argi], "--sha512") == 0) {
            algo = HMAC_SHA512;
        } else if (strcmp(argv[argi], "-x") == 0) {
            fmt = CFX_STR_FMT_HEX;
        } else if (strcmp(argv[argi], "-b64") == 0) {
            fmt = CFX_STR_FMT_BASE64;
        } else if (strcmp(argv[argi], "-f") == 0) {
            show_filename = 1;
        } else if (argv[argi][0] == '-' && argv[argi][1] != '\0') {
            fprintf(stderr, "Unknown option: %s\n", argv[argi]);
            usage(argv[0]);
            return 1;
        } else {
            if (nfiles < 64) files[nfiles++] = argv[argi];
        }
    }

    if (keylen == 0) {
        fprintf(stderr, "Error: -k <key> is required\n");
        return 1;
    }

    for (int i = 0; i < nstrings; i++) {
        rc |= hmac_string(algo, strings[i], key, keylen, fmt);
    }

    for (int i = 0; i < nfiles; i++) {
        FILE *f = fopen(files[i], "rb");
        if (!f) {
            perror(files[i]);
            rc = 1;
        } else {
            rc |= hmac_file(algo, f, show_filename ? files[i] : NULL, key, keylen, fmt);
            fclose(f);
        }
    }

    if (nstrings == 0 && nfiles == 0) {
#ifdef _WIN32
        _setmode(_fileno(stdin), _O_BINARY);
#endif
        rc = hmac_file(algo, stdin, show_filename ? "-" : NULL, key, keylen, fmt);
    }

    return rc;
}
