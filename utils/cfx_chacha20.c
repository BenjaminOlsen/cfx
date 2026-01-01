#include "cfx/chacha20.h"
#include "cfx/base64.h"
#include "misc.h"
#include "cfx_cmd.h"

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#ifdef _WIN32
#include <io.h>
#include <fcntl.h>
#endif



static void usage(const char* name) {
    fprintf(stderr,
        "Usage: %s [options] [text]\n"
        "  Encrypt text with ChaCha20 (RFC8439 layout).\n"
        "  If text is omitted, reads from stdin.\n\n"
        "Options:\n"
        "  -k <key>    Key (auto-detects hex vs ASCII, or '-' for stdin)\n"
        "  -kx         Force key as hex\n"
        "  -ka         Force key as ASCII\n"
        "  -kb         Force key as base64\n"
        "  -c <ctr>    Initial 32-bit counter (default 0)\n"
        "  -n <nonce>  Nonce as 0x24_hex_chars (96-bit, default zeros)\n"
        "  -x          Output as hex (default)\n"
        "  -b64        Output as base64\n"
        "  -v          Verbose output\n"
        "  -h, --help  Show this help\n\n"
        "Examples:\n"
        "  %s -k mykey \"hello world\"\n"
        "  %s -ka deadbeef \"msg\"       Use 'deadbeef' as ASCII\n"
        "  cfx keygen 32 -q | %s -k - \"secret\"\n",
        name, name, name, name);
}

int cfx_chacha20_run(int argc, char** argv) {
    if (argc < 2 || (argc == 2 && (strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0))) {
        usage(argv[0]);
        return argc < 2 ? EXIT_FAILURE : 0;
    }

    const char* key_in   = NULL;
    const char* nonce_in = NULL;
    const char* pt       = NULL;
    uint32_t counter     = 0;
    int verbose = 0;
    enum cfx_str_format fmt = CFX_STR_FMT_HEX;
    enum cfx_str_format key_mode = CFX_STR_FMT_AUTO;

    #define CHECK_ARG(i) if (i >= argc) { usage(argv[0]); return EXIT_FAILURE; }

    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-k") == 0) {
            ++i;
            CHECK_ARG(i);
            key_in = argv[i];
        } else if (strcmp(argv[i], "-kx") == 0) {
            ++i;
            CHECK_ARG(i);
            key_in = argv[i];
            key_mode = CFX_STR_FMT_HEX;
        } else if (strcmp(argv[i], "-ka") == 0) {
            ++i;
            CHECK_ARG(i);
            key_in = argv[i];
            key_mode = CFX_STR_FMT_ASCII;
        } else if (strcmp(argv[i], "-kb") == 0) {
            ++i;
            CHECK_ARG(i);
            key_in = argv[i];
            key_mode = CFX_STR_FMT_BASE64;
        } else if (strcmp(argv[i], "-c") == 0) {
            ++i;
            CHECK_ARG(i);
            errno = 0;
            unsigned long v = strtoul(argv[i], NULL, 0);
            if (errno) { perror("counter"); return EXIT_FAILURE; }
            counter = (uint32_t)v;
        } else if (strcmp(argv[i], "-n") == 0) {
            ++i;
            CHECK_ARG(i);
            nonce_in = argv[i];
        } else if (strcmp(argv[i], "-v") == 0) {
            verbose = 1;
        } else if (strcmp(argv[i], "-x") == 0) {
            fmt = CFX_STR_FMT_HEX;
        } else if (strcmp(argv[i], "-b64") == 0) {
            fmt = CFX_STR_FMT_BASE64;
        } else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            usage(argv[0]);
            return 0;
        } else if (argv[i][0] == '-') {
            fprintf(stderr, "Unknown option: %s\n", argv[i]);
            usage(argv[0]);
            return EXIT_FAILURE;
        } else {
            pt = argv[i];
        }
    }

    uint8_t* pt_buf = NULL;
    size_t len = 0;

    if (!pt) {
#ifdef _WIN32
        _setmode(_fileno(stdin), _O_BINARY);
#endif
        size_t cap = 4096;
        pt_buf = (uint8_t*)malloc(cap);
        if (!pt_buf) { perror("malloc"); return EXIT_FAILURE; }

        int c;
        while ((c = getchar()) != EOF) {
            if (len >= cap) {
                cap *= 2;
                uint8_t* newbuf = (uint8_t*)realloc(pt_buf, cap);
                if (!newbuf) { free(pt_buf); perror("realloc"); return EXIT_FAILURE; }
                pt_buf = newbuf;
            }
            pt_buf[len++] = (uint8_t)c;
        }
    } else {
        len = strlen(pt);
        pt_buf = (uint8_t*)malloc(len);
        if (!pt_buf) { perror("malloc"); return EXIT_FAILURE; }
        memcpy(pt_buf, pt, len);
    }

    int key_from_stdin = (key_in && strcmp(key_in, "-") == 0);

    uint8_t key[32] = {0};

    if (key_from_stdin) {
        if (!pt) {
            fprintf(stderr, "error: cannot read both key and plaintext from stdin\n");
            free(pt_buf);
            return EXIT_FAILURE;
        }
        char* key_str = cfx_read_line_stdin();
        if (!key_str) {
            fprintf(stderr, "error: failed to read key from stdin\n");
            free(pt_buf);
            return EXIT_FAILURE;
        }
        if (cfx_parse_str(key_str, key, sizeof(key), key_mode) < 0) {
            fprintf(stderr, "error: invalid key from stdin\n");
            free(key_str);
            free(pt_buf);
            return EXIT_FAILURE;
        }
        free(key_str);
    } else if (!key_in) {
        const char default_key[] = "1234567890";
        if (verbose) fprintf(stderr, "using default key\n");
        memcpy(key, default_key, strlen(default_key));
    } else {
        if (cfx_parse_str(key_in, key, sizeof(key), key_mode) < 0) {
            fprintf(stderr, "error: invalid key\n");
            free(pt_buf);
            return EXIT_FAILURE;
        }
    }

    if (verbose) {
        fprintf(stderr, "key: ");
        for (size_t i = 0; i < sizeof key; ++i) {
            fprintf(stderr, "%02x ", key[i]);
        }
        fprintf(stderr, "\n");
    }

    uint8_t nonce[12] = {0};
    if (nonce_in) {
        if (strncmp(nonce_in, "0x", 2) != 0 ||
            cfx_parse_hex(nonce_in + 2, nonce, sizeof(nonce)) != 0) {
            fprintf(stderr, "error: -n expects hex: with exactly 24 hex chars, using default\n");
        }
    }

    uint8_t* ct = (uint8_t*)malloc(len ? len : 1);
    if (!ct) { free(pt_buf); perror("malloc"); return EXIT_FAILURE; }

    cfx_chacha20_encrypt(key, counter, nonce, pt_buf, len, ct);
    free(pt_buf);

    if (verbose) printf("ciphertext (%zu bytes):\n", len);

    if (fmt == CFX_STR_FMT_BASE64) {
        size_t b64_len = cfx_base64_enc_len(len);
        char* b64 = (char*)malloc(b64_len + 1);
        if (b64) {
            cfx_base64_encode(b64, &b64_len, ct, len);
            b64[b64_len] = '\0';
            printf("%s\n", b64);
            free(b64);
        }
    } else {
        for (size_t i = 0; i < len; ++i) {
            printf("%02x", ct[i]);
        }
        printf("\n");
    }

    free(ct);
    return 0;
}
