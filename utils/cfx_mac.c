#include "cfx/rand.h"
#include "cfx/poly1305.h"
#include "cfx/base64.h"
#include "cfx_utils_common.h"
#include "cfx_cmd.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#ifdef _WIN32
#include <io.h>
#include <fcntl.h>
#endif


static void usage(const char *prog) {
    fprintf(stderr,
        "Usage: %s [options] [message]\n"
        "  Compute Poly1305 MAC of a message.\n"
        "  If message is omitted, reads from stdin.\n\n"
        "Options:\n"
        "  -k <key>    Key (auto-detects hex vs ASCII, or '-' for stdin)\n"
        "  -kx         Force key as hex\n"
        "  -ka         Force key as ASCII\n"
        "  -kb         Force key as base64\n"
        "  -s <seed>   Seed for random key generation (ignored if -k specified)\n"
        "  -x          Output as hex (default)\n"
        "  -b64        Output as base64\n"
        "  -q          Quiet mode (MAC only, no key output)\n"
        "  -h, --help  Show this help\n\n"
        "Examples:\n"
        "  %s -k mykey \"hello world\"\n"
        "  %s -ka deadbeef \"msg\"       Use 'deadbeef' as ASCII\n"
        "  %s -kx deadbeef \"msg\"       Use 0xdeadbeef as hex\n"
        "  cfx keygen 32 -q | %s -k - \"hello\" -q\n",
        prog, prog, prog, prog, prog);
}

int cfx_mac_run(int argc, char** argv) {
    const char* prog = argv[0];
    if (argc < 2 || (argc == 2 && (strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0))) {
        usage(prog);
        return argc < 2 ? 1 : 0;
    }

    const char* key_in  = NULL;
    const char* seed_in = NULL;
    const char* msg_in  = NULL;
    enum cfx_str_format fmt = CFX_STR_FMT_HEX;
    int quiet = 0;
    enum cfx_str_format key_mode = CFX_STR_FMT_AUTO;

    #define CHECK_ARG(i) \
    do { \
        if (i >= argc) { \
            usage(argv[0]); \
            return EXIT_FAILURE; \
        } \
    } \
    while(0)

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
        } else if (strcmp(argv[i], "-s") == 0) {
            ++i;
            CHECK_ARG(i);
            seed_in = argv[i];
        } else if (strcmp(argv[i], "-x") == 0) {
            fmt = CFX_STR_FMT_HEX;
        } else if (strcmp(argv[i], "-b64") == 0) {
            fmt = CFX_STR_FMT_BASE64;
        } else if (strcmp(argv[i], "-q") == 0) {
            quiet = 1;
        } else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            usage(prog);
            return 0;
        } else if (argv[i][0] == '-') {
            fprintf(stderr, "Unknown option: %s\n", argv[i]);
            usage(prog);
            return EXIT_FAILURE;
        } else {
            msg_in = argv[i];
        }
    }

    uint8_t* msg_buf = NULL;
    size_t len = 0;

    if (!msg_in) {
#ifdef _WIN32
        _setmode(_fileno(stdin), _O_BINARY);
#endif
        size_t cap = 4096;
        msg_buf = (uint8_t*)malloc(cap);
        if (!msg_buf) { perror("malloc"); return EXIT_FAILURE; }

        int c;
        while ((c = getchar()) != EOF) {
            if (len >= cap) {
                cap *= 2;
                uint8_t* newbuf = (uint8_t*)realloc(msg_buf, cap);
                if (!newbuf) { free(msg_buf); perror("realloc"); return EXIT_FAILURE; }
                msg_buf = newbuf;
            }
            msg_buf[len++] = (uint8_t)c;
        }
    } else {
        len = strlen(msg_in);
        msg_buf = (uint8_t*)malloc(len);
        if (!msg_buf) { perror("malloc"); return EXIT_FAILURE; }
        memcpy(msg_buf, msg_in, len);
    }

    const uint8_t* msg = msg_buf;
    uint8_t tag[16];
    uint8_t key[32];
    memset(key, 0, sizeof key);

    int key_from_stdin = (key_in && strcmp(key_in, "-") == 0);

    if (key_from_stdin) {
        if (!msg_in) {
            fprintf(stderr, "error: cannot read both key and message from stdin\n");
            return EXIT_FAILURE;
        }
        char* key_str = cfx_read_line_stdin();
        if (!key_str) {
            fprintf(stderr, "error: failed to read key from stdin\n");
            return EXIT_FAILURE;
        }
        if (cfx_parse_str(key_str, key, sizeof(key), key_mode) < 0) {
            fprintf(stderr, "error: invalid key from stdin\n");
            free(key_str);
            return EXIT_FAILURE;
        }
        free(key_str);
    } else if (!key_in) {
        if (seed_in) {
            unsigned seed = (unsigned)strtoull(seed_in, NULL, 0);
            if (!quiet) printf("got seed: %u\n", seed);
            cfx_srand(seed);
            cfx_rand_bytes((void*)key, sizeof(key));
        } else {
            cfx_rand_bytes_os((void*)key, sizeof(key));
        }
    } else {
        if (cfx_parse_str(key_in, key, sizeof(key), key_mode) < 0) {
            fprintf(stderr, "error: invalid key\n");
            return EXIT_FAILURE;
        }
    }

    cfx_poly1305(tag, msg, len, key);
    free(msg_buf);

    if (!quiet) {
        printf("key: ");
        for (size_t i = 0; i < sizeof(key); ++i) {
            printf("%02x", key[i]);
        }
        printf("\n\nMAC: ");
    }

    if (fmt == CFX_STR_FMT_BASE64) {
        size_t b64_len = cfx_base64_enc_len(sizeof(tag));
        char b64[32];
        cfx_base64_encode(b64, &b64_len, tag, sizeof(tag));
        b64[b64_len] = '\0';
        printf("%s\n", b64);
    } else {
        for (size_t i = 0; i < sizeof(tag); ++i) {
            printf("%02x", tag[i]);
        }
        printf("\n");
    }

    return 0;
}
