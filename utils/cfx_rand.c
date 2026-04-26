#include "cfx/rand.h"
#include "cfx/base64.h"
#include "cfx/compat.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#ifdef _WIN32
#include <io.h>
#include <fcntl.h>
#endif

#include "cfx_cmd.h"
#include "common.h"


static void usage(const char* prog) {
    fprintf(stderr,
        "Usage: %s [options] <num_bytes>\n"
        "  Generate random bytes.\n\n"
        "Options:\n"
        "  -n <num>       Number of bytes (also accepted as positional)\n"
        "  --seed, -s <N> Use integer N as seed (default: random)\n"
        "  --rng <name>   Select RNG (default: %s)\n"
        "  -x             Output as hex (default)\n"
        "  -b64           Output as base64\n"
        "  -b, --bin      Output as raw binary\n"
        "  -v, --verbose  Verbose output\n"
        "  -h, --help     Show this help\n\n"
        "Available RNGs:\n",
        prog,
        g_rand_gens[0].name
    );
    for (size_t i = 0; i < g_rand_gen_cnt; ++i) {
        fprintf(stderr, "  %s\n", g_rand_gens[i].name);
    }
}

static void print_bytes(const uint8_t* bytes, size_t n, enum cfx_str_format fmt) {
    if (fmt == CFX_STR_FMT_BINARY) {
#ifdef _WIN32
        _setmode(_fileno(stdout), _O_BINARY);
#endif
        fwrite(bytes, 1, n, stdout);
    } else if (fmt == CFX_STR_FMT_BASE64) {
        size_t b64_len = cfx_base64_enc_len(n);
        char* b64 = (char*)malloc(b64_len + 1);
        if (b64) {
            cfx_base64_encode(b64, &b64_len, bytes, n);
            b64[b64_len] = '\0';
            printf("%s\n", b64);
            free(b64);
        }
    } else {
        static const char hex_chars[] = "0123456789abcdef";
        char* hex = (char*)malloc(n * 2 + 2);
        if (hex) {
            for (size_t i = 0; i < n; ++i) {
                hex[i * 2]     = hex_chars[bytes[i] >> 4];
                hex[i * 2 + 1] = hex_chars[bytes[i] & 0x0F];
            }
            hex[n * 2] = '\n';
            hex[n * 2 + 1] = '\0';
            fputs(hex, stdout);
            free(hex);
        }
    }
}

int cfx_rand_run(int argc, char** argv) {
    uint64_t ns = cfx_time_ns();
    uint32_t seed = (uint32_t)ns;

    size_t n = 4;
    int verbose = 0;
    enum cfx_str_format fmt = CFX_STR_FMT_HEX;
    const char* prog = argv[0];

    const cfx_rand_desc_t* rand_gen = &g_rand_gens[0]; /* default */

    #define CHECK_NEXT_ARG(i, arg) do { \
        if (i + 1 >= argc) { \
            fprintf(stderr, "%s requires an argument\n\n", arg); \
            usage(prog); \
            return EXIT_FAILURE; \
        } \
    } while (0)

    for (int i = 1; i < argc; i++) {
        const char *arg = argv[i];

        if ((strcmp(arg, "--seed") == 0) || (strcmp(arg, "-s") == 0)) {
            CHECK_NEXT_ARG(i, "--seed");
            arg = argv[++i];
            char* end = NULL;
            seed = strtoull(arg, &end, 0);
            if (end == arg) {
                fprintf(stderr, "Invalid seed: %s\n\n", arg);
                usage(prog);
                return EXIT_FAILURE;
            }
        } else if (strcmp(arg, "--rng") == 0) {
            CHECK_NEXT_ARG(i, "--rng");
            const char *name = argv[++i];
            int found = 0;
            for (size_t j = 0; j < g_rand_gen_cnt; ++j) {
                if (strcmp(name, g_rand_gens[j].name) == 0) {
                    rand_gen = &g_rand_gens[j];
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
            return 0;
        } else if ((strcmp(arg, "-v") == 0) || (strcmp(arg, "--verbose") == 0)) {
            verbose = 1;
        } else if (strcmp(arg, "-x") == 0) {
            fmt = CFX_STR_FMT_HEX;
        } else if (strcmp(arg, "-b64") == 0) {
            fmt = CFX_STR_FMT_BASE64;
        } else if ((strcmp(arg, "-b") == 0) || (strcmp(arg, "--bin") == 0) || (strcmp(arg, "-bin") == 0)) {
            fmt = CFX_STR_FMT_BINARY;
        } else if (strcmp(arg, "-n") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "-n requires an argument\n\n");
                usage(prog);
                return EXIT_FAILURE;
            }
            char* end = NULL;
            errno = 0;
            n = strtoull(argv[++i], &end, 0);
            if (errno || end == argv[i] || *end != '\0') {
                fprintf(stderr, "Invalid num bytes: %s\n\n", argv[i]);
                usage(prog);
                return EXIT_FAILURE;
            }
        } else {
            char* end = NULL;
            errno = 0;
            n = strtoull(arg, &end, 0);
            if (errno || end == arg || *end != '\0') {
                fprintf(stderr, "Invalid num bytes: %s\n\n", arg);
                usage(prog);
                return EXIT_FAILURE;
            }
        }
    }

    if (verbose) printf("selected gen: %s: seed = 0x%08x\n.........\n", rand_gen->name, seed);

    rand_gen->seed(seed);
    uint8_t* bytes = (uint8_t*)malloc(n);
    if (bytes) {
        size_t i = 0;
        while (i < n) {
            uint32_t v = rand_gen->rng32();
            for (size_t vn = 0; vn < sizeof(v) && (i < n); vn++) {
                bytes[i] = v & 0xFF;
                v >>= 8;
                ++i;
            }
        }
    } else {
        return EXIT_FAILURE;
    }

    print_bytes(bytes, n, fmt);
    free(bytes);
    return EXIT_SUCCESS;
}
