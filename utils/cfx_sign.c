/* cfx_sign.c - Sign files with Ed25519 */

#include "cfx/ed25519.h"
#include "cfx/sha512.h"
#include "cfx/base64.h"
#include "cfx/memory.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cfx_cmd.h"
#include "misc.h"

#ifdef _WIN32
#include <io.h>
#include <fcntl.h>
#endif

static void usage(const char* prog) {
    fprintf(stderr,
        "Usage: %s [options] <file>\n"
        "  Sign a file with Ed25519ph (RFC 8032 pre-hashed mode).\n\n"
        "Options:\n"
        "  -k, --seed-file <file>   Read 32-byte seed from file (- for stdin)\n"
        "  -K, --seed <hex>         Provide seed directly as 64 hex chars\n"
        "  -o <sigfile>             Write signature to file (default: stdout)\n"
        "  -x                       Output as hex (default)\n"
        "  -b64                     Output as base64\n"
        "  -h, --help               Show this help\n\n"
        "The seed is your 32-byte private key. The signature is 64 bytes.\n\n"
        "Examples:\n"
        "  %s --seed-file secret.key document.pdf    Sign with key from file\n"
        "  %s --seed <64-hex-chars> document.pdf     Sign with hex seed\n"
        "  %s -k secret.key -o doc.sig doc.pdf       Write signature to file\n",
        prog, prog, prog, prog);
}

/* Read a 32-byte key from file, auto-detecting format (raw binary or hex) */
static int read_key_auto(const char* path, uint8_t* out) {
    FILE* f;
    if (strcmp(path, "-") == 0) {
#ifdef _WIN32
        _setmode(_fileno(stdin), _O_BINARY);
#endif
        f = stdin;
    } else {
        f = fopen(path, "rb");
        if (!f) {
            perror(path);
            return -1;
        }
    }

    char buf[128];
    size_t n = fread(buf, 1, sizeof(buf), f);
    if (f != stdin) fclose(f);

    /* Raw binary: exactly 32 bytes */
    if (n == 32) {
        memcpy(out, buf, 32);
        return 0;
    }

    /* Hex: 64 chars, optionally with trailing newline */
    if (n >= 64 && n <= 66) {
        /* Strip trailing whitespace */
        while (n > 0 && (buf[n-1] == '\n' || buf[n-1] == '\r' || buf[n-1] == ' ')) {
            n--;
        }
        if (n == 64) {
            buf[64] = '\0';
            if (cfx_parse_hex(buf, out, 32) == 0) {
                return 0;
            }
        }
    }

    fprintf(stderr, "error: key file must be 32 raw bytes or 64 hex chars\n");
    return -1;
}

static int hash_file(const char* path, uint8_t* hash) {
    FILE* f = fopen(path, "rb");
    if (!f) {
        perror(path);
        return -1;
    }

    cfx_sha512_ctx_t ctx;
    cfx_sha512_init(&ctx);

    uint8_t buf[262144];  /* 256KB for better throughput on large files */
    size_t n;
    while ((n = fread(buf, 1, sizeof(buf), f)) > 0) {
        cfx_sha512_update(&ctx, buf, n);
    }

    if (ferror(f)) {
        fprintf(stderr, "error reading %s\n", path);
        fclose(f);
        return -1;
    }

    fclose(f);
    cfx_sha512_final(&ctx, hash);
    return 0;
}

int cfx_sign_run(int argc, char** argv) {
    if (argc < 2) {
        usage(argv[0]);
        return 1;
    }

    enum cfx_str_format fmt = CFX_STR_FMT_HEX;
    const char* keyfile = NULL;
    const char* keyhex = NULL;
    const char* outfile = NULL;
    const char* infile = NULL;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            usage(argv[0]);
            return 0;
        } else if (strcmp(argv[i], "-x") == 0) {
            fmt = CFX_STR_FMT_HEX;
        } else if (strcmp(argv[i], "-b64") == 0) {
            fmt = CFX_STR_FMT_BASE64;
        } else if (strcmp(argv[i], "-k") == 0 || strcmp(argv[i], "--seed-file") == 0) {
            if (++i >= argc) {
                fprintf(stderr, "error: %s requires argument\n", argv[i-1]);
                return 1;
            }
            keyfile = argv[i];
        } else if (strcmp(argv[i], "-K") == 0 || strcmp(argv[i], "--seed") == 0) {
            if (++i >= argc) {
                fprintf(stderr, "error: %s requires argument\n", argv[i-1]);
                return 1;
            }
            keyhex = argv[i];
        } else if (strcmp(argv[i], "-o") == 0) {
            if (++i >= argc) {
                fprintf(stderr, "error: -o requires argument\n");
                return 1;
            }
            outfile = argv[i];
        } else if (argv[i][0] == '-') {
            fprintf(stderr, "unknown option: %s\n", argv[i]);
            usage(argv[0]);
            return 1;
        } else {
            if (infile) {
                fprintf(stderr, "error: multiple input files not supported\n");
                return 1;
            }
            infile = argv[i];
        }
    }

    if (!infile) {
        fprintf(stderr, "error: no input file specified\n");
        usage(argv[0]);
        return 1;
    }

    if (!keyfile && !keyhex) {
        fprintf(stderr, "error: must specify -k <keyfile> or -K <hex>\n");
        return 1;
    }

    uint8_t seed[32];
    if (keyhex) {
        if (cfx_parse_hex(keyhex, seed, 32) != 0) {
            fprintf(stderr, "error: seed must be 64 hex characters\n");
            return 1;
        }
    } else {
        if (read_key_auto(keyfile, seed) != 0) {
            return 1;
        }
    }

    uint8_t pk[32], sk[64];
    cfx_ed25519_create_keypair(pk, sk, seed);

    /* Ed25519ph pre-hash */
    uint8_t prehash[64];
    if (hash_file(infile, prehash) != 0) {
        return 1;
    }

    /* sign using Ed25519ph (RFC 8032 pre-hashed mode) */
    uint8_t sig[64];
    cfx_ed25519ph_sign(sig, prehash, sk);

    FILE* out = stdout;
    if (outfile) {
        out = fopen(outfile, "wb");
        if (!out) {
            perror(outfile);
            return 1;
        }
    }

    if (fmt == CFX_STR_FMT_BASE64) {
        size_t b64_len = cfx_base64_enc_len(64);
        char b64[128];
        cfx_base64_encode(b64, &b64_len, sig, 64);
        b64[b64_len] = '\0';
        fprintf(out, "%s\n", b64);
    } else {
        for (int i = 0; i < 64; i++) {
            fprintf(out, "%02x", sig[i]);
        }
        fprintf(out, "\n");
    }

    if (outfile) fclose(out);

    cfx_memzero_s(seed, sizeof(seed));
    cfx_memzero_s(sk, sizeof(sk));

    return 0;
}
