/* cfx_verify.c - Verify Ed25519 signatures */

#include "cfx/ed25519.h"
#include "cfx/sha512.h"
#include "cfx/base64.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "cfx_cmd.h"
#include "misc.h"

#ifdef _WIN32
#include <io.h>
#include <fcntl.h>
#endif

static void usage(const char* prog) {
    fprintf(stderr,
        "Usage: %s [options] <file>\n"
        "  Verify an Ed25519ph signature on a file (RFC 8032 pre-hashed mode).\n\n"
        "Options:\n"
        "  -k, --pubkey-file <file>   Read 32-byte public key from file\n"
        "  -K, --pubkey <hex>         Provide public key as 64 hex chars\n"
        "  -s, --sig-file <file>      Read signature from file\n"
        "  -S, --sig <hex>            Provide signature as 128 hex chars\n"
        "  -q                         Quiet mode (exit code only, no output)\n"
        "  -h, --help                 Show this help\n\n"
        "Exit codes: 0 = valid, 1 = invalid or error\n\n"
        "Examples:\n"
        "  %s --pubkey-file pub.key --sig-file doc.sig document.pdf\n"
        "  %s --pubkey <64-hex> --sig <128-hex> file.txt\n",
        prog, prog, prog);
}

/* Read a 32-byte key from file, auto-detecting format (raw binary or hex) */
static int read_key_auto(const char* path, uint8_t* out) {
    FILE* f = fopen(path, "rb");
    if (!f) {
        perror(path);
        return -1;
    }

    char buf[128];
    size_t n = fread(buf, 1, sizeof(buf), f);
    fclose(f);

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

/* read signature from file - handles hex or base64 text, or raw binary */
static int read_sig_file(const char* path, uint8_t sig[64]) {
    FILE* f = fopen(path, "rb");
    if (!f) {
        perror(path);
        return -1;
    }

    /* read up to 256 bytes to detect format */
    char buf[256];
    size_t n = fread(buf, 1, sizeof(buf) - 1, f);
    fclose(f);

    if (n == 0) {
        fprintf(stderr, "error: empty signature file\n");
        return -1;
    }

    /* raw 64-byte binary? */
    if (n == 64) {
        memcpy(sig, buf, 64);
        return 0;
    }

    buf[n] = '\0';

    /* strip trailing whitespace */
    while (n > 0 && isspace((unsigned char)buf[n-1])) {
        buf[--n] = '\0';
    }

    /* 128 hex chars? */
    if (n == 128) {
        if (cfx_parse_hex(buf, sig, 64) == 0) {
            return 0;
        }
    }

    /* try base64 (~88 chars for 64 bytes) */
    if (n >= 86 && n <= 90) {
        size_t out_len = 64;
        if (cfx_base64_decode(sig, &out_len, buf, n) == 0 && out_len == 64) {
            return 0;
        }
    }

    fprintf(stderr, "error: could not parse signature file (expected 64 raw bytes, 128 hex chars, or base64)\n");
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

int cfx_verify_run(int argc, char** argv) {
    if (argc < 2) {
        usage(argv[0]);
        return 1;
    }

    const char* keyfile = NULL;
    const char* keyhex = NULL;
    const char* sigfile = NULL;
    const char* sighex = NULL;
    const char* infile = NULL;
    int quiet = 0;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            usage(argv[0]);
            return 0;
        } else if (strcmp(argv[i], "-q") == 0) {
            quiet = 1;
        } else if (strcmp(argv[i], "-k") == 0 || strcmp(argv[i], "--pubkey-file") == 0) {
            if (++i >= argc) {
                fprintf(stderr, "error: %s requires argument\n", argv[i-1]);
                return 1;
            }
            keyfile = argv[i];
        } else if (strcmp(argv[i], "-K") == 0 || strcmp(argv[i], "--pubkey") == 0) {
            if (++i >= argc) {
                fprintf(stderr, "error: %s requires argument\n", argv[i-1]);
                return 1;
            }
            keyhex = argv[i];
        } else if (strcmp(argv[i], "-s") == 0 || strcmp(argv[i], "--sig-file") == 0) {
            if (++i >= argc) {
                fprintf(stderr, "error: %s requires argument\n", argv[i-1]);
                return 1;
            }
            sigfile = argv[i];
        } else if (strcmp(argv[i], "-S") == 0 || strcmp(argv[i], "--sig") == 0) {
            if (++i >= argc) {
                fprintf(stderr, "error: %s requires argument\n", argv[i-1]);
                return 1;
            }
            sighex = argv[i];
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

    if (!sigfile && !sighex) {
        fprintf(stderr, "error: must specify -s <sigfile> or -S <hex>\n");
        return 1;
    }

    /* read public key */
    uint8_t pk[32];
    if (keyhex) {
        if (cfx_parse_hex(keyhex, pk, 32) != 0) {
            fprintf(stderr, "error: public key must be 64 hex characters\n");
            return 1;
        }
    } else {
        if (read_key_auto(keyfile, pk) != 0) {
            return 1;
        }
    }

    /* read signature */
    uint8_t sig[64];
    if (sighex) {
        if (cfx_parse_hex(sighex, sig, 64) != 0) {
            fprintf(stderr, "error: signature must be 128 hex characters\n");
            return 1;
        }
    } else {
        if (read_sig_file(sigfile, sig) != 0) {
            return 1;
        }
    }

    /* pre-hash */
    uint8_t prehash[64];
    if (hash_file(infile, prehash) != 0) {
        return 1;
    }

    int ret = cfx_ed25519ph_verify(sig, prehash, pk);

    if (ret == 0) {
        if (!quiet) printf("OK: signature valid\n");
        return 0;
    } else {
        if (!quiet) printf("FAIL: signature invalid\n");
        return 1;
    }
}
