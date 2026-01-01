/* cfx_ed25519.c - Ed25519 digital signature utility */

#include "cfx/ed25519.h"
#include "cfx/rand.h"
#include "cfx/base64.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cfx_cmd.h"
#include "misc.h"


static void usage(const char* prog) {
    fprintf(stderr,
        "Usage: %s <subcommand> [options]\n"
        "  Ed25519 digital signatures (RFC 8032)\n\n"
        "Subcommands:\n"
        "  keygen             Generate keypair (prints seed and public key)\n"
        "  public <seed>      Derive public key from 32-byte seed\n"
        "  sign <seed> <msg>  Sign a message (hex-encoded message)\n"
        "  verify <pk> <sig> <msg>  Verify signature\n\n"
        "Options:\n"
        "  -x          Output as hex (default)\n"
        "  -b64        Output as base64\n"
        "  -q          Quiet mode (output only, no labels)\n"
        "  -h, --help  Show this help\n\n"
        "Examples:\n"
        "  %s keygen                                Generate new keypair\n"
        "  %s public <64-char-hex-seed>             Derive public key\n"
        "  %s sign <seed> <hex-message>             Sign message\n"
        "  %s verify <pk> <sig> <hex-msg>           Verify signature\n",
        prog, prog, prog, prog, prog);
}

static int cmd_keygen(enum cfx_str_format fmt, int quiet) {
    uint8_t seed[32], pk[32], sk[64];

    cfx_srand_os();
    cfx_rand_bytes(seed, 32);

    cfx_ed25519_create_keypair(pk, sk, seed);

    if (quiet) {
        cfx_printf_output(seed, 32, fmt);
        printf("\n");
        cfx_printf_output(pk, 32, fmt);
        printf("\n");
    } else {
        printf("seed:   ");
        cfx_printf_output(seed, 32, fmt);
        printf("\npublic: ");
        cfx_printf_output(pk, 32, fmt);
        printf("\n");
    }

    return 0;
}

static int cmd_public(const char* seed_hex, enum cfx_str_format fmt, int quiet) {
    uint8_t seed[32], pk[32], sk[64];

    if (cfx_parse_hex(seed_hex, seed, 32) != 0) {
        fprintf(stderr, "error: seed must be 64 hex characters\n");
        return 1;
    }

    cfx_ed25519_create_keypair(pk, sk, seed);

    if (!quiet) printf("public: ");
    cfx_printf_output(pk, 32, fmt);
    printf("\n");

    return 0;
}

static int cmd_sign(const char* seed_hex, const char* msg_hex,
                    enum cfx_str_format fmt, int quiet) {
    uint8_t seed[32], pk[32], sk[64], sig[64];
    size_t msg_len;
    uint8_t* msg;

    if (cfx_parse_hex(seed_hex, seed, 32) != 0) {
        fprintf(stderr, "error: seed must be 64 hex characters\n");
        return 1;
    }

    msg_len = strlen(msg_hex) / 2;
    msg = (uint8_t*)malloc(msg_len);
    if (!msg) {
        fprintf(stderr, "error: out of memory\n");
        return 1;
    }

    if (cfx_parse_hex(msg_hex, msg, msg_len) != 0) {
        fprintf(stderr, "error: invalid message hex\n");
        free(msg);
        return 1;
    }

    cfx_ed25519_create_keypair(pk, sk, seed);
    cfx_ed25519_sign(sig, msg, msg_len, sk);

    if (!quiet) printf("signature: ");
    cfx_printf_output(sig, 64, fmt);
    printf("\n");

    free(msg);
    return 0;
}

static int cmd_verify(const char* pk_hex, const char* sig_hex, const char* msg_hex) {
    uint8_t pk[32], sig[64];
    size_t msg_len;
    uint8_t* msg;

    if (cfx_parse_hex(pk_hex, pk, 32) != 0) {
        fprintf(stderr, "error: public key must be 64 hex characters\n");
        return 1;
    }

    if (cfx_parse_hex(sig_hex, sig, 64) != 0) {
        fprintf(stderr, "error: signature must be 128 hex characters\n");
        return 1;
    }

    msg_len = strlen(msg_hex) / 2;
    msg = (uint8_t*)malloc(msg_len > 0 ? msg_len : 1);
    if (!msg) {
        fprintf(stderr, "error: out of memory\n");
        return 1;
    }

    if (msg_len > 0 && cfx_parse_hex(msg_hex, msg, msg_len) != 0) {
        fprintf(stderr, "error: invalid message hex\n");
        free(msg);
        return 1;
    }

    int ret = cfx_ed25519_verify(sig, msg, msg_len, pk);

    if (ret == 0) {
        printf("valid\n");
    } else {
        printf("INVALID\n");
    }

    free(msg);
    return ret == 0 ? 0 : 1;
}

int cfx_ed25519_run(int argc, char** argv) {
    const char* prog = argv[0];

    if (argc < 2) {
        usage(prog);
        return 1;
    }

    enum cfx_str_format fmt = CFX_STR_FMT_HEX;
    int quiet = 0;
    const char* subcmd = NULL;
    const char* arg1 = NULL;
    const char* arg2 = NULL;
    const char* arg3 = NULL;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-x") == 0) {
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
            return 1;
        } else if (subcmd == NULL) {
            subcmd = argv[i];
        } else if (arg1 == NULL) {
            arg1 = argv[i];
        } else if (arg2 == NULL) {
            arg2 = argv[i];
        } else if (arg3 == NULL) {
            arg3 = argv[i];
        } else {
            fprintf(stderr, "Too many arguments\n");
            usage(prog);
            return 1;
        }
    }

    if (subcmd == NULL) {
        usage(prog);
        return 1;
    }

    if (strcmp(subcmd, "keygen") == 0) {
        return cmd_keygen(fmt, quiet);
    } else if (strcmp(subcmd, "public") == 0) {
        if (arg1 == NULL) {
            fprintf(stderr, "error: 'public' requires seed argument\n");
            usage(prog);
            return 1;
        }
        return cmd_public(arg1, fmt, quiet);
    } else if (strcmp(subcmd, "sign") == 0) {
        if (arg1 == NULL || arg2 == NULL) {
            fprintf(stderr, "error: 'sign' requires seed and message arguments\n");
            usage(prog);
            return 1;
        }
        return cmd_sign(arg1, arg2, fmt, quiet);
    } else if (strcmp(subcmd, "verify") == 0) {
        if (arg1 == NULL || arg2 == NULL || arg3 == NULL) {
            fprintf(stderr, "error: 'verify' requires public key, signature, and message\n");
            usage(prog);
            return 1;
        }
        return cmd_verify(arg1, arg2, arg3);
    } else {
        fprintf(stderr, "Unknown subcommand: %s\n", subcmd);
        usage(prog);
        return 1;
    }
}
