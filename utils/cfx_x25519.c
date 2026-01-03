/* cfx_x25519.c - X25519 key exchange utility */

#include "cfx/x25519.h"
#include "cfx/rand.h"
#include "cfx/base64.h"
#include "cfx/memory.h"
#include "misc.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cfx_cmd.h"
#include "misc.h"


static void usage(const char* prog) {
    fprintf(stderr,
        "Usage: %s <subcommand> [options]\n"
        "  X25519 Diffie-Hellman key exchange (RFC 7748)\n\n"
        "Subcommands:\n"
        "  keygen             Generate keypair (prints private and public keys)\n"
        "  public <privkey>   Derive public key from private key\n"
        "  shared <privkey> <pubkey>   Compute shared secret\n\n"
        "Options:\n"
        "  -x          Output as hex (default)\n"
        "  -b64        Output as base64\n"
        "  -q          Quiet mode (keys only, no labels)\n"
        "  -h, --help  Show this help\n\n"
        "Examples:\n"
        "  %s keygen                              Generate new keypair\n"
        "  %s public <64-char-hex-privkey>        Derive public key\n"
        "  %s shared <privkey> <pubkey>           Compute DH shared secret\n",
        prog, prog, prog, prog);
}

static int cmd_keygen(enum cfx_str_format fmt, int quiet) {
    uint8_t privkey[32], pubkey[32];

    /* generate random private key from OS entropy */
    cfx_srand_os();
    cfx_rand_bytes(privkey, 32);

    /* derive public key */
    cfx_x25519_base(pubkey, privkey);

    if (quiet) {
        cfx_printf_output(privkey, 32, fmt);
        printf("\n");
        cfx_printf_output(pubkey, 32, fmt);
        printf("\n");
    } else {
        printf("private: ");
        cfx_printf_output(privkey, 32, fmt);
        printf("\npublic:  ");
        cfx_printf_output(pubkey, 32, fmt);
        printf("\n");
    }

    /* zero sensitive data */
    cfx_memzero_s(privkey, sizeof(privkey));
    return 0;
}

static int cmd_public(const char* privkey_hex, enum cfx_str_format fmt, int quiet) {
    uint8_t privkey[32], pubkey[32];

    if (cfx_parse_hex(privkey_hex, privkey, 32) != 0) {
        fprintf(stderr, "error: private key must be 64 hex characters\n");
        return 1;
    }

    cfx_x25519_base(pubkey, privkey);

    if (!quiet) printf("public: ");
    cfx_printf_output(pubkey, 32, fmt);
    printf("\n");

    /* zero sensitive data */
    cfx_memzero_s(privkey, sizeof(privkey));
    return 0;
}

static int cmd_shared(const char* privkey_hex, const char* pubkey_hex,
                      enum cfx_str_format fmt, int quiet) {
    uint8_t privkey[32], pubkey[32], shared[32];
    int ret;

    if (cfx_parse_hex(privkey_hex, privkey, 32) != 0) {
        fprintf(stderr, "error: private key must be 64 hex characters\n");
        return 1;
    }

    if (cfx_parse_hex(pubkey_hex, pubkey, 32) != 0) {
        cfx_memzero_s(privkey, sizeof(privkey));
        fprintf(stderr, "error: public key must be 64 hex characters\n");
        return 1;
    }

    ret = cfx_x25519(shared, privkey, pubkey);
    if (ret != 0) {
        cfx_memzero_s(privkey, sizeof(privkey));
        fprintf(stderr, "error: invalid public key (result is zero)\n");
        return 1;
    }

    if (!quiet) printf("shared: ");
    cfx_printf_output(shared, 32, fmt);
    printf("\n");

    /* zero sensitive data */
    cfx_memzero_s(privkey, sizeof(privkey));
    cfx_memzero_s(shared, sizeof(shared));
    return 0;
}

int cfx_x25519_run(int argc, char** argv) {
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
            fprintf(stderr, "error: 'public' requires private key argument\n");
            usage(prog);
            return 1;
        }
        return cmd_public(arg1, fmt, quiet);
    } else if (strcmp(subcmd, "shared") == 0) {
        if (arg1 == NULL || arg2 == NULL) {
            fprintf(stderr, "error: 'shared' requires private key and public key arguments\n");
            usage(prog);
            return 1;
        }
        return cmd_shared(arg1, arg2, fmt, quiet);
    } else {
        fprintf(stderr, "Unknown subcommand: %s\n", subcmd);
        usage(prog);
        return 1;
    }
}
