/* cfx_ed25519.c - Ed25519 digital signature utility */

#include "cfx/ed25519.h"
#include "cfx/rand.h"
#include "cfx/base64.h"
#include "cfx/memory.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cfx_cmd.h"
#include "misc.h"

/* Read exactly 'len' bytes from file, return 0 on success */
static int read_key_file(const char* path, uint8_t* out, size_t len) {
    FILE* f = fopen(path, "rb");
    if (!f) {
        perror(path);
        return -1;
    }
    size_t n = fread(out, 1, len, f);
    fclose(f);
    if (n != len) {
        fprintf(stderr, "error: expected %zu bytes, got %zu\n", len, n);
        return -1;
    }
    return 0;
}

/* Load a 32-byte key from either hex string or file */
static int load_key32(const char* hex_or_null, const char* file_or_null, uint8_t out[32]) {
    if (file_or_null) {
        return read_key_file(file_or_null, out, 32);
    } else if (hex_or_null) {
        if (cfx_parse_hex(hex_or_null, out, 32) != 0) {
            fprintf(stderr, "error: key must be 64 hex characters\n");
            return -1;
        }
        return 0;
    }
    fprintf(stderr, "error: no key provided\n");
    return -1;
}

/* Load message from either hex string or file */
static int load_message(const char* hex_or_null, const char* file_or_null,
                        uint8_t** out, size_t* out_len) {
    if (file_or_null) {
        FILE* f = fopen(file_or_null, "rb");
        if (!f) {
            perror(file_or_null);
            return -1;
        }
        int ret = cfx_read_all_file(f, out, out_len);
        fclose(f);
        return ret;
    } else if (hex_or_null) {
        size_t len = strlen(hex_or_null) / 2;
        uint8_t* msg = (uint8_t*)malloc(len > 0 ? len : 1);
        if (!msg) {
            fprintf(stderr, "error: out of memory\n");
            return -1;
        }
        if (len > 0 && cfx_parse_hex(hex_or_null, msg, len) != 0) {
            fprintf(stderr, "error: invalid message hex\n");
            free(msg);
            return -1;
        }
        *out = msg;
        *out_len = len;
        return 0;
    }
    fprintf(stderr, "error: no message provided\n");
    return -1;
}


static void usage(const char* prog) {
    fprintf(stderr,
        "Usage: %s <subcommand> [options]\n"
        "  Ed25519 digital signatures (RFC 8032)\n\n"
        "Subcommands:\n"
        "  keygen                     Generate keypair (prints seed and public key)\n"
        "  public <seed>              Derive public key from 32-byte seed\n"
        "  sign <seed> <msg>          Sign a message\n"
        "  verify <pk> <sig> <msg>    Verify signature\n\n"
        "Options:\n"
        "  -k <file>   Read seed/key from file (32 bytes raw)\n"
        "  -m <file>   Read message from file (instead of hex string)\n"
        "  -x          Output as hex (default)\n"
        "  -b64        Output as base64\n"
        "  -q          Quiet mode (output only, no labels)\n"
        "  -h, --help  Show this help\n\n"
        "Examples:\n"
        "  %s keygen                                Generate new keypair\n"
        "  %s keygen -q > mykey.txt                 Save keypair to file\n"
        "  %s public <64-char-hex-seed>             Derive public key from hex\n"
        "  %s public -k seed.bin                    Derive public key from file\n"
        "  %s sign <seed-hex> <msg-hex>             Sign hex message with hex seed\n"
        "  %s sign -k seed.bin -m message.txt       Sign file with key file\n"
        "  %s verify <pk> <sig> <hex-msg>           Verify signature\n"
        "  %s verify -k pub.bin <sig> -m msg.txt    Verify with key/msg files\n",
        prog, prog, prog, prog, prog, prog, prog, prog, prog);
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

static int cmd_public(const char* seed_hex, const char* seed_file,
                      enum cfx_str_format fmt, int quiet) {
    uint8_t seed[32], pk[32], sk[64];

    if (load_key32(seed_hex, seed_file, seed) != 0) {
        return 1;
    }

    cfx_ed25519_create_keypair(pk, sk, seed);

    if (!quiet) printf("public: ");
    cfx_printf_output(pk, 32, fmt);
    printf("\n");

    cfx_memzero_s(seed, sizeof(seed));
    cfx_memzero_s(sk, sizeof(sk));
    return 0;
}

static int cmd_sign(const char* seed_hex, const char* seed_file,
                    const char* msg_hex, const char* msg_file,
                    enum cfx_str_format fmt, int quiet) {
    uint8_t seed[32], pk[32], sk[64], sig[64];
    size_t msg_len;
    uint8_t* msg = NULL;

    if (load_key32(seed_hex, seed_file, seed) != 0) {
        return 1;
    }

    if (load_message(msg_hex, msg_file, &msg, &msg_len) != 0) {
        cfx_memzero_s(seed, sizeof(seed));
        return 1;
    }

    cfx_ed25519_create_keypair(pk, sk, seed);
    cfx_ed25519_sign(sig, msg, msg_len, sk);

    if (!quiet) printf("signature: ");
    cfx_printf_output(sig, 64, fmt);
    printf("\n");

    free(msg);
    cfx_memzero_s(seed, sizeof(seed));
    cfx_memzero_s(sk, sizeof(sk));
    return 0;
}

static int cmd_verify(const char* pk_hex, const char* pk_file,
                      const char* sig_hex,
                      const char* msg_hex, const char* msg_file) {
    uint8_t pk[32], sig[64];
    size_t msg_len;
    uint8_t* msg = NULL;

    if (load_key32(pk_hex, pk_file, pk) != 0) {
        return 1;
    }

    if (cfx_parse_hex(sig_hex, sig, 64) != 0) {
        fprintf(stderr, "error: signature must be 128 hex characters\n");
        return 1;
    }

    if (load_message(msg_hex, msg_file, &msg, &msg_len) != 0) {
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
    const char* key_file = NULL;   /* -k option */
    const char* msg_file = NULL;   /* -m option */
    const char* arg1 = NULL;       /* positional args */
    const char* arg2 = NULL;
    const char* arg3 = NULL;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-x") == 0) {
            fmt = CFX_STR_FMT_HEX;
        } else if (strcmp(argv[i], "-b64") == 0) {
            fmt = CFX_STR_FMT_BASE64;
        } else if (strcmp(argv[i], "-q") == 0) {
            quiet = 1;
        } else if (strcmp(argv[i], "-k") == 0) {
            if (++i >= argc) {
                fprintf(stderr, "error: -k requires argument\n");
                return 1;
            }
            key_file = argv[i];
        } else if (strcmp(argv[i], "-m") == 0) {
            if (++i >= argc) {
                fprintf(stderr, "error: -m requires argument\n");
                return 1;
            }
            msg_file = argv[i];
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
        if (arg1 == NULL && key_file == NULL) {
            fprintf(stderr, "error: 'public' requires seed (hex arg or -k file)\n");
            usage(prog);
            return 1;
        }
        return cmd_public(arg1, key_file, fmt, quiet);
    } else if (strcmp(subcmd, "sign") == 0) {
        if ((arg1 == NULL && key_file == NULL) || (arg2 == NULL && msg_file == NULL)) {
            fprintf(stderr, "error: 'sign' requires seed and message (hex args or -k/-m files)\n");
            usage(prog);
            return 1;
        }
        /* arg1 is seed hex if no -k, arg2 is msg hex if no -m */
        const char* seed_hex = key_file ? NULL : arg1;
        const char* msg_hex = msg_file ? NULL : (key_file ? arg1 : arg2);
        return cmd_sign(seed_hex, key_file, msg_hex, msg_file, fmt, quiet);
    } else if (strcmp(subcmd, "verify") == 0) {
        /* verify needs: pk, sig, msg - sig is always hex for now */
        const char* pk_hex = key_file ? NULL : arg1;
        const char* sig_hex = key_file ? arg1 : arg2;
        const char* msg_hex_arg = key_file ? arg2 : arg3;
        const char* msg_hex = msg_file ? NULL : msg_hex_arg;

        if ((pk_hex == NULL && key_file == NULL) || sig_hex == NULL ||
            (msg_hex == NULL && msg_file == NULL)) {
            fprintf(stderr, "error: 'verify' requires public key, signature, and message\n");
            usage(prog);
            return 1;
        }
        return cmd_verify(pk_hex, key_file, sig_hex, msg_hex, msg_file);
    } else {
        fprintf(stderr, "Unknown subcommand: %s\n", subcmd);
        usage(prog);
        return 1;
    }
}
