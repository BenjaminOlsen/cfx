/* cfx_ed25519.c - Ed25519 digital signature utility */

#include "cfx/ed25519.h"
#include "cfx/rand.h"
#include "cfx/base64.h"
#include "cfx/memory.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cfx_cmd.h"
#include "common.h"

/* Read a key from file (auto-detects hex/base64/raw binary) */
static int read_key_file(const char* path, uint8_t* out, size_t len) {
    enum cfx_str_format fmt = cfx_detect_file_format(path, NULL);
    if (fmt == CFX_STR_FMT_BINARY)
        return cfx_read_file_bin(path, out, len);
    int n = cfx_read_file_text(path, out, len, fmt);
    if (n < 0 || (size_t)n != len) {
        fprintf(stderr, "error: could not read %zu-byte key from %s\n", len, path);
        return -1;
    }
    return 0;
}

/* Load a 32-byte key from hex string, file, or stdin */
static int load_key32(const char* hex_or_null, const char* file_or_null,
                      uint8_t out[32], int try_stdin) {
    if (file_or_null) {
        return read_key_file(file_or_null, out, 32);
    } else if (hex_or_null) {
        if (cfx_parse_hex(hex_or_null, out, 32) != 0) {
            fprintf(stderr, "error: key must be 64 hex characters\n");
            return -1;
        }
        return 0;
    }
    if (try_stdin) {
        char* line = cfx_read_line_stdin();
        if (line && *line) {
            int ret = cfx_parse_str_exact(line, out, 32);
            if (ret != 0)
                fprintf(stderr, "error: could not parse key from stdin\n");
            cfx_memzero_s(line, strlen(line));
            free(line);
            return ret;
        }
        free(line);
    }
    fprintf(stderr, "error: no key provided (hex arg, -k <file>, or pipe via stdin)\n");
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
    /* fallback: read from stdin */
    return cfx_read_all_file(stdin, out, out_len);
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
        "  -o <file>   Write private key (seed) to file (keygen only)\n"
        "  -p <file>   Write public key to file (keygen only)\n"
        "  -x          Output as hex (default)\n"
        "  -b64        Output as base64\n"
        "  -bin        Output as raw binary (for -o/-p files)\n"
        "  -q          Quiet mode (output only, no labels)\n"
        "  -h, --help  Show this help\n\n"
        "Keys can be provided as hex arguments, read from files (-k),\n"
        "or piped via stdin (one line, hex or b64:... format).\n"
        "Messages can be provided as hex, read from files (-m), or\n"
        "piped via stdin (raw bytes).\n\n"
        "Examples:\n"
        "  %s keygen                                Generate new keypair\n"
        "  %s keygen -o secret.key -p public.key    Save keys to files (hex)\n"
        "  %s keygen -o sec.bin -p pub.bin -bin     Save keys as raw binary\n"
        "  %s public <64-char-hex-seed>             Derive public key from hex\n"
        "  %s public -k seed.bin                    Derive public key from file\n"
        "  echo <hex> | %s public                   Derive public key from stdin\n"
        "  %s sign <seed-hex> <msg-hex>             Sign hex message with hex seed\n"
        "  %s sign -k seed.bin -m message.txt       Sign file with key file\n"
        "  cat msg.txt | %s sign -k seed.bin        Sign stdin with key file\n"
        "  %s verify <pk> <sig> <hex-msg>           Verify signature\n"
        "  %s verify -k pub.bin <sig> -m msg.txt    Verify with key/msg files\n"
        "  cat msg.txt | %s verify -k pk.bin <sig>  Verify message from stdin\n",
        prog, prog, prog, prog, prog, prog, prog, prog, prog, prog, prog, prog, prog);
}

static int write_key_file(const char* path, const uint8_t* data, size_t len, enum cfx_str_format fmt) {
    FILE* f = fopen(path, "wb");
    if (!f) {
        perror(path);
        return -1;
    }

    if (fmt == CFX_STR_FMT_BINARY) {
        if (fwrite(data, 1, len, f) != len) {
            fprintf(stderr, "error: failed to write %s\n", path);
            fclose(f);
            return -1;
        }
    } else if (fmt == CFX_STR_FMT_BASE64) {
        size_t b64_len = cfx_base64_enc_len(len);
        char* b64 = (char*)malloc(b64_len + 1);
        if (!b64) {
            fclose(f);
            return -1;
        }
        cfx_base64_encode(b64, &b64_len, data, len);
        b64[b64_len] = '\0';
        fprintf(f, "%s\n", b64);
        free(b64);
    } else {
        for (size_t i = 0; i < len; i++) {
            fprintf(f, "%02x", data[i]);
        }
        fprintf(f, "\n");
    }

    fclose(f);
    return 0;
}

static int cmd_keygen(enum cfx_str_format fmt, int quiet,
                      const char* seed_file, const char* pub_file) {
    uint8_t seed[32], pk[32], sk[64];

    cfx_srand_os();
    cfx_rand_bytes(seed, 32);

    cfx_ed25519_create_keypair(pk, sk, seed);

    /* Write to files if specified */
    if (seed_file) {
        if (write_key_file(seed_file, seed, 32, fmt) != 0) {
            cfx_memzero_s(seed, sizeof(seed));
            cfx_memzero_s(sk, sizeof(sk));
            return 1;
        }
        if (!quiet) printf("wrote seed to %s\n", seed_file);
    }

    if (pub_file) {
        if (write_key_file(pub_file, pk, 32, fmt) != 0) {
            cfx_memzero_s(seed, sizeof(seed));
            cfx_memzero_s(sk, sizeof(sk));
            return 1;
        }
        if (!quiet) printf("wrote public key to %s\n", pub_file);
    }

    /* Print to stdout if no files specified, or if neither file covers both */
    if (!seed_file || !pub_file) {
        if (quiet) {
            if (!seed_file) {
                cfx_printf_output(seed, 32, fmt);
                printf("\n");
            }
            if (!pub_file) {
                cfx_printf_output(pk, 32, fmt);
                printf("\n");
            }
        } else {
            if (!seed_file) {
                printf("seed:   ");
                cfx_printf_output(seed, 32, fmt);
                printf("\n");
            }
            if (!pub_file) {
                printf("public: ");
                cfx_printf_output(pk, 32, fmt);
                printf("\n");
            }
        }
    }

    cfx_memzero_s(seed, sizeof(seed));
    cfx_memzero_s(sk, sizeof(sk));
    return 0;
}

static int cmd_public(const char* seed_hex, const char* seed_file,
                      enum cfx_str_format fmt, int quiet) {
    uint8_t seed[32], pk[32], sk[64];

    if (load_key32(seed_hex, seed_file, seed, 1) != 0) {
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

    if (load_key32(seed_hex, seed_file, seed, 0) != 0) {
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

    if (load_key32(pk_hex, pk_file, pk, 0) != 0) {
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
    const char* out_seed = NULL;   /* -o option (keygen) */
    const char* out_pub = NULL;    /* -p option (keygen) */
    const char* arg1 = NULL;       /* positional args */
    const char* arg2 = NULL;
    const char* arg3 = NULL;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-x") == 0) {
            fmt = CFX_STR_FMT_HEX;
        } else if (strcmp(argv[i], "-b64") == 0) {
            fmt = CFX_STR_FMT_BASE64;
        } else if (strcmp(argv[i], "-bin") == 0) {
            fmt = CFX_STR_FMT_BINARY;
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
        } else if (strcmp(argv[i], "-o") == 0) {
            if (++i >= argc) {
                fprintf(stderr, "error: -o requires argument\n");
                return 1;
            }
            out_seed = argv[i];
        } else if (strcmp(argv[i], "-p") == 0) {
            if (++i >= argc) {
                fprintf(stderr, "error: -p requires argument\n");
                return 1;
            }
            out_pub = argv[i];
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
        return cmd_keygen(fmt, quiet, out_seed, out_pub);
    } else if (strcmp(subcmd, "public") == 0) {
        /* key from hex arg, -k file, or stdin */
        return cmd_public(arg1, key_file, fmt, quiet);
    } else if (strcmp(subcmd, "sign") == 0) {
        if (arg1 == NULL && key_file == NULL) {
            fprintf(stderr, "error: 'sign' requires seed (hex arg or -k file)\n");
            usage(prog);
            return 1;
        }
        /* arg1 is seed hex if no -k, arg2 is msg hex if no -m
         * message falls back to stdin if neither provided */
        const char* seed_hex = key_file ? NULL : arg1;
        const char* msg_hex = msg_file ? NULL : (key_file ? arg1 : arg2);
        return cmd_sign(seed_hex, key_file, msg_hex, msg_file, fmt, quiet);
    } else if (strcmp(subcmd, "verify") == 0) {
        /* verify needs: pk, sig, msg - sig is always hex for now
         * message falls back to stdin if neither hex arg nor -m file */
        const char* pk_hex = key_file ? NULL : arg1;
        const char* sig_hex = key_file ? arg1 : arg2;
        const char* msg_hex_arg = key_file ? arg2 : arg3;
        const char* msg_hex = msg_file ? NULL : msg_hex_arg;

        if ((pk_hex == NULL && key_file == NULL) || sig_hex == NULL) {
            fprintf(stderr, "error: 'verify' requires public key and signature\n");
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
