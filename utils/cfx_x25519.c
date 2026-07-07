/* cfx_x25519.c - X25519 key exchange utility */

#include "cfx/x25519.h"
#include "cfx/rand.h"
#include "cfx/base64.h"
#include "cfx/memory.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cfx_cmd.h"
#include "cfx_utils_common.h"

/* Read a key from file (auto-detects hex/base64/raw binary) */
static int read_key_file(const char* path, uint8_t* out, size_t len) {
    enum cfx_str_format fmt = cfx_detect_file_format(path, NULL);
    if (fmt == CFX_STR_FMT_BINARY) {
        return cfx_read_file_bin(path, out, len);
    }
    int n = cfx_read_file_text(path, out, len, fmt);
    if (n < 0 || (size_t)n != len) {
        fprintf(stderr, "error: could not read %zu-byte key from %s\n", len, path);
        return -1;
    }
    return 0;
}

/* Load a 32-byte key from hex string, file, or stdin */
static int load_key32(const char* hex_or_null, const char* file_or_null, uint8_t out[32]) {
    if (file_or_null)
        return read_key_file(file_or_null, out, 32);
    if (hex_or_null) {
        if (cfx_parse_hex(hex_or_null, out, 32) != 0) {
            fprintf(stderr, "error: key must be 64 hex characters\n");
            return -1;
        }
        return 0;
    }
    /* fallback: read one line from stdin */
    char* line = cfx_read_line_stdin();
    if (!line || *line == '\0') {
        fprintf(stderr, "error: no key provided (hex arg, -k <file>, or pipe via stdin)\n");
        free(line);
        return -1;
    }
    int ret = cfx_parse_str_exact(line, out, 32);
    if (ret != 0)
        fprintf(stderr, "error: could not parse key from stdin (expected 64 hex chars or b64:...)\n");
    cfx_memzero_s(line, strlen(line));
    free(line);
    return ret;
}

static int write_key_file(const char* path, const uint8_t* data, size_t len,
                          enum cfx_str_format fmt) {
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
        for (size_t i = 0; i < len; i++)
            fprintf(f, "%02x", data[i]);
        fprintf(f, "\n");
    }

    fclose(f);
    return 0;
}

static void usage(const char* prog) {
    fprintf(stderr,
        "Usage: %s <subcommand> [options]\n"
        "  X25519 Diffie-Hellman key exchange (RFC 7748)\n\n"
        "Subcommands:\n"
        "  keygen               Generate keypair\n"
        "  public [<privkey>]   Derive public key from private key\n"
        "  shared [<privkey>] <pubkey>   Compute shared secret\n\n"
        "Options:\n"
        "  -k <file>   Read private key from file (32 bytes raw)\n"
        "  -o <file>   Write private key to file (keygen only)\n"
        "  -p <file>   Write public key to file (keygen only)\n"
        "  -x          Output as hex (default)\n"
        "  -b64        Output as base64\n"
        "  -bin        Output as raw binary (for -o/-p files)\n"
        "  -q          Quiet mode (keys only, no labels)\n"
        "  -h, --help  Show this help\n\n"
        "Keys can be provided as hex arguments, read from files (-k),\n"
        "or piped via stdin (one line, hex or b64:... format).\n\n"
        "Examples:\n"
        "  %s keygen                              Generate new keypair\n"
        "  %s keygen -o priv.key -p pub.key       Save keys to files\n"
        "  %s public <64-hex-privkey>             Derive from hex argument\n"
        "  %s public -k priv.key                  Derive from key file\n"
        "  echo <hex> | %s public                 Derive from stdin pipe\n"
        "  %s shared <privkey-hex> <pubkey-hex>   Shared secret from hex\n"
        "  %s shared -k priv.key <pubkey-hex>     Shared secret with key file\n",
        prog, prog, prog, prog, prog, prog, prog, prog);
}

static int cmd_keygen(enum cfx_str_format fmt, int quiet,
                      const char* priv_file, const char* pub_file) {
    uint8_t privkey[32], pubkey[32];

    cfx_srand_os();
    cfx_rand_bytes(privkey, 32);
    cfx_x25519_base(pubkey, privkey);

    if (priv_file) {
        if (write_key_file(priv_file, privkey, 32, fmt) != 0) {
            cfx_memzero_s(privkey, sizeof(privkey));
            return 1;
        }
        if (!quiet) printf("wrote private key to %s\n", priv_file);
    }

    if (pub_file) {
        if (write_key_file(pub_file, pubkey, 32, fmt) != 0) {
            cfx_memzero_s(privkey, sizeof(privkey));
            return 1;
        }
        if (!quiet) printf("wrote public key to %s\n", pub_file);
    }

    if (!priv_file || !pub_file) {
        if (quiet) {
            if (!priv_file) {
                cfx_printf_output(privkey, 32, fmt);
                printf("\n");
            }
            if (!pub_file) {
                cfx_printf_output(pubkey, 32, fmt);
                printf("\n");
            }
        } else {
            if (!priv_file) {
                printf("private: ");
                cfx_printf_output(privkey, 32, fmt);
                printf("\n");
            }
            if (!pub_file) {
                printf("public:  ");
                cfx_printf_output(pubkey, 32, fmt);
                printf("\n");
            }
        }
    }

    cfx_memzero_s(privkey, sizeof(privkey));
    return 0;
}

static int cmd_public(const char* priv_hex, const char* key_file,
                      enum cfx_str_format fmt, int quiet) {
    uint8_t privkey[32], pubkey[32];

    if (load_key32(priv_hex, key_file, privkey) != 0)
        return 1;

    cfx_x25519_base(pubkey, privkey);

    if (!quiet) printf("public: ");
    cfx_printf_output(pubkey, 32, fmt);
    printf("\n");

    cfx_memzero_s(privkey, sizeof(privkey));
    return 0;
}

static int cmd_shared(const char* priv_hex, const char* priv_file,
                      const char* pub_hex,
                      enum cfx_str_format fmt, int quiet) {
    uint8_t privkey[32], pubkey[32], shared[32];

    if (load_key32(priv_hex, priv_file, privkey) != 0)
        return 1;

    if (cfx_parse_hex(pub_hex, pubkey, 32) != 0) {
        cfx_memzero_s(privkey, sizeof(privkey));
        fprintf(stderr, "error: public key must be 64 hex characters\n");
        return 1;
    }

    int ret = cfx_x25519(shared, privkey, pubkey);
    if (ret != 0) {
        cfx_memzero_s(privkey, sizeof(privkey));
        fprintf(stderr, "error: invalid public key (result is zero)\n");
        return 1;
    }

    if (!quiet) printf("shared: ");
    cfx_printf_output(shared, 32, fmt);
    printf("\n");

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
    const char* key_file = NULL;   /* -k option */
    const char* out_priv = NULL;   /* -o option (keygen) */
    const char* out_pub = NULL;    /* -p option (keygen) */
    const char* arg1 = NULL;       /* positional args */
    const char* arg2 = NULL;

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
        } else if (strcmp(argv[i], "-o") == 0) {
            if (++i >= argc) {
                fprintf(stderr, "error: -o requires argument\n");
                return 1;
            }
            out_priv = argv[i];
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
        return cmd_keygen(fmt, quiet, out_priv, out_pub);
    } else if (strcmp(subcmd, "public") == 0) {
        return cmd_public(arg1, key_file, fmt, quiet);
    } else if (strcmp(subcmd, "shared") == 0) {
        const char* priv_hex = key_file ? NULL : arg1;
        const char* pub_hex = key_file ? arg1 : arg2;
        if (pub_hex == NULL) {
            fprintf(stderr, "error: 'shared' requires peer public key (hex argument)\n");
            usage(prog);
            return 1;
        }
        return cmd_shared(priv_hex, key_file, pub_hex, fmt, quiet);
    } else {
        fprintf(stderr, "Unknown subcommand: %s\n", subcmd);
        usage(prog);
        return 1;
    }
}
