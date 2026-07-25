    
#include "bge.h"
#include "cfx_bge_internal.h"

static int armor_encode(const uint8_t *bin, size_t bin_len,
                        uint8_t **out, size_t *out_len) {
    size_t b64_len = 0;
    cfx_base64_encode(NULL, &b64_len, bin, bin_len);
    char *b64 = malloc(b64_len ? b64_len : 1);
    if (!b64) return -1;
    cfx_base64_encode(b64, &b64_len, bin, bin_len);

    size_t hlen = strlen(BGE_ARMOR_HEADER);
    size_t flen = strlen(BGE_ARMOR_FOOTER);
    size_t nlines = (b64_len + 75) / 76;
    size_t total = hlen + 1 + b64_len + nlines + flen + 1;
    uint8_t *buf = malloc(total + 1);
    if (!buf) {
        free(b64);
        return -1;
    }

    uint8_t *w = buf;
    memcpy(w, BGE_ARMOR_HEADER, hlen);
    w += hlen;
    *w++ = '\n';
    for (size_t i = 0; i < b64_len; i += 76) {
        size_t chunk = b64_len - i;
        if (chunk > 76) chunk = 76;
        memcpy(w, b64 + i, chunk);
        w += chunk;
        *w++ = '\n';
    }
    memcpy(w, BGE_ARMOR_FOOTER, flen);
    w += flen;
    *w++ = '\n';
    free(b64);

    *out = buf;
    *out_len = (size_t)(w - buf);
    return 0;
}

static int read_input(const char *path, uint8_t **data, size_t *len) {
    FILE *file = stdin;
    if (path) {
        file = fopen(path, "rb");
    }
    if (!file) {
        fprintf(stderr, "error: cannot open %s: %s\n", path, strerror(errno));
        return -1;
    }

    int rc = cfx_read_all_file(file, data, len);
    if (path) fclose(file);
    if (rc != 0) {
        fprintf(stderr, "error: cannot read input\n");
        return -1;
    }
    return 0;
}

static int write_output(const char *path, const uint8_t *data, size_t len) {
    FILE *file = stdout;
    if (path) {
        file = fopen(path, "wb");
        if (!file) {
            fprintf(stderr, "error: cannot open %s: %s\n", path, strerror(errno));
            return -1;
        }
    }

    int rc = len == 0 || fwrite(data, 1, len, file) == len ? 0 : -1;
    if (path && fclose(file) != 0) rc = -1;
    if (rc != 0) fprintf(stderr, "error: cannot write output\n");
    return rc;
}

static int parse_file_args(int argc, char **argv, int allow_armor,
                           const char **input, const char **output,
                           int *armor) {
    *input = NULL;
    *output = NULL;
    *armor = 0;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-e") == 0 || strcmp(argv[i], "--encrypt") == 0 ||
            strcmp(argv[i], "-d") == 0 || strcmp(argv[i], "--decrypt") == 0) {
            continue;
        }
        if (allow_armor &&
            (strcmp(argv[i], "-a") == 0 || strcmp(argv[i], "--armor") == 0)) {
            *armor = 1;
        } else if (strcmp(argv[i], "-o") == 0 ||
                   strcmp(argv[i], "--output") == 0) {
            if (++i >= argc) {
                fprintf(stderr, "error: -o requires a path\n");
                return -1;
            }
            *output = argv[i];
        } else if (strcmp(argv[i], "-i") == 0 ||
                   strcmp(argv[i], "--input") == 0) {
            if (++i >= argc) {
                fprintf(stderr, "error: -i requires a path\n");
                return -1;
            }
            *input = argv[i];
        } else if (argv[i][0] == '-') {
            fprintf(stderr, "error: unknown option: %s\n", argv[i]);
            return -1;
        } else if (!*input) {
            *input = argv[i];
        } else {
            fprintf(stderr, "error: unexpected argument: %s\n", argv[i]);
            return -1;
        }
    }
    return 0;
}

int bge_encrypt_file(int argc, char **argv) {
    const char *input_path;
    const char *output_path;
    int armor;
    if (parse_file_args(argc, argv, 1, &input_path, &output_path, &armor) != 0)
        return 1;

    uint8_t *input = NULL;
    size_t input_len = 0;
    if (read_input(input_path, &input, &input_len) != 0) return 1;

    char passphrase[256] = {0};
    int passphrase_len = cfx_prompt_passphrase(passphrase, sizeof(passphrase));
    if (passphrase_len < 0) {
        cfx_bge_free(input, input_len);
        return 1;
    }

    uint8_t *encrypted = NULL;
    size_t encrypted_len = 0;
    int rc = cfx_bge_encrypt(input, input_len,
                             (const uint8_t *)passphrase,
                             (size_t)passphrase_len,
                             &encrypted, &encrypted_len);
    cfx_memzero_s(passphrase, sizeof(passphrase));
    cfx_bge_free(input, input_len);
    if (rc != 0) {
        fprintf(stderr, "error: encryption failed\n");
        return 1;
    }

    uint8_t *result = encrypted;
    size_t result_len = encrypted_len;
    if (armor) {
        if (armor_encode(encrypted, encrypted_len, &result, &result_len) != 0) {
            fprintf(stderr, "error: armor encoding failed\n");
            cfx_bge_free(encrypted, encrypted_len);
            return 1;
        }
    }

    rc = write_output(output_path, result, result_len);
    if (armor) {
        cfx_memzero_s(result, result_len);
        free(result);
        cfx_bge_free(encrypted, encrypted_len);
    } else {
        cfx_bge_free(encrypted, encrypted_len);
    }
    return rc == 0 ? 0 : 1;
}

int bge_decrypt_file(int argc, char **argv) {
    const char *input_path;
    const char *output_path;
    int unused_armor;
    if (parse_file_args(argc, argv, 0, &input_path, &output_path,
                        &unused_armor) != 0)
        return 1;

    uint8_t *input = NULL;
    size_t input_len = 0;
    if (read_input(input_path, &input, &input_len) != 0) return 1;

    char passphrase[256] = {0};
    int passphrase_len = bge_read_passphrase(
        "Enter passphrase: ", passphrase, sizeof(passphrase));
    if (passphrase_len <= 0) {
        cfx_bge_free(input, input_len);
        return 1;
    }

    uint8_t *plaintext = NULL;
    size_t plaintext_len = 0;
    int rc = cfx_bge_decrypt(input, input_len,
                             (const uint8_t *)passphrase,
                             (size_t)passphrase_len,
                             &plaintext, &plaintext_len);
    cfx_memzero_s(passphrase, sizeof(passphrase));
    cfx_bge_free(input, input_len);
    if (rc != 0) {
        fprintf(stderr, rc == -3
            ? "error: authentication failed\n"
            : "error: invalid BGE input\n");
        return 1;
    }

    rc = write_output(output_path, plaintext, plaintext_len);
    cfx_bge_free(plaintext, plaintext_len);
    return rc == 0 ? 0 : 1;
}

static void usage(const char *prog) {
    printf("Usage:\n");
    printf("  %s -e [-a] [-i INPUT] [-o OUTPUT]    Encrypt a file\n", prog);
    printf("  %s -d [-i INPUT] [-o OUTPUT]         Decrypt a file\n\n", prog);
    printf("File encryption using Argon2id key derivation and\n");
    printf("XChaCha20-Poly1305 AEAD.\n\n");
    printf("Options:\n");
    printf("  -e, --encrypt    Encrypt input to output\n");
    printf("  -d, --decrypt    Decrypt input to output\n");
    printf("  -i, --input      Read from file (default: stdin)\n");
    printf("  -o, --output     Write to file (default: stdout)\n");
    printf("  -a, --armor      PEM-encoded base64 (encrypt only; decrypt auto-detects)\n");
    printf("  -p, --passphrase <pw>  Supply passphrase on command line (default: prompt)\n");
    printf("\nExamples:\n");
    printf("  %s -e -i doc.txt -o secret.bge        Encrypt a file\n", prog);
    printf("  %s -e -a -i doc.txt -o secret.asc     Encrypt to armored output\n", prog);
    printf("  %s -d -i secret.bge -o doc.txt        Decrypt a file\n", prog);
    printf("  echo data | %s -e | %s -d             Roundtrip via pipe\n", prog, prog);
}

int cfx_bge_run(int argc, char **argv) {
    g_passphrase_arg = NULL;

    /* extract -p / --passphrase before dispatching */
    for (int i = 1; i < argc; i++) {
        if ((strcmp(argv[i], "-p") == 0 ||
             strcmp(argv[i], "--passphrase") == 0) && i + 1 < argc) {
            g_passphrase_arg = argv[i + 1];
            for (int j = i; j + 2 < argc; j++)
                argv[j] = argv[j + 2];
            argc -= 2;
            i--;
            break;
        }
    }

    if (argc < 2) {
        usage(argv[0]);
        return 1;
    }

    const char *cmd = argv[1];

    if (strcmp(cmd, "-h") == 0 || strcmp(cmd, "--help") == 0 ||
        strcmp(cmd, "help") == 0) {
        usage(argv[0]);
        return 0;
    }

    if (strcmp(cmd, "--version") == 0 || strcmp(cmd, "-V") == 0) {
        printf("%s\n", BGE_VERSION_STR);
        return 0;
    }

    /* scan for -e/-d anywhere in argv */
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-e") == 0 || strcmp(argv[i], "--encrypt") == 0) {
            return bge_encrypt_file(argc, argv);
        }
        if (strcmp(argv[i], "-d") == 0 || strcmp(argv[i], "--decrypt") == 0) {
            return bge_decrypt_file(argc, argv);
        }
    }

    fprintf(stderr, "Unknown command: %s\n", cmd);
    usage(argv[0]);
    return 1;
}
