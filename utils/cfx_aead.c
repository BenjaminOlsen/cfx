
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

#include "cfx/aead_chacha20_poly1305.h"
#include "cfx/aead_chacha20_poly1271.h"
#include "cfx/memory.h"
#include "cfx/rand.h"
#include "cfx_cmd.h"
#include "common.h"

static void usage(const char* prog) {
        printf("Usage: %s [options] <file>\n", prog);
        printf("       %s [options] -           (read stdin)\n\n", prog);
        printf("       %s [options] -s <plaintext string>\n\n", prog);
        printf("Authenticated encryption with ChaCha20-Poly1305 (RFC 8439)\n\n");
        printf("Mode:\n");
        printf("  -e, --encrypt     Encrypt (default)\n");
        printf("  -d, --decrypt     Decrypt and verify tag\n\n");
        printf("Algorithm:\n");
        printf("  -C, --chacha      Use ChaCha20-Poly1305 (12-byte nonce)\n");
        printf("  -X, --xchacha     Use XChaCha20-Poly1305 (24-byte nonce)\n");
        printf("  -P, --p1271       Use ChaCha20-Poly1271 (12-byte nonce)\n");
        printf("  -XP, --xp1271     Use XChaCha20-Poly1271 (24-byte nonce)\n\n");
        printf("Key (required):\n");
        printf("  -k <hex>          Key as 32 byte (ascii or hex or base64) string\n");
        printf("  -K <file>         Read key from file\n\n");
        printf("Nonce:\n");
        printf("  -n <hex>          Nonce as 12 bytes (or 24 with -X)\n");
        printf("  --random-nonce    Generate random nonce, prepend to output\n");
        printf("  --nonce-prefix    Encrypt: prepend nonce to output\n");
        printf("                    Decrypt: read nonce from input\n\n");
        printf("Associated data:\n");
        printf("  -a <string>       AAD as string\n");
        printf("  -A <file>         Read AAD from file\n\n");
        printf("Output:\n");
        printf("  -o <file>         Write to file (default: stdout)\n");
        printf("  -x                Hex output (default for encrypt)\n");
        printf("  -b64              Base64 output\n");
        printf("  -r, --raw         Raw binary output\n");
        printf("  -v, --verbose     Verbose printing\n\n");
        printf("Examples:\n");
        printf("  %s -e -k $(cat key.hex) --random-nonce secret.txt\n", prog);
        printf("  %s -d -k $KEY --nonce-prefix encrypted.bin\n", prog);
        printf("  %s -X --random-nonce -k $KEY secret.txt  (XChaCha20)\n", prog);
    }

typedef enum {
    CHACHA20_POLY1305,
    XCHACHA20_POLY1305,
    CHACHA20_POLY1271,
    XCHACHA20_POLY1271
} algo_type_t;

static int try_reading_file_into(const char* fname, uint8_t* out, size_t outlen, enum cfx_str_format fmt) {
    if (fmt == CFX_STR_FMT_BINARY) {
        if (cfx_read_file_bin(fname, out, outlen) != 0) {
            cfx_memzero_s(out, outlen);
            return -1;
        }
    } else {
        int bytes_read = cfx_read_file_text(fname, out, outlen, CFX_STR_FMT_AUTO);
        if (bytes_read < 0 || (size_t)bytes_read != outlen) {
            cfx_memzero_s(out, outlen);
            return -1;
        }
    }
    return 0;
}


int cfx_aead_run(int argc, char** argv) {

    if (argc < 2) {
        usage(argv[0]);
        return 0;
    }

    int ret = 0;
    int encrypt = 1;  /* == 0 -> decrypt */
    int verbose = 0;
    int nonce_prefix = 0;
    int use_xchacha = 0;
    const char* keystr = NULL;
    const char* noncestr = NULL;
    const char* adstr = NULL;
    enum cfx_str_format outfmt = CFX_STR_FMT_HEX;
    uint8_t key[32] = {0};
    uint8_t nonce[24] = {0};  /* 12 for ChaCha20, 24 for XChaCha20 */
    uint8_t* aad = NULL;
    size_t aad_len = 0;
    uint8_t* pt_in = NULL;  /* owns the pt */
    size_t pt_in_len = 0;
    algo_type_t algo_type = CHACHA20_POLY1305;

    #define CHECK_NEXT_ARG() \
    do { \
        if (argi + 1 >= argc) { \
            usage(argv[0]); \
            return EXIT_FAILURE; \
        } \
    } while(0)

    for (int argi = 1; argi < argc; ++argi) {
        const char* arg = argv[argi];
        const char* prog = argv[0];

        if (strcmp(arg, "-h") == 0) {
            usage(prog);
            return 0;
        } else if (strcmp(arg, "-e") == 0) {
            encrypt = 1;
        } else if (strcmp(arg, "-d") == 0) {
            encrypt = 0;
        } else if (strcmp(arg, "-C") == 0 || strcmp(arg, "--chacha") == 0) {
            algo_type = CHACHA20_POLY1305;
        } else if (strcmp(arg, "-X") == 0 || strcmp(arg, "--xchacha") == 0) {
            algo_type = XCHACHA20_POLY1305;
        } else if (strcmp(arg, "-P") == 0 || strcmp(arg, "--p1271") == 0) {
            algo_type = CHACHA20_POLY1271;
        } else if (strcmp(arg, "-PX") == 0 || strcmp(arg, "--xp1271") == 0) {
            algo_type = XCHACHA20_POLY1271;
        } else if (strcmp(arg, "-k") == 0) {
            CHECK_NEXT_ARG();
            int bytes_read = cfx_parse_str(argv[argi+1], key, sizeof(key), CFX_STR_FMT_AUTO);
            if (bytes_read != 32) {
                fprintf(stderr, "error reading key\n");
                ret = 1;
                goto cleanup;
            }
            argi++;
        } else if (strcmp(arg, "-K") == 0) {
            CHECK_NEXT_ARG();
            const char* fname = argv[argi+1];
            enum cfx_str_format fmt = cfx_detect_file_format(fname, NULL);
            if (try_reading_file_into(fname, key, sizeof(key), fmt) != 0) {
                fprintf(stderr, "error reading key file\n");
                ret = 1;
                goto cleanup;
            }
            argi++;
        } else if (strcmp(arg, "-n") == 0) {
            CHECK_NEXT_ARG();
            noncestr = argv[argi+1];
            argi++;
        } else if (strcmp(arg, "--random-nonce") == 0) {
            nonce_prefix = 2;  /* flag: prepend nonce to output, generate later */
        } else if (strcmp(arg, "--nonce-prefix") == 0) {
            nonce_prefix = 1;
        } else if (strcmp(arg, "-a") == 0) {
            CHECK_NEXT_ARG();
            const char* aadstr = argv[argi+1];
            aad_len = strlen(aadstr);
            aad = (uint8_t*)malloc(aad_len);
            if (!aad) {
                fprintf(stderr, "error allocating aad of length %zu\n", aad_len);
                ret = 1;
                goto cleanup;
            }
            int bytes_read = cfx_parse_str(aadstr, aad, aad_len, CFX_STR_FMT_AUTO);
            if (bytes_read <= 0) {
                fprintf(stderr, "error reading aad\n");
                ret = 1;
                goto cleanup;
            }
            aad_len = (size_t)bytes_read;
            argi++;
        } else if (strcmp(arg, "-A") == 0) {
            CHECK_NEXT_ARG();
            const char* fname = argv[argi+1];
            size_t fsize = 0;
            enum cfx_str_format fmt = cfx_detect_file_format(fname, &fsize);
            aad = (uint8_t*)malloc(fsize);
            if (try_reading_file_into(fname, aad, fsize, fmt) != 0) {
                fprintf(stderr, "error reading aad file, size: %zu\n", fsize);
                ret = 1;
                goto cleanup;
            }
            aad_len = fsize;
            argi++;
        } else if (strcmp(arg, "-x") == 0) {
            outfmt = CFX_STR_FMT_HEX;
        } else if (strcmp(arg, "-b64") == 0) {
            outfmt = CFX_STR_FMT_BASE64;
        } else if ((strcmp(arg, "-r") == 0) || (strcmp(arg, "--raw") == 0)) {
            outfmt = CFX_STR_FMT_BINARY;
        } else if ((strcmp(arg, "-v") == 0) || (strcmp(arg, "--verbose") == 0)) {
            verbose = 1;
        } else if (strcmp(arg, "-s") == 0) {
            CHECK_NEXT_ARG();
            const char* str = argv[argi+1];
            pt_in_len = strlen(str);
            pt_in = (uint8_t*)malloc(pt_in_len);
            if (!pt_in) {
                fprintf(stderr, "error allocating plaintext\n");
                ret = 1;
                goto cleanup;
            }
            memcpy(pt_in, str, pt_in_len);
            argi++;
        } else if (arg[0] == '-' && arg[1] != '\0') {
            fprintf(stderr, "Unknown option: %s\n", argv[argi]);
            ret = 1;
            goto cleanup;
        } else {  /* read plaintext from file or stdin */
            const char *ptfn = arg;
            FILE *f = NULL;

            if (strcmp(ptfn, "-") == 0) {
                f = stdin;
            } else {
                f = fopen(ptfn, "rb");
                if (!f) {
                    perror(ptfn);
                    ret = 1;
                    goto cleanup;
                }
            }

            if (cfx_read_all_file(f, &pt_in, &pt_in_len) != 0) {
                fprintf(stderr, "Failed to read plaintext\n");
                if (f != stdin) fclose(f);
                ret = 1;
                goto cleanup;
            }

            if (f != stdin) fclose(f);
        }
    }
    uint8_t* pt = pt_in;
    size_t pt_len = pt_in_len;
    size_t nonce_len =
        (algo_type == XCHACHA20_POLY1305 || algo_type == XCHACHA20_POLY1271) ? 24 : 12;

    /* parse nonce if provided */
    if (noncestr) {
        int bytes_read = cfx_parse_str(noncestr, nonce, nonce_len, CFX_STR_FMT_AUTO);
        if ((size_t)bytes_read != nonce_len) {
            fprintf(stderr, "error: nonce must be %zu bytes\n", nonce_len);
            ret = 1;
            goto cleanup;
        }
    }

    /* generate random nonce if requested */
    if (nonce_prefix == 2) {
        cfx_srand_os();
        cfx_rand_bytes(nonce, nonce_len);
    }

    if (nonce_prefix == 1 && !encrypt) { /* decrypt: read nonce from input */
        if (pt_len <= nonce_len) {
            fprintf(stderr, "nonce prefix means input must be greater than %zu bytes!\n", nonce_len);
            ret = 1;
            goto cleanup;
        }
        memcpy(nonce, pt, nonce_len);
        pt += nonce_len;
        pt_len -= nonce_len;
    } 

    if (verbose) {
        printf("algorithm: %s\n", use_xchacha ? "XChaCha20-Poly1305" : "ChaCha20-Poly1305");
        printf("nonce: ");
        cfx_print_hex(nonce, nonce_len);
        printf("\n");
        printf("aad: ");
        cfx_print_hex(aad, aad_len);
        printf("\n");
    }

#ifdef _WIN32
    if (outfmt == CFX_STR_FMT_BINARY) {
        _setmode(_fileno(stdout), _O_BINARY);
    }
#endif

    /* encrypt */
    if (encrypt) {
        uint8_t* ct = (uint8_t*)malloc(pt_len);
        uint8_t tag[16];
        if (!ct) {
            ret = 1;
            goto cleanup;
        }
        if (use_xchacha) {
            cfx_xchacha20_poly1305_encrypt(ct, tag, pt, pt_len, aad, aad_len, key, nonce);
        } else {
            cfx_chacha20_poly1305_encrypt(ct, tag, pt, pt_len, aad, aad_len, key, nonce);
        }

        if (outfmt == CFX_STR_FMT_BINARY) {
            if (nonce_prefix) fwrite(nonce, 1, nonce_len, stdout);
            fwrite(ct, 1, pt_len, stdout);
            fwrite(tag, 1, sizeof(tag), stdout);
        } else {
            if (nonce_prefix) cfx_printf_output(nonce, nonce_len, outfmt);
            cfx_printf_output(ct, pt_len, outfmt);
            cfx_printf_output(tag, sizeof(tag), outfmt);
            printf("\n");
        }
        cfx_memzero_s(ct, pt_len);
        cfx_memzero_s(tag, sizeof(tag));
        free(ct);
    } else {
        /* decrypt */
        if (pt_len < 16) {
            fprintf(stderr, "ciphertext too short (need at least 16 bytes for tag)\n");
            ret = 1;
            goto cleanup;
        }
        size_t ct_len = pt_len - 16;
        uint8_t* ct = pt;
        uint8_t* tag = pt + ct_len;
        uint8_t* plaintext = (uint8_t*)malloc(ct_len);
        if (!plaintext) {
            ret = 1;
            goto cleanup;
        }

        int ok;
        if (use_xchacha) {
            ok = cfx_xchacha20_poly1305_decrypt(plaintext, ct, ct_len, aad, aad_len, key, nonce, tag);
        } else {
            ok = cfx_chacha20_poly1305_decrypt(plaintext, ct, ct_len, aad, aad_len, key, nonce, tag);
        }
        if (ok != 0) {
            fprintf(stderr, "authentication failed\n");
            cfx_memzero_s(plaintext, ct_len);
            free(plaintext);
            ret = 1;
            goto cleanup;
        }

        if (outfmt == CFX_STR_FMT_BINARY) {
            fwrite(plaintext, 1, ct_len, stdout);
        } else {
            cfx_printf_output(plaintext, ct_len, outfmt);
            printf("\n");
        }
        cfx_memzero_s(plaintext, ct_len);
        free(plaintext);
    }

cleanup:
    cfx_memzero_s(key, sizeof(key));
    cfx_memzero_s(nonce, 24);
    if (aad) {
        cfx_memzero_s(aad, aad_len);
        free(aad);
    }
    if (pt_in) {
        cfx_memzero_s(pt_in, pt_in_len);
        free(pt_in);
    }
    return ret;
}