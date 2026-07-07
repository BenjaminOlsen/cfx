/* cfx_hash.c - SHA256/SHA3/BLAKE2/SipHash hashing utility */

#include "cfx/sha256.h"
#include "cfx/sha3.h"
#include "cfx/blake2.h"
#include "cfx/siphash.h"
#include "cfx/base64.h"
#include "cfx_utils_common.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
#include <io.h>
#include <fcntl.h>
#endif

#include "cfx_cmd.h"

typedef enum {
    HASH_SHA256,
    HASH_SHA3_224,
    HASH_SHA3_256,
    HASH_SHA3_384,
    HASH_SHA3_512,
    HASH_SHAKE128,
    HASH_SHAKE256,
    HASH_BLAKE2B,
    HASH_BLAKE2S,
    HASH_SIPHASH
} hash_algo_t;

static void usage(const char* prog) {
    printf("Usage: %s [options] [file...]\n", prog);
    printf("  Compute cryptographic hash of files or stdin.\n\n");
    printf("Algorithms:\n");
    printf("  (default)       SHA-256 (32 bytes)\n");
    printf("  --sha3-224      SHA3-224 (28 bytes)\n");
    printf("  --sha3-256      SHA3-256 (32 bytes)\n");
    printf("  --sha3-384      SHA3-384 (48 bytes)\n");
    printf("  --sha3-512      SHA3-512 (64 bytes)\n");
    printf("  --shake128      SHAKE128 XOF (variable, default 32)\n");
    printf("  --shake256      SHAKE256 XOF (variable, default 64)\n");
    printf("  --blake2b       BLAKE2b (1-64 bytes, default 64)\n");
    printf("  --blake2s       BLAKE2s (1-32 bytes, default 32)\n");
    printf("  --siphash       SipHash-2-4 (8 bytes, requires -k)\n\n");
    printf("Options:\n");
    printf("  -s <string>     Hash the given string instead of a file\n");
    printf("  -k <hex>        Key for keyed hashing (BLAKE2: up to 64/32 bytes, SipHash: 16 bytes)\n");
    printf("  -n <len>        Output length in bytes (BLAKE2/SHAKE only)\n");
    printf("  -f              Include filename in output\n");
    printf("  -x              Output as hex (default)\n");
    printf("  -b64            Output as base64\n");
    printf("  -h, --help      Show this help message\n\n");
    printf("Examples:\n");
    printf("  %s file.txt                        SHA256 hash a file\n", prog);
    printf("  %s -s \"hello\"                      SHA256 hash a string\n", prog);
    printf("  %s --sha3-256 file.txt             SHA3-256 hash a file\n", prog);
    printf("  %s --shake256 -n 128 file.txt      SHAKE256 with 128-byte output\n", prog);
    printf("  %s --blake2b file.txt              BLAKE2b (64-byte output)\n", prog);
    printf("  %s --blake2b -n 32 file.txt        BLAKE2b with 32-byte output\n", prog);
    printf("  %s --blake2b -k deadbeef file.txt  BLAKE2b keyed (MAC mode)\n", prog);
    printf("  %s --blake2s -s \"hello\"            BLAKE2s hash a string\n", prog);
    printf("  %s --siphash -k 00010203...0f f    SipHash with 16-byte key\n", prog);
}

static void print_hash(const uint8_t* hash, size_t hash_len, const char* name, enum cfx_str_format fmt) {
    if (fmt == CFX_STR_FMT_BASE64) {
        size_t b64_len = cfx_base64_enc_len(hash_len);
        char b64[256];  /* enough for 128-byte hashes */
        cfx_base64_encode(b64, &b64_len, hash, hash_len);
        b64[b64_len] = '\0';
        printf("%s", b64);
    } else {
        for (size_t i = 0; i < hash_len; i++) {
            printf("%02x", hash[i]);
        }
    }
    if (name) {
        printf("  %s", name);
    }
    printf("\n");
}

static int hash_string_sha256(const char* str, enum cfx_str_format fmt) {
    cfx_sha256_ctx ctx;
    uint8_t hash[32];

    cfx_sha256_init(&ctx);
    cfx_sha256_update(&ctx, (const uint8_t*)str, strlen(str));
    cfx_sha256_final(&ctx, hash);

    print_hash(hash, 32, NULL, fmt);
    return 0;
}

static int hash_string_blake2b(const char* str, size_t outlen,
                               const uint8_t* key, size_t keylen,
                               enum cfx_str_format fmt) {
    uint8_t hash[64];
    if (cfx_blake2b(hash, outlen, str, strlen(str), key, keylen) != 0) {
        fprintf(stderr, "BLAKE2b error\n");
        return 1;
    }
    print_hash(hash, outlen, NULL, fmt);
    return 0;
}

static int hash_string_blake2s(const char* str, size_t outlen,
                               const uint8_t* key, size_t keylen,
                               enum cfx_str_format fmt) {
    uint8_t hash[32];
    if (cfx_blake2s(hash, outlen, str, strlen(str), key, keylen) != 0) {
        fprintf(stderr, "BLAKE2s error\n");
        return 1;
    }
    print_hash(hash, outlen, NULL, fmt);
    return 0;
}

static int hash_string_siphash(const char* str, const uint8_t key[16], enum cfx_str_format fmt) {
    uint64_t h = cfx_siphash((const uint8_t*)str, strlen(str), key);
    uint8_t hash[8];
    for (int i = 0; i < 8; i++) {
        hash[i] = (uint8_t)(h >> (8 * i));
    }
    print_hash(hash, 8, NULL, fmt);
    return 0;
}

static int hash_string_sha3(hash_algo_t algo, const char* str, size_t outlen,
                            enum cfx_str_format fmt) {
    uint8_t hash[64];
    size_t len = strlen(str);

    switch (algo) {
        case HASH_SHA3_224:
            cfx_sha3_224(hash, str, len);
            print_hash(hash, 28, NULL, fmt);
            break;
        case HASH_SHA3_256:
            cfx_sha3_256(hash, str, len);
            print_hash(hash, 32, NULL, fmt);
            break;
        case HASH_SHA3_384:
            cfx_sha3_384(hash, str, len);
            print_hash(hash, 48, NULL, fmt);
            break;
        case HASH_SHA3_512:
            cfx_sha3_512(hash, str, len);
            print_hash(hash, 64, NULL, fmt);
            break;
        case HASH_SHAKE128:
            cfx_shake128(hash, outlen, str, len);
            print_hash(hash, outlen, NULL, fmt);
            break;
        case HASH_SHAKE256:
            cfx_shake256(hash, outlen, str, len);
            print_hash(hash, outlen, NULL, fmt);
            break;
        default:
            return 1;
    }
    return 0;
}

static int hash_file_sha256(FILE* f, const char* filename, enum cfx_str_format fmt) {
    cfx_sha256_ctx ctx;
    uint8_t hash[32];
    uint8_t buf[8192];
    size_t n;

    cfx_sha256_init(&ctx);
    while ((n = fread(buf, 1, sizeof(buf), f)) > 0) {
        cfx_sha256_update(&ctx, buf, n);
    }

    if (ferror(f)) {
        fprintf(stderr, "Error reading %s\n", filename ? filename : "stdin");
        return 1;
    }

    cfx_sha256_final(&ctx, hash);
    print_hash(hash, 32, filename, fmt);
    return 0;
}

static int hash_file_blake2b(FILE* f, const char* filename, size_t outlen,
                             const uint8_t* key, size_t keylen,
                             enum cfx_str_format fmt) {
    cfx_blake2b_ctx_t ctx;
    uint8_t hash[64];
    uint8_t buf[8192];
    size_t n;

    if (key && keylen > 0) {
        if (cfx_blake2b_init_key(&ctx, outlen, key, keylen) != 0) {
            fprintf(stderr, "BLAKE2b init error\n");
            return 1;
        }
    } else {
        if (cfx_blake2b_init(&ctx, outlen) != 0) {
            fprintf(stderr, "BLAKE2b init error\n");
            return 1;
        }
    }

    while ((n = fread(buf, 1, sizeof(buf), f)) > 0) {
        cfx_blake2b_update(&ctx, buf, n);
    }

    if (ferror(f)) {
        fprintf(stderr, "Error reading %s\n", filename ? filename : "stdin");
        return 1;
    }

    cfx_blake2b_final(&ctx, hash);
    print_hash(hash, outlen, filename, fmt);
    return 0;
}

static int hash_file_blake2s(FILE* f, const char* filename, size_t outlen,
                             const uint8_t* key, size_t keylen,
                             enum cfx_str_format fmt) {
    cfx_blake2s_ctx_t ctx;
    uint8_t hash[32];
    uint8_t buf[8192];
    size_t n;

    if (key && keylen > 0) {
        if (cfx_blake2s_init_key(&ctx, outlen, key, keylen) != 0) {
            fprintf(stderr, "BLAKE2s init error\n");
            return 1;
        }
    } else {
        if (cfx_blake2s_init(&ctx, outlen) != 0) {
            fprintf(stderr, "BLAKE2s init error\n");
            return 1;
        }
    }

    while ((n = fread(buf, 1, sizeof(buf), f)) > 0) {
        cfx_blake2s_update(&ctx, buf, n);
    }

    if (ferror(f)) {
        fprintf(stderr, "Error reading %s\n", filename ? filename : "stdin");
        return 1;
    }

    cfx_blake2s_final(&ctx, hash);
    print_hash(hash, outlen, filename, fmt);
    return 0;
}

static int hash_file_siphash(FILE* f, const char* filename, const uint8_t key[16], enum cfx_str_format fmt) {
    cfx_siphash_ctx_t ctx;
    uint8_t buf[8192];
    size_t n;

    cfx_siphash_init(&ctx, key);
    while ((n = fread(buf, 1, sizeof(buf), f)) > 0) {
        cfx_siphash_update(&ctx, buf, n);
    }

    if (ferror(f)) {
        fprintf(stderr, "Error reading %s\n", filename ? filename : "stdin");
        return 1;
    }

    uint64_t h = cfx_siphash_final(&ctx);
    uint8_t hash[8];
    for (int i = 0; i < 8; i++) {
        hash[i] = (uint8_t)(h >> (8 * i));
    }
    print_hash(hash, 8, filename, fmt);
    return 0;
}

static int hash_file_sha3(hash_algo_t algo, FILE* f, const char* filename,
                          size_t outlen, enum cfx_str_format fmt) {
    cfx_sha3_ctx_t ctx;
    uint8_t hash[64];
    uint8_t buf[8192];
    size_t n;

    switch (algo) {
        case HASH_SHA3_224: cfx_sha3_224_init(&ctx); break;
        case HASH_SHA3_256: cfx_sha3_256_init(&ctx); break;
        case HASH_SHA3_384: cfx_sha3_384_init(&ctx); break;
        case HASH_SHA3_512: cfx_sha3_512_init(&ctx); break;
        case HASH_SHAKE128: cfx_shake128_init(&ctx); break;
        case HASH_SHAKE256: cfx_shake256_init(&ctx); break;
        default: return 1;
    }

    while ((n = fread(buf, 1, sizeof(buf), f)) > 0) {
        switch (algo) {
            case HASH_SHA3_224: cfx_sha3_224_update(&ctx, buf, n); break;
            case HASH_SHA3_256: cfx_sha3_256_update(&ctx, buf, n); break;
            case HASH_SHA3_384: cfx_sha3_384_update(&ctx, buf, n); break;
            case HASH_SHA3_512: cfx_sha3_512_update(&ctx, buf, n); break;
            case HASH_SHAKE128: cfx_shake128_update(&ctx, buf, n); break;
            case HASH_SHAKE256: cfx_shake256_update(&ctx, buf, n); break;
            default: break;
        }
    }

    if (ferror(f)) {
        fprintf(stderr, "Error reading %s\n", filename ? filename : "stdin");
        return 1;
    }

    size_t hashlen;
    switch (algo) {
        case HASH_SHA3_224:
            cfx_sha3_224_final(&ctx, hash);
            hashlen = 28;
            break;
        case HASH_SHA3_256:
            cfx_sha3_256_final(&ctx, hash);
            hashlen = 32;
            break;
        case HASH_SHA3_384:
            cfx_sha3_384_final(&ctx, hash);
            hashlen = 48;
            break;
        case HASH_SHA3_512:
            cfx_sha3_512_final(&ctx, hash);
            hashlen = 64;
            break;
        case HASH_SHAKE128:
            cfx_shake128_final(&ctx);
            cfx_shake128_squeeze(&ctx, hash, outlen);
            hashlen = outlen;
            break;
        case HASH_SHAKE256:
            cfx_shake256_final(&ctx);
            cfx_shake256_squeeze(&ctx, hash, outlen);
            hashlen = outlen;
            break;
        default:
            return 1;
    }

    print_hash(hash, hashlen, filename, fmt);
    return 0;
}

static int hash_string(hash_algo_t algo, const char* str, size_t outlen,
                       const uint8_t* key, size_t keylen,
                       enum cfx_str_format fmt) {
    switch (algo) {
        case HASH_SHA256:
            return hash_string_sha256(str, fmt);
        case HASH_SHA3_224:
        case HASH_SHA3_256:
        case HASH_SHA3_384:
        case HASH_SHA3_512:
        case HASH_SHAKE128:
        case HASH_SHAKE256:
            return hash_string_sha3(algo, str, outlen, fmt);
        case HASH_BLAKE2B:
            return hash_string_blake2b(str, outlen, key, keylen, fmt);
        case HASH_BLAKE2S:
            return hash_string_blake2s(str, outlen, key, keylen, fmt);
        case HASH_SIPHASH:
            return hash_string_siphash(str, key, fmt);
    }
    return 1;
}

static int hash_file(hash_algo_t algo, FILE* f, const char* filename, size_t outlen,
                     const uint8_t* key, size_t keylen,
                     enum cfx_str_format fmt) {
    switch (algo) {
        case HASH_SHA256:
            return hash_file_sha256(f, filename, fmt);
        case HASH_SHA3_224:
        case HASH_SHA3_256:
        case HASH_SHA3_384:
        case HASH_SHA3_512:
        case HASH_SHAKE128:
        case HASH_SHAKE256:
            return hash_file_sha3(algo, f, filename, outlen, fmt);
        case HASH_BLAKE2B:
            return hash_file_blake2b(f, filename, outlen, key, keylen, fmt);
        case HASH_BLAKE2S:
            return hash_file_blake2s(f, filename, outlen, key, keylen, fmt);
        case HASH_SIPHASH:
            return hash_file_siphash(f, filename, key, fmt);
    }
    return 1;
}

int cfx_hash_run(int argc, char** argv) {
    if (argc == 2 && (strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0)) {
        usage(argv[0]);
        return 0;
    }

    int rc = 0;
    int show_filename = 0;
    enum cfx_str_format fmt = CFX_STR_FMT_HEX;
    hash_algo_t algo = HASH_SHA256;
    uint8_t key[64] = {0};
    size_t keylen = 0;
    size_t outlen = 0;  /* 0 = use default */

    /* collect string inputs and file inputs */
    const char* strings[64];
    const char* files[64];
    int nstrings = 0;
    int nfiles = 0;

    /* first pass: parse all options */
    for (int argi = 1; argi < argc; argi++) {
        if (strcmp(argv[argi], "-s") == 0) {
            if (argi + 1 >= argc) {
                fprintf(stderr, "Error: -s requires a string argument\n");
                return 1;
            }
            if (nstrings < 64) {
                strings[nstrings++] = argv[argi + 1];
            }
            argi++;
        } else if (strcmp(argv[argi], "-x") == 0) {
            fmt = CFX_STR_FMT_HEX;
        } else if (strcmp(argv[argi], "-b64") == 0) {
            fmt = CFX_STR_FMT_BASE64;
        } else if (strcmp(argv[argi], "-f") == 0) {
            show_filename = 1;
        } else if (strcmp(argv[argi], "--blake2b") == 0) {
            algo = HASH_BLAKE2B;
        } else if (strcmp(argv[argi], "--blake2s") == 0) {
            algo = HASH_BLAKE2S;
        } else if (strcmp(argv[argi], "--siphash") == 0) {
            algo = HASH_SIPHASH;
        } else if (strcmp(argv[argi], "--sha3-224") == 0) {
            algo = HASH_SHA3_224;
        } else if (strcmp(argv[argi], "--sha3-256") == 0) {
            algo = HASH_SHA3_256;
        } else if (strcmp(argv[argi], "--sha3-384") == 0) {
            algo = HASH_SHA3_384;
        } else if (strcmp(argv[argi], "--sha3-512") == 0) {
            algo = HASH_SHA3_512;
        } else if (strcmp(argv[argi], "--shake128") == 0) {
            algo = HASH_SHAKE128;
        } else if (strcmp(argv[argi], "--shake256") == 0) {
            algo = HASH_SHAKE256;
        } else if (strcmp(argv[argi], "-n") == 0) {
            if (argi + 1 >= argc) {
                fprintf(stderr, "Error: -n requires a length argument\n");
                return 1;
            }
            outlen = (size_t)atoi(argv[argi + 1]);
            argi++;
        } else if (strcmp(argv[argi], "-k") == 0) {
            if (argi + 1 >= argc) {
                fprintf(stderr, "Error: -k requires a hex key argument\n");
                return 1;
            }
            int parsed = cfx_parse_hex_auto(argv[argi + 1], key, sizeof(key));
            if (parsed < 0) {
                fprintf(stderr, "Error: invalid hex key\n");
                return 1;
            }
            keylen = (size_t)parsed;
            argi++;
        } else if (argv[argi][0] == '-' && argv[argi][1] != '\0') {
            fprintf(stderr, "Unknown option: %s\n", argv[argi]);
            usage(argv[0]);
            return 1;
        } else {
            /* It's a filename */
            if (nfiles < 64) {
                files[nfiles++] = argv[argi];
            }
        }
    }

    /* validate siphash key requirement */
    if (algo == HASH_SIPHASH && keylen == 0 && (nstrings > 0 || nfiles > 0)) {
        fprintf(stderr, "Error: --siphash requires -k <key>\n");
        return 1;
    }

    /* compute effective output length */
    size_t out = outlen;
    if (out == 0) {
        if (algo == HASH_BLAKE2B || algo == HASH_SHAKE256) out = 64;
        else if (algo == HASH_BLAKE2S || algo == HASH_SHAKE128) out = 32;
    }

    /* second pass: process inputs */
    for (int i = 0; i < nstrings; i++) {
        rc |= hash_string(algo, strings[i], out, key, keylen, fmt);
    }

    for (int i = 0; i < nfiles; i++) {
        FILE* f = fopen(files[i], "rb");
        if (!f) {
            perror(files[i]);
            rc = 1;
        } else {
            rc |= hash_file(algo, f, show_filename ? files[i] : NULL, out, key, keylen, fmt);
            fclose(f);
        }
    }

    /* no input specified - read from stdin */
    if (nstrings == 0 && nfiles == 0) {
        if (algo == HASH_SIPHASH && keylen == 0) {
            fprintf(stderr, "Error: --siphash requires -k <key>\n");
            return 1;
        }
#ifdef _WIN32
        _setmode(_fileno(stdin), _O_BINARY);
#endif
        rc = hash_file(algo, stdin, show_filename ? "-" : NULL, out, key, keylen, fmt);
    }

    return rc;
}
