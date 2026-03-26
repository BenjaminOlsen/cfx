/* cfx_randart.c - drunken bishop randomart for arbitrary data */

#include "cfx_randomart.h"
#include "cfx/sha256.h"
#include "cfx/sha512.h"
#include "cfx/sha3.h"
#include "cfx/blake2.h"
#include "cfx/base64.h"
#include "common.h"
#include "cfx_cmd.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
#include <io.h>
#include <fcntl.h>
#endif

typedef enum {
    RA_HASH_NONE,
    RA_HASH_SHA256,
    RA_HASH_SHA512,
    RA_HASH_SHA3_256,
    RA_HASH_SHA3_512,
    RA_HASH_SHAKE256,
    RA_HASH_BLAKE2B,
    RA_HASH_BLAKE2S
} ra_hash_t;

static void usage(const char *prog) {
    printf("Usage: %s [options] [input]\n", prog);
    printf("  Drunken bishop randomart visualization.\n\n");
    printf("Input modes:\n");
    printf("  -s STRING     Use string as input\n");
    printf("  -f FILE       Read file as input\n");
    printf("  -x HEX        Hex bytes (colons stripped)\n");
    printf("  -x            Read hex from stdin\n");
    printf("  -b64 DATA     Base64-encoded bytes\n");
    printf("  -b64          Read base64 from stdin\n");
    printf("  (default)     Read raw bytes from stdin\n\n");
    printf("Options:\n");
    printf("  -a ALGO       Hash before walking (default: sha256)\n");
    printf("                sha256 sha512 sha3-256 sha3-512\n");
    printf("                shake256 blake2b blake2s none/raw\n");
    printf("  -w WIDTH      Grid width  (odd, default 17)\n");
    printf("  -h HEIGHT     Grid height (odd, default 9)\n");
    printf("  -l LABEL      Frame label\n");
    printf("  --help        Show this help\n\n");
    printf("Examples:\n");
    printf("  %s -s \"hello world\"\n", prog);
    printf("  %s -f key.pub -a sha512 -w 33 -h 17\n", prog);
    printf("  %s -x b69ca37b -a none\n", prog);
    printf("  echo data | %s -a shake256 -w 49 -h 25\n", prog);
}

static ra_hash_t parse_algo(const char *s) {
    if (strcmp(s, "none") == 0 || strcmp(s, "raw") == 0) return RA_HASH_NONE;
    if (strcmp(s, "sha256") == 0) return RA_HASH_SHA256;
    if (strcmp(s, "sha512") == 0) return RA_HASH_SHA512;
    if (strcmp(s, "sha3-256") == 0) return RA_HASH_SHA3_256;
    if (strcmp(s, "sha3-512") == 0) return RA_HASH_SHA3_512;
    if (strcmp(s, "shake256") == 0) return RA_HASH_SHAKE256;
    if (strcmp(s, "blake2b") == 0) return RA_HASH_BLAKE2B;
    if (strcmp(s, "blake2s") == 0) return RA_HASH_BLAKE2S;
    return (ra_hash_t)-1;
}

static const char *algo_name(ra_hash_t a) {
    switch (a) {
    case RA_HASH_NONE:     return "RAW";
    case RA_HASH_SHA256:   return "SHA256";
    case RA_HASH_SHA512:   return "SHA512";
    case RA_HASH_SHA3_256: return "SHA3-256";
    case RA_HASH_SHA3_512: return "SHA3-512";
    case RA_HASH_SHAKE256: return "SHAKE256";
    case RA_HASH_BLAKE2B:  return "BLAKE2b";
    case RA_HASH_BLAKE2S:  return "BLAKE2s";
    }
    return "?";
}

/* strip colons from hex string in-place */
static void strip_colons(char *s) {
    char *dst = s;
    for (char *src = s; *src; src++) {
        if (*src != ':') *dst++ = *src;
    }
    *dst = '\0';
}

/* decode base64 string, return malloc'd bytes or NULL */
static uint8_t *decode_b64(const char *s, size_t *out_len) {
    size_t slen = strlen(s);
    size_t max = slen; /* decoded always <= encoded length */
    uint8_t *buf = malloc(max);
    if (!buf) return NULL;
    size_t n = max;
    if (cfx_base64_decode(buf, &n, s, slen) != 0) {
        free(buf);
        return NULL;
    }
    *out_len = n;
    return buf;
}

/* hash input data, return malloc'd digest and its length */
static uint8_t *hash_data(ra_hash_t algo, const uint8_t *in, size_t in_len,
                          int w, int h, size_t *out_len) {
    uint8_t *out = NULL;
    size_t n = 0;

    switch (algo) {
    case RA_HASH_NONE:
        /* pass through */
        out = malloc(in_len);
        if (out) { memcpy(out, in, in_len); n = in_len; }
        break;
    case RA_HASH_SHA256:
        n = 32;
        out = malloc(n);
        if (out) cfx_sha256(out, in, in_len);
        break;
    case RA_HASH_SHA512:
        n = 64;
        out = malloc(n);
        if (out) cfx_sha512(out, in, in_len);
        break;
    case RA_HASH_SHA3_256:
        n = 32;
        out = malloc(n);
        if (out) cfx_sha3_256(out, in, in_len);
        break;
    case RA_HASH_SHA3_512:
        n = 64;
        out = malloc(n);
        if (out) cfx_sha3_512(out, in, in_len);
        break;
    case RA_HASH_SHAKE256:
        /* auto-size: ~(w*h)/8 bytes for good coverage */
        n = (size_t)(w * h) / 8;
        if (n < 32) n = 32;
        out = malloc(n);
        if (out) cfx_shake256(out, n, in, in_len);
        break;
    case RA_HASH_BLAKE2B:
        n = 64;
        out = malloc(n);
        if (out) cfx_blake2b(out, n, in, in_len, NULL, 0);
        break;
    case RA_HASH_BLAKE2S:
        n = 32;
        out = malloc(n);
        if (out) cfx_blake2s(out, n, in, in_len, NULL, 0);
        break;
    }

    *out_len = n;
    return out;
}

int cfx_randart_run(int argc, char **argv) {
    if (argc == 2 && (strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0)) {
        usage(argv[0]);
        return 0;
    }

    ra_hash_t algo = RA_HASH_SHA256;
    int w = CFX_RA_DEFAULT_WIDTH;
    int h = CFX_RA_DEFAULT_HEIGHT;
    const char *label = NULL;
    const char *str_input = NULL;
    const char *hex_input = NULL;
    const char *file_input = NULL;
    const char *b64_input = NULL;
    int hex_stdin = 0;  /* -x with no arg */
    int b64_stdin = 0;  /* -b64 with no arg */

    /* parse args */
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-s") == 0) {
            if (++i >= argc) { fprintf(stderr, "-s needs argument\n"); return 1; }
            str_input = argv[i];
        } else if (strcmp(argv[i], "-f") == 0) {
            if (++i >= argc) { fprintf(stderr, "-f needs argument\n"); return 1; }
            file_input = argv[i];
        } else if (strcmp(argv[i], "-x") == 0) {
            /* -x HEX or just -x (read hex from stdin) */
            if (i + 1 < argc && argv[i + 1][0] != '-') {
                hex_input = argv[++i];
            } else {
                hex_stdin = 1;
            }
        } else if (strcmp(argv[i], "-b64") == 0) {
            /* -b64 DATA or just -b64 (read from stdin) */
            if (i + 1 < argc && argv[i + 1][0] != '-') {
                b64_input = argv[++i];
            } else {
                b64_stdin = 1;
            }
        } else if (strcmp(argv[i], "-a") == 0) {
            if (++i >= argc) { fprintf(stderr, "-a needs argument\n"); return 1; }
            algo = parse_algo(argv[i]);
            if ((int)algo < 0) {
                fprintf(stderr, "unknown algorithm: %s\n", argv[i]);
                return 1;
            }
        } else if (strcmp(argv[i], "-w") == 0) {
            if (++i >= argc) { fprintf(stderr, "-w needs argument\n"); return 1; }
            w = atoi(argv[i]);
        } else if (strcmp(argv[i], "-h") == 0) {
            if (++i >= argc) { fprintf(stderr, "-h needs argument\n"); return 1; }
            h = atoi(argv[i]);
        } else if (strcmp(argv[i], "-l") == 0) {
            if (++i >= argc) { fprintf(stderr, "-l needs argument\n"); return 1; }
            label = argv[i];
        } else if (strcmp(argv[i], "--help") == 0) {
            usage(argv[0]);
            return 0;
        } else if (argv[i][0] == '-') {
            fprintf(stderr, "unknown option: %s\n", argv[i]);
            return 1;
        } else {
            /* bare arg: treat as string */
            str_input = argv[i];
        }
    }

    /* validate grid */
    if (w < 5 || w > CFX_RA_MAX_WIDTH || (w % 2) == 0) {
        fprintf(stderr, "width must be odd, 5..%d\n", CFX_RA_MAX_WIDTH);
        return 1;
    }
    if (h < 5 || h > CFX_RA_MAX_HEIGHT || (h % 2) == 0) {
        fprintf(stderr, "height must be odd, 5..%d\n", CFX_RA_MAX_HEIGHT);
        return 1;
    }

    /* get input bytes */
    uint8_t *data = NULL;
    size_t data_len = 0;

    if (str_input) {
        data_len = strlen(str_input);
        data = malloc(data_len);
        if (!data) { fprintf(stderr, "alloc failed\n"); return 1; }
        memcpy(data, str_input, data_len);
    } else if (hex_input) {
        /* mutable copy for colon stripping */
        char *hex = strdup(hex_input);
        if (!hex) { fprintf(stderr, "alloc failed\n"); return 1; }
        strip_colons(hex);

        size_t max_bytes = strlen(hex) / 2 + 1;
        data = malloc(max_bytes);
        if (!data) { free(hex); fprintf(stderr, "alloc failed\n"); return 1; }
        int nb = cfx_parse_hex_auto(hex, data, max_bytes);
        free(hex);
        if (nb <= 0) {
            fprintf(stderr, "bad hex\n");
            free(data);
            return 1;
        }
        data_len = (size_t)nb;
    } else if (b64_input) {
        data = decode_b64(b64_input, &data_len);
        if (!data) { fprintf(stderr, "bad base64\n"); return 1; }
    } else if (file_input) {
        FILE *f = fopen(file_input, "rb");
        if (!f) { perror(file_input); return 1; }
        if (cfx_read_all_file(f, &data, &data_len) != 0) {
            fprintf(stderr, "error reading %s\n", file_input);
            fclose(f);
            return 1;
        }
        fclose(f);
    } else if (b64_stdin) {
        char *line = cfx_read_line_stdin();
        if (!line || !*line) {
            fprintf(stderr, "no base64 input on stdin\n");
            free(line);
            return 1;
        }
        data = decode_b64(line, &data_len);
        free(line);
        if (!data) { fprintf(stderr, "bad base64 on stdin\n"); return 1; }
    } else if (hex_stdin) {
        /* read one line of hex from stdin */
        char *line = cfx_read_line_stdin();
        if (!line || !*line) {
            fprintf(stderr, "no hex input on stdin\n");
            free(line);
            return 1;
        }
        strip_colons(line);
        size_t max_bytes = strlen(line) / 2 + 1;
        data = malloc(max_bytes);
        if (!data) { free(line); fprintf(stderr, "alloc failed\n"); return 1; }
        int nb = cfx_parse_hex_auto(line, data, max_bytes);
        free(line);
        if (nb <= 0) {
            fprintf(stderr, "bad hex on stdin\n");
            free(data);
            return 1;
        }
        data_len = (size_t)nb;
    } else {
        /* stdin raw bytes */
#ifdef _WIN32
        _setmode(_fileno(stdin), _O_BINARY);
#endif
        if (cfx_read_all_file(stdin, &data, &data_len) != 0) {
            fprintf(stderr, "error reading stdin\n");
            return 1;
        }
    }

    if (data_len == 0) {
        fprintf(stderr, "no input data\n");
        free(data);
        return 1;
    }

    /* hash */
    size_t walk_len = 0;
    uint8_t *walk_data = hash_data(algo, data, data_len, w, h, &walk_len);
    free(data);
    if (!walk_data) { fprintf(stderr, "alloc failed\n"); return 1; }

    /* walk + render */
    size_t art_sz = (size_t)h * ((size_t)w + 1);
    char *art = malloc(art_sz);
    if (!art) { free(walk_data); fprintf(stderr, "alloc failed\n"); return 1; }

    cfx_randomart(art, w, h, walk_data, walk_len);
    free(walk_data);

    /* default label = algo name */
    const char *lbl = label ? label : algo_name(algo);
    cfx_print_randomart_frame(art, w, h, lbl);

    free(art);
    return 0;
}
