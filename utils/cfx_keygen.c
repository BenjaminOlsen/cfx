/* cfx_keygen.c - cryptographic key generation utility */

#include "cfx/rand.h"
#include "cfx/ed25519.h"
#include "cfx/x25519.h"
#include "cfx/base64.h"
#include "cfx/memory.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>

#ifdef _WIN32
#include <io.h>
#include <fcntl.h>
#include <direct.h>
#include <shlobj.h>
#define mkdir_p(path) _mkdir(path)
#else
#include <unistd.h>
#include <pwd.h>
#define mkdir_p(path) mkdir(path, 0700)
#endif

#include "cfx_cmd.h"
#include "misc.h"

typedef enum {
    KEY_RAW,        /* random bytes */
    KEY_ED25519,
    KEY_X25519
} key_type_t;

/* get path to ~/.cfx directory, creating if needed */
static int get_cfx_dir(char *buf, size_t bufsz) {
#ifdef _WIN32
    char home[MAX_PATH];
    if (SHGetFolderPathA(NULL, CSIDL_PROFILE, NULL, 0, home) != S_OK) {
        fprintf(stderr, "error: cannot determine home directory\n");
        return -1;
    }
#else
    const char *home = getenv("HOME");
    if (!home) {
        struct passwd *pw = getpwuid(getuid());
        if (pw) home = pw->pw_dir;
    }
    if (!home) {
        fprintf(stderr, "error: cannot determine home directory\n");
        return -1;
    }
#endif
    int n = snprintf(buf, bufsz, "%s/.cfx", home);
    if (n < 0 || (size_t)n >= bufsz) {
        fprintf(stderr, "error: path too long\n");
        return -1;
    }
    return 0;
}

static int ensure_cfx_dir(void) {
    char path[1024];
    if (get_cfx_dir(path, sizeof(path)) != 0) return -1;

    struct stat st;
    if (stat(path, &st) == 0) {
        if (S_ISDIR(st.st_mode)) return 0;
        fprintf(stderr, "error: %s exists but is not a directory\n", path);
        return -1;
    }

    if (mkdir_p(path) != 0) {
        fprintf(stderr, "error: cannot create %s: %s\n", path, strerror(errno));
        return -1;
    }
    printf("Created %s\n", path);
    return 0;
}

static int write_key_file(const char *path, const uint8_t *data, size_t len, int is_private) {
    FILE *f = fopen(path, "wb");
    if (!f) {
        fprintf(stderr, "error: cannot create %s: %s\n", path, strerror(errno));
        return -1;
    }

    /* write as hex with newline */
    for (size_t i = 0; i < len; i++) {
        fprintf(f, "%02x", data[i]);
    }
    fprintf(f, "\n");
    fclose(f);

#ifndef _WIN32
    /* set permissions: private keys 0600, public keys 0644 */
    if (chmod(path, is_private ? 0600 : 0644) != 0) {
        fprintf(stderr, "warning: cannot set permissions on %s: %s\n", path, strerror(errno));
    }
#endif

    return 0;
}

static void usage(const char* prog) {
    printf("Usage: %s [options] [bytes]\n", prog);
    printf("  Generate cryptographic keys.\n\n");
    printf("Key types:\n");
    printf("  --ed25519       Generate Ed25519 keypair (default identity key)\n");
    printf("  --x25519        Generate X25519 keypair (encryption key)\n");
    printf("  <bytes>         Generate raw random bytes (for symmetric keys)\n\n");
    printf("Options:\n");
    printf("  -f <name>       Output file basename (default: ~/.cfx/id_<type>)\n");
    printf("  -x              Output as hex (default for raw bytes)\n");
    printf("  -b64            Output as base64\n");
    printf("  -r              Output raw bytes (binary)\n");
    printf("  -h, --help      Show this help\n\n");
    printf("Examples:\n");
    printf("  %s --ed25519              Generate Ed25519 identity key\n", prog);
    printf("  %s --x25519               Generate X25519 encryption key\n", prog);
    printf("  %s --ed25519 -f mykey     Write to mykey and mykey.pub\n", prog);
    printf("  %s 32                     Generate 32-byte random key (hex)\n", prog);
    printf("  %s 32 -r > key            Write raw key to file\n", prog);
}

static int keygen_ed25519(const char *basename) {
    char priv_path[1024], pub_path[1024];

    if (basename) {
        snprintf(priv_path, sizeof(priv_path), "%s", basename);
        snprintf(pub_path, sizeof(pub_path), "%s.pub", basename);
    } else {
        char cfx_dir[1024];
        if (get_cfx_dir(cfx_dir, sizeof(cfx_dir)) != 0) return 1;
        if (ensure_cfx_dir() != 0) return 1;
        snprintf(priv_path, sizeof(priv_path), "%s/id_ed25519", cfx_dir);
        snprintf(pub_path, sizeof(pub_path), "%s/id_ed25519.pub", cfx_dir);
    }

    /* check if files exist */
    struct stat st;
    if (stat(priv_path, &st) == 0) {
        fprintf(stderr, "error: %s already exists (use -f to specify different name)\n", priv_path);
        return 1;
    }

    /* generate keypair */
    uint8_t seed[32], pk[32], sk[64];
    cfx_srand_os();
    cfx_rand_bytes(seed, 32);
    cfx_ed25519_create_keypair(pk, sk, seed);

    /* write files - store seed as private key (32 bytes), not expanded sk */
    if (write_key_file(priv_path, seed, 32, 1) != 0) {
        cfx_memzero_s(seed, sizeof(seed));
        cfx_memzero_s(sk, sizeof(sk));
        return 1;
    }
    if (write_key_file(pub_path, pk, 32, 0) != 0) {
        cfx_memzero_s(seed, sizeof(seed));
        cfx_memzero_s(sk, sizeof(sk));
        return 1;
    }

    printf("Generated Ed25519 keypair:\n");
    printf("  Private: %s\n", priv_path);
    printf("  Public:  %s\n", pub_path);

    cfx_memzero_s(seed, sizeof(seed));
    cfx_memzero_s(sk, sizeof(sk));
    return 0;
}

static int keygen_x25519(const char *basename) {
    char priv_path[1024], pub_path[1024];

    if (basename) {
        snprintf(priv_path, sizeof(priv_path), "%s", basename);
        snprintf(pub_path, sizeof(pub_path), "%s.pub", basename);
    } else {
        char cfx_dir[1024];
        if (get_cfx_dir(cfx_dir, sizeof(cfx_dir)) != 0) return 1;
        if (ensure_cfx_dir() != 0) return 1;
        snprintf(priv_path, sizeof(priv_path), "%s/id_x25519", cfx_dir);
        snprintf(pub_path, sizeof(pub_path), "%s/id_x25519.pub", cfx_dir);
    }

    /* check if files exist */
    struct stat st;
    if (stat(priv_path, &st) == 0) {
        fprintf(stderr, "error: %s already exists (use -f to specify different name)\n", priv_path);
        return 1;
    }

    /* generate keypair */
    uint8_t priv[32], pub[32];
    cfx_srand_os();
    cfx_rand_bytes(priv, 32);
    cfx_x25519_base(pub, priv);

    /* write files */
    if (write_key_file(priv_path, priv, 32, 1) != 0) {
        cfx_memzero_s(priv, sizeof(priv));
        return 1;
    }
    if (write_key_file(pub_path, pub, 32, 0) != 0) {
        cfx_memzero_s(priv, sizeof(priv));
        return 1;
    }

    printf("Generated X25519 keypair:\n");
    printf("  Private: %s\n", priv_path);
    printf("  Public:  %s\n", pub_path);

    cfx_memzero_s(priv, sizeof(priv));
    return 0;
}

static int keygen_raw(long nbytes, enum cfx_str_format fmt) {
    uint8_t* buf = malloc((size_t)nbytes);
    if (!buf) {
        fprintf(stderr, "Allocation failed\n");
        return 1;
    }

    cfx_srand_os();
    cfx_rand_bytes(buf, (size_t)nbytes);

    if (fmt == CFX_STR_FMT_BINARY) {
#ifdef _WIN32
        _setmode(_fileno(stdout), _O_BINARY);
#endif
        fwrite(buf, 1, (size_t)nbytes, stdout);
    } else if (fmt == CFX_STR_FMT_BASE64) {
        size_t b64_len = cfx_base64_enc_len((size_t)nbytes);
        char* b64 = malloc(b64_len + 1);
        if (!b64) {
            fprintf(stderr, "Allocation failed\n");
            free(buf);
            return 1;
        }
        cfx_base64_encode(b64, &b64_len, buf, (size_t)nbytes);
        b64[b64_len] = '\0';
        printf("%s\n", b64);
        free(b64);
    } else {
        for (long i = 0; i < nbytes; i++) {
            printf("%02x", buf[i]);
        }
        printf("\n");
    }

    cfx_memzero_s(buf, (size_t)nbytes);
    free(buf);
    return 0;
}

int cfx_keygen_run(int argc, char** argv) {
    if (argc < 2) {
        usage(argv[0]);
        return 1;
    }

    key_type_t key_type = KEY_RAW;
    long nbytes = -1;
    enum cfx_str_format fmt = CFX_STR_FMT_HEX;
    const char *basename = NULL;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--ed25519") == 0) {
            key_type = KEY_ED25519;
        } else if (strcmp(argv[i], "--x25519") == 0) {
            key_type = KEY_X25519;
        } else if (strcmp(argv[i], "-f") == 0) {
            if (++i >= argc) {
                fprintf(stderr, "error: -f requires argument\n");
                return 1;
            }
            basename = argv[i];
        } else if (strcmp(argv[i], "-x") == 0) {
            fmt = CFX_STR_FMT_HEX;
        } else if (strcmp(argv[i], "-b64") == 0) {
            fmt = CFX_STR_FMT_BASE64;
        } else if (strcmp(argv[i], "-r") == 0) {
            fmt = CFX_STR_FMT_BINARY;
        } else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            usage(argv[0]);
            return 0;
        } else if (argv[i][0] == '-') {
            fprintf(stderr, "Unknown option: %s\n", argv[i]);
            usage(argv[0]);
            return 1;
        } else {
            char* end = NULL;
            nbytes = strtol(argv[i], &end, 10);
            if (*end != '\0' || nbytes <= 0 || nbytes > 1024 * 1024) {
                fprintf(stderr, "Invalid byte count: %s (must be 1-%d)\n", argv[i], 1024 * 1024);
                return 1;
            }
        }
    }

    switch (key_type) {
    case KEY_ED25519:
        return keygen_ed25519(basename);
    case KEY_X25519:
        return keygen_x25519(basename);
    case KEY_RAW:
        if (nbytes <= 0) {
            fprintf(stderr, "Error: specify key type (--ed25519, --x25519) or byte count\n");
            usage(argv[0]);
            return 1;
        }
        return keygen_raw(nbytes, fmt);
    }
    return 1;
}
