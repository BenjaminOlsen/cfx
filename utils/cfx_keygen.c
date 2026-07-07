/* cfx_keygen.c - cryptographic key generation utility */

#include "cfx/rand.h"
#include "cfx/ed25519.h"
#include "cfx/x25519.h"
#include "cfx/sha256.h"
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
#ifndef S_ISDIR
#define S_ISDIR(m) (((m) & _S_IFMT) == _S_IFDIR)
#endif
#else
#include <unistd.h>
#include <pwd.h>
#define mkdir_p(path) mkdir(path, 0700)
#endif

#include "cfx_cmd.h"
#include "cfx_keyfile.h"
#include "cfx_randomart.h"
#include "cfx_utils_common.h"

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
    printf("  --ed25519       Generate Ed25519 keypair (default)\n");
    printf("  --x25519        Generate X25519 keypair (encryption key)\n");
    printf("  <bytes>         Generate raw random bytes (for symmetric keys)\n\n");
    printf("Options:\n");
    printf("  -f <name>       Output file basename (default: ~/.cfx/id_<type>)\n");
    printf("  -x              Output as hex (default for raw bytes)\n");
    printf("  -b64            Output as base64\n");
    printf("  -r              Output raw bytes (binary)\n");
    printf("  --no-pw         Do not prompt for passphrase\n");
    printf("  -h, --help      Show this help\n\n");
    printf("If run with no arguments, generates Ed25519 keypair interactively.\n\n");
    printf("Examples:\n");
    printf("  %s                        Generate Ed25519 key (interactive)\n", prog);
    printf("  %s -f mykey               Generate Ed25519 key to mykey{,.pub}\n", prog);
    printf("  %s --x25519               Generate X25519 encryption key\n", prog);
    printf("  %s 32                     Generate 32-byte random key (hex)\n", prog);
}

/* prompt user for path, return 0 on success */
static int prompt_path(const char *prompt, const char *default_path, char *out, size_t outsz) {
    printf("%s (%s): ", prompt, default_path);
    fflush(stdout);

    if (!fgets(out, (int)outsz, stdin)) {
        out[0] = '\0';
    }
    /* strip newline */
    size_t len = strlen(out);
    while (len > 0 && (out[len-1] == '\n' || out[len-1] == '\r')) {
        out[--len] = '\0';
    }
    /* empty input means use default */
    if (len == 0) {
        snprintf(out, outsz, "%s", default_path);
    }
    return 0;
}

static int ct_pwd_match(const char *pw1, int pw1_len, size_t pw1_bufsz,
                        const char *pw2, int pw2_len, size_t pw2_bufsz) {
    size_t n = pw1_bufsz < pw2_bufsz ? pw1_bufsz : pw2_bufsz;
    int diff = pw1_len ^ pw2_len;
    for (size_t i = 0; i < n; ++i) {
        diff |= pw1[i] ^ pw2[i];
    }
    return diff == 0;
}

/* prompt for passphrase (twice), return length. 0 = no passphrase. */
static int prompt_passphrase(char *pwd, size_t pwdsz) {
    char pwd2[256] = {0};
    int len = cfx_key_read_secret_console("Enter passphrase (empty for no passphrase): ", pwd, pwdsz);
    if (len == 0) return 0;

    int len2 = cfx_key_read_secret_console("Enter same passphrase again: ", pwd2, sizeof(pwd2));
    if (!ct_pwd_match(pwd, len, pwdsz, pwd2, len2, sizeof(pwd2))) {
        fprintf(stderr, "Passphrases do not match.\n");
        cfx_memzero_s(pwd, pwdsz);
        cfx_memzero_s(pwd2, sizeof(pwd2));
        return -1;
    }
    cfx_memzero_s(pwd2, sizeof(pwd2));
    return len;
}

static int keygen_ed25519(const char *basename, int interactive) {
    char priv_path[1024], pub_path[1024], default_priv[1024], cfx_dir[1024];
    char pwd[256] = {0};
    int pwd_len = 0;
    int ret = 1;

    uint8_t seed[32], pk[32], sk[64];
    memset(seed, 0, sizeof(seed));
    memset(sk, 0, sizeof(sk));

    if (get_cfx_dir(cfx_dir, sizeof(cfx_dir)) != 0) goto cleanup;
    snprintf(default_priv, sizeof(default_priv), "%s/id_ed25519", cfx_dir);

    if (basename) {
        snprintf(priv_path, sizeof(priv_path), "%s", basename);
    } else if (interactive) {
        prompt_path("Enter file in which to save the key", default_priv, priv_path, sizeof(priv_path));
    } else {
        snprintf(priv_path, sizeof(priv_path), "%s", default_priv);
    }
    snprintf(pub_path, sizeof(pub_path), "%s.pub", priv_path);

    /* prompt for passphrase if on a tty */
    if (interactive) {
        pwd_len = prompt_passphrase(pwd, sizeof(pwd));
        if (pwd_len < 0) goto cleanup;
    }

    if (strncmp(priv_path, cfx_dir, strlen(cfx_dir)) == 0) {
        if (ensure_cfx_dir() != 0) goto cleanup;
    }

    /* check if files exist */
    struct stat st;
    if (stat(priv_path, &st) == 0) {
        fprintf(stderr, "%s already exists.\n", priv_path);
        if (interactive) {
            printf("Overwrite (y/n)? ");
            fflush(stdout);
            char ans[16];
            if (!fgets(ans, sizeof(ans), stdin) || (ans[0] != 'y' && ans[0] != 'Y')) {
                printf("Aborted.\n");
                goto cleanup;
            }
        } else {
            fprintf(stderr, "Use -f to specify different name, or run interactively.\n");
            goto cleanup;
        }
    }

    printf("Generating public/private ed25519 key pair.\n");

    cfx_srand_os();
    cfx_rand_bytes(seed, 32);
    cfx_ed25519_create_keypair(pk, sk, seed);

    /* write private key — encrypted or plaintext */
    if (pwd_len > 0) {
        if (cfx_key_write_encrypted(priv_path, seed, 32, pwd, (size_t)pwd_len) != 0)
            goto cleanup;
    } else {
        if (write_key_file(priv_path, seed, 32, 1) != 0)
            goto cleanup;
    }

    /* public key is always plaintext */
    if (write_key_file(pub_path, pk, 32, 0) != 0) goto cleanup;

    printf("Your identification has been saved in %s%s\n",
        priv_path, pwd_len > 0 ? " (encrypted)" : "");
    printf("Your public key has been saved in %s\n", pub_path);

    uint8_t fp[32];
    cfx_sha256(fp, pk, 32);
    printf("The key fingerprint is:\nSHA256:");
    for (int i = 0; i < 32; i++) printf("%02x", fp[i]);
    printf("\n");

    cfx_print_randomart(pk, 32, "ED25519");
    ret = 0;

cleanup:
    cfx_memzero_s(seed, sizeof(seed));
    cfx_memzero_s(sk, sizeof(sk));
    cfx_memzero_s(pwd, sizeof(pwd));
    return ret;
}

static int keygen_x25519(const char *basename, int interactive) {
    char priv_path[1024], pub_path[1024], default_priv[1024], cfx_dir[1024];
    char pwd[256] = {0};
    int pwd_len = 0;
    int ret = 1;

    uint8_t priv[32], pub[32];
    memset(priv, 0, sizeof(priv));

    if (get_cfx_dir(cfx_dir, sizeof(cfx_dir)) != 0) goto cleanup;
    snprintf(default_priv, sizeof(default_priv), "%s/id_x25519", cfx_dir);

    if (basename) {
        snprintf(priv_path, sizeof(priv_path), "%s", basename);
    } else if (interactive) {
        prompt_path("Enter file in which to save the key", default_priv, priv_path, sizeof(priv_path));
    } else {
        snprintf(priv_path, sizeof(priv_path), "%s", default_priv);
    }
    snprintf(pub_path, sizeof(pub_path), "%s.pub", priv_path);

    /* prompt for passphrase if on a tty */
    if (interactive) {
        pwd_len = prompt_passphrase(pwd, sizeof(pwd));
        if (pwd_len < 0) goto cleanup;
    }

    if (strncmp(priv_path, cfx_dir, strlen(cfx_dir)) == 0) {
        if (ensure_cfx_dir() != 0) goto cleanup;
    }

    struct stat st;
    if (stat(priv_path, &st) == 0) {
        fprintf(stderr, "%s already exists.\n", priv_path);
        if (interactive) {
            printf("Overwrite (y/n)? ");
            fflush(stdout);
            char ans[16];
            if (!fgets(ans, sizeof(ans), stdin) || (ans[0] != 'y' && ans[0] != 'Y')) {
                printf("Aborted.\n");
                goto cleanup;
            }
        } else {
            fprintf(stderr, "Use -f to specify different name, or run interactively.\n");
            goto cleanup;
        }
    }

    printf("Generating public/private x25519 key pair.\n");

    cfx_srand_os();
    cfx_rand_bytes(priv, 32);
    cfx_x25519_base(pub, priv);

    /* write private key encrypted or plaintext */
    if (pwd_len > 0) {
        if (cfx_key_write_encrypted(priv_path, priv, 32, pwd, (size_t)pwd_len) != 0)
            goto cleanup;
    } else {
        if (write_key_file(priv_path, priv, 32, 1) != 0)
            goto cleanup;
    }

    /* public key is always plaintext */
    if (write_key_file(pub_path, pub, 32, 0) != 0) goto cleanup;

    printf("Your identification has been saved in %s%s\n", priv_path, pwd_len > 0 ? " (encrypted)" : "");
    printf("Your public key has been saved in %s\n", pub_path);

    uint8_t fp[32];
    cfx_sha256(fp, pub, 32);
    printf("The key fingerprint is:\nSHA256:");
    for (int i = 0; i < 32; i++) printf("%02x", fp[i]);
    printf("\n");

    cfx_print_randomart(pub, 32, "X25519");
    ret = 0;

cleanup:
    cfx_memzero_s(priv, sizeof(priv));
    cfx_memzero_s(pwd, sizeof(pwd));
    return ret;
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
        char *b64 = malloc(b64_len + 1);
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
    key_type_t key_type = KEY_ED25519;  /* default to ed25519 */
    long nbytes = -1;
    enum cfx_str_format fmt = CFX_STR_FMT_HEX;
    const char *basename = NULL;
    int no_pw = 0;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--ed25519") == 0) {
            key_type = KEY_ED25519;
        } else if (strcmp(argv[i], "--x25519") == 0) {
            key_type = KEY_X25519;
        } else if (strcmp(argv[i], "--no-pw") == 0) {
            no_pw = 1;
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
            key_type = KEY_RAW;
        }
    }

    int interactive = !no_pw && isatty(fileno(stdin));

    switch (key_type) {
    case KEY_ED25519:
        return keygen_ed25519(basename, interactive);
    case KEY_X25519:
        return keygen_x25519(basename, interactive);
    case KEY_RAW:
        if (nbytes <= 0) {
            fprintf(stderr, "Error: specify byte count for raw key\n");
            usage(argv[0]);
            return 1;
        }
        return keygen_raw(nbytes, fmt);
    }
    return 1;
}
