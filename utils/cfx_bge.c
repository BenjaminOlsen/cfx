/* cfx_bge.c -- passphrase-protected secrets store */

#include "cfx/argon2.h"
#include "cfx/aead_chacha20_poly1305.h"
#include "cfx/rand.h"
#include "cfx/memory.h"
#include "cfx/base64.h"

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
#include <fcntl.h>
#include <sys/file.h>
#include <pwd.h>
#include <termios.h>
#define mkdir_p(path) mkdir(path, 0700)
#endif

#include "cfx_cmd.h"
#include "cfx_keyfile.h"
#include "misc.h"
#include "cfx/macros.h"

/*
 * BGE v2 file layout:
 *
 *   [0..2]    BGE_MAGIC "BGE"
 *   [3]       BGE_VERSION
 *   [4..15]   argon2 params (m, t, p) LE32 each
 *   [16..31]  salt
 *   [32..55]  nonce
 *   [56..71]  key verifier
 *   [72..N]   ciphertext
 *   [last 16] Poly1305 tag
 *
 * argon2 output is 48 bytes: key(32) | verifier(16), one KDF call.
 * AAD = bytes 0-71 (hdr + verifier).  min file = 88 bytes.
 *
 * plaintext is length-prefixed KV pairs:
 *   entry = u16le(name_len) | name | u32le(val_len) | val
 */
#define BGE_MAGIC         "BGE"
#define BGE_VERSION       2       /* KV store format version */
#define BGE_STREAM_VERSION 3      /* streaming file encryption */
#define BGE_VERSION_STR   "2.2.0"
#define BGE_HEADER_LEN    56
#define BGE_VERIFIER_LEN  16
#define BGE_AAD_LEN       (BGE_HEADER_LEN + BGE_VERIFIER_LEN)  /* 72 */
#define BGE_TAG_LEN       16
#define BGE_MIN_FILE      (BGE_AAD_LEN + BGE_TAG_LEN)          /* 88: empty store */

typedef struct {
    uint8_t  magic[3];
    uint8_t  version;
    uint32_t m_cost;    /* little-endian on disk */
    uint32_t t_cost;    /* little-endian on disk */
    uint32_t p_cost;    /* little-endian on disk */
    uint8_t  salt[16];
    uint8_t  nonce[24];
} bge_header;
CFX_STATIC_ASSERT(sizeof(bge_header) == BGE_HEADER_LEN, bge_header_packing);

typedef struct {
    bge_header hdr;
    uint8_t    key[32];
    uint8_t    verifier[BGE_VERIFIER_LEN];
    uint8_t   *file_buf;
    size_t     file_len;
} bge_store;

static void bge_store_wipe(bge_store *s) {
    cfx_memzero_s(s->key, sizeof(s->key));
    cfx_memzero_s(s->verifier, sizeof(s->verifier));
    cfx_memzero_s(&s->hdr, sizeof(s->hdr));
    if (s->file_buf) {
        cfx_memzero_s(s->file_buf, s->file_len);
        free(s->file_buf);
        s->file_buf = NULL;
    }
    s->file_len = 0;
}

#define BGE_DEFAULT_M   0x10000  /* 64 MB (in KB) */
#define BGE_DEFAULT_T   3
#define BGE_DEFAULT_P   4

#define BGE_READ_MAX    (128u * 1024u * 1024u)  /* 16 MB */

#define BGE_ARMOR_BEGIN "-----BEGIN BGE MESSAGE-----\n"
#define BGE_ARMOR_END   "\n-----END BGE MESSAGE-----\n"

static int bge_is_armored(const uint8_t *buf, size_t len) {
    return len >= 28 && memcmp(buf, "-----BEGIN BGE MESSAGE-----", 27) == 0;
}

/* wrap binary blob in PEM-style armor with 76-char lines. caller frees *out. */
static int bge_armor_encode(const uint8_t *bin, size_t bin_len,
                            uint8_t **out, size_t *out_len) {
    size_t b64_len = 0;
    cfx_base64_encode(NULL, &b64_len, bin, bin_len);

    char *b64 = malloc(b64_len);
    if (!b64) return -1;
    cfx_base64_encode(b64, &b64_len, bin, bin_len);

    /* header + base64 + newlines every 76 chars + footer */
    size_t nlines = (b64_len + 75) / 76;
    size_t total = strlen(BGE_ARMOR_BEGIN) + b64_len + nlines + strlen(BGE_ARMOR_END);
    uint8_t *buf = malloc(total + 1);
    if (!buf) { free(b64); return -1; }

    uint8_t *w = buf;
    size_t hlen = strlen(BGE_ARMOR_BEGIN);
    memcpy(w, BGE_ARMOR_BEGIN, hlen); w += hlen;

    for (size_t i = 0; i < b64_len; i += 76) {
        size_t chunk = b64_len - i;
        if (chunk > 76) chunk = 76;
        memcpy(w, b64 + i, chunk); w += chunk;
        *w++ = '\n';
    }

    size_t flen = strlen(BGE_ARMOR_END);
    memcpy(w, BGE_ARMOR_END, flen); w += flen;

    free(b64);
    *out = buf;
    *out_len = (size_t)(w - buf);
    return 0;
}

/* strip PEM header/footer and base64-decode. caller frees *out. */
static int bge_armor_decode(const uint8_t *text, size_t text_len,
                            uint8_t **out, size_t *out_len) {
    /* find start after first newline past header */
    const char *s = (const char *)text;
    const char *body = memchr(s, '\n', text_len);
    if (!body) return -1;
    body++;

    /* find footer */
    const char *footer = strstr(body, "-----END BGE MESSAGE-----");
    if (!footer) return -1;

    size_t b64_len = (size_t)(footer - body);

    /* decode */
    size_t dec_len = cfx_base64_dec_max_len(b64_len);
    uint8_t *dec = malloc(dec_len);
    if (!dec) return -1;

    int rc = cfx_base64_decode(dec, &dec_len, body, b64_len);
    if (rc != 0) { free(dec); return -1; }

    *out = dec;
    *out_len = dec_len;
    return 0;
}

static void usage(const char *prog) {
    printf("Usage:\n");
    printf("  %s <command> [args...]              Secrets store operations\n", prog);
    printf("  %s -e [-a] [-i INPUT] [-o OUTPUT]    Encrypt a file\n", prog);
    printf("  %s -d [-i INPUT] [-o OUTPUT]         Decrypt a file\n\n", prog);
    printf("Passphrase-protected secrets store and file encryption using\n");
    printf("Argon2id key derivation and XChaCha20-Poly1305 AEAD.\n\n");
    printf("Store commands:\n");
    printf("  init   [-s path]                 Create a new encrypted store\n");
    printf("  get    <name> [-s path] [-q]     Retrieve a secret (raw to stdout)\n");
    printf("  set    [name] [value] [-s path]  Set a secret (value from arg, -i, or stdin)\n");
    printf("  rm     <name> [-s path]          Remove a secret\n");
    printf("  ls     [-s path]                 List all secret names\n");
    printf("  info   [-s path]                 Show store location, size, and entry count\n");
    printf("  passwd [-s path]                 Change passphrase\n");
    printf("  dump   [-s path]                 Print all secrets to stdout\n");
    printf("\nFile encryption:\n");
    printf("  -e, --encrypt    Encrypt input to output\n");
    printf("  -d, --decrypt    Decrypt input to output\n");
    printf("  -i, --input      Read from file (default: stdin)\n");
    printf("  -o, --output     Write to file (default: stdout)\n");
    printf("  -a, --armor      PEM-encoded base64 (encrypt only; decrypt auto-detects)\n");
    printf("\nCommon options:\n");
    printf("  -p, --passphrase <pw>  Supply passphrase on command line (default: prompt)\n");
    printf("\nStore options:\n");
    printf("  -s, --store <path>   Path to BGE store (default: ~/.cfx/secrets.bge)\n");
    printf("\nOptions for get:\n");
    printf("  -q, --quoted         Wrap output in single quotes\n");
    printf("\nOptions for set:\n");
    printf("  -i <file>   Read value from file (binary-safe)\n");
    printf("              If no value and no -i, reads from stdin.\n");
    printf("  -f, --force  Overwrite existing entry without asking\n");
    printf("\nExamples:\n");
    printf("  %s -e -i doc.txt -o secret.bge        Encrypt a file\n", prog);
    printf("  %s -e -a -i doc.txt -o secret.asc     Encrypt to armored output\n", prog);
    printf("  %s -d -i secret.bge -o doc.txt        Decrypt a file\n", prog);
    printf("  echo data | %s -e | %s -d             Roundtrip via pipe\n", prog, prog);
    printf("  %s init                              Create store at default path\n", prog);
    printf("  %s set token mytoken                 Set from command line argument\n", prog);
    printf("  %s get token                         Print value (newline on tty)\n", prog);
    printf("  %s ls                                Show all key names\n", prog);
}

#define SUBDIR ".cfx"

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
    int n = snprintf(buf, bufsz, "%s/%s", home, SUBDIR);
    if (n < 0 || (size_t)n >= bufsz) {
        fprintf(stderr, "error: path too long\n");
        return -1;
    }
    return 0;
}

static int bge_default_path(char *buf, size_t bufsz) {
    char cfx_dir[1024];
    if (get_cfx_dir(cfx_dir, sizeof(cfx_dir)) != 0) return -1;
    int n = snprintf(buf, bufsz, "%s/secrets.bge", cfx_dir);
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


/* write entire FILE* into malloc'd buf; works for pipes too. capped at BGE_READ_MAX. */
static int bge_read_all(FILE *f, uint8_t **out, size_t *out_len) {
    uint8_t *buf = NULL;
    size_t cap = 0, len = 0;

    if (!f || !out || !out_len) return -1;
    *out = NULL;
    *out_len = 0;

    for (;;) {
        if (len == cap) {
            size_t newcap = cap ? (cap * 2) : 4096;
            if (newcap > BGE_READ_MAX) newcap = BGE_READ_MAX;
            if (newcap <= cap) {
                if (buf) { cfx_memzero_s(buf, cap); free(buf); }
                fprintf(stderr, "error: input exceeds 16 MB limit\n");
                return -1;
            }
            uint8_t *tmp = realloc(buf, newcap);
            if (!tmp) {
                if (buf) { cfx_memzero_s(buf, cap); free(buf); }
                return -1;
            }
            buf = tmp;
            cap = newcap;
        }

        size_t r = fread(buf + len, 1, cap - len, f);
        len += r;
        if (r == 0) {
            if (feof(f)) break;
            if (buf) {
                cfx_memzero_s(buf, cap);
                free(buf);
            }
            return -1;
        }
    }

    *out = buf;
    *out_len = len;
    return 0;
}

/* read password from /dev/tty (keeps stdin free for pipes).
 * falls back to cfx_key_read_secret_console on win32 or if tty unavailable. */
static int bge_read_secret(const char *prompt, char *buf, size_t bufsz) {
#ifndef _WIN32
    FILE *tty = fopen("/dev/tty", "r+");
    if (tty) {
        fprintf(tty, "%s", prompt);
        fflush(tty);

        /* kill echo */
        int tty_fd = fileno(tty);
        struct termios old, noecho;
        int have_termios = (tcgetattr(tty_fd, &old) == 0);
        if (have_termios) {
            noecho = old;
            noecho.c_lflag &= ~(tcflag_t)ECHO;
            tcsetattr(tty_fd, TCSANOW, &noecho);
        }

        char *r = fgets(buf, (int)bufsz, tty);

        if (have_termios) {
            tcsetattr(tty_fd, TCSANOW, &old);
            fprintf(tty, "\n");
        }

        fclose(tty);
        if (!r) return 0;
        int len = (int)strlen(buf);
        while (len > 0 && (buf[len-1] == '\n' || buf[len-1] == '\r')) {
            buf[--len] = '\0';
        }
        return len;
    }
#endif
    return cfx_key_read_secret_console(prompt, buf, bufsz);
}

static const char *g_passphrase_arg;

static int bge_read_passphrase(const char *prompt, char *buf, size_t bufsz) {
    if (g_passphrase_arg) {
        size_t len = strlen(g_passphrase_arg);
        if (len >= bufsz) len = bufsz - 1;
        memcpy(buf, g_passphrase_arg, len);
        buf[len] = '\0';
        return (int)len;
    }
    return bge_read_secret(prompt, buf, bufsz);
}

/* like bge_read_secret but with echo on (visible input). */
static int bge_read_visible(const char *prompt, char *buf, size_t bufsz) {
#ifndef _WIN32
    FILE *tty = fopen("/dev/tty", "r+");
    if (tty) {
        fprintf(tty, "%s", prompt);
        fflush(tty);

        char *r = fgets(buf, (int)bufsz, tty);
        fclose(tty);
        if (!r) return 0;
        int len = (int)strlen(buf);
        while (len > 0 && (buf[len-1] == '\n' || buf[len-1] == '\r')) {
            buf[--len] = '\0';
        }
        return len;
    }
#endif
    fprintf(stderr, "%s", prompt);
    fflush(stderr);
    char *r = fgets(buf, (int)bufsz, stdin);
    if (!r) return 0;
    int len = (int)strlen(buf);
    while (len > 0 && (buf[len-1] == '\n' || buf[len-1] == '\r')) {
        buf[--len] = '\0';
    }
    return len;
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

/* double-prompt passphrase, returns length or -1 */
static int prompt_passphrase(char *pwd, size_t pwdsz) {
    if (g_passphrase_arg) {
        size_t len = strlen(g_passphrase_arg);
        if (len >= pwdsz) len = pwdsz - 1;
        memcpy(pwd, g_passphrase_arg, len);
        pwd[len] = '\0';
        return (int)len;
    }
    char pwd2[256] = {0};
    int len = bge_read_secret("Enter passphrase: ", pwd, pwdsz);
    if (len <= 0) {
        fprintf(stderr, "error: passphrase required\n");
        cfx_memzero_s(pwd, pwdsz);
        return -1;
    }

    int len2 = bge_read_secret("Enter same passphrase again: ", pwd2, sizeof(pwd2));
    if (!ct_pwd_match(pwd, len, pwdsz, pwd2, len2, sizeof(pwd2))) {
        fprintf(stderr, "Passphrases do not match.\n");
        cfx_memzero_s(pwd, pwdsz);
        cfx_memzero_s(pwd2, sizeof(pwd2));
        return -1;
    }
    cfx_memzero_s(pwd2, sizeof(pwd2));
    return len;
}


/* read file, run argon2 -> key+verifier, check verifier matches.
 * fills store with key, verifier, hdr, raw file buf for later decrypt. */
static int bge_authenticate(const char *path, const char *pwd, size_t pwd_len,
                             bge_store *store) {
    FILE *f = fopen(path, "rb");
    if (!f) {
        fprintf(stderr, "error: cannot open %s: %s\n", path, strerror(errno));
        return -1;
    }

    uint8_t *file_buf = NULL;
    size_t file_len = 0;
    if (bge_read_all(f, &file_buf, &file_len) != 0) {
        fclose(f);
        fprintf(stderr, "error: cannot read %s\n", path);
        return -1;
    }
    fclose(f);

    if (file_len < BGE_MIN_FILE ||
        memcmp(file_buf, BGE_MAGIC, 3) != 0 || file_buf[3] != BGE_VERSION) {
        fprintf(stderr, "error: %s is not a valid BGE store\n", path);
        cfx_memzero_s(file_buf, file_len);
        free(file_buf);
        return -1;
    }

    bge_header hdr;
    memcpy(&hdr, file_buf, sizeof(hdr));

    uint32_t m_cost = cfx_load32_le(&hdr.m_cost);
    uint32_t t_cost = cfx_load32_le(&hdr.t_cost);
    uint32_t p      = cfx_load32_le(&hdr.p_cost);

    if (m_cost < 8 || t_cost < 1 || p < 1 ||
        m_cost > 4194304 || t_cost > 100 || p > 16) {
        fprintf(stderr, "error: unreasonable argon2 parameters in %s\n", path);
        cfx_memzero_s(file_buf, file_len);
        free(file_buf);
        return -1;
    }

    /* one argon2 call -> 48 bytes: key(32) | verifier(16) */
    uint8_t kdf_out[48];
    int rc = cfx_argon2id(kdf_out, 48,
        (const uint8_t *)pwd, pwd_len, hdr.salt, sizeof(hdr.salt),
        m_cost, t_cost, p);
    if (rc != 0) {
        fprintf(stderr, "error: argon2 key derivation failed\n");
        cfx_memzero_s(kdf_out, sizeof(kdf_out));
        cfx_memzero_s(file_buf, file_len);
        free(file_buf);
        return -1;
    }

    /* ct compare derived vs stored verifier (offset 56) */
    const uint8_t *stored_verifier = file_buf + BGE_HEADER_LEN;
    uint8_t diff = 0;
    for (int i = 0; i < BGE_VERIFIER_LEN; i++)
        diff |= kdf_out[32 + i] ^ stored_verifier[i];

    if (diff != 0) {
        fprintf(stderr, "error: wrong passphrase\n");
        cfx_memzero_s(kdf_out, sizeof(kdf_out));
        cfx_memzero_s(file_buf, file_len);
        free(file_buf);
        return -1;
    }

    memcpy(store->key, kdf_out, 32);
    memcpy(store->verifier, kdf_out + 32, BGE_VERIFIER_LEN);
    store->hdr = hdr;
    store->file_buf = file_buf;
    store->file_len = file_len;

    cfx_memzero_s(kdf_out, sizeof(kdf_out));
    return 0;
}

/* AEAD decrypt using already-authenticated store. caller frees *pt_out. */
static int bge_decrypt_store(const bge_store *store,
                              uint8_t **pt_out, size_t *pt_len) {
    const uint8_t *nonce = store->hdr.nonce;
    size_t ct_len = store->file_len - BGE_AAD_LEN - BGE_TAG_LEN;
    const uint8_t *ct    = store->file_buf + BGE_AAD_LEN;
    const uint8_t *tag   = store->file_buf + store->file_len - BGE_TAG_LEN;

    uint8_t *plaintext = malloc(ct_len + 1); /* +1 for NUL convenience */
    if (!plaintext) {
        fprintf(stderr, "error: allocation failed\n");
        return -1;
    }

    /* AAD = first 72 bytes (hdr + verifier) */
    int rc = cfx_xchacha20_poly1305_decrypt(
        plaintext, ct, ct_len,
        store->file_buf, BGE_AAD_LEN,
        store->key, nonce, tag);

    if (rc != 0) {
        fprintf(stderr, "error: store corrupted or tampered\n");
        cfx_memzero_s(plaintext, ct_len + 1);
        free(plaintext);
        return -1;
    }

    plaintext[ct_len] = '\0';

    *pt_out = plaintext;
    *pt_len = ct_len;
    return 0;
}

/* atomic write: flock -> write .tmp (0600) -> fsync -> rename */
static int bge_safe_write(const char *path,
                          const bge_header *header,
                          const uint8_t verifier[BGE_VERIFIER_LEN],
                          const uint8_t *ct, size_t ct_len,
                          const uint8_t *tag) {
    char tmppath[1088];
    int n = snprintf(tmppath, sizeof(tmppath), "%s.tmp", path);
    if (n < 0 || (size_t)n >= sizeof(tmppath)) {
        fprintf(stderr, "error: path too long\n");
        return -1;
    }

    int lock_fd = -1;
    int ret = -1;

#ifndef _WIN32
    /* grab flock if file exists */
    lock_fd = open(path, O_RDONLY);
    if (lock_fd >= 0) {
        if (flock(lock_fd, LOCK_EX) != 0) {
            fprintf(stderr, "error: cannot lock %s: %s\n", path, strerror(errno));
            close(lock_fd);
            return -1;
        }
    }

    /* open .tmp with 0600 from the start -- no chmod race */
    int fd = open(tmppath, O_CREAT | O_WRONLY | O_TRUNC, 0600);
    if (fd < 0) {
        fprintf(stderr, "error: cannot create %s: %s\n", tmppath, strerror(errno));
        if (lock_fd >= 0) close(lock_fd);
        return -1;
    }
    FILE *f = fdopen(fd, "wb");
    if (!f) {
        close(fd);
        unlink(tmppath);
        if (lock_fd >= 0) close(lock_fd);
        return -1;
    }
#else
    FILE *f = fopen(tmppath, "wb");
    if (!f) {
        fprintf(stderr, "error: cannot create %s: %s\n", tmppath, strerror(errno));
        return -1;
    }
#endif

    size_t ok = 1;
    ok = ok && fwrite(header, 1, sizeof(*header), f) == sizeof(*header);
    ok = ok && fwrite(verifier, 1, BGE_VERIFIER_LEN, f) == BGE_VERIFIER_LEN;

    if (ct_len > 0) {
        ok = ok && fwrite(ct, 1, ct_len, f) == ct_len;
    }

    ok = ok && fwrite(tag, 1, BGE_TAG_LEN, f) == BGE_TAG_LEN;

    if (!ok) {
        fprintf(stderr, "error: write failed: %s\n", strerror(errno));
        fclose(f);
        unlink(tmppath);
        if (lock_fd >= 0) close(lock_fd);
        return -1;
    }

#ifndef _WIN32
    fsync(fileno(f));
#else
    fflush(f);
#endif
    fclose(f);

    if (rename(tmppath, path) != 0) {
        fprintf(stderr, "error: rename failed: %s\n", strerror(errno));
        unlink(tmppath);
        if (lock_fd >= 0) close(lock_fd);
        return -1;
    }

    ret = 0;

    if (lock_fd >= 0) close(lock_fd);
    return ret;
}

/* fresh salt + nonce, derive key+verifier, encrypt, write */
static int bge_encrypt_write(const char *path, const uint8_t *pt, size_t pt_len,
                             const char *pwd, size_t pwd_len,
                             uint32_t m, uint32_t t, uint32_t p) {
    bge_header header;
    uint8_t kdf_out[48];
    uint8_t verifier[BGE_VERIFIER_LEN];
    int ret = -1;

    memcpy(header.magic, BGE_MAGIC, 3);
    header.version = BGE_VERSION;
    cfx_store32_le(&header.m_cost, m);
    cfx_store32_le(&header.t_cost, t);
    cfx_store32_le(&header.p_cost, p);

    cfx_srand_os();
    cfx_rand_bytes(header.salt,  sizeof(header.salt));
    cfx_rand_bytes(header.nonce, sizeof(header.nonce));

    /* KDF -> 48 bytes: key | verifier */
    int rc = cfx_argon2id(kdf_out, 48,
        (const uint8_t *)pwd, pwd_len,
        header.salt, sizeof(header.salt), m, t, p);
    if (rc != 0) {
        fprintf(stderr, "error: argon2 key derivation failed\n");
        goto done;
    }
    memcpy(verifier, kdf_out + 32, BGE_VERIFIER_LEN);

    /* AAD = hdr(56) + verifier(16) */
    uint8_t aad[BGE_AAD_LEN];
    memcpy(aad, &header, sizeof(header));
    memcpy(aad + sizeof(header), verifier, BGE_VERIFIER_LEN);

    uint8_t *ct = malloc(pt_len);
    uint8_t tag[BGE_TAG_LEN];
    if (pt_len > 0 && !ct) {
        fprintf(stderr, "error: allocation failed\n");
        goto done;
    }

    rc = cfx_xchacha20_poly1305_encrypt(
        ct, tag, pt, pt_len,
        aad, BGE_AAD_LEN,
        kdf_out, header.nonce);
    if (rc != 0) {
        fprintf(stderr, "error: encryption failed\n");
        free(ct);
        goto done;
    }

    ret = bge_safe_write(path, &header, verifier, ct, pt_len, tag);

    cfx_memzero_s(ct, pt_len);
    free(ct);

done:
    cfx_memzero_s(kdf_out, sizeof(kdf_out));
    cfx_memzero_s(verifier, sizeof(verifier));
    cfx_memzero_s(&header, sizeof(header));
    return ret;
}

/* re-encrypt with same key, fresh nonce. avoids re-running argon2. */
static int bge_write_reusing_key(const char *path, const uint8_t *pt, size_t pt_len,
                                 const bge_store *store) {
    bge_header header = store->hdr;
    int ret = -1;

    /* fresh nonce, keep magic/params/salt from old header */
    cfx_srand_os();
    cfx_rand_bytes(header.nonce, sizeof(header.nonce));

    /* AAD = hdr(56) + verifier(16) */
    uint8_t aad[BGE_AAD_LEN];
    memcpy(aad, &header, sizeof(header));
    memcpy(aad + sizeof(header), store->verifier, BGE_VERIFIER_LEN);

    uint8_t *ct = malloc(pt_len);
    uint8_t tag[BGE_TAG_LEN];
    if (pt_len > 0 && !ct) {
        fprintf(stderr, "error: allocation failed\n");
        return -1;
    }

    int rc = cfx_xchacha20_poly1305_encrypt(
        ct, tag, pt, pt_len,
        aad, BGE_AAD_LEN,
        store->key, header.nonce);

    if (rc != 0) {
        fprintf(stderr, "error: encryption failed\n");
        goto done;
    }

    ret = bge_safe_write(path, &header, store->verifier, ct, pt_len, tag);

done:
    cfx_memzero_s(ct, pt_len);
    free(ct);
    cfx_memzero_s(&header, sizeof(header));
    return ret;
}


/* check if string is a positive integer (for index lookups) */
static int is_numeric(const char *s) {
    if (!s || !*s) return 0;
    for (; *s; s++) {
        if (*s < '0' || *s > '9') return 0;
    }
    return 1;
}

/* resolve 1-based index to entry name. NULL if out of range. */
static const uint8_t *store_name_by_index(const uint8_t *pt, size_t pt_len,
                                           unsigned idx, uint16_t *name_len) {
    const uint8_t *p = pt;
    const uint8_t *end = pt + pt_len;
    unsigned cur = 0;

    while (p + 2 <= end) {
        uint16_t klen = cfx_load16_le(p);
        p += 2;
        if (p + klen > end) break;
        const uint8_t *kptr = p;
        p += klen;
        if (p + 4 > end) break;
        uint32_t vl = cfx_load32_le(p);
        p += 4;
        if (p + vl > end) break;
        p += vl;

        if (++cur == idx) {
            *name_len = klen;
            return kptr;
        }
    }
    return NULL;
}

/* find entry by name, return ptr to value + length. NULL if missing.
    vlen can be null if you don't care about the length */
static const uint8_t *store_get(const uint8_t *pt, size_t pt_len,
                                const char *name, size_t *vlen) {
    size_t nlen = strlen(name);
    const uint8_t *p = pt;
    const uint8_t *end = pt + pt_len;

    while (p + 2 <= end) {
        uint16_t klen = cfx_load16_le(p);
        p += 2;
        if (p + klen > end) break;
        const uint8_t *kptr = p;
        p += klen;

        if (p + 4 > end) break;
        uint32_t val_len = cfx_load32_le(p);
        p += 4;
        if (p + val_len > end) break;
        const uint8_t *vptr = p;
        p += val_len;

        if (klen == nlen && memcmp(kptr, name, nlen) == 0) {
            if (vlen) *vlen = val_len;
            return vptr;
        }
    }
    return NULL;
}

/* upsert entry. returns new malloc'd store buf, sets *new_len. */
static uint8_t *store_set(const uint8_t *pt, size_t pt_len, const char *name, const uint8_t *val, size_t val_len,
                          size_t *new_len) {

    size_t nlen = strlen(name);

    /* size of new entry: 2 + nlen + 4 + val_len */
    size_t entry_len = 2 + nlen + 4 + val_len;

    /* scan for existing entry */
    const uint8_t *p = pt;
    const uint8_t *end = pt + pt_len;
    const uint8_t *found_start = NULL;
    size_t found_entry_len = 0;

    while (p + 2 <= end) {
        const uint8_t *entry_start = p;
        uint16_t klen = cfx_load16_le(p);
        p += 2;
        if (p + klen > end) break;
        const uint8_t *kptr = p;
        p += klen;

        if (p + 4 > end) break;
        uint32_t vl = cfx_load32_le(p);
        p += 4;
        if (p + vl > end) break;
        p += vl;

        if (klen == nlen && memcmp(kptr, name, nlen) == 0) {
            found_start = entry_start;
            found_entry_len = (size_t)(p - entry_start);
            break;
        }
    }

    size_t out_len;
    uint8_t *out;

    if (found_start) {
        /* replace in-place */
        size_t prefix = (size_t)(found_start - pt);
        size_t suffix_off = prefix + found_entry_len;
        size_t suffix_len = pt_len - suffix_off;
        out_len = prefix + entry_len + suffix_len;
        out = malloc(out_len);
        if (!out) return NULL;

        memcpy(out, pt, prefix);
        uint8_t *w = out + prefix;
        cfx_store16_le(w, (uint16_t)nlen); w += 2;
        memcpy(w, name, nlen); w += nlen;
        cfx_store32_le(w, (uint32_t)val_len); w += 4;
        memcpy(w, val, val_len); w += val_len;
        memcpy(w, pt + suffix_off, suffix_len);
    } else {
        /* append */
        out_len = pt_len + entry_len;
        out = malloc(out_len);
        if (!out) return NULL;

        memcpy(out, pt, pt_len);
        uint8_t *w = out + pt_len;
        cfx_store16_le(w, (uint16_t)nlen); w += 2;
        memcpy(w, name, nlen); w += nlen;
        cfx_store32_le(w, (uint32_t)val_len); w += 4;
        memcpy(w, val, val_len);
    }

    *new_len = out_len;
    return out;
}

/* remove entry by name. returns new malloc'd buf or NULL if not found. */
static uint8_t *store_rm(const uint8_t *pt, size_t pt_len, const char *name, size_t *new_len) {
    size_t nlen = strlen(name);
    const uint8_t *p = pt;
    const uint8_t *end = pt + pt_len;
    const uint8_t *found_start = NULL;
    size_t found_entry_len = 0;

    while (p + 2 <= end) {
        const uint8_t *entry_start = p;
        uint16_t klen = cfx_load16_le(p);
        p += 2;
        if (p + klen > end) break;
        const uint8_t *kptr = p;
        p += klen;

        if (p + 4 > end) break;
        uint32_t vl = cfx_load32_le(p);
        p += 4;
        if (p + vl > end) break;
        p += vl;

        if (klen == nlen && memcmp(kptr, name, nlen) == 0) {
            found_start = entry_start;
            found_entry_len = (size_t)(p - entry_start);
            break;
        }
    }

    if (!found_start) {
        fprintf(stderr, "error: '%s' not found\n", name);
        return NULL;
    }

    size_t prefix = (size_t)(found_start - pt);
    size_t suffix_off = prefix + found_entry_len;
    size_t suffix_len = pt_len - suffix_off;
    size_t out_len = prefix + suffix_len;

    uint8_t *out = malloc(out_len + 1);
    if (!out) return NULL;

    memcpy(out, pt, prefix);
    memcpy(out + prefix, pt + suffix_off, suffix_len);
    *new_len = out_len;
    return out;
}

/* for dump - render as "[name]\nvalue\n" pairs. malloc'd, caller frees. */
static uint8_t *store_to_text(const uint8_t *pt, size_t pt_len, size_t *text_len) {
    /* pass 1: size */
    size_t total = 0;
    const uint8_t *p = pt;
    const uint8_t *end = pt + pt_len;

    while (p + 2 <= end) {
        uint16_t klen = cfx_load16_le(p);
        p += 2;
        if (p + klen > end) break;
        p += klen;
        if (p + 4 > end) break;
        uint32_t vl = cfx_load32_le(p);
        p += 4;
        if (p + vl > end) break;
        p += vl;

        /* [name]\n + value + \n */
        total += 1 + klen + 2 + vl + 1;
    }

    uint8_t *out = malloc(total + 1);
    if (!out) return NULL;

    /* pass 2: write */
    uint8_t *w = out;
    p = pt;
    while (p + 2 <= end) {
        uint16_t klen = cfx_load16_le(p);
        p += 2;
        if (p + klen > end) break;
        const uint8_t *kptr = p;
        p += klen;
        if (p + 4 > end) break;
        uint32_t vl = cfx_load32_le(p);
        p += 4;
        if (p + vl > end) break;
        const uint8_t *vptr = p;
        p += vl;

        *w++ = '[';
        memcpy(w, kptr, klen); w += klen;
        *w++ = ']';
        *w++ = '\n';
        memcpy(w, vptr, vl); w += vl;
        *w++ = '\n';
    }

    *w = '\0';
    *text_len = (size_t)(w - out);
    return out;
}

static int bge_init(int argc, char **argv) {
    char path_buf[1024];
    const char *path = NULL;

    for (int i = 2; i < argc; i++) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            printf("Usage: %s init [-s path]\n", argv[0]);
            return 0;
        } else if (strcmp(argv[i], "-s") == 0 || strcmp(argv[i], "--store") == 0) {
            if (i + 1 >= argc) { fprintf(stderr, "error: -s requires a path\n"); return 1; }
            path = argv[++i];
        } else {
            fprintf(stderr, "error: unknown argument: %s\n", argv[i]);
            return 1;
        }
    }

    if (!path) {
        if (bge_default_path(path_buf, sizeof(path_buf)) != 0) return 1;
        path = path_buf;
    }

    /* mkdir if default path */
    char cfx_dir[1024];
    if (get_cfx_dir(cfx_dir, sizeof(cfx_dir)) == 0 &&
        strncmp(path, cfx_dir, strlen(cfx_dir)) == 0) {
        if (ensure_cfx_dir() != 0) return 1;
    }

    /* bail if already exists */
    struct stat st;
    if (stat(path, &st) == 0) {
        fprintf(stderr, "error: %s already exists\n", path);
        return 1;
    }

    char pwd[256] = {0};
    int pwd_len = prompt_passphrase(pwd, sizeof(pwd));
    if (pwd_len < 0) return 1;

    int rc = bge_encrypt_write(path, NULL, 0, pwd, (size_t)pwd_len,
                               BGE_DEFAULT_M, BGE_DEFAULT_T, BGE_DEFAULT_P);
    cfx_memzero_s(pwd, sizeof(pwd));

    if (rc == 0) {
        printf("Initialized encrypted store: %s\n", path);
    }
    return rc != 0;
}

static int bge_get(int argc, char **argv) {
    const char *name = NULL;
    const char *path = NULL;
    int quoted = 0;
    char path_buf[1024];

    for (int i = 2; i < argc; i++) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            printf("Usage: %s get <name> [-s path] [-q]\n", argv[0]);
            return 0;
        } else if (strcmp(argv[i], "-s") == 0 || strcmp(argv[i], "--store") == 0) {
            if (i + 1 >= argc) { fprintf(stderr, "error: -s requires a path\n"); return 1; }
            path = argv[++i];
        } else if (strcmp(argv[i], "-q") == 0 || strcmp(argv[i], "--quoted") == 0) {
            quoted = 1;
        } else if (!name) {
            name = argv[i];
        } else {
            fprintf(stderr, "error: unknown argument: %s\n", argv[i]);
            return 1;
        }
    }


    if (!path) {
        if (bge_default_path(path_buf, sizeof(path_buf)) != 0) return 1;
        path = path_buf;
    }

    char pwd[256] = {0};
    int pwd_len = bge_read_passphrase("Enter passphrase: ", pwd, sizeof(pwd));
    if (pwd_len <= 0) {
        fprintf(stderr, "error: passphrase required\n");
        cfx_memzero_s(pwd, sizeof(pwd));
        return 1;
    }

    bge_store store = {0};
    int rc = bge_authenticate(path, pwd, (size_t)pwd_len, &store);
    cfx_memzero_s(pwd, sizeof(pwd));
    if (rc != 0) return 1;

    uint8_t *pt = NULL;
    size_t pt_len = 0;
    rc = bge_decrypt_store(&store, &pt, &pt_len);
    bge_store_wipe(&store);
    if (rc != 0) return 1;

    char name_buf[256] = {0};
    if (!name) {
        int r = bge_read_visible("Name: ", name_buf, sizeof(name_buf));
        if (r <= 0) {
            fprintf(stderr, "error: name required\n");
            cfx_memzero_s(pt, pt_len);
            free(pt);
            bge_store_wipe(&store);
            return 1;
        }
        name = name_buf;
    }

    /* resolve numeric index to entry name */
    char resolved[256];
    if (is_numeric(name)) {
        unsigned idx = (unsigned)atoi(name);
        uint16_t nlen;
        const uint8_t *nptr = store_name_by_index(pt, pt_len, idx, &nlen);
        if (nptr && nlen < sizeof(resolved)) {
            memcpy(resolved, nptr, nlen);
            resolved[nlen] = '\0';
            name = resolved;
        }
    }

    size_t vlen;
    const uint8_t *val = store_get(pt, pt_len, name, &vlen);
    if (!val) {
        fprintf(stderr, "error: '%s' not found\n", name);
        cfx_memzero_s(pt, pt_len);
        free(pt);
        return 1;
    }

    if (quoted) printf("'");
    fwrite(val, 1, vlen, stdout);
    if (quoted) printf("'");
#ifndef _WIN32
    if (isatty(STDOUT_FILENO)) {
        printf("\n");
    }
#else
    if (_isatty(_fileno(stdout))) {
        printf("\n");
    }
#endif

    cfx_memzero_s(pt, pt_len);
    free(pt);
    return 0;
}

static int bge_set(int argc, char **argv) {
    const char *name = NULL;
    const char *value_str = NULL;
    const char *file_arg = NULL;
    const char *path = NULL;
    int force = 0;
    char path_buf[1024];

    for (int i = 2; i < argc; i++) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            printf("Usage: %s set [name] [value] [-s path] [-f]\n", argv[0]);
            printf("\nOptions:\n");
            printf("  -s, --store <path>   Path to BGE store\n");
            printf("  -i <file>            Read value from file\n");
            printf("                       If no value and no -i, reads from stdin.\n");
            printf("  -f, --force          Overwrite existing entry without asking\n");
            return 0;
        } else if (strcmp(argv[i], "-s") == 0 || strcmp(argv[i], "--store") == 0) {
            if (i + 1 >= argc) { fprintf(stderr, "error: -s requires a path\n"); return 1; }
            path = argv[++i];
        } else if (strcmp(argv[i], "-f") == 0 || strcmp(argv[i], "--force") == 0) {
            force = 1;
        } else if (strcmp(argv[i], "-i") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "error: -i requires a filename argument\n");
                return 1;
            }
            file_arg = argv[++i];
        } else if (!name) {
            name = argv[i];
        } else if (!value_str && !file_arg) {
            value_str = argv[i];
        } else {
            fprintf(stderr, "error: unknown argument: %s\n", argv[i]);
            return 1;
        }
    }

    if (!path) {
        if (bge_default_path(path_buf, sizeof(path_buf)) != 0) return 1;
        path = path_buf;
    }

    /* grab -i value early (before passphrase prompt) */
    uint8_t *val_buf = NULL;
    size_t val_len = 0;
    int val_needs_free = 0;

    if (file_arg) {
        FILE *vf = fopen(file_arg, "rb");
        if (!vf) {
            fprintf(stderr, "error: cannot open %s: %s\n", file_arg, strerror(errno));
            return 1;
        }
        if (bge_read_all(vf, &val_buf, &val_len) != 0) {
            fclose(vf);
            fprintf(stderr, "error: cannot read %s\n", file_arg);
            return 1;
        }
        fclose(vf);
        val_needs_free = 1;
    }

    /* auth before any interactive prompts */
    char pwd[256] = {0};
    int pwd_len = bge_read_passphrase("Enter passphrase: ", pwd, sizeof(pwd));
    if (pwd_len <= 0) {
        fprintf(stderr, "error: passphrase required\n");
        cfx_memzero_s(pwd, sizeof(pwd));
        if (val_needs_free) {
            cfx_memzero_s(val_buf, val_len);
            free(val_buf);
        }
        return 1;
    }

    bge_store store = {0};
    int rc = bge_authenticate(path, pwd, (size_t)pwd_len, &store);
    cfx_memzero_s(pwd, sizeof(pwd));
    if (rc != 0) {
        if (val_needs_free) {
            cfx_memzero_s(val_buf, val_len);
            free(val_buf);
        }
        return 1;
    }

    uint8_t *pt = NULL;
    size_t pt_len = 0;
    rc = bge_decrypt_store(&store, &pt, &pt_len);
    if (rc != 0) {
        bge_store_wipe(&store);
        if (val_needs_free) {
            cfx_memzero_s(val_buf, val_len);
            free(val_buf);
        }
        return 1;
    }

    /* interactive name prompt if not on cmdline */
    char name_buf[256] = {0};
    if (!name) {
        int r = bge_read_visible("Name: ", name_buf, sizeof(name_buf));
        if (r <= 0) {
            fprintf(stderr, "error: name required\n");
            cfx_memzero_s(pt, pt_len);
            free(pt);
            bge_store_wipe(&store);
            if (val_needs_free) {
                cfx_memzero_s(val_buf, val_len);
                free(val_buf);
            }
            return 1;
        }
        name = name_buf;
    }

    /* resolve numeric index to entry name */
    char resolved[256];
    if (is_numeric(name)) {
        unsigned idx = (unsigned)atoi(name);
        uint16_t nlen;
        const uint8_t *nptr = store_name_by_index(pt, pt_len, idx, &nlen);
        if (nptr && nlen < sizeof(resolved)) {
            memcpy(resolved, nptr, nlen);
            resolved[nlen] = '\0';
            name = resolved;
        }
    }

    /* prompt to overwrite unless -f */
    if (!force) {
        if (store_get(pt, pt_len, name, NULL) != NULL) {
            char ans[8] = {0};
            int r = bge_read_visible("Overwrite existing entry? [y/N] ", ans, sizeof(ans));
            if (r <= 0 || (ans[0] != 'y' && ans[0] != 'Y')) {
                fprintf(stderr, "Aborted.\n");
                cfx_memzero_s(pt, pt_len);
                free(pt);
                bge_store_wipe(&store);
                if (val_needs_free) {
                    cfx_memzero_s(val_buf, val_len);
                    free(val_buf);
                }
                return 1;
            }
        }
    }

    /* value: cmdline arg, interactive, or piped stdin */
    if (!val_buf) {
        if (value_str) {
            val_buf = (uint8_t *)value_str;
            val_len = strlen(value_str);
        } else {
#ifndef _WIN32
            if (isatty(STDIN_FILENO)) {
                /* tty: read with echo off */
                char secret_buf[4096] = {0};
                int slen = bge_read_secret("Value: ", secret_buf, sizeof(secret_buf));
                if (slen <= 0) {
                    fprintf(stderr, "error: empty value\n");
                    cfx_memzero_s(secret_buf, sizeof(secret_buf));
                    cfx_memzero_s(pt, pt_len);
                    free(pt);
                    bge_store_wipe(&store);
                    return 1;
                }
                val_buf = malloc((size_t)slen);
                if (!val_buf) {
                    cfx_memzero_s(secret_buf, sizeof(secret_buf));
                    cfx_memzero_s(pt, pt_len);
                    free(pt);
                    bge_store_wipe(&store);
                    return 1;
                }
                memcpy(val_buf, secret_buf, (size_t)slen);
                val_len = (size_t)slen;
                cfx_memzero_s(secret_buf, sizeof(secret_buf));
                val_needs_free = 1;
            } else
#endif
            {
                /* pipe/redirect: slurp stdin */
                if (bge_read_all(stdin, &val_buf, &val_len) != 0) {
                    fprintf(stderr, "error: cannot read from stdin\n");
                    cfx_memzero_s(pt, pt_len);
                    free(pt);
                    bge_store_wipe(&store);
                    return 1;
                }
                val_needs_free = 1;
            }
        }
    }

    size_t new_len;
    uint8_t *new_pt = store_set(pt, pt_len, name, val_buf, val_len, &new_len);
    cfx_memzero_s(pt, pt_len);
    free(pt);
    if (val_needs_free) {
        cfx_memzero_s(val_buf, val_len);
        free(val_buf);
    }

    if (!new_pt) {
        fprintf(stderr, "error: allocation failed\n");
        bge_store_wipe(&store);
        return 1;
    }

    rc = bge_write_reusing_key(path, new_pt, new_len, &store);
    cfx_memzero_s(new_pt, new_len);
    free(new_pt);
    bge_store_wipe(&store);
    if (rc == 0) printf("Ok.\n");
    return rc != 0;
}

static int bge_rm(int argc, char **argv) {
    const char *name = NULL;
    const char *path = NULL;
    int force = 0;
    char path_buf[1024];

    for (int i = 2; i < argc; i++) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            printf("Usage: %s rm <name> [-s path] [-f]\n", argv[0]);
            return 0;
        } else if (strcmp(argv[i], "-s") == 0 || strcmp(argv[i], "--store") == 0) {
            if (i + 1 >= argc) { fprintf(stderr, "error: -s requires a path\n"); return 1; }
            path = argv[++i];
        } else if (strcmp(argv[i], "-f") == 0 || strcmp(argv[i], "--force") == 0) {
            force = 1;
        } else if (!name) {
            name = argv[i];
        } else {
            fprintf(stderr, "error: unknown argument: %s\n", argv[i]);
            return 1;
        }
    }

    if (!path) {
        if (bge_default_path(path_buf, sizeof(path_buf)) != 0) return 1;
        path = path_buf;
    }

    char pwd[256] = {0};
    int pwd_len = bge_read_passphrase("Enter passphrase: ", pwd, sizeof(pwd));
    if (pwd_len <= 0) {
        fprintf(stderr, "error: passphrase required\n");
        cfx_memzero_s(pwd, sizeof(pwd));
        return 1;
    }

    bge_store store = {0};
    int rc = bge_authenticate(path, pwd, (size_t)pwd_len, &store);
    cfx_memzero_s(pwd, sizeof(pwd));
    if (rc != 0) return 1;

    uint8_t *pt = NULL;
    size_t pt_len = 0;
    rc = bge_decrypt_store(&store, &pt, &pt_len);
    if (rc != 0) {
        bge_store_wipe(&store);
        return 1;
    }

    char name_buf[256] = {0};
    if (!name) {
        int r = bge_read_visible("Name: ", name_buf, sizeof(name_buf));
        if (r <= 0) {
            fprintf(stderr, "error: name required\n");
            cfx_memzero_s(pt, pt_len);
            free(pt);
            bge_store_wipe(&store);
            return 1;
        }
        name = name_buf;
    }


    if (!force && store_get(pt, pt_len, name, NULL) != NULL) {
        char prompt[300];
        snprintf(prompt, sizeof(prompt), "Remove '%s' - are you sure? [y/N] ", name);
        char ans[8] = {0};
        int r = bge_read_visible(prompt, ans, sizeof(ans));
        if (r <= 0 || (ans[0] != 'y' && ans[0] != 'Y')) {
            fprintf(stderr, "Aborted.\n");
            cfx_memzero_s(pt, pt_len);
            free(pt);
            bge_store_wipe(&store);
            return 1;
        }
    }

    size_t new_len;
    uint8_t *new_pt = store_rm(pt, pt_len, name, &new_len);
    cfx_memzero_s(pt, pt_len);
    free(pt);

    if (!new_pt) {
        bge_store_wipe(&store);
        return 1;
    }

    rc = bge_write_reusing_key(path, new_pt, new_len, &store);
    cfx_memzero_s(new_pt, new_len);
    free(new_pt);
    bge_store_wipe(&store);
    if (rc == 0) printf("Ok.\n");
    return rc != 0;
}

static int bge_ls(int argc, char **argv) {
    const char *path = NULL;
    char path_buf[1024];

    for (int i = 2; i < argc; i++) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            printf("Usage: %s ls [-s path]\n", argv[0]);
            return 0;
        } else if (strcmp(argv[i], "-s") == 0 || strcmp(argv[i], "--store") == 0) {
            if (i + 1 >= argc) { fprintf(stderr, "error: -s requires a path\n"); return 1; }
            path = argv[++i];
        } else {
            fprintf(stderr, "error: unknown argument: %s\n", argv[i]);
            return 1;
        }
    }

    if (!path) {
        if (bge_default_path(path_buf, sizeof(path_buf)) != 0) return 1;
        path = path_buf;
    }

    char pwd[256] = {0};
    int pwd_len = bge_read_passphrase("Enter passphrase: ", pwd, sizeof(pwd));
    if (pwd_len <= 0) {
        fprintf(stderr, "error: passphrase required\n");
        cfx_memzero_s(pwd, sizeof(pwd));
        return 1;
    }

    bge_store store = {0};
    int rc = bge_authenticate(path, pwd, (size_t)pwd_len, &store);
    cfx_memzero_s(pwd, sizeof(pwd));
    if (rc != 0) return 1;

    uint8_t *pt = NULL;
    size_t pt_len = 0;
    rc = bge_decrypt_store(&store, &pt, &pt_len);
    bge_store_wipe(&store);
    if (rc != 0) return 1;

    /* print entry names */
    const uint8_t *p = pt;
    const uint8_t *end = pt + pt_len;
    unsigned idx = 0;
    while (p + 2 <= end) {
        uint16_t klen = cfx_load16_le(p);
        p += 2;
        if (p + klen > end) break;
        printf("[%u] ", ++idx);
        fwrite(p, 1, klen, stdout);
        printf("\n");
        p += klen;

        if (p + 4 > end) break;
        uint32_t vl = cfx_load32_le(p);
        p += 4;
        if (p + vl > end) break;
        p += vl;
    }

    cfx_memzero_s(pt, pt_len);
    free(pt);
    return 0;
}

static int bge_info(int argc, char **argv) {
    const char *path = NULL;
    char path_buf[1024];

    for (int i = 2; i < argc; i++) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            printf("Usage: %s info [-s path]\n", argv[0]);
            return 0;
        } else if (strcmp(argv[i], "-s") == 0 || strcmp(argv[i], "--store") == 0) {
            if (i + 1 >= argc) { fprintf(stderr, "error: -s requires a path\n"); return 1; }
            path = argv[++i];
        } else {
            fprintf(stderr, "error: unknown argument: %s\n", argv[i]);
            return 1;
        }
    }

    if (!path) {
        if (bge_default_path(path_buf, sizeof(path_buf)) != 0) return 1;
        path = path_buf;
    }

    /* stat doesn't need passphrase */
    struct stat st;
    if (stat(path, &st) != 0) {
        fprintf(stderr, "error: cannot stat %s: %s\n", path, strerror(errno));
        return 1;
    }

    printf("Store:   %s\n", path);
    printf("Size:    %lld bytes\n", (long long)st.st_size);

    /* need to decrypt to count entries */
    char pwd[256] = {0};
    int pwd_len = bge_read_passphrase("Enter passphrase: ", pwd, sizeof(pwd));
    if (pwd_len <= 0) {
        fprintf(stderr, "error: passphrase required\n");
        cfx_memzero_s(pwd, sizeof(pwd));
        return 1;
    }

    bge_store store = {0};
    int rc = bge_authenticate(path, pwd, (size_t)pwd_len, &store);
    cfx_memzero_s(pwd, sizeof(pwd));
    if (rc != 0) return 1;

    uint8_t *pt = NULL;
    size_t pt_len = 0;
    rc = bge_decrypt_store(&store, &pt, &pt_len);
    if (rc != 0) {
        bge_store_wipe(&store);
        return 1;
    }

    /* count entries */
    unsigned count = 0;
    const uint8_t *p = pt;
    const uint8_t *end = pt + pt_len;
    while (p + 2 <= end) {
        uint16_t klen = cfx_load16_le(p);
        p += 2;
        if (p + klen > end) break;
        p += klen;
        if (p + 4 > end) break;
        uint32_t vl = cfx_load32_le(p);
        p += 4;
        if (p + vl > end) break;
        p += vl;
        count++;
    }

    cfx_memzero_s(pt, pt_len);
    free(pt);

    printf("Entries: %u\n", count);

    /* argon2 params */
    uint32_t m_cost = cfx_load32_le(&store.hdr.m_cost);
    uint32_t t_cost = cfx_load32_le(&store.hdr.t_cost);
    uint32_t p_cost = cfx_load32_le(&store.hdr.p_cost);
    printf("Argon2:  m=%u KB, t=%u, p=%u\n", m_cost, t_cost, p_cost);

    bge_store_wipe(&store);
    return 0;
}

static int bge_passwd(int argc, char **argv) {
    const char *path = NULL;
    char path_buf[1024];

    for (int i = 2; i < argc; i++) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            printf("Usage: %s passwd [-s path]\n", argv[0]);
            return 0;
        } else if (strcmp(argv[i], "-s") == 0 || strcmp(argv[i], "--store") == 0) {
            if (i + 1 >= argc) { fprintf(stderr, "error: -s requires a path\n"); return 1; }
            path = argv[++i];
        } else {
            fprintf(stderr, "error: unknown argument: %s\n", argv[i]);
            return 1;
        }
    }

    if (!path) {
        if (bge_default_path(path_buf, sizeof(path_buf)) != 0) return 1;
        path = path_buf;
    }

    /* old passphrase */
    char old_pwd[256] = {0};
    int old_len = bge_read_passphrase("Enter current passphrase: ", old_pwd, sizeof(old_pwd));
    if (old_len <= 0) {
        fprintf(stderr, "error: passphrase required\n");
        cfx_memzero_s(old_pwd, sizeof(old_pwd));
        return 1;
    }

    bge_store store = {0};
    int rc = bge_authenticate(path, old_pwd, (size_t)old_len, &store);
    cfx_memzero_s(old_pwd, sizeof(old_pwd));
    if (rc != 0) return 1;

    uint8_t *pt = NULL;
    size_t pt_len = 0;
    rc = bge_decrypt_store(&store, &pt, &pt_len);
    bge_store_wipe(&store);
    if (rc != 0) return 1;

    /* new passphrase (double prompt) */
    printf("Enter new passphrase.\n");
    char new_pwd[256] = {0};
    int new_len = prompt_passphrase(new_pwd, sizeof(new_pwd));
    if (new_len < 0) {
        cfx_memzero_s(pt, pt_len);
        free(pt);
        return 1;
    }

    rc = bge_encrypt_write(path, pt, pt_len, new_pwd, (size_t)new_len,
                           BGE_DEFAULT_M, BGE_DEFAULT_T, BGE_DEFAULT_P);
    cfx_memzero_s(pt, pt_len);
    free(pt);
    cfx_memzero_s(new_pwd, sizeof(new_pwd));

    if (rc == 0)
        printf("Passphrase changed.\n");
    return rc != 0;
}

static int bge_dump(int argc, char **argv) {
    const char *path = NULL;
    char path_buf[1024];

    for (int i = 2; i < argc; i++) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            printf("Usage: %s dump [-s path]\n", argv[0]);
            return 0;
        } else if (strcmp(argv[i], "-s") == 0 || strcmp(argv[i], "--store") == 0) {
            if (i + 1 >= argc) { fprintf(stderr, "error: -s requires a path\n"); return 1; }
            path = argv[++i];
        } else {
            fprintf(stderr, "error: unknown argument: %s\n", argv[i]);
            return 1;
        }
    }

    if (!path) {
        if (bge_default_path(path_buf, sizeof(path_buf)) != 0) return 1;
        path = path_buf;
    }

    char pwd[256] = {0};
    int pwd_len = bge_read_passphrase("Enter passphrase: ", pwd, sizeof(pwd));
    if (pwd_len <= 0) {
        fprintf(stderr, "error: passphrase required\n");
        cfx_memzero_s(pwd, sizeof(pwd));
        return 1;
    }

    bge_store store = {0};
    int rc = bge_authenticate(path, pwd, (size_t)pwd_len, &store);
    cfx_memzero_s(pwd, sizeof(pwd));
    if (rc != 0) return 1;

    uint8_t *pt = NULL;
    size_t pt_len = 0;
    rc = bge_decrypt_store(&store, &pt, &pt_len);
    bge_store_wipe(&store);
    if (rc != 0) return 1;

    /* render + print */
    size_t text_len = 0;
    uint8_t *text = store_to_text(pt, pt_len, &text_len);
    cfx_memzero_s(pt, pt_len);
    free(pt);

    if (!text) {
        fprintf(stderr, "error: allocation failed\n");
        return 1;
    }

    if (text_len > 0)
        fwrite(text, 1, text_len, stdout);

    cfx_memzero_s(text, text_len);
    free(text);
    return 0;
}

static int bge_encrypt_file(int argc, char **argv) {
    const char *output = NULL;
    const char *input = NULL;
    int armor = 0;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-e") == 0 || strcmp(argv[i], "--encrypt") == 0) {
            continue;
        } else if (strcmp(argv[i], "-a") == 0 || strcmp(argv[i], "--armor") == 0) {
            armor = 1;
        } else if (strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "--output") == 0) {
            if (i + 1 >= argc) { fprintf(stderr, "error: -o requires a path\n"); return 1; }
            output = argv[++i];
        } else if (strcmp(argv[i], "-i") == 0 || strcmp(argv[i], "--input") == 0) {
            if (i + 1 >= argc) { fprintf(stderr, "error: -i requires a path\n"); return 1; }
            input = argv[++i];
        } else if (argv[i][0] == '-') {
            fprintf(stderr, "error: unknown option: %s\n", argv[i]);
            return 1;
        } else if (!input) {
            input = argv[i];
        } else {
            fprintf(stderr, "error: unexpected argument: %s\n", argv[i]);
            return 1;
        }
    }

    FILE *inf = stdin;
    if (input) {
        inf = fopen(input, "rb");
        if (!inf) {
            fprintf(stderr, "error: cannot open %s: %s\n", input, strerror(errno));
            return 1;
        }
    }

    char pwd[256] = {0};
    int pwd_len = prompt_passphrase(pwd, sizeof(pwd));
    if (pwd_len < 0) {
        if (input) fclose(inf);
        return 1;
    }

    /* build v3 header */
    bge_header header;
    uint8_t kdf_out[48], verifier[BGE_VERIFIER_LEN];

    memcpy(header.magic, BGE_MAGIC, 3);
    header.version = BGE_STREAM_VERSION;
    cfx_store32_le(&header.m_cost, BGE_DEFAULT_M);
    cfx_store32_le(&header.t_cost, BGE_DEFAULT_T);
    cfx_store32_le(&header.p_cost, BGE_DEFAULT_P);

    cfx_srand_os();
    cfx_rand_bytes(header.salt,  sizeof(header.salt));
    cfx_rand_bytes(header.nonce, sizeof(header.nonce));

    int rc = cfx_argon2id(kdf_out, 48,
        (const uint8_t *)pwd, (size_t)pwd_len,
        header.salt, sizeof(header.salt),
        BGE_DEFAULT_M, BGE_DEFAULT_T, BGE_DEFAULT_P);
    cfx_memzero_s(pwd, sizeof(pwd));
    if (rc != 0) {
        fprintf(stderr, "error: argon2 key derivation failed\n");
        cfx_memzero_s(kdf_out, sizeof(kdf_out));
        if (input) fclose(inf);
        return 1;
    }
    memcpy(verifier, kdf_out + 32, BGE_VERIFIER_LEN);

    /* open output */
    FILE *outf = stdout;
    if (output) {
        outf = fopen(output, "wb");
        if (!outf) {
            fprintf(stderr, "error: cannot open %s: %s\n", output, strerror(errno));
            cfx_memzero_s(kdf_out, sizeof(kdf_out));
            if (input) fclose(inf);
            return 1;
        }
    }

    /* for armor: accumulate all output in memory, then encode at the end */
    uint8_t *armor_buf = NULL;
    size_t armor_cap = 0, armor_len = 0;
    int ret = 0;

    #define EMIT(data, len) do { \
        if (armor) { \
            size_t _n = (len); \
            if (armor_len + _n > armor_cap) { \
                size_t _nc = armor_cap ? armor_cap * 2 : 4096; \
                while (_nc < armor_len + _n) _nc *= 2; \
                uint8_t *_t = realloc(armor_buf, _nc); \
                if (!_t) { ret = 1; goto done; } \
                armor_buf = _t; armor_cap = _nc; \
            } \
            memcpy(armor_buf + armor_len, (data), _n); \
            armor_len += _n; \
        } else { \
            if (fwrite((data), 1, (len), outf) != (len)) { ret = 1; goto done; } \
        } \
    } while(0)

    /* write header + verifier */
    EMIT(&header, sizeof(header));
    EMIT(verifier, BGE_VERIFIER_LEN);

    /* streaming encrypt loop */
    uint8_t pt_buf[CFX_STREAM_CHUNK_SIZE];
    uint8_t ct_buf[CFX_STREAM_CHUNK_SIZE];
    uint8_t tag[CFX_STREAM_TAG_SIZE];
    uint64_t chunk_counter = 0;

    for (;;) {
        size_t nread = fread(pt_buf, 1, CFX_STREAM_CHUNK_SIZE, inf);
        if (nread == 0 && ferror(inf)) {
            fprintf(stderr, "error: read failed\n");
            ret = 1; goto done;
        }

        /* peek ahead to determine if this is the final chunk */
        int is_final = feof(inf);
        if (!is_final && nread < CFX_STREAM_CHUNK_SIZE)
            is_final = 1;

        rc = cfx_stream_xchacha20_poly1305_encrypt_chunk(
            ct_buf, tag, pt_buf, nread,
            chunk_counter, is_final, kdf_out, header.nonce);
        cfx_memzero_s(pt_buf, sizeof(pt_buf));
        if (rc != 0) {
            fprintf(stderr, "error: encryption failed at chunk %llu\n",
                    (unsigned long long)chunk_counter);
            ret = 1; goto done;
        }

        EMIT(ct_buf, nread);
        EMIT(tag, CFX_STREAM_TAG_SIZE);
        chunk_counter++;

        if (is_final) break;
    }

    #undef EMIT

    /* finalize armor if needed */
    if (armor && ret == 0) {
        uint8_t *armored = NULL;
        size_t armored_len = 0;
        if (bge_armor_encode(armor_buf, armor_len, &armored, &armored_len) != 0) {
            fprintf(stderr, "error: armor encoding failed\n");
            ret = 1;
        } else {
            if (fwrite(armored, 1, armored_len, outf) != armored_len) ret = 1;
            free(armored);
        }
    }

done:
    cfx_memzero_s(kdf_out, sizeof(kdf_out));
    cfx_memzero_s(verifier, sizeof(verifier));
    cfx_memzero_s(pt_buf, sizeof(pt_buf));
    cfx_memzero_s(ct_buf, sizeof(ct_buf));
    if (armor_buf) { free(armor_buf); }
    if (input) fclose(inf);
    if (output && outf) fclose(outf);
    return ret;
}

/* v2 single-shot decrypt from in-memory buffer */
static int bge_decrypt_v2(const uint8_t *file_buf, size_t file_len,
                          const uint8_t key[32], const uint8_t nonce[24],
                          FILE *outf) {
    size_t ct_len = file_len - BGE_AAD_LEN - BGE_TAG_LEN;
    const uint8_t *ct  = file_buf + BGE_AAD_LEN;
    const uint8_t *tag = file_buf + file_len - BGE_TAG_LEN;

    uint8_t *pt = malloc(ct_len + 1);
    if (!pt) {
        fprintf(stderr, "error: allocation failed\n");
        return -1;
    }

    int rc = cfx_xchacha20_poly1305_decrypt(
        pt, ct, ct_len, file_buf, BGE_AAD_LEN, key, nonce, tag);
    if (rc != 0) {
        fprintf(stderr, "error: decryption failed (corrupted or tampered)\n");
        cfx_memzero_s(pt, ct_len + 1); free(pt);
        return -1;
    }

    int ret = 0;
    if (ct_len > 0 && fwrite(pt, 1, ct_len, outf) != ct_len)
        ret = -1;

    cfx_memzero_s(pt, ct_len + 1); free(pt);
    return ret;
}

/* v3 decrypt from in-memory buffer (used for armored v3 input) */
static int bge_decrypt_v3(const uint8_t *file_buf, size_t file_len,
                          const uint8_t key[32], const uint8_t nonce[24],
                          FILE *outf) {
    const uint8_t *p = file_buf + BGE_AAD_LEN;
    const uint8_t *end = file_buf + file_len;
    uint8_t pt_buf[CFX_STREAM_CHUNK_SIZE];
    uint64_t chunk_counter = 0;

    while (p < end) {
        size_t remaining = (size_t)(end - p);
        if (remaining < CFX_STREAM_TAG_SIZE) {
            fprintf(stderr, "error: truncated stream (no tag)\n");
            return -1;
        }

        /* determine chunk data size: remaining minus tag, capped at chunk size */
        size_t chunk_plus_tag;
        int is_final;

        if (remaining <= CFX_STREAM_CHUNK_SIZE + CFX_STREAM_TAG_SIZE) {
            /* this is the last chunk */
            chunk_plus_tag = remaining;
            is_final = 1;
        } else {
            /* full chunk */
            chunk_plus_tag = CFX_STREAM_CHUNK_SIZE + CFX_STREAM_TAG_SIZE;
            is_final = 0;
        }

        size_t ct_len = chunk_plus_tag - CFX_STREAM_TAG_SIZE;
        const uint8_t *ct  = p;
        const uint8_t *tag = p + ct_len;

        int rc = cfx_stream_xchacha20_poly1305_decrypt_chunk(
            pt_buf, ct, ct_len, tag, chunk_counter, is_final, key, nonce);
        if (rc != 0) {
            fprintf(stderr, "error: chunk %llu authentication failed\n",
                    (unsigned long long)chunk_counter);
            cfx_memzero_s(pt_buf, sizeof(pt_buf));
            return -1;
        }

        if (ct_len > 0 && fwrite(pt_buf, 1, ct_len, outf) != ct_len) {
            cfx_memzero_s(pt_buf, sizeof(pt_buf));
            return -1;
        }

        cfx_memzero_s(pt_buf, ct_len);
        p += chunk_plus_tag;
        chunk_counter++;
    }

    return 0;
}

/* v3 streaming decrypt directly from FILE* (no buffering) — poly1305 */
static int bge_decrypt_v3_stream(FILE *inf,
                                  const uint8_t key[32], const uint8_t nonce[24],
                                  FILE *outf) {
    uint8_t chunk_buf[CFX_STREAM_CHUNK_SIZE + CFX_STREAM_TAG_SIZE];
    uint8_t pt_buf[CFX_STREAM_CHUNK_SIZE];
    uint64_t chunk_counter = 0;
    uint8_t lookahead;
    int have_lookahead = 0;

    for (;;) {
        size_t off = 0;
        if (have_lookahead) {
            chunk_buf[0] = lookahead;
            off = 1;
            have_lookahead = 0;
        }

        size_t nread = fread(chunk_buf + off, 1, sizeof(chunk_buf) - off, inf);
        size_t total = off + nread;

        if (nread == 0 && off == 0 && ferror(inf)) {
            fprintf(stderr, "error: read failed\n");
            cfx_memzero_s(chunk_buf, sizeof(chunk_buf));
            return -1;
        }

        if (total == 0) {
            fprintf(stderr, "error: empty stream (missing final chunk)\n");
            return -1;
        }

        if (total < CFX_STREAM_TAG_SIZE) {
            fprintf(stderr, "error: truncated stream (no tag)\n");
            cfx_memzero_s(chunk_buf, sizeof(chunk_buf));
            return -1;
        }

        int is_final;
        if (total < sizeof(chunk_buf)) {
            /* got less than a full chunk+tag: must be final */
            is_final = 1;
        } else {
            /* full chunk+tag: peek one byte to check for more data */
            size_t peeked = fread(&lookahead, 1, 1, inf);
            if (peeked == 0) {
                is_final = 1;   /* exact-size final chunk */
            } else {
                is_final = 0;
                have_lookahead = 1;
            }
        }

        size_t ct_len = total - CFX_STREAM_TAG_SIZE;
        const uint8_t *tag = chunk_buf + ct_len;

        int rc = cfx_stream_xchacha20_poly1305_decrypt_chunk(
            pt_buf, chunk_buf, ct_len, tag, chunk_counter, is_final, key, nonce);
        if (rc != 0) {
            fprintf(stderr, "error: chunk %llu authentication failed\n",
                    (unsigned long long)chunk_counter);
            cfx_memzero_s(pt_buf, sizeof(pt_buf));
            cfx_memzero_s(chunk_buf, sizeof(chunk_buf));
            return -1;
        }

        if (ct_len > 0 && fwrite(pt_buf, 1, ct_len, outf) != ct_len) {
            cfx_memzero_s(pt_buf, sizeof(pt_buf));
            cfx_memzero_s(chunk_buf, sizeof(chunk_buf));
            return -1;
        }

        cfx_memzero_s(pt_buf, ct_len);
        chunk_counter++;

        if (is_final) break;
    }

    cfx_memzero_s(chunk_buf, sizeof(chunk_buf));
    return 0;
}

static int bge_decrypt_file(int argc, char **argv) {
    const char *output = NULL;
    const char *input = NULL;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-d") == 0 || strcmp(argv[i], "--decrypt") == 0)
            continue;
        else if (strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "--output") == 0) {
            if (i + 1 >= argc) { fprintf(stderr, "error: -o requires a path\n"); return 1; }
            output = argv[++i];
        } else if (strcmp(argv[i], "-i") == 0 || strcmp(argv[i], "--input") == 0) {
            if (i + 1 >= argc) { fprintf(stderr, "error: -i requires a path\n"); return 1; }
            input = argv[++i];
        } else if (argv[i][0] == '-') {
            fprintf(stderr, "error: unknown option: %s\n", argv[i]);
            return 1;
        } else if (!input) {
            input = argv[i];
        } else {
            fprintf(stderr, "error: unexpected argument: %s\n", argv[i]);
            return 1;
        }
    }

    FILE *inf = stdin;
    if (input) {
        inf = fopen(input, "rb");
        if (!inf) {
            fprintf(stderr, "error: cannot open %s: %s\n", input, strerror(errno));
            return 1;
        }
    }

    /* peek at the header to decide between streaming and buffered paths */
    uint8_t hdr_peek[BGE_AAD_LEN];
    size_t hdr_n = fread(hdr_peek, 1, BGE_AAD_LEN, inf);

    int is_binary_v3 = (hdr_n == BGE_AAD_LEN &&
                        memcmp(hdr_peek, BGE_MAGIC, 3) == 0 &&
                        hdr_peek[3] == BGE_STREAM_VERSION);

    if (is_binary_v3) {
        /* ── non-armored v3: true streaming decrypt from FILE* ── */
        bge_header hdr;
        memcpy(&hdr, hdr_peek, sizeof(hdr));

        uint32_t m_cost = cfx_load32_le(&hdr.m_cost);
        uint32_t t_cost = cfx_load32_le(&hdr.t_cost);
        uint32_t p      = cfx_load32_le(&hdr.p_cost);

        if (m_cost < 8 || t_cost < 1 || p < 1 ||
            m_cost > 4194304 || t_cost > 100 || p > 16) {
            fprintf(stderr, "error: unreasonable argon2 parameters\n");
            if (input) fclose(inf);
            return 1;
        }

        char pwd[256] = {0};
        int pwd_len = bge_read_passphrase("Enter passphrase: ", pwd, sizeof(pwd));
        if (pwd_len <= 0) {
            fprintf(stderr, "error: passphrase required\n");
            cfx_memzero_s(pwd, sizeof(pwd));
            if (input) fclose(inf);
            return 1;
        }

        uint8_t kdf_out[48];
        int rc = cfx_argon2id(kdf_out, 48,
            (const uint8_t *)pwd, (size_t)pwd_len,
            hdr.salt, sizeof(hdr.salt), m_cost, t_cost, p);
        cfx_memzero_s(pwd, sizeof(pwd));
        if (rc != 0) {
            fprintf(stderr, "error: argon2 key derivation failed\n");
            cfx_memzero_s(kdf_out, sizeof(kdf_out));
            if (input) fclose(inf);
            return 1;
        }

        const uint8_t *stored_verifier = hdr_peek + BGE_HEADER_LEN;
        uint8_t diff = 0;
        for (int i = 0; i < BGE_VERIFIER_LEN; i++)
            diff |= kdf_out[32 + i] ^ stored_verifier[i];

        if (diff != 0) {
            fprintf(stderr, "error: wrong passphrase\n");
            cfx_memzero_s(kdf_out, sizeof(kdf_out));
            if (input) fclose(inf);
            return 1;
        }

        FILE *outf = stdout;
        if (output) {
            outf = fopen(output, "wb");
            if (!outf) {
                fprintf(stderr, "error: cannot open %s: %s\n", output, strerror(errno));
                cfx_memzero_s(kdf_out, sizeof(kdf_out));
                if (input) fclose(inf);
                return 1;
            }
        }

        int ret = bge_decrypt_v3_stream(inf, kdf_out, hdr.nonce, outf);
        cfx_memzero_s(kdf_out, sizeof(kdf_out));
        cfx_memzero_s(hdr_peek, sizeof(hdr_peek));
        if (input) fclose(inf);
        if (output) fclose(outf);
        return ret != 0 ? 1 : 0;
    }

    /* ── armored or v2: read everything into memory ── */
    uint8_t *rest = NULL;
    size_t rest_len = 0;
    int rc = bge_read_all(inf, &rest, &rest_len);
    if (input) fclose(inf);
    if (rc != 0) {
        fprintf(stderr, "error: cannot read input\n");
        return 1;
    }

    /* prepend the bytes we already peeked */
    size_t raw_len = hdr_n + rest_len;
    uint8_t *raw = malloc(raw_len);
    if (!raw) {
        fprintf(stderr, "error: allocation failed\n");
        if (rest) { cfx_memzero_s(rest, rest_len); free(rest); }
        return 1;
    }
    memcpy(raw, hdr_peek, hdr_n);
    if (rest_len > 0) memcpy(raw + hdr_n, rest, rest_len);
    if (rest) { cfx_memzero_s(rest, rest_len); free(rest); }
    cfx_memzero_s(hdr_peek, sizeof(hdr_peek));

    /* auto-detect and strip armor */
    uint8_t *file_buf = raw;
    size_t file_len = raw_len;
    int decoded_armor = 0;

    if (bge_is_armored(raw, raw_len)) {
        uint8_t *dec = NULL;
        size_t dec_len = 0;
        if (bge_armor_decode(raw, raw_len, &dec, &dec_len) != 0) {
            fprintf(stderr, "error: invalid armored input\n");
            cfx_memzero_s(raw, raw_len); free(raw);
            return 1;
        }
        file_buf = dec;
        file_len = dec_len;
        decoded_armor = 1;
    }

    /* validate header */
    if (file_len < BGE_AAD_LEN + BGE_TAG_LEN ||
        memcmp(file_buf, BGE_MAGIC, 3) != 0) {
        fprintf(stderr, "error: not a valid BGE file\n");
        goto fail;
    }

    uint8_t version = file_buf[3];
    if (version != BGE_VERSION && version != BGE_STREAM_VERSION) {
        fprintf(stderr, "error: unsupported BGE version %u\n", version);
        goto fail;
    }

    char pwd[256] = {0};
    int pwd_len = bge_read_passphrase("Enter passphrase: ", pwd, sizeof(pwd));
    if (pwd_len <= 0) {
        fprintf(stderr, "error: passphrase required\n");
        cfx_memzero_s(pwd, sizeof(pwd));
        goto fail;
    }

    /* authenticate via argon2 + verifier */
    bge_header hdr;
    memcpy(&hdr, file_buf, sizeof(hdr));

    uint32_t m_cost = cfx_load32_le(&hdr.m_cost);
    uint32_t t_cost = cfx_load32_le(&hdr.t_cost);
    uint32_t p      = cfx_load32_le(&hdr.p_cost);

    if (m_cost < 8 || t_cost < 1 || p < 1 ||
        m_cost > 4194304 || t_cost > 100 || p > 16) {
        fprintf(stderr, "error: unreasonable argon2 parameters\n");
        cfx_memzero_s(pwd, sizeof(pwd));
        goto fail;
    }

    uint8_t kdf_out[48];
    rc = cfx_argon2id(kdf_out, 48,
        (const uint8_t *)pwd, (size_t)pwd_len,
        hdr.salt, sizeof(hdr.salt), m_cost, t_cost, p);
    cfx_memzero_s(pwd, sizeof(pwd));
    if (rc != 0) {
        fprintf(stderr, "error: argon2 key derivation failed\n");
        cfx_memzero_s(kdf_out, sizeof(kdf_out));
        goto fail;
    }

    const uint8_t *stored_verifier = file_buf + BGE_HEADER_LEN;
    uint8_t diff = 0;
    for (int i = 0; i < BGE_VERIFIER_LEN; i++)
        diff |= kdf_out[32 + i] ^ stored_verifier[i];

    if (diff != 0) {
        fprintf(stderr, "error: wrong passphrase\n");
        cfx_memzero_s(kdf_out, sizeof(kdf_out));
        goto fail;
    }

    /* open output */
    FILE *outf = stdout;
    if (output) {
        outf = fopen(output, "wb");
        if (!outf) {
            fprintf(stderr, "error: cannot open %s: %s\n", output, strerror(errno));
            cfx_memzero_s(kdf_out, sizeof(kdf_out));
            goto fail;
        }
    }

    /* dispatch by version */
    int ret;
    if (version == BGE_VERSION) {
        ret = bge_decrypt_v2(file_buf, file_len, kdf_out, hdr.nonce, outf);
    } else {
        ret = bge_decrypt_v3(file_buf, file_len, kdf_out, hdr.nonce, outf);
    }

    cfx_memzero_s(kdf_out, sizeof(kdf_out));
    if (output) fclose(outf);
    cfx_memzero_s(file_buf, file_len);
    if (decoded_armor) free(file_buf);
    cfx_memzero_s(raw, raw_len); free(raw);
    return ret != 0 ? 1 : 0;

fail:
    cfx_memzero_s(file_buf, file_len);
    if (decoded_armor) free(file_buf);
    cfx_memzero_s(raw, raw_len); free(raw);
    return 1;
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

    /* scan for -e/-d anywhere in argv (flag-style, not subcommand) */
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-e") == 0 || strcmp(argv[i], "--encrypt") == 0) {
            return bge_encrypt_file(argc, argv);
        }
        if (strcmp(argv[i], "-d") == 0 || strcmp(argv[i], "--decrypt") == 0) {
            return bge_decrypt_file(argc, argv);
        }
    }

    if (strcmp(cmd, "init")   == 0) return bge_init(argc, argv);
    if (strcmp(cmd, "get")    == 0) return bge_get(argc, argv);
    if (strcmp(cmd, "set")    == 0) return bge_set(argc, argv);
    if (strcmp(cmd, "rm")     == 0) return bge_rm(argc, argv);
    if (strcmp(cmd, "ls")     == 0 ||
        strcmp(cmd, "list")   == 0) return bge_ls(argc, argv);
    if (strcmp(cmd, "info")   == 0) return bge_info(argc, argv);
    if (strcmp(cmd, "passwd") == 0) return bge_passwd(argc, argv);
    if (strcmp(cmd, "dump")   == 0) return bge_dump(argc, argv);

    fprintf(stderr, "Unknown command: %s\n", cmd);
    usage(argv[0]);
    return 1;
}
