/* cfx_bge.c -- passphrase-protected secrets store */

#include "cfx/argon2.h"
#include "cfx/aead_chacha20_poly1305.h"
#include "cfx/rand.h"
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
#define BGE_VERSION       2       /* file format version (byte [3]) */
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


static void usage(const char *prog) {
    printf("Usage: %s <command> [args...]\n\n", prog);
    printf("Passphrase-protected secrets store. Each store is a single encrypted\n");
    printf("file using Argon2id key derivation and XChaCha20-Poly1305 AEAD.\n");
    printf("Values are arbitrary bytes (passwords, keys, certificates).\n");
    printf("Passphrases and secret values are entered with echo disabled.\n\n");
    printf("Commands:\n");
    printf("  init   [-s path]                 Create a new encrypted store\n");
    printf("  get    <name> [-s path] [-q]     Retrieve a secret (raw to stdout)\n");
    printf("  set    [name] [value] [-s path]  Set a secret (value from arg, -f, or stdin)\n");
    printf("  rm     <name> [-s path]          Remove a secret\n");
    printf("  ls     [-s path]                 List all secret names\n");
    printf("  info   [-s path]                 Show store location, size, and entry count\n");
    printf("  passwd [-s path]                 Change passphrase\n");
    printf("  dump   [-s path]                 Print all secrets to stdout\n");
    printf("\nOptions:\n");
    printf("  -s, --store <path>   Path to BGE store (default: ~/.cfx/secrets.bge)\n");
    printf("\nOptions for get:\n");
    printf("  -q, --quoted         Wrap output in single quotes\n");
    printf("\nOptions for set:\n");
    printf("  -f <file>   Read value from file (binary-safe)\n");
    printf("              If no value and no -f, reads from stdin.\n");
    printf("  -y, --yes   Overwrite existing entry without asking\n");
    printf("\nExamples:\n");
    printf("  %s init                              Create store at default path\n", prog);
    printf("  %s init -s /tmp/secrets.bge          Create store at custom path\n", prog);
    printf("  %s set token mytoken                 Set from command line argument\n", prog);
    printf("  %s set cert -f server.pem            Set from file (binary-safe)\n", prog);
    printf("  echo secret | %s set apikey          Set from piped stdin\n", prog);
#ifdef __APPLE__
    printf("  pbpaste | %s set apikey              Set from clipboard\n", prog);
#endif
    printf("  %s get token                         Print value (newline on tty)\n", prog);
#ifdef __APPLE__
    printf("  %s get token | pbcopy                Copy to clipboard\n", prog);
#endif
    printf("  %s get token > out.txt               Write raw value to file\n", prog);
    printf("  %s ls                                Show all key names\n", prog);
    printf("  %s info                              Show store stats\n", prog);
    printf("  %s dump                              Show all entries as [name]/value\n", prog);
    printf("  %s passwd                            Change the passphrase\n", prog);
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
            if (buf) { cfx_memzero_s(buf, cap); free(buf); }
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
    int pwd_len = bge_read_secret("Enter passphrase: ", pwd, sizeof(pwd));
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
            printf("Usage: %s set [name] [value] [-s path] [-y]\n", argv[0]);
            printf("\nOptions:\n");
            printf("  -s, --store <path>   Path to BGE store\n");
            printf("  -f <file>            Read value from file\n");
            printf("                       If no value and no -f, reads from stdin.\n");
            printf("  -y, --yes            Overwrite existing entry without asking\n");
            return 0;
        } else if (strcmp(argv[i], "-s") == 0 || strcmp(argv[i], "--store") == 0) {
            if (i + 1 >= argc) { fprintf(stderr, "error: -s requires a path\n"); return 1; }
            path = argv[++i];
        } else if (strcmp(argv[i], "-y") == 0 || strcmp(argv[i], "--yes") == 0) {
            force = 1;
        } else if (strcmp(argv[i], "-f") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "error: -f requires a filename argument\n");
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

    /* grab -f value early (before passphrase prompt) */
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
    int pwd_len = bge_read_secret("Enter passphrase: ", pwd, sizeof(pwd));
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

    /* prompt to overwrite unless -y */
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
            printf("Usage: %s rm <name> [-s path] [-y]\n", argv[0]);
            return 0;
        } else if (strcmp(argv[i], "-s") == 0 || strcmp(argv[i], "--store") == 0) {
            if (i + 1 >= argc) { fprintf(stderr, "error: -s requires a path\n"); return 1; }
            path = argv[++i];
        } else if (strcmp(argv[i], "-y") == 0 || strcmp(argv[i], "--yes") == 0) {
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
    int pwd_len = bge_read_secret("Enter passphrase: ", pwd, sizeof(pwd));
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
    int pwd_len = bge_read_secret("Enter passphrase: ", pwd, sizeof(pwd));
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
    while (p + 2 <= end) {
        uint16_t klen = cfx_load16_le(p);
        p += 2;
        if (p + klen > end) break;
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
    int pwd_len = bge_read_secret("Enter passphrase: ", pwd, sizeof(pwd));
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
    int old_len = bge_read_secret("Enter current passphrase: ", old_pwd, sizeof(old_pwd));
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
    int pwd_len = bge_read_secret("Enter passphrase: ", pwd, sizeof(pwd));
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

int cfx_bge_run(int argc, char **argv) {
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
