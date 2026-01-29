/* cfx_bge.c - encrypted secrets store (BGE format) */

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
#include <pwd.h>
#include <termios.h>
#define mkdir_p(path) mkdir(path, 0700)
#endif

#include "cfx_cmd.h"
#include "cfx_keyfile.h"
#include "misc.h"

/*
 * BGE encrypted store layout:
 *
 *   Offset  Size  Field
 *   0       4     magic: "BGE\x01"
 *   4       4     argon2 m_cost  (LE32)
 *   8       4     argon2 t_cost  (LE32)
 *   12      4     argon2 parallelism (LE32)
 *   16      16    salt
 *   32      24    nonce
 *   56      ...   ciphertext (variable length)
 *   -16     16    auth tag (Poly1305)
 *
 * The header (bytes 0-55) is used as AAD for the AEAD.
 *
 * Plaintext is a binary TLV store:
 *   entry = uint16_le(name_len) | name_bytes | uint32_le(val_len) | val_bytes
 *   store = entry entry entry ...
 */
#define BGE_MAGIC       "BGE\x01"
#define BGE_HEADER_LEN  56
#define BGE_TAG_LEN     16
#define BGE_MIN_FILE    (BGE_HEADER_LEN + BGE_TAG_LEN)  /* 72: empty store */

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
    printf("  get    <name> [-s path]          Retrieve a secret (raw to stdout)\n");
    printf("  set    <name> [value] [-s path]  Set a secret (value from arg, -f, or stdin)\n");
    printf("  rm     <name> [-s path]          Remove a secret\n");
    printf("  ls     [-s path]                 List all secret names\n");
    printf("  info   [-s path]                 Show store location, size, and entry count\n");
    printf("  passwd [-s path]                 Change passphrase\n");
    printf("  dump   [-s path]                 Print all secrets to stdout\n");
    printf("\nOptions:\n");
    printf("  -s, --store <path>   Path to BGE store (default: ~/.cfx/secrets.bge)\n");
    printf("\nOptions for set:\n");
    printf("  -f <file>   Read value from file (binary-safe)\n");
    printf("              If no value and no -f, reads from stdin.\n");
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

/* ── path helpers  */
#define SUBDIR "cfx"

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
    int n = snprintf(buf, bufsz, "%s/.%s", home, SUBDIR);
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

/* ── uncapped file/stream reader  */

/*
 * read all data from an open FILE* into a malloc'd buf
 * works for seekable files and non-seekable streams (stdin, pipes).
 * cap: BGE_READ_MAX.  Returns 0 on success.
 */
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

/*
 * Read password from /dev/tty so that stdin remains free for piped data.
 * Falls back to cfx_key_read_password (stdin) on Windows or if /dev/tty
 * cannot be opened.
 */
static int bge_read_password(const char *prompt, char *buf, size_t bufsz) {
#ifndef _WIN32
    FILE *tty = fopen("/dev/tty", "r+");
    if (tty) {
        fprintf(tty, "%s", prompt);
        fflush(tty);

        /* disable echo */
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
        while (len > 0 && (buf[len-1] == '\n' || buf[len-1] == '\r'))
            buf[--len] = '\0';
        return len;
    }
#endif
    return cfx_key_read_password(prompt, buf, bufsz);
}

static int ct_pwd_match(const char *pw1, int pw1_len, size_t pw1_bufsz, const char *pw2, int pw2_len, size_t pw2_bufsz) {
    size_t n = pw1_bufsz < pw2_bufsz ? pw1_bufsz : pw2_bufsz;
    int diff = pw1_len ^ pw2_len;
    for (size_t i = 0; i < n; ++i) {
        diff |= pw1[i] ^ pw2[i];
    }
    return diff == 0;
}

/* prompt for passphrase (twice), return length. -1 on mismatch/error. */
static int prompt_passphrase(char *pwd, size_t pwdsz) {
    char pwd2[256] = {0};
    int len = bge_read_password("Enter passphrase: ", pwd, pwdsz);
    if (len <= 0) {
        fprintf(stderr, "error: passphrase required\n");
        cfx_memzero_s(pwd, pwdsz);
        return -1;
    }

    int len2 = bge_read_password("Enter same passphrase again: ", pwd2, sizeof(pwd2));
    if (!ct_pwd_match(pwd, len, pwdsz, pwd2, len2, sizeof(pwd2))) {
        fprintf(stderr, "Passphrases do not match.\n");
        cfx_memzero_s(pwd, pwdsz);
        cfx_memzero_s(pwd2, sizeof(pwd2));
        return -1;
    }
    cfx_memzero_s(pwd2, sizeof(pwd2));
    return len;
}

/* ── crypto core ─ */

/*
 * Decrypt a BGE file. Returns 0 on success.
 *   *pt_out  = malloc'd plaintext (caller frees after wiping)
 *   *pt_len  = plaintext length
 *   key_out  = derived 32-byte key (for re-encryption without re-running Argon2)
 *   hdr_out  = copy of 56-byte header (for re-encryption)
 */
static int bge_decrypt(const char *path, const char *pwd, size_t pwd_len,
                       uint8_t **pt_out, size_t *pt_len,
                       uint8_t key_out[32], uint8_t hdr_out[BGE_HEADER_LEN]) {
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

    if (file_len < BGE_MIN_FILE || memcmp(file_buf, BGE_MAGIC, 4) != 0) {
        fprintf(stderr, "error: %s is not a valid BGE store\n", path);
        cfx_memzero_s(file_buf, file_len);
        free(file_buf);
        return -1;
    }

    const uint8_t *header = file_buf;
    uint32_t m_cost = cfx_load32_le(header + 4);
    uint32_t t_cost = cfx_load32_le(header + 8);
    uint32_t p      = cfx_load32_le(header + 12);

    if (m_cost < 8 || t_cost < 1 || p < 1 ||
        m_cost > 4194304 || t_cost > 100 || p > 16) {
        fprintf(stderr, "error: unreasonable argon2 parameters in %s\n", path);
        cfx_memzero_s(file_buf, file_len);
        free(file_buf);
        return -1;
    }

    const uint8_t *salt  = header + 16;
    const uint8_t *nonce = header + 32;
    size_t ct_len = file_len - BGE_HEADER_LEN - BGE_TAG_LEN;
    const uint8_t *ct    = file_buf + BGE_HEADER_LEN;
    const uint8_t *tag   = file_buf + file_len - BGE_TAG_LEN;

    /* derive decryption key */
    uint8_t dec_key[32];
    int rc = cfx_argon2id(dec_key, 32,
        (const uint8_t *)pwd, pwd_len, salt, 16,
        m_cost, t_cost, p);
    if (rc != 0) {
        fprintf(stderr, "error: argon2 key derivation failed\n");
        cfx_memzero_s(dec_key, sizeof(dec_key));
        cfx_memzero_s(file_buf, file_len);
        free(file_buf);
        return -1;
    }

    /* decrypt */
    uint8_t *plaintext = malloc(ct_len + 1); /* +1 for NUL convenience */
    if (!plaintext) {
        fprintf(stderr, "error: allocation failed\n");
        cfx_memzero_s(dec_key, sizeof(dec_key));
        cfx_memzero_s(file_buf, file_len);
        free(file_buf);
        return -1;
    }

    rc = cfx_xchacha20_poly1305_decrypt(
        plaintext, ct, ct_len,
        header, BGE_HEADER_LEN,
        dec_key, nonce, tag);

    if (rc != 0) {
        fprintf(stderr, "error: decryption failed (wrong passphrase?)\n");
        cfx_memzero_s(plaintext, ct_len + 1);
        free(plaintext);
        cfx_memzero_s(dec_key, sizeof(dec_key));
        cfx_memzero_s(file_buf, file_len);
        free(file_buf);
        return -1;
    }

    plaintext[ct_len] = '\0';

    /* copy key and header for caller */
    memcpy(key_out, dec_key, 32);
    memcpy(hdr_out, header, BGE_HEADER_LEN);

    *pt_out = plaintext;
    *pt_len = ct_len;

    cfx_memzero_s(dec_key, sizeof(dec_key));
    cfx_memzero_s(file_buf, file_len);
    free(file_buf);
    return 0;
}

/*
 * Encrypt plaintext and write a new BGE file with fresh salt + nonce.
 */
static int bge_encrypt_write(const char *path, const uint8_t *pt, size_t pt_len,
                             const char *pwd, size_t pwd_len,
                             uint32_t m, uint32_t t, uint32_t p) {
    uint8_t header[BGE_HEADER_LEN];
    uint8_t enc_key[32];
    int ret = -1;

    /* build header */
    memcpy(header, BGE_MAGIC, 4);
    cfx_store32_le(header + 4,  m);
    cfx_store32_le(header + 8,  t);
    cfx_store32_le(header + 12, p);

    cfx_srand_os();
    cfx_rand_bytes(header + 16, 16);  /* salt */
    cfx_rand_bytes(header + 32, 24);  /* nonce */

    /* derive key */
    int rc = cfx_argon2id(enc_key, 32,
        (const uint8_t *)pwd, pwd_len, header + 16, 16, m, t, p);
    if (rc != 0) {
        fprintf(stderr, "error: argon2 key derivation failed\n");
        goto done;
    }

    /* allocate ciphertext + tag */
    uint8_t *ct = malloc(pt_len);
    uint8_t tag[BGE_TAG_LEN];
    if (pt_len > 0 && !ct) {
        fprintf(stderr, "error: allocation failed\n");
        goto done;
    }

    rc = cfx_xchacha20_poly1305_encrypt(
        ct, tag, pt, pt_len,
        header, BGE_HEADER_LEN,
        enc_key, header + 32);
    if (rc != 0) {
        fprintf(stderr, "error: encryption failed\n");
        free(ct);
        goto done;
    }

    /* write file */
    FILE *f = fopen(path, "wb");
    if (!f) {
        fprintf(stderr, "error: cannot create %s: %s\n", path, strerror(errno));
        free(ct);
        goto done;
    }

    size_t ok = 1;
    ok = ok && fwrite(header, 1, BGE_HEADER_LEN, f) == BGE_HEADER_LEN;
    if (pt_len > 0)
        ok = ok && fwrite(ct, 1, pt_len, f) == pt_len;
    ok = ok && fwrite(tag, 1, BGE_TAG_LEN, f) == BGE_TAG_LEN;
    fclose(f);

    if (!ok) {
        fprintf(stderr, "error: write failed: %s\n", strerror(errno));
        free(ct);
        goto done;
    }

#ifndef _WIN32
    if (chmod(path, 0600) != 0) {
        fprintf(stderr, "warning: cannot set permissions on %s: %s\n", path, strerror(errno));
    }
#endif

    cfx_memzero_s(ct, pt_len);
    free(ct);
    ret = 0;

done:
    cfx_memzero_s(enc_key, sizeof(enc_key));
    cfx_memzero_s(header, sizeof(header));
    return ret;
}

/*
 * Re-encrypt with existing derived key but a fresh nonce.
 * Used by set/rm/edit to avoid re-running Argon2id.
 */
static int bge_write_reusing_key(const char *path, const uint8_t *pt, size_t pt_len,
                                 const uint8_t key[32],
                                 const uint8_t old_header[BGE_HEADER_LEN]) {
    uint8_t header[BGE_HEADER_LEN];
    int ret = -1;

    /* copy old header (preserves magic, params, salt) */
    memcpy(header, old_header, BGE_HEADER_LEN);

    /* generate fresh nonce only */
    cfx_srand_os();
    cfx_rand_bytes(header + 32, 24);

    /* allocate ciphertext + tag */
    uint8_t *ct = malloc(pt_len);
    uint8_t tag[BGE_TAG_LEN];
    if (pt_len > 0 && !ct) {
        fprintf(stderr, "error: allocation failed\n");
        return -1;
    }

    int rc = cfx_xchacha20_poly1305_encrypt(
        ct, tag, pt, pt_len,
        header, BGE_HEADER_LEN,
        key, header + 32);
    if (rc != 0) {
        fprintf(stderr, "error: encryption failed\n");
        free(ct);
        goto done;
    }

    /* write file */
    FILE *f = fopen(path, "wb");
    if (!f) {
        fprintf(stderr, "error: cannot create %s: %s\n", path, strerror(errno));
        free(ct);
        goto done;
    }

    size_t ok = 1;
    ok = ok && fwrite(header, 1, BGE_HEADER_LEN, f) == BGE_HEADER_LEN;
    if (pt_len > 0)
        ok = ok && fwrite(ct, 1, pt_len, f) == pt_len;
    ok = ok && fwrite(tag, 1, BGE_TAG_LEN, f) == BGE_TAG_LEN;
    fclose(f);

    if (!ok) {
        fprintf(stderr, "error: write failed: %s\n", strerror(errno));
        free(ct);
        goto done;
    }

#ifndef _WIN32
    if (chmod(path, 0600) != 0) {
        fprintf(stderr, "warning: cannot set permissions on %s: %s\n", path, strerror(errno));
    }
#endif

    cfx_memzero_s(ct, pt_len);
    free(ct);
    ret = 0;

done:
    cfx_memzero_s(header, sizeof(header));
    return ret;
}

/* ── TLV store manipulation ───── */

/*
 * TLV entry format:
 *   uint16_le(name_len) | name_bytes | uint32_le(val_len) | val_bytes
 *
 * store = entry entry entry ...
 */

/*
 * Find entry "name" in TLV store.
 * Returns pointer to value start within pt, sets *vlen to value length.
 * Returns NULL if not found.
 */
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
            *vlen = val_len;
            return vptr;
        }
    }
    return NULL;
}

/*
 * Insert or update entry in TLV store.
 * Returns new malloc'd buffer, sets *new_len.
 * val/val_len is an arbitrary byte buffer.
 */
static uint8_t *store_set(const uint8_t *pt, size_t pt_len,
                          const char *name, const uint8_t *val, size_t val_len,
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
        /* replace existing entry */
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
        /* append new entry */
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

/*
 * Remove entry "name" from TLV store.
 * Returns new malloc'd buffer, sets *new_len.
 * Returns NULL if not found or allocation fails.
 */
static uint8_t *store_rm(const uint8_t *pt, size_t pt_len,
                         const char *name, size_t *new_len) {
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

/* ── TLV <-> text conversion for dump/edit ───────── */

/*
 * Convert TLV store to human-readable text:
 *   [keyname]\n
 *   value content\n
 *   [another_key]\n
 *   another value\n
 *
 * Round-trip rule: append \n after each value.
 * Returns malloc'd text buffer, sets *text_len.
 */
static uint8_t *store_to_text(const uint8_t *pt, size_t pt_len, size_t *text_len) {
    /* first pass: compute size */
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

    /* second pass: write */
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

/* ── subcommands ───────────────── */

static int bge_init(int argc, char **argv) {
    char path_buf[1024];
    const char *path = NULL;

    for (int i = 1; i < argc; i++) {
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

    /* ensure directory exists if using default path */
    char cfx_dir[1024];
    if (get_cfx_dir(cfx_dir, sizeof(cfx_dir)) == 0 &&
        strncmp(path, cfx_dir, strlen(cfx_dir)) == 0) {
        if (ensure_cfx_dir() != 0) return 1;
    }

    /* check file doesn't already exist */
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
    char path_buf[1024];

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            printf("Usage: %s get <name> [-s path]\n", argv[0]);
            return 0;
        } else if (strcmp(argv[i], "-s") == 0 || strcmp(argv[i], "--store") == 0) {
            if (i + 1 >= argc) { fprintf(stderr, "error: -s requires a path\n"); return 1; }
            path = argv[++i];
        } else if (!name) {
            name = argv[i];
        } else {
            fprintf(stderr, "error: unknown argument: %s\n", argv[i]);
            return 1;
        }
    }

    if (!name) {
        fprintf(stderr, "Usage: %s get <name> [-s path]\n", argv[0]);
        return 1;
    }

    if (!path) {
        if (bge_default_path(path_buf, sizeof(path_buf)) != 0) return 1;
        path = path_buf;
    }

    char pwd[256] = {0};
    int pwd_len = bge_read_password("Enter passphrase: ", pwd, sizeof(pwd));
    if (pwd_len <= 0) {
        fprintf(stderr, "error: passphrase required\n");
        cfx_memzero_s(pwd, sizeof(pwd));
        return 1;
    }

    uint8_t *pt = NULL;
    size_t pt_len = 0;
    uint8_t key[32], hdr[BGE_HEADER_LEN];
    int rc = bge_decrypt(path, pwd, (size_t)pwd_len, &pt, &pt_len, key, hdr);
    cfx_memzero_s(pwd, sizeof(pwd));
    cfx_memzero_s(key, sizeof(key));
    if (rc != 0) return 1;

    size_t vlen;
    const uint8_t *val = store_get(pt, pt_len, name, &vlen);
    if (!val) {
        fprintf(stderr, "error: '%s' not found\n", name);
        cfx_memzero_s(pt, pt_len);
        free(pt);
        return 1;
    }

    fwrite(val, 1, vlen, stdout);
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
    char path_buf[1024];

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            printf("Usage: %s set <name> [value] [-s path]\n", argv[0]);
            printf("\nOptions:\n");
            printf("  -s, --store <path>   Path to BGE store\n");
            printf("  -f <file>            Read value from file\n");
            printf("                       If no value and no -f, reads from stdin.\n");
            return 0;
        } else if (strcmp(argv[i], "-s") == 0 || strcmp(argv[i], "--store") == 0) {
            if (i + 1 >= argc) { fprintf(stderr, "error: -s requires a path\n"); return 1; }
            path = argv[++i];
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

    if (!name) {
        fprintf(stderr, "Usage: %s set <name> [value] [-s path]\n", argv[0]);
        return 1;
    }

    if (!path) {
        if (bge_default_path(path_buf, sizeof(path_buf)) != 0) return 1;
        path = path_buf;
    }

    /* determine value bytes */
    uint8_t *val_buf = NULL;
    size_t val_len = 0;
    int val_needs_free = 0;

    if (file_arg) {
        /* mode 1: read value from file */
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
    } else if (value_str) {
        /* mode 2: string from argv */
        val_buf = (uint8_t *)value_str;
        val_len = strlen(value_str);
    } else {
        /* mode 3: read from stdin */
#ifndef _WIN32
        if (isatty(STDIN_FILENO)) {
            /* interactive: read secret with echo disabled */
            char secret_buf[4096] = {0};
            int slen = bge_read_password("Enter secret value: ", secret_buf, sizeof(secret_buf));
            if (slen <= 0) {
                fprintf(stderr, "error: empty value\n");
                cfx_memzero_s(secret_buf, sizeof(secret_buf));
                return 1;
            }
            val_buf = malloc((size_t)slen);
            if (!val_buf) {
                cfx_memzero_s(secret_buf, sizeof(secret_buf));
                return 1;
            }
            memcpy(val_buf, secret_buf, (size_t)slen);
            val_len = (size_t)slen;
            cfx_memzero_s(secret_buf, sizeof(secret_buf));
            val_needs_free = 1;
        } else
#endif
        {
            /* pipe/redirect: read all of stdin */
            if (bge_read_all(stdin, &val_buf, &val_len) != 0) {
                fprintf(stderr, "error: cannot read from stdin\n");
                return 1;
            }
            val_needs_free = 1;
        }
    }

    char pwd[256] = {0};
    int pwd_len = bge_read_password("Enter passphrase: ", pwd, sizeof(pwd));
    if (pwd_len <= 0) {
        fprintf(stderr, "error: passphrase required\n");
        cfx_memzero_s(pwd, sizeof(pwd));
        if (val_needs_free) {
            cfx_memzero_s(val_buf, val_len);
            free(val_buf);
        }
        return 1;
    }

    uint8_t *pt = NULL;
    size_t pt_len = 0;
    uint8_t key[32], hdr[BGE_HEADER_LEN];
    int rc = bge_decrypt(path, pwd, (size_t)pwd_len, &pt, &pt_len, key, hdr);
    cfx_memzero_s(pwd, sizeof(pwd));
    if (rc != 0) {
        if (val_needs_free) {
            cfx_memzero_s(val_buf, val_len);
            free(val_buf);
        }
        return 1;
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
        cfx_memzero_s(key, sizeof(key));
        return 1;
    }

    rc = bge_write_reusing_key(path, new_pt, new_len, key, hdr);
    cfx_memzero_s(new_pt, new_len);
    free(new_pt);
    cfx_memzero_s(key, sizeof(key));
    return rc != 0;
}

static int bge_rm(int argc, char **argv) {
    const char *name = NULL;
    const char *path = NULL;
    char path_buf[1024];

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            printf("Usage: %s rm <name> [-s path]\n", argv[0]);
            return 0;
        } else if (strcmp(argv[i], "-s") == 0 || strcmp(argv[i], "--store") == 0) {
            if (i + 1 >= argc) { fprintf(stderr, "error: -s requires a path\n"); return 1; }
            path = argv[++i];
        } else if (!name) {
            name = argv[i];
        } else {
            fprintf(stderr, "error: unknown argument: %s\n", argv[i]);
            return 1;
        }
    }

    if (!name) {
        fprintf(stderr, "Usage: %s rm <name> [-s path]\n", argv[0]);
        return 1;
    }

    if (!path) {
        if (bge_default_path(path_buf, sizeof(path_buf)) != 0) return 1;
        path = path_buf;
    }

    char pwd[256] = {0};
    int pwd_len = bge_read_password("Enter passphrase: ", pwd, sizeof(pwd));
    if (pwd_len <= 0) {
        fprintf(stderr, "error: passphrase required\n");
        cfx_memzero_s(pwd, sizeof(pwd));
        return 1;
    }

    uint8_t *pt = NULL;
    size_t pt_len = 0;
    uint8_t key[32], hdr[BGE_HEADER_LEN];
    int rc = bge_decrypt(path, pwd, (size_t)pwd_len, &pt, &pt_len, key, hdr);
    cfx_memzero_s(pwd, sizeof(pwd));
    if (rc != 0) return 1;

    size_t new_len;
    uint8_t *new_pt = store_rm(pt, pt_len, name, &new_len);
    cfx_memzero_s(pt, pt_len);
    free(pt);

    if (!new_pt) {
        cfx_memzero_s(key, sizeof(key));
        return 1;
    }

    rc = bge_write_reusing_key(path, new_pt, new_len, key, hdr);
    cfx_memzero_s(new_pt, new_len);
    free(new_pt);
    cfx_memzero_s(key, sizeof(key));
    return rc != 0;
}

static int bge_ls(int argc, char **argv) {
    const char *path = NULL;
    char path_buf[1024];

    for (int i = 1; i < argc; i++) {
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
    int pwd_len = bge_read_password("Enter passphrase: ", pwd, sizeof(pwd));
    if (pwd_len <= 0) {
        fprintf(stderr, "error: passphrase required\n");
        cfx_memzero_s(pwd, sizeof(pwd));
        return 1;
    }

    uint8_t *pt = NULL;
    size_t pt_len = 0;
    uint8_t key[32], hdr[BGE_HEADER_LEN];
    int rc = bge_decrypt(path, pwd, (size_t)pwd_len, &pt, &pt_len, key, hdr);
    cfx_memzero_s(pwd, sizeof(pwd));
    cfx_memzero_s(key, sizeof(key));
    if (rc != 0) return 1;

    /* iterate TLV entries, print key names */
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

    for (int i = 1; i < argc; i++) {
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

    /* file metadata (no passphrase needed) */
    struct stat st;
    if (stat(path, &st) != 0) {
        fprintf(stderr, "error: cannot stat %s: %s\n", path, strerror(errno));
        return 1;
    }

    printf("Store:   %s\n", path);
    printf("Size:    %lld bytes\n", (long long)st.st_size);

    /* entry count requires decryption */
    char pwd[256] = {0};
    int pwd_len = bge_read_password("Enter passphrase: ", pwd, sizeof(pwd));
    if (pwd_len <= 0) {
        fprintf(stderr, "error: passphrase required\n");
        cfx_memzero_s(pwd, sizeof(pwd));
        return 1;
    }

    uint8_t *pt = NULL;
    size_t pt_len = 0;
    uint8_t key[32], hdr[BGE_HEADER_LEN];
    int rc = bge_decrypt(path, pwd, (size_t)pwd_len, &pt, &pt_len, key, hdr);
    cfx_memzero_s(pwd, sizeof(pwd));
    cfx_memzero_s(key, sizeof(key));
    if (rc != 0) return 1;

    /* count TLV entries */
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

    /* show argon2 params from header */
    uint32_t m_cost = cfx_load32_le(hdr + 4);
    uint32_t t_cost = cfx_load32_le(hdr + 8);
    uint32_t p_cost = cfx_load32_le(hdr + 12);
    printf("Argon2:  m=%u KB, t=%u, p=%u\n", m_cost, t_cost, p_cost);

    return 0;
}

static int bge_passwd(int argc, char **argv) {
    const char *path = NULL;
    char path_buf[1024];

    for (int i = 1; i < argc; i++) {
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

    /* prompt for old passphrase */
    char old_pwd[256] = {0};
    int old_len = bge_read_password("Enter current passphrase: ", old_pwd, sizeof(old_pwd));
    if (old_len <= 0) {
        fprintf(stderr, "error: passphrase required\n");
        cfx_memzero_s(old_pwd, sizeof(old_pwd));
        return 1;
    }

    uint8_t *pt = NULL;
    size_t pt_len = 0;
    uint8_t key[32], hdr[BGE_HEADER_LEN];
    int rc = bge_decrypt(path, old_pwd, (size_t)old_len, &pt, &pt_len, key, hdr);
    cfx_memzero_s(old_pwd, sizeof(old_pwd));
    cfx_memzero_s(key, sizeof(key));
    if (rc != 0) return 1;

    /* prompt for new passphrase (double prompt) */
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

    for (int i = 1; i < argc; i++) {
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
    int pwd_len = bge_read_password("Enter passphrase: ", pwd, sizeof(pwd));
    if (pwd_len <= 0) {
        fprintf(stderr, "error: passphrase required\n");
        cfx_memzero_s(pwd, sizeof(pwd));
        return 1;
    }

    uint8_t *pt = NULL;
    size_t pt_len = 0;
    uint8_t key[32], hdr[BGE_HEADER_LEN];
    int rc = bge_decrypt(path, pwd, (size_t)pwd_len, &pt, &pt_len, key, hdr);
    cfx_memzero_s(pwd, sizeof(pwd));
    cfx_memzero_s(key, sizeof(key));
    if (rc != 0) return 1;

    /* convert TLV to text and print */
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

/* ── dispatcher ── */

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

    /* shift argc/argv so subcommand sees itself as argv[0] */
    int sub_argc = argc - 1;
    char **sub_argv = argv + 1;

    if (strcmp(cmd, "init")   == 0) return bge_init(sub_argc, sub_argv);
    if (strcmp(cmd, "get")    == 0) return bge_get(sub_argc, sub_argv);
    if (strcmp(cmd, "set")    == 0) return bge_set(sub_argc, sub_argv);
    if (strcmp(cmd, "rm")     == 0) return bge_rm(sub_argc, sub_argv);
    if (strcmp(cmd, "ls")     == 0 ||
        strcmp(cmd, "list")   == 0) return bge_ls(sub_argc, sub_argv);
    if (strcmp(cmd, "info")   == 0) return bge_info(sub_argc, sub_argv);
    if (strcmp(cmd, "passwd") == 0) return bge_passwd(sub_argc, sub_argv);
    if (strcmp(cmd, "dump")   == 0) return bge_dump(sub_argc, sub_argv);

    fprintf(stderr, "Unknown command: %s\n", cmd);
    usage(argv[0]);
    return 1;
}
