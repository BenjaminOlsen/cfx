/* cfx_bge_grace.c -- sudo-style auth grace period for cfx store
 *
 * The cached passphrase is encrypted at rest using XChaCha20-Poly1305.
 * The encryption key is derived from the system boot ID + UID, so:
 *   - the grace file is useless after a reboot (boot ID changes)
 *   - it can't be decrypted on a different machine
 *   - an offline attacker who copies ~/.cfx/.grace gets ciphertext, not plaintext
 */

#include "cfx_bge_internal.h"

#include <time.h>
#include <limits.h>

#ifndef _WIN32
#include <sys/types.h>
#endif

#ifdef __APPLE__
#include <sys/sysctl.h>
#endif

#define GRACE_WINDOW_SEC 60

#ifdef __linux__
  #define GRACE_CLOCK CLOCK_BOOTTIME
#else
  #define GRACE_CLOCK CLOCK_MONOTONIC
#endif

/* forward declaration — grace_check calls grace_delete on invalid files */
void grace_delete(void);

/*
 * Grace file layout (420 bytes total):
 *
 * AAD (authenticated, not encrypted):
 *   [0..7]     mono_ts       monotonic timestamp (8 bytes)
 *   [8..39]    store_hash    SHA-256 of canonical store path (32 bytes)
 *   [40..103]  tty           ttyname, NUL-padded (64 bytes)
 *
 * Encrypted payload:
 *   [104..127] nonce         XChaCha20-Poly1305 nonce (24 bytes)
 *   [128..385] ciphertext    encrypted(pwd_len[2] || pwd[256]) (258 bytes)
 *   [386..401] tag           Poly1305 tag (16 bytes)
 *
 * Encryption key: SHA-256(boot_id || uid_le32 || "cfx-grace")
 */

#define GRACE_AAD_LEN       104   /* ts + hash + tty */
#define GRACE_OFF_TS        0
#define GRACE_OFF_HASH      8
#define GRACE_OFF_TTY       40
#define GRACE_OFF_NONCE     104
#define GRACE_OFF_CT        128
#define GRACE_OFF_TAG       386
#define GRACE_PAYLOAD_LEN   258   /* pwd_len(2) + pwd(256) */
#define GRACE_FILE_LEN      402   /* 104 + 24 + 258 + 16 */

static int grace_path(char *buf, size_t bufsz) {
    char dir[1024];
    if (get_cfx_dir(dir, sizeof(dir)) != 0) return -1;
    int n = snprintf(buf, bufsz, "%s/.grace", dir);
    if (n < 0 || (size_t)n >= bufsz) return -1;
    return 0;
}

static void hash_store_path(const char *store_path, uint8_t out[32]) {
    char resolved[PATH_MAX];
    const char *canon = realpath(store_path, resolved);
    if (!canon) canon = store_path;
    cfx_sha256(out, (const uint8_t *)canon, strlen(canon));
}

static int get_tty(char *buf, size_t bufsz) {
#ifndef _WIN32
    const char *t = ttyname(STDIN_FILENO);
    if (!t) return -1;
    size_t len = strlen(t);
    if (len >= bufsz) return -1;
    memcpy(buf, t, len);
    /* rest is already zeroed by caller */
    return 0;
#else
    (void)buf; (void)bufsz;
    return -1;
#endif
}

static uint64_t mono_now(void) {
    struct timespec ts;
    if (clock_gettime(GRACE_CLOCK, &ts) != 0) return 0;
    return (uint64_t)ts.tv_sec;
}

/* Derive a per-boot, per-user grace encryption key from the system boot ID.
 * Returns 0 on success, -1 if boot ID is unavailable. */
static int grace_derive_key(uint8_t key_out[32]) {
    uint8_t boot_id[64];
    size_t boot_id_len = 0;

#ifdef __APPLE__
    /* macOS: kern.bootsessionuuid */
    size_t len = sizeof(boot_id);
    if (sysctlbyname("kern.bootsessionuuid", boot_id, &len, NULL, 0) != 0)
        return -1;
    boot_id_len = len;
#elif defined(__linux__)
    /* Linux: /proc/sys/kernel/random/boot_id */
    int fd = open("/proc/sys/kernel/random/boot_id", O_RDONLY);
    if (fd < 0) return -1;
    ssize_t rd = read(fd, boot_id, sizeof(boot_id));
    close(fd);
    if (rd <= 0) return -1;
    boot_id_len = (size_t)rd;
#else
    return -1;
#endif

    /* key = SHA-256(boot_id || uid_le32 || "cfx-grace") */
    cfx_sha256_ctx ctx;
    cfx_sha256_init(&ctx);
    cfx_sha256_update(&ctx, boot_id, boot_id_len);
    cfx_memzero_s(boot_id, sizeof(boot_id));

#ifndef _WIN32
    uint32_t uid = (uint32_t)getuid();
    uint8_t uid_le[4];
    uid_le[0] = (uint8_t)(uid);
    uid_le[1] = (uint8_t)(uid >> 8);
    uid_le[2] = (uint8_t)(uid >> 16);
    uid_le[3] = (uint8_t)(uid >> 24);
    cfx_sha256_update(&ctx, uid_le, 4);
#endif

    cfx_sha256_update(&ctx, (const uint8_t *)"cfx-grace", 9);
    cfx_sha256_final(&ctx, key_out);
    return 0;
}

int grace_check(const char *store_path, char *pwd_out, size_t pwd_bufsz) {
    char gpath[1024];
    if (grace_path(gpath, sizeof(gpath)) != 0) return 0;

#ifndef _WIN32
    /* lstat to reject symlinks */
    struct stat st;
    if (lstat(gpath, &st) != 0) return 0;
    if (!S_ISREG(st.st_mode)) { grace_delete(); return 0; }
    if ((st.st_mode & 0777) != 0600) { grace_delete(); return 0; }
    if (st.st_uid != getuid()) { grace_delete(); return 0; }

    int fd = open(gpath, O_RDONLY | O_NOFOLLOW);
    if (fd < 0) return 0;

    uint8_t buf[GRACE_FILE_LEN];
    ssize_t rd = read(fd, buf, sizeof(buf));
    close(fd);
    if (rd != GRACE_FILE_LEN) { grace_delete(); return 0; }
#else
    return 0;  /* no grace on Windows */
#endif

    /* check timestamp */
    uint64_t saved_ts;
    memcpy(&saved_ts, buf + GRACE_OFF_TS, 8);
    uint64_t now = mono_now();
    if (now == 0 || saved_ts == 0) { cfx_memzero_s(buf, sizeof(buf)); return 0; }
    if (now < saved_ts || (now - saved_ts) > GRACE_WINDOW_SEC) {
        cfx_memzero_s(buf, sizeof(buf));
        grace_delete();
        return 0;
    }

    /* check store hash */
    uint8_t cur_hash[32];
    hash_store_path(store_path, cur_hash);
    if (memcmp(cur_hash, buf + GRACE_OFF_HASH, 32) != 0) {
        cfx_memzero_s(buf, sizeof(buf));
        return 0;
    }

    /* check tty */
    char cur_tty[64] = {0};
    if (get_tty(cur_tty, sizeof(cur_tty)) != 0) {
        cfx_memzero_s(buf, sizeof(buf));
        return 0;
    }
    if (memcmp(cur_tty, buf + GRACE_OFF_TTY, 64) != 0) {
        cfx_memzero_s(buf, sizeof(buf));
        return 0;
    }

    /* derive grace encryption key from boot ID */
    uint8_t gkey[32];
    if (grace_derive_key(gkey) != 0) {
        cfx_memzero_s(buf, sizeof(buf));
        return 0;
    }

    /* decrypt the passphrase payload */
    uint8_t plaintext[GRACE_PAYLOAD_LEN];
    int rc = cfx_xchacha20_poly1305_decrypt(
        plaintext,                              /* pt out */
        buf + GRACE_OFF_CT, GRACE_PAYLOAD_LEN,  /* ct, ct_len */
        buf, GRACE_AAD_LEN,                     /* AAD = ts + hash + tty */
        gkey, buf + GRACE_OFF_NONCE,            /* key, nonce */
        buf + GRACE_OFF_TAG);                   /* tag */

    cfx_memzero_s(gkey, sizeof(gkey));
    cfx_memzero_s(buf, sizeof(buf));

    if (rc != 0) {
        /* decryption failed — stale key (reboot?) or tampered */
        grace_delete();
        cfx_memzero_s(plaintext, sizeof(plaintext));
        return 0;
    }

    /* extract passphrase from decrypted payload */
    uint16_t pwd_len;
    memcpy(&pwd_len, plaintext, 2);
    if (pwd_len == 0 || pwd_len > 256 || (size_t)pwd_len >= pwd_bufsz) {
        cfx_memzero_s(plaintext, sizeof(plaintext));
        grace_delete();
        return 0;
    }

    memcpy(pwd_out, plaintext + 2, pwd_len);
    cfx_memzero_s(plaintext, sizeof(plaintext));
    return (int)pwd_len;
}

void grace_stamp(const char *store_path, const char *pwd, size_t pwd_len) {
    if (pwd_len == 0 || pwd_len > 256) return;

#ifdef _WIN32
    return;
#else
    char tty[64] = {0};
    if (get_tty(tty, sizeof(tty)) != 0) return;

    /* derive grace encryption key from boot ID */
    uint8_t gkey[32];
    if (grace_derive_key(gkey) != 0) return;

    uint8_t buf[GRACE_FILE_LEN];
    memset(buf, 0, sizeof(buf));

    /* AAD: timestamp + store hash + tty */
    uint64_t now = mono_now();
    if (now == 0) { cfx_memzero_s(gkey, sizeof(gkey)); return; }
    memcpy(buf + GRACE_OFF_TS, &now, 8);
    hash_store_path(store_path, buf + GRACE_OFF_HASH);
    memcpy(buf + GRACE_OFF_TTY, tty, 64);

    /* build plaintext payload: pwd_len(2) || pwd(256) */
    uint8_t plaintext[GRACE_PAYLOAD_LEN];
    memset(plaintext, 0, sizeof(plaintext));
    uint16_t plen = (uint16_t)pwd_len;
    memcpy(plaintext, &plen, 2);
    memcpy(plaintext + 2, pwd, pwd_len);

    /* generate random nonce */
    cfx_srand_os();
    cfx_rand_bytes(buf + GRACE_OFF_NONCE, 24);

    /* encrypt passphrase payload, AAD = first 104 bytes */
    cfx_xchacha20_poly1305_encrypt(
        buf + GRACE_OFF_CT,                     /* ciphertext out */
        buf + GRACE_OFF_TAG,                    /* tag out */
        plaintext, GRACE_PAYLOAD_LEN,           /* plaintext, pt_len */
        buf, GRACE_AAD_LEN,                     /* AAD = ts + hash + tty */
        gkey, buf + GRACE_OFF_NONCE);           /* key, nonce */

    cfx_memzero_s(plaintext, sizeof(plaintext));
    cfx_memzero_s(gkey, sizeof(gkey));

    /* atomic write via .grace.tmp + fsync + rename */
    char gpath[1024], tmppath[1024];
    if (grace_path(gpath, sizeof(gpath)) != 0) {
        cfx_memzero_s(buf, sizeof(buf));
        return;
    }
    int n = snprintf(tmppath, sizeof(tmppath), "%s.tmp", gpath);
    if (n < 0 || (size_t)n >= sizeof(tmppath)) {
        cfx_memzero_s(buf, sizeof(buf));
        return;
    }

    int fd = open(tmppath, O_WRONLY | O_CREAT | O_TRUNC, 0600);
    if (fd < 0) {
        cfx_memzero_s(buf, sizeof(buf));
        return;
    }

    ssize_t wr = write(fd, buf, GRACE_FILE_LEN);
    cfx_memzero_s(buf, sizeof(buf));

    if (wr != GRACE_FILE_LEN) {
        close(fd);
        unlink(tmppath);
        return;
    }

    fsync(fd);
    close(fd);
    rename(tmppath, gpath);
#endif
}

void grace_delete(void) {
    char gpath[1024];
    if (grace_path(gpath, sizeof(gpath)) != 0) return;
    unlink(gpath);
}
