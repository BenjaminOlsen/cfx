/* cfx_keyfile.c - encrypted key file read/write (CFX\x01 format) */

#include "cfx_keyfile.h"
#include "cfx/argon2.h"
#include "cfx/aead_chacha20_poly1305.h"
#include "cfx/rand.h"
#include "cfx/memory.h"

#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>

#ifdef _WIN32
#include <windows.h>
#else
#include <termios.h>
#include <unistd.h>
#endif


int cfx_key_read_secret_console(const char *prompt, char *buf, size_t bufsz) {
#ifdef _WIN32
    HANDLE h = GetStdHandle(STD_INPUT_HANDLE);
    DWORD mode;
    BOOL have_console = GetConsoleMode(h, &mode);
    if (have_console) {
        SetConsoleMode(h, mode & ~ENABLE_ECHO_INPUT);
    }

    fprintf(stderr, "%s", prompt);

    if (!fgets(buf, (int)bufsz, stdin)) {
        buf[0] = '\0';
    }

    if (have_console) {
        SetConsoleMode(h, mode);
    }
    fprintf(stderr, "\n");

    size_t len = strlen(buf);
    while (len > 0 && (buf[len-1] == '\n' || buf[len-1] == '\r')) {
        buf[--len] = '\0';
    }

    return (int)len;
#else
    /* open the controlling terminal directly, so we still prompt the user even
       when stdin/stdout are pipes or files */
    FILE *tty = fopen("/dev/tty", "r+");
    if (!tty) return -1;

    fprintf(tty, "%s", prompt);
    fflush(tty);

    int tty_fd = fileno(tty);

    /* old = saved attrs to restore; noecho = working copy with ECHO cleared */
    struct termios old, noecho;

    /* if this fails we're not on a real tty, so we just skip echo suppression */
    int have_termios = (tcgetattr(tty_fd, &old) == 0);
    if (have_termios) {
        noecho = old;
        noecho.c_lflag &= ~(tcflag_t)ECHO; /* clear the ECHO bit so we dont see the characters */
        tcsetattr(tty_fd, TCSANOW, &noecho); /* TCSANOW: apply immediately, don't wait for the input queue to drain */
    }

    char *r = fgets(buf, (int)bufsz, tty);

    /* restore the original termios and emit a newline since their Enter keystroke was'nt echoed */
    if (have_termios) {
        tcsetattr(tty_fd, TCSANOW, &old);
        fprintf(tty, "\n");
    }

    fclose(tty);

    if (!r) {
        cfx_memzero_s(buf, bufsz);
        return 0;
    }

    size_t len = strlen(buf);
    while (len > 0 && (buf[len-1] == '\n' || buf[len-1] == '\r')) {
        buf[--len] = '\0';
    }

    return (int)len;

#endif
}

int cfx_key_decrypt(const uint8_t *file_buf, uint8_t *seed_out) {
    if (memcmp(file_buf, CFX_KEY_MAGIC, 4) != 0) {
        fprintf(stderr, "error: not an encrypted key file\n");
        return -1;
    }

    const uint8_t *header = file_buf;
    uint32_t m_cost = cfx_load32_le(header + 4);
    uint32_t t_cost = cfx_load32_le(header + 8);
    uint32_t p      = cfx_load32_le(header + 12);

    if (m_cost < 8 || t_cost < 1 || p < 1 ||
        m_cost > 4194304 || t_cost > 100 || p > 16) {
        fprintf(stderr, "error: unreasonable argon2 parameters in key file\n");
        return -1;
    }

    const uint8_t *salt  = header + 16;
    const uint8_t *nonce = header + 32;
    const uint8_t *ct    = file_buf + 56;
    const uint8_t *tag   = file_buf + 88;

    char pwd[256] = {0};
    int pwd_len = cfx_key_read_secret_console("Enter passphrase: ", pwd, sizeof(pwd));
    if (pwd_len <= 0) {
        fprintf(stderr, "error: encrypted key requires a passphrase\n");
        cfx_memzero_s(pwd, sizeof(pwd));
        return -1;
    }

    /* derive decryption key with argon2id */
    uint8_t dec_key[32];
    int rc = cfx_argon2id(dec_key, 32,
        (const uint8_t *)pwd, (size_t)pwd_len, salt, 16,
        m_cost, t_cost, p);
    cfx_memzero_s(pwd, sizeof(pwd));

    if (rc != 0) {
        fprintf(stderr, "error: argon2 key derivation failed\n");
        cfx_memzero_s(dec_key, sizeof(dec_key));
        return -1;
    }

    /* decrypt with xchacha20-poly1305, AAD = header bytes */
    rc = cfx_xchacha20_poly1305_decrypt(
        seed_out, ct, 32,
        header, CFX_KEY_HEADER_LEN,
        dec_key, nonce, tag);

    cfx_memzero_s(dec_key, sizeof(dec_key));

    if (rc != 0) {
        fprintf(stderr, "error: decryption failed (wrong passphrase?)\n");
        return -1;
    }

    return 0;
}

int cfx_key_write_encrypted(const char *path, const uint8_t *seed, size_t seed_len, const char *pwd, size_t pwd_len) {
    uint8_t file_buf[CFX_KEY_FILE_LEN];
    uint8_t enc_key[32];
    uint8_t *header = file_buf;
    uint8_t *salt   = file_buf + 16;
    uint8_t *nonce  = file_buf + 32;
    uint8_t *ct     = file_buf + 56;
    uint8_t *tag    = file_buf + 88;
    int ret = -1;

    memcpy(header, CFX_KEY_MAGIC, 4);
    cfx_store32_le(header + 4,  CFX_KEY_DEFAULT_M);
    cfx_store32_le(header + 8,  CFX_KEY_DEFAULT_T);
    cfx_store32_le(header + 12, CFX_KEY_DEFAULT_P);

    cfx_rand_bytes(salt, 16);
    cfx_rand_bytes(nonce, 24);

    /* derive encryption key */
    int rc = cfx_argon2id(enc_key, 32,
        (const uint8_t *)pwd, pwd_len, salt, 16,
        CFX_KEY_DEFAULT_M, CFX_KEY_DEFAULT_T, CFX_KEY_DEFAULT_P);
    if (rc != 0) {
        fprintf(stderr, "error: argon2 key derivation failed\n");
        goto done;
    }

    /* encrypt seed, AAD = header (magic + params + salt + nonce) */
    if (cfx_xchacha20_poly1305_encrypt(ct, tag, seed, seed_len, header, CFX_KEY_HEADER_LEN, enc_key, nonce) != 0) {
        fprintf(stderr, "error: encryption failed\n");
        goto done;
    }

    FILE *f = fopen(path, "wb");
    if (!f) {
        fprintf(stderr, "error: cannot create %s: %s\n", path, strerror(errno));
        goto done;
    }
    if (fwrite(file_buf, 1, CFX_KEY_FILE_LEN, f) != CFX_KEY_FILE_LEN) {
        fprintf(stderr, "error: write failed: %s\n", strerror(errno));
        fclose(f);
        goto done;
    }
    fclose(f);

#ifndef _WIN32
    if (chmod(path, 0600) != 0) {
        fprintf(stderr, "warning: cannot set permissions on %s: %s\n", path, strerror(errno));
    }
#endif
    ret = 0;

done:
    cfx_memzero_s(enc_key, sizeof(enc_key));
    cfx_memzero_s(file_buf, sizeof(file_buf));
    return ret;
}
