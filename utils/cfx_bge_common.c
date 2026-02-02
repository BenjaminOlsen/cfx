/* cfx_bge_common.c -- shared I/O, path, and utility helpers for BGE/store */

#include "cfx_bge_internal.h"

const char *g_passphrase_arg;

#define SUBDIR ".cfx"

int get_cfx_dir(char *buf, size_t bufsz) {
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

int bge_default_path(char *buf, size_t bufsz) {
    char cfx_dir[1024];
    if (get_cfx_dir(cfx_dir, sizeof(cfx_dir)) != 0) return -1;
    int n = snprintf(buf, bufsz, "%s/secrets.bge", cfx_dir);
    if (n < 0 || (size_t)n >= bufsz) {
        fprintf(stderr, "error: path too long\n");
        return -1;
    }
    return 0;
}

int ensure_cfx_dir(void) {
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



int bge_read_all(FILE *f, uint8_t **out, size_t *out_len) {
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

int bge_read_secret(const char *prompt, char *buf, size_t bufsz) {
#ifndef _WIN32
    FILE *tty = fopen("/dev/tty", "r+");
    if (tty) {
        fprintf(tty, "%s", prompt);
        fflush(tty);

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

int bge_read_passphrase(const char *prompt, char *buf, size_t bufsz) {
    if (g_passphrase_arg) {
        size_t len = strlen(g_passphrase_arg);
        if (len >= bufsz) len = bufsz - 1;
        memcpy(buf, g_passphrase_arg, len);
        buf[len] = '\0';
        return (int)len;
    }
    return bge_read_secret(prompt, buf, bufsz);
}

int bge_read_visible(const char *prompt, char *buf, size_t bufsz) {
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

int ct_pwd_match(const char *pw1, int pw1_len, size_t pw1_bufsz,
                 const char *pw2, int pw2_len, size_t pw2_bufsz) {
    size_t n = pw1_bufsz < pw2_bufsz ? pw1_bufsz : pw2_bufsz;
    int diff = pw1_len ^ pw2_len;
    for (size_t i = 0; i < n; ++i) {
        diff |= pw1[i] ^ pw2[i];
    }
    return diff == 0;
}

int prompt_passphrase(char *pwd, size_t pwdsz) {
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

/*  store display helpers  */

int is_numeric(const char *s) {
    if (!s || !*s) return 0;
    for (; *s; s++) {
        if (*s < '0' || *s > '9') return 0;
    }
    return 1;
}

void store_print_names(const uint8_t *pt, size_t pt_len) {
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
}
