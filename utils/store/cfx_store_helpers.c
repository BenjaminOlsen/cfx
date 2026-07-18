/* Store-specific paths, name parsing, and clipboard helpers. */

#include "cfx_store_internal.h"

#include <errno.h>

#ifdef _WIN32
#include <direct.h>
#include <shlobj.h>
#define mkdir_p(path) _mkdir(path)
#ifndef S_ISDIR
#define S_ISDIR(m) (((m) & _S_IFMT) == _S_IFDIR)
#endif
#else
#include <pwd.h>
#include <unistd.h>
#define mkdir_p(path) mkdir(path, 0700)
#endif

#define CFX_SUBDIR ".cfx"

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
    int n = snprintf(buf, bufsz, "%s/%s", home, CFX_SUBDIR);
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

int is_numeric(const char *s) {
    if (!s || !*s) return 0;
    for (; *s; s++) {
        if (*s < '0' || *s > '9') return 0;
    }
    return 1;
}

int clip_copy(const uint8_t *data, size_t len) {
    if (!data || len == 0) return -1;

#ifdef _WIN32
    FILE *p = _popen("clip.exe", "wb");
    if (!p) {
        fprintf(stderr, "error: cannot open clip.exe\n");
        return -1;
    }
    fwrite(data, 1, len, p);
    int rc = _pclose(p);
    return rc == 0 ? 0 : -1;
#elif defined(__APPLE__)
    FILE *p = popen("pbcopy", "w");
    if (!p) {
        fprintf(stderr, "error: cannot open pbcopy\n");
        return -1;
    }
    fwrite(data, 1, len, p);
    int rc = pclose(p);
    return rc == 0 ? 0 : -1;
#else
    const char *cmds[] = {
        "wl-copy",
        "xclip -selection clipboard",
        "xsel --clipboard --input",
    };
    const char *checks[] = {
        "command -v wl-copy >/dev/null 2>&1",
        "command -v xclip >/dev/null 2>&1",
        "command -v xsel >/dev/null 2>&1",
    };
    const char *cmd = NULL;
    for (int i = 0; i < 3; i++) {
        if (system(checks[i]) == 0) {
            cmd = cmds[i];
            break;
        }
    }
    if (!cmd) {
        fprintf(stderr, "error: no clipboard tool found (install wl-copy, xclip, or xsel)\n");
        return -1;
    }
    FILE *p = popen(cmd, "w");
    if (!p) {
        fprintf(stderr, "error: cannot run %s\n", cmd);
        return -1;
    }
    fwrite(data, 1, len, p);
    int rc = pclose(p);
    return rc == 0 ? 0 : -1;
#endif
}
