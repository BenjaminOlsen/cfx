/* cfx_store.c -- passphrase-protected secrets store subcommands */

#include "cfx_bge_internal.h"

static void usage(const char *prog) {
    printf("Usage: %s <command> [args...]\n\n", prog);
    printf("Passphrase-protected secrets store using\n");
    printf("Argon2id key derivation and XChaCha20-Poly1305 AEAD.\n\n");
    printf("Commands:\n");
    printf("  init   [-s path]                 Create a new encrypted store\n");
    printf("  get    <name> [-s path] [-q] [-c] Retrieve a secret (raw to stdout)\n");
    printf("  set    [name] [value] [-s path]  Set a secret (value from arg, -i, or stdin)\n");
    printf("  edit   <name> [-s path]          Edit secret in $EDITOR (creates if new)\n");
    printf("  rm     <name> [-s path]          Remove a secret\n");
    printf("  rename <old> <new> [-s path]     Rename a secret entry\n");
    printf("  ls     [-s path]                 List all secret names\n");
    printf("  info   [-s path]                 Show store location, size, and entry count\n");
    printf("  passwd [-s path]                 Change passphrase (current slot)\n");
    printf("  dump   [-s path]                 Print all secrets to stdout\n");
    printf("\nMulti-password (slot) commands:\n");
    printf("  slot ls  [-s path]               List slots and Argon2 params (no passphrase)\n");
    printf("  slot add [-s path]               Add a passphrase slot (migrates v2 to v4)\n");
    printf("  slot rm  [N] [-s path]           Remove slot N (1-based, default: matched)\n");
    printf("  rekey    [-s path]               Generate new DEK, re-wrap all slots\n");
    printf("\nCommon options:\n");
    printf("  -p, --passphrase <pw>  Supply passphrase on command line (default: prompt)\n");
    printf("  -s, --store <path>     Path to BGE store (default: ~/.cfx/secrets.bge)\n");
    printf("\nOptions for get:\n");
    printf("  -q, --quoted           Wrap output in single quotes\n");
    printf("  -c, --clipboard        Copy to system clipboard instead of stdout\n");
    printf("\nOptions for set:\n");
    printf("  -i <file>   Read value from file (binary-safe)\n");
    printf("              If no value and no -i, reads from stdin.\n");
    printf("  -f, --force  Overwrite existing entry without asking\n");
    printf("\nExamples:\n");
    printf("  %s init                              Create store at default path\n", prog);
    printf("  %s set token mytoken                 Set from command line argument\n", prog);
    printf("  %s get token                         Print value (newline on tty)\n", prog);
    printf("  %s ls                                Show all key names\n", prog);
    printf("  %s slot add                          Add a second passphrase\n", prog);
    printf("  %s slot ls                           List passphrase slots\n", prog);
}

static int store_init(int argc, char **argv) {
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

    char cfx_dir[1024];
    if (get_cfx_dir(cfx_dir, sizeof(cfx_dir)) == 0 &&
        strncmp(path, cfx_dir, strlen(cfx_dir)) == 0) {
        if (ensure_cfx_dir() != 0) return 1;
    }

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

static int store_cmd_get(int argc, char **argv) {
    const char *name = NULL;
    const char *path = NULL;
    int quoted = 0;
    int clipboard = 0;
    char path_buf[1024];

    for (int i = 2; i < argc; i++) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            printf("Usage: %s get <name> [-s path] [-q] [-c]\n", argv[0]);
            printf("\nOptions:\n");
            printf("  -s, --store <path>   Path to BGE store\n");
            printf("  -q, --quoted         Wrap output in single quotes\n");
            printf("  -c, --clipboard      Copy to system clipboard instead of stdout\n");
            return 0;
        } else if (strcmp(argv[i], "-s") == 0 || strcmp(argv[i], "--store") == 0) {
            if (i + 1 >= argc) { fprintf(stderr, "error: -s requires a path\n"); return 1; }
            path = argv[++i];
        } else if (strcmp(argv[i], "-q") == 0 || strcmp(argv[i], "--quoted") == 0) {
            quoted = 1;
        } else if (strcmp(argv[i], "-c") == 0 || strcmp(argv[i], "--clipboard") == 0) {
            clipboard = 1;
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

    bge_ustore us = {0};
    int rc = bge_uauthenticate(path, pwd, (size_t)pwd_len, &us);
    cfx_memzero_s(pwd, sizeof(pwd));
    if (rc != 0) return 1;

    uint8_t *pt = NULL;
    size_t pt_len = 0;
    rc = bge_udecrypt(&us, &pt, &pt_len);
    bge_ustore_wipe(&us);
    if (rc != 0) return 1;

    char name_buf[256] = {0};
    if (!name) {
        store_print_names(pt, pt_len);
        int r = bge_read_visible("Name: ", name_buf, sizeof(name_buf));
        if (r <= 0) {
            fprintf(stderr, "error: name required\n");
            cfx_memzero_s(pt, pt_len);
            free(pt);
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

    if (clipboard) {
        if (clip_copy(val, vlen) != 0) {
            fprintf(stderr, "error: clipboard copy failed\n");
            cfx_memzero_s(pt, pt_len);
            free(pt);
            return 1;
        }
        printf("Copied to clipboard.\n");
    } else {
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
    }

    cfx_memzero_s(pt, pt_len);
    free(pt);
    return 0;
}

static int store_cmd_set(int argc, char **argv) {
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

    bge_ustore us = {0};
    int rc = bge_uauthenticate(path, pwd, (size_t)pwd_len, &us);
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
    rc = bge_udecrypt(&us, &pt, &pt_len);
    if (rc != 0) {
        bge_ustore_wipe(&us);
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
            bge_ustore_wipe(&us);
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
                bge_ustore_wipe(&us);
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
            int stdin_is_tty = 0;
#ifndef _WIN32
            stdin_is_tty = isatty(STDIN_FILENO);
#else
            stdin_is_tty = _isatty(_fileno(stdin));
#endif
            if (stdin_is_tty) {
                char secret_buf[4096] = {0};
                int slen = bge_read_secret("Value: ", secret_buf, sizeof(secret_buf));
                if (slen <= 0) {
                    fprintf(stderr, "error: empty value\n");
                    cfx_memzero_s(secret_buf, sizeof(secret_buf));
                    cfx_memzero_s(pt, pt_len);
                    free(pt);
                    bge_ustore_wipe(&us);
                    return 1;
                }
                val_buf = malloc((size_t)slen);
                if (!val_buf) {
                    cfx_memzero_s(secret_buf, sizeof(secret_buf));
                    cfx_memzero_s(pt, pt_len);
                    free(pt);
                    bge_ustore_wipe(&us);
                    return 1;
                }
                memcpy(val_buf, secret_buf, (size_t)slen);
                val_len = (size_t)slen;
                cfx_memzero_s(secret_buf, sizeof(secret_buf));
                val_needs_free = 1;
            } else {
                if (bge_read_all(stdin, &val_buf, &val_len) != 0) {
                    fprintf(stderr, "error: cannot read from stdin\n");
                    cfx_memzero_s(pt, pt_len);
                    free(pt);
                    bge_ustore_wipe(&us);
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
        bge_ustore_wipe(&us);
        return 1;
    }

    rc = bge_uwrite(path, new_pt, new_len, &us);
    cfx_memzero_s(new_pt, new_len);
    free(new_pt);
    bge_ustore_wipe(&us);
    if (rc == 0) printf("Ok.\n");
    return rc != 0;
}

static const char *get_editor(void) {
    const char *ed = getenv("VISUAL");
    if (ed && *ed) return ed;
    ed = getenv("EDITOR");
    if (ed && *ed) return ed;
#ifdef _WIN32
    return "notepad";
#else
    return "vi";
#endif
}

static int create_secure_tmpfile(char *path_out, size_t bufsz) {
#ifdef _WIN32
    char tmpdir[MAX_PATH];
    if (GetTempPathA(MAX_PATH, tmpdir) == 0) return -1;
    if (GetTempFileNameA(tmpdir, "cfx", 0, path_out) == 0) return -1;
    return 0;
#else
    const char *tmpdir = getenv("TMPDIR");
    if (!tmpdir) tmpdir = "/tmp";
    int n = snprintf(path_out, bufsz, "%s/cfx_edit_XXXXXX", tmpdir);
    if (n < 0 || (size_t)n >= bufsz) return -1;
    int fd = mkstemp(path_out);
    if (fd < 0) return -1;
    fchmod(fd, 0600);
    close(fd);
    return 0;
#endif
}

static void wipe_tmpfile(const char *path) {
    FILE *f = fopen(path, "r+b");
    if (f) {
        fseek(f, 0, SEEK_END);
        long sz = ftell(f);
        if (sz > 0) {
            fseek(f, 0, SEEK_SET);
            uint8_t *zeros = calloc(1, (size_t)sz);
            if (zeros) {
                fwrite(zeros, 1, (size_t)sz, f);
                fflush(f);
                free(zeros);
            }
        }
        fclose(f);
    }
    remove(path);
}

static int store_cmd_edit(int argc, char **argv) {
    const char *name = NULL;
    const char *path = NULL;
    char path_buf[1024];

    for (int i = 2; i < argc; i++) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            printf("Usage: %s edit <name> [-s path]\n", argv[0]);
            printf("\nOpens the secret in $VISUAL / $EDITOR (creates if new).\n");
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

    bge_ustore us = {0};
    int rc = bge_uauthenticate(path, pwd, (size_t)pwd_len, &us);
    cfx_memzero_s(pwd, sizeof(pwd));
    if (rc != 0) return 1;

    uint8_t *pt = NULL;
    size_t pt_len = 0;
    rc = bge_udecrypt(&us, &pt, &pt_len);
    if (rc != 0) {
        bge_ustore_wipe(&us);
        return 1;
    }

    /* interactive name prompt if not on cmdline */
    char name_buf[256] = {0};
    if (!name) {
        store_print_names(pt, pt_len);
        int r = bge_read_visible("Name: ", name_buf, sizeof(name_buf));
        if (r <= 0) {
            fprintf(stderr, "error: name required\n");
            cfx_memzero_s(pt, pt_len);
            free(pt);
            bge_ustore_wipe(&us);
            return 1;
        }
        name = name_buf;
    }

    /* resolve numeric index */
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

    /* get current value (NULL means new entry) */
    size_t cur_vlen = 0;
    const uint8_t *cur_val = store_get(pt, pt_len, name, &cur_vlen);

    /* create temp file and write current value */
    char tmppath[1024];
    if (create_secure_tmpfile(tmppath, sizeof(tmppath)) != 0) {
        fprintf(stderr, "error: cannot create temporary file\n");
        cfx_memzero_s(pt, pt_len);
        free(pt);
        bge_ustore_wipe(&us);
        return 1;
    }

    if (cur_val && cur_vlen > 0) {
        FILE *tf = fopen(tmppath, "wb");
        if (!tf) {
            fprintf(stderr, "error: cannot write temporary file\n");
            wipe_tmpfile(tmppath);
            cfx_memzero_s(pt, pt_len);
            free(pt);
            bge_ustore_wipe(&us);
            return 1;
        }
        fwrite(cur_val, 1, cur_vlen, tf);
        fclose(tf);
    }

    /* launch editor */
    const char *editor = get_editor();
    char cmd_buf[2048];
    int n = snprintf(cmd_buf, sizeof(cmd_buf), "%s \"%s\"", editor, tmppath);
    if (n < 0 || (size_t)n >= sizeof(cmd_buf)) {
        fprintf(stderr, "error: editor command too long\n");
        wipe_tmpfile(tmppath);
        cfx_memzero_s(pt, pt_len);
        free(pt);
        bge_ustore_wipe(&us);
        return 1;
    }

    int edit_rc = system(cmd_buf);

    /* read edited file */
    uint8_t *new_val = NULL;
    size_t new_vlen = 0;
    FILE *rf = fopen(tmppath, "rb");
    if (rf) {
        bge_read_all(rf, &new_val, &new_vlen);
        fclose(rf);
    }

    /* wipe temp file immediately */
    wipe_tmpfile(tmppath);

    if (edit_rc != 0) {
        fprintf(stderr, "error: editor exited with status %d, not saving\n", edit_rc);
        if (new_val) { cfx_memzero_s(new_val, new_vlen); free(new_val); }
        cfx_memzero_s(pt, pt_len);
        free(pt);
        bge_ustore_wipe(&us);
        return 1;
    }

    if (!new_val) {
        /* treat as empty value */
        new_val = malloc(1);
        if (!new_val) {
            cfx_memzero_s(pt, pt_len);
            free(pt);
            bge_ustore_wipe(&us);
            return 1;
        }
        new_vlen = 0;
    }

    /* save */
    size_t out_len;
    uint8_t *new_pt = store_set(pt, pt_len, name, new_val, new_vlen, &out_len);
    cfx_memzero_s(pt, pt_len);
    free(pt);
    cfx_memzero_s(new_val, new_vlen);
    free(new_val);

    if (!new_pt) {
        fprintf(stderr, "error: allocation failed\n");
        bge_ustore_wipe(&us);
        return 1;
    }

    rc = bge_uwrite(path, new_pt, out_len, &us);
    cfx_memzero_s(new_pt, out_len);
    free(new_pt);
    bge_ustore_wipe(&us);
    if (rc == 0) printf("Ok.\n");
    return rc != 0;
}

static int store_cmd_rm(int argc, char **argv) {
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

    bge_ustore us = {0};
    int rc = bge_uauthenticate(path, pwd, (size_t)pwd_len, &us);
    cfx_memzero_s(pwd, sizeof(pwd));
    if (rc != 0) return 1;

    uint8_t *pt = NULL;
    size_t pt_len = 0;
    rc = bge_udecrypt(&us, &pt, &pt_len);
    if (rc != 0) {
        bge_ustore_wipe(&us);
        return 1;
    }

    char name_buf[256] = {0};
    if (!name) {
        int r = bge_read_visible("Name: ", name_buf, sizeof(name_buf));
        if (r <= 0) {
            fprintf(stderr, "error: name required\n");
            cfx_memzero_s(pt, pt_len);
            free(pt);
            bge_ustore_wipe(&us);
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
            bge_ustore_wipe(&us);
            return 1;
        }
    }

    size_t new_len;
    uint8_t *new_pt = store_rm(pt, pt_len, name, &new_len);
    cfx_memzero_s(pt, pt_len);
    free(pt);

    if (!new_pt) {
        bge_ustore_wipe(&us);
        return 1;
    }

    rc = bge_uwrite(path, new_pt, new_len, &us);
    cfx_memzero_s(new_pt, new_len);
    free(new_pt);
    bge_ustore_wipe(&us);
    if (rc == 0) printf("Ok.\n");
    return rc != 0;
}

static int store_cmd_rename(int argc, char **argv) {
    const char *old_name = NULL;
    const char *new_name = NULL;
    const char *path = NULL;
    char path_buf[1024];

    for (int i = 2; i < argc; i++) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            printf("Usage: %s rename <old> <new> [-s path]\n", argv[0]);
            return 0;
        } else if (strcmp(argv[i], "-s") == 0 || strcmp(argv[i], "--store") == 0) {
            if (i + 1 >= argc) { fprintf(stderr, "error: -s requires a path\n"); return 1; }
            path = argv[++i];
        } else if (!old_name) {
            old_name = argv[i];
        } else if (!new_name) {
            new_name = argv[i];
        } else {
            fprintf(stderr, "error: unknown argument: %s\n", argv[i]);
            return 1;
        }
    }

    if (!old_name || !new_name) {
        fprintf(stderr, "error: rename requires <old> and <new> arguments\n");
        return 1;
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

    bge_ustore us = {0};
    int rc = bge_uauthenticate(path, pwd, (size_t)pwd_len, &us);
    cfx_memzero_s(pwd, sizeof(pwd));
    if (rc != 0) return 1;

    uint8_t *pt = NULL;
    size_t pt_len = 0;
    rc = bge_udecrypt(&us, &pt, &pt_len);
    if (rc != 0) {
        bge_ustore_wipe(&us);
        return 1;
    }

    /* resolve numeric index for old_name only */
    char resolved[256];
    if (is_numeric(old_name)) {
        unsigned idx = (unsigned)atoi(old_name);
        uint16_t nlen;
        const uint8_t *nptr = store_name_by_index(pt, pt_len, idx, &nlen);
        if (nptr && nlen < sizeof(resolved)) {
            memcpy(resolved, nptr, nlen);
            resolved[nlen] = '\0';
            old_name = resolved;
        }
    }

    if (strcmp(old_name, new_name) == 0) {
        fprintf(stderr, "error: old and new names are the same\n");
        cfx_memzero_s(pt, pt_len);
        free(pt);
        bge_ustore_wipe(&us);
        return 1;
    }

    /* check old exists */
    size_t vlen;
    const uint8_t *val = store_get(pt, pt_len, old_name, &vlen);
    if (!val) {
        fprintf(stderr, "error: '%s' not found\n", old_name);
        cfx_memzero_s(pt, pt_len);
        free(pt);
        bge_ustore_wipe(&us);
        return 1;
    }

    /* check new does not exist */
    if (store_get(pt, pt_len, new_name, NULL) != NULL) {
        fprintf(stderr, "error: '%s' already exists\n", new_name);
        cfx_memzero_s(pt, pt_len);
        free(pt);
        bge_ustore_wipe(&us);
        return 1;
    }

    /* copy value before store_rm (which allocates a new buffer) */
    uint8_t *val_copy = malloc(vlen);
    if (!val_copy) {
        fprintf(stderr, "error: allocation failed\n");
        cfx_memzero_s(pt, pt_len);
        free(pt);
        bge_ustore_wipe(&us);
        return 1;
    }
    memcpy(val_copy, val, vlen);

    /* remove old entry */
    size_t tmp_len;
    uint8_t *tmp_pt = store_rm(pt, pt_len, old_name, &tmp_len);
    cfx_memzero_s(pt, pt_len);
    free(pt);

    if (!tmp_pt) {
        cfx_memzero_s(val_copy, vlen);
        free(val_copy);
        bge_ustore_wipe(&us);
        return 1;
    }

    /* add new entry */
    size_t new_len;
    uint8_t *new_pt = store_set(tmp_pt, tmp_len, new_name, val_copy, vlen, &new_len);
    cfx_memzero_s(tmp_pt, tmp_len);
    free(tmp_pt);
    cfx_memzero_s(val_copy, vlen);
    free(val_copy);

    if (!new_pt) {
        fprintf(stderr, "error: allocation failed\n");
        bge_ustore_wipe(&us);
        return 1;
    }

    rc = bge_uwrite(path, new_pt, new_len, &us);
    cfx_memzero_s(new_pt, new_len);
    free(new_pt);
    bge_ustore_wipe(&us);
    if (rc == 0) printf("Ok.\n");
    return rc != 0;
}

static int store_cmd_ls(int argc, char **argv) {
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

    bge_ustore us = {0};
    int rc = bge_uauthenticate(path, pwd, (size_t)pwd_len, &us);
    cfx_memzero_s(pwd, sizeof(pwd));
    if (rc != 0) return 1;

    uint8_t *pt = NULL;
    size_t pt_len = 0;
    rc = bge_udecrypt(&us, &pt, &pt_len);
    bge_ustore_wipe(&us);
    if (rc != 0) return 1;

    store_print_names(pt, pt_len);

    cfx_memzero_s(pt, pt_len);
    free(pt);
    return 0;
}

static int store_cmd_info(int argc, char **argv) {
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

    struct stat st;
    if (stat(path, &st) != 0) {
        fprintf(stderr, "error: cannot stat %s: %s\n", path, strerror(errno));
        return 1;
    }

    printf("Store:   %s\n", path);
    printf("Size:    %lld bytes\n", (long long)st.st_size);

    char pwd[256] = {0};
    int pwd_len = bge_read_passphrase("Enter passphrase: ", pwd, sizeof(pwd));
    if (pwd_len <= 0) {
        fprintf(stderr, "error: passphrase required\n");
        cfx_memzero_s(pwd, sizeof(pwd));
        return 1;
    }

    bge_ustore us = {0};
    int rc = bge_uauthenticate(path, pwd, (size_t)pwd_len, &us);
    cfx_memzero_s(pwd, sizeof(pwd));
    if (rc != 0) return 1;

    uint8_t *pt = NULL;
    size_t pt_len = 0;
    rc = bge_udecrypt(&us, &pt, &pt_len);
    if (rc != 0) {
        bge_ustore_wipe(&us);
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

    if (us.version == BGE_V4_VERSION) {
        printf("Version: 4 (multi-password)\n");
        int sc = us.u.v4.hdr.slot_count;
        printf("Slots:   %d\n", sc);
        for (int i = 0; i < sc; i++) {
            const bge_v4_slot *s = &us.u.v4.slots[i];
            printf("  [%d] m=%u KB, t=%u, p=%u%s\n",
                   i + 1, s->m_cost, s->t_cost, s->p_cost,
                   (i == us.u.v4.matched_slot) ? " (authenticated)" : "");
        }
    } else {
        printf("Version: 2\n");
        uint32_t m_cost = cfx_load32_le(&us.u.v2.hdr.m_cost);
        uint32_t t_cost = cfx_load32_le(&us.u.v2.hdr.t_cost);
        uint32_t p_cost = cfx_load32_le(&us.u.v2.hdr.p_cost);
        printf("Argon2:  m=%u KB, t=%u, p=%u\n", m_cost, t_cost, p_cost);
    }

    bge_ustore_wipe(&us);
    return 0;
}

static int store_cmd_passwd(int argc, char **argv) {
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

    char old_pwd[256] = {0};
    int old_len = bge_read_passphrase("Enter current passphrase: ", old_pwd, sizeof(old_pwd));
    if (old_len <= 0) {
        fprintf(stderr, "error: passphrase required\n");
        cfx_memzero_s(old_pwd, sizeof(old_pwd));
        return 1;
    }

    bge_ustore us = {0};
    int rc = bge_uauthenticate(path, old_pwd, (size_t)old_len, &us);
    cfx_memzero_s(old_pwd, sizeof(old_pwd));
    if (rc != 0) return 1;

    uint8_t *pt = NULL;
    size_t pt_len = 0;
    rc = bge_udecrypt(&us, &pt, &pt_len);
    if (rc != 0) {
        bge_ustore_wipe(&us);
        return 1;
    }

    printf("Enter new passphrase.\n");
    char new_pwd[256] = {0};
    int new_len = prompt_passphrase(new_pwd, sizeof(new_pwd));
    if (new_len < 0) {
        cfx_memzero_s(pt, pt_len);
        free(pt);
        bge_ustore_wipe(&us);
        return 1;
    }

    if (us.version == BGE_V4_VERSION) {
        /* v4: re-wrap DEK for matched slot, keep other slots as-is */
        bge_v4_store *s4 = &us.u.v4;
        int mi = s4->matched_slot;

        rc = bge_v4_wrap_dek(new_pwd, (size_t)new_len,
                             BGE_DEFAULT_M, BGE_DEFAULT_T, BGE_DEFAULT_P,
                             s4->dek, &s4->slots[mi]);
        cfx_memzero_s(new_pwd, sizeof(new_pwd));
        if (rc != 0) {
            cfx_memzero_s(pt, pt_len);
            free(pt);
            bge_ustore_wipe(&us);
            return 1;
        }

        rc = bge_v4_encrypt_write(path, pt, pt_len,
                                  &s4->hdr, s4->slots, s4->hdr.slot_count, s4->dek);
    } else {
        /* v2: full re-encrypt with new password */
        rc = bge_encrypt_write(path, pt, pt_len, new_pwd, (size_t)new_len,
                               BGE_DEFAULT_M, BGE_DEFAULT_T, BGE_DEFAULT_P);
        cfx_memzero_s(new_pwd, sizeof(new_pwd));
    }

    cfx_memzero_s(pt, pt_len);
    free(pt);
    bge_ustore_wipe(&us);

    if (rc == 0)
        printf("Passphrase changed.\n");
    return rc != 0;
}

static int store_cmd_dump(int argc, char **argv) {
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

    bge_ustore us = {0};
    int rc = bge_uauthenticate(path, pwd, (size_t)pwd_len, &us);
    cfx_memzero_s(pwd, sizeof(pwd));
    if (rc != 0) return 1;

    uint8_t *pt = NULL;
    size_t pt_len = 0;
    rc = bge_udecrypt(&us, &pt, &pt_len);
    bge_ustore_wipe(&us);
    if (rc != 0) return 1;

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

static int store_slot_ls(int argc, char **argv) {
    const char *path = NULL;
    char path_buf[1024];

    for (int i = 3; i < argc; i++) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            printf("Usage: %s slot ls [-s path]\n", argv[0]);
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

    /* read file without passphrase -- just parse header */
    FILE *f = fopen(path, "rb");
    if (!f) {
        fprintf(stderr, "error: cannot open %s: %s\n", path, strerror(errno));
        return 1;
    }

    uint8_t *file_buf = NULL;
    size_t file_len = 0;
    if (bge_read_all(f, &file_buf, &file_len) != 0) {
        fclose(f);
        fprintf(stderr, "error: cannot read %s\n", path);
        return 1;
    }
    fclose(f);

    if (file_len < 4 || memcmp(file_buf, BGE_MAGIC, 3) != 0) {
        fprintf(stderr, "error: %s is not a valid BGE store\n", path);
        free(file_buf);
        return 1;
    }

    uint8_t version = file_buf[3];

    if (version == BGE_V4_VERSION) {
        if (file_len < BGE_V4_FIXED_HDR_LEN) {
            fprintf(stderr, "error: truncated v4 header\n");
            free(file_buf);
            return 1;
        }
        bge_v4_fixed_header hdr;
        memcpy(&hdr, file_buf, BGE_V4_FIXED_HDR_LEN);
        int sc = hdr.slot_count;
        if (sc < 1 || sc > BGE_V4_MAX_SLOTS) {
            fprintf(stderr, "error: invalid slot count %d\n", sc);
            free(file_buf);
            return 1;
        }

        printf("Version: 4 (multi-password)\n");
        printf("Slots:   %d\n", sc);

        const uint8_t *sp = file_buf + BGE_V4_FIXED_HDR_LEN;
        for (int i = 0; i < sc; i++) {
            if ((size_t)(sp - file_buf) + BGE_V4_SLOT_LEN > file_len) {
                fprintf(stderr, "error: truncated slot %d\n", i);
                break;
            }
            bge_v4_slot slot;
            bge_v4_slot_from_buf(&slot, sp);
            printf("  [%d] m=%u KB, t=%u, p=%u\n",
                   i + 1, slot.m_cost, slot.t_cost, slot.p_cost);
            sp += BGE_V4_SLOT_LEN;
        }
    } else if (version == BGE_VERSION) {
        printf("Version: 2 (single password)\n");
        printf("Slots:   1\n");
        if (file_len >= BGE_HEADER_LEN) {
            bge_header hdr;
            memcpy(&hdr, file_buf, sizeof(hdr));
            uint32_t m = cfx_load32_le(&hdr.m_cost);
            uint32_t t = cfx_load32_le(&hdr.t_cost);
            uint32_t p = cfx_load32_le(&hdr.p_cost);
            printf("  [1] m=%u KB, t=%u, p=%u\n", m, t, p);
        }
    } else {
        printf("Version: %u (unknown)\n", version);
    }

    free(file_buf);
    return 0;
}

static int store_slot_add(int argc, char **argv) {
    const char *path = NULL;
    char path_buf[1024];

    for (int i = 3; i < argc; i++) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            printf("Usage: %s slot add [-s path]\n", argv[0]);
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

    /* authenticate with existing password */
    char pwd[256] = {0};
    int pwd_len = bge_read_passphrase("Enter existing passphrase: ", pwd, sizeof(pwd));
    if (pwd_len <= 0) {
        fprintf(stderr, "error: passphrase required\n");
        cfx_memzero_s(pwd, sizeof(pwd));
        return 1;
    }

    bge_ustore us = {0};
    int rc = bge_uauthenticate(path, pwd, (size_t)pwd_len, &us);
    if (rc != 0) {
        cfx_memzero_s(pwd, sizeof(pwd));
        return 1;
    }

    /* check max slots before doing any more work */
    if (us.version == BGE_V4_VERSION &&
        us.u.v4.hdr.slot_count >= BGE_V4_MAX_SLOTS) {
        fprintf(stderr, "error: maximum %d slots reached\n", BGE_V4_MAX_SLOTS);
        cfx_memzero_s(pwd, sizeof(pwd));
        bge_ustore_wipe(&us);
        return 1;
    }

    /* decrypt store data */
    uint8_t *pt = NULL;
    size_t pt_len = 0;
    rc = bge_udecrypt(&us, &pt, &pt_len);
    if (rc != 0) {
        cfx_memzero_s(pwd, sizeof(pwd));
        bge_ustore_wipe(&us);
        return 1;
    }

    /* prompt for new passphrase */
    printf("Enter new passphrase for additional slot.\n");
    char new_pwd[256] = {0};
    int new_len = prompt_passphrase(new_pwd, sizeof(new_pwd));
    if (new_len < 0) {
        cfx_memzero_s(pwd, sizeof(pwd));
        cfx_memzero_s(pt, pt_len);
        free(pt);
        bge_ustore_wipe(&us);
        return 1;
    }

    if (us.version == BGE_V4_VERSION) {
        /* already v4: append a slot */
        bge_v4_store *s4 = &us.u.v4;
        int sc = s4->hdr.slot_count;

        /* wrap DEK with new password */
        rc = bge_v4_wrap_dek(new_pwd, (size_t)new_len,
                             BGE_DEFAULT_M, BGE_DEFAULT_T, BGE_DEFAULT_P,
                             s4->dek, &s4->slots[sc]);
        cfx_memzero_s(new_pwd, sizeof(new_pwd));
        cfx_memzero_s(pwd, sizeof(pwd));
        if (rc != 0) {
            cfx_memzero_s(pt, pt_len);
            free(pt);
            bge_ustore_wipe(&us);
            return 1;
        }

        s4->hdr.slot_count = (uint8_t)(sc + 1);
        rc = bge_v4_encrypt_write(path, pt, pt_len,
                                  &s4->hdr, s4->slots, sc + 1, s4->dek);
    } else {
        /* v2 -> v4 migration */
        bge_v4_fixed_header hdr;
        memcpy(hdr.magic, BGE_MAGIC, 3);
        hdr.version = BGE_V4_VERSION;
        hdr.slot_count = 2;
        memset(hdr.reserved, 0, sizeof(hdr.reserved));

        /* generate random DEK */
        uint8_t dek[BGE_V4_DEK_LEN];
        cfx_srand_os();
        cfx_rand_bytes(dek, sizeof(dek));

        bge_v4_slot slots[2];
        memset(slots, 0, sizeof(slots));

        /* slot 0: wrap DEK under existing password */
        rc = bge_v4_wrap_dek(pwd, (size_t)pwd_len,
                             BGE_DEFAULT_M, BGE_DEFAULT_T, BGE_DEFAULT_P,
                             dek, &slots[0]);
        cfx_memzero_s(pwd, sizeof(pwd));
        if (rc != 0) {
            cfx_memzero_s(new_pwd, sizeof(new_pwd));
            cfx_memzero_s(dek, sizeof(dek));
            cfx_memzero_s(pt, pt_len);
            free(pt);
            bge_ustore_wipe(&us);
            return 1;
        }

        /* slot 1: wrap DEK under new password */
        rc = bge_v4_wrap_dek(new_pwd, (size_t)new_len,
                             BGE_DEFAULT_M, BGE_DEFAULT_T, BGE_DEFAULT_P,
                             dek, &slots[1]);
        cfx_memzero_s(new_pwd, sizeof(new_pwd));
        if (rc != 0) {
            cfx_memzero_s(dek, sizeof(dek));
            cfx_memzero_s(pt, pt_len);
            free(pt);
            bge_ustore_wipe(&us);
            return 1;
        }

        rc = bge_v4_encrypt_write(path, pt, pt_len, &hdr, slots, 2, dek);
        cfx_memzero_s(dek, sizeof(dek));
        cfx_memzero_s(slots, sizeof(slots));
    }

    int was_v4 = (us.version == BGE_V4_VERSION);

    cfx_memzero_s(pt, pt_len);
    free(pt);
    bge_ustore_wipe(&us);

    if (rc == 0) {
        printf("Slot added.\n");
        if (!was_v4)
            printf("Store migrated from v2 to v4.\n");
    }
    return rc != 0;
}

static int store_slot_rm(int argc, char **argv) {
    const char *path = NULL;
    char path_buf[1024];
    int target_slot = -1;  /* -1 means: use matched slot */

    for (int i = 3; i < argc; i++) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            printf("Usage: %s slot rm [N] [-s path]\n", argv[0]);
            printf("  N is 1-based slot index (default: authenticated slot)\n");
            return 0;
        } else if (strcmp(argv[i], "-s") == 0 || strcmp(argv[i], "--store") == 0) {
            if (i + 1 >= argc) { fprintf(stderr, "error: -s requires a path\n"); return 1; }
            path = argv[++i];
        } else if (target_slot < 0 && is_numeric(argv[i])) {
            target_slot = atoi(argv[i]);
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

    bge_ustore us = {0};
    int rc = bge_uauthenticate(path, pwd, (size_t)pwd_len, &us);
    cfx_memzero_s(pwd, sizeof(pwd));
    if (rc != 0) return 1;

    if (us.version != BGE_V4_VERSION) {
        fprintf(stderr, "error: slot rm requires a v4 store (use 'slot add' first)\n");
        bge_ustore_wipe(&us);
        return 1;
    }

    bge_v4_store *s4 = &us.u.v4;
    int sc = s4->hdr.slot_count;

    if (sc <= 1) {
        fprintf(stderr, "error: cannot remove the only slot\n");
        bge_ustore_wipe(&us);
        return 1;
    }

    /* resolve target slot (convert 1-based to 0-based) */
    int rm_idx;
    if (target_slot > 0) {
        rm_idx = target_slot - 1;
        if (rm_idx >= sc) {
            fprintf(stderr, "error: slot %d does not exist (store has %d slots)\n",
                    target_slot, sc);
            bge_ustore_wipe(&us);
            return 1;
        }
    } else {
        rm_idx = s4->matched_slot;
    }

    /* confirm */
    char prompt[128];
    snprintf(prompt, sizeof(prompt), "Remove slot %d? [y/N] ", rm_idx + 1);
    char ans[8] = {0};
    int r = bge_read_visible(prompt, ans, sizeof(ans));
    if (r <= 0 || (ans[0] != 'y' && ans[0] != 'Y')) {
        fprintf(stderr, "Aborted.\n");
        bge_ustore_wipe(&us);
        return 1;
    }

    /* decrypt store data */
    uint8_t *pt = NULL;
    size_t pt_len = 0;
    rc = bge_udecrypt(&us, &pt, &pt_len);
    if (rc != 0) {
        bge_ustore_wipe(&us);
        return 1;
    }

    /* remove slot from array */
    for (int i = rm_idx; i < sc - 1; i++) {
        s4->slots[i] = s4->slots[i + 1];
    }
    memset(&s4->slots[sc - 1], 0, sizeof(bge_v4_slot));
    s4->hdr.slot_count = (uint8_t)(sc - 1);

    /* re-encrypt */
    rc = bge_v4_encrypt_write(path, pt, pt_len,
                              &s4->hdr, s4->slots, sc - 1, s4->dek);
    cfx_memzero_s(pt, pt_len);
    free(pt);
    bge_ustore_wipe(&us);

    if (rc == 0)
        printf("Slot %d removed.\n", rm_idx + 1);
    return rc != 0;
}

static int store_slot(int argc, char **argv) {
    if (argc < 3) {
        printf("Usage: %s slot <add|rm|ls> [-s path]\n", argv[0]);
        return 1;
    }

    const char *sub = argv[2];
    if (strcmp(sub, "ls")   == 0 || strcmp(sub, "list") == 0) return store_slot_ls(argc, argv);
    if (strcmp(sub, "add")  == 0) return store_slot_add(argc, argv);
    if (strcmp(sub, "rm")   == 0 || strcmp(sub, "remove") == 0) return store_slot_rm(argc, argv);

    fprintf(stderr, "Unknown slot subcommand: %s\n", sub);
    printf("Usage: %s slot <add|rm|ls> [-s path]\n", argv[0]);
    return 1;
}

static int store_cmd_rekey(int argc, char **argv) {
    const char *path = NULL;
    char path_buf[1024];

    for (int i = 2; i < argc; i++) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            printf("Usage: %s rekey [-s path]\n", argv[0]);
            printf("Generate a new DEK and re-wrap all slots. V4 only.\n");
            printf("You will be prompted for each slot's passphrase.\n");
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

    /* first, authenticate to decrypt the store */
    char pwd[256] = {0};
    int pwd_len = bge_read_passphrase("Enter any passphrase to unlock: ", pwd, sizeof(pwd));
    if (pwd_len <= 0) {
        fprintf(stderr, "error: passphrase required\n");
        cfx_memzero_s(pwd, sizeof(pwd));
        return 1;
    }

    bge_ustore us = {0};
    int rc = bge_uauthenticate(path, pwd, (size_t)pwd_len, &us);
    cfx_memzero_s(pwd, sizeof(pwd));
    if (rc != 0) return 1;

    if (us.version != BGE_V4_VERSION) {
        fprintf(stderr, "error: rekey requires a v4 store (use 'slot add' first)\n");
        bge_ustore_wipe(&us);
        return 1;
    }

    uint8_t *pt = NULL;
    size_t pt_len = 0;
    rc = bge_udecrypt(&us, &pt, &pt_len);
    if (rc != 0) {
        bge_ustore_wipe(&us);
        return 1;
    }

    bge_v4_store *s4 = &us.u.v4;
    int sc = s4->hdr.slot_count;

    /* generate new DEK */
    uint8_t new_dek[BGE_V4_DEK_LEN];
    cfx_srand_os();
    cfx_rand_bytes(new_dek, sizeof(new_dek));

    /* for each slot, prompt for that slot's passphrase, verify, re-wrap */
    bge_v4_slot new_slots[BGE_V4_MAX_SLOTS];
    memset(new_slots, 0, sizeof(new_slots));

    for (int i = 0; i < sc; i++) {
        char prompt[128];
        snprintf(prompt, sizeof(prompt), "Enter passphrase for slot %d: ", i + 1);
        char slot_pwd[256] = {0};
        int slot_pwd_len = bge_read_secret(prompt, slot_pwd, sizeof(slot_pwd));
        if (slot_pwd_len <= 0) {
            fprintf(stderr, "error: passphrase required for slot %d\n", i + 1);
            cfx_memzero_s(slot_pwd, sizeof(slot_pwd));
            cfx_memzero_s(new_dek, sizeof(new_dek));
            cfx_memzero_s(pt, pt_len);
            free(pt);
            bge_ustore_wipe(&us);
            return 1;
        }

        /* verify against slot's existing verifier */
        uint8_t kdf_out[48];
        rc = cfx_argon2id(kdf_out, 48,
            (const uint8_t *)slot_pwd, (size_t)slot_pwd_len,
            s4->slots[i].salt, sizeof(s4->slots[i].salt),
            s4->slots[i].m_cost, s4->slots[i].t_cost, s4->slots[i].p_cost);
        if (rc != 0) {
            fprintf(stderr, "error: argon2 failed for slot %d\n", i + 1);
            cfx_memzero_s(kdf_out, sizeof(kdf_out));
            cfx_memzero_s(slot_pwd, sizeof(slot_pwd));
            cfx_memzero_s(new_dek, sizeof(new_dek));
            cfx_memzero_s(pt, pt_len);
            free(pt);
            bge_ustore_wipe(&us);
            return 1;
        }

        uint8_t diff = 0;
        for (int j = 0; j < 16; j++)
            diff |= kdf_out[32 + j] ^ s4->slots[i].verifier[j];

        cfx_memzero_s(kdf_out, sizeof(kdf_out));

        if (diff != 0) {
            fprintf(stderr, "error: wrong passphrase for slot %d\n", i + 1);
            cfx_memzero_s(slot_pwd, sizeof(slot_pwd));
            cfx_memzero_s(new_dek, sizeof(new_dek));
            cfx_memzero_s(pt, pt_len);
            free(pt);
            bge_ustore_wipe(&us);
            return 1;
        }

        /* wrap new DEK with this slot's password */
        rc = bge_v4_wrap_dek(slot_pwd, (size_t)slot_pwd_len,
                             s4->slots[i].m_cost,
                             s4->slots[i].t_cost,
                             s4->slots[i].p_cost,
                             new_dek, &new_slots[i]);
        cfx_memzero_s(slot_pwd, sizeof(slot_pwd));
        if (rc != 0) {
            cfx_memzero_s(new_dek, sizeof(new_dek));
            cfx_memzero_s(pt, pt_len);
            free(pt);
            bge_ustore_wipe(&us);
            return 1;
        }
    }

    /* re-encrypt store with new DEK */
    rc = bge_v4_encrypt_write(path, pt, pt_len, &s4->hdr, new_slots, sc, new_dek);

    cfx_memzero_s(new_dek, sizeof(new_dek));
    cfx_memzero_s(new_slots, sizeof(new_slots));
    cfx_memzero_s(pt, pt_len);
    free(pt);
    bge_ustore_wipe(&us);

    if (rc == 0)
        printf("Rekey complete. All %d slots re-wrapped with new data key.\n", sc);
    return rc != 0;
}

int cfx_store_run(int argc, char **argv) {
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

    if (strcmp(cmd, "init")   == 0) return store_init(argc, argv);
    if (strcmp(cmd, "get")    == 0) return store_cmd_get(argc, argv);
    if (strcmp(cmd, "set")    == 0) return store_cmd_set(argc, argv);
    if (strcmp(cmd, "edit")   == 0) return store_cmd_edit(argc, argv);
    if (strcmp(cmd, "rm")     == 0) return store_cmd_rm(argc, argv);
    if (strcmp(cmd, "rename") == 0 ||
        strcmp(cmd, "mv")    == 0) return store_cmd_rename(argc, argv);
    if (strcmp(cmd, "ls")     == 0 ||
        strcmp(cmd, "list")   == 0) return store_cmd_ls(argc, argv);
    if (strcmp(cmd, "info")   == 0) return store_cmd_info(argc, argv);
    if (strcmp(cmd, "passwd") == 0) return store_cmd_passwd(argc, argv);
    if (strcmp(cmd, "dump")   == 0) return store_cmd_dump(argc, argv);
    if (strcmp(cmd, "slot")   == 0) return store_slot(argc, argv);
    if (strcmp(cmd, "rekey")  == 0) return store_cmd_rekey(argc, argv);

    fprintf(stderr, "Unknown command: %s\n", cmd);
    usage(argv[0]);
    return 1;
}
