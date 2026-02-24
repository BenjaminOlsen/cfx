/* cfx_argon2.c - Password hashing with Argon2id */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cfx/argon2.h"
#include "cfx/rand.h"
#include "cfx/memory.h"
#include "misc.h"

static void print_help(const char* prog) {
    printf("Usage: %s [options] [password]\n", prog);
    printf("       %s -v <encoded> [password]\n\n", prog);
    printf("Hash passwords with Argon2id (RFC 9106)\n\n");
    printf("Options:\n");
    printf("  -m <KB>      Memory cost in KB (default: 65536 = 64MB)\n");
    printf("  -t <N>       Time cost / iterations (default: 3)\n");
    printf("  -p <N>       Parallelism / lanes (default: 4)\n");
    printf("  -v <encoded> Verify password against encoded hash\n");
    printf("  -h, --help   Show this help\n\n");
    printf("If password is omitted, reads from stdin (recommended).\n\n");
    printf("Examples:\n");
    printf("  echo -n 'password' | %s          Read from stdin (secure)\n", prog);
    printf("  %s -m 32768 -t 2 < pw.txt        Use 32MB, 2 iterations\n", prog);
    printf("  %s -v '$argon2id$...' < pw.txt   Verify password\n", prog);
    printf("\n");
    printf("WARNING: Passing password on command line exposes it in process\n");
    printf("listings and shell history. Prefer stdin or file redirection.\n");
}

static char* read_password_stdin(void) {
    char* line = cfx_read_line_stdin();
    if (!line) {
        fprintf(stderr, "Error reading password from stdin\n");
        return NULL;
    }
    return line;
}

int cfx_argon2_run(int argc, char** argv) {
    uint32_t m_cost = 65536;  /* 64 MB */
    uint32_t t_cost = 3;
    uint32_t p_cost = 4;
    const char* verify_encoded = NULL;
    const char* password = NULL;

    /* parse arguments */
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            print_help(argv[0]);
            return 0;
        } else if (strcmp(argv[i], "-m") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "Error: -m requires an argument\n");
                return 1;
            }
            m_cost = (uint32_t)atoi(argv[++i]);
            if (m_cost < 8) {
                fprintf(stderr, "Error: memory cost must be at least 8 KB\n");
                return 1;
            }
        } else if (strcmp(argv[i], "-t") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "Error: -t requires an argument\n");
                return 1;
            }
            t_cost = (uint32_t)atoi(argv[++i]);
            if (t_cost < 1) {
                fprintf(stderr, "Error: time cost must be at least 1\n");
                return 1;
            }
        } else if (strcmp(argv[i], "-p") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "Error: -p requires an argument\n");
                return 1;
            }
            p_cost = (uint32_t)atoi(argv[++i]);
            if (p_cost < 1) {
                fprintf(stderr, "Error: parallelism must be at least 1\n");
                return 1;
            }
        } else if (strcmp(argv[i], "-v") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "Error: -v requires an encoded hash argument\n");
                return 1;
            }
            verify_encoded = argv[++i];
        } else if (argv[i][0] == '-') {
            fprintf(stderr, "Unknown option: %s\n", argv[i]);
            return 1;
        } else {
            if (password == NULL) {
                password = argv[i];
            } else {
                fprintf(stderr, "Error: too many arguments\n");
                return 1;
            }
        }
    }

    /* read password from stdin if not provided */
    char* stdin_password = NULL;
    if (password == NULL) {
        stdin_password = read_password_stdin();
        if (!stdin_password) {
            return 1;
        }
        password = stdin_password;
    }

    int result;

    if (verify_encoded) {
        /* verify mode */
        result = cfx_argon2_verify(verify_encoded,
                                   (const uint8_t*)password, strlen(password));
        if (result == 0) {
            printf("OK\n");
            result = 0;
        } else {
            printf("FAIL\n");
            result = 1;
        }
    } else {
        /* hash mode */
        uint8_t salt[16];
        cfx_rand_bytes(salt, sizeof(salt));

        uint8_t hash[32];
        int rc = cfx_argon2id(hash, sizeof(hash),
                              (const uint8_t*)password, strlen(password),
                              salt, sizeof(salt),
                              m_cost, t_cost, p_cost);

        if (rc != 0) {
            fprintf(stderr, "Error: argon2 failed\n");
            result = 1;
        } else {
            char encoded[256];
            int len = cfx_argon2_encode(encoded, sizeof(encoded),
                                        hash, sizeof(hash),
                                        salt, sizeof(salt),
                                        m_cost, t_cost, p_cost,
                                        CFX_ARGON2ID);
            if (len < 0) {
                fprintf(stderr, "Error: encoding failed\n");
                result = 1;
            } else {
                printf("%s\n", encoded);
                result = 0;
            }

            cfx_memzero_s(hash, sizeof(hash));
        }

        cfx_memzero_s(salt, sizeof(salt));
    }

    if (stdin_password) {
        cfx_memzero_s(stdin_password, strlen(stdin_password));
        free(stdin_password);
    }

    return result;
}
