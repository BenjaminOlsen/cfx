/* cfx_bge.c -- file encryption dispatch (XChaCha20-Poly1305) */

#include "cfx_bge_internal.h"

static void usage(const char *prog) {
    printf("Usage:\n");
    printf("  %s -e [-a] [-i INPUT] [-o OUTPUT]    Encrypt a file\n", prog);
    printf("  %s -d [-i INPUT] [-o OUTPUT]         Decrypt a file\n\n", prog);
    printf("File encryption using Argon2id key derivation and\n");
    printf("XChaCha20-Poly1305 AEAD.\n\n");
    printf("Options:\n");
    printf("  -e, --encrypt    Encrypt input to output\n");
    printf("  -d, --decrypt    Decrypt input to output\n");
    printf("  -i, --input      Read from file (default: stdin)\n");
    printf("  -o, --output     Write to file (default: stdout)\n");
    printf("  -a, --armor      PEM-encoded base64 (encrypt only; decrypt auto-detects)\n");
    printf("  -p, --passphrase <pw>  Supply passphrase on command line (default: prompt)\n");
    printf("\nExamples:\n");
    printf("  %s -e -i doc.txt -o secret.bge        Encrypt a file\n", prog);
    printf("  %s -e -a -i doc.txt -o secret.asc     Encrypt to armored output\n", prog);
    printf("  %s -d -i secret.bge -o doc.txt        Decrypt a file\n", prog);
    printf("  echo data | %s -e | %s -d             Roundtrip via pipe\n", prog, prog);
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

    /* scan for -e/-d anywhere in argv */
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-e") == 0 || strcmp(argv[i], "--encrypt") == 0) {
            return bge_encrypt_file(argc, argv);
        }
        if (strcmp(argv[i], "-d") == 0 || strcmp(argv[i], "--decrypt") == 0) {
            return bge_decrypt_file(argc, argv);
        }
    }

    fprintf(stderr, "Unknown command: %s\n", cmd);
    usage(argv[0]);
    return 1;
}
