#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "cfx/big.h"
#include "cfx_utils.h"

static void usage(const char* prog) {
    fprintf(stderr, "Usage: %s [-x|-d] [-ix|-id|ib] file\n", prog);
    fprintf(stderr, "   -x : print in hex\n");
    fprintf(stderr, "   -b : print in binary\n");
    fprintf(stderr, "   -d : print in decimal (default)\n");
    fprintf(stderr, "   -ix : interpret input as hex\n");
    fprintf(stderr, "   -id : interpret input as decimal\n");
    fprintf(stderr, "   -ib : interpret input as binary\n");
}

int cfx_fread_run(int argc, char* argv[]) {
    int cfx_print_hex = 0;
    int print_bin = 0;
    const char* fname = NULL;

    if (argc < 2 || (argc == 2 && (strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0))) {
        usage(argv[0]);
        return argc < 2 ? 1 : 0;
    }


    /* Parse flags */
    int argi = 1;
    int base = 10; /* default */

    if (strcmp(argv[argi], "-x") == 0) {
        cfx_print_hex = 1; argi++;
    } else if (strcmp(argv[argi], "-d") == 0) {
        cfx_print_hex = 0; argi++;
    } else if (strcmp(argv[argi], "-b") == 0) {
        print_bin = 1; argi++;
    }
    if (strcmp(argv[argi], "-ix") == 0) {
        base = 16; argi++;
    } else if (strcmp(argv[argi], "-id") == 0) {
        base = 10; argi++;
    } else if (strcmp(argv[argi], "-ib") == 0) {
        base = 2; argi++;
    }

    if (argi >= argc) {
        usage(argv[0]);
        return 1;
    }

    fname = argv[argi];

    FILE* fp = fopen(fname, "rb");
    if (!fp) {
        perror("fopen");
        return 1;
    }

    cfx_big_t big;
    cfx_big_init(&big);

    if (cfx_big_from_file(&big, fp, base) != 0) {
        perror("cfx_big_from_file");
        fclose(fp);
        return 1;
    }
    fclose(fp);

    size_t sz = 0;
    char* s = NULL;

    if (cfx_print_hex) {
        s = cfx_big_to_hex(&big, &sz);
    } else if (print_bin) {
        s = cfx_big_to_bin(&big, &sz);
    } else {
        s = cfx_big_to_str(&big, &sz);
    }

    if (!s) {
        fprintf(stderr, "conversion failed\n");
        cfx_big_free(&big);
        return 1;
    }

    printf("%s\n", s);

    printf("digits: %zu, limbs: %zu\n", strlen(s), big.n);
    free(s);
    cfx_big_free(&big);
    return 0;
}
