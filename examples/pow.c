#include "cfx/big.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

static void usage(const char* prog) {
    fprintf(stderr, "Usage: %s [-x|-d] [-ix|-id] n p\n", prog);
    fprintf(stderr, "    calculate n^p\n");
    fprintf(stderr, "   -x : print in hex\n");
    fprintf(stderr, "   -b : print in binary\n");
    fprintf(stderr, "   -d : print in decimal (default)\n");
    fprintf(stderr, "   -ix : interpret input as hex\n");
    fprintf(stderr, "   -id : interpret input as decimal\n");
}



int main(int argc, char* argv[]) {
    int print_hex = 0;
    int print_bin = 0;

    if (argc < 2) {
        usage(argv[0]);
        return 1;
    }

    /* Parse flags */
    int argi = 1;
    int base = 10; /* default */

    if (strcmp(argv[argi], "-x") == 0) {
        print_hex = 1; argi++;
    } else if (strcmp(argv[argi], "-d") == 0) {
        print_hex = 0; argi++;
    } else if (strcmp(argv[argi], "-b") == 0) {
        print_bin = 1; argi++;
    } 
    if (strcmp(argv[argi], "-ix") == 0) {
        base = 16; argi++;
    } else if (strcmp(argv[argi], "-id") == 0) {
        base = 10; argi++;
    } 

    if (argi >= argc - 1) {
        usage(argv[0]);
        return 1;
    }
    char* nstr = argv[argi++];
    char* pstr = argv[argi];

    cfx_big_t n, p;
    cfx_big_init(&n);
    cfx_big_init(&p);
    if (base == 16) {
        cfx_big_from_hex(&n, nstr);
        cfx_big_from_hex(&p, pstr);
    } else if (base == 10) {
        cfx_big_from_str(&n, nstr);
        cfx_big_from_str(&p, pstr);
    }

    cfx_big_t np;
    cfx_big_init(&np);
    cfx_big_pow(&np, &n, &p);

    char* npstr;
    size_t sz;
    if (print_hex) {
        npstr = cfx_big_to_hex(&np, &sz);
    } else if (print_bin) {
        npstr = cfx_big_to_bin(&np, &sz);
    } else {
        npstr = cfx_big_to_str(&np, &sz);
    }

    printf("n: %s\np: %s\nn^p: %s\ndigits: %zu\n", nstr, pstr, npstr, sz);

    cfx_big_free(&n);
    cfx_big_free(&p);
    // cfx_big_free(&acc);
    cfx_big_free(&np);
    free(npstr);
    
    return 0;
}
