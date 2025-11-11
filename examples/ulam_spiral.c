#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cfx/algo.h"

static void usage(const char* prog) {
    fprintf(stderr,
        "Usage: %s [-c] [-o <outputfilename>] [-s <scale>] size \n"
        "  Creates a PBM (P1 ASCII) image of an Ulam spiral of width & height 'size'.\n"
        "  -c (optional flag) console output\n"
        "  -o (optional) specify filename\n"
        "  -s scale (optional) >=1 scales each cell to SxS pixels (default 1)\n",
        prog);
}

enum direction { UP, DOWN, LEFT, RIGHT };

void turn(enum direction* d) {
    switch(*d) {
        case UP: *d = LEFT; break;
        case DOWN: *d = RIGHT; break;
        case LEFT: *d = DOWN; break;
        case RIGHT: *d = UP; break;
        default: printf("turn problem\n"); exit(-11);
    }
}

int main(int argc, char** argv) {
    if (argc < 2) {
        usage(argv[0]);
        return EXIT_FAILURE;
    }
    int console = 0;
    const char* outfile = "ulam.pbm";
    int scale = 1;
    int w = 0;

    /* cmd line */
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-c") == 0) {
            console = 1;
        } else if (strcmp(argv[i], "-o") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "Error: -o requires a filename\n");
                usage(argv[0]);
                return EXIT_FAILURE;
            }
            outfile = argv[++i];
        } else if (strcmp(argv[i], "-s") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "Error: -s requires a scale value\n");
                usage(argv[0]);
                return EXIT_FAILURE;
            }
            scale = atoi(argv[++i]);
            if (scale <= 0) scale = 1;
        } else if (argv[i][0] == '-') {
            fprintf(stderr, "Unknown option: %s\n", argv[i]);
            usage(argv[0]);
            return EXIT_FAILURE;
        } else {

            w = atoi(argv[i]);
        }
    }

    printf("writing a %d x %d image to %s, scale %d\n", w, w, outfile, scale);
    
    if (w <= 0 || (w % 2) == 0) {
        fprintf(stderr, "Size must be a positive odd integer.\n");
        return EXIT_FAILURE;
    }
    
    int* values = (int*)calloc((size_t)w*w, sizeof(*values));
    unsigned char* prime = (unsigned char*)calloc((size_t)w*w, 1);

    /*
    "go step times, turn, go step times, turn, increase step by one, repeat"
    dir = right;
    turn: right->up->left->down->right (D4 :P)
    */

    int row = (int)(w/2);
    int col = row;
    int curstep = 1;
    enum direction d = RIGHT;
    int sides_before_inc = 2;  /* every 2 sides of the spiral we increment the segments per side */
    int segs_left = curstep;

    for (int k = 1; k <= w*w; ++k) {
        values[col+w*row] = k;
        prime[col+w*row] = cfx_is_prime_u64(k) ? (unsigned char)1 : (unsigned char)0;
        // printf("values[%d][%d] = %d%c\n", col, row, k,
        //     cfx_is_prime_u64(k) ? '*': ' ');

        switch (d) {
            case UP: --row; break;
            case DOWN: ++row; break;
            case LEFT: --col; break;
            case RIGHT: ++col; break;
            default: exit(-1);
        }
        --segs_left;
        if (!segs_left) {
            turn(&d);
            // printf("turn!\n");
            segs_left = curstep;
            --sides_before_inc;
        }
        if (!sides_before_inc) {
            curstep++;
            // printf("curstep: %d\n", curstep);
            segs_left++;
            sides_before_inc = 2;
        }
    }

    // console out:
    if (console) {
        for (int row = 0; row < w; ++row) {
            for (int col = 0; col < w; ++col) {
                cfx_limb_t val = (cfx_limb_t)values[col+w*row];
                if (val == 1) printf("O");
                else if (cfx_is_prime_u64(val)) printf("%3d ", values[col+w*row]);
                else printf("%c ", ' ');
            }
            printf("\n");
        }
    }

    // file out: PBM (P1) ascii old school
    FILE* f = fopen(outfile, "w");
    if (!f) {
        perror("fopen");
        free(prime);
        return EXIT_FAILURE;
    }

    int W = w * scale;
    int H = w * scale;
    fprintf(f, "P1\n%d %d\n", W, H);

    /* make the .pbm pixels*/
    for (int row = 0; row < w; ++row) {
        for (int vrep = 0; vrep < scale; ++vrep) {  /* vertical */
            for (int col = 0; col < w; ++col) {
                int bit = prime[(size_t)row * (size_t)w + (size_t)col] ? 1 : 0;
                for (int hrep = 0; hrep < scale; ++hrep) { /* horizontal */
                    fputc(bit ? '1' : '0', f);
                    fputc(' ', f);
                }
            }
            fputc('\n', f);
        }
    }

    fclose(f);
    free(values);
    free(prime);
    return EXIT_SUCCESS;
}
