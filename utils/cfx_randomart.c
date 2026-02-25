/* cfx_randomart.c - drunken bishop ASCII art
   ref: https://www.dirk-loss.de/sshvis/drunken_bishop.pdf */

#include "cfx_randomart.h"
#include "cfx/sha256.h"

#include <stdio.h>
#include <string.h>

static const char glyphs[] = " .o+=*BOX@%&#/^SE";

void cfx_randomart(char *out, int w, int h,
                   const uint8_t *data, size_t len) {
    /* zero visit counts */
    uint8_t grid[CFX_RA_MAX_HEIGHT][CFX_RA_MAX_WIDTH];
    memset(grid, 0, sizeof(grid));

    int x = w / 2;
    int y = h / 2;
    int sx = x, sy = y;

    /* walk */
    for (size_t i = 0; i < len; i++) {
        uint8_t b = data[i];
        for (int j = 0; j < 4; j++) {
            int dx = (b & 1) ? 1 : -1;
            int dy = (b & 2) ? 1 : -1;
            b >>= 2;

            x += dx;
            y += dy;
            if (x < 0) x = 0;
            if (x >= w) x = w - 1;
            if (y < 0) y = 0;
            if (y >= h) y = h - 1;

            if (grid[y][x] < 14) grid[y][x]++;
        }
    }

    /* render to output buffer: out[row * (w+1) + col] */
    for (int row = 0; row < h; row++) {
        for (int col = 0; col < w; col++) {
            char c;
            if (row == sy && col == sx)
                c = 'S';
            else if (row == y && col == x)
                c = 'E';
            else
                c = glyphs[grid[row][col]];
            out[row * (w + 1) + col] = c;
        }
        out[row * (w + 1) + w] = '\0';
    }
}

void cfx_print_randomart_frame(const char *art, int w, int h,
                               const char *label) {
    /* top border: +---[label]---+ */
    int llen = label ? (int)strlen(label) : 0;
    int pad = w - 2 - llen;  /* space for dashes around label */
    int left = 2;  /* at least "+-" before "[" */
    int right = pad - left;
    if (right < 1) { left = pad - 1; right = 1; }

    printf("+");
    for (int i = 0; i < left; i++) printf("-");
    if (label) printf("[%s]", label);
    for (int i = 0; i < right; i++) printf("-");
    printf("+\n");

    for (int row = 0; row < h; row++) {
        printf("|%s|\n", art + row * (w + 1));
    }

    printf("+");
    for (int i = 0; i < w; i++) printf("-");
    printf("+\n");
}

void cfx_print_randomart(const uint8_t *data, size_t len,
                         const char *label) {
    uint8_t hash[32];
    cfx_sha256(hash, data, len);

    int w = CFX_RA_DEFAULT_WIDTH;
    int h = CFX_RA_DEFAULT_HEIGHT;
    char art[CFX_RA_DEFAULT_HEIGHT][CFX_RA_DEFAULT_WIDTH + 1];
    cfx_randomart((char *)art, w, h, hash, 32);
    cfx_print_randomart_frame((const char *)art, w, h, label);
}
