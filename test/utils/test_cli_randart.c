/* test_randart.c - tests for drunken bishop randomart */

#include "test_common.h"
#include "../../utils/cfx_randomart.h"

#define W CFX_RA_DEFAULT_WIDTH
#define H CFX_RA_DEFAULT_HEIGHT

/* count occurrences of char c in art buffer */
static int count_char(const char *art, int w, int h, char c) {
    int n = 0;
    for (int r = 0; r < h; r++)
        for (int col = 0; col < w; col++)
            if (art[r * (w + 1) + col] == c) n++;
    return n;
}

/* all zeros -> bishop goes NW repeatedly, pins at (0,0) */
static void test_all_zeros(void) {
    uint8_t data[32];
    memset(data, 0x00, sizeof(data));

    char art[H][W + 1];
    cfx_randomart((char *)art, W, H, data, 32);

    CFX_ASSERT(art[0][0] == 'E' || art[0][0] == 'S');
    CFX_ASSERT(count_char((char *)art, W, H, 'S') == 1);
    CFX_ASSERT(count_char((char *)art, W, H, 'E') == 1);
}

/* all 0xFF -> bishop goes SE repeatedly, pins at bottom-right */
static void test_all_ff(void) {
    uint8_t data[32];
    memset(data, 0xFF, sizeof(data));

    char art[H][W + 1];
    cfx_randomart((char *)art, W, H, data, 32);

    CFX_ASSERT(art[H - 1][W - 1] == 'E');
}

/* start and end markers present exactly once */
static void test_markers(void) {
    uint8_t data[] = { 0xab, 0xcd, 0xef, 0x01, 0x23, 0x45, 0x67, 0x89 };
    char art[H][W + 1];
    cfx_randomart((char *)art, W, H, data, sizeof(data));

    CFX_ASSERT(count_char((char *)art, W, H, 'S') == 1);
    CFX_ASSERT(count_char((char *)art, W, H, 'E') == 1);
}

/* empty input: S and E at center, same cell */
static void test_empty(void) {
    char art[H][W + 1];
    cfx_randomart((char *)art, W, H, NULL, 0);

    int cx = W / 2, cy = H / 2;
    char c = art[cy][cx];
    CFX_ASSERT(c == 'S' || c == 'E');
}

/* custom grid dimensions */
static void test_custom_dims(void) {
    int cw = 9, ch = 5;
    uint8_t data[] = { 0xde, 0xad, 0xbe, 0xef };
    char art[5][10];
    cfx_randomart((char *)art, cw, ch, data, sizeof(data));

    for (int r = 0; r < ch; r++) {
        CFX_ASSERT(art[r][cw] == '\0');
        CFX_ASSERT((int)strlen(art[r]) == cw);
    }
    CFX_ASSERT(count_char((char *)art, cw, ch, 'S') == 1);
    CFX_ASSERT(count_char((char *)art, cw, ch, 'E') == 1);
}

/* deterministic: same input always gives same output */
static void test_deterministic(void) {
    uint8_t data[] = { 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08 };
    char art1[H][W + 1], art2[H][W + 1];

    cfx_randomart((char *)art1, W, H, data, sizeof(data));
    cfx_randomart((char *)art2, W, H, data, sizeof(data));

    CFX_ASSERT(memcmp(art1, art2, sizeof(art1)) == 0);
}

int main(void) {
    CFX_TEST(test_all_zeros);
    CFX_TEST(test_all_ff);
    CFX_TEST(test_markers);
    CFX_TEST(test_empty);
    CFX_TEST(test_custom_dims);
    CFX_TEST(test_deterministic);
    return 0;
}
