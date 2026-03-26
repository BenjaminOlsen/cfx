#include "cfx/chacha20.h"
#include "cfx/macros.h"

#if CFX_HAVE_AVX2
#include <immintrin.h>

/* This are some tests of some tricky transpose logic that's
 * used in the avx2 mm256 chacha20 implementation...
 */

#define LANE32(v, idx) (uint32_t)_mm_extract_epi32(v, idx)

#define DBG_P printf

static void print_m128(const char *label, __m128i v) {
    uint32_t a = LANE32(v, 0);
    uint32_t b = LANE32(v, 1);
    uint32_t c = LANE32(v, 2);
    uint32_t d = LANE32(v, 3);

    DBG_P("%s = [%u %u %u %u]\n", label, a, b, c, d);
}

static void print_m256(const char *label, __m256i v) {
    uint32_t tmp[8];
    _mm256_storeu_si256((__m256i *)tmp, v);

    DBG_P("%s = [%u %u %u %u %u %u %u %u]\n",
        label,
        tmp[0], tmp[1], tmp[2], tmp[3],
        tmp[4], tmp[5], tmp[6], tmp[7]);
}


/* This function takes 16 __m256i vectors x[0..15], each representing:
 *
 *   x[w] = [in[0][w], in[1][w], ..., in[7][w]]
 *
 * and reconstructs per-block AoS layout into out[8][16]:
 *
 *   out[b][w] = original in[b][w]
 *
 * using AVX2 + SSE transposes.
 */

/* x[w] = [ in[0][w], in[1][w], ..., in[7][w] ] */
static inline void transpose_16x8_to_blocks(const __m256i x[16], uint32_t out[8][16]) {
    __m128i lo[16], hi[16];

    DBG_P("\n---- Splitting x[w] into lo[w] and hi[w] ----\n");

    for (int w = 0; w < 16; ++w) {
        /* lo[w] = _mm256_castsi256_si128(x[w]);     */
        lo[w] = _mm256_extracti128_si256(x[w], 0); /* lanes 0..3 */
        hi[w] = _mm256_extracti128_si256(x[w], 1);  /* lanes 4..7 */

        DBG_P("w=%02d: ", w);
        print_m256("x[w]", x[w]);
        print_m128("lo[w]", lo[w]);
        print_m128("hi[w]", hi[w]);
    }

    uint32_t *dst;

#define EXTRACT_BLOCK(ARR, lane_idx, out_idx)                           \
        do {                                                                \
            DBG_P("\n-- Extract block %d from %s lane %d --\n",            \
        (out_idx), #ARR, (lane_idx));                            \
            dst = out[(out_idx)];                                           \
            for (int w = 0; w < 16; ++w) {                                  \
                dst[w] = LANE32((ARR)[w], (lane_idx));                      \
                DBG_P("out[%d][%02d] = %u\n", out_idx, w, dst[w]);         \
            }                                                               \
        } while (0)

    /* Blocks 0..3 from lo */
    EXTRACT_BLOCK(lo, 0, 0);
    EXTRACT_BLOCK(lo, 1, 1);
    EXTRACT_BLOCK(lo, 2, 2);
    EXTRACT_BLOCK(lo, 3, 3);

    /* Blocks 4..7 from hi */
    EXTRACT_BLOCK(hi, 0, 4);
    EXTRACT_BLOCK(hi, 1, 5);
    EXTRACT_BLOCK(hi, 2, 6);
    EXTRACT_BLOCK(hi, 3, 7);

#undef EXTRACT_BLOCK
#undef LANE32
}


static void test_avx2_transpose(void) {
    /* Original AoS data: 8 blocks, 16 words each.
       Pattern: in[blk][word] = blk*1000 + word */
    uint32_t in[8][16];
    for (int blk = 0; blk < 8; ++blk) {
        for (int w = 0; w < 16; ++w) {
            in[blk][w] = (uint32_t)(blk * 1000 + w);
        }
    }

    /* Build SoA: x[w] = vector of 8 blocks' word w */
    __m256i x[16];
    for (int w = 0; w < 16; ++w) {
        x[w] = _mm256_setr_epi32(
            (int)in[0][w], (int)in[1][w], (int)in[2][w], (int)in[3][w],
            (int)in[4][w], (int)in[5][w], (int)in[6][w], (int)in[7][w]
            );
    }

    uint32_t out[8][16] = {0};
    transpose_16x8_to_blocks(x, out);

    int errors = 0;
    for (int blk = 0; blk < 8; ++blk) {
        for (int w = 0; w < 16; ++w) {
            uint32_t expect = in[blk][w];
            uint32_t got    = out[blk][w];
            if (expect != got) {
                DBG_P("Mismatch at block %d, word %d: expect %u, got %u\n",
                    blk, w, expect, got);
                ++errors;
            }
        }
    }

    if (errors == 0) {
        printf("-- avx2 transpose test: OK\n");
    } else {
        printf("-- avx2 transpose test: %d mismatches\n", errors);
    }

    CFX_ASSERT(errors == 0);
}
#endif

int main(void) {

#if CFX_HAVE_AVX2
    CFX_TEST(test_avx2_transpose);
#endif

    return 0;
}
