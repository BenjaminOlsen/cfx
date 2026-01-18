/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

/*
 * ARM Cortex-M4 optimized ChaCha20 block generation.
 *
 * ChaCha20 is inherently 32-bit friendly. The Cortex-M4's barrel shifter
 * makes the rotation operations essentially free when combined with other
 * operations (XOR + rotate is one instruction cycle).
 *
 * Requires: -mcpu=cortex-m4 -mthumb
 *
 * Self-guarding: only compiles when CFX_TARGET_ARM_CORTEX_M4 is defined.
 */

#ifdef CFX_TARGET_ARM_CORTEX_M4

#include "chacha20_backend.h"
#include "cfx/chacha20.h"
#include "cfx/arch.h"

#include <string.h>

/* Inline rotation macros using barrel shifter */
#define ROTL16(x) ({ \
    uint32_t _r; \
    __asm__ ("ror %[out], %[in], #16" : [out] "=r" (_r) : [in] "r" (x)); \
    _r; \
})

#define ROTL12(x) ({ \
    uint32_t _r; \
    __asm__ ("ror %[out], %[in], #20" : [out] "=r" (_r) : [in] "r" (x)); \
    _r; \
})

#define ROTL8(x) ({ \
    uint32_t _r; \
    __asm__ ("ror %[out], %[in], #24" : [out] "=r" (_r) : [in] "r" (x)); \
    _r; \
})

#define ROTL7(x) ({ \
    uint32_t _r; \
    __asm__ ("ror %[out], %[in], #25" : [out] "=r" (_r) : [in] "r" (x)); \
    _r; \
})

/*
 * Optimized quarter-round macro.
 */
#define QR_M4(a, b, c, d) do { \
    (a) += (b); \
    (d) ^= (a); \
    (d) = ROTL16(d); \
    (c) += (d); \
    (b) ^= (c); \
    (b) = ROTL12(b); \
    (a) += (b); \
    (d) ^= (a); \
    (d) = ROTL8(d); \
    (c) += (d); \
    (b) ^= (c); \
    (b) = ROTL7(b); \
} while (0)

/*
 * Generate a single ChaCha20 block with ARM-optimized quarter-rounds.
 */
static void chacha20_block_m4(const uint32_t s[16], uint32_t counter, uint32_t w[16])
{
    /* Copy state to working array */
    for (int i = 0; i < 16; ++i) {
        w[i] = s[i];
    }
    w[12] = counter;

    /* 20 rounds = 10 double-rounds */
    for (int i = 0; i < 10; ++i) {
        /* Column rounds */
        QR_M4(w[0], w[4], w[8],  w[12]);
        QR_M4(w[1], w[5], w[9],  w[13]);
        QR_M4(w[2], w[6], w[10], w[14]);
        QR_M4(w[3], w[7], w[11], w[15]);

        /* Diagonal rounds */
        QR_M4(w[0], w[5], w[10], w[15]);
        QR_M4(w[1], w[6], w[11], w[12]);
        QR_M4(w[2], w[7], w[8],  w[13]);
        QR_M4(w[3], w[4], w[9],  w[14]);
    }

    /* Add original state */
    for (int i = 0; i < 12; ++i) {
        w[i] += s[i];
    }
    w[12] += counter;
    w[13] += s[13];
    w[14] += s[14];
    w[15] += s[15];
}

/*
 * Generate 8 ChaCha20 blocks.
 *
 * For Cortex-M4 without SIMD, we generate blocks sequentially.
 * The optimization comes from the efficient quarter-round implementation.
 */
void cfx_chacha20_block8_impl(const uint8_t key[32],
                              uint32_t counter,
                              const uint8_t nonce[12],
                              uint8_t out[8][64])
{
    /* Initialize state once */
    uint32_t s[16];

    /* Constants "expand 32-byte k" */
    s[0]  = 0x61707865u;
    s[1]  = 0x3320646eu;
    s[2]  = 0x79622d32u;
    s[3]  = 0x6b206574u;

    /* Key (8 words) */
    s[4]  = CFX_LOAD32_LE(key + 0);
    s[5]  = CFX_LOAD32_LE(key + 4);
    s[6]  = CFX_LOAD32_LE(key + 8);
    s[7]  = CFX_LOAD32_LE(key + 12);
    s[8]  = CFX_LOAD32_LE(key + 16);
    s[9]  = CFX_LOAD32_LE(key + 20);
    s[10] = CFX_LOAD32_LE(key + 24);
    s[11] = CFX_LOAD32_LE(key + 28);

    /* Counter placeholder (set per block) */
    s[12] = 0;

    /* Nonce (3 words) */
    s[13] = CFX_LOAD32_LE(nonce + 0);
    s[14] = CFX_LOAD32_LE(nonce + 4);
    s[15] = CFX_LOAD32_LE(nonce + 8);

    /* Generate 8 blocks sequentially */
    for (int block = 0; block < 8; ++block) {
        uint32_t w[16];

        chacha20_block_m4(s, counter + (uint32_t)block, w);

        /* Serialize output in little-endian */
        for (int i = 0; i < 16; ++i) {
            CFX_STORE32_LE(out[block] + 4 * i, w[i]);
        }
    }
}

#endif /* CFX_TARGET_ARM_CORTEX_M4 */
