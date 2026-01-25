/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

/*
 * ARM Cortex-M4 optimized Poly1305 block processing.
 *
 * Uses UMULL/UMLAL instructions for the 5-term multiply-accumulate chains
 * that form the core of Poly1305 field multiplication.
 *
 * Each accumulator d[i] requires 5 multiply-accumulate operations.
 * On Cortex-M4, UMLAL executes in 1 cycle, so each accumulator takes
 * ~5 cycles vs ~20+ cycles for portable emulated 64-bit arithmetic.
 *
 * Requires: -mcpu=cortex-m4 -mthumb
 */

#include "../poly1305_backend.h"

/*
 * Poly1305 field multiply: h = (h + t) * r mod (2^130 - 5)
 *
 * ARM Cortex-M4 optimized using UMULL/UMLAL instruction chains.
 * Each 5-term accumulator is computed with a UMULL followed by 4 UMLALs.
 */
void cfx_poly1305_block_impl(
    uint32_t *h0, uint32_t *h1, uint32_t *h2, uint32_t *h3, uint32_t *h4,
    uint32_t t0, uint32_t t1, uint32_t t2, uint32_t t3, uint32_t t4,
    uint32_t r0, uint32_t r1, uint32_t r2, uint32_t r3, uint32_t r4,
    uint32_t s1, uint32_t s2, uint32_t s3, uint32_t s4){
    uint32_t c;

    /* h += t */
    uint32_t H0 = *h0 + t0;
    uint32_t H1 = *h1 + t1;
    uint32_t H2 = *h2 + t2;
    uint32_t H3 = *h3 + t3;
    uint32_t H4 = *h4 + t4;

    /*
     * Compute d[i] = 5-term multiply-accumulate using UMULL/UMLAL chains.
     *
     * d0 = H0*r0 + H1*s4 + H2*s3 + H3*s2 + H4*s1
     * d1 = H0*r1 + H1*r0 + H2*s4 + H3*s3 + H4*s2
     * d2 = H0*r2 + H1*r1 + H2*r0 + H3*s4 + H4*s3
     * d3 = H0*r3 + H1*r2 + H2*r1 + H3*r0 + H4*s4
     * d4 = H0*r4 + H1*r3 + H2*r2 + H3*r1 + H4*r0
     */

    uint32_t d0_lo, d0_hi;
    uint32_t d1_lo, d1_hi;
    uint32_t d2_lo, d2_hi;
    uint32_t d3_lo, d3_hi;
    uint32_t d4_lo, d4_hi;

    /* d0 = H0*r0 + H1*s4 + H2*s3 + H3*s2 + H4*s1 */
    __asm__ volatile (
        "umull %[lo], %[hi], %[a0], %[b0]\n\t"
        "umlal %[lo], %[hi], %[a1], %[b1]\n\t"
        "umlal %[lo], %[hi], %[a2], %[b2]\n\t"
        "umlal %[lo], %[hi], %[a3], %[b3]\n\t"
        "umlal %[lo], %[hi], %[a4], %[b4]"
        : [lo] "=&r" (d0_lo), [hi] "=&r" (d0_hi)
        : [a0] "r" (H0), [b0] "r" (r0),
        [a1] "r" (H1), [b1] "r" (s4),
        [a2] "r" (H2), [b2] "r" (s3),
        [a3] "r" (H3), [b3] "r" (s2),
        [a4] "r" (H4), [b4] "r" (s1)
        );

    /* d1 = H0*r1 + H1*r0 + H2*s4 + H3*s3 + H4*s2 */
    __asm__ volatile (
        "umull %[lo], %[hi], %[a0], %[b0]\n\t"
        "umlal %[lo], %[hi], %[a1], %[b1]\n\t"
        "umlal %[lo], %[hi], %[a2], %[b2]\n\t"
        "umlal %[lo], %[hi], %[a3], %[b3]\n\t"
        "umlal %[lo], %[hi], %[a4], %[b4]"
        : [lo] "=&r" (d1_lo), [hi] "=&r" (d1_hi)
        : [a0] "r" (H0), [b0] "r" (r1),
        [a1] "r" (H1), [b1] "r" (r0),
        [a2] "r" (H2), [b2] "r" (s4),
        [a3] "r" (H3), [b3] "r" (s3),
        [a4] "r" (H4), [b4] "r" (s2)
        );

    /* d2 = H0*r2 + H1*r1 + H2*r0 + H3*s4 + H4*s3 */
    __asm__ volatile (
        "umull %[lo], %[hi], %[a0], %[b0]\n\t"
        "umlal %[lo], %[hi], %[a1], %[b1]\n\t"
        "umlal %[lo], %[hi], %[a2], %[b2]\n\t"
        "umlal %[lo], %[hi], %[a3], %[b3]\n\t"
        "umlal %[lo], %[hi], %[a4], %[b4]"
        : [lo] "=&r" (d2_lo), [hi] "=&r" (d2_hi)
        : [a0] "r" (H0), [b0] "r" (r2),
        [a1] "r" (H1), [b1] "r" (r1),
        [a2] "r" (H2), [b2] "r" (r0),
        [a3] "r" (H3), [b3] "r" (s4),
        [a4] "r" (H4), [b4] "r" (s3)
        );

    /* d3 = H0*r3 + H1*r2 + H2*r1 + H3*r0 + H4*s4 */
    __asm__ volatile (
        "umull %[lo], %[hi], %[a0], %[b0]\n\t"
        "umlal %[lo], %[hi], %[a1], %[b1]\n\t"
        "umlal %[lo], %[hi], %[a2], %[b2]\n\t"
        "umlal %[lo], %[hi], %[a3], %[b3]\n\t"
        "umlal %[lo], %[hi], %[a4], %[b4]"
        : [lo] "=&r" (d3_lo), [hi] "=&r" (d3_hi)
        : [a0] "r" (H0), [b0] "r" (r3),
        [a1] "r" (H1), [b1] "r" (r2),
        [a2] "r" (H2), [b2] "r" (r1),
        [a3] "r" (H3), [b3] "r" (r0),
        [a4] "r" (H4), [b4] "r" (s4)
        );

    /* d4 = H0*r4 + H1*r3 + H2*r2 + H3*r1 + H4*r0 */
    __asm__ volatile (
        "umull %[lo], %[hi], %[a0], %[b0]\n\t"
        "umlal %[lo], %[hi], %[a1], %[b1]\n\t"
        "umlal %[lo], %[hi], %[a2], %[b2]\n\t"
        "umlal %[lo], %[hi], %[a3], %[b3]\n\t"
        "umlal %[lo], %[hi], %[a4], %[b4]"
        : [lo] "=&r" (d4_lo), [hi] "=&r" (d4_hi)
        : [a0] "r" (H0), [b0] "r" (r4),
        [a1] "r" (H1), [b1] "r" (r3),
        [a2] "r" (H2), [b2] "r" (r2),
        [a3] "r" (H3), [b3] "r" (r1),
        [a4] "r" (H4), [b4] "r" (r0)
        );

    /*
     * Carry chain: reduce to 26-bit limbs.
     *
     * Each d[i] is a 64-bit value. We extract the low 26 bits and
     * propagate the high bits as carry to the next accumulator.
     */

    /* c = d0 >> 26; H0 = d0 & 0x3ffffff; d1 += c */
    c = (d0_lo >> 26) | (d0_hi << 6);
    H0 = d0_lo & 0x3ffffffu;

    /* Add carry to d1 */
    __asm__ volatile (
        "adds %[lo], %[lo], %[c]\n\t"
        "adc  %[hi], %[hi], #0"
        : [lo] "+r" (d1_lo), [hi] "+r" (d1_hi)
        : [c] "r" (c)
        : "cc"
        );

    /* c = d1 >> 26; H1 = d1 & 0x3ffffff; d2 += c */
    c = (d1_lo >> 26) | (d1_hi << 6);
    H1 = d1_lo & 0x3ffffffu;

    __asm__ volatile (
        "adds %[lo], %[lo], %[c]\n\t"
        "adc  %[hi], %[hi], #0"
        : [lo] "+r" (d2_lo), [hi] "+r" (d2_hi)
        : [c] "r" (c)
        : "cc"
        );

    /* c = d2 >> 26; H2 = d2 & 0x3ffffff; d3 += c */
    c = (d2_lo >> 26) | (d2_hi << 6);
    H2 = d2_lo & 0x3ffffffu;

    __asm__ volatile (
        "adds %[lo], %[lo], %[c]\n\t"
        "adc  %[hi], %[hi], #0"
        : [lo] "+r" (d3_lo), [hi] "+r" (d3_hi)
        : [c] "r" (c)
        : "cc"
        );

    /* c = d3 >> 26; H3 = d3 & 0x3ffffff; d4 += c */
    c = (d3_lo >> 26) | (d3_hi << 6);
    H3 = d3_lo & 0x3ffffffu;

    __asm__ volatile (
        "adds %[lo], %[lo], %[c]\n\t"
        "adc  %[hi], %[hi], #0"
        : [lo] "+r" (d4_lo), [hi] "+r" (d4_hi)
        : [c] "r" (c)
        : "cc"
        );

    /* c = d4 >> 26; H4 = d4 & 0x3ffffff */
    c = (d4_lo >> 26) | (d4_hi << 6);
    H4 = d4_lo & 0x3ffffffu;

    /* Final reduction: c * 2^130 = c * 5 (mod 2^130 - 5) */
    H0 += c * 5u;
    c = H0 >> 26;
    H0 &= 0x3ffffffu;
    H1 += c;

    *h0 = H0;
    *h1 = H1;
    *h2 = H2;
    *h3 = H3;
    *h4 = H4;
}
