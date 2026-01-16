/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

/*
 * Portable Poly1305 block processing.
 *
 * This is the reference implementation using standard C arithmetic.
 * The 26-bit radix representation allows using 32-bit multiplies
 * that fit in 64-bit accumulators without overflow.
 */

#include "../poly1305_backend.h"

/*
 * Poly1305 field multiply: h = (h + t) * r mod (2^130 - 5)
 *
 * Uses 26-bit limbs for both operands.
 * The magic of 2^130 - 5: overflow terms get multiplied by 5 and folded back.
 */
void cfx_poly1305_block_impl(
    uint32_t *h0, uint32_t *h1, uint32_t *h2, uint32_t *h3, uint32_t *h4,
    uint32_t  t0, uint32_t  t1, uint32_t  t2, uint32_t  t3, uint32_t  t4,
    uint32_t  r0, uint32_t  r1, uint32_t  r2, uint32_t  r3, uint32_t  r4,
    uint32_t  s1, uint32_t  s2, uint32_t  s3, uint32_t  s4)
{
    uint64_t d0, d1, d2, d3, d4;
    uint32_t c;

    /* h += t */
    uint32_t H0 = *h0 + t0;
    uint32_t H1 = *h1 + t1;
    uint32_t H2 = *h2 + t2;
    uint32_t H3 = *h3 + t3;
    uint32_t H4 = *h4 + t4;

    /* h *= r (with reduction via s[i] = 5*r[i]) */
    d0 = (uint64_t)H0 * r0
       + (uint64_t)H1 * s4
       + (uint64_t)H2 * s3
       + (uint64_t)H3 * s2
       + (uint64_t)H4 * s1;

    d1 = (uint64_t)H0 * r1
       + (uint64_t)H1 * r0
       + (uint64_t)H2 * s4
       + (uint64_t)H3 * s3
       + (uint64_t)H4 * s2;

    d2 = (uint64_t)H0 * r2
       + (uint64_t)H1 * r1
       + (uint64_t)H2 * r0
       + (uint64_t)H3 * s4
       + (uint64_t)H4 * s3;

    d3 = (uint64_t)H0 * r3
       + (uint64_t)H1 * r2
       + (uint64_t)H2 * r1
       + (uint64_t)H3 * r0
       + (uint64_t)H4 * s4;

    d4 = (uint64_t)H0 * r4
       + (uint64_t)H1 * r3
       + (uint64_t)H2 * r2
       + (uint64_t)H3 * r1
       + (uint64_t)H4 * r0;

    /* Carry chain: keep limbs < 2^26 */
    c  = (uint32_t)(d0 >> 26); H0 = (uint32_t)d0 & 0x3ffffffu; d1 += c;
    c  = (uint32_t)(d1 >> 26); H1 = (uint32_t)d1 & 0x3ffffffu; d2 += c;
    c  = (uint32_t)(d2 >> 26); H2 = (uint32_t)d2 & 0x3ffffffu; d3 += c;
    c  = (uint32_t)(d3 >> 26); H3 = (uint32_t)d3 & 0x3ffffffu; d4 += c;
    c  = (uint32_t)(d4 >> 26); H4 = (uint32_t)d4 & 0x3ffffffu;

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
