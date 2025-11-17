#include "cfx/poly1305.h"
#include "cfx/memory.h"

/* Watch out for macro collisions: */
#ifdef S
#message ""
#define OLD_S S
#undef S
#endif

#ifdef X
#define OLD_X X
#undef X
#endif

/* sum: */
#define SUM(a,b,c,d,e) ((a)+(b)+(c)+(d)+(e))

/* cast to accumulator, multiply */
#define MUL(a,b) ((uint64_t)(a)*(uint64_t)(b))


/* ---------------------------------------------------------------------------------------------------- */
/* The amazing poly1305 message authenticator algorithm - shamelessly plagiarized from libsodium but copied here
 * to feel the masters' hands on my keyboard...
 * ----- */

/* tag = Poly1305(key, m) per RFC8439 A.3
 *     = m1*r^n + m2*r(n-1) + ... + mn*r + k mod p with p = 2^130 - 5
 * key: 32 bytes (r || s)
 * m  : message
 * mlen: length in bytes
 * tag: 16 byte MAC
 */
void cfx_poly1305_mac(const uint8_t key[32], const uint8_t *m, size_t mlen, uint8_t tag[16]) {
    /*  clamp(r): r &= 0x0ffffffc0ffffffc0ffffffc0fffffff */
    uint32_t r0 =  cfx_load32_le(&key[0])  & 0x3ffffffu;
    uint32_t r1 = (cfx_load32_le(&key[3])  >> 2) & 0x3ffff03u;
    uint32_t r2 = (cfx_load32_le(&key[6])  >> 4) & 0x3ffc0ffu;
    uint32_t r3 = (cfx_load32_le(&key[9])  >> 6) & 0x3f03fffu;
    uint32_t r4 = (cfx_load32_le(&key[12]) >> 8) & 0x00fffffu;

    /* precompute the 5*r[i] terms: */
    uint64_t s1 = (uint64_t)r1 * 5u;
    uint64_t s2 = (uint64_t)r2 * 5u;
    uint64_t s3 = (uint64_t)r3 * 5u;
    uint64_t s4 = (uint64_t)r4 * 5u;

    /* pad s = second half of the key */
    uint32_t pad0 = cfx_load32_le(&key[16]);
    uint32_t pad1 = cfx_load32_le(&key[20]);
    uint32_t pad2 = cfx_load32_le(&key[24]);
    uint32_t pad3 = cfx_load32_le(&key[28]);

    /* accumulator h in 26-bit limbs */
    uint32_t h0 = 0, h1 = 0, h2 = 0, h3 = 0, h4 = 0;

    const uint8_t *p = m;
    size_t bytes = mlen;

    while (bytes >= 16) { /* full 16-byte blocks (each with hibit set) */
        uint32_t t0, t1, t2, t3, t4;
        uint64_t d0, d1, d2, d3, d4;
        uint32_t c;
        const uint32_t hibit = 1u << 24;  /* 2^128 term for full blocks */

        /* h += m[i] (decoded into 26-bit limbs) */
        t0 =  cfx_load32_le(p + 0)  & 0x3ffffffu;
        t1 = (cfx_load32_le(p + 3)  >> 2) & 0x3ffffffu;
        t2 = (cfx_load32_le(p + 6)  >> 4) & 0x3ffffffu;
        t3 = (cfx_load32_le(p + 9)  >> 6) & 0x3ffffffu;
        t4 = (cfx_load32_le(p + 12) >> 8) | hibit;

        h0 += t0;
        h1 += t1;
        h2 += t2;
        h3 += t3;
        h4 += t4;

        /* h *= r  (mod 2^130-5), with 32×32 -> 64 bit multiplies */
        d0 = SUM(MUL(h0,r0), MUL(h1,s4), MUL(h2,s3), MUL(h3,s2), MUL(h4,s1));
        d1 = SUM(MUL(h0,r1), MUL(h1,r0), MUL(h2,s4), MUL(h3,s3), MUL(h4,s2));
        d2 = SUM(MUL(h0,r2), MUL(h1,r1), MUL(h2,r0), MUL(h3,s4), MUL(h4,s3));
        d3 = SUM(MUL(h0,r3), MUL(h1,r2), MUL(h2,r1), MUL(h3,r0), MUL(h4,s4));
        d4 = SUM(MUL(h0,r4), MUL(h1,r3), MUL(h2,r2), MUL(h3,r1), MUL(h4,r0));

        /* (partial) h %= p : carry chain, keep limbs < 2^26 */
        c  = (uint32_t)(d0 >> 26); h0 = (uint32_t)d0 & 0x3ffffffu; d1 += c;
        c  = (uint32_t)(d1 >> 26); h1 = (uint32_t)d1 & 0x3ffffffu; d2 += c;
        c  = (uint32_t)(d2 >> 26); h2 = (uint32_t)d2 & 0x3ffffffu; d3 += c;
        c  = (uint32_t)(d3 >> 26); h3 = (uint32_t)d3 & 0x3ffffffu; d4 += c;
        c  = (uint32_t)(d4 >> 26); h4 = (uint32_t)d4 & 0x3ffffffu;

        h0 += c * 5u;
        c   = h0 >> 26;
        h0 &= 0x3ffffffu;
        h1 += c;

        p     += 16;
        bytes -= 16;
    }

    /* ----- final partial block (if any, with explicit 0x01 pad)  */
    if (bytes > 0) {
        uint8_t  buf[16] = {0};
        uint32_t t0, t1, t2, t3, t4;
        uint64_t d0, d1, d2, d3, d4;
        uint32_t c;
        const uint32_t hibit = 0u;  /* no 2^128 term for partial block */

        /* copy leftover bytes, append 0x01, zero-pad */
        for (size_t i = 0; i < bytes; i++) {
            buf[i] = p[i];
        }
        buf[bytes] = 1;

        /* h += m_last */
        t0 =  cfx_load32_le(buf + 0)  & 0x3ffffffu;
        t1 = (cfx_load32_le(buf + 3)  >> 2) & 0x3ffffffu;
        t2 = (cfx_load32_le(buf + 6)  >> 4) & 0x3ffffffu;
        t3 = (cfx_load32_le(buf + 9)  >> 6) & 0x3ffffffu;
        t4 = (cfx_load32_le(buf + 12) >> 8) | hibit;

        h0 += t0;
        h1 += t1;
        h2 += t2;
        h3 += t3;
        h4 += t4;

        /* h *= r (same as above) */
        d0 = SUM(MUL(h0,r0), MUL(h1,s4), MUL(h2,s3), MUL(h3,s2), MUL(h4,s1));
        d1 = SUM(MUL(h0,r1), MUL(h1,r0), MUL(h2,s4), MUL(h3,s3), MUL(h4,s2));
        d2 = SUM(MUL(h0,r2), MUL(h1,r1), MUL(h2,r0), MUL(h3,s4), MUL(h4,s3));
        d3 = SUM(MUL(h0,r3), MUL(h1,r2), MUL(h2,r1), MUL(h3,r0), MUL(h4,s4));
        d4 = SUM(MUL(h0,r4), MUL(h1,r3), MUL(h2,r2), MUL(h3,r1), MUL(h4,r0));

        /* carry */
        c  = (uint32_t)(d0 >> 26); h0 = (uint32_t)d0 & 0x3ffffffu; d1 += c;
        c  = (uint32_t)(d1 >> 26); h1 = (uint32_t)d1 & 0x3ffffffu; d2 += c;
        c  = (uint32_t)(d2 >> 26); h2 = (uint32_t)d2 & 0x3ffffffu; d3 += c;
        c  = (uint32_t)(d3 >> 26); h3 = (uint32_t)d3 & 0x3ffffffu; d4 += c;
        c  = (uint32_t)(d4 >> 26); h4 = (uint32_t)d4 & 0x3ffffffu;

        /* last carry*/
        h0 += c * 5u;
        c   = h0 >> 26;
        h0 &= 0x3ffffffu;
        h1 += c;

        /* scrub buf */
        cfx_memzero_s(buf, sizeof(buf));
    }

    /* final reduction: carry fully, then conditional subtract p */
    uint32_t c;

    /* carry h into canonical 26-bit limbs */
    c  = h1 >> 26; h1 &= 0x3ffffffu; h2 += c;
    c  = h2 >> 26; h2 &= 0x3ffffffu; h3 += c;
    c  = h3 >> 26; h3 &= 0x3ffffffu; h4 += c;
    c  = h4 >> 26; h4 &= 0x3ffffffu; h0 += c * 5u;
    c  = h0 >> 26; h0 &= 0x3ffffffu; h1 += c;

    /* compute g = h + 5 - 2^130 */
    uint32_t g0 = h0 + 5u;
    c  = g0 >> 26; g0 &= 0x3ffffffu;
    uint32_t g1 = h1 + c;
    c  = g1 >> 26; g1 &= 0x3ffffffu;
    uint32_t g2 = h2 + c;
    c  = g2 >> 26; g2 &= 0x3ffffffu;
    uint32_t g3 = h3 + c;
    c  = g3 >> 26; g3 &= 0x3ffffffu;
    uint32_t g4 = h4 + c - (1u << 26);

    /* mask = all ones if h >= p, else 0 */
    uint32_t mask  = (g4 >> 31) - 1u;
    uint32_t nmask = ~mask;

    /* select h if h < p, else g */
    h0 = (h0 & nmask) | (g0 & mask);
    h1 = (h1 & nmask) | (g1 & mask);
    h2 = (h2 & nmask) | (g2 & mask);
    h3 = (h3 & nmask) | (g3 & mask);
    h4 = (h4 & nmask) | (g4 & mask);

    /* pack h into 128 bits and add pad (s) */
    /* h = h % (2^128) into t: four 32 bit words */
    uint32_t t0 = ( h0        | (h1 << 26))        & 0xffffffffu;
    uint32_t t1 = ((h1 >> 6)  | (h2 << 20))        & 0xffffffffu;
    uint32_t t2 = ((h2 >> 12) | (h3 << 14))        & 0xffffffffu;
    uint32_t t3 = ((h3 >> 18) | (h4 << 8))         & 0xffffffffu;

    /* mac = (h + pad) % 2^128 */
    uint64_t f;

    f  = (uint64_t)t0 + pad0;
    t0 = (uint32_t)f;
    f  = (uint64_t)t1 + pad1 + (f >> 32);
    t1 = (uint32_t)f;
    f  = (uint64_t)t2 + pad2 + (f >> 32);
    t2 = (uint32_t)f;
    f  = (uint64_t)t3 + pad3 + (f >> 32);
    t3 = (uint32_t)f;

    cfx_store32_le(tag + 0, t0);
    cfx_store32_le(tag + 4, t1);
    cfx_store32_le(tag + 8, t2);
    cfx_store32_le(tag + 12, t3);

    /* scrub sensitive locals: */
    cfx_memzero_s(&r0, sizeof(r0));
    cfx_memzero_s(&r1, sizeof(r1));
    cfx_memzero_s(&r2, sizeof(r2));
    cfx_memzero_s(&r3, sizeof(r3));
    cfx_memzero_s(&r4, sizeof(r4));

    cfx_memzero_s(&s1, sizeof(s1));
    cfx_memzero_s(&s2, sizeof(s2));
    cfx_memzero_s(&s3, sizeof(s3));
    cfx_memzero_s(&s4, sizeof(s4));

    cfx_memzero_s(&h0, sizeof(h0));
    cfx_memzero_s(&h1, sizeof(h1));
    cfx_memzero_s(&h2, sizeof(h2));
    cfx_memzero_s(&h3, sizeof(h3));
    cfx_memzero_s(&h4, sizeof(h4));

    cfx_memzero_s(&pad0, sizeof(pad0));
    cfx_memzero_s(&pad1, sizeof(pad1));
    cfx_memzero_s(&pad2, sizeof(pad2));
    cfx_memzero_s(&pad3, sizeof(pad3));
    cfx_memzero_s(&t0, sizeof(t0));
    cfx_memzero_s(&t1, sizeof(t1));
    cfx_memzero_s(&t2, sizeof(t2));
    cfx_memzero_s(&t3, sizeof(t3));
    cfx_memzero_s(&f,  sizeof(f));
}

#ifdef OLD_S
    #define S OLD_S
    #undef OLD_S
#endif
#ifdef OLD_X
    #define X OLD_X
    #undef OLD_X
#endif
