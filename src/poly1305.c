#include "cfx/poly1305.h"
#include "cfx/memory.h"
#include "cfx/arch.h"

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
/* The amazing poly1305 message authenticator algorithm - shamelessly plagiarized from libsodium.

the '1305' comes from  p = 2^130 - 5, so we can do some tricks:

if a and b are expressed in 5 x 26  bit limbs, then
we can take advantage of the fact that the top column (in parens below)
which are in the (2^26)^5 = 2^130 place. So if we have a quantity n in that
term it's:
  n * 2^130 mod p
= n*2^130 mod (2^130 - 5)
= (n mod p) * (2^130 mod p)
= (n mod p) * 5

= n * 5     - (if 0 <= n < p, then n mod p = n.)

This makes the multiplication of two 5 limb numbers look like this:

                                   a4      a3      a2      a1      a0
                              ×    b4      b3      b2      b1      b0
                              ---------------------------------------
                                a4×b0   a3×b0   a2×b0   a1×b0   a0×b0
                     (a4*b1)  + a3×b1   a2×b1   a1×b1   a0×b1 5×a4×b1  <--
              ...  + (a3*b2)  + a2×b2   a1×b2   a0×b2 5×a4×b2 5×a3×b2  <--
        ... + ...  + (a2*b3)  + a1×b3   a0×b3 5×a4×b3 5×a3×b3 5×a2×b3
 ... + .... + ...  + (a1*b4)  + a0×b4 5×a4×b4 5×a3×b4 5×a2×b4 5×a1×b4

 * ----- */

/* tag = Poly1305(key, m) per RFC8439 A.3
 *     = m1*r^n + m2*r(n-1) + ... + mn*r + k mod p with p = 2^130 - 5
 * key: 32 bytes (r || s)
 * m  : message
 * mlen: length in bytes
 * tag: 16 byte MAC
 */


void cfx_poly1305_init(cfx_poly1305_state_t* s, const uint8_t key[32]) {

    /*  clamp(r): r &= 0x0ffffffc0ffffffc0ffffffc0fffffff */
    s->r0 =  CFX_LOAD32_LE(&key[0])  & 0x3ffffffu;
    s->r1 = (CFX_LOAD32_LE(&key[3])  >> 2) & 0x3ffff03u;
    s->r2 = (CFX_LOAD32_LE(&key[6])  >> 4) & 0x3ffc0ffu;
    s->r3 = (CFX_LOAD32_LE(&key[9])  >> 6) & 0x3f03fffu;
    s->r4 = (CFX_LOAD32_LE(&key[12]) >> 8) & 0x00fffffu;
    s->s1 = (uint64_t)s->r1 * 5u;
    s->s2 = (uint64_t)s->r2 * 5u;
    s->s3 = (uint64_t)s->r3 * 5u;
    s->s4 = (uint64_t)s->r4 * 5u;
    s->pad0 = CFX_LOAD32_LE(&key[16]);
    s->pad1 = CFX_LOAD32_LE(&key[20]);
    s->pad2 = CFX_LOAD32_LE(&key[24]);
    s->pad3 = CFX_LOAD32_LE(&key[28]);
    s->h0 = 0;
    s->h1 = 0;
    s->h2 = 0;
    s->h3 = 0;
    s->h4 = 0;

    memset(s->buffer, 0, sizeof s->buffer);
    s->buflen = 0;
    s->done = 0;
}

CFX_INLINE void cfx_poly1305_block(
    uint32_t *h0, uint32_t *h1, uint32_t *h2, uint32_t *h3, uint32_t *h4,
    uint32_t  r0, uint32_t  r1, uint32_t  r2, uint32_t  r3, uint32_t  r4,
    uint64_t  s1, uint64_t  s2, uint64_t  s3, uint64_t  s4,
    const uint8_t *p)
{
    uint32_t t0, t1, t2, t3, t4;
    uint64_t d0, d1, d2, d3, d4;
    uint32_t c;
    const uint32_t hibit = 1u << 24;  /* 2^128 term for full blocks */

    uint32_t H0 = *h0;
    uint32_t H1 = *h1;
    uint32_t H2 = *h2;
    uint32_t H3 = *h3;
    uint32_t H4 = *h4;

    /* h += m[i] (decoded into 26-bit limbs) */
    t0 =  CFX_LOAD32_LE(p + 0)  & 0x3ffffffu;
    t1 = (CFX_LOAD32_LE(p + 3)  >> 2) & 0x3ffffffu;
    t2 = (CFX_LOAD32_LE(p + 6)  >> 4) & 0x3ffffffu;
    t3 = (CFX_LOAD32_LE(p + 9)  >> 6) & 0x3ffffffu;
    t4 = (CFX_LOAD32_LE(p + 12) >> 8) | hibit;

    H0 += t0;
    H1 += t1;
    H2 += t2;
    H3 += t3;
    H4 += t4;

    /* h *= r  (mod 2^130-5), with 32×32 -> 64 bit multiplies */
    d0 = SUM(MUL(H0,r0), MUL(H1,s4), MUL(H2,s3), MUL(H3,s2), MUL(H4,s1));
    d1 = SUM(MUL(H0,r1), MUL(H1,r0), MUL(H2,s4), MUL(H3,s3), MUL(H4,s2));
    d2 = SUM(MUL(H0,r2), MUL(H1,r1), MUL(H2,r0), MUL(H3,s4), MUL(H4,s3));
    d3 = SUM(MUL(H0,r3), MUL(H1,r2), MUL(H2,r1), MUL(H3,r0), MUL(H4,s4));
    d4 = SUM(MUL(H0,r4), MUL(H1,r3), MUL(H2,r2), MUL(H3,r1), MUL(H4,r0));

    /* (partial) h %= p : carry chain, keep limbs < 2^26 */
    c  = (uint32_t)(d0 >> 26); H0 = (uint32_t)d0 & 0x3ffffffu; d1 += c;
    c  = (uint32_t)(d1 >> 26); H1 = (uint32_t)d1 & 0x3ffffffu; d2 += c;
    c  = (uint32_t)(d2 >> 26); H2 = (uint32_t)d2 & 0x3ffffffu; d3 += c;
    c  = (uint32_t)(d3 >> 26); H3 = (uint32_t)d3 & 0x3ffffffu; d4 += c;
    c  = (uint32_t)(d4 >> 26); H4 = (uint32_t)d4 & 0x3ffffffu;

    H0 += c * 5u;
    c   = H0 >> 26;
    H0 &= 0x3ffffffu;
    H1 += c;

    /* write back */
    *h0 = H0;
    *h1 = H1;
    *h2 = H2;
    *h3 = H3;
    *h4 = H4;
}

void cfx_poly1305_update(cfx_poly1305_state_t* s, const uint8_t* msg, size_t mlen) {
    if (s->done || mlen == 0) { return; }
    uint8_t buflen = s->buflen;

    uint32_t h0 = s->h0, h1 = s->h1, h2 = s->h2, h3 = s->h3, h4 = s->h4;
    const uint32_t r0 = s->r0, r1 = s->r1, r2 = s->r2, r3 = s->r3, r4 = s->r4;
    const uint64_t s1 = s->s1, s2 = s->s2, s3 = s->s3, s4 = s->s4;

    if (buflen) {
        size_t remain = 16u - buflen;
        if (mlen < remain) {
            memcpy(s->buffer + buflen, msg, mlen);
            s->buflen += mlen;
            return;
        }

        memcpy(s->buffer + buflen, msg, remain);
        cfx_poly1305_block(&h0, &h1, &h2, &h3, &h4,
                            r0, r1, r2, r3, r4,
                            s1, s2, s3, s4,
                            s->buffer);
        msg += remain;
        mlen -= remain;
        s->buflen = 0;
    }

    while (mlen >= 16) { /* full 16-byte blocks  */
        cfx_poly1305_block(&h0, &h1, &h2, &h3, &h4,
                            r0, r1, r2, r3, r4,
                            s1, s2, s3, s4,
                            msg);
        msg  += 16;
        mlen -= 16;
    }

    /* the rest into s->buffer */
    if (mlen > 0) {
        memcpy(s->buffer, msg, mlen);
        s->buflen = mlen;
    }

    s->h0 = h0;
    s->h1 = h1;
    s->h2 = h2;
    s->h3 = h3;
    s->h4 = h4;
}

void cfx_poly1305_finish(cfx_poly1305_state_t* s, uint8_t tag[16]) {
    if (s->done) return;

    uint32_t r0 = s->r0, r1 = s->r1, r2 = s->r2, r3 = s->r3, r4 = s->r4;
    uint64_t s1 = s->s1, s2 = s->s2, s3 = s->s3, s4 = s->s4;
    uint32_t pad0 = s->pad0, pad1 = s->pad1, pad2 = s->pad2, pad3 = s->pad3;
    uint32_t h0 = s->h0, h1 = s->h1, h2 = s->h2, h3 = s->h3, h4 = s->h4;

    if (s->buflen) {
        uint8_t buf[16] = {0};
        uint32_t t0, t1, t2, t3, t4;
        uint64_t d0, d1, d2, d3, d4;
        uint32_t c;
        const uint32_t hibit = 0u;  /* no 2^128 term for partial block */

        /* copy leftover bytes, append 0x01, zero-pad */
        for (size_t i = 0; i < s->buflen; i++) {
            buf[i] = s->buffer[i];
        }
        buf[s->buflen] = 1;

        /* h += m_last */
        t0 =  CFX_LOAD32_LE(buf + 0)  & 0x3ffffffu;
        t1 = (CFX_LOAD32_LE(buf + 3)  >> 2) & 0x3ffffffu;
        t2 = (CFX_LOAD32_LE(buf + 6)  >> 4) & 0x3ffffffu;
        t3 = (CFX_LOAD32_LE(buf + 9)  >> 6) & 0x3ffffffu;
        t4 = (CFX_LOAD32_LE(buf + 12) >> 8) | hibit;

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
        CFX_MEMZERO_S(buf, sizeof(buf));
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

    CFX_STORE32_LE(tag + 0, t0);
    CFX_STORE32_LE(tag + 4, t1);
    CFX_STORE32_LE(tag + 8, t2);
    CFX_STORE32_LE(tag + 12, t3);

    /* scrub sensitive locals: */
    CFX_MEMZERO_S(&r0, sizeof(r0));
    CFX_MEMZERO_S(&r1, sizeof(r1));
    CFX_MEMZERO_S(&r2, sizeof(r2));
    CFX_MEMZERO_S(&r3, sizeof(r3));
    CFX_MEMZERO_S(&r4, sizeof(r4));

    CFX_MEMZERO_S(&s1, sizeof(s1));
    CFX_MEMZERO_S(&s2, sizeof(s2));
    CFX_MEMZERO_S(&s3, sizeof(s3));
    CFX_MEMZERO_S(&s4, sizeof(s4));

    CFX_MEMZERO_S(&h0, sizeof(h0));
    CFX_MEMZERO_S(&h1, sizeof(h1));
    CFX_MEMZERO_S(&h2, sizeof(h2));
    CFX_MEMZERO_S(&h3, sizeof(h3));
    CFX_MEMZERO_S(&h4, sizeof(h4));

    CFX_MEMZERO_S(&pad0, sizeof(pad0));
    CFX_MEMZERO_S(&pad1, sizeof(pad1));
    CFX_MEMZERO_S(&pad2, sizeof(pad2));
    CFX_MEMZERO_S(&pad3, sizeof(pad3));
    CFX_MEMZERO_S(&t0, sizeof(t0));
    CFX_MEMZERO_S(&t1, sizeof(t1));
    CFX_MEMZERO_S(&t2, sizeof(t2));
    CFX_MEMZERO_S(&t3, sizeof(t3));
    CFX_MEMZERO_S(&f,  sizeof(f));

    CFX_MEMZERO_S(s,  sizeof(*s));
    s->done = 1;
}

void cfx_poly1305_mac_2(const uint8_t key[32], const uint8_t* m, size_t mlen, uint8_t tag[16]) {
    cfx_poly1305_state_t st;
    cfx_poly1305_init(&st, key);
    cfx_poly1305_update(&st, m, mlen);
    cfx_poly1305_finish(&st, tag);
}

void cfx_poly1305_mac(const uint8_t key[32], const uint8_t* m, size_t mlen, uint8_t tag[16]) {
    /*  clamp(r): r &= 0x0ffffffc0ffffffc0ffffffc0fffffff */
    uint32_t r0 =  CFX_LOAD32_LE(&key[0])  & 0x3ffffffu;
    uint32_t r1 = (CFX_LOAD32_LE(&key[3])  >> 2) & 0x3ffff03u;
    uint32_t r2 = (CFX_LOAD32_LE(&key[6])  >> 4) & 0x3ffc0ffu;
    uint32_t r3 = (CFX_LOAD32_LE(&key[9])  >> 6) & 0x3f03fffu;
    uint32_t r4 = (CFX_LOAD32_LE(&key[12]) >> 8) & 0x00fffffu;

    /* precompute the 5*r[i] terms: */
    uint64_t s1 = (uint64_t)r1 * 5u;
    uint64_t s2 = (uint64_t)r2 * 5u;
    uint64_t s3 = (uint64_t)r3 * 5u;
    uint64_t s4 = (uint64_t)r4 * 5u;

    /* pad s = second half of the key */
    uint32_t pad0 = CFX_LOAD32_LE(&key[16]);
    uint32_t pad1 = CFX_LOAD32_LE(&key[20]);
    uint32_t pad2 = CFX_LOAD32_LE(&key[24]);
    uint32_t pad3 = CFX_LOAD32_LE(&key[28]);

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
        t0 =  CFX_LOAD32_LE(p + 0)  & 0x3ffffffu;
        t1 = (CFX_LOAD32_LE(p + 3)  >> 2) & 0x3ffffffu;
        t2 = (CFX_LOAD32_LE(p + 6)  >> 4) & 0x3ffffffu;
        t3 = (CFX_LOAD32_LE(p + 9)  >> 6) & 0x3ffffffu;
        t4 = (CFX_LOAD32_LE(p + 12) >> 8) | hibit;

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
        t0 =  CFX_LOAD32_LE(buf + 0)  & 0x3ffffffu;
        t1 = (CFX_LOAD32_LE(buf + 3)  >> 2) & 0x3ffffffu;
        t2 = (CFX_LOAD32_LE(buf + 6)  >> 4) & 0x3ffffffu;
        t3 = (CFX_LOAD32_LE(buf + 9)  >> 6) & 0x3ffffffu;
        t4 = (CFX_LOAD32_LE(buf + 12) >> 8) | hibit;

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
        CFX_MEMZERO_S(buf, sizeof(buf));
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

    CFX_STORE32_LE(tag + 0, t0);
    CFX_STORE32_LE(tag + 4, t1);
    CFX_STORE32_LE(tag + 8, t2);
    CFX_STORE32_LE(tag + 12, t3);

    /* scrub sensitive locals: */
    CFX_MEMZERO_S(&r0, sizeof(r0));
    CFX_MEMZERO_S(&r1, sizeof(r1));
    CFX_MEMZERO_S(&r2, sizeof(r2));
    CFX_MEMZERO_S(&r3, sizeof(r3));
    CFX_MEMZERO_S(&r4, sizeof(r4));

    CFX_MEMZERO_S(&s1, sizeof(s1));
    CFX_MEMZERO_S(&s2, sizeof(s2));
    CFX_MEMZERO_S(&s3, sizeof(s3));
    CFX_MEMZERO_S(&s4, sizeof(s4));

    CFX_MEMZERO_S(&h0, sizeof(h0));
    CFX_MEMZERO_S(&h1, sizeof(h1));
    CFX_MEMZERO_S(&h2, sizeof(h2));
    CFX_MEMZERO_S(&h3, sizeof(h3));
    CFX_MEMZERO_S(&h4, sizeof(h4));

    CFX_MEMZERO_S(&pad0, sizeof(pad0));
    CFX_MEMZERO_S(&pad1, sizeof(pad1));
    CFX_MEMZERO_S(&pad2, sizeof(pad2));
    CFX_MEMZERO_S(&pad3, sizeof(pad3));
    CFX_MEMZERO_S(&t0, sizeof(t0));
    CFX_MEMZERO_S(&t1, sizeof(t1));
    CFX_MEMZERO_S(&t2, sizeof(t2));
    CFX_MEMZERO_S(&t3, sizeof(t3));
    CFX_MEMZERO_S(&f,  sizeof(f));
}

#ifdef OLD_S
    #define S OLD_S
    #undef OLD_S
#endif
#ifdef OLD_X
    #define X OLD_X
    #undef OLD_X
#endif
