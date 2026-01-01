/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

/*
 * ge25519.c - Group element operations on Ed25519
 *
 * Extended coordinate formulas from Hisil-Wong-Carter-Dawson 2008.
 * Twisted Edwards curve: -x² + y² = 1 + d·x²·y² (mod p)
 *
 * Constant-time implementation: no branches on secret data.
 */

#include "cfx/ge25519.h"
#include <string.h>

/* d = -121665/121666 mod p */
const fe25519_t cfx_ge25519_d = {{
    0x00034dca135978a3ULL,
    0x0001a8283b156ebdULL,
    0x0005e7a26001c029ULL,
    0x000739c663a03cbbULL,
    0x00052036cee2b6ffULL
}};

/* 2*d */
const fe25519_t cfx_ge25519_2d = {{
    0x00069b9426b2f159ULL,
    0x00035050762add7aULL,
    0x0003cf44c0038052ULL,
    0x0006738cc7407977ULL,
    0x0002406d9dc56dffULL
}};

/* sqrt(-1) mod p */
const fe25519_t cfx_ge25519_sqrtm1 = {{
    0x00061b274a0ea0b0ULL,
    0x0000d5a5fc8f189dULL,
    0x0007ef5e9cbd0c60ULL,
    0x00078595a6804c9eULL,
    0x0002b8324804fc1dULL
}};

/*
 * Basepoint B:
 *   y = 4/5 mod p = 0x6666...6658
 *   x = positive sqrt satisfying the curve equation
 */
const ge25519_t cfx_ge25519_base = {
    /* X */
    {{
        0x00062d608f25d51aULL,
        0x000412a4b4f6592aULL,
        0x00075b7171a4b31dULL,
        0x0001ff60527118feULL,
        0x000216936d3cd6e5ULL
    }},
    /* Y = 4/5 mod p */
    {{
        0x0006666666666658ULL,
        0x0004ccccccccccccULL,
        0x0001999999999999ULL,
        0x0003333333333333ULL,
        0x0006666666666666ULL
    }},
    /* Z = 1 */
    {{
        0x0000000000000001ULL,
        0x0000000000000000ULL,
        0x0000000000000000ULL,
        0x0000000000000000ULL,
        0x0000000000000000ULL
    }},
    /* T = X*Y */
    {{
        0x00068ab3a5b7dda3ULL,
        0x00000eea2a5eadbbULL,
        0x0002af8df483c27eULL,
        0x000332b375274732ULL,
        0x00067875f0fd78b7ULL
    }}
};

void cfx_ge25519_0(ge25519_t* r) {
    cfx_fe25519_0(&r->X);
    cfx_fe25519_1(&r->Y);
    cfx_fe25519_1(&r->Z);
    cfx_fe25519_0(&r->T);
}

void cfx_ge25519_copy(ge25519_t* r, const ge25519_t* p) {
    cfx_fe25519_copy(&r->X, &p->X);
    cfx_fe25519_copy(&r->Y, &p->Y);
    cfx_fe25519_copy(&r->Z, &p->Z);
    cfx_fe25519_copy(&r->T, &p->T);
}

void cfx_ge25519_to_cached(ge25519_cached_t* r, const ge25519_t* p) {
    cfx_fe25519_add(&r->YplusX, &p->Y, &p->X);
    cfx_fe25519_sub(&r->YminusX, &p->Y, &p->X);
    cfx_fe25519_copy(&r->Z, &p->Z);
    cfx_fe25519_mul(&r->T2d, &p->T, &cfx_ge25519_2d);
}

/*
 * Point doubling using dedicated formula (faster than generic add).
 * Cost: 4M + 4S + 1D
 *
 * A = X1²
 * B = Y1²
 * C = 2·Z1²
 * D = -A  (a = -1 for our twisted Edwards curve)
 * E = (X1 + Y1)² - A - B
 * G = D + B
 * F = G - C
 * H = D - B
 * X3 = E·F
 * Y3 = G·H
 * Z3 = F·G
 * T3 = E·H
 */
void cfx_ge25519_double(ge25519_t* r, const ge25519_t* p) {
    fe25519_t A, B, C, D, E, F, G, H, t0;

    cfx_fe25519_sqr(&A, &p->X);              /* A = X² */
    cfx_fe25519_sqr(&B, &p->Y);              /* B = Y² */
    cfx_fe25519_sqr(&C, &p->Z);              /* C = Z² */
    cfx_fe25519_add(&C, &C, &C);             /* C = 2·Z² */
    cfx_fe25519_neg(&D, &A);                 /* D = -A (since a = -1) */
    cfx_fe25519_add(&t0, &p->X, &p->Y);      /* t0 = X + Y */
    cfx_fe25519_sqr(&E, &t0);                /* E = (X + Y)² */
    cfx_fe25519_sub(&E, &E, &A);             /* E = E - A */
    cfx_fe25519_sub(&E, &E, &B);             /* E = E - B */
    cfx_fe25519_add(&G, &D, &B);             /* G = D + B */
    cfx_fe25519_sub(&F, &G, &C);             /* F = G - C */
    cfx_fe25519_sub(&H, &D, &B);             /* H = D - B */
    cfx_fe25519_mul(&r->X, &E, &F);          /* X3 = E·F */
    cfx_fe25519_mul(&r->Y, &G, &H);          /* Y3 = G·H */
    cfx_fe25519_mul(&r->Z, &F, &G);          /* Z3 = F·G */
    cfx_fe25519_mul(&r->T, &E, &H);          /* T3 = E·H */
}

/*
 * Point addition with cached representation.
 * Cost: 8M + 1D
 *
 * Unified formula (works for doubling too, but cfx_ge25519_double is faster).
 */
void cfx_ge25519_add_cached(ge25519_t* r, const ge25519_t* p, const ge25519_cached_t* q) {
    fe25519_t A, B, C, D, E, F, G, H, t0;

    cfx_fe25519_add(&t0, &p->Y, &p->X);      /* t0 = Y1 + X1 */
    cfx_fe25519_sub(&A, &p->Y, &p->X);       /* A = Y1 - X1 */
    cfx_fe25519_mul(&A, &A, &q->YminusX);    /* A = (Y1 - X1)·(Y2 - X2) */
    cfx_fe25519_mul(&B, &t0, &q->YplusX);    /* B = (Y1 + X1)·(Y2 + X2) */
    cfx_fe25519_mul(&C, &p->T, &q->T2d);     /* C = T1·(2·d·T2) */
    cfx_fe25519_mul(&D, &p->Z, &q->Z);       /* D = Z1·Z2 */
    cfx_fe25519_add(&D, &D, &D);             /* D = 2·Z1·Z2 */
    cfx_fe25519_sub(&E, &B, &A);             /* E = B - A */
    cfx_fe25519_sub(&F, &D, &C);             /* F = D - C */
    cfx_fe25519_add(&G, &D, &C);             /* G = D + C */
    cfx_fe25519_add(&H, &B, &A);             /* H = B + A */
    cfx_fe25519_mul(&r->X, &E, &F);          /* X3 = E·F */
    cfx_fe25519_mul(&r->Y, &G, &H);          /* Y3 = G·H */
    cfx_fe25519_mul(&r->Z, &F, &G);          /* Z3 = F·G */
    cfx_fe25519_mul(&r->T, &E, &H);          /* T3 = E·H */
}

void cfx_ge25519_sub_cached(ge25519_t* r, const ge25519_t* p, const ge25519_cached_t* q) {
    fe25519_t A, B, C, D, E, F, G, H, t0;

    /* same as add but swap YplusX/YminusX and negate C */
    cfx_fe25519_add(&t0, &p->Y, &p->X);
    cfx_fe25519_sub(&A, &p->Y, &p->X);
    cfx_fe25519_mul(&A, &A, &q->YplusX);     /* swapped: use YplusX */
    cfx_fe25519_mul(&B, &t0, &q->YminusX);   /* swapped: use YminusX */
    cfx_fe25519_mul(&C, &p->T, &q->T2d);
    cfx_fe25519_mul(&D, &p->Z, &q->Z);
    cfx_fe25519_add(&D, &D, &D);
    cfx_fe25519_sub(&E, &B, &A);
    cfx_fe25519_add(&F, &D, &C);             /* swapped: D + C instead of D - C */
    cfx_fe25519_sub(&G, &D, &C);             /* swapped: D - C instead of D + C */
    cfx_fe25519_add(&H, &B, &A);
    cfx_fe25519_mul(&r->X, &E, &F);
    cfx_fe25519_mul(&r->Y, &G, &H);
    cfx_fe25519_mul(&r->Z, &F, &G);
    cfx_fe25519_mul(&r->T, &E, &H);
}

void cfx_ge25519_add_precomp(ge25519_t* r, const ge25519_t* p, const ge25519_precomp_t* q) {
    fe25519_t A, B, C, D, E, F, G, H, t0;

    cfx_fe25519_add(&t0, &p->Y, &p->X);
    cfx_fe25519_sub(&A, &p->Y, &p->X);
    cfx_fe25519_mul(&A, &A, &q->yminusx);
    cfx_fe25519_mul(&B, &t0, &q->yplusx);
    cfx_fe25519_mul(&C, &q->xy2d, &p->T);
    cfx_fe25519_add(&D, &p->Z, &p->Z);
    cfx_fe25519_sub(&E, &B, &A);
    cfx_fe25519_sub(&F, &D, &C);
    cfx_fe25519_add(&G, &D, &C);
    cfx_fe25519_add(&H, &B, &A);
    cfx_fe25519_mul(&r->X, &E, &F);
    cfx_fe25519_mul(&r->Y, &G, &H);
    cfx_fe25519_mul(&r->Z, &F, &G);
    cfx_fe25519_mul(&r->T, &E, &H);
}

void cfx_ge25519_sub_precomp(ge25519_t* r, const ge25519_t* p, const ge25519_precomp_t* q) {
    fe25519_t A, B, C, D, E, F, G, H, t0;

    cfx_fe25519_add(&t0, &p->Y, &p->X);
    cfx_fe25519_sub(&A, &p->Y, &p->X);
    cfx_fe25519_mul(&A, &A, &q->yplusx);
    cfx_fe25519_mul(&B, &t0, &q->yminusx);
    cfx_fe25519_mul(&C, &q->xy2d, &p->T);
    cfx_fe25519_add(&D, &p->Z, &p->Z);
    cfx_fe25519_sub(&E, &B, &A);
    cfx_fe25519_add(&F, &D, &C);
    cfx_fe25519_sub(&G, &D, &C);
    cfx_fe25519_add(&H, &B, &A);
    cfx_fe25519_mul(&r->X, &E, &F);
    cfx_fe25519_mul(&r->Y, &G, &H);
    cfx_fe25519_mul(&r->Z, &F, &G);
    cfx_fe25519_mul(&r->T, &E, &H);
}

void cfx_ge25519_neg(ge25519_t* r, const ge25519_t* p) {
    cfx_fe25519_neg(&r->X, &p->X);
    cfx_fe25519_copy(&r->Y, &p->Y);
    cfx_fe25519_copy(&r->Z, &p->Z);
    cfx_fe25519_neg(&r->T, &p->T);
}

void cfx_ge25519_cmov(ge25519_t* r, const ge25519_t* p, unsigned int b) {
    cfx_fe25519_cmov(&r->X, &p->X, b);
    cfx_fe25519_cmov(&r->Y, &p->Y, b);
    cfx_fe25519_cmov(&r->Z, &p->Z, b);
    cfx_fe25519_cmov(&r->T, &p->T, b);
}

void cfx_ge25519_cneg(ge25519_t* p, unsigned int b) {
    fe25519_t negX, negT;
    cfx_fe25519_neg(&negX, &p->X);
    cfx_fe25519_neg(&negT, &p->T);
    cfx_fe25519_cmov(&p->X, &negX, b);
    cfx_fe25519_cmov(&p->T, &negT, b);
}

/*
 * Scalar multiplication: r = [s]p
 * Uses constant-time double-and-add with Montgomery ladder style.
 */
void cfx_ge25519_scalarmult(ge25519_t* r, const uint8_t s[32], const ge25519_t* p) {
    ge25519_t R0, R1;
    ge25519_cached_t p_cached;
    int i;

    cfx_ge25519_0(&R0);                      /* R0 = identity */
    cfx_ge25519_copy(&R1, p);                /* R1 = P */
    cfx_ge25519_to_cached(&p_cached, p);

    /* process bits from most significant to least */
    for (i = 254; i >= 0; i--) {
        unsigned int bit = (s[i >> 3] >> (i & 7)) & 1;
        ge25519_cached_t R0_cached, R1_cached;

        /* constant-time swap based on bit */
        cfx_ge25519_to_cached(&R0_cached, &R0);
        cfx_ge25519_to_cached(&R1_cached, &R1);

        /* if bit: R0 = R0 + R1, R1 = 2*R1 */
        /* else:   R0 = 2*R0, R1 = R0 + R1 */
        ge25519_t R0_new, R1_new;

        cfx_ge25519_double(&R0_new, &R0);
        cfx_ge25519_add_cached(&R1_new, &R0, &R1_cached);

        ge25519_t R0_alt, R1_alt;
        cfx_ge25519_add_cached(&R0_alt, &R0, &R1_cached);
        cfx_ge25519_double(&R1_alt, &R1);

        /* select based on bit */
        cfx_ge25519_cmov(&R0_new, &R0_alt, bit);
        cfx_ge25519_cmov(&R1_new, &R1_alt, bit);

        cfx_ge25519_copy(&R0, &R0_new);
        cfx_ge25519_copy(&R1, &R1_new);
    }

    cfx_ge25519_copy(r, &R0);
}

/*
 * Scalar multiplication with basepoint: r = [s]B
 * For now, same as cfx_ge25519_scalarmult with basepoint.
 * TODO: precomputed table for faster basepoint multiplication.
 */
void cfx_ge25519_scalarmult_base(ge25519_t* r, const uint8_t s[32]) {
    cfx_ge25519_scalarmult(r, s, &cfx_ge25519_base);
}

/*
 * Pack point to 32 bytes.
 * Encoding: y coordinate with sign of x in the top bit.
 */
void cfx_ge25519_pack(uint8_t out[32], const ge25519_t* p) {
    fe25519_t x, y, z_inv;

    /* convert to affine: x = X/Z, y = Y/Z */
    cfx_fe25519_inv(&z_inv, &p->Z);
    cfx_fe25519_mul(&x, &p->X, &z_inv);
    cfx_fe25519_mul(&y, &p->Y, &z_inv);

    /* encode y */
    cfx_fe25519_tobytes(out, &y);

    /* set top bit to sign of x */
    out[31] |= (uint8_t)(cfx_fe25519_isnegative(&x) << 7);
}

/*
 * Unpack 32 bytes to point.
 * Returns 0 on success, -1 if not on curve.
 */
int cfx_ge25519_unpack(ge25519_t* p, const uint8_t in[32]) {
    fe25519_t u, v, v3, vxx, check;
    int x_sign;

    /* extract y, clearing top bit */
    uint8_t y_bytes[32];
    memcpy(y_bytes, in, 32);
    x_sign = (y_bytes[31] >> 7) & 1;
    y_bytes[31] &= 0x7f;

    cfx_fe25519_frombytes(&p->Y, y_bytes);

    /* compute x² = (y² - 1) / (d·y² + 1) */
    cfx_fe25519_sqr(&u, &p->Y);              /* u = y² */
    cfx_fe25519_mul(&v, &u, &cfx_ge25519_d);     /* v = d·y² */
    cfx_fe25519_sub(&u, &u, &cfx_fe25519_one);   /* u = y² - 1 */
    cfx_fe25519_add(&v, &v, &cfx_fe25519_one);   /* v = d·y² + 1 */

    /* compute v³ for the sqrt computation */
    cfx_fe25519_sqr(&v3, &v);                /* v² */
    cfx_fe25519_mul(&v3, &v3, &v);           /* v³ */

    /* x = u·v³·(u·v⁷)^((p-5)/8) */
    cfx_fe25519_sqr(&p->X, &v3);             /* v⁶ */
    cfx_fe25519_mul(&p->X, &p->X, &v);       /* v⁷ */
    cfx_fe25519_mul(&p->X, &p->X, &u);       /* u·v⁷ */
    cfx_fe25519_pow22523(&p->X, &p->X);      /* (u·v⁷)^((p-5)/8) */
    cfx_fe25519_mul(&p->X, &p->X, &v3);      /* v³·(u·v⁷)^((p-5)/8) */
    cfx_fe25519_mul(&p->X, &p->X, &u);       /* u·v³·(u·v⁷)^((p-5)/8) */

    /* check: v·x² == u? */
    cfx_fe25519_sqr(&vxx, &p->X);
    cfx_fe25519_mul(&vxx, &vxx, &v);
    cfx_fe25519_sub(&check, &vxx, &u);

    if (cfx_fe25519_iszero(&check) == 0) {
        /* try x = x·sqrt(-1) */
        cfx_fe25519_mul(&p->X, &p->X, &cfx_ge25519_sqrtm1);
        cfx_fe25519_sqr(&vxx, &p->X);
        cfx_fe25519_mul(&vxx, &vxx, &v);
        cfx_fe25519_sub(&check, &vxx, &u);
        if (cfx_fe25519_iszero(&check) == 0) {
            return -1;  /* not on curve */
        }
    }

    /* fix sign of x */
    if (cfx_fe25519_isnegative(&p->X) != x_sign) {
        cfx_fe25519_neg(&p->X, &p->X);
    }

    /* set Z = 1 and T = X*Y */
    cfx_fe25519_1(&p->Z);
    cfx_fe25519_mul(&p->T, &p->X, &p->Y);

    return 0;
}

int cfx_ge25519_is_on_curve(const ge25519_t* p) {
    fe25519_t x2, y2, z2, z4, lhs, rhs, tmp;

    /* compute in projective: check -X²Z² + Y²Z² = Z⁴ + d·X²Y² */
    cfx_fe25519_sqr(&x2, &p->X);
    cfx_fe25519_sqr(&y2, &p->Y);
    cfx_fe25519_sqr(&z2, &p->Z);
    cfx_fe25519_sqr(&z4, &z2);

    /* lhs = -X²·Z² + Y²·Z² = Z²·(Y² - X²) */
    cfx_fe25519_sub(&lhs, &y2, &x2);
    cfx_fe25519_mul(&lhs, &lhs, &z2);

    /* rhs = Z⁴ + d·X²·Y² */
    cfx_fe25519_mul(&tmp, &x2, &y2);
    cfx_fe25519_mul(&rhs, &tmp, &cfx_ge25519_d);
    cfx_fe25519_add(&rhs, &rhs, &z4);

    cfx_fe25519_sub(&tmp, &lhs, &rhs);
    return cfx_fe25519_iszero(&tmp);
}

int cfx_ge25519_is_identity(const ge25519_t* p) {
    /* identity is (0 : c : c : 0) for any c, or in extended: X = 0, Y = Z */
    fe25519_t tmp;
    int x_zero = cfx_fe25519_iszero(&p->X);
    cfx_fe25519_sub(&tmp, &p->Y, &p->Z);
    int y_eq_z = cfx_fe25519_iszero(&tmp);
    return x_zero && y_eq_z;
}
