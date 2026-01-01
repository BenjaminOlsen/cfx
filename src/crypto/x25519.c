/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

/*
 * x25519.c - X25519 Diffie-Hellman (RFC 7748)
 *
 * montgomery ladder on curve25519, x-coordinate only.
 * constant-time: no branches or memory access patterns based on secret data.
 */

#include <cfx/x25519.h>
#include <cfx/fe25519.h>
#include <string.h>

/* basepoint x = 9 (little-endian) */
const uint8_t cfx_x25519_basepoint[32] = {
    9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};

/*
 * scalar clamping (RFC 7748 section 5):
 * - clear bottom 3 bits (multiple of 8 for cofactor)
 * - clear bit 255 (always < 2^255)
 * - set bit 254 (fixed position for timing)
 */
static void clamp_scalar(uint8_t e[32]) {
    e[0]  &= 0xF8;   /* clear bits 0,1,2 */
    e[31] &= 0x7F;   /* clear bit 255 */
    e[31] |= 0x40;   /* set bit 254 */
}

/*
 * montgomery ladder step
 *
 * given projective points (x2:z2) and (x3:z3) where x3 - x2 = base,
 * computes (x2:z2) = 2*(x2:z2) and (x3:z3) = (x2:z2) + (x3:z3)
 *
 * formulas from https://cr.yp.to/ecdh.html
 * a24 = (A+2)/4 = 121666 where A = 486662
 */
static void ladder_step(fe25519_t* x2, fe25519_t* z2,
                        fe25519_t* x3, fe25519_t* z3,
                        const fe25519_t* x1) {
    fe25519_t a, aa, b, bb, c, d, da, cb, e;

    /* doubling side */
    cfx_fe25519_add(&a, x2, z2);      /* a = x2 + z2 */
    cfx_fe25519_sqr(&aa, &a);         /* aa = a² */
    cfx_fe25519_sub(&b, x2, z2);      /* b = x2 - z2 */
    cfx_fe25519_sqr(&bb, &b);         /* bb = b² */
    cfx_fe25519_sub(&e, &aa, &bb);    /* e = aa - bb */

    /* differential addition side */
    cfx_fe25519_add(&c, x3, z3);      /* c = x3 + z3 */
    cfx_fe25519_sub(&d, x3, z3);      /* d = x3 - z3 */
    cfx_fe25519_mul(&da, &d, &a);     /* da = d * a */
    cfx_fe25519_mul(&cb, &c, &b);     /* cb = c * b */

    /* new x3, z3 (addition result) */
    cfx_fe25519_add(x3, &da, &cb);    /* x3 = da + cb */
    cfx_fe25519_sqr(x3, x3);          /* x3 = (da + cb)² */
    cfx_fe25519_sub(z3, &da, &cb);    /* z3 = da - cb */
    cfx_fe25519_sqr(z3, z3);          /* z3 = (da - cb)² */
    cfx_fe25519_mul(z3, z3, x1);      /* z3 = x1 * (da - cb)² */

    /* new x2, z2 (doubling result) */
    cfx_fe25519_mul(x2, &aa, &bb);    /* x2 = aa * bb */
    cfx_fe25519_mul121666(z2, &e);    /* z2 = 121666 * e */
    cfx_fe25519_add(z2, z2, &bb);     /* z2 = bb + 121666*e */
    cfx_fe25519_mul(z2, z2, &e);      /* z2 = e * (bb + 121666*e) */
}

/*
 * montgomery ladder for scalar multiplication
 *
 * computes [scalar] * point using constant-time x-only ladder.
 * processes bits top-down, maintaining x2 = [n]*P and x3 = [n+1]*P.
 */
static void ladder(fe25519_t* result, const uint8_t scalar[32], const fe25519_t* point) {
    fe25519_t x1, x2, z2, x3, z3;
    int i, bit, swap, prev_swap;

    cfx_fe25519_copy(&x1, point);     /* x1 = point (the difference, never changes) */
    cfx_fe25519_1(&x2);               /* x2 = 1 (projective for [0]*P, really infinity) */
    cfx_fe25519_0(&z2);               /* z2 = 0 */
    cfx_fe25519_copy(&x3, point);     /* x3 = point ([1]*P) */
    cfx_fe25519_1(&z3);               /* z3 = 1 */

    prev_swap = 0;

    /* process 255 bits, top down (bit 254 is highest due to clamping) */
    for (i = 254; i >= 0; i--) {
        bit = (scalar[i >> 3] >> (i & 7)) & 1;

        /* constant-time conditional swap based on bit XOR previous bit */
        swap = bit ^ prev_swap;
        cfx_fe25519_cswap(&x2, &x3, swap);
        cfx_fe25519_cswap(&z2, &z3, swap);
        prev_swap = bit;

        ladder_step(&x2, &z2, &x3, &z3, &x1);
    }

    /* final swap to ensure x2 has the result */
    cfx_fe25519_cswap(&x2, &x3, prev_swap);
    cfx_fe25519_cswap(&z2, &z3, prev_swap);

    /* result = x2 / z2 = x2 * z2^(-1) */
    cfx_fe25519_inv(&z2, &z2);
    cfx_fe25519_mul(result, &x2, &z2);
}

/* raw scalar multiplication (no clamping) */
void cfx_x25519_scalarmult(uint8_t out[32], const uint8_t scalar[32], const uint8_t point[32]) {
    fe25519_t p, result;

    cfx_fe25519_frombytes(&p, point);
    ladder(&result, scalar, &p);
    cfx_fe25519_tobytes(out, &result);
}

/* compute public key from private key */
void cfx_x25519_base(uint8_t public_key[32], const uint8_t private_key[32]) {
    uint8_t e[32];

    memcpy(e, private_key, 32);
    clamp_scalar(e);
    cfx_x25519_scalarmult(public_key, e, cfx_x25519_basepoint);
}

/* full DH function with clamping and zero check */
int cfx_x25519(uint8_t shared[32], const uint8_t scalar[32], const uint8_t point[32]) {
    uint8_t e[32];
    int i;
    unsigned int zero_check;

    memcpy(e, scalar, 32);
    clamp_scalar(e);
    cfx_x25519_scalarmult(shared, e, point);

    /* constant-time check for all-zero result (indicates bad input) */
    zero_check = 0;
    for (i = 0; i < 32; i++) {
        zero_check |= shared[i];
    }

    /* return -1 if result is all zeros, 0 otherwise */
    return (int)((zero_check - 1) >> 31) * -1;
}
