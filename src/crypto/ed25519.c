/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

/*
 * ed25519.c - Ed25519 digital signatures (RFC 8032)
 *
 * Deterministic Schnorr signatures on the twisted Edwards curve.
 * Uses SHA-512 for hashing, ge25519 for point operations.
 */

#include "cfx/ed25519.h"
#include "cfx/ge25519.h"
#include "cfx/sc25519.h"
#include "cfx/sha512.h"
#include <string.h>

/*
 * Clamp a 32-byte scalar for Ed25519 key derivation.
 * Clear bottom 3 bits (cofactor clearing)
 * Clear top bit (bit 255)
 * Set bit 254 (ensure fixed position)
 */
static void clamp_scalar(uint8_t s[32]) {
    s[0] &= 248;
    s[31] &= 127;
    s[31] |= 64;
}

void cfx_ed25519_create_keypair(uint8_t pk[32], uint8_t sk[64], const uint8_t seed[32]) {
    uint8_t hash[64];
    ge25519_t A;

    /* hash seed to get expanded secret key */
    cfx_sha512(hash, seed, 32);

    /* clamp first half for scalar */
    clamp_scalar(hash);

    /* compute public key A = [a]B */
    cfx_ge25519_scalarmult_base(&A, hash);
    cfx_ge25519_pack(pk, &A);

    /* secret key = seed || public_key */
    memcpy(sk, seed, 32);
    memcpy(sk + 32, pk, 32);
}

void cfx_ed25519_get_public_key(uint8_t pk[32], const uint8_t sk[64]) {
    memcpy(pk, sk + 32, 32);
}

void cfx_ed25519_sign(uint8_t sig[64], const uint8_t* msg, size_t msg_len, const uint8_t sk[64]) {
    uint8_t hash[64];
    uint8_t r_scalar[64];
    uint8_t k_scalar[64];
    uint8_t a[32];
    ge25519_t R;
    cfx_sha512_ctx_t ctx;

    /* hash(seed) -> (scalar a, prefix) */
    cfx_sha512(hash, sk, 32);
    memcpy(a, hash, 32);
    clamp_scalar(a);

    /* r = SHA-512(prefix || msg) mod L (deterministic nonce) */
    cfx_sha512_init(&ctx);
    cfx_sha512_update(&ctx, hash + 32, 32);  /* prefix = second half of hash */
    cfx_sha512_update(&ctx, msg, msg_len);
    cfx_sha512_final(&ctx, r_scalar);

    /* reduce r mod L */
    uint8_t r[32];
    cfx_sc25519_reduce(r, r_scalar);

    /* R = [r]B */
    cfx_ge25519_scalarmult_base(&R, r);
    cfx_ge25519_pack(sig, &R);  /* first 32 bytes of sig = R */

    /* k = SHA-512(R || A || msg) mod L */
    cfx_sha512_init(&ctx);
    cfx_sha512_update(&ctx, sig, 32);       /* R */
    cfx_sha512_update(&ctx, sk + 32, 32);   /* A (public key) */
    cfx_sha512_update(&ctx, msg, msg_len);
    cfx_sha512_final(&ctx, k_scalar);

    uint8_t k[32];
    cfx_sc25519_reduce(k, k_scalar);

    /* s = (r + k*a) mod L */
    cfx_sc25519_muladd(sig + 32, k, a, r);  /* second 32 bytes of sig = s */
}

int cfx_ed25519_verify(const uint8_t sig[64], const uint8_t* msg, size_t msg_len, const uint8_t pk[32]) {
    ge25519_t A, R, sB, kA, check;
    ge25519_cached_t R_cached, kA_cached;
    uint8_t k_scalar[64];
    uint8_t k[32];
    uint8_t check_bytes[32];
    cfx_sha512_ctx_t ctx;

    /* unpack public key */
    if (cfx_ge25519_unpack(&A, pk) != 0) {
        return -1;
    }

    /* unpack R from signature */
    if (cfx_ge25519_unpack(&R, sig) != 0) {
        return -1;
    }

    /* check s < L */
    /* simplified check: top byte should be < 0x10 for valid s */
    if (sig[63] & 0xf0) {
        return -1;
    }

    /* k = SHA-512(R || A || msg) mod L */
    cfx_sha512_init(&ctx);
    cfx_sha512_update(&ctx, sig, 32);     /* R */
    cfx_sha512_update(&ctx, pk, 32);      /* A */
    cfx_sha512_update(&ctx, msg, msg_len);
    cfx_sha512_final(&ctx, k_scalar);
    cfx_sc25519_reduce(k, k_scalar);

    /* verify: [s]B == R + [k]A */
    /* compute [s]B */
    cfx_ge25519_scalarmult_base(&sB, sig + 32);

    /* compute [k]A */
    cfx_ge25519_scalarmult(&kA, k, &A);

    /* compute R + [k]A */
    cfx_ge25519_to_cached(&R_cached, &R);
    cfx_ge25519_to_cached(&kA_cached, &kA);
    cfx_ge25519_add_cached(&check, &R, &kA_cached);

    /* compare [s]B with R + [k]A */
    cfx_ge25519_pack(check_bytes, &sB);
    uint8_t expected_bytes[32];
    cfx_ge25519_pack(expected_bytes, &check);

    /* constant-time comparison */
    uint8_t diff = 0;
    for (int i = 0; i < 32; i++) {
        diff |= check_bytes[i] ^ expected_bytes[i];
    }

    return diff == 0 ? 0 : -1;
}
