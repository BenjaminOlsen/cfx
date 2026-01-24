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
#include "cfx/memory.h"
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

#ifdef CFX_ED25519_PARANOID
/*
 * Multiply point by cofactor (8 = 2^3).
 * Used for paranoid small-subgroup and cofactored verification checks.
 */
static void ge25519_mul_cofactor(ge25519_t *r, const ge25519_t *p) {
    cfx_ge25519_double(r, p);   /* 2P */
    cfx_ge25519_double(r, r);   /* 4P */
    cfx_ge25519_double(r, r);   /* 8P */
}
#endif

/*
 * Check if scalar s is in canonical form (s < L).
 * L = 2^252 + 27742317777372353535851937790883648493
 *   = 0x1000000000000000000000000000000014def9dea2f79cd65812631a5cf5d3ed
 * Returns 1 if s < L (valid), 0 otherwise.
 * Constant-time: computes s - L and checks for borrow.
 */

/*************************************
simple example in base 10:
........................................
ex 1:
 s = 100 L = 99
-> s = {0, 0, 1}, L = {9, 9, 0}

i = 0:
diff = 0 - 9 - 0 = -9 -> borrow = 1
i = 1
diff = 0 -9 -1 = -10 -> borrow = 1;
i = 2
diff = 1 - 0 - 1 = 0 -> borrow = 0; -> L < s

........................................
ex 2:
s = 99, L = 100:
i = 0
diff = 9 - 0 - 0 -> borrow = 0;
diff = 9 - 0 - 0 -> borrow = 0;
diff = 1 - 0 - 0 -> borrow = 1; -> s < L
*/
static int sc25519_is_canonical(const uint8_t s[32]) {
    static const uint8_t L[32] = {
        0xed, 0xd3, 0xf5, 0x5c, 0x1a, 0x63, 0x12, 0x58,
        0xd6, 0x9c, 0xf7, 0xa2, 0xde, 0xf9, 0xde, 0x14,
        0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x10
    };

    int borrow = 0;
    for (int i = 0; i < 32; i++) {
        int diff = (int)s[i] - (int)L[i] - borrow;
        borrow = (diff >> 8) & 1;
    }
    return borrow;  /* 1 if s < L, 0 otherwise */
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

    CFX_MEMZERO_S(hash, sizeof(hash));
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

    CFX_MEMZERO_S(hash, sizeof(hash));
    CFX_MEMZERO_S(r_scalar, sizeof(r_scalar));
    CFX_MEMZERO_S(k_scalar, sizeof(k_scalar));
    CFX_MEMZERO_S(a, sizeof(a));
    CFX_MEMZERO_S(r, sizeof(r));
    CFX_MEMZERO_S(k, sizeof(k));
}

int cfx_ed25519_verify(const uint8_t sig[64], const uint8_t* msg, size_t msg_len, const uint8_t pk[32]) {
    ge25519_t A, R, sB, kA, check;
    ge25519_cached_t kA_cached;
    uint8_t k_scalar[64];
    uint8_t k[32];
    uint8_t check_bytes[32];
    uint8_t expected_bytes[32];
    cfx_sha512_ctx_t ctx;

    /* unpack public key */
    if (cfx_ge25519_unpack(&A, pk) != 0) {
        return -1;
    }

    /* unpack R from signature */
    if (cfx_ge25519_unpack(&R, sig) != 0) {
        return -1;
    }

    /* check s < L (group order) - reject non-canonical signatures */
    if (!sc25519_is_canonical(sig + 32)) {
        return -1;
    }

#ifdef CFX_ED25519_PARANOID
    /* Paranoid: reject small-order public keys.
     * If [8]A == identity, then A is in the small subgroup (order dividing 8).
     * Legitimate keys are A = [a]G which always have order L. */
    {
        ge25519_t A8;
        ge25519_mul_cofactor(&A8, &A);
        if (cfx_ge25519_is_identity(&A8)) {
            return -1;
        }
    }
#endif

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
    cfx_ge25519_to_cached(&kA_cached, &kA);
    cfx_ge25519_add_cached(&check, &R, &kA_cached);

#ifdef CFX_ED25519_PARANOID
    /* Paranoid: cofactored verification.
     * Check [8][s]B == [8](R + [k]A) instead of [s]B == R + [k]A.
     * This clears any small-order components that might allow forgery. */
    ge25519_mul_cofactor(&sB, &sB);
    ge25519_mul_cofactor(&check, &check);
#endif

    /* compare [s]B with R + [k]A (or their cofactor multiples) */
    cfx_ge25519_pack(check_bytes, &sB);
    cfx_ge25519_pack(expected_bytes, &check);

    uint8_t diff = 0;
    for (int i = 0; i < 32; i++) {
        diff |= check_bytes[i] ^ expected_bytes[i];
    }

    return diff == 0 ? 0 : -1;
}

/*
 * Ed25519ph - pre-hashed variant (RFC 8032 section 5.1)
 *
 * dom2(x, y) = "SigEd25519 no Ed25519 collisions" || octet(x) || octet(len(y)) || y
 *
 * For Ed25519ph with empty context:
 *   x = 1 (pre-hash flag)
 *   y = "" (empty context)
 *   dom2(1, "") = "SigEd25519 no Ed25519 collisions" || 0x01 || 0x00
 */
static const uint8_t DOM2_PH_PREFIX[34] = {
    /* "SigEd25519 no Ed25519 collisions" (32 bytes) */
    'S', 'i', 'g', 'E', 'd', '2', '5', '5', '1', '9', ' ',
    'n', 'o', ' ',
    'E', 'd', '2', '5', '5', '1', '9', ' ',
    'c', 'o', 'l', 'l', 'i', 's', 'i', 'o', 'n', 's',
    0x01,  /* phflag = 1 for pre-hash */
    0x00   /* context length = 0 */
};

void cfx_ed25519ph_sign(uint8_t sig[64], const uint8_t prehash[64], const uint8_t sk[64]) {
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

    /* r = SHA-512(dom2(1,"") || prefix || PH(M)) mod L */
    cfx_sha512_init(&ctx);
    cfx_sha512_update(&ctx, DOM2_PH_PREFIX, sizeof(DOM2_PH_PREFIX));
    cfx_sha512_update(&ctx, hash + 32, 32);  /* prefix */
    cfx_sha512_update(&ctx, prehash, 64);    /* PH(M) */
    cfx_sha512_final(&ctx, r_scalar);

    uint8_t r[32];
    cfx_sc25519_reduce(r, r_scalar);

    /* R = [r]B */
    cfx_ge25519_scalarmult_base(&R, r);
    cfx_ge25519_pack(sig, &R);

    /* k = SHA-512(dom2(1,"") || R || A || PH(M)) mod L */
    cfx_sha512_init(&ctx);
    cfx_sha512_update(&ctx, DOM2_PH_PREFIX, sizeof(DOM2_PH_PREFIX));
    cfx_sha512_update(&ctx, sig, 32);        /* R */
    cfx_sha512_update(&ctx, sk + 32, 32);    /* A (public key) */
    cfx_sha512_update(&ctx, prehash, 64);    /* PH(M) */
    cfx_sha512_final(&ctx, k_scalar);

    uint8_t k[32];
    cfx_sc25519_reduce(k, k_scalar);

    cfx_sc25519_muladd(sig + 32, k, a, r);  /* s = (r + k*a) mod L */

    CFX_MEMZERO_S(hash, sizeof(hash));
    CFX_MEMZERO_S(r_scalar, sizeof(r_scalar));
    CFX_MEMZERO_S(k_scalar, sizeof(k_scalar));
    CFX_MEMZERO_S(a, sizeof(a));
    CFX_MEMZERO_S(r, sizeof(r));
    CFX_MEMZERO_S(k, sizeof(k));
}

int cfx_ed25519ph_verify(const uint8_t sig[64], const uint8_t prehash[64], const uint8_t pk[32]) {
    ge25519_t A, R, sB, kA, check;
    ge25519_cached_t kA_cached;
    uint8_t k_scalar[64];
    uint8_t k[32];
    uint8_t check_bytes[32];
    uint8_t expected_bytes[32];
    cfx_sha512_ctx_t ctx;

    /* unpack public key */
    if (cfx_ge25519_unpack(&A, pk) != 0) {
        return -1;
    }

    /* unpack R from signature */
    if (cfx_ge25519_unpack(&R, sig) != 0) {
        return -1;
    }

    /* check s < L (group order) - reject non-canonical signatures */
    if (!sc25519_is_canonical(sig + 32)) {
        return -1;
    }

#ifdef CFX_ED25519_PARANOID
    /* Paranoid: reject small-order public keys */
    {
        ge25519_t A8;
        ge25519_mul_cofactor(&A8, &A);
        if (cfx_ge25519_is_identity(&A8)) {
            return -1;
        }
    }
#endif

    /* k = SHA-512(dom2(1,"") || R || A || PH(M)) mod L */
    cfx_sha512_init(&ctx);
    cfx_sha512_update(&ctx, DOM2_PH_PREFIX, sizeof(DOM2_PH_PREFIX));
    cfx_sha512_update(&ctx, sig, 32);      /* R */
    cfx_sha512_update(&ctx, pk, 32);       /* A */
    cfx_sha512_update(&ctx, prehash, 64);  /* PH(M) */
    cfx_sha512_final(&ctx, k_scalar);
    cfx_sc25519_reduce(k, k_scalar);

    /* verify: [s]B == R + [k]A */
    cfx_ge25519_scalarmult_base(&sB, sig + 32);
    cfx_ge25519_scalarmult(&kA, k, &A);
    cfx_ge25519_to_cached(&kA_cached, &kA);
    cfx_ge25519_add_cached(&check, &R, &kA_cached);

#ifdef CFX_ED25519_PARANOID
    /* Paranoid: cofactored verification */
    ge25519_mul_cofactor(&sB, &sB);
    ge25519_mul_cofactor(&check, &check);
#endif

    /* compare [s]B with R + [k]A (or their cofactor multiples) */
    cfx_ge25519_pack(check_bytes, &sB);
    cfx_ge25519_pack(expected_bytes, &check);

    uint8_t diff = 0;
    for (int i = 0; i < 32; i++) {
        diff |= check_bytes[i] ^ expected_bytes[i];
    }

    return diff == 0 ? 0 : -1;
}
