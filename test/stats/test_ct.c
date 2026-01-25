/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */
/*
 * test_ct.c - Constant-time testing using Welch's t-test
 *
 * Compares timing distributions of two input classes:
 *   Class 0: Fixed/predictable inputs
 *   Class 1: Random inputs
 *
 * If |t| > threshold, the implementation likely has timing leakage.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#include "cfx/ed25519.h"
#include "cfx/x25519.h"
#include "cfx/poly1305.h"
#include "cfx/chacha20.h"
#include "cfx/sha512.h"
#include "cfx/rand.h"

#ifdef _WIN32
#include <windows.h>
#else
#include <time.h>
#endif

/* Test parameters */
#define WARMUP_ITERATIONS    1000
#define DEFAULT_MEASUREMENTS 500000
#define THRESHOLD            4.5

static int verbose = 1;
static int fail_fast = 1;
static int num_measurements = DEFAULT_MEASUREMENTS;

/* Get a random bit for class selection */
static int rand_bit(void) {
    uint8_t b;
    cfx_rand_bytes(&b, 1);
    return b & 1;
}

static uint64_t get_time(void) {
#ifdef _WIN32
    LARGE_INTEGER t;
    QueryPerformanceCounter(&t);
    return (uint64_t)t.QuadPart;
#elif defined(__APPLE__)
    return clock_gettime_nsec_np(CLOCK_MONOTONIC);
#else
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (uint64_t)ts.tv_sec * 1000000000ULL + ts.tv_nsec;
#endif
}

/* Welford's online algorithm for mean/variance */
typedef struct {
    double n, mean, m2;
} ttest_ctx;

static void ttest_init(ttest_ctx *ctx) {
    ctx->n = ctx->mean = ctx->m2 = 0;
}

static void ttest_push(ttest_ctx *ctx, double x) {
    ctx->n++;
    double delta = x - ctx->mean;
    ctx->mean += delta / ctx->n;
    ctx->m2 += delta * (x - ctx->mean);
}

static double ttest_compute(ttest_ctx *a, ttest_ctx *b) {
    if (a->n < 2 || b->n < 2) return 0;
    double var_a = a->m2 / (a->n - 1);
    double var_b = b->m2 / (b->n - 1);
    double se = sqrt(var_a / a->n + var_b / b->n);
    if (se < 1e-10) return 0;
    return (a->mean - b->mean) / se;
}

/* Prevent compiler from optimizing away results */
static volatile uint8_t sink;

/*============================================================================
 * Ed25519 tests
 *============================================================================*/

static int test_ed25519_sign(void) {
    printf("Testing cfx_ed25519_sign()...\n");

    ttest_ctx class0, class1;
    ttest_init(&class0);
    ttest_init(&class1);

    /* Fixed inputs for class 0 */
    uint8_t seed_fixed[32], msg_fixed[64];
    memset(seed_fixed, 0x42, 32);
    memset(msg_fixed, 0x00, 64);

    uint8_t seed[32], pk[32], sk[64], msg[64], sig[64];

    /* Warmup */
    for (int i = 0; i < WARMUP_ITERATIONS; i++) {
        cfx_rand_bytes(seed, 32);
        cfx_ed25519_create_keypair(pk, sk, seed);
        cfx_ed25519_sign(sig, msg_fixed, 64, sk);
    }

    for (int i = 0; i < num_measurements; i++) {
        int class_id = rand_bit();

        if (class_id) {
            cfx_rand_bytes(seed, 32);
            cfx_rand_bytes(msg, 64);
        } else {
            memcpy(seed, seed_fixed, 32);
            memcpy(msg, msg_fixed, 64);
        }

        cfx_ed25519_create_keypair(pk, sk, seed);

        uint64_t t0 = get_time();
        cfx_ed25519_sign(sig, msg, 64, sk);
        uint64_t t1 = get_time();

        sink = sig[0];

        ttest_push(class_id ? &class1 : &class0, (double)(t1 - t0));

        if (i > (num_measurements / 10) && (i % (num_measurements / 10)) == 0) {
            double t = ttest_compute(&class0, &class1);
            if (verbose) printf("  %d: t = %.3f\n", i, t);
            if (fail_fast && fabs(t) > THRESHOLD) {
                printf("  FAIL (timing leakage detected)\n");
                return 1;
            }
        }
    }

    double t = ttest_compute(&class0, &class1);
    printf("  final: t = %.3f %s\n", t, fabs(t) > THRESHOLD ? "FAIL" : "PASS");
    return fabs(t) > THRESHOLD ? 1 : 0;
}

static int test_ed25519_verify(void) {
    printf("Testing cfx_ed25519_verify() (valid vs invalid sig)...\n");

    ttest_ctx class0, class1;
    ttest_init(&class0);
    ttest_init(&class1);

    /* Generate a keypair and valid signature */
    uint8_t seed[32], pk[32], sk[64], msg[64], sig_valid[64], sig_invalid[64];
    cfx_rand_bytes(seed, 32);
    cfx_rand_bytes(msg, 64);
    cfx_ed25519_create_keypair(pk, sk, seed);
    cfx_ed25519_sign(sig_valid, msg, 64, sk);
    memcpy(sig_invalid, sig_valid, 64);
    sig_invalid[0] ^= 0x01;  /* Flip one bit */

    uint8_t sig_test[64];

    /* Warmup */
    for (int i = 0; i < WARMUP_ITERATIONS; i++) {
        cfx_ed25519_verify(sig_valid, msg, 64, pk);
    }

    for (int i = 0; i < num_measurements; i++) {
        int class_id = rand_bit();
        memcpy(sig_test, class_id ? sig_invalid : sig_valid, 64);

        uint64_t t0 = get_time();
        volatile int r = cfx_ed25519_verify(sig_test, msg, 64, pk);
        uint64_t t1 = get_time();
        (void)r;

        ttest_push(class_id ? &class1 : &class0, (double)(t1 - t0));

        if (i > (num_measurements / 10) && (i % (num_measurements / 10)) == 0) {
            double t = ttest_compute(&class0, &class1);
            if (verbose) printf("  %d: t = %.3f\n", i, t);
            if (fail_fast && fabs(t) > THRESHOLD) {
                printf("  FAIL (timing leakage detected)\n");
                return 1;
            }
        }
    }

    double t = ttest_compute(&class0, &class1);
    printf("  final: t = %.3f %s\n", t, fabs(t) > THRESHOLD ? "FAIL" : "PASS");
    return fabs(t) > THRESHOLD ? 1 : 0;
}

/*============================================================================
 * X25519 tests
 *============================================================================*/

static int test_x25519(void) {
    printf("Testing cfx_x25519()...\n");

    ttest_ctx class0, class1;
    ttest_init(&class0);
    ttest_init(&class1);

    /* Fixed scalar for class 0 */
    uint8_t scalar_fixed[32], point_fixed[32];
    memset(scalar_fixed, 0x42, 32);
    memset(point_fixed, 0x09, 32);  /* Basepoint-ish */

    uint8_t scalar[32], point[32], out[32];

    /* Warmup */
    for (int i = 0; i < WARMUP_ITERATIONS; i++) {
        cfx_rand_bytes(scalar, 32);
        cfx_x25519(out, scalar, point_fixed);
    }

    for (int i = 0; i < num_measurements; i++) {
        int class_id = rand_bit();

        if (class_id) {
            cfx_rand_bytes(scalar, 32);
            cfx_rand_bytes(point, 32);
        } else {
            memcpy(scalar, scalar_fixed, 32);
            memcpy(point, point_fixed, 32);
        }

        uint64_t t0 = get_time();
        cfx_x25519(out, scalar, point);
        uint64_t t1 = get_time();

        sink = out[0];

        ttest_push(class_id ? &class1 : &class0, (double)(t1 - t0));

        if (i > (num_measurements / 10) && (i % (num_measurements / 10)) == 0) {
            double t = ttest_compute(&class0, &class1);
            if (verbose) printf("  %d: t = %.3f\n", i, t);
            if (fail_fast && fabs(t) > THRESHOLD) {
                printf("  FAIL (timing leakage detected)\n");
                return 1;
            }
        }
    }

    double t = ttest_compute(&class0, &class1);
    printf("  final: t = %.3f %s\n", t, fabs(t) > THRESHOLD ? "FAIL" : "PASS");
    return fabs(t) > THRESHOLD ? 1 : 0;
}

/*============================================================================
 * Poly1305 tests
 *============================================================================*/

static int test_poly1305(void) {
    printf("Testing cfx_poly1305()...\n");

    ttest_ctx class0, class1;
    ttest_init(&class0);
    ttest_init(&class1);

    uint8_t key_fixed[32], msg_fixed[256];
    memset(key_fixed, 0x42, 32);
    memset(msg_fixed, 0x00, 256);

    uint8_t key[32], msg[256], tag[16];

    /* Warmup */
    for (int i = 0; i < WARMUP_ITERATIONS; i++) {
        cfx_rand_bytes(key, 32);
        cfx_poly1305(tag, msg_fixed, 256, key);
    }

    for (int i = 0; i < num_measurements; i++) {
        int class_id = rand_bit();

        if (class_id) {
            cfx_rand_bytes(key, 32);
            cfx_rand_bytes(msg, 256);
        } else {
            memcpy(key, key_fixed, 32);
            memcpy(msg, msg_fixed, 256);
        }

        uint64_t t0 = get_time();
        cfx_poly1305(tag, msg, 256, key);
        uint64_t t1 = get_time();

        sink = tag[0];

        ttest_push(class_id ? &class1 : &class0, (double)(t1 - t0));

        if (i > (num_measurements / 10) && (i % (num_measurements / 10)) == 0) {
            double t = ttest_compute(&class0, &class1);
            if (verbose) printf("  %d: t = %.3f\n", i, t);
            if (fail_fast && fabs(t) > THRESHOLD) {
                printf("  FAIL (timing leakage detected)\n");
                return 1;
            }
        }
    }

    double t = ttest_compute(&class0, &class1);
    printf("  final: t = %.3f %s\n", t, fabs(t) > THRESHOLD ? "FAIL" : "PASS");
    return fabs(t) > THRESHOLD ? 1 : 0;
}

/* Constant-time compare (for verify test) */
static int ct_compare(const uint8_t *a, const uint8_t *b, size_t len) {
    uint8_t diff = 0;
    for (size_t i = 0; i < len; i++) {
        diff |= a[i] ^ b[i];
    }
    return diff == 0 ? 0 : -1;
}

static int test_poly1305_verify(void) {
    printf("Testing poly1305 verify (ct_compare on valid vs invalid tag)...\n");

    ttest_ctx class0, class1;
    ttest_init(&class0);
    ttest_init(&class1);

    uint8_t key[32], msg[256], tag_valid[16], tag_invalid[16], tag_computed[16];
    cfx_rand_bytes(key, 32);
    cfx_rand_bytes(msg, 256);
    cfx_poly1305(tag_valid, msg, 256, key);
    memcpy(tag_invalid, tag_valid, 16);
    tag_invalid[0] ^= 0x01;

    uint8_t tag_test[16];

    /* Warmup */
    for (int i = 0; i < WARMUP_ITERATIONS; i++) {
        cfx_poly1305(tag_computed, msg, 256, key);
        ct_compare(tag_computed, tag_valid, 16);
    }

    for (int i = 0; i < num_measurements; i++) {
        int class_id = rand_bit();
        memcpy(tag_test, class_id ? tag_invalid : tag_valid, 16);

        uint64_t t0 = get_time();
        cfx_poly1305(tag_computed, msg, 256, key);
        volatile int r = ct_compare(tag_computed, tag_test, 16);
        uint64_t t1 = get_time();
        (void)r;

        ttest_push(class_id ? &class1 : &class0, (double)(t1 - t0));

        if (i > (num_measurements / 10) && (i % (num_measurements / 10)) == 0) {
            double t = ttest_compute(&class0, &class1);
            if (verbose) printf("  %d: t = %.3f\n", i, t);
            if (fail_fast && fabs(t) > THRESHOLD) {
                printf("  FAIL (timing leakage detected)\n");
                return 1;
            }
        }
    }

    double t = ttest_compute(&class0, &class1);
    printf("  final: t = %.3f %s\n", t, fabs(t) > THRESHOLD ? "FAIL" : "PASS");
    return fabs(t) > THRESHOLD ? 1 : 0;
}

/*============================================================================
 * ChaCha20 tests
 *============================================================================*/

static int test_chacha20(void) {
    printf("Testing cfx_chacha20_block()...\n");

    ttest_ctx class0, class1;
    ttest_init(&class0);
    ttest_init(&class1);

    uint8_t key_fixed[32], nonce_fixed[12];
    memset(key_fixed, 0x42, 32);
    memset(nonce_fixed, 0x00, 12);

    uint8_t key[32], nonce[12], out[64];
    cfx_chacha20_ctx_t ctx;

    /* Warmup */
    for (int i = 0; i < WARMUP_ITERATIONS; i++) {
        cfx_rand_bytes(key, 32);
        cfx_chacha20_ctx_init(&ctx, key, nonce_fixed);
        cfx_chacha20_block(&ctx, 0, out);
    }

    for (int i = 0; i < num_measurements; i++) {
        int class_id = rand_bit();

        if (class_id) {
            cfx_rand_bytes(key, 32);
            cfx_rand_bytes(nonce, 12);
        } else {
            memcpy(key, key_fixed, 32);
            memcpy(nonce, nonce_fixed, 12);
        }

        cfx_chacha20_ctx_init(&ctx, key, nonce);

        uint64_t t0 = get_time();
        cfx_chacha20_block(&ctx, 0, out);
        uint64_t t1 = get_time();

        sink = out[0];

        ttest_push(class_id ? &class1 : &class0, (double)(t1 - t0));

        if (i > (num_measurements / 10) && (i % (num_measurements / 10)) == 0) {
            double t = ttest_compute(&class0, &class1);
            if (verbose) printf("  %d: t = %.3f\n", i, t);
            if (fail_fast && fabs(t) > THRESHOLD) {
                printf("  FAIL (timing leakage detected)\n");
                return 1;
            }
        }
    }

    double t = ttest_compute(&class0, &class1);
    printf("  final: t = %.3f %s\n", t, fabs(t) > THRESHOLD ? "FAIL" : "PASS");
    return fabs(t) > THRESHOLD ? 1 : 0;
}

/*============================================================================
 * SHA-512 tests
 *============================================================================*/

static int test_sha512(void) {
    printf("Testing cfx_sha512()...\n");

    ttest_ctx class0, class1;
    ttest_init(&class0);
    ttest_init(&class1);

    uint8_t msg_fixed[256];
    memset(msg_fixed, 0x00, 256);

    uint8_t msg[256], hash[64];

    /* Warmup */
    for (int i = 0; i < WARMUP_ITERATIONS; i++) {
        cfx_rand_bytes(msg, 256);
        cfx_sha512(hash, msg, 256);
    }

    for (int i = 0; i < num_measurements; i++) {
        int class_id = rand_bit();

        if (class_id) {
            cfx_rand_bytes(msg, 256);
        } else {
            memcpy(msg, msg_fixed, 256);
        }

        uint64_t t0 = get_time();
        cfx_sha512(hash, msg, 256);
        uint64_t t1 = get_time();

        sink = hash[0];

        ttest_push(class_id ? &class1 : &class0, (double)(t1 - t0));

        if (i > (num_measurements / 10) && (i % (num_measurements / 10)) == 0) {
            double t = ttest_compute(&class0, &class1);
            if (verbose) printf("  %d: t = %.3f\n", i, t);
            if (fail_fast && fabs(t) > THRESHOLD) {
                printf("  FAIL (timing leakage detected)\n");
                return 1;
            }
        }
    }

    double t = ttest_compute(&class0, &class1);
    printf("  final: t = %.3f %s\n", t, fabs(t) > THRESHOLD ? "FAIL" : "PASS");
    return fabs(t) > THRESHOLD ? 1 : 0;
}

/*============================================================================
 * Edge case tests
 *============================================================================*/

/* Ed25519 sign with edge-case seeds: leading zeros vs random
 * Seeds with leading zero bytes produce scalars with leading zeros,
 * which could cause timing variations in scalar multiplication. */
static int test_ed25519_sign_edge(void) {
    printf("Testing cfx_ed25519_sign() [edge: leading-zero seeds]...\n");

    ttest_ctx class0, class1;
    ttest_init(&class0);
    ttest_init(&class1);

    /* Class 0: seed with many leading zeros */
    uint8_t seed_zeros[32], msg[64];
    memset(seed_zeros, 0x00, 32);
    seed_zeros[31] = 0x01;  /* Minimal non-zero seed */

    uint8_t seed[32], pk[32], sk[64], sig[64];

    /* Warmup */
    for (int i = 0; i < WARMUP_ITERATIONS; i++) {
        cfx_rand_bytes(seed, 32);
        cfx_ed25519_create_keypair(pk, sk, seed);
        cfx_ed25519_sign(sig, msg, 64, sk);
    }

    for (int i = 0; i < num_measurements; i++) {
        int class_id = rand_bit();
        cfx_rand_bytes(msg, 64);

        if (class_id) {
            cfx_rand_bytes(seed, 32);
        } else {
            memcpy(seed, seed_zeros, 32);
        }

        cfx_ed25519_create_keypair(pk, sk, seed);

        uint64_t t0 = get_time();
        cfx_ed25519_sign(sig, msg, 64, sk);
        uint64_t t1 = get_time();

        sink = sig[0];

        ttest_push(class_id ? &class1 : &class0, (double)(t1 - t0));

        if (i > (num_measurements / 10) && (i % (num_measurements / 10)) == 0) {
            double t = ttest_compute(&class0, &class1);
            if (verbose) printf("  %d: t = %.3f\n", i, t);
            if (fail_fast && fabs(t) > THRESHOLD) {
                printf("  FAIL (timing leakage detected)\n");
                return 1;
            }
        }
    }

    double t = ttest_compute(&class0, &class1);
    printf("  final: t = %.3f %s\n", t, fabs(t) > THRESHOLD ? "FAIL" : "PASS");
    return fabs(t) > THRESHOLD ? 1 : 0;
}

/* X25519 with low-order points
 * These are the 8 points of small order on Curve25519.
 * A non-constant-time implementation might handle these differently. */
static int test_x25519_loworder(void) {
    printf("Testing cfx_x25519() [edge: low-order points]...\n");

    ttest_ctx class0, class1;
    ttest_init(&class0);
    ttest_init(&class1);

    /* Low-order point of order 8 (one of the torsion points) */
    static const uint8_t low_order_point[32] = {
        0xe0, 0xeb, 0x7a, 0x7c, 0x3b, 0x41, 0xb8, 0xae,
        0x16, 0x56, 0xe3, 0xfa, 0xf1, 0x9f, 0xc4, 0x6a,
        0xda, 0x09, 0x8d, 0xeb, 0x9c, 0x32, 0xb1, 0xfd,
        0x86, 0x62, 0x05, 0x16, 0x5f, 0x49, 0xb8, 0x00
    };

    /* Normal basepoint */
    static const uint8_t basepoint[32] = {
        0x09, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00
    };

    uint8_t scalar[32], point[32], out[32];

    /* Warmup */
    for (int i = 0; i < WARMUP_ITERATIONS; i++) {
        cfx_rand_bytes(scalar, 32);
        cfx_x25519(out, scalar, basepoint);
    }

    for (int i = 0; i < num_measurements; i++) {
        int class_id = rand_bit();
        cfx_rand_bytes(scalar, 32);

        if (class_id) {
            memcpy(point, basepoint, 32);
        } else {
            memcpy(point, low_order_point, 32);
        }

        uint64_t t0 = get_time();
        cfx_x25519(out, scalar, point);
        uint64_t t1 = get_time();

        sink = out[0];

        ttest_push(class_id ? &class1 : &class0, (double)(t1 - t0));

        if (i > (num_measurements / 10) && (i % (num_measurements / 10)) == 0) {
            double t = ttest_compute(&class0, &class1);
            if (verbose) printf("  %d: t = %.3f\n", i, t);
            if (fail_fast && fabs(t) > THRESHOLD) {
                printf("  FAIL (timing leakage detected)\n");
                return 1;
            }
        }
    }

    double t = ttest_compute(&class0, &class1);
    printf("  final: t = %.3f %s\n", t, fabs(t) > THRESHOLD ? "FAIL" : "PASS");
    return fabs(t) > THRESHOLD ? 1 : 0;
}

/* X25519 with edge-case scalars: many leading zeros vs random */
static int test_x25519_scalar_edge(void) {
    printf("Testing cfx_x25519() [edge: leading-zero scalars]...\n");

    ttest_ctx class0, class1;
    ttest_init(&class0);
    ttest_init(&class1);

    /* Scalar with many leading zeros */
    uint8_t scalar_zeros[32];
    memset(scalar_zeros, 0x00, 32);
    scalar_zeros[0] = 0x08;  /* Minimal after clamping */

    static const uint8_t basepoint[32] = { 9 };

    uint8_t scalar[32], out[32];

    /* Warmup */
    for (int i = 0; i < WARMUP_ITERATIONS; i++) {
        cfx_rand_bytes(scalar, 32);
        cfx_x25519(out, scalar, basepoint);
    }

    for (int i = 0; i < num_measurements; i++) {
        int class_id = rand_bit();

        if (class_id) {
            cfx_rand_bytes(scalar, 32);
        } else {
            memcpy(scalar, scalar_zeros, 32);
        }

        uint64_t t0 = get_time();
        cfx_x25519(out, scalar, basepoint);
        uint64_t t1 = get_time();

        sink = out[0];

        ttest_push(class_id ? &class1 : &class0, (double)(t1 - t0));

        if (i > (num_measurements / 10) && (i % (num_measurements / 10)) == 0) {
            double t = ttest_compute(&class0, &class1);
            if (verbose) printf("  %d: t = %.3f\n", i, t);
            if (fail_fast && fabs(t) > THRESHOLD) {
                printf("  FAIL (timing leakage detected)\n");
                return 1;
            }
        }
    }

    double t = ttest_compute(&class0, &class1);
    printf("  final: t = %.3f %s\n", t, fabs(t) > THRESHOLD ? "FAIL" : "PASS");
    return fabs(t) > THRESHOLD ? 1 : 0;
}

/* Poly1305 with edge-case keys: r=0 vs random
 * The r portion of the key affects multiplication; r=0 is a degenerate case. */
static int test_poly1305_edge(void) {
    printf("Testing cfx_poly1305() [edge: r=0 key]...\n");

    ttest_ctx class0, class1;
    ttest_init(&class0);
    ttest_init(&class1);

    /* Key with r=0 (first 16 bytes = 0, s = random) */
    uint8_t key_rzero[32];
    memset(key_rzero, 0x00, 16);      /* r = 0 */
    memset(key_rzero + 16, 0x42, 16); /* s = fixed */

    uint8_t key[32], msg[256], tag[16];

    /* Warmup */
    for (int i = 0; i < WARMUP_ITERATIONS; i++) {
        cfx_rand_bytes(key, 32);
        cfx_rand_bytes(msg, 256);
        cfx_poly1305(tag, msg, 256, key);
    }

    for (int i = 0; i < num_measurements; i++) {
        int class_id = rand_bit();
        cfx_rand_bytes(msg, 256);

        if (class_id) {
            cfx_rand_bytes(key, 32);
        } else {
            memcpy(key, key_rzero, 32);
        }

        uint64_t t0 = get_time();
        cfx_poly1305(tag, msg, 256, key);
        uint64_t t1 = get_time();

        sink = tag[0];

        ttest_push(class_id ? &class1 : &class0, (double)(t1 - t0));

        if (i > (num_measurements / 10) && (i % (num_measurements / 10)) == 0) {
            double t = ttest_compute(&class0, &class1);
            if (verbose) printf("  %d: t = %.3f\n", i, t);
            if (fail_fast && fabs(t) > THRESHOLD) {
                printf("  FAIL (timing leakage detected)\n");
                return 1;
            }
        }
    }

    double t = ttest_compute(&class0, &class1);
    printf("  final: t = %.3f %s\n", t, fabs(t) > THRESHOLD ? "FAIL" : "PASS");
    return fabs(t) > THRESHOLD ? 1 : 0;
}

/*============================================================================
 * Main
 *============================================================================*/

static void usage(const char *prog) {
    printf("Usage: %s [options] [test...]\n", prog);
    printf("  Constant-time testing using Welch's t-test.\n\n");
    printf("Options:\n");
    printf("  -n <count>       Number of measurements (default: %d)\n", DEFAULT_MEASUREMENTS);
    printf("  -q, --quiet      Less verbose output\n");
    printf("  -c, --continue   Don't stop on first failure\n");
    printf("  -h, --help       Show this help\n\n");
    printf("Tests:\n");
    printf("  ed25519          Ed25519 sign/verify (incl. edge cases)\n");
    printf("  x25519           X25519 (incl. low-order, edge cases)\n");
    printf("  poly1305         Poly1305 MAC/verify (incl. edge cases)\n");
    printf("  chacha20         ChaCha20 block\n");
    printf("  sha512           SHA-512 hash\n");
    printf("  all              Run all tests (default)\n\n");
    printf("A test passes if |t| < %.1f after %d measurements.\n",
        THRESHOLD, num_measurements);
}

int main(int argc, char **argv) {
    int run_all = 1;
    int run_ed25519 = 0;
    int run_x25519 = 0;
    int run_poly1305 = 0;
    int run_chacha20 = 0;
    int run_sha512 = 0;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            usage(argv[0]);
            return 0;
        } else if (strcmp(argv[i], "-q") == 0 || strcmp(argv[i], "--quiet") == 0) {
            verbose = 0;
        } else if (strcmp(argv[i], "-c") == 0 || strcmp(argv[i], "--continue") == 0) {
            fail_fast = 0;
        } else if (strcmp(argv[i], "-n") == 0) {
            if (++i >= argc) {
                fprintf(stderr, "error: -n requires argument\n");
                return 1;
            }
            num_measurements = atoi(argv[i]);
            if (num_measurements < 1000) {
                fprintf(stderr, "error: need at least 1000 measurements\n");
                return 1;
            }
        } else if (strcmp(argv[i], "all") == 0) {
            run_all = 1;
        } else if (strcmp(argv[i], "ed25519") == 0) {
            run_ed25519 = 1; run_all = 0;
        } else if (strcmp(argv[i], "x25519") == 0) {
            run_x25519 = 1; run_all = 0;
        } else if (strcmp(argv[i], "poly1305") == 0) {
            run_poly1305 = 1; run_all = 0;
        } else if (strcmp(argv[i], "chacha20") == 0) {
            run_chacha20 = 1; run_all = 0;
        } else if (strcmp(argv[i], "sha512") == 0) {
            run_sha512 = 1; run_all = 0;
        } else {
            fprintf(stderr, "Unknown option or test: %s\n", argv[i]);
            usage(argv[0]);
            return 1;
        }
    }

    printf("Constant-time test (Welch's t-test)\n");
    printf("Measurements: %d, threshold: |t| > %.1f\n\n",
        num_measurements, THRESHOLD);

    /* Seed RNG from OS entropy */
    cfx_srand_os();

    int failures = 0;

    if (run_all || run_ed25519) {
        failures += test_ed25519_sign();
        failures += test_ed25519_verify();
        failures += test_ed25519_sign_edge();
        printf("\n");
    }
    if (run_all || run_x25519) {
        failures += test_x25519();
        failures += test_x25519_loworder();
        failures += test_x25519_scalar_edge();
        printf("\n");
    }
    if (run_all || run_poly1305) {
        failures += test_poly1305();
        failures += test_poly1305_verify();
        failures += test_poly1305_edge();
        printf("\n");
    }
    if (run_all || run_chacha20) {
        failures += test_chacha20();
        printf("\n");
    }
    if (run_all || run_sha512) {
        failures += test_sha512();
        printf("\n");
    }

    if (failures > 0) {
        printf("=== %d test(s) FAILED ===\n", failures);
        return 1;
    } else {
        printf("=== All tests PASSED ===\n");
        return 0;
    }
}
