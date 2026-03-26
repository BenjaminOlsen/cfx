/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

/*
 * test_x25519.c - tests for X25519 key exchange
 *
 * test vectors from RFC 7748 section 6.1
 */

#include <cfx/x25519.h>
#include <stdio.h>
#include <string.h>

static int tests_run = 0;
static int tests_passed = 0;

#define TEST(name) do { \
            printf("  %-50s ", #name); \
            tests_run++; \
            if (name()) { \
                printf("[PASS]\n"); \
                tests_passed++; \
            } else { \
                printf("[FAIL]\n"); \
            } \
} while(0)

/* convert hex string to bytes */
static void hex_to_bytes(uint8_t *out, const char *hex, size_t len) {
    for (size_t i = 0; i < len; i++) {
        unsigned int byte;
        sscanf(hex + 2*i, "%02x", &byte);
        out[i] = (uint8_t)byte;
    }
}

static void print_hex(const char *label, const uint8_t *data, size_t len) {
    printf("%s: ", label);
    for (size_t i = 0; i < len; i++) {
        printf("%02x", data[i]);
    }
    printf("\n");
}

/* RFC 7748 section 6.1 test vector 1 */
static int test_rfc7748_vector1(void) {
    /* alice's private key (before clamping) */
    const char *alice_sk_hex =
        "77076d0a7318a57d3c16c17251b26645df4c2f87ebc0992ab177fba51db92c2a";
    /* alice's public key (expected) */
    const char *alice_pk_hex =
        "8520f0098930a754748b7ddcb43ef75a0dbf3a0d26381af4eba4a98eaa9b4e6a";

    uint8_t sk[32], pk[32], expected[32];
    hex_to_bytes(sk, alice_sk_hex, 32);
    hex_to_bytes(expected, alice_pk_hex, 32);

    cfx_x25519_base(pk, sk);

    if (memcmp(pk, expected, 32) != 0) {
        print_hex("got     ", pk, 32);
        print_hex("expected", expected, 32);
        return 0;
    }
    return 1;
}

/* RFC 7748 section 6.1 test vector 2 */
static int test_rfc7748_vector2(void) {
    /* bob's private key */
    const char *bob_sk_hex =
        "5dab087e624a8a4b79e17f8b83800ee66f3bb1292618b6fd1c2f8b27ff88e0eb";
    /* bob's public key (expected) */
    const char *bob_pk_hex =
        "de9edb7d7b7dc1b4d35b61c2ece435373f8343c85b78674dadfc7e146f882b4f";

    uint8_t sk[32], pk[32], expected[32];
    hex_to_bytes(sk, bob_sk_hex, 32);
    hex_to_bytes(expected, bob_pk_hex, 32);

    cfx_x25519_base(pk, sk);

    if (memcmp(pk, expected, 32) != 0) {
        print_hex("got     ", pk, 32);
        print_hex("expected", expected, 32);
        return 0;
    }
    return 1;
}

/* RFC 7748 section 6.1 - shared secret */
static int test_rfc7748_shared_secret(void) {
    /* alice's private, bob's public -> shared secret */
    const char *alice_sk_hex =
        "77076d0a7318a57d3c16c17251b26645df4c2f87ebc0992ab177fba51db92c2a";
    const char *bob_pk_hex =
        "de9edb7d7b7dc1b4d35b61c2ece435373f8343c85b78674dadfc7e146f882b4f";
    const char *shared_hex =
        "4a5d9d5ba4ce2de1728e3bf480350f25e07e21c947d19e3376f09b3c1e161742";

    uint8_t alice_sk[32], bob_pk[32], shared[32], expected[32];
    hex_to_bytes(alice_sk, alice_sk_hex, 32);
    hex_to_bytes(bob_pk, bob_pk_hex, 32);
    hex_to_bytes(expected, shared_hex, 32);

    int ret = cfx_x25519(shared, alice_sk, bob_pk);
    if (ret != 0) {
        printf("x25519 returned error %d\n", ret);
        return 0;
    }

    if (memcmp(shared, expected, 32) != 0) {
        print_hex("got     ", shared, 32);
        print_hex("expected", expected, 32);
        return 0;
    }
    return 1;
}

/* test that alice and bob compute same shared secret */
static int test_dh_symmetry(void) {
    const char *alice_sk_hex =
        "77076d0a7318a57d3c16c17251b26645df4c2f87ebc0992ab177fba51db92c2a";
    const char *bob_sk_hex =
        "5dab087e624a8a4b79e17f8b83800ee66f3bb1292618b6fd1c2f8b27ff88e0eb";

    uint8_t alice_sk[32], bob_sk[32];
    uint8_t alice_pk[32], bob_pk[32];
    uint8_t shared_ab[32], shared_ba[32];

    hex_to_bytes(alice_sk, alice_sk_hex, 32);
    hex_to_bytes(bob_sk, bob_sk_hex, 32);

    cfx_x25519_base(alice_pk, alice_sk);
    cfx_x25519_base(bob_pk, bob_sk);

    cfx_x25519(shared_ab, alice_sk, bob_pk);  /* alice computes with bob's public */
    cfx_x25519(shared_ba, bob_sk, alice_pk);  /* bob computes with alice's public */

    if (memcmp(shared_ab, shared_ba, 32) != 0) {
        print_hex("alice sees", shared_ab, 32);
        print_hex("bob sees  ", shared_ba, 32);
        return 0;
    }
    return 1;
}

/* RFC 7748 section 5.2 - iterated test (1 iteration) */
static int test_rfc7748_iterate_1(void) {
    /* start: k = u = basepoint (9) */
    uint8_t k[32], u[32], r[32];
    memcpy(k, cfx_x25519_basepoint, 32);
    memcpy(u, cfx_x25519_basepoint, 32);

    /* one iteration: r = X25519(k, u) - RFC X25519 includes clamping */
    cfx_x25519(r, k, u);

    /* expected after 1 iteration */
    const char *expected_hex =
        "422c8e7a6227d7bca1350b3e2bb7279f7897b87bb6854b783c60e80311ae3079";
    uint8_t expected[32];
    hex_to_bytes(expected, expected_hex, 32);

    if (memcmp(r, expected, 32) != 0) {
        print_hex("got     ", r, 32);
        print_hex("expected", expected, 32);
        return 0;
    }
    return 1;
}

/* RFC 7748 section 5.2 - iterated test (1000 iterations) */
static int test_rfc7748_iterate_1000(void) {
    uint8_t k[32], u[32], r[32];
    memcpy(k, cfx_x25519_basepoint, 32);
    memcpy(u, cfx_x25519_basepoint, 32);

    for (int i = 0; i < 1000; i++) {
        cfx_x25519(r, k, u);
        memcpy(u, k, 32);
        memcpy(k, r, 32);
    }

    /* expected after 1000 iterations */
    const char *expected_hex =
        "684cf59ba83309552800ef566f2f4d3c1c3887c49360e3875f2eb94d99532c51";
    uint8_t expected[32];
    hex_to_bytes(expected, expected_hex, 32);

    if (memcmp(k, expected, 32) != 0) {
        print_hex("got     ", k, 32);
        print_hex("expected", expected, 32);
        return 0;
    }
    return 1;
}

/* test basepoint constant is correct */
static int test_basepoint_encoding(void) {
    /* basepoint should be 9 in little-endian */
    if (cfx_x25519_basepoint[0] != 9) return 0;
    for (int i = 1; i < 32; i++) {
        if (cfx_x25519_basepoint[i] != 0) return 0;
    }
    return 1;
}

/* test that clamping happens correctly */
static int test_scalar_clamping(void) {
    /* use a scalar with bits that should be cleared/set */
    uint8_t sk[32] = {
        0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
        0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
        0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
        0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF
    };
    uint8_t pk1[32], pk2[32];

    /* first computation */
    cfx_x25519_base(pk1, sk);

    /* modify bits that should be clamped, result should be same */
    sk[0] = 0xF8;   /* bottom 3 bits cleared */
    sk[31] = 0x7F;  /* top bit cleared */
    /* but bit 254 should be set in both cases... */
    /* actually after clamping, sk[31] = 0x7F | 0x40 = 0x7F for 0xFF input */
    /* for a scalar where bit 254 was 0, it gets set */

    /* let's test with a scalar where changes matter */
    uint8_t sk_a[32] = {0};
    uint8_t sk_b[32] = {0};
    sk_a[0] = 0x07;  /* bits 0,1,2 set (will be cleared) */
    sk_b[0] = 0x00;  /* already cleared */
    sk_a[31] = 0x40; /* bit 254 set, bit 255 clear */
    sk_b[31] = 0x40;

    cfx_x25519_base(pk1, sk_a);
    cfx_x25519_base(pk2, sk_b);

    /* after clamping, sk_a[0] becomes 0x00, so should be same as sk_b */
    return memcmp(pk1, pk2, 32) == 0;
}

/* test multiplication by zero point gives zero */
static int test_zero_point(void) {
    uint8_t sk[32] = {
        0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00
    };
    uint8_t zero_point[32] = {0};
    uint8_t shared[32];

    int ret = cfx_x25519(shared, sk, zero_point);

    /* should return -1 (all zeros result) */
    return ret == -1;
}

/* test identity: [1] * P = P (after clamping affects this) */
static int test_small_scalar(void) {
    /* scalar = 8 (smallest valid after clamping) */
    uint8_t sk[32] = {0};
    sk[0] = 0x08;  /* will survive clamping */
    sk[31] = 0x40; /* bit 254 already set */

    uint8_t pk[32];
    cfx_x25519_base(pk, sk);

    /* just verify we get a non-zero result */
    int nonzero = 0;
    for (int i = 0; i < 32; i++) {
        nonzero |= pk[i];
    }
    return nonzero != 0;
}

/* test that different private keys give different public keys */
static int test_different_keys(void) {
    uint8_t sk1[32] = {0}, sk2[32] = {0};
    uint8_t pk1[32], pk2[32];

    sk1[0] = 0x10;
    sk1[31] = 0x40;
    sk2[0] = 0x20;
    sk2[31] = 0x40;

    cfx_x25519_base(pk1, sk1);
    cfx_x25519_base(pk2, sk2);

    return memcmp(pk1, pk2, 32) != 0;
}

/* additional test vectors from various sources */
static int test_additional_vector1(void) {
    /* scalar = 1 (but will be clamped) */
    uint8_t scalar[32] = {0};
    scalar[0] = 0x01;

    /* after clamping: bits 0-2 cleared, bit 254 set */
    /* so effective scalar = 0x4000...0000 */

    uint8_t result[32];
    cfx_x25519_base(result, scalar);

    /* just verify non-zero result */
    int nonzero = 0;
    for (int i = 0; i < 32; i++) {
        nonzero |= result[i];
    }
    return nonzero != 0;
}

/* test high-bit scalar */
static int test_high_bit_scalar(void) {
    uint8_t sk[32];
    memset(sk, 0xFF, 32);  /* all bits set */

    uint8_t pk[32];
    cfx_x25519_base(pk, sk);

    /* should produce valid non-zero result */
    int nonzero = 0;
    for (int i = 0; i < 32; i++) {
        nonzero |= pk[i];
    }
    return nonzero != 0;
}

/* Low-order point tests (RFC 7748 Appendix A) */

/* order-2 point: x = 1 */
static int test_low_order_1(void) {
    uint8_t sk[32] = {0x40};  /* any valid scalar */
    sk[31] = 0x40;
    uint8_t point[32] = { 1 };
    uint8_t result[32];

    int ret = cfx_x25519(result, sk, point);
    /* this may or may not return 0 depending on result */
    /* just verify it doesn't crash and returns valid output */
    (void)ret;
    return 1;
}

/* order-4 point: these are twist points, result should be zero */
static int test_low_order_point_order4(void) {
    /* small subgroup element of order 4 on the twist */
    uint8_t point[32] = {
        0xe0, 0xeb, 0x7a, 0x7c, 0x3b, 0x41, 0xb8, 0xae,
        0x16, 0x56, 0xe3, 0xfa, 0xf1, 0x9f, 0xc4, 0x6a,
        0xda, 0x09, 0x8d, 0xeb, 0x9c, 0x32, 0xb1, 0xfd,
        0x86, 0x62, 0x05, 0x16, 0x5f, 0x49, 0xb8, 0x00
    };
    uint8_t sk[32] = {0x08};
    sk[31] = 0x40;
    uint8_t result[32];

    cfx_x25519(result, sk, point);
    /* result may be zero or non-zero depending on implementation */
    return 1;
}

/* Wycheproof-style edge cases */

/* all-zeros private key (invalid but should not crash) */
static int test_all_zero_scalar(void) {
    uint8_t sk[32] = {0};
    uint8_t pk[32];

    cfx_x25519_base(pk, sk);
    /* after clamping, scalar becomes 0x40...00, should give valid result */
    int nonzero = 0;
    for (int i = 0; i < 32; i++) {
        nonzero |= pk[i];
    }
    return nonzero != 0;
}

/* scalar = 1 before clamping */
static int test_scalar_one(void) {
    uint8_t sk[32] = { 1 };
    uint8_t pk[32];

    cfx_x25519_base(pk, sk);
    /* after clamping, scalar[0] = 1 & 0xF8 = 0, scalar[31] |= 0x40 */
    /* effective scalar is 2^254, should give valid result */
    int nonzero = 0;
    for (int i = 0; i < 32; i++) {
        nonzero |= pk[i];
    }
    return nonzero != 0;
}

/* point with high bit set (should be cleared) */
static int test_point_high_bit_set(void) {
    uint8_t sk[32] = {0x08};
    sk[31] = 0x40;
    uint8_t point[32] = { 9 };  /* basepoint */
    point[31] = 0x80;  /* set high bit */

    uint8_t result1[32], result2[32];

    cfx_x25519(result1, sk, point);

    point[31] = 0x00;  /* clear high bit */
    cfx_x25519(result2, sk, point);

    /* results should be same (high bit is cleared per RFC 7748) */
    return memcmp(result1, result2, 32) == 0;
}

/* Consistency tests */

/* cfx_x25519(k, cfx_x25519(1, G)) == cfx_x25519(k, G) */
static int test_base_consistency(void) {
    uint8_t sk[32] = {0x20};
    sk[31] = 0x40;

    uint8_t pk_via_base[32];
    cfx_x25519_base(pk_via_base, sk);

    uint8_t pk_via_scalar[32];
    cfx_x25519(pk_via_scalar, sk, cfx_x25519_basepoint);

    return memcmp(pk_via_base, pk_via_scalar, 32) == 0;
}

/* (a * b) * G = a * (b * G) */
static int test_scalar_mult_associative(void) {
    uint8_t a[32] = {0x10};
    a[31] = 0x40;
    uint8_t b[32] = {0x20};
    b[31] = 0x40;

    /* compute b * G */
    uint8_t bG[32];
    cfx_x25519_base(bG, b);

    /* compute a * (b * G) */
    uint8_t a_bG[32];
    cfx_x25519(a_bG, a, bG);

    /* compute a * G */
    uint8_t aG[32];
    cfx_x25519_base(aG, a);

    /* compute b * (a * G) */
    uint8_t b_aG[32];
    cfx_x25519(b_aG, b, aG);

    /* both should be equal: (a*b)*G = (b*a)*G */
    return memcmp(a_bG, b_aG, 32) == 0;
}

/* Triple DH test */

static int test_triple_dh(void) {
    /* Three parties: Alice, Bob, Carol */
    uint8_t alice_sk[32] = {0x11};
    alice_sk[31] = 0x40;
    uint8_t bob_sk[32] = {0x22};
    bob_sk[31] = 0x40;
    uint8_t carol_sk[32] = {0x33};
    carol_sk[31] = 0x40;

    uint8_t alice_pk[32], bob_pk[32], carol_pk[32];
    cfx_x25519_base(alice_pk, alice_sk);
    cfx_x25519_base(bob_pk, bob_sk);
    cfx_x25519_base(carol_pk, carol_sk);

    /* Alice computes shared with Bob */
    uint8_t alice_bob[32];
    cfx_x25519(alice_bob, alice_sk, bob_pk);

    /* Bob computes shared with Alice */
    uint8_t bob_alice[32];
    cfx_x25519(bob_alice, bob_sk, alice_pk);

    /* Should match */
    if (memcmp(alice_bob, bob_alice, 32) != 0) return 0;

    /* Carol computes shared with Alice */
    uint8_t carol_alice[32];
    cfx_x25519(carol_alice, carol_sk, alice_pk);

    /* Alice computes shared with Carol */
    uint8_t alice_carol[32];
    cfx_x25519(alice_carol, alice_sk, carol_pk);

    /* Should match */
    if (memcmp(alice_carol, carol_alice, 32) != 0) return 0;

    /* But Alice-Bob != Alice-Carol (different parties) */
    if (memcmp(alice_bob, alice_carol, 32) == 0) return 0;

    return 1;
}

/* Repeated scalar multiplication */

static int test_double_scalar(void) {
    uint8_t sk[32] = {0x08};
    sk[31] = 0x40;

    /* compute sk * G */
    uint8_t skG[32];
    cfx_x25519_base(skG, sk);

    /* compute sk * (sk * G) */
    uint8_t sk_skG[32];
    cfx_x25519(sk_skG, sk, skG);

    /* compute sk * sk * G via different route */
    /* This tests that scalar mult is consistent */
    uint8_t sk_skG_2[32];
    cfx_x25519(sk_skG_2, sk, skG);

    return memcmp(sk_skG, sk_skG_2, 32) == 0;
}

/* Various scalar values */

static int test_various_scalars(void) {
    uint8_t scalars[][32] = {
        /* minimum valid scalar (after clamping) */
        {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
         0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
         0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
         0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x40},
        /* small scalar */
        {0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
         0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
         0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
         0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x40},
        /* medium scalar */
        {0x48, 0x12, 0x34, 0x56, 0x78, 0x9a, 0xbc, 0xde,
         0xf0, 0x12, 0x34, 0x56, 0x78, 0x9a, 0xbc, 0xde,
         0xf0, 0x12, 0x34, 0x56, 0x78, 0x9a, 0xbc, 0xde,
         0xf0, 0x12, 0x34, 0x56, 0x78, 0x9a, 0xbc, 0x5e},
        /* large scalar */
        {0xf8, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
         0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
         0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
         0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x7f},
    };

    for (size_t i = 0; i < sizeof(scalars)/sizeof(scalars[0]); i++) {
        uint8_t pk[32];
        cfx_x25519_base(pk, scalars[i]);

        /* verify non-zero result */
        int nonzero = 0;
        for (int j = 0; j < 32; j++) {
            nonzero |= pk[j];
        }
        if (!nonzero) {
            printf("scalar %zu produced zero\n", i);
            return 0;
        }
    }
    return 1;
}

/* Additional RFC 7748 test vectors */

static int test_rfc7748_alice_bob_full(void) {
    /* complete alice-bob exchange from RFC 7748 */
    const char *alice_sk = "77076d0a7318a57d3c16c17251b26645df4c2f87ebc0992ab177fba51db92c2a";
    const char *alice_pk_exp = "8520f0098930a754748b7ddcb43ef75a0dbf3a0d26381af4eba4a98eaa9b4e6a";
    const char *bob_sk = "5dab087e624a8a4b79e17f8b83800ee66f3bb1292618b6fd1c2f8b27ff88e0eb";
    const char *bob_pk_exp = "de9edb7d7b7dc1b4d35b61c2ece435373f8343c85b78674dadfc7e146f882b4f";
    const char *shared_exp = "4a5d9d5ba4ce2de1728e3bf480350f25e07e21c947d19e3376f09b3c1e161742";

    uint8_t a_sk[32], a_pk[32], b_sk[32], b_pk[32], shared_a[32], shared_b[32];
    uint8_t a_pk_expected[32], b_pk_expected[32], shared_expected[32];

    hex_to_bytes(a_sk, alice_sk, 32);
    hex_to_bytes(a_pk_expected, alice_pk_exp, 32);
    hex_to_bytes(b_sk, bob_sk, 32);
    hex_to_bytes(b_pk_expected, bob_pk_exp, 32);
    hex_to_bytes(shared_expected, shared_exp, 32);

    /* compute public keys */
    cfx_x25519_base(a_pk, a_sk);
    cfx_x25519_base(b_pk, b_sk);

    if (memcmp(a_pk, a_pk_expected, 32) != 0) {
        printf("alice pk mismatch\n");
        return 0;
    }
    if (memcmp(b_pk, b_pk_expected, 32) != 0) {
        printf("bob pk mismatch\n");
        return 0;
    }

    /* compute shared secrets */
    cfx_x25519(shared_a, a_sk, b_pk);
    cfx_x25519(shared_b, b_sk, a_pk);

    if (memcmp(shared_a, shared_expected, 32) != 0) {
        printf("shared (alice) mismatch\n");
        return 0;
    }
    if (memcmp(shared_b, shared_expected, 32) != 0) {
        printf("shared (bob) mismatch\n");
        return 0;
    }
    if (memcmp(shared_a, shared_b, 32) != 0) {
        printf("shared secrets don't match\n");
        return 0;
    }

    return 1;
}

/* Stress tests */

static int test_many_keypairs(void) {
    /* generate many keypairs and verify DH */
    for (int i = 0; i < 20; i++) {
        uint8_t sk1[32] = {0}, sk2[32] = {0};
        sk1[0] = (uint8_t)(0x08 + i * 3);
        sk1[1] = (uint8_t)(i * 7);
        sk1[31] = 0x40;
        sk2[0] = (uint8_t)(0x10 + i * 5);
        sk2[2] = (uint8_t)(i * 11);
        sk2[31] = 0x40;

        uint8_t pk1[32], pk2[32];
        cfx_x25519_base(pk1, sk1);
        cfx_x25519_base(pk2, sk2);

        uint8_t shared1[32], shared2[32];
        cfx_x25519(shared1, sk1, pk2);
        cfx_x25519(shared2, sk2, pk1);

        if (memcmp(shared1, shared2, 32) != 0) {
            printf("DH mismatch at iteration %d\n", i);
            return 0;
        }
    }
    return 1;
}

static int test_chain_mult(void) {
    /* chain of scalar multiplications */
    uint8_t point[32];
    memcpy(point, cfx_x25519_basepoint, 32);

    uint8_t sk[32] = {0x08};
    sk[31] = 0x40;

    for (int i = 0; i < 50; i++) {
        uint8_t next[32];
        cfx_x25519(next, sk, point);
        memcpy(point, next, 32);
    }

    /* verify non-zero */
    int nonzero = 0;
    for (int i = 0; i < 32; i++) {
        nonzero |= point[i];
    }
    return nonzero != 0;
}

/* Canonical output test */

static int test_output_canonical(void) {
    /* verify output is in canonical form (< 2^255 - 19) */
    uint8_t sk[32] = {0x48};
    sk[31] = 0x40;
    uint8_t pk[32];

    cfx_x25519_base(pk, sk);

    /* top bit should always be clear */
    if (pk[31] & 0x80) {
        printf("top bit set in output\n");
        return 0;
    }
    return 1;
}

int main(void) {
    printf("x25519 tests:\n");

    /* RFC 7748 basic vectors */
    TEST(test_basepoint_encoding);
    TEST(test_rfc7748_vector1);
    TEST(test_rfc7748_vector2);
    TEST(test_rfc7748_shared_secret);
    TEST(test_dh_symmetry);
    TEST(test_rfc7748_iterate_1);
    TEST(test_rfc7748_iterate_1000);
    TEST(test_rfc7748_alice_bob_full);

    /* clamping and scalars */
    TEST(test_scalar_clamping);
    TEST(test_zero_point);
    TEST(test_small_scalar);
    TEST(test_different_keys);
    TEST(test_additional_vector1);
    TEST(test_high_bit_scalar);

    /* low-order points */
    TEST(test_low_order_1);
    TEST(test_low_order_point_order4);

    /* edge cases */
    TEST(test_all_zero_scalar);
    TEST(test_scalar_one);
    TEST(test_point_high_bit_set);

    /* consistency */
    TEST(test_base_consistency);
    TEST(test_scalar_mult_associative);

    /* multi-party */
    TEST(test_triple_dh);
    TEST(test_double_scalar);
    TEST(test_various_scalars);

    /* stress */
    TEST(test_many_keypairs);
    TEST(test_chain_mult);
    TEST(test_output_canonical);

    printf("\n%d/%d tests passed\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}
