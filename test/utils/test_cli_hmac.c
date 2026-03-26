/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#include "test_cli_common.h"

/* RFC 4231 test case 1: key = 20 bytes of 0x0b */
static void test_hmac_sha256_rfc4231_1(void) {
    char *argv[] = {"cfx_hmac", "-k",
        "0b0b0b0b0b0b0b0b0b0b0b0b0b0b0b0b0b0b0b0b",
        "-s", "Hi There", NULL};
    char buf[8192];
    int rc = run_capture(cfx_hmac_run, 5, argv, buf, sizeof(buf));
    CFX_ASSERT(rc == 0);
    assert_output_starts_with(buf,
        "b0344c61d8db38535ca8afceaf0bf12b881dc200c9833da726e9376c2e32cff7");
}

/* RFC 4231 test case 2: key = "Jefe" */
static void test_hmac_sha256_rfc4231_2(void) {
    char *argv[] = {"cfx_hmac", "-k", "4a656665",
        "-s", "what do ya want for nothing?", NULL};
    char buf[8192];
    int rc = run_capture(cfx_hmac_run, 5, argv, buf, sizeof(buf));
    CFX_ASSERT(rc == 0);
    assert_output_starts_with(buf,
        "5bdcc146bf60754e6a042426089575c75a003f089d2739839dec58b964ec3843");
}

static void test_hmac_sha512_rfc4231_1(void) {
    char *argv[] = {"cfx_hmac", "--sha512", "-k",
        "0b0b0b0b0b0b0b0b0b0b0b0b0b0b0b0b0b0b0b0b",
        "-s", "Hi There", NULL};
    char buf[8192];
    int rc = run_capture(cfx_hmac_run, 6, argv, buf, sizeof(buf));
    CFX_ASSERT(rc == 0);
    assert_output_starts_with(buf,
        "87aa7cdea5ef619d4ff0b4241a1d6cb0"
        "2379f4e2ce4ec2787ad0b30545e17cde"
        "daa833b7d6b8a702038b274eaea3f4e4"
        "be9d914eeb61f1702e696c203a126854");
}

static void test_hmac_no_key_fails(void) {
    char *argv[] = {"cfx_hmac", "-s", "hello", NULL};
    char buf[8192];
    int rc = run_capture(cfx_hmac_run, 3, argv, buf, sizeof(buf));
    CFX_ASSERT(rc != 0);
}

/* key > block size (65 bytes of 0xaa) — forces key hashing in HMAC-SHA256 */
static void test_hmac_sha256_long_key(void) {
    char *argv[] = {"cfx_hmac", "-k",
        "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
        "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
        "aaaa",  /* 66 bytes = 132 hex chars, > SHA-256 block size of 64 */
        "-s", "test", NULL};
    char buf[8192];
    int rc = run_capture(cfx_hmac_run, 5, argv, buf, sizeof(buf));
    CFX_ASSERT(rc == 0);
    /* just verify it succeeds — the long key path is exercised */
    CFX_ASSERT(strlen(buf) == 65);  /* 64 hex + newline */
}

/* 0x-prefixed key */
static void test_hmac_0x_prefix_key(void) {
    char *argv[] = {"cfx_hmac", "-k", "0x4a656665",
        "-s", "what do ya want for nothing?", NULL};
    char buf[8192];
    int rc = run_capture(cfx_hmac_run, 5, argv, buf, sizeof(buf));
    CFX_ASSERT(rc == 0);
    /* same as RFC 4231 case 2 — 0x prefix should be stripped */
    assert_output_starts_with(buf,
        "5bdcc146bf60754e6a042426089575c75a003f089d2739839dec58b964ec3843");
}

/* base64 output */
static void test_hmac_base64_output(void) {
    char *argv[] = {"cfx_hmac", "-b64", "-k", "4a656665",
        "-s", "what do ya want for nothing?", NULL};
    char buf[8192];
    int rc = run_capture(cfx_hmac_run, 6, argv, buf, sizeof(buf));
    CFX_ASSERT(rc == 0);
    /* should be base64, not hex */
    CFX_ASSERT(strncmp(buf, "5bdcc146", 8) != 0);
}

static void test_help(void) {
    char *argv[] = {"cfx_hmac", "--help", NULL};
    char buf[8192];
    int rc = run_capture(cfx_hmac_run, 2, argv, buf, sizeof(buf));
    CFX_ASSERT(rc == 0);
    assert_output_contains(buf, "Usage:");
}

int main(void) {
    CFX_TEST(test_hmac_sha256_rfc4231_1);
    CFX_TEST(test_hmac_sha256_rfc4231_2);
    CFX_TEST(test_hmac_sha512_rfc4231_1);
    CFX_TEST(test_hmac_sha256_long_key);
    CFX_TEST(test_hmac_0x_prefix_key);
    CFX_TEST(test_hmac_base64_output);
    CFX_TEST(test_hmac_no_key_fails);
    CFX_TEST(test_help);
    printf("All cfx_hmac CLI tests passed.\n");
    return 0;
}
