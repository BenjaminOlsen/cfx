/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#include "test_cli_common.h"

static void test_sha256_string(void) {
    char *argv[] = {"cfx_hash", "-s", "hello", NULL};
    char buf[8192];
    int rc = run_capture(cfx_hash_run, 3, argv, buf, sizeof(buf));
    CFX_ASSERT(rc == 0);
    assert_output_starts_with(buf,
        "2cf24dba5fb0a30e26e83b2ac5b9e29e1b161e5c1fa7425e73043362938b9824");
}

static void test_sha256_empty(void) {
    char *argv[] = {"cfx_hash", "-s", "", NULL};
    char buf[8192];
    int rc = run_capture(cfx_hash_run, 3, argv, buf, sizeof(buf));
    CFX_ASSERT(rc == 0);
    assert_output_starts_with(buf,
        "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855");
}

static void test_sha3_256_string(void) {
    char *argv[] = {"cfx_hash", "--sha3-256", "-s", "hello", NULL};
    char buf[8192];
    int rc = run_capture(cfx_hash_run, 4, argv, buf, sizeof(buf));
    CFX_ASSERT(rc == 0);
    assert_output_starts_with(buf,
        "3338be694f50c5f338814986cdf0686453a888b84f424d792af4b9202398f392");
}

static void test_blake2b_string(void) {
    /* BLAKE2b-256 of "hello" */
    char *argv[] = {"cfx_hash", "--blake2b", "-n", "32", "-s", "hello", NULL};
    char buf[8192];
    int rc = run_capture(cfx_hash_run, 6, argv, buf, sizeof(buf));
    CFX_ASSERT(rc == 0);
    assert_output_starts_with(buf,
        "324dcf027dd4a30a932c441f365a25e86b173defa4b8e58948253471b81b72cf");
}

/* base64 output mode */
static void test_sha256_base64(void) {
    char *argv[] = {"cfx_hash", "-b64", "-s", "hello", NULL};
    char buf[8192];
    int rc = run_capture(cfx_hash_run, 4, argv, buf, sizeof(buf));
    CFX_ASSERT(rc == 0);
    /* base64 of the SHA-256 digest */
    assert_output_starts_with(buf, "LPJNul+wow4m6DsqxbninhsWHlwfp0JecwQzYpOLmCQ=");
}

/* SHAKE256 with explicit output length */
static void test_shake256_custom_len(void) {
    char *argv[] = {"cfx_hash", "--shake256", "-n", "16", "-s", "hello", NULL};
    char buf[8192];
    int rc = run_capture(cfx_hash_run, 6, argv, buf, sizeof(buf));
    CFX_ASSERT(rc == 0);
    /* output should be exactly 32 hex chars (16 bytes) + newline */
    CFX_ASSERT(strlen(buf) == 33);
}

/* multiple -s flags */
static void test_multiple_strings(void) {
    char *argv[] = {"cfx_hash", "-s", "hello", "-s", "world", NULL};
    char buf[8192];
    int rc = run_capture(cfx_hash_run, 5, argv, buf, sizeof(buf));
    CFX_ASSERT(rc == 0);
    /* should produce two lines of output */
    assert_output_contains(buf,
        "2cf24dba5fb0a30e26e83b2ac5b9e29e1b161e5c1fa7425e73043362938b9824");
    assert_output_contains(buf,
        "486ea46224d1bb4fb680f34f7c9ad96a8f24ec88be73ea8e5a6c65260e9cb8a7");
}

/* blake2b keyed hashing */
static void test_blake2b_keyed(void) {
    char *argv[] = {"cfx_hash", "--blake2b", "-n", "32",
        "-k", "deadbeef", "-s", "hello", NULL};
    char buf[8192];
    int rc = run_capture(cfx_hash_run, 8, argv, buf, sizeof(buf));
    CFX_ASSERT(rc == 0);
    /* keyed output should differ from unkeyed */
    CFX_ASSERT(strncmp(buf, "324dcf027dd4a30a", 16) != 0);
}

static void test_help(void) {
    char *argv[] = {"cfx_hash", "--help", NULL};
    char buf[8192];
    int rc = run_capture(cfx_hash_run, 2, argv, buf, sizeof(buf));
    CFX_ASSERT(rc == 0);
    assert_output_contains(buf, "Usage:");
}

int main(void) {
    CFX_TEST(test_sha256_string);
    CFX_TEST(test_sha256_empty);
    CFX_TEST(test_sha3_256_string);
    CFX_TEST(test_blake2b_string);
    CFX_TEST(test_sha256_base64);
    CFX_TEST(test_shake256_custom_len);
    CFX_TEST(test_multiple_strings);
    CFX_TEST(test_blake2b_keyed);
    CFX_TEST(test_help);
    printf("All cfx_hash CLI tests passed.\n");
    return 0;
}
