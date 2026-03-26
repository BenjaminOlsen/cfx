/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#include "test_cli_common.h"

static void test_prime_7(void) {
    char *argv[] = {"cfx_isprime", "7", NULL};
    char buf[8192];
    int rc = run_capture(cfx_isprime_run, 2, argv, buf, sizeof(buf));
    CFX_ASSERT(rc == 0);
    assert_output_contains(buf, "prime");
}

static void test_prime_97(void) {
    char *argv[] = {"cfx_isprime", "97", NULL};
    char buf[8192];
    int rc = run_capture(cfx_isprime_run, 2, argv, buf, sizeof(buf));
    CFX_ASSERT(rc == 0);
    assert_output_contains(buf, "prime");
}

static void test_composite_4(void) {
    char *argv[] = {"cfx_isprime", "4", NULL};
    char buf[8192];
    int rc = run_capture(cfx_isprime_run, 2, argv, buf, sizeof(buf));
    CFX_ASSERT(rc == 1);
    assert_output_contains(buf, "composite");
}

static void test_composite_100(void) {
    char *argv[] = {"cfx_isprime", "100", NULL};
    char buf[8192];
    int rc = run_capture(cfx_isprime_run, 2, argv, buf, sizeof(buf));
    CFX_ASSERT(rc == 1);
    assert_output_contains(buf, "composite");
}

/* 2 is prime (smallest prime) */
static void test_prime_2(void) {
    char *argv[] = {"cfx_isprime", "2", NULL};
    char buf[8192];
    int rc = run_capture(cfx_isprime_run, 2, argv, buf, sizeof(buf));
    CFX_ASSERT(rc == 0);
}

/* 1 is not prime */
static void test_not_prime_1(void) {
    char *argv[] = {"cfx_isprime", "1", NULL};
    char buf[8192];
    int rc = run_capture(cfx_isprime_run, 2, argv, buf, sizeof(buf));
    CFX_ASSERT(rc == 1);
}

/* quiet mode — no output, just return code */
static void test_quiet_mode(void) {
    char *argv[] = {"cfx_isprime", "-q", "7", NULL};
    char buf[8192];
    int rc = run_capture(cfx_isprime_run, 3, argv, buf, sizeof(buf));
    CFX_ASSERT(rc == 0);
    CFX_ASSERT(buf[0] == '\0');  /* no output in quiet mode */
}

/* hex input */
static void test_hex_input(void) {
    /* 0x61 = 97, which is prime */
    char *argv[] = {"cfx_isprime", "-x", "61", NULL};
    char buf[8192];
    int rc = run_capture(cfx_isprime_run, 3, argv, buf, sizeof(buf));
    CFX_ASSERT(rc == 0);
}

int main(void) {
    CFX_TEST(test_prime_7);
    CFX_TEST(test_prime_97);
    CFX_TEST(test_composite_4);
    CFX_TEST(test_composite_100);
    CFX_TEST(test_prime_2);
    CFX_TEST(test_not_prime_1);
    CFX_TEST(test_quiet_mode);
    CFX_TEST(test_hex_input);
    printf("All cfx_isprime CLI tests passed.\n");
    return 0;
}
