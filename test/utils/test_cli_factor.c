/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#include "test_cli_common.h"

static void test_factor_12(void) {
    char *argv[] = {"cfx_factor", "12", NULL};
    char buf[8192];
    int rc = run_capture(cfx_factor_run, 2, argv, buf, sizeof(buf));
    CFX_ASSERT(rc == 0);
    assert_output_contains(buf, "2");
    assert_output_contains(buf, "3");
}

static void test_factor_prime(void) {
    char *argv[] = {"cfx_factor", "97", NULL};
    char buf[8192];
    int rc = run_capture(cfx_factor_run, 2, argv, buf, sizeof(buf));
    CFX_ASSERT(rc == 0);
    assert_output_contains(buf, "97");
}

static void test_factor_1001(void) {
    /* 1001 = 7 * 11 * 13 */
    char *argv[] = {"cfx_factor", "1001", NULL};
    char buf[8192];
    int rc = run_capture(cfx_factor_run, 2, argv, buf, sizeof(buf));
    CFX_ASSERT(rc == 0);
    assert_output_contains(buf, "7");
    assert_output_contains(buf, "11");
    assert_output_contains(buf, "13");
}

static void test_help(void) {
    char *argv[] = {"cfx_factor", "--help", NULL};
    char buf[8192];
    int rc = run_capture(cfx_factor_run, 2, argv, buf, sizeof(buf));
    CFX_ASSERT(rc == 0);
}

/* power of 2 */
static void test_factor_power_of_2(void) {
    char *argv[] = {"cfx_factor", "256", NULL};
    char buf[8192];
    int rc = run_capture(cfx_factor_run, 2, argv, buf, sizeof(buf));
    CFX_ASSERT(rc == 0);
    assert_output_contains(buf, "2");
}

/* factor 2 (smallest prime) */
static void test_factor_2(void) {
    char *argv[] = {"cfx_factor", "2", NULL};
    char buf[8192];
    int rc = run_capture(cfx_factor_run, 2, argv, buf, sizeof(buf));
    CFX_ASSERT(rc == 0);
    assert_output_contains(buf, "2");
}

/* hex input */
static void test_factor_hex(void) {
    /* 0xff = 255 = 3 * 5 * 17 */
    char *argv[] = {"cfx_factor", "-x", "ff", NULL};
    char buf[8192];
    int rc = run_capture(cfx_factor_run, 3, argv, buf, sizeof(buf));
    CFX_ASSERT(rc == 0);
    assert_output_contains(buf, "3");
    assert_output_contains(buf, "5");
    assert_output_contains(buf, "17");
}

int main(void) {
    CFX_TEST(test_factor_12);
    CFX_TEST(test_factor_prime);
    CFX_TEST(test_factor_1001);
    CFX_TEST(test_factor_power_of_2);
    CFX_TEST(test_factor_2);
    CFX_TEST(test_factor_hex);
    CFX_TEST(test_help);
    printf("All cfx_factor CLI tests passed.\n");
    return 0;
}
