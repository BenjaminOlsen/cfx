/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#include "test_cli_common.h"

static void test_gcd_35_15(void) {
    char *argv[] = {"cfx_xgcd", "35", "15", NULL};
    char buf[8192];
    int rc = run_capture(cfx_xgcd_run, 3, argv, buf, sizeof(buf));
    CFX_ASSERT(rc == 0);
    assert_output_contains(buf, "5");
    assert_output_contains(buf, "OK");
}

static void test_gcd_coprime(void) {
    /* gcd(17, 13) = 1 */
    char *argv[] = {"cfx_xgcd", "17", "13", NULL};
    char buf[8192];
    int rc = run_capture(cfx_xgcd_run, 3, argv, buf, sizeof(buf));
    CFX_ASSERT(rc == 0);
    assert_output_contains(buf, "gcd");
    assert_output_contains(buf, "1");
    assert_output_contains(buf, "OK");
}

int main(void) {
    CFX_TEST(test_gcd_35_15);
    CFX_TEST(test_gcd_coprime);
    printf("All cfx_xgcd CLI tests passed.\n");
    return 0;
}
