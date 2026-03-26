/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#include "test_cli_common.h"

static void test_100_div_7(void) {
    /* 100 / 7 = 14 remainder 2 */
    char *argv[] = {"cfx_div", "100", "7", NULL};
    char buf[8192];
    int rc = run_capture(cfx_div_run, 3, argv, buf, sizeof(buf));
    CFX_ASSERT(rc == 0);
    assert_output_contains(buf, "14");
    assert_output_contains(buf, "2");
}

static void test_exact_division(void) {
    /* 144 / 12 = 12 remainder 0 */
    char *argv[] = {"cfx_div", "144", "12", NULL};
    char buf[8192];
    int rc = run_capture(cfx_div_run, 3, argv, buf, sizeof(buf));
    CFX_ASSERT(rc == 0);
    assert_output_contains(buf, "12");
}

int main(void) {
    CFX_TEST(test_100_div_7);
    CFX_TEST(test_exact_division);
    printf("All cfx_div CLI tests passed.\n");
    return 0;
}
