/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#include "test_cli_common.h"

static void test_2_pow_10(void) {
    char *argv[] = {"cfx_pow", "2", "10", NULL};
    char buf[8192];
    int rc = run_capture(cfx_pow_run, 3, argv, buf, sizeof(buf));
    CFX_ASSERT(rc == 0);
    assert_output_contains(buf, "1024");
}

static void test_3_pow_5(void) {
    char *argv[] = {"cfx_pow", "3", "5", NULL};
    char buf[8192];
    int rc = run_capture(cfx_pow_run, 3, argv, buf, sizeof(buf));
    CFX_ASSERT(rc == 0);
    assert_output_contains(buf, "243");
}

static void test_anything_pow_0(void) {
    char *argv[] = {"cfx_pow", "42", "0", NULL};
    char buf[8192];
    int rc = run_capture(cfx_pow_run, 3, argv, buf, sizeof(buf));
    CFX_ASSERT(rc == 0);
    assert_output_contains(buf, "1");
}

int main(void) {
    CFX_TEST(test_2_pow_10);
    CFX_TEST(test_3_pow_5);
    CFX_TEST(test_anything_pow_0);
    printf("All cfx_pow CLI tests passed.\n");
    return 0;
}
