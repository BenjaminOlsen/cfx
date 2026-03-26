/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#include "test_cli_common.h"

static void test_c_10_3(void) {
    char *argv[] = {"cfx_choose", "10", "3", NULL};
    char buf[8192];
    int rc = run_capture(cfx_choose_run, 3, argv, buf, sizeof(buf));
    CFX_ASSERT(rc == 0);
    assert_output_contains(buf, "120");
}

static void test_c_20_10(void) {
    char *argv[] = {"cfx_choose", "20", "10", NULL};
    char buf[8192];
    int rc = run_capture(cfx_choose_run, 3, argv, buf, sizeof(buf));
    CFX_ASSERT(rc == 0);
    assert_output_contains(buf, "184756");
}

static void test_c_5_0(void) {
    char *argv[] = {"cfx_choose", "5", "0", NULL};
    char buf[8192];
    int rc = run_capture(cfx_choose_run, 3, argv, buf, sizeof(buf));
    CFX_ASSERT(rc == 0);
    assert_output_contains(buf, "1");
}

/* C(n,n) = 1 */
static void test_c_n_n(void) {
    char *argv[] = {"cfx_choose", "7", "7", NULL};
    char buf[8192];
    int rc = run_capture(cfx_choose_run, 3, argv, buf, sizeof(buf));
    CFX_ASSERT(rc == 0);
    assert_output_contains(buf, "1");
}

/* symmetry: C(10,3) == C(10,7) */
static void test_symmetry(void) {
    char *argv1[] = {"cfx_choose", "10", "3", NULL};
    char *argv2[] = {"cfx_choose", "10", "7", NULL};
    char buf1[8192], buf2[8192];
    int rc1 = run_capture(cfx_choose_run, 3, argv1, buf1, sizeof(buf1));
    int rc2 = run_capture(cfx_choose_run, 3, argv2, buf2, sizeof(buf2));
    CFX_ASSERT(rc1 == 0 && rc2 == 0);
    assert_output_contains(buf1, "120");
    assert_output_contains(buf2, "120");
}

int main(void) {
    CFX_TEST(test_c_10_3);
    CFX_TEST(test_c_20_10);
    CFX_TEST(test_c_5_0);
    CFX_TEST(test_c_n_n);
    CFX_TEST(test_symmetry);
    printf("All cfx_choose CLI tests passed.\n");
    return 0;
}
