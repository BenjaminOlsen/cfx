/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#include "test_cli_common.h"

static void test_addition(void) {
    char *argv[] = {"cfx_dc", "2", "+", "3", NULL};
    char buf[8192];
    int rc = run_capture(cfx_dc_run, 4, argv, buf, sizeof(buf));
    CFX_ASSERT(rc == 0);
    assert_output_contains(buf, "5");
}

static void test_power(void) {
    char *argv[] = {"cfx_dc", "2", "**", "10", NULL};
    char buf[8192];
    int rc = run_capture(cfx_dc_run, 4, argv, buf, sizeof(buf));
    CFX_ASSERT(rc == 0);
    assert_output_contains(buf, "1024");
}

static void test_modulo(void) {
    char *argv[] = {"cfx_dc", "100", "%", "7", NULL};
    char buf[8192];
    int rc = run_capture(cfx_dc_run, 4, argv, buf, sizeof(buf));
    CFX_ASSERT(rc == 0);
    assert_output_contains(buf, "2");
}

static void test_multiplication(void) {
    char *argv[] = {"cfx_dc", "12", "*", "12", NULL};
    char buf[8192];
    int rc = run_capture(cfx_dc_run, 4, argv, buf, sizeof(buf));
    CFX_ASSERT(rc == 0);
    assert_output_contains(buf, "144");
}

static void test_hex_output(void) {
    char *argv[] = {"cfx_dc", "-ox", "255", NULL};
    char buf[8192];
    int rc = run_capture(cfx_dc_run, 3, argv, buf, sizeof(buf));
    CFX_ASSERT(rc == 0);
    assert_output_contains(buf, "ff");
}

/* subtraction producing zero */
static void test_subtraction_zero(void) {
    char *argv[] = {"cfx_dc", "42", "-", "42", NULL};
    char buf[8192];
    int rc = run_capture(cfx_dc_run, 4, argv, buf, sizeof(buf));
    CFX_ASSERT(rc == 0);
    assert_output_contains(buf, "0");
}

/* parentheses for precedence */
static void test_parentheses(void) {
    /* (2 + 3) * 4 = 20, not 2 + 12 = 14 */
    char *argv[] = {"cfx_dc", "(", "2", "+", "3", ")", "*", "4", NULL};
    char buf[8192];
    int rc = run_capture(cfx_dc_run, 8, argv, buf, sizeof(buf));
    CFX_ASSERT(rc == 0);
    assert_output_contains(buf, "20");
}

/* hex input + hex output */
static void test_hex_roundtrip(void) {
    char *argv[] = {"cfx_dc", "-ix", "-ox", "0xff", "+", "0x1", NULL};
    char buf[8192];
    int rc = run_capture(cfx_dc_run, 6, argv, buf, sizeof(buf));
    CFX_ASSERT(rc == 0);
    assert_output_contains(buf, "100");
}

/* RPN mode */
static void test_rpn(void) {
    char *argv[] = {"cfx_dc", "-rpn", "3", "4", "+", "5", "*", NULL};
    char buf[8192];
    int rc = run_capture(cfx_dc_run, 7, argv, buf, sizeof(buf));
    CFX_ASSERT(rc == 0);
    assert_output_contains(buf, "35");
}

/* bitwise shift */
static void test_shift_left(void) {
    char *argv[] = {"cfx_dc", "1", "<<", "10", NULL};
    char buf[8192];
    int rc = run_capture(cfx_dc_run, 4, argv, buf, sizeof(buf));
    CFX_ASSERT(rc == 0);
    assert_output_contains(buf, "1024");
}

int main(void) {
    CFX_TEST(test_addition);
    CFX_TEST(test_power);
    CFX_TEST(test_modulo);
    CFX_TEST(test_multiplication);
    CFX_TEST(test_hex_output);
    CFX_TEST(test_subtraction_zero);
    CFX_TEST(test_parentheses);
    CFX_TEST(test_hex_roundtrip);
    CFX_TEST(test_rpn);
    CFX_TEST(test_shift_left);
    printf("All cfx_dc CLI tests passed.\n");
    return 0;
}
