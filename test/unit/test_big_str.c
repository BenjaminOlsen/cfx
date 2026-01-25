/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */
/* test_big_str.c - Tests for big integer string parsing and conversion */

#include "test_common.h"

static void test_limb1(void) {
    cfx_limb_t L[] = { 123456789u };
    check_str_conversion("single limb", L, 1, "123456789");
    PRINT_TEST(1);
}

static void test_limb2(void) {
    cfx_limb_t L[] = { 123456789u, 4200u };
    check_str_conversion("two limbs pad", L, 2, "4200123456789");
    PRINT_TEST(1);
}

static void test_limb3(void) {
    cfx_limb_t L[] = { 1u, 1u };
    check_str_conversion("two limbs tiny low", L, 2, "1000000001");
    PRINT_TEST(1);
}

static void test_limb4(void) {
    cfx_limb_t L[] = { 999999999u, 999999999u };
    check_str_conversion("two limbs max", L, 2, "999999999999999999");
    PRINT_TEST(1);
}

static void test_limb5(void) {
    cfx_limb_t L[] = { 5u, 42u, 7u, 1u };
    check_str_conversion("four limbs pad", L, 4, "1000000007000000042000000005");
    PRINT_TEST(1);
}

static void test_limb6(void) {
    cfx_limb_t L[] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    check_str_conversion("9 limbs pad", L, 9,
        "900000000"
        "800000000"
        "700000000"
        "600000000"
        "500000000"
        "400000000"
        "300000000"
        "200000000"
        "1"
        );
    PRINT_TEST(1);
}

static void test_limb7(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    cfx_big_from_limb(&b, 1);
    for (int i = 0; i < 10; ++i) cfx_big_mul_sm_eq(&b, 1000000000u - 1u);
    char *s = cfx_big_dec_alloc(&b, NULL);
    CFX_ASSERT(s[0] == '9');
    CFX_ASSERT(strlen(s) >= 9);
    char *expect = "999999990000000044999999880000000209999999"
                   "748000000209999999880000000044999999990000000001";
    CFX_ASSERT(strcmp(s, expect) == 0);
    free(s);
    cfx_big_free(&b);
    PRINT_TEST(1);
}

static void test_str1(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    const char *sin =   "99911231231238761239876981273469128374169283476129"
                      "38471629384761250389172603459812630498672312387123"
                      "87123981723918273912891238719248719238719248169723"
                      "00091203901290909090911100091231283761000101023882";
    cfx_big_from_dec(&b, sin);
    char *sout = cfx_big_dec_alloc(&b, NULL);
    CFX_ASSERT(strcmp(sin, sout) == 0);
    free(sout);
    cfx_big_free(&b);
    PRINT_TEST(1);
}

static void test_str2(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    const char *sin = "9218";
    cfx_big_from_dec(&b, sin);
    char *sout = cfx_big_dec_alloc(&b, NULL);
    CFX_ASSERT(strcmp(sin, sout) == 0);
    free(sout);
    cfx_big_free(&b);
    PRINT_TEST(1);
}

static void test_from_str_matches1(void) {
    cfx_big_t b1, b2;
    cfx_big_init(&b1);
    cfx_big_init(&b2);

    CFX_ASSERT(cfx_big_from_str(&b1, "0xFF") == 0);
    CFX_ASSERT(cfx_big_from_str(&b2, "255") == 0);
    CFX_ASSERT(cfx_big_eq(&b1, &b2));

    cfx_big_free(&b1);
    cfx_big_free(&b2);
}

static void test_from_str_matches2(void) {
    cfx_big_t b1, b2, b3;
    cfx_big_init(&b1);
    cfx_big_init(&b2);
    cfx_big_init(&b3);

    CFX_ASSERT(cfx_big_from_str(&b1, "0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF") == 0);
    CFX_ASSERT(cfx_big_from_str(&b2, "1393796574908163946345982392040522594123775") == 0);
    CFX_ASSERT(cfx_big_from_str(&b3,
        "0b11111111111111111111111111111111111111111"
        "11111111111111111111111111111111111111111111"
        "1111111111111111111111111111111111111111111111111111111") == 0);

    CFX_ASSERT(cfx_big_eq(&b1, &b2));
    CFX_ASSERT(cfx_big_eq(&b1, &b3));

    cfx_big_free(&b1);
    cfx_big_free(&b2);
    cfx_big_free(&b3);
}

static void test_from_str_zero_forms(void) {
    cfx_big_t z1, z2, z3, z4;
    cfx_big_init(&z1);
    cfx_big_init(&z2);
    cfx_big_init(&z3);
    cfx_big_init(&z4);

    CFX_ASSERT(cfx_big_from_str(&z1, "0") == 0);
    CFX_ASSERT(cfx_big_from_str(&z2, "0000000") == 0);
    CFX_ASSERT(cfx_big_from_str(&z3, "0x0") == 0);
    CFX_ASSERT(cfx_big_from_str(&z4, "0b0") == 0);

    CFX_ASSERT(cfx_big_eq(&z1, &z2));
    CFX_ASSERT(cfx_big_eq(&z1, &z3));
    CFX_ASSERT(cfx_big_eq(&z1, &z4));

    cfx_big_free(&z1);
    cfx_big_free(&z2);
    cfx_big_free(&z3);
    cfx_big_free(&z4);
}

static void test_from_str_whitespace_ok(void) {
    cfx_big_t a, b, c;
    cfx_big_init(&a);
    cfx_big_init(&b);
    cfx_big_init(&c);

    CFX_ASSERT(cfx_big_from_str(&a, "   255") == 0);
    CFX_ASSERT(cfx_big_from_str(&b, "255   \n\t") == 0);
    CFX_ASSERT(cfx_big_from_str(&c, " \t  0xFF \r\n") == 0);

    CFX_ASSERT(cfx_big_eq(&a, &b));
    CFX_ASSERT(cfx_big_eq(&a, &c));

    cfx_big_free(&a);
    cfx_big_free(&b);
    cfx_big_free(&c);
}

static void test_from_str_hex_prefix_case(void) {
    cfx_big_t a, b, c;
    cfx_big_init(&a);
    cfx_big_init(&b);
    cfx_big_init(&c);

    CFX_ASSERT(cfx_big_from_str(&a, "0xdeadBEEF") == 0);
    CFX_ASSERT(cfx_big_from_str(&b, "0XDEADBEEF") == 0);
    CFX_ASSERT(cfx_big_from_str(&c, "3735928559") == 0);

    CFX_ASSERT(cfx_big_eq(&a, &b));
    CFX_ASSERT(cfx_big_eq(&a, &c));

    cfx_big_free(&a);
    cfx_big_free(&b);
    cfx_big_free(&c);
}

static void test_from_str_bin_prefix_case(void) {
    cfx_big_t a, b;
    cfx_big_init(&a);
    cfx_big_init(&b);

    CFX_ASSERT(cfx_big_from_str(&a, "0b101010") == 0);
    CFX_ASSERT(cfx_big_from_str(&b, "0B101010") == 0);

    CFX_ASSERT(cfx_big_eq(&a, &b));

    cfx_big_free(&a);
    cfx_big_free(&b);
}

static void test_from_str_rejects_trailing_junk(void) {
    cfx_big_t a;
    cfx_big_init(&a);

    CFX_ASSERT(cfx_big_from_str(&a, "255x") < 0);
    CFX_ASSERT(cfx_big_from_str(&a, "0xFFgg") < 0);
    CFX_ASSERT(cfx_big_from_str(&a, "0b10102") < 0);

    cfx_big_free(&a);
}

static void test_from_str_rejects_empty(void) {
    cfx_big_t a;
    cfx_big_init(&a);

    CFX_ASSERT(cfx_big_from_str(&a, "") < 0);
    CFX_ASSERT(cfx_big_from_str(&a, "   \t\r\n") < 0);

    cfx_big_free(&a);
}

static void test_from_str_rejects_prefix_only(void) {
    cfx_big_t a;
    cfx_big_init(&a);

    CFX_ASSERT(cfx_big_from_str(&a, "0x") < 0);
    CFX_ASSERT(cfx_big_from_str(&a, "0b") < 0);
    CFX_ASSERT(cfx_big_from_str(&a, "0X   ") < 0);

    cfx_big_free(&a);
}

static void test_from_str_no_legacy_octal(void) {
    cfx_big_t a, b;
    cfx_big_init(&a);
    cfx_big_init(&b);

    CFX_ASSERT(cfx_big_from_str(&a, "010") == 0);
    CFX_ASSERT(cfx_big_from_str(&b, "10") == 0);
    CFX_ASSERT(cfx_big_eq(&a, &b));

    cfx_big_free(&a);
    cfx_big_free(&b);
}

static void test_from_str_limb_boundary_all_ones(void) {
    cfx_big_t a, b;
    cfx_big_init(&a);
    cfx_big_init(&b);

    CFX_ASSERT(cfx_big_from_str(&a, "0xFFFFFFFFFFFFFFFF") == 0);
    CFX_ASSERT(cfx_big_from_str(&b, "18446744073709551615") == 0);
    CFX_ASSERT(cfx_big_eq(&a, &b));

    cfx_big_free(&a);
    cfx_big_free(&b);
}

static void test_from_str_limb_boundary_power_of_two(void) {
    cfx_big_t a, b;
    cfx_big_init(&a);
    cfx_big_init(&b);

    CFX_ASSERT(cfx_big_from_str(&a, "0x10000000000000000") == 0);
    CFX_ASSERT(cfx_big_from_str(&b, "18446744073709551616") == 0);
    CFX_ASSERT(cfx_big_eq(&a, &b));

    cfx_big_free(&a);
    cfx_big_free(&b);
}

static void test_scan_num_hex_basic(void) {
    cfx_big_t b, ref;
    cfx_big_init(&b);
    cfx_big_init(&ref);

    const char *s = "0xFFF+1";
    size_t n = 12345;

    CFX_ASSERT(cfx_big_scan_num_n(&b, (const uint8_t *)s, strlen(s), &n) == 0);
    CFX_ASSERT(n == 5);

    CFX_ASSERT(cfx_big_from_str(&ref, "4095") == 0);
    CFX_ASSERT(cfx_big_eq(&b, &ref));

    cfx_big_free(&b);
    cfx_big_free(&ref);
}

static void test_scan_num_bin_basic(void) {
    cfx_big_t b, ref;
    cfx_big_init(&b);
    cfx_big_init(&ref);

    const char *s = "0b0000000001+7";
    size_t n = 0;

    CFX_ASSERT(cfx_big_scan_num_n(&b, (const uint8_t *)s, strlen(s), &n) == 0);
    CFX_ASSERT(n == 12);

    CFX_ASSERT(cfx_big_from_str(&ref, "1") == 0);
    CFX_ASSERT(cfx_big_eq(&b, &ref));

    cfx_big_free(&b);
    cfx_big_free(&ref);
}

static void test_scan_num_b64_basic(void) {
    cfx_big_t b, ref;
    cfx_big_init(&b);
    cfx_big_init(&ref);

    const char *s = "b64:AQ==+7";
    size_t n = 0;

    CFX_ASSERT(cfx_big_scan_num_n(&b, (const uint8_t *)s, strlen(s), &n) == 0);
    CFX_ASSERT(n == 8);

    CFX_ASSERT(cfx_big_from_str(&ref, "1") == 0);
    CFX_ASSERT(cfx_big_eq(&b, &ref));

    cfx_big_free(&b);
    cfx_big_free(&ref);
}

static void test_scan_num_b64_whitespace_inside(void) {
    cfx_big_t b, ref;
    cfx_big_init(&b);
    cfx_big_init(&ref);

    const char *s = "b64: A Q = =   +1";
    size_t n = 12345;

    CFX_ASSERT(cfx_big_scan_num_n(&b, (const uint8_t *)s, strlen(s), &n) == 0);
    CFX_ASSERT(n == 15);
    CFX_ASSERT(cfx_big_from_str(&ref, "1") == 0);
    CFX_ASSERT(cfx_big_eq(&b, &ref));

    cfx_big_free(&b);
    cfx_big_free(&ref);
}

static void test_scan_num_b64_multibyte_value(void) {
    cfx_big_t b, ref;
    cfx_big_init(&b);
    cfx_big_init(&ref);

    const char *s = "b64:AQA=+1";
    size_t n = 0;

    CFX_ASSERT(cfx_big_scan_num_n(&b, (const uint8_t *)s, strlen(s), &n) == 0);
    CFX_ASSERT(n == 8);
    CFX_ASSERT(cfx_big_from_str(&ref, "256") == 0);
    CFX_ASSERT(cfx_big_eq(&b, &ref));

    cfx_big_free(&b);
    cfx_big_free(&ref);
}

static void test_scan_num_b64_stops_before_junk(void) {
    cfx_big_t b, ref;
    cfx_big_init(&b);
    cfx_big_init(&ref);

    const char *s = "b64:AQ==!+7";
    size_t n = 0;

    CFX_ASSERT(cfx_big_scan_num_n(&b, (const uint8_t *)s, strlen(s), &n) == 0);
    CFX_ASSERT(n == 8);

    CFX_ASSERT(cfx_big_from_str(&ref, "1") == 0);
    CFX_ASSERT(cfx_big_eq(&b, &ref));

    cfx_big_free(&b);
    cfx_big_free(&ref);
}

static void test_from_str_b64_allows_trailing_spaces_only(void) {
    cfx_big_t b, ref;
    cfx_big_init(&b);
    cfx_big_init(&ref);

    CFX_ASSERT(cfx_big_from_str(&b, "b64:AQ==   \r\n\t") == 0);
    CFX_ASSERT(cfx_big_from_str(&ref, "1") == 0);
    CFX_ASSERT(cfx_big_eq(&b, &ref));

    CFX_ASSERT(cfx_big_from_str(&b, "b64:AQ==XYZ") != 0);

    cfx_big_free(&b);
    cfx_big_free(&ref);
}

static void test_from_str_b64_rejects_invalid(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    CFX_ASSERT(cfx_big_from_str(&b, "b64:Zm9v*") != 0);
    CFX_ASSERT(cfx_big_from_str(&b, "b64:Zg=") != 0);
    CFX_ASSERT(cfx_big_from_str(&b, "b64:Zg=A") != 0);
    CFX_ASSERT(cfx_big_from_str(&b, "b64:   \n\t") != 0);
    cfx_big_free(&b);
}

static void test_scan_num_oct_basic(void) {
    /* octal scanning not fully implemented - placeholder */
    cfx_big_t b, ref;
    cfx_big_init(&b);
    cfx_big_init(&ref);
    cfx_big_free(&b);
    cfx_big_free(&ref);
}

static void test_scan_num_dec_basic(void) {
    cfx_big_t b, ref;
    cfx_big_init(&b);
    cfx_big_init(&ref);

    const char *s = "12345*9";
    size_t n = 0;

    CFX_ASSERT(cfx_big_scan_num_n(&b, (const uint8_t *)s, strlen(s), &n) == 0);
    CFX_ASSERT(n == 5);

    CFX_ASSERT(cfx_big_from_str(&ref, "12345") == 0);
    CFX_ASSERT(cfx_big_eq(&b, &ref));

    cfx_big_free(&b);
    cfx_big_free(&ref);
}

static void test_scan_num_stops_at_invalid_digit(void) {
    cfx_big_t b, ref;
    cfx_big_init(&b);
    cfx_big_init(&ref);

    const char *s = "0x12zz";
    size_t n = 0;

    CFX_ASSERT(cfx_big_scan_num_n(&b, (const uint8_t *)s, strlen(s), &n) == 0);
    CFX_ASSERT(n == 4);

    CFX_ASSERT(cfx_big_from_str(&ref, "18") == 0);
    CFX_ASSERT(cfx_big_eq(&b, &ref));

    cfx_big_free(&b);
    cfx_big_free(&ref);
}

static void test_scan_num_prefix_only_fails(void) {
    cfx_big_t b;
    cfx_big_init(&b);

    size_t n = 777;
    CFX_ASSERT(cfx_big_scan_num_n(&b, (const uint8_t *)"0x", 2, &n) < 0);
    CFX_ASSERT(n == 0);

    n = 777;
    CFX_ASSERT(cfx_big_scan_num_n(&b, (const uint8_t *)"0b", 2, &n) < 0);
    CFX_ASSERT(n == 0);

    cfx_big_free(&b);
}

static void test_scan_num_non_number_fails(void) {
    cfx_big_t b;
    cfx_big_init(&b);

    size_t n = 777;
    CFX_ASSERT(cfx_big_scan_num_n(&b, (const uint8_t *)"+123", 4, &n) < 0);
    CFX_ASSERT(n == 0);

    n = 777;
    CFX_ASSERT(cfx_big_scan_num_n(&b, (const uint8_t *)"(", 1, &n) < 0);
    CFX_ASSERT(n == 0);

    n = 777;
    CFX_ASSERT(cfx_big_scan_num_n(&b, (const uint8_t *)" 123", 4, &n) < 0);
    CFX_ASSERT(n == 0);

    cfx_big_free(&b);
}

static void test_scan_num_prefix_case_insensitive(void) {
    cfx_big_t a, b;
    cfx_big_init(&a);
    cfx_big_init(&b);

    size_t na = 0, nb = 0;

    CFX_ASSERT(cfx_big_scan_num_n(&a, (const uint8_t *)"0XFF ", 5, &na) == 0);
    CFX_ASSERT(na == 4);

    CFX_ASSERT(cfx_big_scan_num_n(&b, (const uint8_t *)"0B1010", 6, &nb) == 0);
    CFX_ASSERT(nb == 6);

    cfx_big_free(&a);
    cfx_big_free(&b);
}

static void test_scan_num_respects_in_len(void) {
    cfx_big_t b, ref;
    cfx_big_init(&b);
    cfx_big_init(&ref);

    const char *s = "0x1234+999";
    size_t n = 0;

    CFX_ASSERT(cfx_big_scan_num_n(&b, (const uint8_t *)s, 4, &n) == 0);
    CFX_ASSERT(n == 4);

    CFX_ASSERT(cfx_big_from_str(&ref, "18") == 0);
    CFX_ASSERT(cfx_big_eq(&b, &ref));

    cfx_big_free(&b);
    cfx_big_free(&ref);
}

static void test_to_str_buf_basic(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    cfx_big_from_str(&b, "12345678901234567890");

    char buf[64];
    size_t len = cfx_big_snprint_dec(&b, buf, sizeof(buf));
    CFX_ASSERT(len == 20);
    CFX_ASSERT(strcmp(buf, "12345678901234567890") == 0);

    cfx_big_free(&b);
}

static void test_to_str_buf_zero(void) {
    cfx_big_t b;
    cfx_big_init(&b);

    char buf[16];
    size_t len = cfx_big_snprint_dec(&b, buf, sizeof(buf));
    CFX_ASSERT(len == 1);
    CFX_ASSERT(strcmp(buf, "0") == 0);

    cfx_big_free(&b);
}

static void test_to_str_buf_truncate(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    cfx_big_from_str(&b, "12345678901234567890");

    char buf[10];
    size_t len = cfx_big_snprint_dec(&b, buf, sizeof(buf));
    CFX_ASSERT(len == 20);
    CFX_ASSERT(strlen(buf) == 9);
    CFX_ASSERT(strcmp(buf, "123456789") == 0);

    cfx_big_free(&b);
}

static void test_to_str_buf_null(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    cfx_big_from_str(&b, "999");

    size_t len = cfx_big_snprint_dec(&b, NULL, 0);
    CFX_ASSERT(len == 3);

    cfx_big_free(&b);
}

static void test_to_hex_buf_basic(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    cfx_big_from_hex(&b, "DEADBEEF");

    char buf[32];
    size_t len = cfx_big_snprint_hex(&b, buf, sizeof(buf));
    CFX_ASSERT(len == 8);
    CFX_ASSERT(strcmp(buf, "DEADBEEF") == 0 || strcmp(buf, "deadbeef") == 0);

    cfx_big_free(&b);
}

static void test_to_hex_buf_truncate(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    cfx_big_from_hex(&b, "DEADBEEFCAFE");

    char buf[5];
    size_t len = cfx_big_snprint_hex(&b, buf, sizeof(buf));
    CFX_ASSERT(len == 12);
    CFX_ASSERT(strlen(buf) == 4);

    cfx_big_free(&b);
}

static void test_to_bin_buf_basic(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    cfx_big_from_limb(&b, 13);

    char buf[32];
    size_t len = cfx_big_snprint_bin(&b, buf, sizeof(buf));
    CFX_ASSERT(len == 4);
    CFX_ASSERT(strcmp(buf, "1101") == 0);

    cfx_big_free(&b);
}

static void test_to_bin_buf_truncate(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    cfx_big_from_limb(&b, 255);

    char buf[5];
    size_t len = cfx_big_snprint_bin(&b, buf, sizeof(buf));
    CFX_ASSERT(len == 8);
    CFX_ASSERT(strlen(buf) == 4);
    CFX_ASSERT(strcmp(buf, "1111") == 0);

    cfx_big_free(&b);
}

int main(void) {
    CFX_TEST(test_limb1);
    CFX_TEST(test_limb2);
    CFX_TEST(test_limb3);
    CFX_TEST(test_limb4);
    CFX_TEST(test_limb5);
    CFX_TEST(test_limb6);
    CFX_TEST(test_limb7);
    CFX_TEST(test_str1);
    CFX_TEST(test_str2);
    CFX_TEST(test_from_str_matches1);
    CFX_TEST(test_from_str_matches2);
    CFX_TEST(test_from_str_zero_forms);
    CFX_TEST(test_from_str_whitespace_ok);
    CFX_TEST(test_from_str_hex_prefix_case);
    CFX_TEST(test_from_str_bin_prefix_case);
    CFX_TEST(test_from_str_rejects_trailing_junk);
    CFX_TEST(test_from_str_rejects_empty);
    CFX_TEST(test_from_str_rejects_prefix_only);
    CFX_TEST(test_from_str_no_legacy_octal);
    CFX_TEST(test_from_str_limb_boundary_all_ones);
    CFX_TEST(test_from_str_limb_boundary_power_of_two);
    CFX_TEST(test_scan_num_hex_basic);
    CFX_TEST(test_scan_num_bin_basic);
    CFX_TEST(test_scan_num_oct_basic);
    CFX_TEST(test_scan_num_dec_basic);
    CFX_TEST(test_scan_num_stops_at_invalid_digit);
    CFX_TEST(test_scan_num_prefix_only_fails);
    CFX_TEST(test_scan_num_non_number_fails);
    CFX_TEST(test_scan_num_prefix_case_insensitive);
    CFX_TEST(test_scan_num_respects_in_len);
    CFX_TEST(test_scan_num_b64_basic);
    CFX_TEST(test_scan_num_b64_whitespace_inside);
    CFX_TEST(test_scan_num_b64_multibyte_value);
    CFX_TEST(test_scan_num_b64_stops_before_junk);
    CFX_TEST(test_from_str_b64_allows_trailing_spaces_only);
    CFX_TEST(test_from_str_b64_rejects_invalid);
    CFX_TEST(test_to_str_buf_basic);
    CFX_TEST(test_to_str_buf_zero);
    CFX_TEST(test_to_str_buf_truncate);
    CFX_TEST(test_to_str_buf_null);
    CFX_TEST(test_to_hex_buf_basic);
    CFX_TEST(test_to_hex_buf_truncate);
    CFX_TEST(test_to_bin_buf_basic);
    CFX_TEST(test_to_bin_buf_truncate);
    puts("OK");
    return 0;
}
