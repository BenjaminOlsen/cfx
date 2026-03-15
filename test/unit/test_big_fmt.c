#include "test_common.h"

#include <math.h>

/* helper: wrap a string as a FILE* for cfx_big_from_file */
static FILE *str_to_file(const char *s) {
    FILE *f = tmpfile();
    if (!f) return NULL;
    fwrite(s, 1, strlen(s), f);
    rewind(f);
    return f;
}

/* ---- cfx_big_from_file -------------------------------------------------- */

static void test_from_file_decimal(void) {
    FILE *f = str_to_file("12345");
    cfx_big_t b;
    cfx_big_init(&b);
    CFX_ASSERT(cfx_big_from_file(&b, f, 10) == 0);

    cfx_big_t expected = make_u64(12345);
    CFX_ASSERT(cfx_big_eq(&b, &expected));

    cfx_big_free(&b);
    cfx_big_free(&expected);
    fclose(f);
}

static void test_from_file_hex_prefix(void) {
    FILE *f = str_to_file("0xFF");
    cfx_big_t b;
    cfx_big_init(&b);
    CFX_ASSERT(cfx_big_from_file(&b, f, 16) == 0);

    cfx_big_t expected = make_u64(255);
    CFX_ASSERT(cfx_big_eq(&b, &expected));

    cfx_big_free(&b);
    cfx_big_free(&expected);
    fclose(f);
}

static void test_from_file_hex_no_prefix(void) {
    FILE *f = str_to_file("DEADBEEF");
    cfx_big_t b;
    cfx_big_init(&b);
    CFX_ASSERT(cfx_big_from_file(&b, f, 16) == 0);

    cfx_big_t expected = make_u64(0xDEADBEEF);
    CFX_ASSERT(cfx_big_eq(&b, &expected));

    cfx_big_free(&b);
    cfx_big_free(&expected);
    fclose(f);
}

static void test_from_file_binary(void) {
    FILE *f = str_to_file("0b1010");
    cfx_big_t b;
    cfx_big_init(&b);
    CFX_ASSERT(cfx_big_from_file(&b, f, 2) == 0);

    cfx_big_t expected = make_u64(10);
    CFX_ASSERT(cfx_big_eq(&b, &expected));

    cfx_big_free(&b);
    cfx_big_free(&expected);
    fclose(f);
}

static void test_from_file_separators(void) {
    FILE *f = str_to_file("1_000_000");
    cfx_big_t b;
    cfx_big_init(&b);
    CFX_ASSERT(cfx_big_from_file(&b, f, 10) == 0);

    cfx_big_t expected = make_u64(1000000);
    CFX_ASSERT(cfx_big_eq(&b, &expected));

    cfx_big_free(&b);
    cfx_big_free(&expected);
    fclose(f);
}

static void test_from_file_hex_separators(void) {
    FILE *f = str_to_file("DEAD_BEEF");
    cfx_big_t b;
    cfx_big_init(&b);
    CFX_ASSERT(cfx_big_from_file(&b, f, 16) == 0);

    cfx_big_t expected = make_u64(0xDEADBEEF);
    CFX_ASSERT(cfx_big_eq(&b, &expected));

    cfx_big_free(&b);
    cfx_big_free(&expected);
    fclose(f);
}

static void test_from_file_leading_space(void) {
    FILE *f = str_to_file("   42");
    cfx_big_t b;
    cfx_big_init(&b);
    CFX_ASSERT(cfx_big_from_file(&b, f, 10) == 0);

    cfx_big_t expected = make_u64(42);
    CFX_ASSERT(cfx_big_eq(&b, &expected));

    cfx_big_free(&b);
    cfx_big_free(&expected);
    fclose(f);
}

static void test_from_file_plus_sign(void) {
    FILE *f = str_to_file("+99");
    cfx_big_t b;
    cfx_big_init(&b);
    CFX_ASSERT(cfx_big_from_file(&b, f, 10) == 0);

    cfx_big_t expected = make_u64(99);
    CFX_ASSERT(cfx_big_eq(&b, &expected));

    cfx_big_free(&b);
    cfx_big_free(&expected);
    fclose(f);
}

static void test_from_file_zero(void) {
    FILE *f = str_to_file("0");
    cfx_big_t b;
    cfx_big_init(&b);
    CFX_ASSERT(cfx_big_from_file(&b, f, 10) == 0);

    cfx_big_t expected = make_u64(0);
    CFX_ASSERT(cfx_big_eq(&b, &expected));

    cfx_big_free(&b);
    cfx_big_free(&expected);
    fclose(f);
}

static void test_from_file_negative(void) {
    FILE *f = str_to_file("-42");
    cfx_big_t b;
    cfx_big_init(&b);
    CFX_ASSERT(cfx_big_from_file(&b, f, 10) != 0);
    cfx_big_free(&b);
    fclose(f);
}

static void test_from_file_invalid_char(void) {
    FILE *f = str_to_file("12z34");
    cfx_big_t b;
    cfx_big_init(&b);
    CFX_ASSERT(cfx_big_from_file(&b, f, 10) != 0);
    cfx_big_free(&b);
    fclose(f);
}

static void test_from_file_empty(void) {
    FILE *f = str_to_file("");
    cfx_big_t b;
    cfx_big_init(&b);
    CFX_ASSERT(cfx_big_from_file(&b, f, 10) != 0);
    cfx_big_free(&b);
    fclose(f);
}

static void test_from_file_large_decimal(void) {
    /* 2^64 = 18446744073709551616 (20 digits, spans multiple chunks) */
    FILE *f = str_to_file("18446744073709551616");
    cfx_big_t b;
    cfx_big_init(&b);
    CFX_ASSERT(cfx_big_from_file(&b, f, 10) == 0);

    cfx_big_t expected;
    cfx_big_init(&expected);
    cfx_big_t one;
    cfx_big_init(&one);
    cfx_big_from_limb(&one, 1);
    cfx_big_shl_bits(&expected, &one, 64);

    CFX_ASSERT(cfx_big_eq(&b, &expected));

    cfx_big_free(&b);
    cfx_big_free(&expected);
    cfx_big_free(&one);
    fclose(f);
}

static void test_from_file_hex_large(void) {
    FILE *f = str_to_file("0x1FFFFFFFFFFFFFFFF");
    cfx_big_t b;
    cfx_big_init(&b);
    CFX_ASSERT(cfx_big_from_file(&b, f, 16) == 0);
    ASSERT_BIG_HEX(&b, "1FFFFFFFFFFFFFFFF");
    cfx_big_free(&b);
    fclose(f);
}

static void test_from_file_binary_large(void) {
    /* 0b11111111 = 255 */
    FILE *f = str_to_file("0b11111111");
    cfx_big_t b;
    cfx_big_init(&b);
    CFX_ASSERT(cfx_big_from_file(&b, f, 2) == 0);

    cfx_big_t expected = make_u64(255);
    CFX_ASSERT(cfx_big_eq(&b, &expected));

    cfx_big_free(&b);
    cfx_big_free(&expected);
    fclose(f);
}

/* ---- cfx_big_log -------------------------------------------------------- */

static void test_big_log_base10(void) {
    cfx_big_t b = make_u64(1000);
    double l = cfx_big_log(&b, 10.0);
    CFX_ASSERT(fabs(l - 3.0) < 0.001);
    cfx_big_free(&b);
}

static void test_big_log_base2(void) {
    cfx_big_t b = make_u64(1024);
    double l = cfx_big_log(&b, 2.0);
    CFX_ASSERT(fabs(l - 10.0) < 0.001);
    cfx_big_free(&b);
}

static void test_big_log_large(void) {
    /* 2^128 */
    cfx_big_t b;
    cfx_big_init(&b);
    cfx_big_t one;
    cfx_big_init(&one);
    cfx_big_from_limb(&one, 1);
    cfx_big_shl_bits(&b, &one, 128);

    double l = cfx_big_log(&b, 2.0);
    CFX_ASSERT(fabs(l - 128.0) < 0.1);

    cfx_big_free(&b);
    cfx_big_free(&one);
}

static void test_big_log_one(void) {
    cfx_big_t b = make_u64(1);
    double l = cfx_big_log(&b, 10.0);
    CFX_ASSERT(fabs(l) < 0.001);
    cfx_big_free(&b);
}

/* ---- cfx_big_to_sci ----------------------------------------------------- */

static void test_to_sci_zero(void) {
    cfx_big_t b = make_u64(0);
    char buf[64];
    int rc = cfx_big_to_sci(&b, 10, 3, buf, sizeof(buf));
    CFX_ASSERT(rc == 1);
    CFX_ASSERT(strcmp(buf, "0") == 0);
    cfx_big_free(&b);
}

static void test_to_sci_small(void) {
    cfx_big_t b = make_u64(12345);
    char buf[64];
    int rc = cfx_big_to_sci(&b, 10, 5, buf, sizeof(buf));
    CFX_ASSERT(rc == 1);
    /* ~1.2345e4 */
    CFX_ASSERT(strncmp(buf, "1.234", 5) == 0);
    CFX_ASSERT(strstr(buf, "e4") != NULL);
    cfx_big_free(&b);
}

static void test_to_sci_large(void) {
    /* 10^100 (googol) — fp rounding may give 9.999e99 or 1.000e100 */
    cfx_big_t b;
    cfx_big_init(&b);
    cfx_big_pow_sm(&b, 10, 100);

    char buf[128];
    int rc = cfx_big_to_sci(&b, 10, 4, buf, sizeof(buf));
    CFX_ASSERT(rc == 1);
    /* accept either rounding */
    CFX_ASSERT(strstr(buf, "e99") != NULL || strstr(buf, "e100") != NULL);

    cfx_big_free(&b);
}

static void test_to_sci_base2(void) {
    /* 2^64 */
    cfx_big_t b;
    cfx_big_init(&b);
    cfx_big_t one;
    cfx_big_init(&one);
    cfx_big_from_limb(&one, 1);
    cfx_big_shl_bits(&b, &one, 64);

    char buf[128];
    int rc = cfx_big_to_sci(&b, 2, 3, buf, sizeof(buf));
    CFX_ASSERT(rc == 1);
    /* non-base-10 format: "X.XX x 2^64" */
    CFX_ASSERT(strstr(buf, "x 2^64") != NULL);

    cfx_big_free(&b);
    cfx_big_free(&one);
}

static void test_to_sci_invalid_args(void) {
    char buf[64];
    /* NULL x */
    CFX_ASSERT(cfx_big_to_sci(NULL, 10, 3, buf, sizeof(buf)) == 0);
    /* NULL out */
    cfx_big_t b = make_u64(42);
    CFX_ASSERT(cfx_big_to_sci(&b, 10, 3, NULL, 64) == 0);
    /* zero outsz */
    CFX_ASSERT(cfx_big_to_sci(&b, 10, 3, buf, 0) == 0);
    /* base < 2 */
    CFX_ASSERT(cfx_big_to_sci(&b, 1, 3, buf, sizeof(buf)) == 0);
    cfx_big_free(&b);
}

static void test_to_sci_sig_digits_1(void) {
    cfx_big_t b = make_u64(99999);
    char buf[64];
    int rc = cfx_big_to_sci(&b, 10, 1, buf, sizeof(buf));
    CFX_ASSERT(rc == 1);
    /* ~1e5 (one significant digit) */
    CFX_ASSERT(strstr(buf, "e") != NULL);
    cfx_big_free(&b);
}

/* ---- main --------------------------------------------------------------- */

int main(void) {
    CFX_TEST(test_from_file_decimal);
    CFX_TEST(test_from_file_hex_prefix);
    CFX_TEST(test_from_file_hex_no_prefix);
    CFX_TEST(test_from_file_binary);
    CFX_TEST(test_from_file_separators);
    CFX_TEST(test_from_file_hex_separators);
    CFX_TEST(test_from_file_leading_space);
    CFX_TEST(test_from_file_plus_sign);
    CFX_TEST(test_from_file_zero);
    CFX_TEST(test_from_file_negative);
    CFX_TEST(test_from_file_invalid_char);
    CFX_TEST(test_from_file_empty);
    CFX_TEST(test_from_file_large_decimal);
    CFX_TEST(test_from_file_hex_large);
    CFX_TEST(test_from_file_binary_large);

    CFX_TEST(test_big_log_base10);
    CFX_TEST(test_big_log_base2);
    CFX_TEST(test_big_log_large);
    CFX_TEST(test_big_log_one);

    CFX_TEST(test_to_sci_zero);
    CFX_TEST(test_to_sci_small);
    CFX_TEST(test_to_sci_large);
    CFX_TEST(test_to_sci_base2);
    CFX_TEST(test_to_sci_invalid_args);
    CFX_TEST(test_to_sci_sig_digits_1);

    puts("OK");
    return 0;
}
