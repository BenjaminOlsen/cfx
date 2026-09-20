/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#include "test_common.h"
#include "cfx_utils_common.h"

static void test_parse_u32(void) {
    const struct {
        const char *text;
        uint32_t value;
    } valid[] = {
        {"0", 0},
        {"+0", 0},
        {"42", 42},
        {" \t+42", 42},
        {"052", 42},
        {"0x2a", 42},
        {"0X2A", 42},
        {"4294967294", UINT32_MAX - 1},
        {"4294967295", UINT32_MAX},
        {"0xffffffff", UINT32_MAX},
        {"037777777777", UINT32_MAX},
        {"00000000000000000001", 1},
        {"-0", 0}, {"-1", UINT32_MAX},
        {"-2", UINT32_MAX - 1},
        {" \t-0x1", UINT32_MAX},
        {"-01", UINT32_MAX},
        {"-2147483648", UINT32_C(2147483648)},
        {"-4294967295", 1},
        {"-0xffffffff", 1},
        {"-037777777777", 1}
    };

    const char *invalid[] = {
        "", " \t", "+", "-", "0x", "+0X", "08",
        "0xg", "12a", "1 ", "1\n", "++1", "+ 1",
        "4294967296", "0x100000000", "040000000000",
        "-4294967296", "-0x100000000", "-040000000000",
        "-18446744073709551616", "--1", "-+1", "- 1", "-0x", "-08",
        "18446744073709551616", "9999999999999999999999999999999999", "blah",
        NULL
    };

    for (size_t i = 0; i < sizeof(valid) / sizeof(valid[0]); ++i) {
        uint32_t counter = 123;
        CFX_ASSERT(cfx_parse_u32(valid[i].text, &counter) == 0);
        CFX_ASSERT(counter == valid[i].value);
    }
    for (size_t i = 0; i < sizeof(invalid) / sizeof(invalid[0]); ++i) {
        uint32_t counter = 123;
        CFX_ASSERT(cfx_parse_u32(invalid[i], &counter) == -1);
        CFX_ASSERT(counter == 123);
    }
    CFX_ASSERT(cfx_parse_u32("1", NULL) == -1);
}

int main(void) {
    CFX_TEST(test_parse_u32);
    return 0;
}

