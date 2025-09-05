#include "cfx/big.h"
#include "cfx/macros.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

static void test_cfx_big_init(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    assert(b.limb == NULL);
    assert(b.n == 0);
    assert(b.cap == 0);
}

static void test_cfx_big_reserve() {
    cfx_big_t b;
    size_t rcap1 = 55;
    cfx_big_reserve(&b, rcap1);
    assert(b.cap == rcap1);
    assert(b.n == 0);
    assert(b.limb != NULL);

    size_t rcap2 = rcap1 / 2;
    cfx_big_reserve(&b, rcap2);  /* shouldn't reserve less space. */
    assert(b.cap == rcap1);
}

static void big_init_from_limbs(cfx_big_t* b, const uint64_t *limbs, size_t n) {
    cfx_big_init(b);
    if (n == 0) return;
    b->limb = (uint64_t*)realloc(b->limb, n * sizeof(uint64_t));
    memcpy(b->limb, limbs, n * sizeof(uint64_t));
    b->n = n;
    b->cap = n;
}

static void test_mul_by_zero(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    cfx_big_set_val(&b, 123);
    cfx_big_mul_sm(&b, 2838);
    cfx_big_mul_sm(&b, 1928);
    cfx_big_mul_sm(&b, 9);
    cfx_big_mul_sm(&b, 123765);
    cfx_big_mul_sm(&b, 0);
    size_t sz = 0;
    char* s = cfx_big_to_str(&b, &sz);
    PRINT_DBG("s: %s, sz: %zu\n", s, sz);
    assert(sz == 1);
    assert(strcmp(s, "0") == 0);
    free(s);
}

static void test_add_sm(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    char *s = cfx_big_to_str(&b, NULL);
    PRINT_DBG("init: %s\n", s);
    
    cfx_big_set_val(&b, 123);
    s = cfx_big_to_str(&b, NULL);
    PRINT_DBG("after setting val: %s\n", s);

    cfx_big_add_sm(&b, 321);
    s = cfx_big_to_str(&b, NULL);
    PRINT_DBG("after add: %s\n",s);
    assert(b.limb[0] == 444);
    cfx_big_set_val(&b, CFX_BIG_BASE-1);
    s = cfx_big_to_str(&b, NULL);
    PRINT_DBG("after set: %s\n", s);
    free(s);
    assert(b.limb[0] == CFX_BIG_BASE-1);
    cfx_big_add_sm(&b, 1);
    assert(b.limb[0] == 0);
    assert(b.limb[1] == 1);
    s = cfx_big_to_str(&b, NULL);
    PRINT_DBG("after carry: %s\n", s);
    free(s);
}


static void test_sub_sm(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    PRINTF_BIG(&b, "init: b.n:%ld; ", b.n);

    cfx_big_set_val(&b, 1);
    cfx_big_sub_sm(&b, 1);
    assert(b.limb[0] == 0);
    cfx_big_sub_sm(&b, 1);
    assert(b.limb[0] == 0);
    cfx_big_set_val(&b, CFX_BIG_BASE-1);
    PRINTF_BIG(&b, "set to CFX_BIG_BASE-1: ");
    cfx_big_add_sm(&b, 1);
    PRINTF_BIG(&b, "add 1: ");
    cfx_big_sub_sm(&b, 1);
    PRINTF_BIG(&b, "sub 1: ");
    cfx_big_sub_sm(&b, 1);
    PRINTF_BIG(&b, "sub 1: ");
    cfx_big_sub_sm(&b, 1);
    PRINTF_BIG(&b, "sub 1: ");
    cfx_big_sub_sm(&b, 1);
    PRINTF_BIG(&b, "sub 1: ");
    
    cfx_big_add_sm(&b, 100);
    PRINTF_BIG(&b, "add 100: ");
    cfx_big_sub_sm(&b, 100);
    PRINTF_BIG(&b, "sub 100: ");
    cfx_big_add_sm(&b, 100);
    PRINTF_BIG(&b, "add 100: ");
    cfx_big_sub_sm(&b, 100);
    PRINTF_BIG(&b, "sub 100: ");
    cfx_big_add_sm(&b, 100);
    PRINTF_BIG(&b, "add 100: ");
    cfx_big_sub_sm(&b, 100);
    PRINTF_BIG(&b, "sub 100: ");

}

/* helper to run one test */
static void check(const char *label, const uint64_t *limbs, size_t n, const char *expect) {
    cfx_big_t b;
    big_init_from_limbs(&b, limbs, n);
    size_t len = 0;
    char *s = cfx_big_to_str(&b, &len);
    // PRINTF_BIG(&b, "str is:\n");
    // PRINT_DBG("str should be:\n%s\n", expect);
    assert(strcmp(s, expect) == 0);
    assert(len == strlen(expect));
    assert(s[len] == '\0');
    PRINT_DBG("[ok] %s -> %s\n", label, s);
    free(s);
    cfx_big_free(&b);
}


// Single small limb
static void test_limb1(void) {
    uint64_t L[] = { 123456789u };
    check("single limb", L, 1, "123456789"); 
}

// -->) Single 0 limb but n>0 (should normally be normalized away; if you allow it, behavior is “000…000”?)
// Prefer: represent zero as n==0. You can skip this if you enforce normalization.

// Two limbs with inner zero-padding needed:
// value = limb[1]*1e9 + limb[0] = 42*1e9 + 123456789
static void test_limb2(void) {
    uint64_t L[] = { 123456789u, 4200u };
    check("two limbs pad", L, 2, "4200123456789");
}

// Inner limb exact zero-padding boundary: limb[0] has fewer than 9 digits
static void test_limb3(void) {
    uint64_t L[] = { 1u, 1u };
    check("two limbs tiny low", L, 2, "1000000001");
}

// Max limb values
static void test_limb4(void) {
    uint64_t L[] = { 999999999u, 999999999u };
    check("two limbs max", L, 2, "999999999999999999");
}

// Four limbs mixed
// value = 1*1e27 + 7*1e18 + 42*1e9 + 5 → "1 000000007 000000042 000000005"
static void test_limb5(void) {
    uint64_t L[] = { 5u, 42u, 7u, 1u };
    check("four limbs pad", L, 4, "1000000007000000042000000005");
}

static void test_limb6(void) {
    uint64_t L[] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    check("9 limbs pad", L, 9,
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
}

// Large ndigits sanity: build via mul to exercise carry
static void test_limb7(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    cfx_big_set_val(&b, 1);
    for (int i = 0; i < 10; ++i) cfx_big_mul_sm(&b, 1000000000u - 1u); // (1e9-1)^10
    char *s = cfx_big_to_str(&b, NULL);
    // spot checks: starts with '9' and length >= 9
    PRINT_DBG("%s\n", s);
    assert(s[0] == '9');
    assert(strlen(s) >= 9);
    char* expect = "999999990000000044999999880000000209999999"
        "748000000209999999880000000044999999990000000001";
    assert(strcmp(s, expect) == 0);
    free(s);
    cfx_big_free(&b);
}

static void test_str1(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    const char *sin = "99911231231238761239876981273469128374169283476129384716293847612503891726034598126304986723123871238791238719248719238719248169723";
    cfx_big_from_str(&b, sin);
    char *sout = cfx_big_to_str(&b, NULL);
    int ok = (strcmp(sin, sout) == 0); 
    PRINT_DBG("test str1: in:\n%s \n%s\nout.. %s\n", sin, sout,
        ok ? "ok":"NOT ok");
    assert(ok);
}

static void test_str2(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    const char *sin = "9218";
    cfx_big_from_str(&b, sin);
    char *sout = cfx_big_to_str(&b, NULL);
    int ok = (strcmp(sin, sout) == 0); 
    PRINT_DBG("test str: in:\n%s \n%s\nout.. %s\n", sin, sout,
        ok ? "ok":"NOT ok");
    assert(ok);
}

int main(void) {
    test_cfx_big_init();
    test_cfx_big_reserve();
    test_mul_by_zero();
    test_limb1();
    test_limb2();
    test_limb3();
    test_limb4();
    test_limb5();
    test_limb6();
    test_limb7();
    test_str1();
    test_str2();
    test_add_sm();
    test_sub_sm();
    return 0;
}
