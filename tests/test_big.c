#include "cfx/big.h"
#include <assert.h>

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

int main(void) {
    test_cfx_big_init();
    test_cfx_big_reserve();
    return 0;
}