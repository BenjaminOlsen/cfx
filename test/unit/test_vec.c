/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */
/* test_vec.c - Tests for cfx_vec dynamic array */

#include "cfx/vec.h"
#include "cfx/macros.h"

#include <stdio.h>
#include <string.h>

static void test_vec_init(void) {
    cfx_vec_t v;
    cfx_vec_init(&v);

    CFX_ASSERT(v.data == NULL);
    CFX_ASSERT(v.size == 0);
    CFX_ASSERT(v.cap == 0);

    cfx_vec_free(&v);
}

static void test_vec_push_basic(void) {
    cfx_vec_t v;
    cfx_vec_init(&v);

    CFX_ASSERT(cfx_vec_push(&v, 42) == 0);
    CFX_ASSERT(v.size == 1);
    CFX_ASSERT(v.data[0] == 42);
    CFX_ASSERT(v.cap >= 1);

    CFX_ASSERT(cfx_vec_push(&v, 99) == 0);
    CFX_ASSERT(v.size == 2);
    CFX_ASSERT(v.data[0] == 42);
    CFX_ASSERT(v.data[1] == 99);

    cfx_vec_free(&v);
}

static void test_vec_push_many(void) {
    cfx_vec_t v;
    cfx_vec_init(&v);

    for (cfx_limb_t i = 0; i < 1000; ++i) {
        CFX_ASSERT(cfx_vec_push(&v, i * 2) == 0);
    }

    CFX_ASSERT(v.size == 1000);
    CFX_ASSERT(v.cap >= 1000);

    for (size_t i = 0; i < 1000; ++i) {
        CFX_ASSERT(v.data[i] == i * 2);
    }

    cfx_vec_free(&v);
}

static void test_vec_reserve(void) {
    cfx_vec_t v;
    cfx_vec_init(&v);

    CFX_ASSERT(cfx_vec_reserve(&v, 100) == 0);
    CFX_ASSERT(v.cap >= 100);
    CFX_ASSERT(v.size == 0);

    size_t old_cap = v.cap;
    CFX_ASSERT(cfx_vec_reserve(&v, 50) == 0);
    CFX_ASSERT(v.cap == old_cap);

    CFX_ASSERT(cfx_vec_reserve(&v, 1000) == 0);
    CFX_ASSERT(v.cap >= 1000);

    cfx_vec_free(&v);
}

static void test_vec_reserve_min_cap(void) {
    cfx_vec_t v;
    cfx_vec_init(&v);

    CFX_ASSERT(cfx_vec_reserve(&v, 1) == 0);
    CFX_ASSERT(v.cap >= CFX_VEC_MIN_CAP);

    cfx_vec_free(&v);
}

static void test_vec_resize_grow(void) {
    cfx_vec_t v;
    cfx_vec_init(&v);

    CFX_ASSERT(cfx_vec_resize(&v, 10, 0xFF) == 0);
    CFX_ASSERT(v.size == 10);
    CFX_ASSERT(v.cap >= 10);

    for (size_t i = 0; i < 10; ++i) {
        CFX_ASSERT(v.data[i] == 0xFF);
    }

    cfx_vec_free(&v);
}

static void test_vec_resize_shrink(void) {
    cfx_vec_t v;
    cfx_vec_init(&v);

    for (cfx_limb_t i = 0; i < 20; ++i) {
        cfx_vec_push(&v, i);
    }
    CFX_ASSERT(v.size == 20);

    CFX_ASSERT(cfx_vec_resize(&v, 5, 0) == 0);
    CFX_ASSERT(v.size == 5);

    for (size_t i = 0; i < 5; ++i) {
        CFX_ASSERT(v.data[i] == i);
    }

    cfx_vec_free(&v);
}

static void test_vec_resize_grow_preserves(void) {
    cfx_vec_t v;
    cfx_vec_init(&v);

    cfx_vec_push(&v, 1);
    cfx_vec_push(&v, 2);
    cfx_vec_push(&v, 3);

    CFX_ASSERT(cfx_vec_resize(&v, 6, 99) == 0);
    CFX_ASSERT(v.size == 6);

    CFX_ASSERT(v.data[0] == 1);
    CFX_ASSERT(v.data[1] == 2);
    CFX_ASSERT(v.data[2] == 3);
    CFX_ASSERT(v.data[3] == 99);
    CFX_ASSERT(v.data[4] == 99);
    CFX_ASSERT(v.data[5] == 99);

    cfx_vec_free(&v);
}

static void test_vec_clear(void) {
    cfx_vec_t v;
    cfx_vec_init(&v);

    for (cfx_limb_t i = 0; i < 50; ++i) {
        cfx_vec_push(&v, i);
    }

    size_t old_cap = v.cap;
    cfx_vec_clear(&v);

    CFX_ASSERT(v.size == 0);
    CFX_ASSERT(v.cap == old_cap);
    CFX_ASSERT(v.data != NULL);

    cfx_vec_free(&v);
}

static void test_vec_free(void) {
    cfx_vec_t v;
    cfx_vec_init(&v);

    for (cfx_limb_t i = 0; i < 100; ++i) {
        cfx_vec_push(&v, i);
    }

    cfx_vec_free(&v);

    CFX_ASSERT(v.data == NULL);
    CFX_ASSERT(v.size == 0);
    CFX_ASSERT(v.cap == 0);
}

static void test_vec_free_empty(void) {
    cfx_vec_t v;
    cfx_vec_init(&v);

    cfx_vec_free(&v);
    CFX_ASSERT(v.data == NULL);

    cfx_vec_free(&v);
    CFX_ASSERT(v.data == NULL);
}

static void test_vec_resize_zero(void) {
    cfx_vec_t v;
    cfx_vec_init(&v);

    CFX_ASSERT(cfx_vec_resize(&v, 0, 0) == 0);
    CFX_ASSERT(v.size == 0);

    cfx_vec_push(&v, 1);
    cfx_vec_push(&v, 2);
    CFX_ASSERT(cfx_vec_resize(&v, 0, 0) == 0);
    CFX_ASSERT(v.size == 0);
    CFX_ASSERT(v.cap > 0);

    cfx_vec_free(&v);
}

static void test_vec_capacity_doubling(void) {
    cfx_vec_t v;
    cfx_vec_init(&v);

    cfx_vec_reserve(&v, 1);
    size_t cap1 = v.cap;
    CFX_ASSERT(cap1 >= CFX_VEC_MIN_CAP);

    cfx_vec_reserve(&v, cap1 + 1);
    CFX_ASSERT(v.cap >= cap1 * 2);

    cfx_vec_free(&v);
}

int main(void) {
    CFX_TEST(test_vec_init);
    CFX_TEST(test_vec_push_basic);
    CFX_TEST(test_vec_push_many);
    CFX_TEST(test_vec_reserve);
    CFX_TEST(test_vec_reserve_min_cap);
    CFX_TEST(test_vec_resize_grow);
    CFX_TEST(test_vec_resize_shrink);
    CFX_TEST(test_vec_resize_grow_preserves);
    CFX_TEST(test_vec_clear);
    CFX_TEST(test_vec_free);
    CFX_TEST(test_vec_free_empty);
    CFX_TEST(test_vec_resize_zero);
    CFX_TEST(test_vec_capacity_doubling);
    puts("OK");
    return 0;
}
