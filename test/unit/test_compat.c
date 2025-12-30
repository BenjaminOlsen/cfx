/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */
/* test_compat.c - Tests for cross-platform compatibility layer */

#include "cfx/compat.h"
#include "cfx/macros.h"

#include <stdio.h>
#include <string.h>

/* ---- Pure function tests ---- */

static void test_clz32(void) {
    CFX_ASSERT(cfx_clz32(0) == 32);
    CFX_ASSERT(cfx_clz32(1) == 31);
    CFX_ASSERT(cfx_clz32(2) == 30);
    CFX_ASSERT(cfx_clz32(3) == 30);
    CFX_ASSERT(cfx_clz32(0x80000000u) == 0);
    CFX_ASSERT(cfx_clz32(0x40000000u) == 1);
    CFX_ASSERT(cfx_clz32(0x00010000u) == 15);
    CFX_ASSERT(cfx_clz32(0x0000FFFFu) == 16);
    CFX_ASSERT(cfx_clz32(0xFFFFFFFFu) == 0);

    for (int i = 0; i < 32; ++i) {
        uint32_t v = 1u << i;
        CFX_ASSERT(cfx_clz32(v) == 31 - i);
    }
}

static void test_clz64(void) {
    CFX_ASSERT(cfx_clz64(0) == 64);
    CFX_ASSERT(cfx_clz64(1) == 63);
    CFX_ASSERT(cfx_clz64(2) == 62);
    CFX_ASSERT(cfx_clz64(0x8000000000000000ull) == 0);
    CFX_ASSERT(cfx_clz64(0x4000000000000000ull) == 1);
    CFX_ASSERT(cfx_clz64(0x0000000100000000ull) == 31);
    CFX_ASSERT(cfx_clz64(0x00000000FFFFFFFFull) == 32);
    CFX_ASSERT(cfx_clz64(0xFFFFFFFFFFFFFFFFull) == 0);

    for (int i = 0; i < 64; ++i) {
        uint64_t v = 1ull << i;
        CFX_ASSERT(cfx_clz64(v) == 63 - i);
    }
}

static void test_strndup_basic(void) {
    char* s = cfx_strndup("hello", 10);
    CFX_ASSERT(s != NULL);
    CFX_ASSERT(strcmp(s, "hello") == 0);
    free(s);
}

static void test_strndup_truncate(void) {
    char* s = cfx_strndup("hello world", 5);
    CFX_ASSERT(s != NULL);
    CFX_ASSERT(strcmp(s, "hello") == 0);
    CFX_ASSERT(strlen(s) == 5);
    free(s);
}

static void test_strndup_zero(void) {
    char* s = cfx_strndup("hello", 0);
    CFX_ASSERT(s != NULL);
    CFX_ASSERT(s[0] == '\0');
    CFX_ASSERT(strlen(s) == 0);
    free(s);
}

static void test_strndup_exact(void) {
    char* s = cfx_strndup("abc", 3);
    CFX_ASSERT(s != NULL);
    CFX_ASSERT(strcmp(s, "abc") == 0);
    free(s);
}

static void test_strndup_embedded_null(void) {
    /* strndup stops at first null regardless of n */
    char* s = cfx_strndup("ab\0cd", 10);
    CFX_ASSERT(s != NULL);
    CFX_ASSERT(strcmp(s, "ab") == 0);
    CFX_ASSERT(strlen(s) == 2);
    free(s);
}

/* ---- Sanity check tests ---- */

static void test_cpu_count(void) {
    int n = cfx_cpu_count();
    CFX_ASSERT(n >= 1);
    CFX_ASSERT(n < 10000);  /* sanity upper bound */
}

static void test_time_ns_monotonic(void) {
    uint64_t t1 = cfx_time_ns();
    uint64_t t2 = cfx_time_ns();
    uint64_t t3 = cfx_time_ns();

    CFX_ASSERT(t2 >= t1);
    CFX_ASSERT(t3 >= t2);
    CFX_ASSERT(t1 > 0);
}

static void test_time_ns_advances(void) {
    uint64_t t1 = cfx_time_ns();

    /* busy wait a tiny bit */
    volatile int x = 0;
    for (int i = 0; i < 100000; ++i) x += i;
    (void)x;

    uint64_t t2 = cfx_time_ns();
    CFX_ASSERT(t2 > t1);
}

/* ---- Threading tests ---- */

static void* thread_increment(void* arg) {
    int* counter = (int*)arg;
    (*counter) += 1;
    return NULL;
}

static void test_thread_create_join(void) {
    int counter = 0;
    cfx_thread_t thread;

    int rc = cfx_thread_create(&thread, thread_increment, &counter);
    CFX_ASSERT(rc == 0);

    rc = cfx_thread_join(thread, NULL);
    CFX_ASSERT(rc == 0);
    CFX_ASSERT(counter == 1);
}

static void test_thread_multiple(void) {
    int counters[4] = {0, 0, 0, 0};
    cfx_thread_t threads[4];

    for (int i = 0; i < 4; ++i) {
        int rc = cfx_thread_create(&threads[i], thread_increment, &counters[i]);
        CFX_ASSERT(rc == 0);
    }

    for (int i = 0; i < 4; ++i) {
        int rc = cfx_thread_join(threads[i], NULL);
        CFX_ASSERT(rc == 0);
        CFX_ASSERT(counters[i] == 1);
    }
}

/* ---- Mutex tests ---- */

typedef struct {
    cfx_mutex_t* mtx;
    int* counter;
    int iterations;
} mutex_test_arg_t;

static void* thread_mutex_increment(void* arg) {
    mutex_test_arg_t* mta = (mutex_test_arg_t*)arg;
    for (int i = 0; i < mta->iterations; ++i) {
        cfx_mutex_lock(mta->mtx);
        (*mta->counter)++;
        cfx_mutex_unlock(mta->mtx);
    }
    return NULL;
}

static void test_mutex_basic(void) {
    cfx_mutex_t mtx;
    cfx_mutex_init(&mtx);

    cfx_mutex_lock(&mtx);
    cfx_mutex_unlock(&mtx);

    cfx_mutex_destroy(&mtx);
}

static void test_mutex_contention(void) {
    cfx_mutex_t mtx;
    cfx_mutex_init(&mtx);

    int counter = 0;
    const int iterations = 10000;
    const int num_threads = 4;

    mutex_test_arg_t args[4];
    cfx_thread_t threads[4];

    for (int i = 0; i < num_threads; ++i) {
        args[i].mtx = &mtx;
        args[i].counter = &counter;
        args[i].iterations = iterations;
        cfx_thread_create(&threads[i], thread_mutex_increment, &args[i]);
    }

    for (int i = 0; i < num_threads; ++i) {
        cfx_thread_join(threads[i], NULL);
    }

    CFX_ASSERT(counter == num_threads * iterations);

    cfx_mutex_destroy(&mtx);
}

/* ---- Atomic tests ---- */

static void test_atomic_basic(void) {
    cfx_atomic_int val = 0;

    cfx_atomic_store(&val, 42);
    CFX_ASSERT(cfx_atomic_load(&val) == 42);

    cfx_atomic_store(&val, 0);
    CFX_ASSERT(cfx_atomic_load(&val) == 0);

    cfx_atomic_store(&val, -1);
    CFX_ASSERT(cfx_atomic_load(&val) == -1);
}

typedef struct {
    cfx_atomic_int* flag;
    int* shared_data;
} atomic_test_arg_t;

static void* thread_wait_for_flag(void* arg) {
    atomic_test_arg_t* ata = (atomic_test_arg_t*)arg;

    /* spin until flag is set */
    while (cfx_atomic_load(ata->flag) == 0) {
        /* busy wait */
    }

    /* flag was set, read shared data */
    return (void*)(intptr_t)(*ata->shared_data);
}

static void test_atomic_synchronization(void) {
    cfx_atomic_int flag = 0;
    int shared_data = 0;

    atomic_test_arg_t arg;
    arg.flag = &flag;
    arg.shared_data = &shared_data;

    cfx_thread_t thread;
    cfx_thread_create(&thread, thread_wait_for_flag, &arg);

    /* write data then set flag */
    shared_data = 12345;
    cfx_atomic_store(&flag, 1);

    cfx_thread_join(thread, NULL);

    /* thread should have seen shared_data = 12345 */
    CFX_ASSERT(shared_data == 12345);
}

int main(void) {
    /* pure functions */
    CFX_TEST(test_clz32);
    CFX_TEST(test_clz64);
    CFX_TEST(test_strndup_basic);
    CFX_TEST(test_strndup_truncate);
    CFX_TEST(test_strndup_zero);
    CFX_TEST(test_strndup_exact);
    CFX_TEST(test_strndup_embedded_null);

    CFX_TEST(test_cpu_count);
    CFX_TEST(test_time_ns_monotonic);
    CFX_TEST(test_time_ns_advances);

    CFX_TEST(test_thread_create_join);
    CFX_TEST(test_thread_multiple);

    CFX_TEST(test_mutex_basic);
    CFX_TEST(test_mutex_contention);

    CFX_TEST(test_atomic_basic);
    CFX_TEST(test_atomic_synchronization);

    puts("OK");
    return 0;
}
