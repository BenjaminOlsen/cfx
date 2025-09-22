#ifndef CFX_MACROS_H
#define CFX_MACROS_H

#include "cfx/fmt.h"

#include <stdio.h>
#include <string.h>
#include <assert.h>
/* removes folders from __FILE__ path */
#define __FILENAME__ \
    (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__)

#ifdef CFX_DEBUG
#define CFX_PRINT_DBG(fmt, ...) \
    do { \
        fprintf(stderr, "[%s:%d:%s()]: " fmt, \
        __FILENAME__, __LINE__, __func__, ##__VA_ARGS__); \
    } while(0)
#else
#define CFX_PRINT_DBG(...) do {} while(0)
#endif

#define CFX_PRINT_ERR(fmt, ...) do {                       \
    char buf[128];                                               \
    snprintf(buf, sizeof(buf), fmt, ##__VA_ARGS__);              \
    fprintf(stderr, "[%s] - %s\n", __func__, buf);               \
} while (0)

#define U128_HI(x) (uint64_t)((x>>64)&0xFFFFFFFFFFFFFFFF)
#define U128_LO(x) (uint64_t)((x)&0xFFFFFFFFFFFFFFFF)

#define PRINT_TEST(ok) \
    do { \
        printf("[%s] ---- %s\n", __func__, ok ? "ok" : "NOT OK"); \
    } while (0)

#define PRINT_ARR(A, na) do { \
    printf("["); \
    for (size_t i = na; i--;) { \
        printf(""U64F"", A[i]); \
        if (i == 0) printf("]\n"); \
        else printf(", "); \
    } \
} while (0)

#define PRINT_BIG(s, b) \
do { \
    printf(".................\n"); \
    printf("%s, cap: %zu, n: %zu, limb %p\n", s, (b)->cap, (b)->n, (void*)(b)->limb); \
    for(size_t i=(b)->n; i--;) {printf("    limb[%zu]=0x%016llx\n", i, (b)->limb[i]);} \
    if((b)->n) printf("\n"); \
} while(0)

#define STR(x) #x
#define CFX_TEST(f) f(); printf(STR(f) "() - OK\n")

#define CFX_PRINT_TEST_COND(cond) printf("[%s] cond: '%s' - %s\n", __func__,  STR(cond), cond ? "OK" : ">>>> NOT OK <<<<")

#ifdef CFX_DEBUG
#define CFX_ASSERT_PRINT(cond) \
    do { \
        CFX_PRINT_TEST_COND(cond); \
        assert(cond); \
    } while (0)
#else
#define CFX_ASSERT_PRINT(cond) CFX_PRINT_TEST_COND(cond)
#endif 

#ifdef CFX_DEBUG
  #define PRINT_DBG(...) printf(__VA_ARGS__)
#else
  #define PRINT_DBG(...) do { } while (0)
#endif


#define CFX_ASSERT(cond) assert(cond)

#endif /* CFX_MACROS_H */
