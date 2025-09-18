#ifndef CFX_MACROS_H
#define CFX_MACROS_H

#include <stdio.h>
#include <string.h>

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
    for (ssize_t i = na-1; i >= 0; --i) { \
        printf("%llu", A[i]); \
        if (i == 0) printf("]\n"); \
        else printf(", "); \
    } \
} while (0)

#endif /* CFX_MACROS_H */
