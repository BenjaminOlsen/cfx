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
        fprintf(stderr, "[%s:%d:%s():]: " fmt, \
        __FILENAME__, __LINE__, __func__, ##__VA_ARGS__); \
    } while(0)
#else
#define CFX_PRINT_DBG(...) do {} while(0)
#endif

#endif /* CFX_MACROS_H */
