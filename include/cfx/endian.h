#ifndef CFX_ENDIAN_H
#define CFX_ENDIAN_H

/* GCC / Clang */
#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
#define CFX_LITTLE_ENDIAN 1

/* Microsoft: always LE */
#elif defined(_MSC_VER)
#define CFX_LITTLE_ENDIAN 1

/* else use system headers if available */
#elif defined(__linux__) || defined(__ANDROID__)
#include <endian.h>
#if __BYTE_ORDER == __LITTLE_ENDIAN
#define CFX_LITTLE_ENDIAN 1
#endif

#elif defined(__APPLE__)
#include <machine/endian.h>
#if _BYTE_ORDER == _LITTLE_ENDIAN
#define CFX_LITTLE_ENDIAN 1
#endif

#else
#define CFX_LITTLE_ENDIAN 0
#endif

#endif /* CFX_ENDIAN_H */
