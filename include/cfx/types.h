#ifndef CFX_TYPES_H
#define CFX_TYPES_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef __uint128_t uint128_t;

typedef struct {
    uint64_t hi, lo;
} cfx_u128;

#ifdef __cplusplus
}
#endif
#endif
