#ifndef CFX_TYPES_H
#define CFX_TYPES_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define CFX_BIG_BASE 1000000000u
#define CFX_DIGITS_PER_LIMB 9

typedef struct {
    uint64_t hi, lo;
} cfx_u128;

#ifdef __cplusplus
}
#endif
#endif
