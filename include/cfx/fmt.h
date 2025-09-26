#ifndef CFX_FMT_H
#define CFX_FMT_H

#include <inttypes.h>

/* 8-bit */
#define D8F "%" PRId8
#define U8F "%" PRIu8
#define X8F "%" PRIx8
#define X8F_UP "%" PRIX8

/* 16-bit */
#define D16F "%" PRId16
#define U16F "%" PRIu16
#define X16F "%" PRIx16
#define X16F_UP "%" PRIX16

/* 32-bit */
#define D32F "%" PRId32
#define U32F "%" PRIu32
#define X32F "%" PRIx32
#define X32F_UP "%" PRIX32

/* 64-bit */
#define D64F "%" PRId64
#define CFX_PRIuLIMB "%" PRIu64
#define X64F "%" PRIx64
#define X64F_UP "%" PRIX64

/* intptr/uintptr_t */
#define DPF "%" PRIdPTR
#define UPF "%" PRIuPTR
#define XPF "%" PRIxPTR
#define XPF_UP "%" PRIXPTR

#define ZF "%zu"   /* size_t */
#define TF "%td"   /* ptrdiff_t */

#endif  /* CFX_FMT_H */
