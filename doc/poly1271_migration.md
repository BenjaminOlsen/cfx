# Poly1271 Migration from Reference Implementation

This document logs the migration of poly1271 code from `../poly1271` reference implementation to cfx.

## Date: 2026-01-16

## Overview

Poly1271 is a polynomial MAC using the Mersenne prime 2^127 - 1. It's similar to Poly1305 but uses:
- 15-byte blocks (vs 16 for Poly1305)
- Mersenne prime reduction (just addition, no Montgomery reduction)
- ~123 bits security (vs ~103 for Poly1305)

## Source Files

### Reference Implementation
- `../poly1271/src/poly1271.c` - Scalar implementation (~709 lines)
- `../poly1271/src/poly1271_avx2.c` - AVX2 intrinsics implementation (~753 lines)
- `../poly1271/src/poly1271_avx2_asm_v2_macos.c` - AVX2 ASM for macOS (~545 lines)
- `../poly1271/src/poly1271_avx2_asm.c` - AVX2 ASM for Linux
- `../poly1271/test/test_poly1271.c` - Test suite (~795 lines)

### CFX Files Created/Modified
- `src/crypto/poly1271/poly1271.c` - Scalar implementation (adapted)
- `src/crypto/poly1271/poly1271_avx2.c` - AVX2 intrinsics implementation (adapted)
- `include/cfx/poly1271.h` - Header file (updated with verify functions)
- `test/unit/test_poly1271.c` - Tests (enhanced with AVX2 tests)
- `src/crypto/CMakeLists.txt` - Build config (updated)

## Changes Made

### Naming Convention
All functions renamed from `poly1271_*` to `cfx_poly1271_*`:
- `poly1271_init` -> `cfx_poly1271_init`
- `poly1271_update` -> `cfx_poly1271_update`
- `poly1271_finish` -> `cfx_poly1271_finish`
- `poly1271` -> `cfx_poly1271`
- `poly1271_verify` -> `cfx_poly1271_verify`
- `poly1271_avx2_*` -> `cfx_poly1271_avx2_*`

### Types
- `poly1271_ctx_t` -> `cfx_poly1271_ctx_t`
- `poly1271_avx2_ctx_t` -> `cfx_poly1271_avx2_ctx_t`

### Memory Zeroing
- `secure_zero()` -> `CFX_MEMZERO_S()`

### Constants
- `POLY1271_BLOCK_SIZE` -> `CFX_POLY1271_BLOCK_SIZE`
- `POLY1271_KEY_SIZE` -> `CFX_POLY1271_KEY_SIZE`
- `POLY1271_TAG_SIZE` -> `CFX_POLY1271_TAG_SIZE`

### AVX2 Detection
- `#ifdef __AVX2__` -> `#if CFX_HAVE_AVX2`

## Implementation Details

### Scalar Implementation (`poly1271.c`)
- Uses radix-2^64 representation (2 limbs)
- 64-bit multiplication with `__int128` or MSVC intrinsics
- Precomputes r, r^2, r^3, r^4 for multi-block processing
- Lazy reduction every 8 blocks
- Constant-time finalization

### AVX2 Implementation (`poly1271_avx2.c`)
- Uses radix-2^26 representation (5 limbs)
- Processes 4 blocks in parallel
- Uses `_mm256_i64gather_epi64` for efficient block loading
- Schoolbook multiplication with folding (2^130 = 8 mod p)
- Horizontal sum for combining parallel lanes

## Test Vectors

All test vectors from reference implementation preserved:
- Empty message
- Single byte (0x0d)
- Exactly 1 block (15 bytes)
- 1 block + 1 byte (16 bytes)
- Exactly 2 blocks (30 bytes)
- 100 bytes
- 256 bytes
- 1000 bytes
- Walt Whitman quote (120 bytes)

## Status

**Migration Complete** - 2026-01-16

All 41 cfx tests pass including:
- `test_poly1271` with scalar and AVX2 tests
- `test_aead_chacha20_poly1271` (uses poly1271 for authentication)

## TODO

- [ ] Consider adding ASM implementations for performance
- [ ] Benchmark AVX2 vs scalar
- [ ] Add runtime dispatch (auto-select best implementation)

## Files Not Migrated

The following files from the reference were not migrated:
- `poly1271_opt.c` - Optimized scalar (functionality merged into poly1271.c)
- `poly1271_dispatch.c` - Runtime dispatch (not needed for compile-time selection)
- `poly1271_avx2_asm*.c` - ASM implementations (intrinsics sufficient for now)
