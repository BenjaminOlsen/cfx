/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#ifndef CFX_TARGET_H
#define CFX_TARGET_H

/*
 * Target-specific build information.
 *
 * CFX_TARGET_NAME is a string literal set at compile time by CMake,
 * identifying which optimized backend was selected.
 *
 * Example values: "base", "x86_64_bmi2", "x86_64_avx2", "arm_cortex_m4"
 *
 * Usage:
 *   printf("cfx built for: %s\n", CFX_TARGET_NAME);
 */

#ifndef CFX_TARGET_NAME
#define CFX_TARGET_NAME "base"
#endif

/*
 * Capability macros (defined by CMake when the target supports them):
 *
 *   CFX_CAP_BMI2      - x86 BMI2 intrinsics (_mulx_u64, _addcarry_u64)
 *   CFX_CAP_AVX2      - x86 AVX2 SIMD
 *   CFX_CAP_AVX512    - x86 AVX-512 SIMD
 *   CFX_CAP_NEON      - ARM NEON SIMD
 *   CFX_CAP_DSP       - ARM DSP instructions (Cortex-M4 UMULL/UMLAL)
 *   CFX_CAP_RVV       - RISC-V Vector extension
 *   CFX_CAP_WASM_SIMD - WebAssembly SIMD128
 */

/* Target hierarchy (compile-time only, for reference):
 *
 * base                         <- root: implements EVERYTHING
 * ├── x86_64                      <- inherits from base
 * │   └── x86_64_bmi2             <- inherits from x86_64
 * │       └── x86_64_avx2         <- inherits from x86_64_bmi2
 * │           └── x86_64_avx512   <- inherits from x86_64_avx2
 * ├── arm_cortex_m4               <- inherits from base (no NEON)
 * ├── arm_neon                    <- inherits from base
 * │   └── aarch64_neon            <- inherits from arm_neon
 * ├── riscv64                     <- inherits from base
 * │   └── riscv64_v               <- inherits from riscv64
 * └── wasm32                      <- inherits from base
 *     └── wasm_simd128            <- inherits from wasm32
 */

#ifdef __cplusplus
extern "C" {
#endif

/* Returns the target name string (same as CFX_TARGET_NAME macro) */
static inline const char* cfx_target_name(void) {
    return CFX_TARGET_NAME;
}

#ifdef __cplusplus
}
#endif

#endif /* CFX_TARGET_H */
