/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#ifndef CFX_TARGET_H
#define CFX_TARGET_H

/*
 * Target-specific build information.
 *
 * CFX_TARGET_NAME is a string literal set at compile time by CMake,
 * identifying which optimized backend was selected.
 *
 * Example values: "portable", "x86_64_bmi2", "x86_64_avx2", "arm_cortex_m4"
 *
 * Usage:
 *   printf("cfx built for: %s\n", CFX_TARGET_NAME);
 */

#ifndef CFX_TARGET_NAME
#define CFX_TARGET_NAME "portable"
#endif

/*
 * Capability macros (defined by CMake when the target supports them):
 *
 *   CFX_CAP_AVX2      - x86 AVX2 SIMD
 */

/* Target hierarchy (compile-time only, for reference):
 *
 * portable                 <- root: implements EVERYTHING
 * ├── x86_64_bmi2          <- inherits from portable
 * │   └── x86_64_avx2      <- inherits from x86_64_bmi2
 * ├── arm_cortex_m4        <- inherits from portable
 */

#ifdef __cplusplus
extern "C" {
#endif

/* Returns the target name string (same as CFX_TARGET_NAME macro) */
static inline const char * cfx_target_name(void) {
    return CFX_TARGET_NAME;
}

#ifdef __cplusplus
}
#endif

#endif /* CFX_TARGET_H */
