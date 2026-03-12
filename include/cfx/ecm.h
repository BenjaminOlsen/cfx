/*
 * ecm.h - Elliptic Curve Method (ECM) for integer factorization
 *
 * ECM is efficient for finding "small" factors (up to ~60-70 bits) in large
 * composites. Its complexity depends on the size of the smallest factor,
 * not the size of n.
 *
 * Algorithm overview:
 * 1. Pick a random elliptic curve E and point P on it
 * 2. Compute Q = k*P where k = product of prime powers up to B1
 * 3. If gcd(Z_Q, n) is non-trivial, we found a factor!
 * 4. If not, try a different curve
 *
 * We use Montgomery curves By² = x³ + Ax² + x with projective coordinates
 * (X:Z) where x = X/Z. This allows efficient arithmetic using only the
 * x-coordinate (no y needed).
 */

#ifndef CFX_ECM_H
#define CFX_ECM_H

#include "cfx/big.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Point on Montgomery curve in projective coordinates.
 * Represents the x coordinate as X/Z.
 * The y-coordinate is not stored (Montgomery ladder doesn't need it).
 */
typedef struct {
    cfx_big_t X;
    cfx_big_t Z;
} cfx_ecm_point_t;

/* must call before use - allocates */
void cfx_ecm_point_init(cfx_ecm_point_t *P);

/* frees allocated resources */
void cfx_ecm_point_free(cfx_ecm_point_t *P);

/* Copy point: dst = src */
void cfx_ecm_point_copy(cfx_ecm_point_t *dst, const cfx_ecm_point_t *src);

/*
 * Try to find a factor of n using ECM.
 *
 * Parameters:
 *   factor  - Output: a non-trivial factor if found
 *   n       - The number to factor (should be odd and composite)
 *   B1      - Stage 1 bound. Larger = finds larger factors but slower.
 *             Rough guide: B1=50000 for ~40-bit, B1=1000000 for ~50-bit factors
 *   curves  - Number of curves to try. More curves = higher success probability.
 *
 * Returns:
 *   1 if a non-trivial factor was found (stored in 'factor')
 *   0 if no factor found after trying all curves
 */
int cfx_ecm_factor(cfx_big_t *factor, const cfx_big_t *n,
    uint64_t B1, unsigned curves);

/*
 * Auto-tuned ECM: picks B1 and curve count based on n's size.
 * Good default for most uses.
 *
 * Returns:
 *   1 if a non-trivial factor was found
 *   0 if no factor found
 */
int cfx_ecm_factor_auto(cfx_big_t *factor, const cfx_big_t *n);

#ifdef __cplusplus
}
#endif

#endif /* CFX_ECM_H */
