/*
 * cfx_pi.c - Compute digits of pi
 *
 * Chudnovsky algorithm (1989):
 *   1/π = 12 · Σ ((-1)^k · (6k)! · (13591409 + 545140134k))
 *               / ((3k)! · (k!)^3 · 640320^(3k+3/2))
 *
 * About 14 digits per term. Uses binary splitting for efficiency.
 */

#include "cfx/big.h"
#include "cfx/sbig.h"
#include "cfx_cmd.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void usage(const char* prog) {
    fprintf(stderr, "Usage: %s <digits>\n\n", prog);
    fprintf(stderr, "Compute digits of pi (Chudnovsky algorithm)\n\n");
    fprintf(stderr, "Output options:\n");
    fprintf(stderr, "  -w <cols>         Wrap at <cols> columns (default: 80)\n");
    fprintf(stderr, "  -w 0              No wrapping\n");
    fprintf(stderr, "  -v                Verbose printing\n");
    fprintf(stderr, "  --raw             No formatting, just digits\n\n");
    fprintf(stderr, "Examples:\n");
    fprintf(stderr, "  %s 100                       # first 100 decimal digits\n", prog);
    fprintf(stderr, "  %s 1000                      # ~0.01 sec\n", prog);
    fprintf(stderr, "  %s 10000                     # ~0.5 sec\n", prog);
}


static void big_exp_u64(cfx_big_t* out, const cfx_big_t* n, cfx_limb_t p, int verbose) {
    if (p == 0) {
        cfx_big_from_limb(out, 1);
        return;
    }
    if (cfx_big_is_zero(n)) {
        cfx_big_from_limb(out, 0);
        return;
    }
    if (cfx_big_eq_u64(n, 1)) {
        cfx_big_from_limb(out, 1);
        return;
    }

    cfx_big_t acc, np; /* accumulator, p copy, n^p*/
    cfx_big_init(&acc);
    cfx_big_init(&np);
    cfx_big_from_limb(&np, 1);
    cfx_big_copy(&acc, n);
    cfx_limb_t p_orig = p;
    while (p) {
        if (p & 1) {
            cfx_big_mul_auto(&np, &acc);
        }
        p >>= 1;
        if (p) cfx_big_sq_eq(&acc);
        if (verbose) { printf("pow " CFX_PRIuLIMB " / " CFX_PRIuLIMB "\n", p, p_orig); }
    }
    cfx_big_move(out, &np);
    cfx_big_free(&acc);
}

/*
 * Chudnovsky algorithm using binary splitting
 *
 * 1/π = 12 · Σ ((-1)^k · (6k)! · (13591409 + 545140134k))
 *            / ((3k)! · (k!)^3 · 640320^(3k+3/2))
 *
 * We compute P(0,n), Q(0,n), T(0,n) where:
 *   P(a,b) = product of numerator factors (always positive)
 *   Q(a,b) = product of denominator factors (always positive)
 *   T(a,b) = signed sum (can be negative due to (-1)^k)
 *
 * Then π = Q * 426880 * sqrt(10005) / T
 */

#define C 640320ULL
#define C3_OVER_24 (C * C * C / 24)  /* 10939058860032000 */
#define A 13591409LL
#define B 545140134LL

/*
 * Binary splitting for Chudnovsky
 * P, Q are unsigned (always positive)
 * T is signed (alternates based on (-1)^k)
 */
static void bs(size_t a, size_t b, cfx_big_t* P, cfx_big_t* Q, cfx_sbig_t* T, int verbose) {

    if (verbose) printf("bs(%zu, %zu)\n", a, b);
    if (b - a == 1) {
        /* Base case: compute for single term k = a */
        cfx_big_t tmp1;
        cfx_big_init(&tmp1);

        if (a == 0) {
            cfx_big_from_u64(P, 1);
            cfx_big_from_u64(Q, 1);
        } else {
            size_t k = a;

            /* P_k = (6k-5)(2k-1)(6k-1) */
            cfx_big_from_u64(P, 6*k - 5);
            cfx_big_from_u64(&tmp1, 2*k - 1);
            cfx_big_mul_eq(P, &tmp1);
            cfx_big_from_u64(&tmp1, 6*k - 1);
            cfx_big_mul_eq(P, &tmp1);

            /* Q_k = k^3 * C^3/24 */
            cfx_big_from_u64(Q, k);
            cfx_big_from_u64(&tmp1, k);
            cfx_big_mul_eq(Q, &tmp1);
            cfx_big_mul_eq(Q, &tmp1);  /* k^3 */

            cfx_big_from_u64(&tmp1, (uint64_t)C3_OVER_24);
            cfx_big_mul_eq(Q, &tmp1);
        }

        /* T_k = P_k * (A + B*k) * (-1)^k */
        cfx_big_t t_mag;
        cfx_big_init(&t_mag);

        cfx_big_from_u64(&t_mag, (uint64_t)B);
        cfx_big_from_u64(&tmp1, (uint64_t)a);
        cfx_big_mul_eq(&t_mag, &tmp1);
        cfx_big_from_u64(&tmp1, (uint64_t)A);
        cfx_big_add_eq(&t_mag, &tmp1);
        cfx_big_mul_eq(&t_mag, P);

        /* Set T with sign based on (-1)^a */
        int8_t sign = (a % 2 == 0) ? 1 : -1;
        cfx_sbig_assign_big(T, &t_mag, sign);

        cfx_big_free(&t_mag);
        cfx_big_free(&tmp1);
    } else {
        /* Recursive case: split in the middle */
        size_t m = (a + b) / 2;

        cfx_big_t P1, Q1, P2, Q2;
        cfx_sbig_t T1, T2;

        cfx_big_init(&P1); cfx_big_init(&Q1);
        cfx_big_init(&P2); cfx_big_init(&Q2);
        cfx_sbig_init(&T1); cfx_sbig_init(&T2);

        bs(a, m, &P1, &Q1, &T1, verbose);
        bs(m, b, &P2, &Q2, &T2, verbose);

        /* Combine:
         * P = P1 * P2
         * Q = Q1 * Q2
         * T = T1 * Q2 + T2 * P1  (signed arithmetic)
         */
        cfx_big_mul(P, &P1, &P2);
        cfx_big_mul(Q, &Q1, &Q2);

        /* T = T1 * Q2 + T2 * P1 */
        cfx_sbig_t term1, term2;
        cfx_sbig_init(&term1);
        cfx_sbig_init(&term2);

        /* term1 = T1 * Q2 */
        cfx_big_t prod1;
        cfx_big_init(&prod1);
        cfx_big_mul(&prod1, &T1.mag, &Q2);
        cfx_sbig_assign_big(&term1, &prod1, T1.sign);
        cfx_big_free(&prod1);

        /* term2 = T2 * P1 */
        cfx_big_t prod2;
        cfx_big_init(&prod2);
        cfx_big_mul(&prod2, &T2.mag, &P1);
        cfx_sbig_assign_big(&term2, &prod2, T2.sign);
        cfx_big_free(&prod2);

        /* T = term1 + term2 */
        cfx_sbig_add(T, &term1, &term2);

        cfx_sbig_free(&term1);
        cfx_sbig_free(&term2);
        cfx_big_free(&P1); cfx_big_free(&Q1);
        cfx_big_free(&P2); cfx_big_free(&Q2);
        cfx_sbig_free(&T1); cfx_sbig_free(&T2);
    }
}

/*
 * Integer square root using Newton's method
 */
static void big_isqrt(cfx_big_t* result, const cfx_big_t* n) {
    if (cfx_big_is_zero(n)) {
        cfx_big_from_u64(result, 0);
        return;
    }

    size_t bits = cfx_big_bitlen(n);
    cfx_big_t x, x_new, tmp, rem;
    cfx_big_init(&x);
    cfx_big_init(&x_new);
    cfx_big_init(&tmp);
    cfx_big_init(&rem);

    cfx_big_from_u64(&x, 1);
    cfx_big_shl_bits_eq(&x, (unsigned)((bits + 1) / 2));

    while (1) {
        cfx_big_divrem(&tmp, &rem, n, &x);
        cfx_big_add(&x_new, &x, &tmp);
        cfx_big_shr_bits_eq(&x_new, 1);

        if (cfx_big_cmp(&x_new, &x) >= 0) {
            break;
        }
        cfx_big_assign(&x, &x_new);
    }

    cfx_big_assign(result, &x);

    cfx_big_free(&x);
    cfx_big_free(&x_new);
    cfx_big_free(&tmp);
    cfx_big_free(&rem);
}


static void compute_pi_chudnovsky(cfx_big_t* pi, size_t digits, int verbose) {
    size_t n_terms = digits / 14 + 2;
    size_t precision = digits + 50;

    fprintf(stderr, "Computing %zu terms via binary splitting...\n", n_terms);

    cfx_big_t P, Q;
    cfx_sbig_t T;
    cfx_big_init(&P);
    cfx_big_init(&Q);
    cfx_sbig_init(&T);

    bs(0, n_terms, &P, &Q, &T, verbose);

    fprintf(stderr, "Computing sqrt(10005)...\n");

    /* π = Q * 426880 * sqrt(10005) / T
     * Note: 426880 = 53360 * 8, and 53360 = 640320/12, and sqrt(640320) = 8*sqrt(10005)
     * So this is equivalent to: Q * 640320 * sqrt(640320) / (12 * T)
     */
    cfx_big_t scale, k_scaled, sqrt_k, tmp, rem;
    cfx_big_init(&scale);
    cfx_big_init(&k_scaled);
    cfx_big_init(&sqrt_k);
    cfx_big_init(&tmp);
    cfx_big_init(&rem);

    cfx_big_t ten;
    cfx_big_init(&ten);
    cfx_big_from_u64(&ten, 10);
    big_exp_u64(&scale, &ten, (cfx_limb_t)precision, verbose);

    cfx_big_from_u64(&k_scaled, 10005);
    cfx_big_t scale_sq;
    cfx_big_init(&scale_sq);
    cfx_big_mul(&scale_sq, &scale, &scale);
    cfx_big_mul_eq(&k_scaled, &scale_sq);

    big_isqrt(&sqrt_k, &k_scaled);

    /* numerator = Q * sqrt(10005) * 426880 */
    cfx_big_mul(&tmp, &Q, &sqrt_k);
    cfx_big_mul_sm_eq(&tmp, 426880);

    /* π = numerator / |T| */
    cfx_big_divrem(pi, &rem, &tmp, &T.mag);

    /* Truncate extra digits */
    cfx_big_t extra;
    cfx_big_init(&extra);
    big_exp_u64(&extra, &ten, (cfx_limb_t)(precision - digits), verbose);
    cfx_big_divrem(pi, &rem, pi, &extra);

    cfx_big_free(&extra);
    cfx_big_free(&ten);
    cfx_big_free(&scale);
    cfx_big_free(&scale_sq);
    cfx_big_free(&k_scaled);
    cfx_big_free(&sqrt_k);
    cfx_big_free(&tmp);
    cfx_big_free(&rem);
    cfx_big_free(&P);
    cfx_big_free(&Q);
    cfx_sbig_free(&T);
}

static void print_pi_formatted(const char* s, size_t len, int wrap_cols, int raw) {
    if (raw) {
        printf("%s\n", s);
        return;
    }

    if (len == 0) {
        printf("3\n");
        return;
    }

    printf("3.");

    const char* digits = s + 1;
    size_t remaining = len - 1;

    if (wrap_cols == 0) {
        printf("%s\n", digits);
        return;
    }

    size_t col = 2;
    for (size_t i = 0; i < remaining; i++) {
        putchar(digits[i]);
        col++;
        if (col >= (size_t)wrap_cols && i + 1 < remaining) {
            putchar('\n');
            col = 0;
        }
    }
    putchar('\n');
}

int cfx_pi_run(int argc, char* argv[]) {
    if (argc < 2 ||
        (argc == 2 && (strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0))) {
        usage(argv[0]);
        return argc < 2 ? 1 : 0;
    }

    int wrap_cols = 80;
    int raw = 0;
    int verbose = 0;
    size_t digits = 0;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-w") == 0 && i + 1 < argc) {
            wrap_cols = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-v") == 0) {
            verbose = 1;
        } else if (strcmp(argv[i], "--raw") == 0) {
            raw = 1;
        } else if (argv[i][0] != '-') {
            digits = (size_t)atol(argv[i]);
        } else {
            fprintf(stderr, "Unknown option: %s\n", argv[i]);
            usage(argv[0]);
            return 1;
        }
    }

    if (digits == 0) {
        fprintf(stderr, "Error: specify number of digits\n");
        usage(argv[0]);
        return 1;
    }

    if (digits > 10000000) {
        fprintf(stderr, "Warning: %zu digits will take a while...\n", digits);
    }

    cfx_big_t pi;
    cfx_big_init(&pi);

    fprintf(stderr, "Computing %zu digits of pi...\n", digits);
    compute_pi_chudnovsky(&pi, digits, verbose);

    fprintf(stderr, "Converting to decimal...\n");

    size_t len;
    char* s = cfx_big_dec_alloc(&pi, &len);

    print_pi_formatted(s, len, wrap_cols, raw);

    free(s);
    cfx_big_free(&pi);

    return 0;
}
