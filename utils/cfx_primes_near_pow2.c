/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cfx/big.h"
#include "cfx/compat.h"
#include "cfx_utils.h"

#define MAX_DELTA 1000u
static cfx_mutex_t g_print_mu = CFX_MUTEX_INITIALIZER;

static void big_set_pow2(cfx_big_t *x, size_t k) {
    cfx_big_from_u64(x, 1);
    cfx_big_shl_bits(x, x, (unsigned)k);
}

/* dist = |a - b| (not used in the threaded version but keeping it here) */
static void big_sub_abs(const cfx_big_t *a, const cfx_big_t *b, cfx_big_t *dist_out) {
    if (cfx_big_cmp(a, b) >= 0) {
        cfx_big_copy(dist_out, a);
        cfx_big_sub(dist_out, b);
    } else {
        cfx_big_copy(dist_out, b);
        cfx_big_sub(dist_out, a);
    }
}

static void print_big_str(const char *s) {
    printf("%s", s);
}

/* Search symmetrically around target:
 * check target-1, target+1, target-3, target+3, ... up to max_delta.
 * On success:
 *   - prime_out = nearest prime
 *   - dist_out  = distance as a small positive big-int
 *   - *dir_out  = -1 for "below", +1 for "above"
 * Returns 1 on success, 0 if no prime found within max_delta.
 */
static int nearest_prime_near(const cfx_big_t *target,
                              cfx_big_t *prime_out,
                              cfx_big_t *dist_out,
                              int *dir_out,
                              cfx_limb_t max_delta)
{
    cfx_big_t n, delta;
    cfx_big_init(&delta);
    cfx_big_init(&n);
    cfx_big_from_limb(&delta, 1);

    int found = 0;

    while (cfx_big_cmp_sm(&delta, max_delta) <= 0) {

        /* below */
        if (cfx_big_cmp(&delta, target) < 0) {
            cfx_big_copy(&n, target);
            cfx_big_sub(&n, &delta);   /* n = target - delta */

            if (cfx_big_cmp_sm(&n, 2) > 0) {
                if (cfx_big_is_prime(&n)) {
                    cfx_big_copy(prime_out, &n);
                    cfx_big_copy(dist_out, &delta);
                    if (dir_out) *dir_out = -1;
                    found = 1;
                    break;
                }
            }
        }

        /* above */
        cfx_big_copy(&n, target);
        cfx_big_add(&n, &delta);       /* n = target + delta */

        if (cfx_big_is_prime(&n)) {
            cfx_big_copy(prime_out, &n);
            cfx_big_copy(dist_out, &delta);
            if (dir_out) *dir_out = +1;
            found = 1;
            break;
        }

        cfx_big_add_sm(&delta, 2);
    }

    cfx_big_free(&delta);
    cfx_big_free(&n);
    return found;
}

/* ---- Parallelization structs ------------------------------------- */

struct prime_result {
    int ok;
    int dir;        /* -1 = below, +1 = above */
    char *prime_dec;
    char *dist_dec;
};

struct worker_args {
    long k_start;
    long k_end;
    long global_k_start;
    cfx_limb_t max_delta;
    struct prime_result *results; /* array indexed by (k - global_k_start) */
    int print_immediately;
    int print_hex;
};

static void* worker_func(void *arg) {
    struct worker_args *wa = (struct worker_args*)arg;

    cfx_big_t N, p, dist;
    cfx_big_init(&N);
    cfx_big_init(&p);
    cfx_big_init(&dist);

    for (long k = wa->k_start; k <= wa->k_end; ++k) {
        long idx = k - wa->global_k_start;
        struct prime_result *res = &wa->results[idx];

        big_set_pow2(&N, (size_t)k);

        int dir = 0;
        int ok = nearest_prime_near(&N, &p, &dist, &dir, wa->max_delta);

        if (wa->print_immediately) {
            /* Build strings OUTSIDE the lock */
            char *p_str = NULL;
            char *d_str = cfx_big_to_str(&dist, NULL);;

            if (ok) {
                if (wa->print_hex) {
                    p_str = cfx_big_to_hex(&p, NULL);
                } else {
                    p_str = cfx_big_to_str(&p, NULL);
                }
            }

            cfx_mutex_lock(&g_print_mu);

            printf("2^%ld", k);
            if (!ok) {
                printf(" -\n");
            } else {
                printf("%c%s: %s\n", (dir < 0) ? '-' : '+', d_str, p_str);
            }
            fflush(stdout);

            cfx_mutex_unlock(&g_print_mu);

            free(p_str);
            free(d_str);
        } else {
            res->ok = ok;
            if (!ok) {
                res->dir = 0;
                res->prime_dec = NULL;
                res->dist_dec = NULL;
            } else {
                res->dir = dir;
                char* p_str;
                if (wa->print_hex) {
                    p_str = cfx_big_to_hex(&p, NULL);
                } else {
                    p_str = cfx_big_to_str(&p, NULL);
                }
                res->prime_dec = p_str;
                res->dist_dec = cfx_big_to_str(&dist, NULL);
            }
        }
    }

    cfx_big_free(&N);
    cfx_big_free(&p);
    cfx_big_free(&dist);

    return NULL;
}


void usage(const char* prog) {
    fprintf(stderr,
        "Usage: %s k_start [k_end] [max_delta] [options]\n"
        "\n"
        "For each k in [k_start, k_end], compute 2^k and\n"
        "find the closest prime above or below.\n"
        "\n"
        "Positional arguments:\n"
        "  k_start            starting exponent (>= 1)\n"
        "  k_end              ending exponent (default = k_start)\n"
        "  max_delta          search distance limit (default = %d)\n"
        "\n"
        "Options:\n"
        "  -t, --threads N    number of worker threads (default = 1)\n"
        "  -p                print results in order (default: out-of-order)\n"
        "  -x               print results in hex\n",
        prog, MAX_DELTA);
}


int cfx_primes_near_pow2_run(int argc, char **argv) {
    if (argc < 2 || (argc == 2 && (strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0))) {
        usage(argv[0]);
        return argc < 2 ? 1 : 0;
    }

    char *endp = NULL;
    int argi = 1;

    /* ---------------- positional args */

    /* k_start */
    long k_start = strtol(argv[argi++], &endp, 10);
    if (*endp != '\0') {
        fprintf(stderr, "Invalid k_start: '%s'\n", argv[argi-1]);
        return 1;
    }

    /* k_end (optional) */
    long k_end = k_start;
    if (argi < argc && argv[argi][0] != '-') {
        endp = NULL;
        k_end = strtol(argv[argi], &endp, 10);
        if (*endp != '\0') {
            fprintf(stderr, "Invalid k_end: '%s'\n", argv[argi]);
            return 1;
        }
        argi++;
    }

    if (k_start < 1 || k_end < k_start) {
        fprintf(stderr, "Require: 1 <= k_start <= k_end\n");
        return 1;
    }

    /* max_delta (optional) */
    cfx_limb_t max_delta = MAX_DELTA;
    if (argi < argc && argv[argi][0] != '-') {
        endp = NULL;
        long tmp = strtol(argv[argi], &endp, 10);
        if (*endp != '\0' || tmp <= 0) {
            fprintf(stderr, "Invalid max_delta: '%s'\n", argv[argi]);
            return 1;
        }
        max_delta = (cfx_limb_t)tmp;
        argi++;
    }

    /* ---------------- options  */
    int num_threads = 1;
    int print_in_order = 0;
    int print_hex = 0;

    while (argi < argc) {
        if (strcmp(argv[argi], "-p") == 0) {
            print_in_order = 1;
            argi++;
        } else if (strcmp(argv[argi], "-t") == 0 ||
                 strcmp(argv[argi], "--threads") == 0) {

            if (argi + 1 >= argc) {
                fprintf(stderr, "Option %s requires an argument\n", argv[argi]);
                return 1;
            }

            endp = NULL;
            long tmp = strtol(argv[argi + 1], &endp, 10);
            if (*endp != '\0' || tmp <= 0) {
                fprintf(stderr, "Invalid thread count: '%s'\n", argv[argi + 1]);
                return 1;
            }

            num_threads = (int)tmp;
            argi += 2;
        } else if (strcmp(argv[argi], "-x") == 0) {
            print_hex = 1;
            argi++;
        } else {
            fprintf(stderr, "Unknown option: '%s'\n", argv[argi]);
            usage(argv[0]);
            return 1;
        }
    }

    long total_k = k_end - k_start + 1;

    int cpu_count = cfx_cpu_count();
    const int max_threads = cpu_count < total_k ? cpu_count : (int)total_k;
    num_threads = num_threads < max_threads ? num_threads : max_threads;

    printf("# k_start=%ld  k_end=%ld  max_delta=%llu  threads=%d\n",
           k_start, k_end, (unsigned long long)max_delta, num_threads);

    struct prime_result *results = (struct prime_result*)calloc((size_t)total_k, sizeof(*results));
    
    if (!results) {
        fprintf(stderr, "Out of memory for results\n");
        return 1;
    }

    cfx_thread_t *threads = (cfx_thread_t*)malloc((size_t) num_threads * sizeof(cfx_thread_t));
    struct worker_args *wargs = (struct worker_args*)malloc((size_t) num_threads * sizeof(*wargs));

    if (!threads || !wargs) {
        fprintf(stderr, "Out of memory for thread structures\n");
        free(results);
        free(threads);
        free(wargs);
        return 1;
    }

    /* split [k_start, k_end] into  num_threads chunks */
    long remaining = total_k;
    long current_k = k_start;

    for (int i = 0; i <  num_threads; ++i) {
        long chunk = remaining / ( num_threads - i);
        if (chunk <= 0) chunk = 1;

        long ks = current_k;
        long ke = ks + chunk - 1;

        wargs[i].k_start        = ks;
        wargs[i].k_end          = ke;
        wargs[i].global_k_start = k_start;
        wargs[i].max_delta      = max_delta;
        wargs[i].results        = results;
        wargs[i].print_immediately = !print_in_order;
        wargs[i].print_hex = print_hex;

        current_k = ke + 1;
        remaining -= chunk;

        cfx_thread_create(&threads[i], worker_func, &wargs[i]);
    }

    for (int i = 0; i < num_threads; ++i) {
        cfx_thread_join(threads[i], NULL);
    }

    /* print results in order */
    if (print_in_order) {
        for (long k = k_start; k <= k_end; ++k) {
            long idx = k - k_start;
            struct prime_result *res = &results[idx];

            printf("2^%ld", k);

            if (!res->ok) {
                printf(" -\n");
                continue;
            }

            printf("%c", (res->dir < 0) ? '-' : '+');
            print_big_str(res->dist_dec);
            printf(": ");
            print_big_str(res->prime_dec);
            printf("\n");
        }
    }

    /* Cleanup */
    for (long i = 0; i < total_k; ++i) {
        free(results[i].prime_dec);
        free(results[i].dist_dec);
    }
    free(results);
    free(threads);
    free(wargs);

    return 0;
}
