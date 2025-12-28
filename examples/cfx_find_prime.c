#include "cfx/big.h"
#include "cfx/rand.h"
#include "cfx/primes.h"
#include "cfx/macros.h"
#include "cfx/compat.h"

#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>

#define SMALL_PRIME_CNT 1024u
#define MAX_THREAD_CNT 16

static void usage(const char* prog) {
    fprintf(stderr, "Usage: %s <bits> [--safe] [--top2] [--seed=<int>] [--verbose|v] [--threads=<int>]\n", prog);
    fprintf(stderr, "  Generate an N-bit prime.\n");
    fprintf(stderr, "  --safe: safe (Sophie Germain) prime\n");
    fprintf(stderr, "  --top2: make sure top 2 bits are set\n");
    fprintf(stderr, "  --seed: provide rng seed (default random seed by cfx_rand_os)\n");
    fprintf(stderr, "  -v or --verbose: verbose printing\n");
    fprintf(stderr, "  --threads: thread cnt (default 1)\n");
}

enum {
    CFX_PRIME_FLAG_NONE      = 0,
    CFX_PRIME_FLAG_SAFE      = 1 << 0,  /* generate a safe prime p with (p-1)/2 also prime */
    CFX_PRIME_FLAG_TOP2      = 1 << 1,  /* set the top two bits (FIPS-ish) */
};

/* Generate an N-bit prime into `p`. */
/* - `bits`  : number of bits (e.g., 1024, 2048, 4096). */
/* - `rng`   : CSPRNG callback. */
/* - `rngctx`: user context for RNG. */
/* - `flags` : OR of CFX_PRIME_FLAG_*. */
/* Returns 0 on success, nonzero on error. */
int cfx_big_gen_nbit_prime(cfx_big_t* p, size_t bits,
                           cfx_rng_fn rng, void* rngctx,
                           unsigned flags, volatile int *stop_flag,
                           const int threadid, int verbose);

static const uint32_t* SMALL_PRIMES = cfx_primes;

/* --- make a random n bit odd number ---------------------------------------------- */
static int rand_nbit_odd(cfx_big_t* out, size_t bits, cfx_rng_fn rng, void* ctx, int flags) {
    if (bits < 2) return -1;
    size_t nbytes = (bits + 7) / 8;

    uint8_t buf[8192];                  /* buffer for up to 8192*8 = 65536 bits */
    if (nbytes > sizeof(buf)) return -2;

    if (rng(ctx, buf, nbytes) != 0) return -3;

    /* mask the top partial byte so value < 2^bits */
    unsigned excess = (unsigned)(8*nbytes - bits);
    if (excess) buf[0] &= (uint8_t)(0xFFu >> excess);

    /* force top bit to ensure exactly 'bits' bits */
    buf[0] |= (uint8_t)(1u << (7 - excess));

    /* (TODO? Rethink this) force the second highest bit */
    if ((flags & CFX_PRIME_FLAG_TOP2) && bits >= 2) {
        unsigned pos = (7 - excess) - 1;
        if ((int)pos >= 0) buf[0] |= (uint8_t)(1u << pos);
    }

    buf[nbytes - 1] |= 1u;  /* force odd */

    /* every safe prime is congruent to 3 mod 4, proof:
    if p > 3 is prime, then p = 1 or 3 mod 4, 
    then 2p + 1 is either:
    - 2*1 + 1 = 0 mod 4 (not prime)
    - 2*3 + 1 = 3 mod 4 (safe)
    then q = 2p+1 = 3 mod 4 
    */
    if (flags & CFX_PRIME_FLAG_SAFE) {
        uint8_t r = buf[nbytes - 1] & 3u;
        if (r != 3u) {
            buf[nbytes - 1] += 3u - r;
        }
    }

    int rc = cfx_big_from_bytes_be(out, buf, nbytes);
    if (rc != 0) return rc;

    /* sanity: todo */
    /* assert(cfx_big_bitlen(out) == bits); */
    return 0;
}

/* --- Update remainders after adding 2 or 4 for safe primes -------------------- */
static inline void remainders_add_step(uint32_t* rem, uint32_t step, const uint32_t* primes, const size_t prime_cnt) {
    for (size_t i = 0; i < prime_cnt; ++i) {
        uint32_t p = primes[i];
        uint32_t r = rem[i] + step;
        rem[i] = (uint32_t)(r >= p ? (r - p) : r);
    }
}

/* --- Initialize remainder table rem[i] = n % primes[i] ------------------------ */
static void remainders_init(uint32_t* rem, const cfx_big_t* n, const uint32_t* primes, const size_t prime_cnt) {
    for (size_t i = 0; i < prime_cnt; ++i) rem[i] = (uint32_t)cfx_big_mod_sm(n, primes[i]);
}

/* --- Try small-prime rejection; return 0 if divisible by any ----------------- */
static int passes_small_trial(const uint32_t* rem, const uint32_t* primes, const size_t prime_cnt) {
    (void)primes;
    for (size_t i = 0; i < prime_cnt; ++i) {
        if (rem[i] == 0) return 0;
    }
    return 1;
}

/* A safe prime is a prime p that can be written p = 2q + 1 where q is also prime.
AKA Sophie Germain prime. It's called 'safe' because the multiplicative group Zp
(integers mod p without 0) behaves in a way that avoids very small subgroups */
static int is_safe_prime(const cfx_big_t* p) {
    cfx_big_t q;
    cfx_big_init(&q);
    cfx_big_copy(&q, p);
    cfx_big_sub_sm(&q, 1);
    cfx_big_shr_bits_eq(&q, 1);  /* q = (p-1)/2 */
    int is_safe = cfx_big_is_prime(&q);
    cfx_big_free(&q);
    return is_safe;
}


static void print_big_hex(const cfx_big_t* b) {
    if (b->n == 0) { puts("0x0"); return; }
    size_t i = b->n - 1;
    printf("0x%llx", (unsigned long long)b->limb[i]);

    /* remaining limbs zero-padded to full limb width */
    int w = (int)(sizeof(cfx_limb_t) * 2); /* hex digits per limb */
    while (i-- > 0) {
        printf("%0*llx", w, (unsigned long long)b->limb[i]);
    }
    putchar('\n');
}

#define ATOMIC_LOAD(p) cfx_atomic_load(p)
#define ATOMIC_STORE(p, v) cfx_atomic_store(p, v)

/* prime search summary: 
1) Choose a random 'bits' bit odd number
2) fill a list of remainders after division by each SMALL_PRIMES
3) loop:
*/
int cfx_big_gen_nbit_prime(cfx_big_t* p, size_t bits,
                           cfx_rng_fn rng, void* rngctx,
                           unsigned flags, volatile int *stop_flag,
                           const int threadid, int verbose) {
    if (bits < 2 || !p || !rng) return -1;

    /* Draw initial random N-bit odd candidate. */
    cfx_big_t n;
    cfx_big_init(&n);
    int rc = rand_nbit_odd(&n, bits, rng, rngctx, flags);
    if (rc != 0) goto end;

    /* prepare small prime remainders. */
    uint32_t rem[SMALL_PRIME_CNT] = {0};
    remainders_init(rem, &n, SMALL_PRIMES, SMALL_PRIME_CNT);

    uint32_t step = (flags & CFX_PRIME_FLAG_SAFE) ? 4 : 2;

    /* Ensure we’re not accidentally equal to a small prime (only possible if bits are tiny). */
    for (size_t i = 0; i < SMALL_PRIME_CNT; ++i) {
        if (cfx_big_cmp_sm(&n, SMALL_PRIMES[i]) == 0) {
            /* bump by 2 to avoid returning small primes */
            cfx_big_add_sm(&n, 2);
            remainders_add_step(rem, step, SMALL_PRIMES, SMALL_PRIME_CNT);
            break;
        }
    }

    unsigned reject_cnt = 0;
    for (;;) {
         if (stop_flag && ATOMIC_LOAD(stop_flag)) {
            rc = -2; /* cancelled */
            goto end;
        }
        /* Small trial division */
        if (passes_small_trial(rem, SMALL_PRIMES, SMALL_PRIME_CNT)) {
            /* Probable-prime test (Miller–Rabin in cfx_big_is_prime) */
            if (cfx_big_is_prime(&n)) {
                if (flags & CFX_PRIME_FLAG_SAFE) {
                    if (is_safe_prime(&n)) {
                        cfx_big_copy(p, &n);
                        goto end;
                    } else if (verbose) {
                        printf("[thread %d] - REJECTED prime #%u for not being safe:\n", threadid, ++reject_cnt);
                        print_big_hex(&n);
                    }
                } else {
                    cfx_big_copy(p, &n);
                    goto end;
                }
            }
        }

        /* Next odd candidate: n += step */
        cfx_big_add_sm(&n, step);
        remainders_add_step(rem, step, SMALL_PRIMES, SMALL_PRIME_CNT);

        /* Keep within N bits: if we overflowed (rare), redraw. */
        if (cfx_big_bitlen(&n) > bits) {
            rc = rand_nbit_odd(&n, bits, rng, rngctx, (flags & CFX_PRIME_FLAG_TOP2) != 0);
            if (rc != 0) goto end;
            remainders_init(rem, &n, SMALL_PRIMES, SMALL_PRIME_CNT);
        }
    }

end:
    cfx_big_free(&n);
    return rc;
}


typedef struct {
    size_t bits;
    unsigned flags;
    uint32_t base_seed;
    int thread_id;
    cfx_big_t *result;
    cfx_atomic_int *stop_flag;       /* 0 = keep going, 1 = someone won */
    int verbose;
    cfx_mutex_t *result_mutex;
} prime_worker_args_t;


static void* worker(void* arg) {
    prime_worker_args_t* w = (prime_worker_args_t*)arg;

    /* Quick early exit if someone already won before we even start */
    if (ATOMIC_LOAD(w->stop_flag)) {
        return NULL;
    }

    cfx_rng_ctx_t rng_ctx;
    uint32_t u32 = (uint32_t)w->thread_id;
    uint32_t seed = w->base_seed ^ cfx_splitmix32(&u32);
    cfx_chacha20_rng_init(&rng_ctx, seed);

    cfx_big_t local;
    cfx_big_init(&local);

    int rc = cfx_big_gen_nbit_prime(&local, w->bits,
                                       cfx_chacha20_rng, &rng_ctx,
                                       w->flags,
                                       w->stop_flag,
                                       w->thread_id,
                                       w->verbose);
    if (rc != 0) {
        /* rc > 0 : cancelled, rc < 0 : error; either way, we don't claim victory */
        if (w->verbose) {
            printf("[thread %d] cancelled\n", w->thread_id);
        }
        cfx_big_free(&local);
        return NULL;
    }

    /* We found a prime. Try to publish it if we're first. */
    cfx_mutex_lock(w->result_mutex);
    if (!ATOMIC_LOAD(w->stop_flag)) {
        /* We're the first winner */
        cfx_big_copy(w->result, &local);
        ATOMIC_STORE(w->stop_flag, 1);
        if (w->verbose) {
            printf("[thread %d] - I WIN!\n", w->thread_id);
            print_big_hex(&local);
        }
    } else if (w->verbose) {
        printf("[thread %d] found a prime but not the first...\n", w->thread_id);
        print_big_hex(&local);
    }
    cfx_mutex_unlock(w->result_mutex);

    cfx_big_free(&local);
    return NULL;
}

int main(int argc, char* argv[]) {
    const char* prog = argv[0];
    if (argc < 2) {
        usage(prog);
        return 1;
    }

    if (SMALL_PRIME_CNT > cfx_primes_len) {
        printf("have to compile with CFX_PRIMES_COUNT >= %u\n", SMALL_PRIME_CNT);
        return 1;
    }

    unsigned flags = CFX_PRIME_FLAG_NONE;
    size_t bits = 0;
    int verbose = 0;
    unsigned threadcnt = 1;
    cfx_srand_os();
    unsigned seed = cfx_rand_os();

    for (int argi = 1; argi < argc; ++argi) {
        const char* arg = argv[argi];
        if (strcmp(argv[argi], "--top2") == 0) {
            flags |= CFX_PRIME_FLAG_TOP2;
        } else if (strcmp(argv[argi], "--safe") == 0) {
            flags |= CFX_PRIME_FLAG_SAFE;
        } else if ((strcmp(argv[argi], "--help") == 0) || (strcmp(argv[argi], "-h") == 0)) {
            usage(argv[0]);
            return 1;
        } else if (strncmp(arg, "--seed=", 7) == 0) {
            char* end = NULL;
            seed = strtoull(arg + 7, &end, 0);  /* base 0: deduces base from a */
            if (end == arg + 7) {
                fprintf(stderr, "Invalid seed: %s\n\n", arg + 7);
                usage(prog);
                return EXIT_FAILURE;
            }
        } else if (strncmp(arg, "--threads=", 10) == 0) {
            char* end = NULL;
            threadcnt = strtoul(arg + 10, &end, 0);
            if (end == arg + 10) {
                fprintf(stderr, "Invalid threadcnt: %s\n\n", arg + 10);
                usage(prog);
                return EXIT_FAILURE;
            }
        } else if ((strcmp(argv[argi], "--verbose") == 0) || (strcmp(argv[argi], "-v") == 0)) {
            verbose = 1;
        } else {
            char* endp = NULL;
            bits = (size_t)strtol(argv[argi], &endp, 10);
            if (*argv[argi] == '\0' || (endp && *endp != '\0')){
                return 1;
            }
        }

        if (argi > argc - 1) {
            usage(argv[0]);
            return 1;
        }
    }

    if (bits < 2) {
        fprintf(stderr, "error: <bits> must be an integer >= 2\n");
        return 1;
    }

    if (threadcnt < 1 || threadcnt > MAX_THREAD_CNT) {
        fprintf(stderr, "threadcnt should be > 1 and < %u\n", MAX_THREAD_CNT);
        return 1;
    }

    cfx_big_t result;
    cfx_big_init(&result);
    cfx_atomic_int stop_flag = 0;
    cfx_mutex_t result_mutex = CFX_MUTEX_INITIALIZER;
    cfx_thread_t threads[MAX_THREAD_CNT];
    prime_worker_args_t args[MAX_THREAD_CNT];

    for (int i = 0; i < threadcnt; ++i) {
        args[i].bits            = bits;
        args[i].flags           = flags;
        args[i].base_seed       = seed;
        args[i].thread_id       = i;
        args[i].result          = &result;
        args[i].stop_flag       = &stop_flag;
        args[i].result_mutex    = &result_mutex;
        args[i].verbose         = verbose;
        

        int rc = cfx_thread_create(&threads[i], worker, &args[i]);
        if (rc != 0) {
            fprintf(stderr, "thread_create failed: %d\n", rc);
            ATOMIC_STORE(&stop_flag, 1);  /* don't start more if the first failed.... */
            /* you might want to join already-started threads and abort */
        }
    }

    for (int i = 0; i < threadcnt; ++i) {
        cfx_thread_join(threads[i], NULL);
    }

    if (!ATOMIC_LOAD(&stop_flag)) {
        fprintf(stderr, "No prime found (all threads aborted?)\n");
        cfx_big_free(&result);
        return 1;
    }

    puts("---------------");
    print_big_hex(&result);
    printf(" bits = %zu\n", cfx_big_bitlen(&result));

    cfx_big_free(&result);
    return 0;
}
