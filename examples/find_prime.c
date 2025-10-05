#include "cfx/big.h"

#include <fcntl.h>
#include <unistd.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>

// RNG: fill `out[len]` with cryptographic random bytes. Return 0 on success.
typedef int (*cfx_rng_fn)(void* ctx, uint8_t* out, size_t len);

// Flags
enum {
    CFX_PRIME_FLAG_NONE      = 0,
    CFX_PRIME_FLAG_SAFE      = 1 << 0,  // generate a safe prime p with (p-1)/2 also prime
    CFX_PRIME_FLAG_TOP2      = 1 << 1,  // set the top two bits (FIPS-ish)
};

// Generate an N-bit prime into `p`.
// - `bits`  : number of bits (e.g., 1024, 2048, 4096).
// - `rng`   : CSPRNG callback.
// - `rngctx`: user context for RNG.
// - `flags` : OR of CFX_PRIME_FLAG_*.
// Returns 0 on success, nonzero on error.
int cfx_big_gen_nbit_prime(cfx_big_t* p, size_t bits,
                           cfx_rng_fn rng, void* rngctx,
                           unsigned flags);

static const uint16_t SMALL_PRIMES[] = {
    3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,
    101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,
    191,193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,
    281,283,293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,
    389,397,401,409,419,421,431,433,439,443,449,457,461,463,467,479,487,
    491,499,503,509,521,523,541,0
};

// --- make a random n bit odd number ----------------------------------------------
static int rand_nbit_odd(cfx_big_t* out, size_t bits, cfx_rng_fn rng, void* ctx, int set_top2) {
    if (bits < 2) return -1;
    size_t nbytes = (bits + 7) / 8;
    uint8_t buf[8192];                  // stack buffer for up to 65536 bits; adjust if needed
    if (nbytes > sizeof(buf)) return -2;

    if (rng(ctx, buf, nbytes) != 0) return -3;

    // mask the top partial byte so value < 2^bits
    unsigned excess = (unsigned)(8*nbytes - bits);
    if (excess) buf[0] &= (uint8_t)(0xFFu >> excess);

    // force top bit to ensure exactly 'bits' bits
    buf[0] |= (uint8_t)(1u << (7 - excess));

    // (TODO? Rethink this) force the second highest bit
    if (set_top2 && bits >= 2) {
        unsigned pos = (7 - excess) - 1;
        if ((int)pos >= 0) buf[0] |= (uint8_t)(1u << pos);
    }

    buf[nbytes - 1] |= 1u;  // force odd

    int rc = cfx_big_from_bytes_be(out, buf, nbytes);
    if (rc != 0) return rc;

    // sanity: todo
    assert(cfx_big_bitlen(out) == bits);
    return 0;
}

// --- Update remainders after adding +2 ---------------------------------------
static inline void remainders_add2(uint16_t* rem, const uint16_t* primes) {
    for (size_t i = 0; primes[i]; ++i) {
        uint16_t p = primes[i];
        uint16_t r = rem[i] + 2;
        rem[i] = (uint16_t)(r >= p ? (r - p) : r);
    }
}

// --- Initialize remainder table rem[i] = n % primes[i] ------------------------
static void remainders_init(uint16_t* rem, const cfx_big_t* n, const uint16_t* primes) {
    for (size_t i = 0; primes[i]; ++i) rem[i] = (uint16_t)cfx_big_mod_sm(n, primes[i]);
}

// --- Try small-prime rejection; return 0 if divisible by any -----------------
static int passes_small_trial(const uint16_t* rem, const uint16_t* primes) {
    for (size_t i = 0; primes[i]; ++i) {
        if (rem[i] == 0) return 0;
    }
    return 1;
}

static int is_safe_prime(const cfx_big_t* p) {
    cfx_big_t q, t;
    cfx_big_copy(&t, p);
    cfx_big_sub_sm(&t, 1);
    cfx_big_copy(&q, &t);
    cfx_big_shr_bits_eq(&q, 1);  // q = (p-1)/2
    return cfx_big_is_prime(&q);
}

// --- Main generator -----------------------------------------------------------
int cfx_big_gen_nbit_prime(cfx_big_t* p, size_t bits,
                           cfx_rng_fn rng, void* rngctx,
                           unsigned flags) {
    if (bits < 2 || !p || !rng) return -1;

    // 1) Draw initial random N-bit odd candidate.
    cfx_big_t n;
    cfx_big_init(&n);
    int rc = rand_nbit_odd(&n, bits, rng, rngctx, (flags & CFX_PRIME_FLAG_TOP2) != 0);
    if (rc != 0) return rc;

    // 2) Prepare small-prime remainders.
    uint16_t rem[sizeof(SMALL_PRIMES)/sizeof(SMALL_PRIMES[0])] = {0};
    remainders_init(rem, &n, SMALL_PRIMES);

    // Ensure we’re not accidentally equal to a small prime (only possible if bits are tiny).
    for (size_t i = 0; SMALL_PRIMES[i]; ++i) {
        if (cfx_big_cmp_sm(&n, SMALL_PRIMES[i]) == 0) {
            // bump by 2 to avoid returning small primes
            cfx_big_add_sm(&n, 2);
            remainders_add2(rem, SMALL_PRIMES);
            break;
        }
    }

    // 3) Search loop
    for (;;) {
        // Small trial division
        if (passes_small_trial(rem, SMALL_PRIMES)) {
            // Probable-prime test (Miller–Rabin in your cfx_big_is_prime)
            if (cfx_big_is_prime(&n)) {
                if (flags & CFX_PRIME_FLAG_SAFE) {
                    if (is_safe_prime(&n)) {
                        cfx_big_copy(p, &n);
                        return 0;
                    }
                } else {
                    cfx_big_copy(p, &n);
                    return 0;
                }
            }
        }

        // Next odd candidate: n += 2
        cfx_big_add_sm(&n, 2);
        remainders_add2(rem, SMALL_PRIMES);

        // Keep within N bits: if we overflowed (rare), redraw.
        if (cfx_big_bitlen(&n) > bits) {
            rc = rand_nbit_odd(&n, bits, rng, rngctx, (flags & CFX_PRIME_FLAG_TOP2) != 0);
            if (rc != 0) return rc;
            remainders_init(rem, &n, SMALL_PRIMES);
        }
    }
}


static int urandom_rng(void* ctx, uint8_t* out, size_t len) {
    (void)ctx;
    int fd = open("/dev/urandom", O_RDONLY);
    if (fd < 0) return -1;
    ssize_t got = read(fd, out, len);
    close(fd);
    return (got == (ssize_t)len) ? 0 : -1;
}



// optional: update usage text to match this tool
static void usage(const char* prog) {
    fprintf(stderr, "Usage: %s <bits> [--safe] [--top2]\n", prog);
    fprintf(stderr, "  Generate an N-bit prime (optionally safe) and print it in hex.\n");
}

// tiny hex printer (assumes cfx_big_t uses b->limb[] little-endian limbs)
static void print_big_hex(const cfx_big_t* b) {
    if (b->n == 0) { puts("0x0"); return; }

    // print most-significant limb without zero padding
    size_t i = b->n - 1;
    printf("0x%llx", (unsigned long long)b->limb[i]);

    // remaining limbs zero-padded to full limb width
    int w = (int)(sizeof(cfx_limb_t) * 2); // hex digits per limb
    while (i-- > 0) {
        printf("%0*llx", w, (unsigned long long)b->limb[i]);
    }
    putchar('\n');
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        usage(argv[0]);
        return 1;
    }

    char* endp = NULL;
    long nbits_long = strtol(argv[1], &endp, 10);
    if (*argv[1] == '\0' || (endp && *endp != '\0') || nbits_long < 2) {
        fprintf(stderr, "error: <bits> must be an integer >= 2\n");
        return 1;
    }
    size_t bits = (size_t)nbits_long;

    unsigned flags = CFX_PRIME_FLAG_NONE;
    for (int i = 2; i < argc; ++i) {
        if (strcmp(argv[i], "--safe") == 0) flags |= CFX_PRIME_FLAG_SAFE;
        else if (strcmp(argv[i], "--top2") == 0) flags |= CFX_PRIME_FLAG_TOP2;
        else {
            fprintf(stderr, "error: unknown option '%s'\n", argv[i]);
            usage(argv[0]);
            return 1;
        }
    }

    cfx_big_t p; // assuming your cfx_big_copy() allocates dest; if not, call your init
    int rc = cfx_big_gen_nbit_prime(&p, bits, urandom_rng, NULL, flags);
    if (rc != 0) {
        fprintf(stderr, "error: cfx_big_gen_nbit_prime failed (rc=%d)\n", rc);
        return 2;
    }

    // print results
    print_big_hex(&p);
    // optional: also show bit length
    printf("// bits = %zu\n", cfx_big_bitlen(&p));

    // if your API requires an explicit free/clear, call it here (e.g., cfx_big_free(&p);)
    return 0;
}
