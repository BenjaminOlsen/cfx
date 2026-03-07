#include "cfx/rand.h"
#include "cfx/rng_splitmix.h"
#include "cfx/chacha20.h"
#include "cfx/memory.h"
#include "cfx/macros.h"
#include "cfx/arith.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#if defined(_WIN32) || defined(_WIN64)
#  define WIN32_LEAN_AND_MEAN
#  include <windows.h>
#  include <bcrypt.h>
#  ifdef _MSC_VER
#    pragma comment(lib, "bcrypt.lib")
#  endif
#else
#  include <unistd.h>
#  include <fcntl.h>
#  include <errno.h>
#endif

/* ---------------------------------------------------------------------------------------------- */
/* chacha20 based rng */

#define CHACHA20_BLOCK_BYTES 64
#if CFX_HAVE_AVX2
#define CHACHA20_LANE_CNT 8
#elif CFX_SIMD
#define CHACHA20_LANE_CNT 4
#else
#define CHACHA20_LANE_CNT 1
#endif

#define CHACHA20_BUF_BYTES (CHACHA20_BLOCK_BYTES * CHACHA20_LANE_CNT)

typedef struct {
    CFX_ALIGNAS(CFX_CHACHA20_CTX_ALIGN) uint8_t buf[CHACHA20_BUF_BYTES];
    cfx_chacha20_ctx_t ctx;
    uint32_t counter;
    size_t idx;
    int seeded;
} cfx_chacha20_rng_t;

CFX_STATIC_ASSERT(sizeof(cfx_chacha20_rng_t) <= CFX_CHACHA_RNG_CTX_SIZE,
    chacha_rng_ctx_too_small);

static inline void chacha20_advance_counter(cfx_chacha20_rng_t *st, uint32_t blocks) {
    uint32_t old = st->counter;
    st->counter += blocks;
    if (st->counter < old) {
        /* counter wrapped - increment nonce to avoid reuse */
        cfx_chacha20_ctx_inc_nonce(&st->ctx);
    }
}


static cfx_chacha20_rng_t G = {0};

CFX_INLINE void cfx_chacha20_refill(cfx_chacha20_rng_t *st) {
    if (!st->seeded) return;

#if CFX_HAVE_AVX2
    cfx_chacha20_block8(&st->ctx, st->counter, (uint8_t (*)[64]) st->buf);
    chacha20_advance_counter(st, 8);
#elif CFX_SIMD
    cfx_chacha20_block4(&st->ctx, st->counter, (uint8_t (*)[64]) st->buf);
    chacha20_advance_counter(st, 4);
#else
    cfx_chacha20_block(&st->ctx, st->counter, st->buf);
    chacha20_advance_counter(st, 1);
#endif
    st->idx = 0;
}


CFX_INLINE void cfx_chacha20_generate(cfx_chacha20_rng_t *st, uint8_t *out, size_t len) {
    assert(st->seeded);
    if (!st->seeded) return;

    while (len && st->idx < CHACHA20_BUF_BYTES) {
        size_t avail = CHACHA20_BUF_BYTES - st->idx;
        size_t n = (len < avail ? len : avail);
        memcpy(out, st->buf + st->idx, n);
        st->idx += n;
        out += n;
        len -= n;
    }

#if CFX_HAVE_AVX2
    while (len >= CHACHA20_BUF_BYTES) {
        cfx_chacha20_block8(&st->ctx, st->counter, (uint8_t (*)[64]) out);
        chacha20_advance_counter(st, 8);
        out += CHACHA20_BUF_BYTES;
        len -= CHACHA20_BUF_BYTES;
    }
#elif CFX_SIMD
    while (len >= CHACHA20_BUF_BYTES) {
        cfx_chacha20_block4(&st->ctx, st->counter, (uint8_t (*)[64]) out);
        chacha20_advance_counter(st, 4);
        out += CHACHA20_BUF_BYTES;
        len -= CHACHA20_BUF_BYTES;
    }
#else
    while (len >= 64) {
        cfx_chacha20_block(&st->ctx, st->counter, out);
        chacha20_advance_counter(st, 1);
        out += 64;
        len -= 64;
    }
#endif

    if (len) {
        cfx_chacha20_refill(st);
        memcpy(out, st->buf, len);
        st->idx = len;
    }
}

static void splitmix64_bytes(uint64_t seed, uint8_t *out, size_t len) {
    uint64_t s = seed ? seed : 1;
    while (len >= 8) {
        uint64_t x = cfx_splitmix64(&s);
        memcpy(out, &x, 8);
        out += 8;
        len -= 8;
    }
    if (len > 0) {
        uint64_t x = cfx_splitmix64(&s);
        memcpy(out, &x, len);
    }
}

void cfx_chacha20_rng_init(cfx_rng_ctx_t *st, uint32_t seed) {
    cfx_chacha20_rng_t *s = (cfx_chacha20_rng_t *)st;
    uint8_t key[32], nonce[12];
    uint64_t sm_seed = seed ? (uint64_t)seed : 1;
    splitmix64_bytes(sm_seed, key, 32);
    sm_seed = cfx_splitmix64(&sm_seed);
    splitmix64_bytes(sm_seed, nonce, 12);
    cfx_chacha20_ctx_init(&s->ctx, key, nonce);
    s->counter = 0;
    s->idx = CHACHA20_BUF_BYTES;
    s->seeded = 1;
    CFX_MEMZERO_S(key, 32);
    CFX_MEMZERO_S(nonce, 12);
}

int cfx_chacha20_rng(void *ctx, uint8_t *out, size_t len) {
    cfx_chacha20_rng_t *s = (cfx_chacha20_rng_t *)ctx;
    if (!s->seeded) return 1;
    cfx_chacha20_generate(s, (uint8_t *)out, len);
    return 0;
}

void cfx_chacha20_seed(uint32_t seed) {
    uint8_t key[32], nonce[12];
    uint64_t sm_seed = seed ? (uint64_t)seed : 1;
    splitmix64_bytes(sm_seed, key, 32);
    sm_seed = cfx_splitmix64(&sm_seed);
    splitmix64_bytes(sm_seed, nonce, 12);
    cfx_chacha20_ctx_init(&G.ctx, key, nonce);
    G.counter = 0;
    G.idx = CHACHA20_BUF_BYTES;
    G.seeded = 1;
    CFX_MEMZERO_S(key, 32);
    CFX_MEMZERO_S(nonce, 12);
}

void cfx_chacha20_bytes(void *buf, size_t len) {
    cfx_chacha20_generate(&G, buf, len);
}

uint32_t cfx_chacha20_gen32(void) {
    uint8_t b[4];
    cfx_chacha20_generate(&G, b, 4);
    return (uint32_t)b[0]
           | ((uint32_t)b[1] << 8)
           | ((uint32_t)b[2] << 16)
           | ((uint32_t)b[3] << 24);
}


/* ---------------------------------------------------------------------------------------------- */
/* cfx's chosen rand */
void cfx_srand(uint32_t seed) {
    cfx_chacha20_seed(seed);
}

uint32_t cfx_urand(void) {
    return cfx_chacha20_gen32();
}

int cfx_rand(void) {
    return (int)(cfx_urand() & 0x7fffffff);
}

void cfx_rand_bytes(void *buf, size_t len) {
    cfx_chacha20_bytes(buf, len);
}

int cfx_rng(void *ctx, uint8_t *out, size_t len) {
    return cfx_chacha20_rng(ctx, out, len);
}

/* ---------------------------------------------- */
/* OS-specific cryptographic RNG */
#if defined(_WIN32) || defined(_WIN64)
/* Windows: BCryptGenRandom */
static int os_getrandom(void *out, size_t len) {
    NTSTATUS status = BCryptGenRandom(NULL, (PUCHAR)out, (ULONG)len,
        BCRYPT_USE_SYSTEM_PREFERRED_RNG);
    return BCRYPT_SUCCESS(status) ? 0 : -1;
}
#else
/* POSIX: /dev/urandom */
static int os_getrandom(void *out, size_t len) {
    int fd = open("/dev/urandom", O_RDONLY);
    if (fd < 0) return -1;
    size_t got = 0;
    while (got < len) {
        ssize_t r = read(fd, (unsigned char *)out + got, len - got);
        if (r < 0) {
            if (errno == EINTR) continue;
            close(fd);
            return -1;
        }
        if (r == 0) {
            close(fd); return -1;
        }
        got += (size_t)r;
    }
    close(fd);
    return 0;
}
#endif

void cfx_srand_os(void) {
    uint8_t tmp[32 + 12 + 4];
    if (os_getrandom(tmp, sizeof(tmp)) != 0) abort();
    cfx_chacha20_ctx_init(&G.ctx, tmp, tmp + 32);
    memcpy(&G.counter, tmp + 44, sizeof G.counter);
    cfx_memzero_s(tmp, sizeof(tmp));
    G.idx = CHACHA20_BUF_BYTES;
    G.seeded = 1;
}

uint32_t cfx_rand_os(void) {
    uint32_t x = 0;
    if (os_getrandom(&x, sizeof x) != 0) {
        abort();
    }
    return x;
}

void cfx_rand_bytes_os(void *buf, size_t len) {
    if (os_getrandom(buf, len) != 0) {
        abort();
    }
}

cfx_limb_t cfx_rand_limb(void) {
#if CFX_LIMB_BITS == 64
    uint64_t lo = cfx_urand();
    uint64_t hi = cfx_urand();
    return (cfx_limb_t)((hi << 32) | lo);
#else
    return (cfx_limb_t)cfx_urand();
#endif
}
