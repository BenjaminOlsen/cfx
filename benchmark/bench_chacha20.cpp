#include <benchmark/benchmark.h>
#include <cstdint>
#include <cstring>

#include "cfx/chacha20.h"

static void init_key_nonce(uint8_t key[32], uint8_t nonce1[12], uint8_t nonce4[4][12]) {
    // real values don’t matter
    for (int i = 0; i < 32; ++i) {
        key[i] = static_cast<uint8_t>(i);
    }

    for (int i = 0; i < 12; ++i) {
        nonce1[i] = static_cast<uint8_t>(100 + i);
    }

    for (int b = 0; b < 4; ++b) {
        for (int i = 0; i < 12; ++i) {
            // Use same nonce for all 4 lanes to keep behavior identical
            nonce4[b][i] = nonce1[i];
        }
    }
}

/********************************************************************************************
 Copied here for benchmarking comparison only - alternative chacha20 block functions
 ********************************************************************************************/
#define ROTL32(x, n) CFX_CHACHA20_ROTL32(x, n)
#define QR(a, b, c, d) CFX_CHACHA20_QR(a, b, c, d)

#define _EXPA 0x61707865u
#define _ND_3 0x3320646eu
#define _2_BY 0x79622d32u
#define _TE_K 0x6b206574u

void cfx_chacha20_block_rfc8439_2(const uint8_t key[32], uint32_t counter, const uint8_t nonce[12], uint8_t out[64]) {

    uint32_t s0,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15;
    uint32_t w0,w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13,w14,w15;

    s0  = 0x61707865u;  /* "expa" */
    s1  = 0x3320646eu;  /* "nd 3" */
    s2  = 0x79622d32u;  /* "2-by" */
    s3  = 0x6b206574u;  /* "te k" */

    s4  = CFX_LOAD32_LE(key + 0);
    s5  = CFX_LOAD32_LE(key + 4);
    s6  = CFX_LOAD32_LE(key + 8);
    s7  = CFX_LOAD32_LE(key + 12);
    s8  = CFX_LOAD32_LE(key + 16);
    s9  = CFX_LOAD32_LE(key + 20);
    s10 = CFX_LOAD32_LE(key + 24);
    s11 = CFX_LOAD32_LE(key + 28);

    s12 = counter;                            /* 32-bit block counter */

    s13 = CFX_LOAD32_LE(nonce + 0);
    s14 = CFX_LOAD32_LE(nonce + 4);
    s15 = CFX_LOAD32_LE(nonce + 8);           /* 96-bit nonce */

    w0  = s0 ;
    w1  = s1 ;
    w2  = s2 ;
    w3  = s3 ;
    w4  = s4 ;
    w5  = s5 ;
    w6  = s6 ;
    w7  = s7 ;
    w8  = s8 ;
    w9  = s9 ;
    w10 = s10;
    w11 = s11;
    w12 = s12;
    w13 = s13;
    w14 = s14;
    w15 = s15;

    for (size_t i = 20; i; i-=2){
        /* column rounds */
        QR(w0, w4, w8, w12)
        QR(w1, w5, w9, w13)
        QR(w2, w6, w10, w14)
        QR(w3, w7, w11, w15)
        /* diagonal rounds */
        QR(w0, w5, w10, w15)
        QR(w1, w6, w11, w12)
        QR(w2, w7, w8, w13)
        QR(w3, w4, w9, w14)
    }

    w0  += s0 ;
    w1  += s1 ;
    w2  += s2 ;
    w3  += s3 ;
    w4  += s4 ;
    w5  += s5 ;
    w6  += s6 ;
    w7  += s7 ;
    w8  += s8 ;
    w9  += s9 ;
    w10 += s10;
    w11 += s11;
    w12 += s12;
    w13 += s13;
    w14 += s14;
    w15 += s15;

    CFX_STORE32_LE(out + 0, w0);
    CFX_STORE32_LE(out + 4, w1);
    CFX_STORE32_LE(out + 8, w2);
    CFX_STORE32_LE(out + 12, w3);
    CFX_STORE32_LE(out + 16, w4);
    CFX_STORE32_LE(out + 20, w5);
    CFX_STORE32_LE(out + 24, w6);
    CFX_STORE32_LE(out + 28, w7);
    CFX_STORE32_LE(out + 32, w8);
    CFX_STORE32_LE(out + 36, w9);
    CFX_STORE32_LE(out + 40, w10);
    CFX_STORE32_LE(out + 44, w11);
    CFX_STORE32_LE(out + 48, w12);
    CFX_STORE32_LE(out + 52, w13);
    CFX_STORE32_LE(out + 56, w14);
    CFX_STORE32_LE(out + 60, w15);
}

void cfx_chacha20_block_rfc8439_3(const uint8_t *CFX_RESTRICT key,
                                  uint32_t counter,
                                  const uint8_t *CFX_RESTRICT nonce,
                                  uint8_t *CFX_RESTRICT out) {

    uint32_t x0  = _EXPA;
    uint32_t x1  = _ND_3;
    uint32_t x2  = _2_BY;
    uint32_t x3  = _TE_K;
    uint32_t x4  = CFX_LOAD32_LE(key + 0);
    uint32_t x5  = CFX_LOAD32_LE(key + 4);
    uint32_t x6  = CFX_LOAD32_LE(key + 8);
    uint32_t x7  = CFX_LOAD32_LE(key + 12);
    uint32_t x8  = CFX_LOAD32_LE(key + 16);
    uint32_t x9  = CFX_LOAD32_LE(key + 20);
    uint32_t x10 = CFX_LOAD32_LE(key + 24);
    uint32_t x11 = CFX_LOAD32_LE(key + 28);
    uint32_t x12 = counter;
    uint32_t x13 = CFX_LOAD32_LE(nonce + 0);
    uint32_t x14 = CFX_LOAD32_LE(nonce + 4);
    uint32_t x15 = CFX_LOAD32_LE(nonce + 8);

    uint32_t o0 = x0,  o1 = x1,  o2 = x2,  o3 = x3;
    uint32_t o4 = x4,  o5 = x5,  o6 = x6,  o7 = x7;
    uint32_t o8 = x8,  o9 = x9,  o10 = x10, o11 = x11;
    uint32_t o12 = x12, o13 = x13, o14 = x14, o15 = x15;

    for (unsigned i = 20; i; i -= 2) {
        QR(x0, x4, x8,  x12);
        QR(x1, x5, x9,  x13);
        QR(x2, x6, x10, x14);
        QR(x3, x7, x11, x15);
        QR(x0, x5, x10, x15);
        QR(x1, x6, x11, x12);
        QR(x2, x7, x8,  x13);
        QR(x3, x4, x9,  x14);
    }

    x0  += o0;  x1  += o1;  x2  += o2;  x3  += o3;
    x4  += o4;  x5  += o5;  x6  += o6;  x7  += o7;
    x8  += o8;  x9  += o9;  x10 += o10; x11 += o11;
    x12 += o12; x13 += o13; x14 += o14; x15 += o15;

    CFX_STORE32_LE(out +  0, x0);
    CFX_STORE32_LE(out +  4, x1);
    CFX_STORE32_LE(out +  8, x2);
    CFX_STORE32_LE(out + 12, x3);
    CFX_STORE32_LE(out + 16, x4);
    CFX_STORE32_LE(out + 20, x5);
    CFX_STORE32_LE(out + 24, x6);
    CFX_STORE32_LE(out + 28, x7);
    CFX_STORE32_LE(out + 32, x8);
    CFX_STORE32_LE(out + 36, x9);
    CFX_STORE32_LE(out + 40, x10);
    CFX_STORE32_LE(out + 44, x11);
    CFX_STORE32_LE(out + 48, x12);
    CFX_STORE32_LE(out + 52, x13);
    CFX_STORE32_LE(out + 56, x14);
    CFX_STORE32_LE(out + 60, x15);
}


// ---------------------------------------------------------------------
// scalar 4 x single-block
// ---------------------------------------------------------------------

static void BM_Chacha20_Block_Scalar4(benchmark::State& state) {
    uint8_t key[32];
    uint8_t nonce[12];
    uint8_t out[4][64];
    uint32_t counter = 0;

    uint8_t nonce4_dummy[4][12];   // not used here, but reused init helper
    init_key_nonce(key, nonce, nonce4_dummy);

    for (auto _ : state) {
        for (size_t k = 2; k--;) {
            cfx_chacha20_block_rfc8439(key, counter + 0, nonce, out[0]);
            cfx_chacha20_block_rfc8439(key, counter + 1, nonce, out[1]);
            cfx_chacha20_block_rfc8439(key, counter + 2, nonce, out[2]);
            cfx_chacha20_block_rfc8439(key, counter + 3, nonce, out[3]);

            benchmark::DoNotOptimize(out);
            benchmark::ClobberMemory();

            counter += 4;
        }
    }
    const auto total_bytes = static_cast<int64_t>(state.iterations()) * 4 * 64 * 2;
    state.SetBytesProcessed(total_bytes);
    state.counters["total_bytes"] = static_cast<double>(total_bytes);
}

// ---------------------------------------------------------------------
// same as above, but different
// ---------------------------------------------------------------------
static void BM_Chacha20_Block_Scalar4_2(benchmark::State& state) {
    uint8_t key[32];
    uint8_t nonce[12];
    uint8_t out[4][64];
    uint32_t counter = 0;

    uint8_t nonce4_dummy[4][12];
    init_key_nonce(key, nonce, nonce4_dummy);

    for (auto _ : state) {
        for (size_t k = 2; k--;) {
            cfx_chacha20_block_rfc8439_2(key, counter + 0, nonce, out[0]);
            cfx_chacha20_block_rfc8439_2(key, counter + 1, nonce, out[1]);
            cfx_chacha20_block_rfc8439_2(key, counter + 2, nonce, out[2]);
            cfx_chacha20_block_rfc8439_2(key, counter + 3, nonce, out[3]);

            benchmark::DoNotOptimize(out);
            benchmark::ClobberMemory();

            counter += 4;
        }
    }
    const auto total_bytes = static_cast<int64_t>(state.iterations()) * 4 * 64 * 2;
    state.SetBytesProcessed(total_bytes);
    state.counters["total_bytes"] = static_cast<double>(total_bytes);
}

// ---------------------------------------------------------------------
// same as above, but different
// ---------------------------------------------------------------------
static void BM_Chacha20_Block_Scalar4_3(benchmark::State& state) {
    uint8_t key[32];
    uint8_t nonce[12];
    uint8_t out[4][64];
    uint32_t counter = 0;

    uint8_t nonce4_dummy[4][12];
    init_key_nonce(key, nonce, nonce4_dummy);

    for (auto _ : state) {
        for (size_t k = 2; k--;) {
            cfx_chacha20_block_rfc8439_3(key, counter + 0, nonce, out[0]);
            cfx_chacha20_block_rfc8439_3(key, counter + 1, nonce, out[1]);
            cfx_chacha20_block_rfc8439_3(key, counter + 2, nonce, out[2]);
            cfx_chacha20_block_rfc8439_3(key, counter + 3, nonce, out[3]);

            benchmark::DoNotOptimize(out);
            benchmark::ClobberMemory();

            counter += 4;
        }
    }
    const auto total_bytes = static_cast<int64_t>(state.iterations()) * 4 * 64 * 2;
    state.SetBytesProcessed(total_bytes);
    state.counters["total_bytes"] = static_cast<double>(total_bytes);
}

static void BM_Chacha20_Block_Ctx(benchmark::State& state) {
    uint8_t key[32];
    uint8_t nonce[12];
    uint8_t out[4][64];
    uint32_t counter = 0;

    uint8_t nonce4_dummy[4][12];
    init_key_nonce(key, nonce, nonce4_dummy);

    cfx_chacha20_ctx_t s;
    cfx_chacha20_ctx_init(&s, key, nonce);

    for (auto _ : state) {
        for (size_t k = 2; k--;) {
            cfx_chacha20_block(&s, counter, out[0]);
            cfx_chacha20_block(&s, counter+1, out[1]);
            cfx_chacha20_block(&s, counter+2, out[2]);
            cfx_chacha20_block(&s, counter+3, out[3]);

            benchmark::DoNotOptimize(out);
            benchmark::ClobberMemory();

            counter += 4;
        }
    }
    const auto total_bytes = static_cast<int64_t>(state.iterations()) * 4 * 64 * 2;
    state.SetBytesProcessed(total_bytes);
    state.counters["total_bytes"] = static_cast<double>(total_bytes);
}

// ---------------------------------------------------------------------
// block4 SIMD
// ---------------------------------------------------------------------
static void BM_Chacha20_Block4_Simd(benchmark::State& state) {
    uint8_t key[32];
    uint8_t nonce1[12];
    uint8_t nonce4[4][12];
    uint8_t out[4][64];

    uint32_t counters[4] = {0, 1, 2, 3};

    init_key_nonce(key, nonce1, nonce4);

    for (auto _ : state) {
        for (size_t k = 2; k--;) {
            cfx_chacha20_block4_simd(key, counters, nonce4, out);

            benchmark::DoNotOptimize(out);
            benchmark::ClobberMemory();

            counters[0] += 4;
            counters[1] += 4;
            counters[2] += 4;
            counters[3] += 4;
        }
    }

    const auto total_bytes = static_cast<int64_t>(state.iterations()) * 4 * 64 * 2;
    state.SetBytesProcessed(total_bytes);
    state.counters["total_bytes"] = static_cast<double>(total_bytes);
}

static void BM_Chacha20_Block4_Simd_2(benchmark::State& state) {
    uint8_t key[32];
    uint8_t nonce1[12];
    uint8_t nonce4[4][12];
    uint8_t out[4][64];

    uint32_t counters[4] = {0, 1, 2, 3};

    init_key_nonce(key, nonce1, nonce4);
    cfx_chacha20_ctx4_t s;
    cfx_chacha20_state_init4(&s, key, nonce4);

    for (auto _ : state) {
        for (size_t k = 2; k--;) {
            cfx_chacha20_block4(&s, counters, out);

            benchmark::DoNotOptimize(out);
            benchmark::ClobberMemory();

            counters[0] += 4;
            counters[1] += 4;
            counters[2] += 4;
            counters[3] += 4;
        }
    }

    const auto total_bytes = static_cast<int64_t>(state.iterations()) * 4 * 64 * 2;
    state.SetBytesProcessed(total_bytes);
    state.counters["total_bytes"] = static_cast<double>(total_bytes);
}

static void BM_Chacha20_Block8_avx2(benchmark::State& state) {

    uint8_t key[32];
    uint8_t nonce1[12];
    uint8_t nonce4[4][12];  // wont use this
    uint8_t out[8][64];

    init_key_nonce(key, nonce1, nonce4);

    uint32_t counter = 0;
    for (auto _ : state) {
        cfx_chacha20_block8_avx2(key, counter, nonce1, out);
        benchmark::DoNotOptimize(out);
        benchmark::ClobberMemory();

        counter += 8;
    }

    const auto total_bytes = static_cast<int64_t>(state.iterations()) * 8 * 64;
    state.SetBytesProcessed(total_bytes);
    state.counters["total_bytes"] = static_cast<double>(total_bytes);

}

BENCHMARK(BM_Chacha20_Block_Scalar4)->Arg(1 << 16 );
BENCHMARK(BM_Chacha20_Block_Scalar4_2)->Arg(1 << 16 );
BENCHMARK(BM_Chacha20_Block_Scalar4_3)->Arg(1 << 16 );
BENCHMARK(BM_Chacha20_Block_Ctx)->Arg(1 << 16 );
BENCHMARK(BM_Chacha20_Block4_Simd)->Arg(1 << 16 );
BENCHMARK(BM_Chacha20_Block4_Simd_2)->Arg(1 << 16 );
BENCHMARK(BM_Chacha20_Block8_avx2)->Arg(1 << 16 );

BENCHMARK_MAIN();