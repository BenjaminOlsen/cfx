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

// ---------------------------------------------------------------------
// Benchmark: scalar 4× single-block
// ---------------------------------------------------------------------

static void BM_Chacha20_Block_Scalar4(benchmark::State& state) {
    uint8_t key[32];
    uint8_t nonce[12];
    uint8_t out[4][64];
    uint32_t counter = 0;

    uint8_t nonce4_dummy[4][12];  // not used here, but reused init helper
    init_key_nonce(key, nonce, nonce4_dummy);

    for (auto _ : state) {
        // Process 4 blocks, same total work as one block4 call
        cfx_chacha20_block_rfc8439(key, counter + 0, nonce, out[0]);
        cfx_chacha20_block_rfc8439(key, counter + 1, nonce, out[1]);
        cfx_chacha20_block_rfc8439(key, counter + 2, nonce, out[2]);
        cfx_chacha20_block_rfc8439(key, counter + 3, nonce, out[3]);

        benchmark::DoNotOptimize(out);
        benchmark::ClobberMemory();

        counter += 4;
    }
    const auto total_bytes = static_cast<int64_t>(state.iterations()) * 4 * 64;
    state.SetBytesProcessed(total_bytes);
    state.counters["total_bytes"] = static_cast<double>(total_bytes);
}

// ---------------------------------------------------------------------
// Benchmark: block4 “SIMD” version
// ---------------------------------------------------------------------

static void BM_Chacha20_Block4_Simd(benchmark::State& state) {
    uint8_t key[32];
    uint8_t nonce1[12];
    uint8_t nonce4[4][12];
    uint8_t out[4][64];

    uint32_t counters[4] = {0, 1, 2, 3};

    init_key_nonce(key, nonce1, nonce4);

    for (auto _ : state) {
        cfx_chacha20_block4_simd(key, counters, nonce4, out);

        benchmark::DoNotOptimize(out);
        benchmark::ClobberMemory();

        // so we don’t keep encrypting the same block values
        counters[0] += 4;
        counters[1] += 4;
        counters[2] += 4;
        counters[3] += 4;
    }

    const auto total_bytes = static_cast<int64_t>(state.iterations()) * 4 * 64;
    state.SetBytesProcessed(total_bytes);
    state.counters["total_bytes"] = static_cast<double>(total_bytes);
}

BENCHMARK(BM_Chacha20_Block_Scalar4);
BENCHMARK(BM_Chacha20_Block4_Simd);

BENCHMARK_MAIN();