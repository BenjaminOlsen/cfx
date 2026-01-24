#include <benchmark/benchmark.h>
#include <cstdint>
#include <cstring>

#include "cfx/chacha20.h"

static void init_key_nonce(uint8_t key[32], uint8_t nonce[12]) {
    for (int i = 0; i < 32; ++i)
        key[i] = static_cast<uint8_t>(i);
    for (int i = 0; i < 12; ++i)
        nonce[i] = static_cast<uint8_t>(100 + i);
}

// ---------------------------------------------------------------------
// Single block using ctx API
// ---------------------------------------------------------------------
static void BM_Chacha20_Block(benchmark::State& state) {
    uint8_t key[32];
    uint8_t nonce[12];
    uint8_t out[64];
    uint32_t counter = 0;

    init_key_nonce(key, nonce);

    cfx_chacha20_ctx_t ctx;
    cfx_chacha20_ctx_init(&ctx, key, nonce);

    for (auto _ : state) {
        cfx_chacha20_block(&ctx, counter, out);
        benchmark::DoNotOptimize(out);
        benchmark::ClobberMemory();
        ++counter;
    }

    state.SetBytesProcessed(static_cast<int64_t>(state.iterations()) * 64);
}

// ---------------------------------------------------------------------
// 4 blocks using ctx API (calls block 4 times)
// ---------------------------------------------------------------------
static void BM_Chacha20_Block_x4(benchmark::State& state) {
    uint8_t key[32];
    uint8_t nonce[12];
    uint8_t out[4][64];
    uint32_t counter = 0;

    init_key_nonce(key, nonce);

    cfx_chacha20_ctx_t ctx;
    cfx_chacha20_ctx_init(&ctx, key, nonce);

    for (auto _ : state) {
        cfx_chacha20_block(&ctx, counter + 0, out[0]);
        cfx_chacha20_block(&ctx, counter + 1, out[1]);
        cfx_chacha20_block(&ctx, counter + 2, out[2]);
        cfx_chacha20_block(&ctx, counter + 3, out[3]);

        benchmark::DoNotOptimize(out);
        benchmark::ClobberMemory();
        counter += 4;
    }

    state.SetBytesProcessed(static_cast<int64_t>(state.iterations()) * 4 * 64);
}

// ---------------------------------------------------------------------
// block4 (SIMD on AVX2, scalar fallback otherwise)
// ---------------------------------------------------------------------
static void BM_Chacha20_Block4(benchmark::State& state) {
    uint8_t key[32];
    uint8_t nonce[12];
    uint8_t out[4][64];
    uint32_t counter = 0;

    init_key_nonce(key, nonce);

    cfx_chacha20_ctx_t ctx;
    cfx_chacha20_ctx_init(&ctx, key, nonce);

    for (auto _ : state) {
        cfx_chacha20_block4(&ctx, counter, out);

        benchmark::DoNotOptimize(out);
        benchmark::ClobberMemory();
        counter += 4;
    }

    state.SetBytesProcessed(static_cast<int64_t>(state.iterations()) * 4 * 64);
}

// ---------------------------------------------------------------------
// block8 (AVX2 8-lane SIMD, or 8x scalar fallback)
// ---------------------------------------------------------------------
static void BM_Chacha20_Block8(benchmark::State& state) {
    uint8_t key[32];
    uint8_t nonce[12];
    uint8_t out[8][64];
    uint32_t counter = 0;

    init_key_nonce(key, nonce);

    cfx_chacha20_ctx_t ctx;
    cfx_chacha20_ctx_init(&ctx, key, nonce);

    for (auto _ : state) {
        cfx_chacha20_block8(&ctx, counter, out);

        benchmark::DoNotOptimize(out);
        benchmark::ClobberMemory();
        counter += 8;
    }

    state.SetBytesProcessed(static_cast<int64_t>(state.iterations()) * 8 * 64);
}

// ---------------------------------------------------------------------
// encrypt (streaming, uses block internally)
// ---------------------------------------------------------------------
static void BM_Chacha20_Encrypt(benchmark::State& state) {
    const size_t len = static_cast<size_t>(state.range(0));
    uint8_t key[32];
    uint8_t nonce[12];

    init_key_nonce(key, nonce);

    std::vector<uint8_t> pt(len, 0x42);
    std::vector<uint8_t> ct(len);

    for (auto _ : state) {
        cfx_chacha20_encrypt(key, 0, nonce, pt.data(), len, ct.data());

        benchmark::DoNotOptimize(ct.data());
        benchmark::ClobberMemory();
    }

    state.SetBytesProcessed(static_cast<int64_t>(state.iterations()) * static_cast<int64_t>(len));
}

// ---------------------------------------------------------------------
// encrypt_ctx (streaming with pre-initialized ctx)
// ---------------------------------------------------------------------
static void BM_Chacha20_Encrypt_Ctx(benchmark::State& state) {
    const size_t len = static_cast<size_t>(state.range(0));
    uint8_t key[32];
    uint8_t nonce[12];

    init_key_nonce(key, nonce);

    cfx_chacha20_ctx_t ctx;
    cfx_chacha20_ctx_init(&ctx, key, nonce);

    std::vector<uint8_t> pt(len, 0x42);
    std::vector<uint8_t> ct(len);

    for (auto _ : state) {
        uint32_t counter = 0;
        cfx_chacha20_encrypt_ctx(&ctx, &counter, pt.data(), len, ct.data());

        benchmark::DoNotOptimize(ct.data());
        benchmark::ClobberMemory();
    }

    state.SetBytesProcessed(static_cast<int64_t>(state.iterations()) * static_cast<int64_t>(len));
}

BENCHMARK(BM_Chacha20_Block);
BENCHMARK(BM_Chacha20_Block_x4);
BENCHMARK(BM_Chacha20_Block4);
BENCHMARK(BM_Chacha20_Block8);
BENCHMARK(BM_Chacha20_Encrypt)->RangeMultiplier(4)->Range(64, 1 << 16);
BENCHMARK(BM_Chacha20_Encrypt_Ctx)->RangeMultiplier(4)->Range(64, 1 << 16);

BENCHMARK_MAIN();
