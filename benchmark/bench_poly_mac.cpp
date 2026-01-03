#include <benchmark/benchmark.h>
#include <cstdint>
#include <cstring>
#include <vector>

#include "cfx/poly1305.h"
#include "cfx/poly1271.h"

static void fill_data(std::vector<uint8_t>& data) {
    for (size_t i = 0; i < data.size(); i++) {
        data[i] = static_cast<uint8_t>(i & 0xff);
    }
}

static void BM_Poly1305(benchmark::State& state) {
    const size_t size = static_cast<size_t>(state.range(0));
    std::vector<uint8_t> data(size);
    fill_data(data);

    uint8_t key[32];
    for (int i = 0; i < 32; i++) key[i] = static_cast<uint8_t>(i);

    uint8_t tag[16];

    for (auto _ : state) {
        cfx_poly1305_ctx_t ctx;
        cfx_poly1305_init(&ctx, key);
        cfx_poly1305_update(&ctx, data.data(), data.size());
        cfx_poly1305_finish(&ctx, tag);

        benchmark::DoNotOptimize(tag);
        benchmark::ClobberMemory();
    }

    state.SetBytesProcessed(static_cast<int64_t>(state.iterations()) * size);
}

static void BM_Poly1271(benchmark::State& state) {
    const size_t size = static_cast<size_t>(state.range(0));
    std::vector<uint8_t> data(size);
    fill_data(data);

    uint8_t key[32];
    for (int i = 0; i < 32; i++) key[i] = static_cast<uint8_t>(i);

    uint8_t tag[16];

    for (auto _ : state) {
        cfx_poly1271_ctx_t ctx;
        cfx_poly1271_init(&ctx, key);
        cfx_poly1271_update(&ctx, data.data(), data.size());
        cfx_poly1271_finish(&ctx, tag);

        benchmark::DoNotOptimize(tag);
        benchmark::ClobberMemory();
    }

    state.SetBytesProcessed(static_cast<int64_t>(state.iterations()) * size);
}

/* one-shot versions */
static void BM_Poly1305_OneShot(benchmark::State& state) {
    const size_t size = static_cast<size_t>(state.range(0));
    std::vector<uint8_t> data(size);
    fill_data(data);

    uint8_t key[32];
    for (int i = 0; i < 32; i++) key[i] = static_cast<uint8_t>(i);

    uint8_t tag[16];

    for (auto _ : state) {
        cfx_poly1305(key, data.data(), data.size(), tag);

        benchmark::DoNotOptimize(tag);
        benchmark::ClobberMemory();
    }

    state.SetBytesProcessed(static_cast<int64_t>(state.iterations()) * size);
}

static void BM_Poly1271_OneShot(benchmark::State& state) {
    const size_t size = static_cast<size_t>(state.range(0));
    std::vector<uint8_t> data(size);
    fill_data(data);

    uint8_t key[32];
    for (int i = 0; i < 32; i++) key[i] = static_cast<uint8_t>(i);

    uint8_t tag[16];

    for (auto _ : state) {
        cfx_poly1271(tag, data.data(), data.size(), key);

        benchmark::DoNotOptimize(tag);
        benchmark::ClobberMemory();
    }

    state.SetBytesProcessed(static_cast<int64_t>(state.iterations()) * size);
}

/* data sizes: 64B, 256B, 1KB, 4KB, 16KB, 64KB */
#define MAC_BENCH_SIZES \
    ->Arg(64) \
    ->Arg(256) \
    ->Arg(1 << 10) \
    ->Arg(4 << 10) \
    ->Arg(16 << 10) \
    ->Arg(64 << 10)

BENCHMARK(BM_Poly1305) MAC_BENCH_SIZES;
BENCHMARK(BM_Poly1271) MAC_BENCH_SIZES;
BENCHMARK(BM_Poly1305_OneShot) MAC_BENCH_SIZES;
BENCHMARK(BM_Poly1271_OneShot) MAC_BENCH_SIZES;

BENCHMARK_MAIN();
