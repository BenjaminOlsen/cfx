/* bench_poly_mac.cpp — throughput comparison: Poly1305 vs Poly1271 (portable & AVX2) */

#include <benchmark/benchmark.h>
#include <cstdint>
#include <cstring>
#include <vector>

#include "cfx/poly1305.h"
#include "cfx/poly1271.h"

#if defined(_MSC_VER)
#include <intrin.h>
#else
#include <x86intrin.h>
#endif

static void fill_data(std::vector<uint8_t>& data) {
    for (size_t i = 0; i < data.size(); i++) {
        data[i] = static_cast<uint8_t>(i & 0xff);
    }
}

/* ── Poly1305 (streaming) ── */
static void BM_Poly1305(benchmark::State& state) {
    const size_t size = static_cast<size_t>(state.range(0));
    std::vector<uint8_t> data(size);
    fill_data(data);

    uint8_t key[32];
    for (int i = 0; i < 32; i++) key[i] = static_cast<uint8_t>(i);

    uint8_t tag[16];
    uint64_t tsc = 0;

    for (auto _ : state) {
        unsigned aux;
        uint64_t t0 = __rdtscp(&aux);

        cfx_poly1305_ctx_t ctx;
        cfx_poly1305_init(&ctx, key);
        cfx_poly1305_update(&ctx, data.data(), data.size());
        cfx_poly1305_finish(&ctx, tag);

        tsc += __rdtscp(&aux) - t0;
        benchmark::DoNotOptimize(tag);
        benchmark::ClobberMemory();
    }

    double bytes = static_cast<double>(state.iterations()) * static_cast<double>(size);
    state.SetBytesProcessed(static_cast<int64_t>(bytes));
    state.counters["cpb"] = static_cast<double>(tsc) / bytes;
}

/* ── Poly1271 portable (streaming) ── */
static void BM_Poly1271(benchmark::State& state) {
    const size_t size = static_cast<size_t>(state.range(0));
    std::vector<uint8_t> data(size);
    fill_data(data);

    uint8_t key[32];
    for (int i = 0; i < 32; i++) key[i] = static_cast<uint8_t>(i);

    uint8_t tag[16];
    uint64_t tsc = 0;

    for (auto _ : state) {
        unsigned aux;
        uint64_t t0 = __rdtscp(&aux);

        cfx_poly1271_ctx_t ctx;
        cfx_poly1271_init(&ctx, key);
        cfx_poly1271_update(&ctx, data.data(), data.size());
        cfx_poly1271_finish(&ctx, tag);

        tsc += __rdtscp(&aux) - t0;
        benchmark::DoNotOptimize(tag);
        benchmark::ClobberMemory();
    }

    double bytes = static_cast<double>(state.iterations()) * static_cast<double>(size);
    state.SetBytesProcessed(static_cast<int64_t>(bytes));
    state.counters["cpb"] = static_cast<double>(tsc) / bytes;
}

#if CFX_HAVE_AVX2
/* ── Poly1271 AVX2 (streaming) ── */
static void BM_Poly1271_AVX2(benchmark::State& state) {
    const size_t size = static_cast<size_t>(state.range(0));
    std::vector<uint8_t> data(size);
    fill_data(data);

    uint8_t key[32];
    for (int i = 0; i < 32; i++) key[i] = static_cast<uint8_t>(i);

    uint8_t tag[16];
    uint64_t tsc = 0;

    for (auto _ : state) {
        unsigned aux;
        uint64_t t0 = __rdtscp(&aux);

        cfx_poly1271_avx2_ctx_t ctx;
        cfx_poly1271_avx2_init(&ctx, key);
        cfx_poly1271_avx2_update(&ctx, data.data(), data.size());
        cfx_poly1271_avx2_finish(&ctx, tag);

        tsc += __rdtscp(&aux) - t0;
        benchmark::DoNotOptimize(tag);
        benchmark::ClobberMemory();
    }

    double bytes = static_cast<double>(state.iterations()) * static_cast<double>(size);
    state.SetBytesProcessed(static_cast<int64_t>(bytes));
    state.counters["cpb"] = static_cast<double>(tsc) / bytes;
}
#endif

/* ── One-shot variants ── */
static void BM_Poly1305_OneShot(benchmark::State& state) {
    const size_t size = static_cast<size_t>(state.range(0));
    std::vector<uint8_t> data(size);
    fill_data(data);

    uint8_t key[32];
    for (int i = 0; i < 32; i++) key[i] = static_cast<uint8_t>(i);

    uint8_t tag[16];
    uint64_t tsc = 0;

    for (auto _ : state) {
        unsigned aux;
        uint64_t t0 = __rdtscp(&aux);
        cfx_poly1305(tag, data.data(), data.size(), key);
        tsc += __rdtscp(&aux) - t0;
        benchmark::DoNotOptimize(tag);
        benchmark::ClobberMemory();
    }

    double bytes = static_cast<double>(state.iterations()) * static_cast<double>(size);
    state.SetBytesProcessed(static_cast<int64_t>(bytes));
    state.counters["cpb"] = static_cast<double>(tsc) / bytes;
}

static void BM_Poly1271_OneShot(benchmark::State& state) {
    const size_t size = static_cast<size_t>(state.range(0));
    std::vector<uint8_t> data(size);
    fill_data(data);

    uint8_t key[32];
    for (int i = 0; i < 32; i++) key[i] = static_cast<uint8_t>(i);

    uint8_t tag[16];
    uint64_t tsc = 0;

    for (auto _ : state) {
        unsigned aux;
        uint64_t t0 = __rdtscp(&aux);
        cfx_poly1271(tag, data.data(), data.size(), key);
        tsc += __rdtscp(&aux) - t0;
        benchmark::DoNotOptimize(tag);
        benchmark::ClobberMemory();
    }

    double bytes = static_cast<double>(state.iterations()) * static_cast<double>(size);
    state.SetBytesProcessed(static_cast<int64_t>(bytes));
    state.counters["cpb"] = static_cast<double>(tsc) / bytes;
}

#if CFX_HAVE_AVX2
static void BM_Poly1271_AVX2_OneShot(benchmark::State& state) {
    const size_t size = static_cast<size_t>(state.range(0));
    std::vector<uint8_t> data(size);
    fill_data(data);

    uint8_t key[32];
    for (int i = 0; i < 32; i++) key[i] = static_cast<uint8_t>(i);

    uint8_t tag[16];
    uint64_t tsc = 0;

    for (auto _ : state) {
        unsigned aux;
        uint64_t t0 = __rdtscp(&aux);
        cfx_poly1271_avx2(tag, data.data(), data.size(), key);
        tsc += __rdtscp(&aux) - t0;
        benchmark::DoNotOptimize(tag);
        benchmark::ClobberMemory();
    }

    double bytes = static_cast<double>(state.iterations()) * static_cast<double>(size);
    state.SetBytesProcessed(static_cast<int64_t>(bytes));
    state.counters["cpb"] = static_cast<double>(tsc) / bytes;
}
#endif

/* data sizes: 64B, 256B, 1KB, 4KB, 16KB, 64KB, 256KB, 1MB */
#define MAC_BENCH_SIZES \
    ->Arg(64) \
    ->Arg(256) \
    ->Arg(1 << 10) \
    ->Arg(4 << 10) \
    ->Arg(16 << 10) \
    ->Arg(64 << 10) \
    ->Arg(256 << 10) \
    ->Arg(1 << 20)

BENCHMARK(BM_Poly1305) MAC_BENCH_SIZES;
BENCHMARK(BM_Poly1271) MAC_BENCH_SIZES;
#if CFX_HAVE_AVX2
BENCHMARK(BM_Poly1271_AVX2) MAC_BENCH_SIZES;
#endif

BENCHMARK(BM_Poly1305_OneShot) MAC_BENCH_SIZES;
BENCHMARK(BM_Poly1271_OneShot) MAC_BENCH_SIZES;
#if CFX_HAVE_AVX2
BENCHMARK(BM_Poly1271_AVX2_OneShot) MAC_BENCH_SIZES;
#endif

BENCHMARK_MAIN();
