#include <benchmark/benchmark.h>
#include <cstdint>
#include <cstring>
#include <vector>

#include "cfx/sha256.h"
#include "cfx/sha3.h"
#include "cfx/blake2.h"

static void fill_data(std::vector<uint8_t>& data) {
    for (size_t i = 0; i < data.size(); i++) {
        data[i] = static_cast<uint8_t>(i & 0xff);
    }
}

static void BM_SHA256(benchmark::State& state) {
    const size_t size = static_cast<size_t>(state.range(0));
    std::vector<uint8_t> data(size);
    fill_data(data);
    uint8_t hash[32];

    for (auto _ : state) {
        cfx_sha256_ctx ctx;
        cfx_sha256_init(&ctx);
        cfx_sha256_update(&ctx, data.data(), data.size());
        cfx_sha256_final(&ctx, hash);

        benchmark::DoNotOptimize(hash);
        benchmark::ClobberMemory();
    }

    state.SetBytesProcessed(static_cast<int64_t>(state.iterations()) * size);
}

static void BM_SHA3_256(benchmark::State& state) {
    const size_t size = static_cast<size_t>(state.range(0));
    std::vector<uint8_t> data(size);
    fill_data(data);
    uint8_t hash[32];

    for (auto _ : state) {
        cfx_sha3_256(hash, data.data(), data.size());

        benchmark::DoNotOptimize(hash);
        benchmark::ClobberMemory();
    }

    state.SetBytesProcessed(static_cast<int64_t>(state.iterations()) * size);
}

static void BM_SHA3_512(benchmark::State& state) {
    const size_t size = static_cast<size_t>(state.range(0));
    std::vector<uint8_t> data(size);
    fill_data(data);
    uint8_t hash[64];

    for (auto _ : state) {
        cfx_sha3_512(hash, data.data(), data.size());

        benchmark::DoNotOptimize(hash);
        benchmark::ClobberMemory();
    }

    state.SetBytesProcessed(static_cast<int64_t>(state.iterations()) * size);
}

static void BM_SHAKE128(benchmark::State& state) {
    const size_t size = static_cast<size_t>(state.range(0));
    std::vector<uint8_t> data(size);
    fill_data(data);
    uint8_t hash[32];

    for (auto _ : state) {
        cfx_shake128(hash, 32, data.data(), data.size());

        benchmark::DoNotOptimize(hash);
        benchmark::ClobberMemory();
    }

    state.SetBytesProcessed(static_cast<int64_t>(state.iterations()) * size);
}

static void BM_SHAKE256(benchmark::State& state) {
    const size_t size = static_cast<size_t>(state.range(0));
    std::vector<uint8_t> data(size);
    fill_data(data);
    uint8_t hash[64];

    for (auto _ : state) {
        cfx_shake256(hash, 64, data.data(), data.size());

        benchmark::DoNotOptimize(hash);
        benchmark::ClobberMemory();
    }

    state.SetBytesProcessed(static_cast<int64_t>(state.iterations()) * size);
}

static void BM_BLAKE2b(benchmark::State& state) {
    const size_t size = static_cast<size_t>(state.range(0));
    std::vector<uint8_t> data(size);
    fill_data(data);
    uint8_t hash[64];

    for (auto _ : state) {
        cfx_blake2b(hash, 64, data.data(), data.size(), nullptr, 0);

        benchmark::DoNotOptimize(hash);
        benchmark::ClobberMemory();
    }

    state.SetBytesProcessed(static_cast<int64_t>(state.iterations()) * size);
}

static void BM_BLAKE2s(benchmark::State& state) {
    const size_t size = static_cast<size_t>(state.range(0));
    std::vector<uint8_t> data(size);
    fill_data(data);
    uint8_t hash[32];

    for (auto _ : state) {
        cfx_blake2s(hash, 32, data.data(), data.size(), nullptr, 0);

        benchmark::DoNotOptimize(hash);
        benchmark::ClobberMemory();
    }

    state.SetBytesProcessed(static_cast<int64_t>(state.iterations()) * size);
}

static void BM_BLAKE2b_Stream(benchmark::State& state) {
    const size_t size = static_cast<size_t>(state.range(0));
    std::vector<uint8_t> data(size);
    fill_data(data);
    uint8_t hash[64];

    for (auto _ : state) {
        cfx_blake2b_ctx_t ctx;
        cfx_blake2b_init(&ctx, 64);
        cfx_blake2b_update(&ctx, data.data(), data.size());
        cfx_blake2b_final(&ctx, hash);

        benchmark::DoNotOptimize(hash);
        benchmark::ClobberMemory();
    }

    state.SetBytesProcessed(static_cast<int64_t>(state.iterations()) * size);
}

static void BM_SHA3_256_Stream(benchmark::State& state) {
    const size_t size = static_cast<size_t>(state.range(0));
    std::vector<uint8_t> data(size);
    fill_data(data);
    uint8_t hash[32];

    for (auto _ : state) {
        cfx_sha3_ctx_t ctx;
        cfx_sha3_256_init(&ctx);
        cfx_sha3_256_update(&ctx, data.data(), data.size());
        cfx_sha3_256_final(&ctx, hash);

        benchmark::DoNotOptimize(hash);
        benchmark::ClobberMemory();
    }

    state.SetBytesProcessed(static_cast<int64_t>(state.iterations()) * size);
}

/* data sizes: 64B, 256B, 1KB, 4KB, 16KB, 64KB, 256KB, 1MB */
#define HASH_BENCH_SIZES \
    ->Arg(64) \
    ->Arg(256) \
    ->Arg(1 << 10) \
    ->Arg(4 << 10) \
    ->Arg(16 << 10) \
    ->Arg(64 << 10) \
    ->Arg(256 << 10) \
    ->Arg(1 << 20)

BENCHMARK(BM_SHA256) HASH_BENCH_SIZES;
BENCHMARK(BM_SHA3_256) HASH_BENCH_SIZES;
BENCHMARK(BM_SHA3_256_Stream) HASH_BENCH_SIZES;
BENCHMARK(BM_SHA3_512) HASH_BENCH_SIZES;
BENCHMARK(BM_SHAKE128) HASH_BENCH_SIZES;
BENCHMARK(BM_SHAKE256) HASH_BENCH_SIZES;
BENCHMARK(BM_BLAKE2b) HASH_BENCH_SIZES;
BENCHMARK(BM_BLAKE2b_Stream) HASH_BENCH_SIZES;
BENCHMARK(BM_BLAKE2s) HASH_BENCH_SIZES;

BENCHMARK_MAIN();
