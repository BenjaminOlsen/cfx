// bm_mul_rows.cc
// Benchmark for cfx_big_mul_rows_pthreads (row-parallel, single carry pass).
// Build example (adjust include/lib paths):
//   c++ -O3 -DNDEBUG -std=c++17 bm_mul_rows.cc -I./include -L./build -lcfx -lbenchmark -lpthread -o bm_mul_rows

#include "cfx/big.h"

#include <benchmark/benchmark.h>
#include <cstdint>
#include <cstring>
#include <string>


// --- Helpers -----------------------------------------------------------------
static inline uint64_t splitmix64(uint64_t& s) {
    uint64_t z = (s += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

static void big_rand_limbs(cfx_big_t* x, size_t n, uint64_t seed0) {
    cfx_big_reserve(x, n);
    uint64_t s = seed0;
    for (size_t i = 0; i < n; ++i) x->limb[i] = splitmix64(s);
    x->n = n;
    if (n && x->limb[n-1] == 0) x->limb[n-1] = 1; // avoid trim-to-zero edge
}

static void big_copy_from(const cfx_big_t* src, cfx_big_t* dst) {
    cfx_big_reserve(dst, src->n);
    if (src->n) std::memcpy(dst->limb, src->limb, src->n * sizeof(uint64_t));
    dst->n = src->n;
}

// --- Benchmark ---------------------------------------------------------------
// state.range(0) = na, range(1) = nb, range(2) = threads (<=0 => auto)
static void BM_MulRows(benchmark::State& state) {
    const size_t na = static_cast<size_t>(state.range(0));
    const size_t nb = static_cast<size_t>(state.range(1));
    const int threads = static_cast<int>(state.range(2));

    // Prepare deterministic operands once (outside timing loop)
    cfx_big_t a, m;
    cfx_big_init(&a); cfx_big_init(&m);
    big_rand_limbs(&a, na, 0xA5A5A5A5DEADBEEFULL);
    big_rand_limbs(&m, nb, 0x0123456789ABCDEFULL);

    // Label & counters
    state.counters["na"] = static_cast<double>(na);
    state.counters["nb"] = static_cast<double>(nb);
    state.counters["thr"] = static_cast<double>(threads);
    state.SetItemsProcessed(static_cast<int64_t>(state.iterations()) * static_cast<int64_t>(na) * static_cast<int64_t>(nb));

    for (auto _ : state) {
        cfx_big_t tmp;
        cfx_big_init(&tmp);
        big_copy_from(&a, &tmp);

        cfx_big_mul_rows_pthreads(&tmp, &m, threads);

        // Prevent the optimizer from discarding work
        if (tmp.n) benchmark::DoNotOptimize(tmp.limb[tmp.n - 1]);
        benchmark::ClobberMemory();

        cfx_big_free(&tmp);
    }

    cfx_big_free(&a);
    cfx_big_free(&m);
}
BENCHMARK(BM_MulRows)
    ->Args({64,   64,   1})
    ->Args({64,   64,   8})
    ->Args({256,  256,  8})
    ->Args({512,  512,  8})
    ->Args({1024, 1024, 8})
    ->Args({2048, 2048, 8})
    ->Args({4096, 4096, 8})
    ->Args({8192, 8192, 8})
    ->Args({8192, 4096, 8})
    ->Args({16384, 16384, 1})
    ->Args({16384, 16384, 2})
    ->Args({16384, 16384, 4})
    ->Args({16384, 16384, 8})
    ->Args({16384, 16384, 16})
    // ->Unit(benchmark::kMillisecond)
    ->UseRealTime(); // wall clock is more representative for internal threading

// Optional: explore thread auto-detect (-1) vs fixed
static void BM_MulRows_AutoThreads(benchmark::State& state) {
    const size_t na = static_cast<size_t>(state.range(0));
    const size_t nb = static_cast<size_t>(state.range(1));
    const int threads = -1; // auto

    cfx_big_t a, m;
    cfx_big_init(&a); cfx_big_init(&m);
    big_rand_limbs(&a, na, 0xCAFEBABE'12345678ULL);
    big_rand_limbs(&m, nb, 0x0F0F0F0F'F0F0F0F0ULL);

    state.counters["na"] = static_cast<double>(na);
    state.counters["nb"] = static_cast<double>(nb);
    state.counters["thr"] = -1.0;
    state.SetItemsProcessed(static_cast<int64_t>(state.iterations()) * static_cast<int64_t>(na) * static_cast<int64_t>(nb));

    for (auto _ : state) {
        cfx_big_t tmp; cfx_big_init(&tmp);
        big_copy_from(&a, &tmp);
        cfx_big_mul_rows_pthreads(&tmp, &m, threads);
        if (tmp.n) benchmark::DoNotOptimize(tmp.limb[tmp.n - 1]);
        benchmark::ClobberMemory();
        cfx_big_free(&tmp);
    }

    cfx_big_free(&a); cfx_big_free(&m);
}
BENCHMARK(BM_MulRows_AutoThreads)
    ->Args({1024, 1024})
    ->Args({4096, 4096})
    ->Args({16384, 16384})
    // ->Unit(benchmark::kMillisecond)
    ->UseRealTime();

BENCHMARK_MAIN();