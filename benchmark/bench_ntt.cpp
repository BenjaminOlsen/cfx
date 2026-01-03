
#include "cfx/big.h"
#include "cfx/ntt.h"

#include <benchmark/benchmark.h>
#include <cstdint>
#include <cstring>

// --- Helpers -----------------------------------------------------------------
static inline cfx_limb_t splitmix64(uint64_t& s) {
    uint64_t z = (s += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return static_cast<cfx_limb_t>(z ^ (z >> 31));
}

static void big_rand_limbs(cfx_big_t* x, size_t n, uint64_t seed0) {
    cfx_big_reserve(x, n);
    uint64_t s = seed0;
    for (size_t i = 0; i < n; ++i) x->limb[i] = splitmix64(s);
    x->n = n;
    if (n && x->limb[n-1] == 0) x->limb[n-1] = 1;
}

// --- Schoolbook multiplication -----------------------------------------------
static void BM_MulSchoolbook(benchmark::State& state) {
    const size_t n = static_cast<size_t>(state.range(0));

    cfx_big_t a, b;
    cfx_big_init(&a);
    cfx_big_init(&b);
    big_rand_limbs(&a, n, 0xDEADBEEF12345678ULL);
    big_rand_limbs(&b, n, 0xCAFEBABE87654321ULL);

    state.counters["limbs"] = static_cast<double>(n);

    for (auto _ : state) {
        cfx_big_t tmp;
        cfx_big_init(&tmp);
        cfx_big_copy(&tmp, &a);

        cfx_big_mul_eq(&tmp, &b);

        if (tmp.n) benchmark::DoNotOptimize(tmp.limb[tmp.n - 1]);
        benchmark::ClobberMemory();
        cfx_big_free(&tmp);
    }

    cfx_big_free(&a);
    cfx_big_free(&b);
}

// --- NTT multiplication ------------------------------------------------------
static void BM_MulNTT(benchmark::State& state) {
    const size_t n = static_cast<size_t>(state.range(0));

    cfx_big_t a, b;
    cfx_big_init(&a);
    cfx_big_init(&b);
    big_rand_limbs(&a, n, 0xDEADBEEF12345678ULL);
    big_rand_limbs(&b, n, 0xCAFEBABE87654321ULL);

    state.counters["limbs"] = static_cast<double>(n);

    for (auto _ : state) {
        cfx_big_t out;
        cfx_big_init(&out);

        cfx_big_mul_fft(&out, &a, &b);

        if (out.n) benchmark::DoNotOptimize(out.limb[out.n - 1]);
        benchmark::ClobberMemory();
        cfx_big_free(&out);
    }

    cfx_big_free(&a);
    cfx_big_free(&b);
}

// --- Auto multiplication (picks best method) ---------------------------------
static void BM_MulAuto(benchmark::State& state) {
    const size_t n = static_cast<size_t>(state.range(0));

    cfx_big_t a, b;
    cfx_big_init(&a);
    cfx_big_init(&b);
    big_rand_limbs(&a, n, 0xDEADBEEF12345678ULL);
    big_rand_limbs(&b, n, 0xCAFEBABE87654321ULL);

    state.counters["limbs"] = static_cast<double>(n);

    for (auto _ : state) {
        cfx_big_t tmp;
        cfx_big_init(&tmp);
        cfx_big_copy(&tmp, &a);

        cfx_big_mul_auto(&tmp, &b);

        if (tmp.n) benchmark::DoNotOptimize(tmp.limb[tmp.n - 1]);
        benchmark::ClobberMemory();
        cfx_big_free(&tmp);
    }

    cfx_big_free(&a);
    cfx_big_free(&b);
}

// Register benchmarks at various sizes
// Small sizes where schoolbook should win
BENCHMARK(BM_MulSchoolbook)->Arg(32)->Arg(64)->Arg(128)->Arg(256);
BENCHMARK(BM_MulNTT)->Arg(32)->Arg(64)->Arg(128)->Arg(256);
BENCHMARK(BM_MulAuto)->Arg(32)->Arg(64)->Arg(128)->Arg(256);

BENCHMARK(BM_MulSchoolbook)->Arg(4096)->Arg(6200)->Arg(8192)->Arg(16384)->Arg(32768)->Arg(65536)->Arg(131072);
BENCHMARK(BM_MulNTT)->Arg(4096)->Arg(6200)->Arg(8192)->Arg(16384)->Arg(32768)->Arg(65536)->Arg(131072);
BENCHMARK(BM_MulAuto)->Arg(4096)->Arg(6200)->Arg(8192)->Arg(16384)->Arg(32768)->Arg(65536)->Arg(131072);

BENCHMARK(BM_MulSchoolbook)->Arg(4096)->Arg(6000)->Arg(7000)->Arg(7200)->Arg(7400)->Arg(7600)->Arg(7800)->Arg(8192);
BENCHMARK(BM_MulNTT)->Arg(4096)->Arg(6000)->Arg(7000)->Arg(7200)->Arg(7400)->Arg(7600)->Arg(7800)->Arg(8192);
BENCHMARK(BM_MulAuto)->Arg(4096)->Arg(6000)->Arg(7000)->Arg(7200)->Arg(7400)->Arg(7600)->Arg(7800)->Arg(8192);

BENCHMARK_MAIN();
