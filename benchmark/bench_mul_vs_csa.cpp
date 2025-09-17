// bench_mul_vs_csa.cpp

#include <cfx/big.h>

#include <benchmark/benchmark.h>

#include <cstdint>
#include <vector>
#include <random>
#include <cstring>
#include <cassert>

// ---- helpers ---------------------------------------------------------------
static void ensure_cap(cfx_big_t* b, size_t need) {
    if (b->cap < need) cfx_big_reserve(b, need);
    if (b->n > need)   b->n = need; // caller manages n exactly
}

static void fill_random(cfx_big_t* b, size_t limbs, uint64_t seed) {
    cfx_big_reserve(b, limbs);
    std::mt19937_64 rng(seed);
    for (size_t i = 0; i < limbs; ++i) b->limb[i] = rng();
    if (limbs) b->limb[limbs-1] |= (1ULL << 63); // avoid leading zero
    b->n = limbs;
}

static bool eq_limbs(const cfx_big_t& x, const cfx_big_t& y) {
    if (x.n != y.n) return false;
    for (size_t i = 0; i < x.n; ++i) if (x.limb[i] != y.limb[i]) return false;
    return true;
}

struct InplaceInputs {
    cfx_big_t b0{}, m{}, bt{}, out_ref{};
    size_t nb{}, nm{};
    InplaceInputs(size_t nb_, size_t nm_, uint64_t seed_base) : nb(nb_), nm(nm_) {
        cfx_big_init(&b0); cfx_big_init(&m); cfx_big_init(&bt); cfx_big_init(&out_ref);
        fill_random(&b0, nb, seed_base ^ 0xA55A5AA5u);
        fill_random(&m,  nm, seed_base ^ 0xC3C3C3C3u);
        // Work buffers large enough for product
        cfx_big_reserve(&bt,      nb + nm);
        cfx_big_reserve(&out_ref, nb + nm);
    }
    ~InplaceInputs() { cfx_big_free(&b0); cfx_big_free(&m); cfx_big_free(&bt); cfx_big_free(&out_ref); }
};

// One-time correctness check: both in-place paths should match on (b0, m)
static void verify_correctness(InplaceInputs& I) {
    // out_ref = b0; out_ref *= m using CSA
    cfx_big_copy(&I.out_ref, &I.b0);
    cfx_big_mul_csa(&I.out_ref, &I.m);

    // bt = b0; bt *= m using legacy
    cfx_big_copy(&I.bt, &I.b0);
    cfx_big_mul(&I.bt, &I.m);

    if (!eq_limbs(I.out_ref, I.bt)) {
        fprintf(stderr, "[verify] mismatch nb=%zu nm=%zu\n", I.nb, I.nm);
        std::abort();
    }
}

// ---- Bench templates (in-place) --------------------------------------------
template <typename F>
static void bench_inplace(benchmark::State& state, F&& inplace_mul) {
    const size_t limbs = static_cast<size_t>(state.range(0));
    InplaceInputs I(limbs, limbs, /*seed*/0xDEADBEEF);
    verify_correctness(I);  // guard once

    for (auto _ : state) {
        // Reset bt to original multiplicand b0 (no reallocs)
        cfx_big_copy(&I.bt, &I.b0);
        inplace_mul(I.bt, I.m);   // bt *= m

        benchmark::DoNotOptimize(I.bt.limb[0]);
        benchmark::ClobberMemory();
    }

    state.SetItemsProcessed(state.iterations());
    state.counters["limbs"] = static_cast<double>(limbs);
}

static void BM_cfx_big_mul(benchmark::State& state) {
    bench_inplace(state, [](cfx_big_t& b, const cfx_big_t& m){
        cfx_big_mul(&b, &m);
    });
}

static void BM_cfx_big_mul_csa(benchmark::State& state) {
    bench_inplace(state, [](cfx_big_t& b, const cfx_big_t& m){
        cfx_big_mul_csa(&b, &m);
    });
}

// ---- Register sizes --------------------------------------------------------
#define REG(sz) \
    BENCHMARK(BM_cfx_big_mul)->Arg(sz); \
    BENCHMARK(BM_cfx_big_mul_csa)->Arg(sz)

REG(4);     // 256 bits
REG(8);     // 512 bits
REG(16);    // 1024 bits
REG(32);    // 2048 bits
REG(64);    // 4096 bits
REG(128);   // 8192 bits
REG(256);   // 16384 bits

BENCHMARK_MAIN();