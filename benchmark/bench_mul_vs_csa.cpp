// bench_mul_vs_csa.cpp

#include <cfx/big.h>
#include <cfx/algo.h>

#include <benchmark/benchmark.h>

#include <cstdint>
#include <vector>
#include <random>
#include <cstring>
#include <cassert>
#include <cstdio>

// ---- helpers ---------------------------------------------------------------
static void ensure_cap(cfx_big_t* b, size_t need) {
    if (b->cap < need) cfx_big_reserve(b, need);
    if (b->n > need)   b->n = need; // caller manages n exactly
}

static void fill_random(cfx_big_t* b, size_t limbs, cfx_limb_t seed) {
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
    InplaceInputs(size_t nb_, size_t nm_, cfx_limb_t seed_base) : nb(nb_), nm(nm_) {
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

    // bt = b0; bt *= m using schoolbook
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

// --- CSA with scratch benchmark (in-place) ---
static void BM_cfx_big_mul_csa_scratch(benchmark::State& state) {
    const size_t limbs = static_cast<size_t>(state.range(0));
    InplaceInputs I(limbs, limbs, /*seed*/0xDEADBEEF);

    // one-time correctness check against your legacy mul
    {
        cfx_big_t tmp; cfx_big_init(&tmp);
        cfx_big_copy(&tmp, &I.b0);
        cfx_big_mul(&tmp, &I.m);  // legacy
        cfx_mul_scratch_t s{};    // temp scratch just for verify
        cfx_mul_scratch_alloc(&s, I.nb + I.nm);
        cfx_mul_scratch_zero(&s, I.nb + I.nm);
        cfx_big_t ref; cfx_big_init(&ref);
        cfx_big_copy(&ref, &I.b0);
        cfx_big_mul_csa_scratch(&ref, &I.m, &s);
        assert(eq_limbs(tmp, ref));
        cfx_mul_scratch_free(&s);
        cfx_big_free(&tmp);
        cfx_big_free(&ref);
    }

    // persistent scratch sized to nout = nb + nm
    cfx_mul_scratch_t scratch{};
    const size_t nout = I.nb + I.nm;
    cfx_mul_scratch_alloc(&scratch, nout);

    for (auto _ : state) {
        // reset inputs and scratch
        cfx_big_copy(&I.bt, &I.b0);
        cfx_mul_scratch_zero(&scratch, nout);

        // in-place multiply using CSA + scratch
        cfx_big_mul_csa_scratch(&I.bt, &I.m, &scratch);

        benchmark::DoNotOptimize(I.bt.limb[0]);
        benchmark::ClobberMemory();
    }

    cfx_mul_scratch_free(&scratch);
    state.SetItemsProcessed(state.iterations());
    state.counters["limbs"] = static_cast<double>(limbs);
}


// ---- Row multiplying  -------------------------------------------------------
// --- Helpers -----------------------------------------------------------------
static inline cfx_limb_t splitmix64(cfx_limb_t& s) {
    cfx_limb_t z = (s += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

static void big_rand_limbs(cfx_big_t* x, size_t n, cfx_limb_t seed0) {
    cfx_big_reserve(x, n);
    cfx_limb_t s = seed0;
    for (size_t i = 0; i < n; ++i) x->limb[i] = splitmix64(s);
    x->n = n;
    if (n && x->limb[n-1] == 0) x->limb[n-1] = 1; // avoid trim-to-zero edge
}

static void big_copy_from(const cfx_big_t* src, cfx_big_t* dst) {
    cfx_big_reserve(dst, src->n);
    if (src->n) std::memcpy(dst->limb, src->limb, src->n * sizeof(cfx_limb_t));
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

static void BM_MulAuto(benchmark::State& state) {
    const size_t na = static_cast<size_t>(state.range(0));
    const size_t nb = static_cast<size_t>(state.range(1));

    // Prepare deterministic operands once (outside timing loop)
    cfx_big_t a, m;
    cfx_big_init(&a); cfx_big_init(&m);
    big_rand_limbs(&a, na, 0xA5A5A5A5DEADBEEFULL);
    big_rand_limbs(&m, nb, 0x0123456789ABCDEFULL);

    // Label & counters
    state.counters["na"] = static_cast<double>(na);
    state.counters["nb"] = static_cast<double>(nb);
    state.SetItemsProcessed(static_cast<int64_t>(state.iterations()) * static_cast<int64_t>(na) * static_cast<int64_t>(nb));

    for (auto _ : state) {
        cfx_big_t tmp;
        cfx_big_init(&tmp);
        big_copy_from(&a, &tmp);

        cfx_big_mul_auto(&tmp, &m);

        // Prevent the optimizer from discarding work
        if (tmp.n) benchmark::DoNotOptimize(tmp.limb[tmp.n - 1]);
        benchmark::ClobberMemory();

        cfx_big_free(&tmp);
    }

    cfx_big_free(&a);
    cfx_big_free(&m);
}

// BENCHMARK(BM_MulRows)
//     ->Args({64,   64,   1})
//     ->Args({64,   64,   8})
//     ->Args({256,  256,  8})
//     ->Args({512,  512,  8})
//     ->Args({1024, 1024, 8})
//     ->Args({2048, 2048, 8})
//     ->Args({4096, 4096, 8})
//     ->Args({8192, 8192, 8})
//     ->Args({8192, 4096, 8})
//     ->Args({16384, 16384, 1})
//     ->Args({16384, 16384, 2})
//     ->Args({16384, 16384, 4})
//     ->Args({16384, 16384, 8})
//     ->Args({16384, 16384, 16})
//     // ->Unit(benchmark::kMillisecond)
//     ->UseRealTime(); // wall clock is more representative for internal threading


// ---- Register sizes --------------------------------------------------------
#define REG(sz) \
    BENCHMARK(BM_cfx_big_mul)->Arg(sz); \
    BENCHMARK(BM_cfx_big_mul_csa)->Arg(sz); \
    BENCHMARK(BM_cfx_big_mul_csa_scratch)->Arg(sz); \
    BENCHMARK(BM_MulRows)->Args({sz, sz, 1}); \
    BENCHMARK(BM_MulRows)->Args({sz, sz, 2}); \
    BENCHMARK(BM_MulRows)->Args({sz, sz, 4}); \
    BENCHMARK(BM_MulRows)->Args({sz, sz, 6}); \
    BENCHMARK(BM_MulRows)->Args({sz, sz, 8}); \
    BENCHMARK(BM_MulRows)->Args({sz, sz, 16}); \
    BENCHMARK(BM_MulAuto)->Args({sz, sz})
   

REG(4);     // 256 bits
REG(8);     // 512 bits
REG(16);    // 1024 bits
REG(32);    // 2048 bits
REG(64);    // 4096 bits
REG(128);   // 8192 bits
REG(256);   // 16384 bits
REG(512);   // 32768 bits
REG(1024);   // 65536 bits
REG(2048);   // 131072 bits
REG(4096);   // 262144 bits ---> segfault!
REG(8192);
// REG(16384);
// REG(16384*2);

BENCHMARK_MAIN();