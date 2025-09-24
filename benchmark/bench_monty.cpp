// bench_monty.cpp
// Compare Montgomery modular multiplication vs. plain multiply + reduction.

#include <benchmark/benchmark.h>
#include <cstdint>
#include <random>
#include <cstring>
#include <cassert>
#include <algorithm>
#include <iostream>

#include "cfx/cfx.h"


// ----------------------- helpers -----------------------

static inline void random_fill_limbs(uint64_t* p, size_t n, std::mt19937_64& rng) {
    for (size_t i = 0; i < n; ++i) p[i] = rng();
}

static void make_random_big(cfx_big_t* out, size_t limbs, std::mt19937_64& rng) {
    std::vector<uint64_t> vec(limbs, 0);
    std::generate(vec.begin(), vec.end(), [&rng]() -> uint64_t {return rng();});
    cfx_big_from_limbs(out, vec.data(), vec.size());
}

static void make_random_odd_modulus(cfx_big_t* n, size_t limbs, std::mt19937_64& rng) {
    make_random_big(n, limbs, rng);
    // ensure odd
    n->limb[0] |= 1ull;
    // avoid degenerate tiny n
    if (limbs == 1 && n->limb[0] <= 3) n->limb[0] = 5;
}


static void modmul_plain(cfx_big_t* out, const cfx_big_t* a, const cfx_big_t* b, const cfx_big_t* n) {
    cfx_big_t prod;
    cfx_big_init(&prod);
    cfx_big_copy(&prod, a);
    cfx_big_mul(&prod, b);
    cfx_big_mod(out, &prod, n);
    cfx_big_free(&prod);
}

static void modmul_auto(cfx_big_t* out, const cfx_big_t* a, const cfx_big_t* b, const cfx_big_t* n) {
    cfx_big_t prod;
    cfx_big_init(&prod);
    cfx_big_copy(&prod, a);
    cfx_big_mul_auto(&prod, b);
    cfx_big_mod(out, &prod, n);
    cfx_big_free(&prod);
}

// ----------------------- benchmarks -----------------------

static void BM_MontMul(benchmark::State& state) {
    const size_t limbs = static_cast<size_t>(state.range(0));

    std::mt19937_64 rng(0xC0FFEEu);
    
    cfx_big_mont_ctx_t ctx;
    cfx_big_t n;
    cfx_big_init(&n);
    make_random_odd_modulus(&n, limbs, rng);
    int ok = cfx_big_mont_ctx_init(&ctx, &n);
    assert(ok && "mont ctx init failed");

    // random a,b and their Montgomery forms
    cfx_big_t a, b, aM, bM, outM;
    cfx_big_init(&a);
    cfx_big_init(&b);
    cfx_big_init(&aM);
    cfx_big_init(&bM);
    cfx_big_init(&outM);

    make_random_big(&a, limbs, rng);
    make_random_big(&b, limbs, rng);

    ok = cfx_big_mont_to(&aM, &a, &ctx); assert(ok);
    ok = cfx_big_mont_to(&bM, &b, &ctx); assert(ok);

    // Correctness sanity (once, outside timed loop): out = a*b mod n
    {
        cfx_big_t out_check, out_from_mont;
        cfx_big_init(&out_check);
        cfx_big_init(&out_from_mont);

        modmul_plain(&out_check, &a, &b, &n);
        cfx_big_mont_mul(&outM, &aM, &bM, &ctx);
        cfx_big_mont_from(&out_from_mont, &outM, &ctx);

        int c = cfx_big_cmp(&out_check, &out_from_mont);
        if (c != 0) {
            std::cout << "oye!" << std::endl;
        }

        cfx_big_free(&out_check);
        cfx_big_free(&out_from_mont);
    }

    for (auto _ : state) {
        benchmark::DoNotOptimize(outM.n);
        cfx_big_mont_mul(&outM, &aM, &bM, &ctx);
        benchmark::ClobberMemory();
    }

    state.SetItemsProcessed(state.iterations());
    state.counters["limbs"] = static_cast<double>(limbs);
    state.counters["bits"] = static_cast<double>(limbs*64);

    // cleanup
    cfx_big_free(&outM);
    cfx_big_free(&aM);
    cfx_big_free(&bM);
    cfx_big_free(&a);
    cfx_big_free(&b);
    // If you have a ctx_free, call it; else free fields manually:
    cfx_big_mont_ctx_free(&ctx);
    cfx_big_free(&n);
}

static void BM_ModMul_Plain(benchmark::State& state) {

    const size_t limbs = static_cast<size_t>(state.range(0));

    std::mt19937_64 rng(0xBADC0DEu);

    cfx_big_t n, a, b, out;
    cfx_big_init(&n); cfx_big_init(&a); cfx_big_init(&b); cfx_big_init(&out);

    make_random_odd_modulus(&n, limbs, rng);
    make_random_big(&a, limbs, rng);
    make_random_big(&b, limbs, rng);

    for (auto _ : state) {
        benchmark::DoNotOptimize(out.n);
        modmul_plain(&out, &a, &b, &n);
        benchmark::ClobberMemory();
    }

    state.SetItemsProcessed(state.iterations());
    state.counters["limbs"] = static_cast<double>(limbs);
    state.counters["bits"] = static_cast<double>(limbs*64);

    cfx_big_free(&out);
    cfx_big_free(&a);
    cfx_big_free(&b);
    cfx_big_free(&n);
}

static void BM_ModMul_Auto(benchmark::State& state) {

    const size_t limbs = static_cast<size_t>(state.range(0));

    std::mt19937_64 rng(0xBADC0DEu);

    cfx_big_t n, a, b, out;
    cfx_big_init(&n); cfx_big_init(&a); cfx_big_init(&b); cfx_big_init(&out);

    make_random_odd_modulus(&n, limbs, rng);
    make_random_big(&a, limbs, rng);
    make_random_big(&b, limbs, rng);

    for (auto _ : state) {
        benchmark::DoNotOptimize(out.n);
        modmul_auto(&out, &a, &b, &n);
        benchmark::ClobberMemory();
    }

    state.SetItemsProcessed(state.iterations());
    state.counters["limbs"] = static_cast<double>(limbs);
    state.counters["bits"] = static_cast<double>(limbs*64);

    cfx_big_free(&out);
    cfx_big_free(&a);
    cfx_big_free(&b);
    cfx_big_free(&n);
}

// ----------------------- registration -----------------------

BENCHMARK(BM_ModMul_Plain)
    ->Arg(4)->Arg(8)->Arg(16)->Arg(32)->Arg(64)->Arg(128)->Arg(256)->Arg(512)->Arg(1024)->Arg(2048)->Arg(4096);

BENCHMARK(BM_ModMul_Auto)
    ->Arg(4)->Arg(8)->Arg(16)->Arg(32)->Arg(64)->Arg(128)->Arg(256)->Arg(512)->Arg(1024)->Arg(2048)->Arg(4096);

BENCHMARK(BM_MontMul)
    ->Arg(4)->Arg(8)->Arg(16)->Arg(32)->Arg(64)->Arg(128)->Arg(256)->Arg(512)->Arg(1024)->Arg(2048)->Arg(4096);;

BENCHMARK_MAIN();
