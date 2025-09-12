
#include <benchmark/benchmark.h>
#include <cstdint>
#include <cfx/big.h>

static void BM_InitFree(benchmark::State& state) {
    for (auto _ : state) {
        cfx_big_t b;
        cfx_big_init(&b);
        cfx_big_set_val(&b, 42);
        cfx_big_free(&b);
        benchmark::DoNotOptimize(b);
    }
}
BENCHMARK(BM_InitFree);

static void BM_AddSmall(benchmark::State& state) {
    cfx_big_t b;
    cfx_big_init(&b);
    cfx_big_set_val(&b, 0);
    for (auto _ : state) {
        cfx_big_add_sm(&b, 123456789ULL);
        benchmark::DoNotOptimize(b);
    }
    cfx_big_free(&b);
}
BENCHMARK(BM_AddSmall);

static void BM_MulSmall(benchmark::State& state) {
    cfx_big_t b;
    cfx_big_init(&b);
    for (auto _ : state) {
        cfx_big_set_val(&b, 123456789ULL);
        cfx_big_mul_sm(&b, 7);
    }
    benchmark::DoNotOptimize(b);
    cfx_big_free(&b);
}
BENCHMARK(BM_MulSmall);

static void BM_MulBig(benchmark::State& state) {
    cfx_big_t a, b;
    cfx_big_init(&a);
    cfx_big_init(&b);
    
    cfx_big_set_val(&b, 987654321ULL);
    for (auto _ : state) {
        cfx_big_set_val(&a, 123456789ULL);
        cfx_big_mul_eq(&a, &b);
        benchmark::DoNotOptimize(a);
        benchmark::DoNotOptimize(b);
    }
    cfx_big_free(&a);
    cfx_big_free(&b);
}
BENCHMARK(BM_MulBig);

static void BM_DivSm(benchmark::State& state) {
    cfx_big_t a;
    cfx_big_init(&a);
    const char* s = "1239487123698471263948712630489172364019823740192384701237864012837461923876501293874619283746192348751602398476";
    cfx_big_from_str(&a, s);
    
    for (auto _ : state) {
        cfx_big_t tmp;
        cfx_big_init(&tmp);
        cfx_big_copy(&tmp, &a);
        uint64_t rem = cfx_big_div_eq_sm(&tmp, 7);
        benchmark::DoNotOptimize(rem);
        cfx_big_free(&tmp);
    }
    cfx_big_free(&a);
}

BENCHMARK(BM_DivSm);
BENCHMARK_MAIN();
