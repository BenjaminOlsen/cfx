
#include "cfx/arith.h"
#include "cfx/algo.h"
#include <benchmark/benchmark.h>

extern "C" {
static inline cfx_limb_t isqrt_u64(cfx_limb_t n) {
    /* Integer sqrt via double for speed, then adjust */
    double d = (double)n;
    cfx_limb_t x = (cfx_limb_t)(d > 0 ? __builtin_floor(__builtin_sqrt(d)) : 0);
    while ((x + 1) > 0 && (x + 1) <= n / (x+1)) ++x;
    while (x > 0 && x > n / x) --x;
    return x;
}
}

#define BM_SQRT(func)                               \
static void BM_##func(benchmark::State& state) {    \
    cfx_limb_t n = 0xABCD;                          \
    cfx_limb_t sqrt;                                \
    for (auto _ : state) {                          \
        benchmark::DoNotOptimize(sqrt);             \
        sqrt = func(n);                             \
        benchmark::ClobberMemory();                 \
    }                                               \
    state.SetItemsProcessed(state.iterations());    \
}                                                   \
BENCHMARK(BM_##func)

BM_SQRT(isqrt_u64);
#if !defined(CFX_NO_FP)
BM_SQRT(cfx_isqrt_fp);
#endif
BM_SQRT(cfx_isqrt_nr);

BENCHMARK_MAIN();
