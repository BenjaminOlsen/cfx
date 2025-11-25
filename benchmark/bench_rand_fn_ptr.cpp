#include <benchmark/benchmark.h>
#include <cstdint>

//duplicated from cfx
static inline uint32_t pcg32_core(uint64_t* s) {
    const uint64_t pcg_inc = UINT64_C(0xda3e39cb94b95bdb);
    uint64_t oldstate = *s;
    *s = oldstate * UINT64_C(6364136223846793005) + (pcg_inc | 1);
    uint32_t xorshifted = (uint32_t)(((oldstate >> 18u) ^ oldstate) >> 27u);
    uint32_t rot = oldstate >> 59u;
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}

static uint64_t g_pcg_state_inline  = 0x123456789abcdef0ULL;
static uint64_t g_pcg_state_fnptr   = 0x123456789abcdef0ULL;

// ---------------------------------------------------------------------
// Variant 1: direct / inline call
static void BM_PCG32_Inline(benchmark::State& state) {
    const int inner_loop = static_cast<int>(state.range(0));
    uint64_t s = g_pcg_state_inline;
    uint32_t sink = 0;

    for (auto _ : state) {
        for (int i = 0; i < inner_loop; ++i) {
            sink ^= pcg32_core(&s);   // direct, inlinable call
        }
    }

    // write back state so optimizer can't assume it is unused
    g_pcg_state_inline = s;
    benchmark::DoNotOptimize(sink);
    state.SetItemsProcessed(static_cast<int64_t>(state.iterations()) * inner_loop);
}

// ---------------------------------------------------------------------
// Variant 2: call via function pointer
using rng32_fn = uint32_t(*)(void*);

// Make the core accessible via a function taking void* state.
// Mark noinline to avoid compiler "helpfully" inlining anyway.
#if defined(__GNUC__) || defined(__clang__)
__attribute__((noinline))
#endif
static uint32_t pcg32_core_void(void* p) {
    uint64_t* s = static_cast<uint64_t*>(p);
    return pcg32_core(s);
}

static void BM_PCG32_FnPtr(benchmark::State& state) {
    const int inner_loop = static_cast<int>(state.range(0));
    uint64_t s = g_pcg_state_fnptr;
    uint32_t sink = 0;

    rng32_fn gen32 = &pcg32_core_void;

    for (auto _ : state) {
        for (int i = 0; i < inner_loop; ++i) {
            sink ^= gen32(&s);   // indirect call through function pointer
        }
    }

    g_pcg_state_fnptr = s;
    benchmark::DoNotOptimize(sink);
    state.SetItemsProcessed(static_cast<int64_t>(state.iterations()) * inner_loop);
}

BENCHMARK(BM_PCG32_Inline)->Arg(1 << 10)->Arg(1 << 15);
BENCHMARK(BM_PCG32_FnPtr)->Arg(1 << 10)->Arg(1 << 15);

BENCHMARK_MAIN();
