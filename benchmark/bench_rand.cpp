#include <cstdint>
#include <cstdlib>
#include <vector>

#include <benchmark/benchmark.h>

#include "cfx/rand.h"

#if CFX_HAVE_OPENSSL

#include <openssl/rand.h>

static inline void openssl_bytes(void *buf, size_t len) {
    RAND_bytes(reinterpret_cast<unsigned char*>(buf), (int)len);
}

static const cfx_rand_desc_t g_openssl_desc = {
    "OpenSSL RAND_bytes",
    /* gen32 = */ nullptr,
    /* seed  = */ nullptr,
    /* bytes = */ openssl_bytes
};

#endif


static void BM_RNG(benchmark::State& state, const cfx_rand_desc_t* desc) {
    const int inner_loop = static_cast<int>(state.range(0));
    const size_t bytes_per_iter = static_cast<size_t>(inner_loop) * sizeof(uint32_t);

    std::vector<uint32_t> buf(inner_loop);

    // Fixed seed for reproducibility
    if (desc->seed) {
        desc->seed(0x12345678u);
    }

    for (auto _ : state) {
        desc->bytes(buf.data(), bytes_per_iter);
        benchmark::DoNotOptimize(buf);
    }

    const double iters = static_cast<double>(state.iterations());
    const double bytes = iters * bytes_per_iter;
    state.SetBytesProcessed(static_cast<int64_t>(bytes));
    state.SetItemsProcessed(static_cast<int64_t>(iters * inner_loop));
}


#define BENCH_ARGS(b) b->Arg(1<<6)->Arg(1<<8)->Arg(1<<10)->Arg(1<<12)->Arg(1<<14)

int main(int argc, char** argv) {

    std::vector<const cfx_rand_desc_t*> rngs;
    rngs.reserve(g_rand_gen_cnt + 1);

    for (size_t i = 0; i < g_rand_gen_cnt; ++i) {
        rngs.push_back(&g_rand_gens[i]);
    }

    for (const cfx_rand_desc_t* desc : rngs) {
        // Capture `desc` by value so that each lambda keeps its own pointer.
        auto* b = benchmark::RegisterBenchmark(
            desc->name,
            [desc](benchmark::State& state) { BM_RNG(state, desc); });

        // Generate N words per iteration; tune this to your liking.
        // 1<<10 = 1024 words ~ 4 KiB per iteration
        BENCH_ARGS(b);
    }

    #if CFX_HAVE_OPENSSL
    auto* b = benchmark::RegisterBenchmark(
        "OpenSSL RAND_bytes",
        [](benchmark::State& state) { BM_RNG(state, &g_openssl_desc); } );
    BENCH_ARGS(b);
    #endif


    // .-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
    benchmark::Initialize(&argc, argv);
    if (benchmark::ReportUnrecognizedArguments(argc, argv)) return 1;
    benchmark::RunSpecifiedBenchmarks();
    return 0;
}
