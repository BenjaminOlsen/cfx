#include <cstdint>
#include <cstdlib>
#include <vector>
#include <cstring>

#include <benchmark/benchmark.h>

#include "cfx/aead_chacha20_poly1305.h"
#include "cfx/macros.h"
#include "cfx/rand.h"

#if defined(_MSC_VER)
#include <intrin.h>
#else
#include <x86intrin.h>
#endif

/* Encrypt benchmark: AEAD( key, nonce, pt, aad ) -> ct, tag */
static void BM_AEAD_Encrypt(benchmark::State& state, size_t aad_len) {
    const size_t pt_len = static_cast<size_t>(state.range(0));

    std::vector<uint8_t> key(32);
    std::vector<uint8_t> nonce(12);
    std::vector<uint8_t> aad(aad_len);
    std::vector<uint8_t> pt(pt_len);
    std::vector<uint8_t> ct(pt_len);
    std::vector<uint8_t> pt_out(pt_len);
    uint8_t tag[16];

    cfx_rand_bytes(key.data(),   key.size());
    cfx_rand_bytes(nonce.data(), nonce.size());
    cfx_rand_bytes(pt.data(),    pt.size());
    if (aad_len) {
        cfx_rand_bytes(aad.data(), aad.size());
    }

    const uint8_t* aad_ptr = aad_len ? aad.data() : nullptr;

    /* quick correctness sanity check before timing */
    int rc = cfx_chacha20_poly1305_encrypt(
        ct.data(), tag,
        pt.data(), pt_len,
        aad_ptr, aad_len,
        key.data(), nonce.data()
    );
    CFX_ASSERT(rc == 0);

    rc = cfx_chacha20_poly1305_decrypt(
        pt_out.data(),
        ct.data(), pt_len,
        aad_ptr, aad_len,
        key.data(), nonce.data(),
        tag
    );
    CFX_ASSERT(rc == 0);
    CFX_ASSERT(std::memcmp(pt_out.data(), pt.data(), pt_len) == 0);

    // measure CPU cycles per byte!
    uint64_t tsc_ticks = 0;

    for (auto _ : state) {
        unsigned aux;
        uint64_t start_cycles = __rdtscp(&aux);

        rc = cfx_chacha20_poly1305_encrypt(
            ct.data(), tag,
            pt.data(), pt_len,
            aad_ptr, aad_len,
            key.data(), nonce.data()
        );

        uint64_t end_cycles = __rdtscp(&aux);
        CFX_ASSERT(rc == 0);

        tsc_ticks += (end_cycles - start_cycles);

        benchmark::DoNotOptimize(ct.data());
        benchmark::DoNotOptimize(tag);
    }

    const double iters = static_cast<double>(state.iterations());
    const double bytes = iters * static_cast<double>(pt_len);
    const double avg_cycles = static_cast<double>(tsc_ticks) / iters;
    const double cpb = static_cast<double>(tsc_ticks) / bytes;

    state.SetBytesProcessed(static_cast<int64_t>(bytes));
    state.SetItemsProcessed(static_cast<int64_t>(iters));

    state.counters["tsc_ticks_per_iter"] = avg_cycles;
    state.counters["tsc_ticks_per_byte"] = cpb;
}

/* Decrypt benchmark: AEAD-Decrypt( key, nonce, ct, tag, aad ) -> pt */
static void BM_AEAD_Decrypt(benchmark::State& state, size_t aad_len) {
    const size_t pt_len = static_cast<size_t>(state.range(0));

    std::vector<uint8_t> key(32);
    std::vector<uint8_t> nonce(12);
    std::vector<uint8_t> aad(aad_len);
    std::vector<uint8_t> pt(pt_len);
    std::vector<uint8_t> ct(pt_len);
    std::vector<uint8_t> pt_out(pt_len);
    uint8_t tag[16];

    cfx_rand_bytes(key.data(),   key.size());
    cfx_rand_bytes(nonce.data(), nonce.size());
    cfx_rand_bytes(pt.data(),    pt.size());
    if (aad_len) {
        cfx_rand_bytes(aad.data(), aad.size());
    }

    const uint8_t* aad_ptr = aad_len ? aad.data() : nullptr;

    /* sanity check */
    int rc = cfx_chacha20_poly1305_encrypt(
        ct.data(), tag,
        pt.data(), pt_len,
        aad_ptr, aad_len,
        key.data(), nonce.data()
    );
    CFX_ASSERT(rc == 0);

    rc = cfx_chacha20_poly1305_decrypt(
        pt_out.data(),
        ct.data(), pt_len,
        aad_ptr, aad_len,
        key.data(), nonce.data(),
        tag
    );
    CFX_ASSERT(rc == 0);
    CFX_ASSERT(std::memcmp(pt_out.data(), pt.data(), pt_len) == 0);

    uint64_t tsc_ticks = 0;
    for (auto _ : state) {
        unsigned aux;
        uint64_t start_cycles = __rdtscp(&aux);
        rc = cfx_chacha20_poly1305_decrypt(
            pt_out.data(),
            ct.data(), pt_len,
            aad_ptr, aad_len,
            key.data(), nonce.data(),
            tag
        );
        uint64_t end_cycles = __rdtscp(&aux);
        tsc_ticks += (end_cycles - start_cycles);
        CFX_ASSERT(rc == 0);
        benchmark::DoNotOptimize(pt_out);
    }

    const double iters = static_cast<double>(state.iterations());
    const double bytes = iters * static_cast<double>(pt_len);
    const double avg_cycles = static_cast<double>(tsc_ticks) / iters;
    const double cpb = static_cast<double>(tsc_ticks) / bytes;

    state.SetBytesProcessed(static_cast<int64_t>(bytes));
    state.SetItemsProcessed(static_cast<int64_t>(iters));
    state.counters["tsc_ticks_per_iter"] = avg_cycles;
    state.counters["tsc_ticks_per_byte"] = cpb;
}

#define BENCH_ARGS(b) \
    b->Arg(1<<6)      /*   64 B */ \
     ->Arg(1<<8)      /*  256 B */ \
     ->Arg(1<<10)     /* 1024 B */ \
     ->Arg(1<<12)     /* 4096 B */ \
     ->Arg(1<<14)     /* 16384 B */ \
     ->Arg(1<<16)     /* 65536 B */

int main(int argc, char** argv) {

    cfx_srand(0xBEBEBEBE);

    /* Encrypt, no AAD */
    {
        auto* b = benchmark::RegisterBenchmark(
            "AEAD_Encrypt_NoAAD",
            [](benchmark::State& state) {
                BM_AEAD_Encrypt(state, /*aad_len=*/0);
            });
        BENCH_ARGS(b);
    }

    /* Decrypt, no AAD */
    {
        auto* b = benchmark::RegisterBenchmark(
            "AEAD_Decrypt_NoAAD",
            [](benchmark::State& state) {
                BM_AEAD_Decrypt(state, /*aad_len=*/0);
            });
        BENCH_ARGS(b);
    }

    /* Encrypt, 32-byte AAD */
    {
        auto* b = benchmark::RegisterBenchmark(
            "AEAD_Encrypt_AAD32",
            [](benchmark::State& state) {
                BM_AEAD_Encrypt(state, /*aad_len=*/32);
            });
        BENCH_ARGS(b);
    }

    /* Decrypt, 32-byte AAD */
    {
        auto* b = benchmark::RegisterBenchmark(
            "AEAD_Decrypt_AAD32",
            [](benchmark::State& state) {
                BM_AEAD_Decrypt(state, /*aad_len=*/32);
            });
        BENCH_ARGS(b);
    }

    benchmark::Initialize(&argc, argv);
    if (benchmark::ReportUnrecognizedArguments(argc, argv)) {
        return 1;
    }
    benchmark::RunSpecifiedBenchmarks();
    return 0;
}
