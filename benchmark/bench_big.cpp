
#include <benchmark/benchmark.h>
#include <cstdint>
#include <cfx/big.h>

static void BM_InitFree(benchmark::State& state) {
    for (auto _ : state) {
        cfx_big_t b;
        cfx_big_init(&b);
        cfx_big_from_u64(&b, 42);
        cfx_big_free(&b);
        benchmark::DoNotOptimize(b);
    }
}
BENCHMARK(BM_InitFree);

static void BM_AddSmall(benchmark::State& state) {
    cfx_big_t b;
    cfx_big_init(&b);
    cfx_big_from_u64(&b, 0);
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
        cfx_big_from_u64(&b, 123456789ULL);
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
    
    cfx_big_from_u64(&b, 987654321ULL);
    for (auto _ : state) {
        cfx_big_from_u64(&a, 123456789ULL);
        cfx_big_mul(&a, &b);
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
        uint64_t rem = cfx_big_div_sm(&tmp, 7);
        benchmark::DoNotOptimize(rem);
        cfx_big_free(&tmp);
    }
    cfx_big_free(&a);
}
BENCHMARK(BM_DivSm);


static void BM_DivBig(benchmark::State& state) {
    cfx_big_t u;
    cfx_big_init(&u);
    const char* s = "1239487123691239872398127309183701928370918237019283709870708711037507120306120489389261001010192743768471263948712630489172364019823740192384701237864012837461923876501293874619283746192348751602398476";
    cfx_big_from_str(&u, s);
    
    for (auto _ : state) {
        cfx_big_t v;
        cfx_big_init(&v);
        cfx_big_from_hex(&v, "fe01023987fe08abc102870a987d87e981296300ea");
        cfx_big_t q, r;
        cfx_big_init(&q);
        cfx_big_init(&r);
        cfx_big_divrem(&q, &r, &u, &v);
        
        benchmark::DoNotOptimize(r);
        benchmark::DoNotOptimize(q);
        cfx_big_free(&q);
        cfx_big_free(&r);
        cfx_big_free(&v);

    }
    cfx_big_free(&u);
}
BENCHMARK(BM_DivBig);
// ----------------------------------------------------------------------------
// different squaring functions
cfx_big_t sq1(const cfx_big_t* b) {
    cfx_big_t ret;
    cfx_big_init(&ret);
    if (b->n == 0) return ret;
    size_t n = b->n;
    size_t szout = 2*n;
    cfx_big_reserve(&ret, szout);
    ret.n = szout;

    for (size_t i = 0; i < n; ++i) {
        uint128_t carry = 0;
        uint128_t bi = b->limb[i];
        /* terms ij == terms ji -> iterate over half the terms but double them*/
        for (size_t j = i+1; j < n; ++j) {
            uint128_t prod = bi * b->limb[j];
            prod <<= 1; /* double */
            prod |= carry;
            prod += ret.limb[i+j];
            ret.limb[i+j] = (uint64_t)prod;
            carry = prod >> 64;
            // printf("doubling term i: %zu, j: %zu; prod: %llu, carry: %llu\n", i, j, prod, (uint64_t)carry);
        }
        uint128_t sq = bi*bi;
        sq += ret.limb[2*i];
        uint64_t lo = (uint64_t)sq;
        uint128_t c2 = sq >> 64;
        // printf("squaring term i: %zu, lo: %llu, carry: %llu\n", i, lo, (uint64_t)carry);

        // propagate carry from cross terms into next limb
        uint128_t u = (uint128_t)ret.limb[i + n] + carry + c2;
        ret.limb[i + n] = (uint64_t)u;
        // store fixed lower
        ret.limb[2 * i] = lo;
        // carry continues if u overflowed 64 bits
        uint128_t k = u >> 64;
        if (k) {
            size_t idx = i + n + 1;
            while (k) {
                printf("carry continues %llu\n", (uint64_t)k);
                if (idx >= szout) break;
                uint128_t w = (uint128_t)ret.limb[idx] + (uint64_t)k;
                ret.limb[idx] = (uint64_t)w;
                k = w >> 64;
                ++idx;
            }
        }
    }
    // _trim_leading_zeros(&ret);
    return ret;
}
cfx_big_t sq2(const cfx_big_t* b) {
    cfx_big_t ret;
    cfx_big_init(&ret);
    if (b->n == 0) return ret;

    size_t n = b->n;
    size_t out_n = 2 * n;
    // temp buffer initialized to zero
    uint64_t* tmp = (uint64_t*)calloc(out_n, sizeof(uint64_t));
    if (!tmp) { /* OOM handling as per your policy */ return ret; }

    for (size_t i = 0; i < n; ++i) {
        uint128_t carry = 0;

        // cross terms j > i counted twice
        size_t j = i + 1;
        for (; j < n; ++j) {
            uint128_t t = (uint128_t)b->limb[i] * (uint128_t)b->limb[j];
            t <<= 1; // times 2
            t += tmp[i + j];
            t += carry;
            tmp[i + j] = (uint64_t)t;
            carry = t >> 64;
        }

        // add the square term (i,i)
        uint128_t t2 = (uint128_t)b->limb[i] * (uint128_t)b->limb[i];
        t2 += tmp[2 * i];
        uint64_t lo = (uint64_t)t2;
        uint128_t c2 = t2 >> 64;

        // propagate carry from cross terms into next limb
        uint128_t u = (uint128_t)tmp[i + n] + carry + c2;
        tmp[i + n] = (uint64_t)u;
        // store fixed lower
        tmp[2 * i] = lo;
        // carry continues if u overflowed 64 bits
        uint128_t k = u >> 64;
        if (k) {
            size_t idx = i + n + 1;
            while (k) {
                if (idx >= out_n) break;
                uint128_t w = (uint128_t)tmp[idx] + (uint64_t)k;
                tmp[idx] = (uint64_t)w;
                k = w >> 64;
                ++idx;
            }
        }
    }

    // move tmp into ret
    cfx_big_reserve(&ret, out_n);
    memcpy(ret.limb, tmp, out_n * sizeof(uint64_t));
    ret.n = out_n;
    // _trim_leading_zeros(ret);
    free(tmp);
    return ret;
}

// -------------
cfx_big_t sq3(const cfx_big_t* b) {
    const size_t n = b->n;
    cfx_big_t ret;
    cfx_big_init(&ret);

    if (n == 0) {
        return ret;
    }
    
    size_t cnt = 0;
    
    cfx_big_reserve(&ret, 2*n);
    memset(ret.limb, 0, 2*n * sizeof(uint64_t));
    ret.n = 2*n;

    // 1) Cross terms: for i < j, add 2*b[i]*b[j] into ret[i+j]
    for (size_t i = 0; i < n; ++i) {
        uint128_t carry = 0;
        for (size_t j = i + 1; j < n; ++j) {
            uint128_t p = (uint128_t)b->limb[i] * b->limb[j];

            // add p once
            uint128_t t = (uint128_t)ret.limb[i + j]
                            + (uint64_t)p
                            + (uint64_t)carry;
            ret.limb[i + j] = (uint64_t)t;
            carry = (carry >> 64) + (t >> 64) + (p >> 64);

            // add p again (to double) -- same carry rule
            t = (uint128_t)ret.limb[i + j]
                + (uint64_t)p
                + (uint64_t)carry;
            ret.limb[i + j] = (uint64_t)t;
            carry = (carry >> 64) + (t >> 64) + (p >> 64);
            cnt++;
        }

        // propagate whatever is left in 'carry'
        size_t k = i + n; // next column after the last updated (i + (n-1))
        while (carry) {
            uint128_t t = (uint128_t)ret.limb[k] + (uint64_t)carry;
            ret.limb[k] = (uint64_t)t;
            carry = (carry >> 64) + (t >> 64);
            ++k;
            cnt++;
        }
    }

    // 2) Diagonals: add b[i]^2 once at ret[2*i]
    for (size_t i = 0; i < n; ++i) {
        uint128_t sq = (uint128_t)b->limb[i] * b->limb[i];

        uint128_t t = (uint128_t)ret.limb[2*i] + (uint64_t)sq;
        ret.limb[2*i] = (uint64_t)t;
        uint128_t c = (t >> 64) + (sq >> 64);

        size_t k = 2*i + 1;
        cnt++;
        while (c) {
            t = (uint128_t)ret.limb[k] + (uint64_t)c;
            ret.limb[k] = (uint64_t)t;
            c = (c >> 64) + (t >> 64);
            ++k;
        }
    }

    // trim
    while (ret.n && ret.limb[ret.n - 1] == 0) --ret.n;

    // printf("%s loop cnt: %zu\n", __func__, cnt);
    return ret;
}


static void BM_sq1(benchmark::State& state) {
    cfx_big_t a, b;
    cfx_big_init(&a);
    cfx_big_init(&b);

    for (auto _ : state) {
        cfx_big_from_u64(&a, 123456789ULL);
        b = sq1(&a);
        benchmark::DoNotOptimize(a);
        benchmark::DoNotOptimize(b);
    }
    cfx_big_free(&a);
    cfx_big_free(&b);
}
BENCHMARK(BM_sq1);

static void BM_sq2(benchmark::State& state) {
    cfx_big_t a, b;
    cfx_big_init(&a);
    cfx_big_init(&b);

    for (auto _ : state) {
        cfx_big_from_u64(&a, 123456789ULL);
        b = sq1(&a);
        benchmark::DoNotOptimize(a);
        benchmark::DoNotOptimize(b);
    }
    cfx_big_free(&a);
    cfx_big_free(&b);
}
BENCHMARK(BM_sq2);

static void BM_sq3(benchmark::State& state) {
    cfx_big_t a, b;
    cfx_big_init(&a);
    cfx_big_init(&b);

    for (auto _ : state) {
        cfx_big_from_u64(&a, 123456789ULL);
        b = sq1(&a);
        benchmark::DoNotOptimize(a);
        benchmark::DoNotOptimize(b);
    }
    cfx_big_free(&a);
    cfx_big_free(&b);
}
BENCHMARK(BM_sq3);

BENCHMARK_MAIN();
