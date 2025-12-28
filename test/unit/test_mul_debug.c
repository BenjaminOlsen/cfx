/*
 * Debug test for multiplication issue
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>

#include "cfx/big.h"
#include "cfx/arith.h"

/* Local debug multiplication to trace what's happening */
static void debug_mul(cfx_big_t* out, const cfx_big_t* a, const cfx_big_t* b) {
    printf("\n=== DEBUG MUL ===\n");
    printf("a: n=%zu, cap=%zu, limb=%p\n", a->n, a->cap, (void*)a->limb);
    printf("b: n=%zu, cap=%zu, limb=%p\n", b->n, b->cap, (void*)b->limb);
    printf("out: n=%zu, cap=%zu, limb=%p\n", out->n, out->cap, (void*)out->limb);

    for (size_t i = 0; i < a->n; i++)
        printf("  a[%zu] = 0x%08" PRIx32 "\n", i, (uint32_t)a->limb[i]);
    for (size_t i = 0; i < b->n; i++)
        printf("  b[%zu] = 0x%08" PRIx32 "\n", i, (uint32_t)b->limb[i]);

    size_t na = a->n;
    size_t nb = b->n;

    cfx_big_reserve(out, na + nb);
    printf("After reserve: out.cap=%zu, out.limb=%p\n", out->cap, (void*)out->limb);

    out->n = na + nb;
    memset(out->limb, 0, out->n * sizeof(cfx_limb_t));

    printf("\nStarting schoolbook multiplication...\n");

    for (size_t i = 0; i < nb; ++i) {
        cfx_acc_t carry = 0;
        cfx_acc_t bi = (cfx_acc_t)b->limb[i];
        printf("\nOuter i=%zu, bi=0x%" PRIx64 "\n", i, (uint64_t)bi);

        size_t k = i;
        for (size_t j = 0; j < na; ++j, ++k) {
            cfx_limb_t aj_raw = a->limb[j];
            cfx_acc_t aj = (cfx_acc_t)aj_raw;
            cfx_acc_t outk = (cfx_acc_t)out->limb[k];
            cfx_acc_t prod = bi * aj;
            cfx_acc_t s = outk + prod + carry;

            printf("  j=%zu k=%zu:\n", j, k);
            printf("    a[j]_raw = 0x%08" PRIx32 " (%u)\n", (uint32_t)aj_raw, (unsigned)aj_raw);
            printf("    aj (cast to acc) = 0x%" PRIx64 "\n", (uint64_t)aj);
            printf("    bi = 0x%" PRIx64 "\n", (uint64_t)bi);
            printf("    prod = bi * aj = 0x%" PRIx64 "\n", (uint64_t)prod);
            printf("    outk = 0x%" PRIx64 ", carry = 0x%" PRIx64 "\n", (uint64_t)outk, (uint64_t)carry);
            printf("    s = outk + prod + carry = 0x%" PRIx64 "\n", (uint64_t)s);

            out->limb[k] = (cfx_limb_t)s;
            carry = s >> CFX_LIMB_BITS;
            printf("    -> out[%zu]=0x%08" PRIx32 ", new carry=0x%" PRIx64 "\n", k, (uint32_t)out->limb[k], (uint64_t)carry);
        }

        while (carry) {
            cfx_acc_t s = (cfx_acc_t)out->limb[k] + carry;
            printf("  carry propagate: k=%zu, out[k]=0x%" PRIx64 " + carry=0x%" PRIx64 " = s=0x%" PRIx64 "\n",
                   k, (uint64_t)out->limb[k], (uint64_t)carry, (uint64_t)s);
            out->limb[k] = (cfx_limb_t)s;
            carry = s >> CFX_LIMB_BITS;
            printf("     -> out[%zu]=0x%08" PRIx32 ", carry=0x%" PRIx64 "\n", k, (uint32_t)out->limb[k], (uint64_t)carry);
            ++k;
        }

        printf("  After i=%zu: out = {", i);
        for (size_t x = 0; x < out->n; x++) {
            printf("0x%08" PRIx32 "%s", (uint32_t)out->limb[x], x+1<out->n ? ", " : "");
        }
        printf("}\n");
    }

    /* Trim */
    while (out->n && out->limb[out->n-1] == 0) --out->n;
    printf("\nFinal (trimmed) out: n=%zu, limbs = {", out->n);
    for (size_t x = 0; x < out->n; x++) {
        printf("0x%08" PRIx32 "%s", (uint32_t)out->limb[x], x+1<out->n ? ", " : "");
    }
    printf("}\n");
}

int main(void) {
    printf("=== Multiplication Debug Test ===\n");
    printf("CFX_LIMB_BITS = %d\n", CFX_LIMB_BITS);
    printf("sizeof(cfx_limb_t) = %zu\n", sizeof(cfx_limb_t));
    printf("sizeof(cfx_acc_t) = %zu\n\n", sizeof(cfx_acc_t));

    /* Test raw 64-bit multiplication */
    printf("--- Testing raw 64-bit multiplication ---\n");
    {
        cfx_limb_t a_raw = 0xFFFFFFC5u;
        cfx_limb_t b_raw = 0xFFFFFFADu;

        printf("a_raw = 0x%08" PRIx32 " (%u)\n", a_raw, a_raw);
        printf("b_raw = 0x%08" PRIx32 " (%u)\n", b_raw, b_raw);

        cfx_acc_t a_64 = (cfx_acc_t)a_raw;
        cfx_acc_t b_64 = (cfx_acc_t)b_raw;

        printf("a_64 = 0x%" PRIx64 "\n", (uint64_t)a_64);
        printf("b_64 = 0x%" PRIx64 "\n", (uint64_t)b_64);

        cfx_acc_t prod = a_64 * b_64;
        printf("prod = a_64 * b_64 = 0x%" PRIx64 "\n", (uint64_t)prod);
        printf("Expected: 0xffffff7200001321\n");

        /* Also try with explicit uint64_t */
        uint64_t aa = (uint64_t)a_raw;
        uint64_t bb = (uint64_t)b_raw;
        uint64_t pp = aa * bb;
        printf("uint64_t version: pp = 0x%" PRIx64 "\n", pp);

        /* Also try direct - no casting */
        uint64_t direct = (uint64_t)0xFFFFFFC5u * (uint64_t)0xFFFFFFADu;
        printf("Direct literal: 0x%" PRIx64 "\n", direct);

        /* Try with ULL suffix */
        uint64_t ull = 0xFFFFFFC5ull * 0xFFFFFFADull;
        printf("With ULL: 0x%" PRIx64 "\n", ull);
    }
    printf("\n");

    cfx_big_t p, q, n;
    cfx_big_init(&p);
    cfx_big_init(&q);
    cfx_big_init(&n);

    /* Parse p */
    cfx_big_from_str(&p, "18446744073709551557");
    printf("p = 18446744073709551557 (hex: 0xFFFFFFFFFFFFFFC5)\n");
    printf("  p.n = %zu limbs\n", p.n);
    for (size_t i = 0; i < p.n; i++) {
        printf("  p.limb[%zu] = 0x%08" PRIx32 "\n", i, (uint32_t)p.limb[i]);
    }

    /* Verify p limbs are correct */
    int p_ok = (p.n == 2 && p.limb[0] == 0xFFFFFFC5u && p.limb[1] == 0xFFFFFFFFu);
    printf("  p limbs correct: %s\n", p_ok ? "YES" : "NO");

    /* Parse q */
    cfx_big_from_str(&q, "18446744073709551533");
    printf("\nq = 18446744073709551533 (hex: 0xFFFFFFFFFFFFFFAD)\n");
    printf("  q.n = %zu limbs\n", q.n);
    for (size_t i = 0; i < q.n; i++) {
        printf("  q.limb[%zu] = 0x%08" PRIx32 "\n", i, (uint32_t)q.limb[i]);
    }

    /* Verify q limbs are correct */
    int q_ok = (q.n == 2 && q.limb[0] == 0xFFFFFFADu && q.limb[1] == 0xFFFFFFFFu);
    printf("  q limbs correct: %s\n", q_ok ? "YES" : "NO");

    /* Simple test: 3 * 5 = 15 */
    printf("\n--- Simple test: 3 * 5 ---\n");
    cfx_big_t a, b, c;
    cfx_big_init(&a);
    cfx_big_init(&b);
    cfx_big_init(&c);
    cfx_big_from_limb(&a, 3);
    cfx_big_from_limb(&b, 5);
    cfx_big_mul(&c, &a, &b);
    printf("3 * 5 = %u (expected 15)\n", (unsigned)c.limb[0]);
    cfx_big_free(&a);
    cfx_big_free(&b);
    cfx_big_free(&c);

    /* Test: 0xFFFFFFFF * 0xFFFFFFFF */
    printf("\n--- Test: 0xFFFFFFFF * 0xFFFFFFFF ---\n");
    cfx_big_init(&a);
    cfx_big_init(&b);
    cfx_big_init(&c);
    cfx_big_from_limb(&a, 0xFFFFFFFFu);
    cfx_big_from_limb(&b, 0xFFFFFFFFu);
    cfx_big_mul(&c, &a, &b);
    printf("Result: n=%zu, limb[0]=0x%08" PRIx32 ", limb[1]=0x%08" PRIx32 "\n",
           c.n, (uint32_t)c.limb[0], c.n > 1 ? (uint32_t)c.limb[1] : 0);
    printf("Expected: 0xFFFFFFFE00000001 -> limb[0]=0x00000001, limb[1]=0xFFFFFFFE\n");
    cfx_big_free(&a);
    cfx_big_free(&b);
    cfx_big_free(&c);

    /* Now multiply p * q using debug version */
    printf("\n--- Main test: p * q ---\n");
    debug_mul(&n, &p, &q);

    printf("n = p * q:\n");
    printf("  n.n = %zu limbs\n", n.n);
    for (size_t i = 0; i < n.n && i < 8; i++) {
        printf("  n.limb[%zu] = 0x%08" PRIx32 "\n", i, (uint32_t)n.limb[i]);
    }

    /* Expected (from Python):
     * p * q = 340282366920938460843936948965011886881
     * hex(p*q) = 0xffffffffffffff720000000000001321
     *
     * In 32-bit little-endian limbs:
     *   limb[0] = 0x00001321
     *   limb[1] = 0x00000000
     *   limb[2] = 0xffffff72
     *   limb[3] = 0xffffffff
     */
    printf("\nExpected limbs (from Python hex(p*q) = 0xffffffffffffff720000000000001321):\n");
    printf("  limb[0] = 0x00001321\n");
    printf("  limb[1] = 0x00000000\n");
    printf("  limb[2] = 0xffffff72\n");
    printf("  limb[3] = 0xffffffff\n");

    int n_ok = (n.n == 4 &&
                n.limb[0] == 0x00001321u &&
                n.limb[1] == 0x00000000u &&
                n.limb[2] == 0xffffff72u &&
                n.limb[3] == 0xffffffffu);
    printf("\nn limbs correct: %s\n", n_ok ? "YES" : "NO");

    char* str = cfx_big_to_str(&n, NULL);
    printf("\nn as decimal: %s\n", str);
    printf("Expected:     340282366920938460843936948965011886881\n");
    free(str);

    cfx_big_free(&p);
    cfx_big_free(&q);
    cfx_big_free(&n);

    printf("\nDone.\n");
    return n_ok ? 0 : 1;
}
