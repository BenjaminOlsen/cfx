/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */
/* test_big_mul.c - Tests for big integer multiplication */

#include "test_common.h"

static void test_mul1(void) {
    cfx_big_t b, m;
    cfx_big_init(&b);
    cfx_big_from_limb(&b, 5);
    cfx_big_init(&m);
    cfx_big_from_limb(&m, 7);

    cfx_big_mul_eq(&b, &m);
    cfx_limb_t expect[] = {35};
    big_expect_limbs(__func__, &b, expect, 1);

    cfx_big_free(&b);
    cfx_big_free(&m);
    PRINT_TEST(1);
}

static void test_mul_adduiv(void) {
    cfx_big_t a, m;
    cfx_big_init(&a);
    cfx_big_from_limb(&a, 0);
    cfx_big_init(&m);
    cfx_big_from_limb(&m, 1);
    cfx_limb_t K = 0x1B2CDE;
    cfx_big_mul_sm(&m, K);

    for (cfx_limb_t i = 0; i < K; ++i) {
        cfx_big_add_sm(&a, 1);
    }
    CFX_ASSERT(cfx_big_eq(&a, &m));

    cfx_big_free(&a);
    cfx_big_free(&m);
    PRINT_TEST(1);
}

static void test_carry_two_limbs_times_2(void) {
    cfx_big_t b, m;
    cfx_big_init(&b);
#if CFX_LIMB_BITS == 64
    cfx_limb_t limbs_b[] = {UINT64_C(0xFFFFFFFFFFFFFFFF), UINT64_C(0xFFFFFFFFFFFFFFFF)};
#elif CFX_LIMB_BITS == 32
    cfx_limb_t limbs_b[] = {UINT32_C(0xFFFFFFFF), UINT32_C(0xFFFFFFFF)};
#endif
    cfx_big_from_limbs(&b, limbs_b, 2);
    cfx_big_init(&m);
    cfx_big_from_limb(&m, 2);

    cfx_big_mul_eq(&b, &m);
#if CFX_LIMB_BITS == 64
    cfx_limb_t expect[] = {UINT64_C(0xFFFFFFFFFFFFFFFE), UINT64_C(0xFFFFFFFFFFFFFFFF), UINT64_C(1)};
#elif CFX_LIMB_BITS == 32
    cfx_limb_t expect[] = {UINT32_C(0xFFFFFFFE), UINT32_C(0xFFFFFFFF), UINT32_C(1)};
#endif
    big_expect_limbs(__func__, &b, expect, 3);

    cfx_big_free(&b);
    cfx_big_free(&m);
    PRINT_TEST(1);
}

static void test_mul_by_base_2_64_shift(void) {
    cfx_big_t b, m;
    cfx_big_init(&b);
    cfx_big_init(&m);
#if CFX_LIMB_BITS == 64
    cfx_limb_t limbs_b[] = {UINT64_C(0x0123456789ABCDEF), UINT64_C(0x0FEDCBA987654321), UINT64_C(0x0000000000000001)};
#elif CFX_LIMB_BITS == 32
    cfx_limb_t limbs_b[] = {UINT32_C(0x01234567), UINT32_C(0x0FEDCBA9), UINT32_C(0x00000001)};
#endif
    size_t sz0 = sizeof(limbs_b)/sizeof(limbs_b[0]);
    cfx_big_from_limbs(&b, limbs_b, sz0);
#if CFX_LIMB_BITS == 64
    cfx_limb_t limbs_m[] = {UINT64_C(0), UINT64_C(1)};
#elif CFX_LIMB_BITS == 32
    cfx_limb_t limbs_m[] = {UINT32_C(0), UINT32_C(1)};
#endif
    cfx_big_from_limbs(&m, limbs_m, 2);

    const size_t N = 100;  /* reduced from 1111 for faster tests */

    for (size_t n = 1; n < N; ++n) {
        cfx_big_mul_eq(&b, &m);
        cfx_limb_t* expect = (cfx_limb_t*)malloc((n + sz0)*sizeof(cfx_limb_t));
        for (size_t k = 0; k < n; ++k) {
            expect[k] = 0;
        }
#if CFX_LIMB_BITS == 64
        expect[n+sz0-3] = UINT64_C(0x0123456789ABCDEF);
        expect[n+sz0-2] = UINT64_C(0x0FEDCBA987654321);
        expect[n+sz0-1] = UINT64_C(0x0000000000000001);
#elif CFX_LIMB_BITS == 32
        expect[n+sz0-3] = UINT32_C(0x01234567);
        expect[n+sz0-2] = UINT32_C(0x0FEDCBA9);
        expect[n+sz0-1] = UINT32_C(0x00000001);
#endif

        big_expect_limbs(__func__, &b, expect, sz0 + n);
        free(expect);
    }

    cfx_big_free(&b);
    cfx_big_free(&m);
    PRINT_TEST(1);
}

static void test_self_multiply_square(void) {
    cfx_big_t b;
    cfx_big_init(&b);
#if CFX_LIMB_BITS == 64
    cfx_limb_t limbs[] = {UINT64_C(0xFFFFFFFFFFFFFFFF)};
#elif CFX_LIMB_BITS == 32
    cfx_limb_t limbs[] = {UINT32_C(0xFFFFFFFF)};
#endif

    cfx_big_from_limbs(&b, limbs, 1);
    big_expect_limbs(__func__, &b, limbs, 1);
    cfx_big_mul_eq(&b, &b);
#if CFX_LIMB_BITS == 64
    cfx_limb_t expect[] = {UINT64_C(1), UINT64_C(0xFFFFFFFFFFFFFFFE)};
#elif CFX_LIMB_BITS == 32
    cfx_limb_t expect[] = {UINT32_C(1), UINT32_C(0xFFFFFFFE)};
#endif

    big_expect_limbs(__func__, &b, expect, 2);
    cfx_big_free(&b);
    PRINT_TEST(1);
}

static void test_self_multiply_big(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    size_t N = 10;
    cfx_limb_t* limbs = (cfx_limb_t*)malloc(N * sizeof(cfx_limb_t));
    for (size_t i = 0; i < N; ++i) {
        limbs[i] = 2;
    }
    cfx_big_from_limbs(&b, limbs, N);
    cfx_big_mul_eq(&b, &b);
    cfx_big_free(&b);
    free(limbs);
    PRINT_TEST(1);
}

static void test_known_squares(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    cfx_big_from_dec(&b,
        "12554203470773361528352143580257209"
        "759168353591939024551938");
    cfx_big_sq(&b);
    char* s = cfx_big_to_str(&b, NULL);
    char* expect =
        "15760802478557791686620405668794173"
        "58808686358020659979337980952437065"
        "51032029871143396883518477623727858"
        "361659555844";
    CFX_ASSERT(strcmp(s, expect) == 0);
    free(s);

    cfx_big_from_dec(&b,
        "14536774485912137811774516281687980"
        "27136112775646765168338161504493023"
        "20618275753867907499968765052767290"
        "9553511701551770336148695547906");
    cfx_big_sq(&b);
    s = cfx_big_to_str(&b, NULL);

    expect =
        "211317812454266098563847037101386611"
        "572936979272149330366075109813640724"
        "474345860994392765176988025566104455"
        "837520931631539188377309778866973982"
        "401217516000816262383671909413419933"
        "179663682199874639617901251522391335"
        "258229693927038827878426679857709370"
        "5798799065540984836";
    CFX_ASSERT(strcmp(s, expect) == 0);
    free(s);

    cfx_big_from_dec(&b,
        "788040123927889584288300542721290477751811"
        "061317446299559616920266928982601233965513"
        "84510299169195903266945438318594");
    cfx_big_sq(&b);
    s = cfx_big_to_str(&b, NULL);

    expect =
        "621007236920283574126921534924820203000914"
        "006088690925637430299255754323320524243949"
        "820727721413359346225814236510862646127210"
        "854935998980345500359712802179223868607115"
        "902506409479615101839926583004283044656249"
        "4079870273849846136836";
    CFX_ASSERT(strcmp(s, expect) == 0);
    free(s);

    cfx_big_free(&b);
    PRINT_TEST(1);
}

static void test_known_squares_2(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    cfx_big_from_dec(&b,
        "4946608029462090681478206578991795742708644"
        "2564742658586426229545514803499564697000372"
        "3095350971345437292114654548843072761868784"
        "674125049315509629339381027496416991005114370");
    cfx_big_mul_eq_csa(&b, &b);
    char* s = cfx_big_to_str(&b, NULL);

    char* expect =
        "24468930997138827791465884302231872460739875"
        "68207399045598121707346936565872814004700782"
        "89425671088010093095401304027496760620656589"
        "55442387983480255066626620918677531232865400"
        "98771942569641370987258548756784517935582567"
        "13069390644656846533223333372022161629347849"
        "16188372769675748466345955863374240885700087"
        "1949664098243244687781332547496780496900";
    CFX_ASSERT(strcmp(s, expect) == 0);
    free(s);

    cfx_big_mul_eq_csa(&b, &b);
    s = cfx_big_to_str(&b, NULL);

    expect =
        "598728584142741349308708530097477655730195121"
        "8921561456478992373941874807025827989117005224"
        "4544159576837861157818585384355732464079298824"
        "8840307403854743209208295687980560227400969404"
        "2099334329107992954561531802255073931559398056"
        "4137316892818829069898815609946436650172187960"
        "4903063742870707747793924033318816999047826082"
        "6483886794388751538092273101337790920556525050"
        "8202721113875166122496565273190156790054685080"
        "2591709257787294752819636451016941276946114924"
        "0557203469761158225922311367092258610357178141"
        "4512816029237243505436725602625497437379829908"
        "1001420696314530701776191559963115206965790843"
        "4600574609622913370900822575527296202868593414"
        "9525207968313075153130098691732396070700210909"
        "610000";
    CFX_ASSERT(strcmp(s, expect) == 0);
    free(s);

    cfx_big_free(&b);
    PRINT_TEST(1);
}

static void test_sq_zero(void) {
    cfx_big_t a;
    cfx_big_init(&a);
    cfx_big_sq(&a);
    CFX_ASSERT(cfx_big_is_zero(&a));
    cfx_big_free(&a);
}

static void test_sq_one(void) {
    cfx_big_t a;
    cfx_big_init(&a);
    cfx_big_from_limb(&a, 1);
    cfx_big_sq(&a);
    CFX_ASSERT(cfx_big_eq_u64(&a, 1));
    cfx_big_free(&a);
}

static void test_sq_small(void) {
    cfx_big_t a;
    cfx_big_init(&a);
    cfx_big_from_limb(&a, 12);
    cfx_big_sq(&a);
    CFX_ASSERT(cfx_big_eq_u64(&a, 144));
    cfx_big_free(&a);
}

int main(void) {
    CFX_TEST(test_mul1);
    CFX_TEST(test_mul_adduiv);
    CFX_TEST(test_carry_two_limbs_times_2);
    CFX_TEST(test_mul_by_base_2_64_shift);
    CFX_TEST(test_self_multiply_square);
    CFX_TEST(test_self_multiply_big);
    CFX_TEST(test_known_squares);
    CFX_TEST(test_known_squares_2);
    CFX_TEST(test_sq_zero);
    CFX_TEST(test_sq_one);
    CFX_TEST(test_sq_small);
    puts("OK");
    return 0;
}
