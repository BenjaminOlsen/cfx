#include "cfx/big.h"
#include "cfx/macros.h"
#include "cfx/error.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define PRINT_TEST(ok) \
    do { \
        printf("%s() ---- %s\n", __func__, ok ? "ok" : "NOT OK"); \
    } while (0)

static void test_cfx_big_init(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    assert(b.limb == NULL);
    assert(b.n == 0);
    assert(b.cap == 0);
    PRINT_TEST(1);
}

static void test_cfx_big_reserve(void) {
    cfx_big_t b;
    size_t rcap1 = 55;
    cfx_big_reserve(&b, rcap1);
    assert(b.cap == rcap1);
    assert(b.n == 0);
    assert(b.limb != NULL);
    size_t rcap2 = rcap1 / 2;
    cfx_big_reserve(&b, rcap2);  /* shouldn't reserve less space. */
    assert(b.cap == rcap1);
    PRINT_TEST(1);
}

static void test_copy_swap(void) {
    cfx_big_t a, b;
    uint64_t la[] = {1, 2, 3, 4};
    uint64_t lb[] = {UINT64_MAX-1, UINT64_MAX-2, UINT64_MAX-3, UINT64_MAX-4};

    cfx_big_from_limbs(&a, la, 4);
    cfx_big_from_limbs(&b, lb, 4);
    assert(a.n == b.n);

    for (size_t i = 0; i < sizeof(la)/sizeof(uint64_t); ++i) {
        assert(a.limb[i] == la[i]);
        assert(b.limb[i] == lb[i]);
    }
    cfx_big_t aa, bb;
    cfx_big_copy(&aa, &a);
    cfx_big_copy(&bb, &b);

    assert(cfx_big_eq(&aa, &a));
    assert(cfx_big_eq(&bb, &b));

    cfx_big_swap(&a, &b);

    assert(cfx_big_eq(&aa, &b));
    assert(cfx_big_eq(&bb, &a));
    PRINT_TEST(1);
}

static void big_init_from_limbs_base_1e9(cfx_big_t* b, const uint64_t *limbs, size_t n) {
    cfx_big_init(b);
    if (n == 0) return;
    for (ssize_t i = n - 1; i >= 0; --i) {
        cfx_big_mul_sm(b, 1000000000ull);
        cfx_big_add_sm(b, limbs[i]);
    }
    PRINT_TEST(1);
}

static void test_mul_by_zero(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    cfx_big_set_val(&b, 123);
    cfx_big_mul_sm(&b, 2838);
    cfx_big_mul_sm(&b, 1928);
    cfx_big_mul_sm(&b, 9);
    cfx_big_mul_sm(&b, 123765);
    cfx_big_mul_sm(&b, 0);
    size_t sz = 0;
    char* s = cfx_big_to_str(&b, &sz);
    CFX_PRINT_DBG("s: %s, sz: %zu\n", s, sz);
    assert(sz == 1);
    assert(strcmp(s, "0") == 0);
    free(s);
    PRINT_TEST(1);
}

static void test_add_sm(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    CFX_BIG_PRINTF(&b, "init: b.n:%ld; ", b.n);
    
    cfx_big_set_val(&b, 123);
    CFX_BIG_PRINTF(&b, "after setting val: ");

    cfx_big_add_sm(&b, 321);
    CFX_BIG_PRINTF(&b, "after add: ");
    assert(b.limb[0] == 444);
    cfx_big_set_val(&b, UINT64_MAX);
    CFX_BIG_PRINTF(&b, "after set:");
    
    assert(b.limb[0] == UINT64_MAX);
    cfx_big_add_sm(&b, 1);
    CFX_BIG_PRINTF(&b, "after carry: ");
    assert(b.limb[0] == 0);
    assert(b.limb[1] == 1);
    PRINT_TEST(1);
}


static void test_sub_sm(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    CFX_BIG_PRINTF(&b, "init: b.n:%ld; ", b.n);

    cfx_big_set_val(&b, 1);
    cfx_big_sub_sm(&b, 1);
    assert(b.limb[0] == 0);
    cfx_big_sub_sm(&b, 1);
    assert(b.limb[0] == 0);
    
    cfx_big_set_val(&b, UINT64_MAX);
    CFX_BIG_PRINTF(&b, "set to UINT64_MAX: ");
    assert(b.limb[0] == UINT64_MAX);
    cfx_big_add_sm(&b, 1);
    CFX_BIG_PRINTF(&b, "add 1: ");
    assert(b.limb[0] == 0);
    assert(b.limb[1] == 1);

    cfx_big_sub_sm(&b, 1);
    CFX_BIG_PRINTF(&b, "sub 1: ");
    assert(b.limb[0] == UINT64_MAX);
    assert(b.n == 1);

    const uint64_t N = 12;
    for (uint64_t n = 0; n < N; ++n) {
        cfx_big_sub_sm(&b, 1);
        CFX_BIG_PRINTF(&b, "sub 1: ");
    }
    assert(b.limb[0] == UINT64_MAX - N);

    uint64_t q = 100;
    uint64_t orig = b.limb[0];
    for (uint64_t n = 0; n < 2; ++n) {
        uint128_t s = (uint128_t)b.limb[0] + q;
        cfx_big_add_sm(&b, q);
        CFX_BIG_PRINTF(&b, "add %llu: ", q);
        assert(b.limb[0] == (uint64_t)s);
        assert(b.n > 1);
        assert(b.limb[1] == (uint64_t)(s >> 64));
        cfx_big_sub_sm(&b, q);
        CFX_BIG_PRINTF(&b, "sub %llu: ", q);
        assert(b.limb[0] == orig);
    }
    PRINT_TEST(1);
}

/* helper to run one limb test */
static void check(const char *label, const uint64_t *limbs, size_t n, const char *expect) {
    cfx_big_t b;
    big_init_from_limbs_base_1e9(&b, limbs, n);
    size_t len = 0;
    char *s = cfx_big_to_str(&b, &len);
    // CFX_BIG_PRINTF(&b, "str is:\n");
    // CFX_PRINT_DBG("str should be:\n%s\n", expect);
    assert(strcmp(s, expect) == 0);
    assert(len == strlen(expect));
    assert(s[len] == '\0');
    CFX_PRINT_DBG("[ok] %s -> %s\n", label, s);
    free(s);
    cfx_big_free(&b);
}


// Single small limb
static void test_limb1(void) {
    uint64_t L[] = { 123456789u };
    check("single limb", L, 1, "123456789");
    PRINT_TEST(1);
}

// -->) Single 0 limb but n>0 (should normally be normalized away; if you allow it, behavior is “000…000”?)
// Prefer: represent zero as n==0. You can skip this if you enforce normalization.

// Two limbs with inner zero-padding needed:
// value = limb[1]*1e9 + limb[0] = 42*1e9 + 123456789
static void test_limb2(void) {
    uint64_t L[] = { 123456789u, 4200u };
    check("two limbs pad", L, 2, "4200123456789");
    PRINT_TEST(1);
}

// Inner limb exact zero-padding boundary: limb[0] has fewer than 9 digits
static void test_limb3(void) {
    uint64_t L[] = { 1u, 1u };
    check("two limbs tiny low", L, 2, "1000000001");
    PRINT_TEST(1);
}

// Max limb values
static void test_limb4(void) {
    uint64_t L[] = { 999999999u, 999999999u };
    check("two limbs max", L, 2, "999999999999999999");
    PRINT_TEST(1);
}

// Four limbs mixed
// value = 1*1e27 + 7*1e18 + 42*1e9 + 5 → "1 000000007 000000042 000000005"
static void test_limb5(void) {
    uint64_t L[] = { 5u, 42u, 7u, 1u };
    check("four limbs pad", L, 4, "1000000007000000042000000005");
    PRINT_TEST(1);
}

static void test_limb6(void) {
    uint64_t L[] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    check("9 limbs pad", L, 9,
        "900000000"
        "800000000"
        "700000000"
        "600000000"
        "500000000"
        "400000000"
        "300000000"
        "200000000"
        "1"
    );
    PRINT_TEST(1);
}

// Large ndigits sanity: build via mul to exercise carry
static void test_limb7(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    cfx_big_set_val(&b, 1);
    for (int i = 0; i < 10; ++i) cfx_big_mul_sm(&b, 1000000000u - 1u); // (1e9-1)^10
    char *s = cfx_big_to_str(&b, NULL);
    // spot checks: starts with '9' and length >= 9
    CFX_PRINT_DBG("%s\n", s);
    assert(s[0] == '9');
    assert(strlen(s) >= 9);
    char* expect = "999999990000000044999999880000000209999999"
        "748000000209999999880000000044999999990000000001";
    assert(strcmp(s, expect) == 0);
    free(s);
    cfx_big_free(&b);
    PRINT_TEST(1);
}

static void test_str1(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    const char *sin =   "99911231231238761239876981273469128374169283476129"
                        "38471629384761250389172603459812630498672312387123"
                        "87123981723918273912891238719248719238719248169723"
                        "00091203901290909090911100091231283761000101023882";
    cfx_big_from_str(&b, sin);
    char *sout = cfx_big_to_str(&b, NULL);
    int ok = (strcmp(sin, sout) == 0); 
    CFX_PRINT_DBG("test str1: in:\n%s \n%s\nout.. %s\n", sin, sout,
        ok ? "ok":"NOT ok");
    assert(ok);
    PRINT_TEST(1);
}

static void test_str2(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    const char *sin = "9218";
    cfx_big_from_str(&b, sin);
    char *sout = cfx_big_to_str(&b, NULL);
    int ok = (strcmp(sin, sout) == 0); 
    CFX_PRINT_DBG("test str: in:\n%s \n%s\nout.. %s\n", sin, sout,
        ok ? "ok":"NOT ok");
    assert(ok);
    PRINT_TEST(1);
}

static void test_cache(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    assert(b.cache == NULL);
    cfx_big_enable_fac(&b);
    assert(b.cache != NULL);

    cfx_big_set_val(&b, 1);
    // assert(b.cache->primes.data == NULL);
    // assert(b.cache->state == CFX_FAC_FULL); todo
    PRINT_TEST(1);
}

static void test_zero_right(void) {
    cfx_big_t b, m;
    cfx_big_init(&b);
    cfx_big_set_val(&b, 123);
    cfx_big_init(&m);
    cfx_big_set_val(&m, 0);

    cfx_big_mul(&b, &m);
    assert(cfx_big_is_zero(&b));

    cfx_big_free(&b);
    cfx_big_free(&m);
    PRINT_TEST(1);
}

static void test_zero_left(void) {
    cfx_big_t b, m;
    cfx_big_init(&b);
    cfx_big_set_val(&b, 0);
    cfx_big_init(&m);
    cfx_big_set_val(&m, 987);

    cfx_big_mul(&b, &m);
    
    assert(cfx_big_is_zero(&b));

    cfx_big_free(&b);
    cfx_big_free(&m);
    PRINT_TEST(1);
}

static void big_expect_limbs(const char* s, const cfx_big_t* b, const uint64_t* limbs, size_t n) {
    if (b->n != n) {
        printf("[%s]: size mismatch! n: %zu, b->n: %zu!\n", s, n, b->n);
    }
    assert(b->n == n);
    int ok = 1;
    for (size_t i = 0; i < n; ++i) {
        if (b->limb[i] != limbs[i]) {
            printf("[%s]: limb mismatch at idx %zu!\n", s, i);
            ok = 0;
            break;
        }
    }
    assert(ok);
    PRINT_TEST(1);
}

static void test_mul1(void) {
    cfx_big_t b, m;
    cfx_big_init(&b);
    cfx_big_set_val(&b, 5);
    cfx_big_init(&m);
    cfx_big_set_val(&m, 7);

    cfx_big_mul(&b, &m);
    uint64_t expect[] = {35};
    CFX_BIG_PRINT(b);
    big_expect_limbs(__func__, &b, expect, 1);

    cfx_big_free(&b);
    cfx_big_free(&m);
    PRINT_TEST(1);
}

static void test_mul_adduiv(void) {
    cfx_big_t a, m;
    cfx_big_init(&a);
    cfx_big_set_val(&a, 0);
    cfx_big_init(&m);
    cfx_big_set_val(&m, 1);
    uint64_t K = 0x1B2CDE;
    cfx_big_mul_sm(&m, K);

    for (uint64_t i = 0; i < K; ++i) {
        cfx_big_add_sm(&a, 1);
    }
    assert(cfx_big_eq(&a, &m));
    PRINT_TEST(1);
}

static void test_carry_two_limbs_times_2(void) {
    cfx_big_t b, m;
    uint64_t limbs_b[] = {0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFFFFFFull};
    cfx_big_from_limbs(&b, limbs_b, 2);
    cfx_big_init(&m);
    cfx_big_set_val(&m, 2);

    cfx_big_mul(&b, &m);
    uint64_t expect[] = {0xFFFFFFFFFFFFFFFEull, 0xFFFFFFFFFFFFFFFFull, 1ull};
    big_expect_limbs(__func__, &b, expect, 3);

    cfx_big_free(&b);
    cfx_big_free(&m);
    PRINT_TEST(1);
}


static void test_mul_by_base_2_64_shift(void) {
    // Multiply by 2^64 (limbs = [0,1]) should shift by one limb
    cfx_big_t b, m;
    uint64_t limbs_b[] = {0x0123456789ABCDEFull, 0x0FEDCBA987654321ull, 0x0000000000000001ull};
    size_t sz0 = sizeof(limbs_b)/sizeof(limbs_b[0]);
    cfx_big_from_limbs(&b, limbs_b, sz0);
    uint64_t limbs_m[] = {0ull, 1ull};
    cfx_big_from_limbs(&m, limbs_m, 2);

    const size_t N = 1111;

    for (size_t n = 1; n < N; ++n) {

        cfx_big_mul(&b, &m);
        uint64_t* expect = (uint64_t*)malloc((n + sz0)*sizeof(uint64_t));
        for (size_t k = 0; k < n; ++k) {
            expect[k] = 0;
        }
        expect[n+sz0-3] = 0x0123456789ABCDEFull,
        expect[n+sz0-2] = 0x0FEDCBA987654321ull;
        expect[n+sz0-1] = 0x0000000000000001ull;

        big_expect_limbs(__func__, &b, expect, sz0 + n);
        // printf("[test_mul_by_base_2_64_shift]: tested shift of %zu OK\n", n);
        free(expect);
        expect = NULL;
    }

    cfx_big_free(&b);
    cfx_big_free(&m);
    PRINT_TEST(1);
}

static void test_self_multiply_square(void) {
    // (2^64 - 1)^2 = [0xFFFFFFFFFFFFFFFE, 1] in base 2^64
    cfx_big_t b;
    uint64_t limbs[] = {0xFFFFFFFFFFFFFFFFull};
    cfx_big_from_limbs(&b, limbs, 1);
    big_expect_limbs(__func__, &b, limbs, 1);
    cfx_big_mul(&b, &b); // self-mul path
    uint64_t expect[] = {1ull, 0xFFFFFFFFFFFFFFFEull};
    big_expect_limbs(__func__, &b, expect, 2);
    cfx_big_free(&b);
    PRINT_TEST(1);
}

static void test_self_multiply_big(void) {
    cfx_big_t b;
    size_t N = 10;
    uint64_t* limbs = (uint64_t*)malloc(N * sizeof(uint64_t));
    for (size_t i = 0; i < N; ++i) {
        limbs[i] = 2;
    }
    cfx_big_from_limbs(&b, limbs, N);
    printf("before: "); CFX_BIG_PRINT(b);
    cfx_big_mul(&b, &b);
    printf("after: "); CFX_BIG_PRINT(b);
    cfx_big_free(&b);
    free(limbs);
    limbs = NULL;
    PRINT_TEST(1);
}


void test_mul(cfx_big_t* b, const cfx_big_t* m) {
    if(cfx_big_is_zero(m)) {
        printf("multiplying b by zero!\n");
        cfx_big_free(b);
        cfx_big_set_val(b, 0);
        return;
    }
    // if (cfx_big_eq(b, m)) {
    //     cfx_big_sq(b);
    //     return;
    // }
    PRINT_TEST(1);
    
}

static void test_known_squares(void) {
    cfx_big_t b;
    cfx_big_from_str(&b, 
        "12554203470773361528352143580257209"
        "759168353591939024551938");
    cfx_big_sq(&b);
    char* s = cfx_big_to_str(&b, NULL);
    char* expect = 
        "15760802478557791686620405668794173"
        "58808686358020659979337980952437065"
        "51032029871143396883518477623727858"
        "361659555844";
    assert(strcmp(s, expect) == 0);
    
    /* ----------------------------------- */
    cfx_big_from_str(&b,
        "14536774485912137811774516281687980"
        "27136112775646765168338161504493023"
        "20618275753867907499968765052767290"
        "9553511701551770336148695547906");
    cfx_big_sq(&b);
    free(s);
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
    assert(strcmp(s, expect) == 0);

    /* ----------------------------------- */
    cfx_big_from_str(&b,
        "788040123927889584288300542721290477751811"
        "061317446299559616920266928982601233965513"
        "84510299169195903266945438318594");
    cfx_big_sq(&b);
    free(s);
    s = cfx_big_to_str(&b, NULL);
    
    expect = 
        "621007236920283574126921534924820203000914"
        "006088690925637430299255754323320524243949"
        "820727721413359346225814236510862646127210"
        "854935998980345500359712802179223868607115"
        "902506409479615101839926583004283044656249"
        "4079870273849846136836";
    assert(strcmp(s, expect) == 0);

    /* ----------------------------------- */
    cfx_big_from_str(&b,
        "4946608029462090681478206578991795742708644"
        "2564742658586426229545514803499564697000372"
        "3095350971345437292114654548843072761868784"
        "674125049315509629339381027496416991005114370");
    cfx_big_sq(&b); // 1
    free(s);
    s = cfx_big_to_str(&b, NULL);
    
    expect = 
        "24468930997138827791465884302231872460739875"
        "68207399045598121707346936565872814004700782"
        "89425671088010093095401304027496760620656589"
        "55442387983480255066626620918677531232865400"
        "98771942569641370987258548756784517935582567"
        "13069390644656846533223333372022161629347849"
        "16188372769675748466345955863374240885700087"
        "1949664098243244687781332547496780496900";
    assert(strcmp(s, expect) == 0);

    /* ----------------------------------- */
    cfx_big_sq(&b); // 2
    free(s);
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
    assert(strcmp(s, expect) == 0);

    cfx_big_sq(&b); // 4
    free(s);
    s = cfx_big_to_str(&b, NULL);
    
    expect = 
        "3584759174695717079200799669904391027371015034"
        "5591263067167166620789779901562032648556807601"
        "2950877610090545985260292606958082794183173841"
        "0531271568492233387532134452461086273878119171"
        "5847641737220633768273029588942704704129884200"
        "9958682084547837574688649786501681990292232657"
        "7415652116159576158769762214147356000386767514"
        "9511670598205070085582557818076333629441697422"
        "9875224640887417011148169060064563559950573932"
        "5226205227402913380194309571812791018912125575"
        "1791447739605363127819124149182487530334381706"
        "8896663768263012311430826689280351346703702625"
        "8078023280801979632626881628873401907467070917"
        "5815446662382885870340476544847226508933219513"
        "1908572998226578819696551060038441941135034737"
        "0068773357966050317170281108465947743312221154"
        "8985206867725840086507000464147209225780870815"
        "2470140579958636771833023777601016674168494135"
        "5848715099748171254374285813092185256574867981"
        "4922088922467495708073886349364048519341765660"
        "2078730664244522169366402136000455068825485572"
        "3466936983629748500464140600538157767988061984"
        "0563990463715200572971418532092893645704169342"
        "4426635419183083905555687209236852205715503369"
        "2336623040698494499300396883453818464450499457"
        "8627966191480124089611256089345895606330526380"
        "0963369509588825439132963159814190607878092811"
        "6978547138214917414204956687046025729039394600"
        "8260306269245606621871407784590216117503519011"
        "6801526199652563056546897824927378333686359035"
        "2100000000";
    printf("\n\n%s\n\n", s);
    assert(strcmp(s, expect) == 0);
    // CFX_BIG_PRINT_LIMBS(b);

    
    // cfx_big_sq(&b); // 8
    cfx_big_mul(&b, &b);
    free(s);
    s = cfx_big_to_str(&b, NULL);
    //// sanity check:
    cfx_big_t B;
    cfx_big_from_str(&B, s);
    assert(cfx_big_eq(&B, &b));
    char* sanity = cfx_big_to_str(&B, NULL);
    assert(strcmp(s, sanity) == 0);
    expect = "12850498340565118640831124663943566637163477270761965804775357209550765230785703949035729472531170991846950765062231843931603991577837122702202578619413163343802199724298959199688594791894654764450428039375530903550022460719636135025802643105174385371702288121193336984280156229260813405468551910887207681430170841269088143925314707610300410631962624300469472858167288519605269872627946462373021700902442150967139315116172545474219540016600664971616358393885751109046842533050183826430384246617880667515303122705936214635228967259835005388450384900687056571344933537245392429971062926387006758181202634581358151550436165392630946855812938662698021351517729227827971635271649805374176176940281753018739171849301949232131565817532850243174625134450254366859930465467865935871669065946539235916390634180117949435696071838514048795111101267436279103470592457465410839430058374274997666026789906916220271637971862818408213637491993426120588003613013514121426918778425609885983132067309940493824037246336471087139347566932289527147579507444510833470558987620105627392491106342931918229603332167519231845701934759523730028587564778522689440283636622636742434431091107839566802176801151089094414637433966684982281939928036568083598873318392485777034090015016265841172493029035363177944812390618909768667573763207436642271055957361540865636837910531318607009603088734768184385262471743568922118461465472820244894124358177759866070832621634448723670657830089717393358547747703735801523850299292208462216521430657525373415183422463571048104152842920726565698681208502981595141456253796350916896570590072884174607249064969441514018160468583314653443352671437479825779995459466312201730735736110098653787331865116531446616641987026900818215864623614884767160893264417488979936319828793408610211926232514931040704633236494245315339349224851289392387332474248568479216841581124860328079900988652092392614618200827173162666704740632419460637566175022280383967839231370086128119069239566679344622785639260673870367122861560404871658960189532546990946727895958575281330314215816959430714079369740156022375856316371677709024505475091778006456248457155666092229065835357913477405119329372494342721686545800319902824829662799906780886313825440342423658177641941252419690985063319979924638260305722020541317833802446921285816583681614997935512759860955018042317788791616381896041300997716291956236555872471276333548522003017903198166576927518763140985360706323427487549508484530837381768891367647567673106372095645829146242426938921111540219157093806251149323044120727645913143734374450245354473392643178616466237265916832878586939031607485456329408461025836997224109513396266781143113993671870210102616033541932236431301110708951036119790977069093536028698906465137616493507440201974410000000000000000";
    cfx_big_from_str(&B, expect);
    free(sanity);
    sanity = cfx_big_to_str(&B, NULL);
    assert(strcmp(sanity, expect) == 0);

    for (size_t i = 0; i < b.n; ++i) {
        printf("calculated: b.limb[%zu]: %llu; correct: B.limb[%zu]: %llu: %s: diff %d\n",
        i, b.limb[i], i, B.limb[i], b.limb[i] == B.limb[i] ? "ok" : "--- NOT OK",
        (int)(b.limb[i] & 0xFFFF) - (int)(B.limb[i] & 0xFFFF));
    }
    printf("\n\n");
    assert(cfx_big_eq(&B, &b));
    assert(strcmp(s, expect) == 0);
    int cnt = 0;
    cfx_big_mul(&b, &b); // 16
    printf("mul %d len: %zu \n", ++cnt, b.n);
    cfx_big_mul(&b, &b); // 32
    printf("mul %d len: %zu \n", ++cnt, b.n);
    cfx_big_mul(&b, &b); // 64
    printf("mul %d len: %zu \n", ++cnt, b.n);
    cfx_big_mul(&b, &b); // 128
    printf("mul %d len: %zu \n", ++cnt, b.n);
    cfx_big_mul(&b, &b); // 256
    printf("mul %d len: %zu \n", ++cnt, b.n);
    // cfx_big_mul(&b, &b); // 512
    // printf("mul %d len: %zu \n", ++cnt, b.n);
    // cfx_big_mul(&b, &b); // 1024
    // printf("mul %d len: %zu \n", ++cnt, b.n);
    // cfx_big_mul(&b, &b); // 2048
    // printf("mul %d len: %zu \n", ++cnt, b.n);
    // cfx_big_mul(&b, &b); // 4096
    // printf("mul %d len: %zu \n", ++cnt, b.n);

    for (size_t i = 0; i < b.n; ++i) {
        printf("calculated: b.limb[%zu]: %llu (0x%llx)\n", i, b.limb[i], b.limb[i]);
    }
    
    char* huge = cfx_big_to_str(&b, NULL);
    printf("%s\n", huge);
    printf("digits: %zu\n", strlen(huge));
    free(huge);
    free(s);
    cfx_big_free(&b);
    PRINT_TEST(1);
}

int main(void) {
    test_copy_swap();
    test_cfx_big_init();
    test_cfx_big_reserve();
    test_mul_by_zero();
    test_limb1();
    test_limb2();
    test_limb3();
    test_limb4();
    test_limb5();
    test_limb6();
    test_limb7();
    test_str1();
    test_str2();
    test_add_sm();
    test_sub_sm();
    test_cache();
    test_zero_right();
    test_zero_left();
    test_mul1();
    test_carry_two_limbs_times_2();
    test_mul_by_base_2_64_shift();
    test_self_multiply_square();
    test_self_multiply_big();
    test_known_squares();
    test_mul_adduiv();
    return 0;
}
