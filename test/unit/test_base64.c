#include "cfx/base64.h"
#include "cfx/macros.h"
#include "cfx/rand.h"

#include <string.h>
#include <stdint.h>
#include <stddef.h>


static void test_base64_rfc4648_vectors(void) {
    struct {
        const char *plain;
        const char *b64;
    } v[] = {
        { "",       ""          },
        { "f",      "Zg=="      },
        { "fo",     "Zm8="      },
        { "foo",    "Zm9v"      },
        { "foob",   "Zm9vYg=="  },
        { "fooba",  "Zm9vYmE="  },
        { "foobar", "Zm9vYmFy"  },
    };

    for (size_t i = 0; i < sizeof(v)/sizeof(v[0]); i++) {
        const uint8_t *in = (const uint8_t *)v[i].plain;
        size_t in_len = strlen(v[i].plain);

        /* encode */
        {
            char out[128];
            size_t out_len = sizeof(out);
            int rc = cfx_base64_encode(out, &out_len, in, in_len);
            CFX_ASSERT(rc == 0);
            CFX_ASSERT(out_len == strlen(v[i].b64));
            CFX_ASSERT(memcmp(out, v[i].b64, out_len) == 0);

            /* length query */
            size_t need = 0;
            rc = cfx_base64_encode(NULL, &need, in, in_len);
            CFX_ASSERT(rc == 0);
            CFX_ASSERT(need == out_len);
        }

        /* decode */
        {
            uint8_t out[128];
            size_t out_len = sizeof(out);
            int rc = cfx_base64_decode(out, &out_len, v[i].b64, strlen(v[i].b64));
            CFX_ASSERT(rc == 0);
            CFX_ASSERT(out_len == in_len);
            CFX_ASSERT(memcmp(out, in, in_len) == 0);

            /* length query */
            size_t need = 0;
            rc = cfx_base64_decode(NULL, &need, v[i].b64, strlen(v[i].b64));
            CFX_ASSERT(rc == 0);
            CFX_ASSERT(need == in_len);
        }
    }
}

static void test_base64_binary_known_cases(void) {
    /* 0x00 */
    {
        const uint8_t in[] = { 0x00 };
        const char *exp = "AA==";

        char b64[16];
        size_t b64_len = sizeof(b64);
        int rc = cfx_base64_encode(b64, &b64_len, in, sizeof(in));
        CFX_ASSERT(rc == 0);
        CFX_ASSERT(b64_len == strlen(exp));
        CFX_ASSERT(memcmp(b64, exp, b64_len) == 0);

        uint8_t out[16];
        size_t out_len = sizeof(out);
        rc = cfx_base64_decode(out, &out_len, exp, strlen(exp));
        CFX_ASSERT(rc == 0);
        CFX_ASSERT(out_len == sizeof(in));
        CFX_ASSERT(memcmp(out, in, sizeof(in)) == 0);
    }

    /* 0x00 0x00 */
    {
        const uint8_t in[] = { 0x00, 0x00 };
        const char *exp = "AAA=";

        char b64[16];
        size_t b64_len = sizeof(b64);
        int rc = cfx_base64_encode(b64, &b64_len, in, sizeof(in));
        CFX_ASSERT(rc == 0);
        CFX_ASSERT(b64_len == strlen(exp));
        CFX_ASSERT(memcmp(b64, exp, b64_len) == 0);

        uint8_t out[16];
        size_t out_len = sizeof(out);
        rc = cfx_base64_decode(out, &out_len, exp, strlen(exp));
        CFX_ASSERT(rc == 0);
        CFX_ASSERT(out_len == sizeof(in));
        CFX_ASSERT(memcmp(out, in, sizeof(in)) == 0);
    }

    /* 0x00 0x00 0x00 */
    {
        const uint8_t in[] = { 0x00, 0x00, 0x00 };
        const char *exp = "AAAA";

        char b64[16];
        size_t b64_len = sizeof(b64);
        int rc = cfx_base64_encode(b64, &b64_len, in, sizeof(in));
        CFX_ASSERT(rc == 0);
        CFX_ASSERT(b64_len == strlen(exp));
        CFX_ASSERT(memcmp(b64, exp, b64_len) == 0);

        uint8_t out[16];
        size_t out_len = sizeof(out);
        rc = cfx_base64_decode(out, &out_len, exp, strlen(exp));
        CFX_ASSERT(rc == 0);
        CFX_ASSERT(out_len == sizeof(in));
        CFX_ASSERT(memcmp(out, in, sizeof(in)) == 0);
    }

    /* roundtrip some high bytes */
    {
        const uint8_t in[] = { 0xFF, 0xEE, 0xDD, 0xCC, 0xBB, 0xAA };

        char b64[64];
        size_t b64_len = sizeof(b64);
        int rc = cfx_base64_encode(b64, &b64_len, in, sizeof(in));
        CFX_ASSERT(rc == 0);

        uint8_t out[64];
        size_t out_len = sizeof(out);
        rc = cfx_base64_decode(out, &out_len, b64, b64_len);
        CFX_ASSERT(rc == 0);
        CFX_ASSERT(out_len == sizeof(in));
        CFX_ASSERT(memcmp(out, in, sizeof(in)) == 0);
    }
}

static void test_base64_decode_ignores_whitespace(void) {
    /* "foobar" -> Zm9vYmFy */
    {
        const uint8_t exp[] = { 'f','o','o','b','a','r' };
        const char *b64_ws = " Z m 9 v \n Y m F y \r";

        uint8_t out[32];
        size_t out_len = sizeof(out);
        int rc = cfx_base64_decode(out, &out_len, b64_ws, strlen(b64_ws));
        CFX_ASSERT(rc == 0);
        CFX_ASSERT(out_len == sizeof(exp));
        CFX_ASSERT(memcmp(out, exp, sizeof(exp)) == 0);
    }

    /* "fooba" -> Zm9vYmE= with whitespace */
    {
        const uint8_t exp[] = { 'f','o','o','b','a' };
        const char *b64_ws = "\n Zm9vYmE= \t\r\n";

        uint8_t out[32];
        size_t out_len = sizeof(out);
        int rc = cfx_base64_decode(out, &out_len, b64_ws, strlen(b64_ws));
        CFX_ASSERT(rc == 0);
        CFX_ASSERT(out_len == sizeof(exp));
        CFX_ASSERT(memcmp(out, exp, sizeof(exp)) == 0);
    }
}

static void test_base64_rejects_invalid_inputs(void) {
    /* length not multiple of 4 (non-ws chars) */
    {
        const char *bad[] = { "Zg=", "Z", "Zm9vYg=" }; /* 3,1,7 non-ws */
        for (size_t i = 0; i < sizeof(bad)/sizeof(bad[0]); i++) {
            uint8_t out[32];
            size_t out_len = sizeof(out);
            int rc = cfx_base64_decode(out, &out_len, bad[i], strlen(bad[i]));
            CFX_ASSERT(rc != 0);
        }
    }

    /* invalid characters */
    {
        const char *bad[] = { "Zm9v*", "Zm9v_===" };
        for (size_t i = 0; i < sizeof(bad)/sizeof(bad[0]); i++) {
            uint8_t out[32];
            size_t out_len = sizeof(out);
            int rc = cfx_base64_decode(out, &out_len, bad[i], strlen(bad[i]));
            CFX_ASSERT(rc != 0);
        }
    }

    /* '=' placement rules and trailing junk after padding */
    {
        const char *bad[] = {
            "=m9v",
            "Zm=9",
            "Zg===",
            "Zg==Zg==",
            "Zg==\nZg=="
        };
        for (size_t i = 0; i < sizeof(bad)/sizeof(bad[0]); i++) {
            uint8_t out[32];
            size_t out_len = sizeof(out);
            int rc = cfx_base64_decode(out, &out_len, bad[i], strlen(bad[i]));
            CFX_ASSERT(rc != 0);;
        }
    }

    /* non-canonical tail bits (your decoder enforces this) */
    {
        const char *bad = "Zg=A";
        uint8_t out[32];
        size_t out_len = sizeof(out);
        int rc = cfx_base64_decode(out, &out_len, bad, strlen(bad));
        CFX_ASSERT(rc != 0);;

        /* control: canonical "Zg==" must work */
        const char *good = "Zg==";
        out_len = sizeof(out);
        rc = cfx_base64_decode(out, &out_len, good, strlen(good));
        CFX_ASSERT(rc == 0);
        CFX_ASSERT(out_len == 1);
        CFX_ASSERT(out[0] == (uint8_t)'f');
    }
}

static void test_base64_small_buffer_errors(void) {
    /* encode needs 8 for "foobar" -> 8 chars */
    {
        const uint8_t in[] = { 'f','o','o','b','a','r' };
        char out[8];
        size_t out_len = 7; /* too small */
        int rc = cfx_base64_encode(out, &out_len, in, sizeof(in));
        CFX_ASSERT(rc != 0);;
    }

    /* decode needs 6 for "Zm9vYmFy" */
    {
        const char *b64 = "Zm9vYmFy";
        uint8_t out[2];
        size_t out_len = sizeof(out);
        int rc = cfx_base64_decode(out, &out_len, b64, strlen(b64));
        CFX_ASSERT(rc != 0);;
    }
}

static void test_base64_roundtrip_many_lengths(void) {
    uint8_t bin[256];
    char b64[512];
    uint8_t out[256];

    for (size_t n = 0; n < sizeof(bin); n++) {
        cfx_pcg32_seed(n);
        cfx_pcg32_bytes(bin,n);
        for (size_t k = 0; k < n; ++k) {
            printf("%02x", bin[k]);
        }
        printf("\n");
        size_t b64_len = sizeof(b64);
        int rc = cfx_base64_encode(b64, &b64_len, bin, n);
        CFX_ASSERT(rc == 0);

        size_t out_len = sizeof(out);
        rc = cfx_base64_decode(out, &out_len, b64, b64_len);
        CFX_ASSERT(rc == 0);

        CFX_ASSERT(out_len == n);
        CFX_ASSERT((n == 0) || (memcmp(out, bin, n) == 0));
    }
}


int main(void) {
    CFX_TEST(test_base64_rfc4648_vectors);
    CFX_TEST(test_base64_binary_known_cases);
    CFX_TEST(test_base64_decode_ignores_whitespace);
    CFX_TEST(test_base64_rejects_invalid_inputs);
    CFX_TEST(test_base64_small_buffer_errors);
    CFX_TEST(test_base64_roundtrip_many_lengths);
    return 0;
}
