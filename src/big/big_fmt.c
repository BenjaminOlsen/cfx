/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#include "big_internal.h"

#include "cfx/algo.h"
#include "cfx/base64.h"

#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <inttypes.h>

static inline size_t hex_digits_limb(cfx_limb_t v) {
    if (!v) return 1;
#if defined(__GNUC__) || defined(__clang__)
    #if CFX_LIMB_BITS == 64
    unsigned lead = (unsigned)__builtin_clzll(v);
    #else
    unsigned lead = (unsigned)__builtin_clz(v);
    #endif
    unsigned bits = CFX_LIMB_BITS - lead;
    return (bits + 3u) / 4u; /* ceil(bits/4) */
#else
    size_t d = 0;
    while (v) {
        v >>= 4; ++d;
    }
    return d;
#endif
}
static size_t bits_in_limb(cfx_limb_t x) {
    if (!x) return 0;
#if defined(__GNUC__) || defined(__clang__)
    #if CFX_LIMB_BITS == 64
    return CFX_LIMB_BITS - (size_t)__builtin_clzll(x);
    #else
    return CFX_LIMB_BITS - (size_t)__builtin_clz(x);
    #endif
#else
    size_t n = 0;
    while (x) {
        ++n; x >>= 1;
    }
    return n;
#endif
}

char * cfx_big_b64_alloc(const cfx_big_t *src, size_t *sz_out) {
    if (!src || src->n == 0) {
        char *s = (char *)malloc(2);
        if (!s) return NULL;
        s[0] = '0';
        s[1] = '\0';
        if (sz_out) *sz_out = 1;
        return s;
    }

    size_t bytecnt = 0;
    cfx_big_to_bytes_be(NULL, &bytecnt, src);
    uint8_t *bytes = (uint8_t *)malloc(bytecnt);
    cfx_big_to_bytes_be(bytes, &bytecnt, src);
    size_t charcnt = 0;
    cfx_base64_encode(NULL, &charcnt, bytes, bytecnt);
    char *s = (char *)malloc(charcnt + 1);
    cfx_base64_encode(s, &charcnt, bytes, bytecnt);
    s[charcnt] = '\0';
    free(bytes);
    if (sz_out) *sz_out = charcnt;
    return s;
}

char * cfx_big_bin_alloc(const cfx_big_t *src, size_t *sz_out) {
    if (!src || src->n == 0) {
        char *s = (char *)malloc(2);
        if (!s) return NULL;
        s[0] = '0';
        s[1] = '\0';
        if (sz_out) *sz_out = 1;
        return s;
    }

    const size_t ms_idx  = src->n - 1;              /* most-significant limb index */
    const cfx_limb_t msval  = src->limb[ms_idx];
    const size_t ms_bits  = bits_in_limb(msval);      /* 1..CFX_LIMB_BITS */
    const size_t total_len = ms_bits + (size_t)CFX_LIMB_BITS * ms_idx; /* total bits as characters */

    char *s = (char *)malloc(total_len + 1);
    if (!s) return NULL;

    size_t pos = 0;
    const char bch[2] = {'0', '1'};

    for (size_t b = ms_bits; b-- > 0; ) {
        s[pos++] = bch[(size_t)((msval >> b) & 1u)];
    }


    for (size_t i = ms_idx; i-- > 0; ) {
        cfx_limb_t limb = src->limb[i];
        for (size_t b = CFX_LIMB_BITS; b-- > 0; ) {
            s[pos++] = bch[(size_t)((limb >> b) & 1u)];
        }
    }

    s[total_len] = '\0';
    if (sz_out) *sz_out = total_len;
    return s;
}

char * cfx_big_hex_alloc(const cfx_big_t *src, size_t *sz_out) {
    /* Treat empty/zero as "0" */
    if (!src || src->n == 0) {
        char *s = (char *)malloc(2);
        if (!s) return NULL;
        s[0] = '0'; s[1] = '\0';
        if (sz_out) *sz_out = 1;
        return s;
    }

    /* trim leading zeros */
    size_t ms = src->n;

    const cfx_limb_t ms_val = src->limb[ms - 1];
    const size_t ms_digits = hex_digits_limb(ms_val);
    const size_t hex_per_limb = CFX_LIMB_BITS / 4;  /* 8 for 32-bit, 16 for 64-bit */
    const size_t total_len = ms_digits + (ms - 1) * hex_per_limb;

    char *s = (char *)malloc(total_len + 1); /* +1 for NUL */
    if (!s) return NULL;

    char *p = s;
    size_t rem = total_len + 1;

    /* Most-significant limb without leading zeros */
    int written = snprintf(p, rem, CFX_PRIxLIMB, (cfx_limb_t)ms_val);
    assert(written > 0 && (size_t)written == ms_digits);
    p   += written;
    rem -= (size_t)written;

    /* Remaining limbs, zero-padded to hex_per_limb hex chars each */
    #ifdef CFX_DEBUG
    size_t k = 0;
    size_t cnt = 0;
    #endif
    for (size_t i = ms-1; i--;) {
        written = snprintf(p, rem, "" CFX_PRI0xLIMB, (cfx_limb_t)src->limb[i]);
        assert(written == (int)hex_per_limb);
        p   += written;
        #ifdef CFX_DEBUG
        cnt += written;
        #endif
        rem -= (size_t)written;

        #ifdef CFX_DEBUG
        if (!(cnt % 7)) {  /* some number.. */
            const char spinner[] = "|/-\\";
            ++k;
            CFX_PRINT_DBG("%zu hex digits done... %zu/%zu limbs remain... %c        \r",
                cnt, rem, total_len, spinner[k % 4]);
            fflush(stdout);
        }
        #endif
    }
    CFX_PRINT_DBG("\n");
    /* `snprintf` already wrote the final '\0' on the last call */
    if (sz_out) *sz_out = total_len;
    return s;
}

/* Convert cfx_big_t to decimal string */
char * cfx_big_dec_alloc(const cfx_big_t *src, size_t *sz_out) {
    if (src->n == 0) {
        char *s = (char *)malloc(2);
        s[0]='0';
        s[1]='\0';
        if (sz_out) *sz_out = 1;
        return s;
    }

    cfx_big_t tmp = *src;
    tmp.cap = tmp.n;
    tmp.limb = (cfx_limb_t *)malloc(tmp.n * sizeof(cfx_limb_t));
    memcpy(tmp.limb, src->limb, tmp.n * sizeof(cfx_limb_t));

    enum { CHUNK_BASE = 1000000000u, CHUNK_DIGS = 9 };
    size_t maxdig = src->n * 20; /* log10(2^64) == 19.2659... */
    if (maxdig / 20 != src->n) {
        /* overflow in src->n * 20 */
        if (sz_out) *sz_out = 0;
        return NULL;
    }
    uint32_t *chunks = (uint32_t *)malloc(maxdig);
    size_t k = 0;

    /* printf("[cfx_big_dec_alloc] building base %u representation... max digits: %zu\n",CHUNK_BASE, maxdig); */
    #ifdef CFX_DEBUG
    int cnt = 0;
    const size_t n0 = tmp.n;
    #endif
    while (tmp.n) {

        #ifdef CFX_DEBUG
        if (!(tmp.n % 100)) {
            const char spinner[] = "|/-\\";
            CFX_PRINT_DBG("%zu decimal digits done... %zu/%zu limbs remain... %c        \r",
                k*9, tmp.n, n0, spinner[cnt++ % 4]);
            fflush(stdout);
        }
        #endif
        chunks[k++] = cfx_big_div_sm_u32_eq(&tmp, CHUNK_BASE);
    }
    /* printf("\n"); */
    cfx_big_free(&tmp);

    /* build string */
    /* first chunk has no zero-padding, others are 9-digit padded */
    size_t len = 0;
    {
        uint32_t first = chunks[k-1];
        char buf[16];
        int l = snprintf(buf, sizeof buf, "%" PRIu32, first);
        len += l + (k-1) * CHUNK_DIGS;
    }
    char *s = (char *)malloc(len+1);
    char *p = s;
    size_t rem = len + 1;

    /* write first */
    int w = snprintf(p, rem, "%" PRIu32, chunks[k-1]);
    p += w;
    rem -= (size_t)w;
    /* write rest padded (indices k-2 -> 0) */
    for (size_t i = k-1; i--;) {
        w = snprintf(p, rem, "%09" PRIu32, chunks[i]);
        p += w; rem -= (size_t)w;
    }
    s[len] = '\0';
    if (sz_out) *sz_out = len;
    free(chunks);
    return s;
}

size_t cfx_big_snprint_dec(const cfx_big_t *b, char *out, size_t outlen) {
    size_t len;
    char *s = cfx_big_dec_alloc(b, &len);
    if (!s) {
        if (out && outlen > 0) out[0] = '\0';
        return 0;
    }
    if (out && outlen > 0) {
        size_t copy = (len < outlen) ? len : outlen - 1;
        memcpy(out, s, copy);
        out[copy] = '\0';
    }
    free(s);
    return len;
}

size_t cfx_big_snprint_hex(const cfx_big_t *b, char *out, size_t outlen) {
    size_t len;
    char *s = cfx_big_hex_alloc(b, &len);
    if (!s) {
        if (out && outlen > 0) out[0] = '\0';
        return 0;
    }
    if (out && outlen > 0) {
        size_t copy = (len < outlen) ? len : outlen - 1;
        memcpy(out, s, copy);
        out[copy] = '\0';
    }
    free(s);
    return len;
}

size_t cfx_big_snprint_bin(const cfx_big_t *b, char *out, size_t outlen) {
    size_t len;
    char *s = cfx_big_bin_alloc(b, &len);
    if (!s) {
        if (out && outlen > 0) out[0] = '\0';
        return 0;
    }
    if (out && outlen > 0) {
        size_t copy = (len < outlen) ? len : outlen - 1;
        memcpy(out, s, copy);
        out[copy] = '\0';
    }
    free(s);
    return len;
}

/* Scan a single numeric literal token at in[0..in_len).
   Accepts optional prefixes: 0x/0X, 0b/0B, 0o/0O.
   No internal whitespace, no sign.
   On success: parses into out, sets *consumed to total chars consumed, returns 0.
   On failure (not a number at start): sets *consumed=0, returns -1. */
int cfx_big_scan_num_n(cfx_big_t *out, const uint8_t *in, size_t in_len, size_t *consumed) {
    if (consumed) *consumed = 0;
    if (!out || (!in && in_len)) return -1;
    if (in_len == 0) return -1;

    size_t prefix = 0;
    int ret = 0;
    size_t digits = 0;

    /* base prefix detection */
    if (in_len >= 2 && in[0] == (uint8_t)'0') {
        uint8_t p = in[1];
        if (p == (uint8_t)'x' || p == (uint8_t)'X' ||
            p == (uint8_t)'b' || p == (uint8_t)'B' ||
            p == (uint8_t)'o' || p == (uint8_t)'O') {
            prefix = 2;
        }
    } else if (in_len >= 4 && strncmp((const char *)in, "b64:", 4) == 0) {
        prefix = 4;
    }

    // printf("SCAN_NUM prefix=0%c, calling bin on '%c' '%c' '%c'...\n",
    //     (char)in[1],
    //     (char)in[prefix + 0],
    //     (char)in[prefix + 1],
    //     (char)in[prefix + 2]);

    if (prefix) {
        if (prefix == 2) {
            uint8_t p = in[1];
            if (p == (uint8_t)'x' || p == (uint8_t)'X') {
                ret = cfx_big_scan_hex_n(out, in + prefix, in_len - prefix, &digits);
            } else if (p == (uint8_t)'b' || p == (uint8_t)'B') {
                ret = cfx_big_scan_bin_n(out, in + prefix, in_len - prefix, &digits);
            } else { /* o/O */
                ret = cfx_big_from_oct_n(out, in + prefix, in_len - prefix, &digits);
            }
        } else if (prefix == 4) {
            ret = cfx_big_scan_b64_n(out, in + prefix, in_len - prefix, &digits);
        }

        if (ret != 0) return -1;          /* no digits after prefix */
        if (consumed) *consumed = prefix + digits;
        return 0;
    }

    /* no prefix => decimal */
    ret = cfx_big_scan_dec_n(out, in, in_len, &digits);
    if (ret != 0) return -1;
    if (consumed) *consumed = digits;
    return 0;
}

int cfx_big_from_str(cfx_big_t *out, const char *s) {
    if (!out || !s) return -1;

    size_t len = strlen(s);
    if (len == 0) return -1;

    while (len > 0 && isspace((unsigned char)s[0])) {
        ++s;
        --len;
    }

    size_t consumed = 0;  /* digits consumed */
    size_t prefix_len = 0;
    int ret = 0;

    if (len >= 2 && (s[0]=='0' && (s[1]=='x' || s[1]=='X'))) {
        prefix_len = 2;
        ret = cfx_big_scan_hex_n(out, (const uint8_t *)(s+prefix_len), len-prefix_len, &consumed);
    } else if (len >=2 && (s[0]=='0' && (s[1]=='b' || s[1]=='B'))) {
        prefix_len = 2;
        ret = cfx_big_scan_bin_n(out, (const uint8_t *)(s+prefix_len), len-prefix_len, &consumed);
    } else if (len >= 4 && (s[0]=='b' || s[0]=='B') && s[1]=='6' && s[2]=='4' && s[3]==':') {
        prefix_len = 4;
        ret = cfx_big_scan_b64_n(out, (const uint8_t *)(s + prefix_len), len - prefix_len, &consumed);
    } else {
        prefix_len = 0;
        ret = cfx_big_scan_dec_n(out, (const uint8_t *)s, len, &consumed);
    }

    s += (prefix_len + consumed);

    while (*s && isspace((unsigned char)*s)) ++s;   /* allow trailing spaces only */
    if (ret != 0 || *s != '\0') return -1;          /* no digits or junk trailing */

    return 0;
}

static const int8_t hex_table[256] = {
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, -1,-1,-1,-1,-1,-1,
    -1,10,11,12,13,14,15,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,10,11,12,13,14,15,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1
};

int cfx_big_scan_hex_n(cfx_big_t *out, const uint8_t *in, size_t in_len, size_t *consumed) {
    cfx_big_assign_zero(out);
    if (consumed) *consumed = 0;
    size_t pos = 0;
    for (; pos < in_len; ++pos) {
        const uint8_t c = in[pos];
        int d = hex_table[c];
        if (d < 0) break;
        cfx_big_shl_bits(out, out, 4);
        cfx_big_add_sm_eq(out, (cfx_limb_t)d);
    }
    if (pos == 0) return -1;
    if (consumed) *consumed = pos;
    return 0;
}

int cfx_big_from_hex(cfx_big_t *out, const char *s) {
    while (*s && isspace((unsigned char)*s)) ++s;

    if (s[0]=='0' && (s[1]=='x' || s[1]=='X')) s += 2;
    size_t len = strlen(s);
    size_t digits = 0;
    int ret = cfx_big_scan_hex_n(out, (const uint8_t *)s, len, &digits);
    s += digits;
    while (isspace((unsigned char)*s)) ++s;
    if (ret != 0 || *s != '\0') return -1;
    return 0;
}

int cfx_big_scan_bin_n(cfx_big_t *out, const uint8_t *in, size_t in_len, size_t *consumed) {
    size_t ndig = 0;

    if (consumed) *consumed = 0;
    if (!out || (!in && in_len)) return -1;

    while (ndig < in_len) {
        uint8_t c = in[ndig];
        if (c != (uint8_t)'0' && c != (uint8_t)'1') break;
        ++ndig;
    }

    if (consumed) *consumed = ndig;
    if (ndig == 0) {
        cfx_big_assign_zero(out); return -1;
    }

    cfx_big_assign_zero(out);

    size_t limb_cnt = (ndig + CFX_LIMB_BITS - 1) / CFX_LIMB_BITS;
    cfx_big_reserve(out, limb_cnt);
    out->n = limb_cnt;

    for (size_t j = 0; j < ndig; ++j) {
        if (in[ndig - 1 - j] == (uint8_t)'1') {
            out->limb[j / CFX_LIMB_BITS] |= (cfx_limb_t)1 << (j % CFX_LIMB_BITS);
        }
    }

    while (out->n > 0 && out->limb[out->n - 1] == 0) --out->n;
    return 0;
}



int cfx_big_from_oct_n(cfx_big_t *out, const uint8_t *in, size_t in_len, size_t *consumed) {
    (void)out;
    (void)in;
    (void)in_len;
    (void)consumed;
    /* TODO: implement octal parsing */
    return -1;
}

int cfx_big_from_bin(cfx_big_t *b, const char *s) {
    size_t len = strlen(s);
    if ((len > 2) && (s[0] == '0' && ((s[1] == 'b') || (s[1] == 'B')))) {
        s += 2;
        len -= 2;
    }
    size_t digits = 0;
    int ret = cfx_big_scan_bin_n(b, (const uint8_t *)s, len, &digits);
    return ret;
}

int cfx_big_scan_dec_n(cfx_big_t *out, const uint8_t *in, size_t in_len, size_t *consumed) {
    cfx_big_assign_zero(out);
    if (consumed) *consumed = 0;
    size_t pos = 0;
    while (pos < in_len) {
        if (!isdigit((unsigned char)in[pos])) break;
        uint32_t digit = (uint32_t)(in[pos] - '0');
        cfx_big_mul_sm_eq(out, 10); /* out = out * 10 */
        cfx_big_add_sm_eq(out, digit); /* out = out + digit */
        ++pos;
    }
    if (pos == 0) return -1;
    if (consumed) *consumed = pos;
    return 0;
}

int cfx_big_from_dec(cfx_big_t *out, const char *s) {
    while (isspace((unsigned char)*s)) s++;
    size_t len = strlen(s);
    size_t digits = 0;
    int ret = cfx_big_scan_dec_n(out, (const uint8_t *)s, len, &digits);
    return ret;
}

static int is_b64_char(unsigned char c) {
    if (c == '+' || c == '/' || c == '=') return 1;
    if (c >= 'A' && c <= 'Z') return 1;
    if (c >= 'a' && c <= 'z') return 1;
    if (c >= '0' && c <= '9') return 1;
    /*if (isspace(c)) return 1;*/
    return 0;
}

/* Decodes 'in' as base 64 into 'out', up until the end of token: either '=' padding or first invalid character.
   Writes consumed characters to *consumed. */
int cfx_big_scan_b64_n(cfx_big_t *out, const uint8_t *in, size_t in_len, size_t *consumed) {
    if (!out || !in || !consumed) return -1;

    size_t n = 0;        /* bytes consumed from input */
    size_t nonws = 0;    /* count of non-ws chars consumed */
    size_t eq = 0;
    int seen_eq = 0;

    /* consume one base64 token */
    while (n < in_len) {
        unsigned char c = (unsigned char)in[n];

        if (isspace(c)) {
            n++; continue;
        }                                  /* allow internal whitespace */

        if (c == '=') {
            seen_eq = 1;
            eq++;
            nonws++;
            n++;
            continue;
        }

        if (is_b64_char(c)) {
            if (seen_eq) break;      /* data after '=' ends token */
            nonws++;
            n++;
            continue;
        }

        break; /* non-b64, non-ws terminates token */
    }

    if (nonws == 0) return -1;
    if ((nonws % 4) != 0) return -1;
    if (eq > 2) return -1;

    /* now decode n bytes of base64: */
    size_t bin_len = 0;
    if (cfx_base64_decode(NULL, &bin_len, (const char *)in, n) != 0) return -1;

    uint8_t stack_buf[256];
    uint8_t *bin = stack_buf;

    if (bin_len > sizeof(stack_buf)) {
        bin = (uint8_t *)malloc(bin_len);
        if (!bin) return -1;
    }

    size_t got = bin_len;
    if (cfx_base64_decode(bin, &got, (const char *)in, n) != 0 || got != bin_len) {
        if (bin != stack_buf) free(bin);
        return -1;
    }

    if (cfx_big_from_bytes_be(out, bin, bin_len) != 0) {
        if (bin != stack_buf) free(bin);
        return -1;
    }

    if (bin != stack_buf) free(bin);

    *consumed = n;
    return 0;
}


/* Precomputed 10^k for k=0..18 */
static const cfx_limb_t POW10U64[CFX_LIMB_DIGITS_DEC + 1] = {
    1ULL,
    10ULL,
    100ULL,
    1000ULL,
    10000ULL,
    100000ULL,
    1000000ULL,
    10000000ULL,
    100000000ULL,
    1000000000ULL,
    #if (CFX_LIMB_BITS == 64)
    10000000000ULL,
    100000000000ULL,
    1000000000000ULL,
    10000000000000ULL,
    100000000000000ULL,
    1000000000000000ULL,
    10000000000000000ULL,
    100000000000000000ULL,
    1000000000000000000ULL,
    /*10000000000000000000ULL */
    #endif
};

static void flush_chunk(cfx_big_t *out, unsigned base, cfx_limb_t *chunk_val, unsigned *chunk_len) {
    if (*chunk_len == 0) return;

    if (base == 10) {
        cfx_big_mul_sm_eq(out, POW10U64[*chunk_len]);
        cfx_big_add_sm_eq(out, *chunk_val);
    } else if (base == 16) {
        unsigned bits = 4u * (*chunk_len);
        cfx_big_shl_bits_eq(out, bits);
        cfx_big_add_sm_eq(out, *chunk_val);
    } else { /* base 2 */
        unsigned bits = *chunk_len;
        cfx_big_shl_bits_eq(out, bits);
        cfx_big_add_sm_eq(out, *chunk_val);
    }

    *chunk_val = 0;
    *chunk_len = 0;
}

/* Return 0 on success, nonzero on parse error. */
int cfx_big_from_file(cfx_big_t *out, FILE *fp, int base) {
    cfx_big_from_limb(out, 0);

    enum { BASE_DEC=10, BASE_HEX=16, BASE_BIN=2 };
    int detected_base;
    if ((base != BASE_DEC) && (base != BASE_HEX) && (base != BASE_BIN)) {
        detected_base = BASE_DEC; /* default */
    } else {
        detected_base = base;
    }

    int saw_digit = 0;
    int negative = 0;
    int in_prefix = 1; /* we're skipping leading ws, sign, 0x */

    /* chunk accumulators */
    cfx_limb_t chunk_val = 0;
    unsigned chunk_len = 0; /* digits in current chunk */

    unsigned dec_chunk_max = CFX_LIMB_DIGITS_DEC;
    unsigned hex_chunk_max = CFX_LIMB_DIGITS_HEX;

    unsigned char buf[64 * 1024];
    size_t nread;
    while ((nread = fread(buf, 1, sizeof(buf), fp)) > 0) {
        for (size_t i = 0; i < nread; ++i) {
            unsigned char c = buf[i];

            /* allow underscores, quotes, spaces, newlines, tabs as visual separators */
            if ( (c == '\n') ||(c == '_') || (c == '"') || isspace(c) || (c == '\t') || (c == '\r') )  continue;

            if (in_prefix) {
                if (isspace(c)) continue;
                if (c == '+') {
                    in_prefix = 0; continue;
                }
                if (c == '-') {
                    negative = 1; in_prefix = 0; continue;
                }

                /* Base detection: 0x / 0X (hex), 0b / 0B (bin) */
                if (c == '0') {
                    /* Peek ahead safely: if at buffer end, we'll detect on next loop */
                    unsigned char next = (i + 1 < nread) ? buf[i + 1] : 0;
                    if (next == 'x' || next == 'X') {
                        if (detected_base != BASE_HEX) {
                            CFX_PRINT_ERR("hex '0x' prefix in file,"
                                " but different base (%d) specified!", detected_base);

                        }
                        detected_base = BASE_HEX;
                        i++;
                        in_prefix = 0;
                        continue;
                    }
                    if (next == 'b' || next == 'B') {
                        if (detected_base != BASE_BIN) {
                            CFX_PRINT_ERR("binary '0b' prefix in file,"
                                " but different base (%d) specified!", detected_base);
                            return -1;
                        }
                        detected_base = BASE_BIN;
                        i++;
                        in_prefix = 0;
                        continue;
                    }
                    /* Otherwise, treat as decimal 0 and fall through as a digit */
                }
                /* We've seen a non-space, non-sign */
                in_prefix = 0;
            }

            /* Digit handling by base */
            if (detected_base == BASE_DEC) {
                if (c >= '0' && c <= '9') {
                    saw_digit = 1;
                    if (chunk_len == dec_chunk_max) flush_chunk(out, detected_base, &chunk_val, &chunk_len);
                    chunk_val = chunk_val * 10u + (cfx_limb_t)(c - '0');
                    chunk_len++;
                    continue;
                }
            } else if (detected_base == BASE_HEX) {
                int v = hex_table[c];

                if (v != -1) {
                    saw_digit = 1;
                    if (chunk_len == hex_chunk_max) flush_chunk(out, detected_base, &chunk_val, &chunk_len);
                    chunk_val = (chunk_val << 4) | (cfx_limb_t)v;
                    chunk_len++;
                    continue;
                }
            } else { /* BASE_BIN */
                if (c == '0' || c == '1') {
                    saw_digit = 1;
                    if (chunk_len == 60u) flush_chunk(out, detected_base, &chunk_val, &chunk_len); /* keep under 64 bits */
                    chunk_val = (chunk_val << 1) | (cfx_limb_t)(c - '0');
                    chunk_len++;
                    continue;
                }
            }

            /* Any other non-space char terminates the number (or is an error). */
            if (isspace(c)) {
                /* End of the number */
                goto done_reading;
            } else {
                /* Invalid character in the number */
                /* errno = EINVAL; */
                CFX_PRINT_ERR("Invalid character found: '%c' (0x%x)!", c, c);
                return -1;
            }
        }
    }

done_reading:
    /* Flush any pending chunk */
    flush_chunk(out, detected_base, &chunk_val, &chunk_len);

    if (!saw_digit) {
        CFX_PRINT_ERR("didn't find any digit!");
        return -1;
    }

    /* TODO: track sign */
    if (negative) {
        /* reject negatives for now */
        /* errno = ERANGE; -> or implement signed handling */
        CFX_PRINT_ERR("negative number found");
        return -1;
    }

    return 0;
}

double cfx_big_log(const cfx_big_t *b, double base) {
    size_t l = b->n;
    cfx_limb_t hi = b->limb[l-1];
    double ln_base = log(base);
    double ln_hi = log(hi);
    double ln_B = CFX_LIMB_BITS * log(2.0);
    return ((l - 1) * ln_B + ln_hi) / ln_base;
}

int cfx_big_to_sci(const cfx_big_t *x, unsigned base,
    int sig_digits, char *out, size_t outsz) {

    if (!x || !out || outsz == 0 || base < 2) return 0;

    /* zero? */
    if (x->n == 0) {
        snprintf(out, outsz, "0");
        return 1;
    }

    /* find top nonzero limb (defensive) */
    size_t k = x->n;
    while (k > 0 && x->limb[k-1] == 0) --k;
    if (k == 0) {
        snprintf(out, outsz, "0"); return 1;
    }

    cfx_limb_t hi = x->limb[k-1];

    const long double lnB   = (long double)CFX_LIMB_BITS * logl(2.0L);
    const long double lnb   = logl((long double)base);

    /* ln(n) ≈ (k-1)*ln(B) + ln(hi) */
    long double ln_n = ((long double)(k - 1)) * lnB + logl((long double)hi);

    /* e = floor( ln(n) / ln(b) ) */
    long double logb_n = ln_n / lnb;
    long long e = (long long) floorl(logb_n);

    /* mantissa m = b^( fractional part ) */
    long double frac = logb_n - (long double)e;
    long double m = expl(frac * lnb);   /* == powl(base, frac);  guarantees 1 <= m < b (up to rounding) */

    /* rounding guard: if m rounds to exactly b, bump e and renormalize */
    if (!(m < (long double)base)) {
        m /= (long double)base; ++e;
    }

    /* format: decimal mantissa * base^e */
    if (sig_digits < 1) sig_digits = 1;
    int after_decimal = sig_digits - 1;

    if (base == 10) {
        /* 3.236e8930 */
        snprintf(out, outsz, "%.*Lf" "e%lld", after_decimal, m, e);
    } else {
        /* 3.236 * 7^12345 */
        snprintf(out, outsz, "%.*Lf x %u^%lld", after_decimal, m, base, e);
    }
    return 1;
}
