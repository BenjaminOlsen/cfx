#include "cfx/base64.h"

static int cfx_is_ws(unsigned char c) {
    return (c == ' ' || c == '\t' || c == '\r' || c == '\n' || c == '\f' || c == '\v');
}

size_t cfx_base64_enc_len(size_t bin_len) {
    /* 4 * ceil(n/3) */
    return (bin_len / 3) * 4 + ((bin_len % 3) ? 4 : 0);
}

size_t cfx_base64_dec_max_len(size_t b64_len) {
    /* Upper bound assuming no whitespace and minimal padding checks: 3*(n/4) */
    return (b64_len / 4) * 3 + 3;
}

int cfx_base64_encode(char *out, size_t *out_len, const uint8_t *in, size_t in_len) {
    static const char tbl[65] =
        "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

    size_t need = cfx_base64_enc_len(in_len);
    if (out_len == NULL) return -1;

    if (out == NULL) {
        *out_len = need;
        return 0;
    }
    if (*out_len < need) return -1;

    size_t i = 0, o = 0;

    while (i + 3 <= in_len) {
        uint32_t v = ((uint32_t)in[i] << 16) |
                     ((uint32_t)in[i + 1] << 8) |
                     ((uint32_t)in[i + 2]);
        out[o + 0] = tbl[(v >> 18) & 0x3F];
        out[o + 1] = tbl[(v >> 12) & 0x3F];
        out[o + 2] = tbl[(v >>  6) & 0x3F];
        out[o + 3] = tbl[(v >>  0) & 0x3F];
        i += 3;
        o += 4;
    }

    /* tail */
    if (i < in_len) {
        uint32_t v = (uint32_t)in[i] << 16;
        out[o + 0] = tbl[(v >> 18) & 0x3F];

        if (i + 1 < in_len) {
            v |= (uint32_t)in[i + 1] << 8;
            out[o + 1] = tbl[(v >> 12) & 0x3F];
            out[o + 2] = tbl[(v >>  6) & 0x3F];
            out[o + 3] = '=';
        } else {
            out[o + 1] = tbl[(v >> 12) & 0x3F];
            out[o + 2] = '=';
            out[o + 3] = '=';
        }
        o += 4;
    }

    *out_len = o;
    return 0;
}

static int8_t cfx_b64_val(unsigned char c) {
    if (c >= 'A' && c <= 'Z') return (int8_t)(c - 'A');
    if (c >= 'a' && c <= 'z') return (int8_t)(c - 'a' + 26);
    if (c >= '0' && c <= '9') return (int8_t)(c - '0' + 52);
    if (c == '+') return 62;
    if (c == '/') return 63;
    return -1;
}

int cfx_base64_decode(uint8_t *out, size_t *out_len, const char *in, size_t in_len) {
    if (out_len == NULL) return -2;

    /* First pass: count non-ws chars and validate charset/padding placement roughly. */
    size_t n = 0;          /* non-ws chars */
    size_t eq = 0;         /* number of '=' seen (must be 0,1,2 and only at end) */
    int seen_pad = 0;

    for (size_t i = 0; i < in_len; i++) {
        unsigned char c = (unsigned char)in[i];
        if (cfx_is_ws(c)) continue;

        if (c == '=') {
            seen_pad = 1;
            eq++;
            n++;
            continue;
        }

        if (seen_pad) {
            /* After '=', only whitespace is allowed. */
            return -3;
        }

        if (cfx_b64_val(c) < 0) return -4;
        n++;
    }

    if (n == 0) {
        if (out == NULL) {
            *out_len = 0; return 0;
        }
        *out_len = 0;
        return 0;
    }

    if ((n % 4) != 0) return -5;
    if (eq > 2) return -6;

    /* Determine exact output length from padding count. */
    size_t quads = n / 4;
    size_t dec_len = quads * 3;
    if (eq) dec_len -= eq;

    if (out == NULL) {
        *out_len = dec_len;
        return 0;
    }
    if (*out_len < dec_len) {
        return -7;
    }

    /* Second pass: decode. */
    size_t o = 0;
    int8_t q[4];
    size_t qi = 0;
    int in_padding_quad = 0;

    for (size_t i = 0; i < in_len; i++) {
        unsigned char c = (unsigned char)in[i];
        if (cfx_is_ws(c)) continue;

        if (c == '=') {
            q[qi++] = -2; /* sentinel for pad */
        } else {
            q[qi++] = cfx_b64_val(c);
            if (q[qi - 1] < 0) return -8;
        }

        if (qi < 4) continue;
        qi = 0;

        /*
         * Padding rules (canonical):
         *  - '=' can appear only in positions 3 and/or 4 of the last quad.
         *  - If q[2] is pad, q[3] must be pad (==).
         *  - If q[3] is pad, q[2] must not be pad (single '=').
         *  - No non-pad quads after a padded quad.
         */
        int pad2 = (q[2] == -2);
        int pad3 = (q[3] == -2);

        if (in_padding_quad) {
            /* Once we've decoded a padded quad, there must be nothing else. */
            return -9;
        }

        if (!pad2 && !pad3) {
            uint32_t v = ((uint32_t)q[0] << 18) |
                         ((uint32_t)q[1] << 12) |
                         ((uint32_t)q[2] <<  6) |
                         ((uint32_t)q[3] <<  0);
            out[o++] = (uint8_t)((v >> 16) & 0xFF);
            out[o++] = (uint8_t)((v >>  8) & 0xFF);
            out[o++] = (uint8_t)((v >>  0) & 0xFF);
            continue;
        }

        /* padded quad: must be last */
        in_padding_quad = 1;

        if (pad2 && !pad3) return -10;          /* "= in 3rd but not 4th" impossible */
        if (q[0] < 0 || q[1] < 0) return -11;   /* first two must be real */

        if (pad2 && pad3) {
            /* 'xx==' -> 1 byte output; enforce canonical unused bits: low 4 bits of q1 must be 0 */
            if ((q[1] & 0x0F) != 0) return -12;
            uint32_t v = ((uint32_t)q[0] << 18) | ((uint32_t)q[1] << 12);
            out[o++] = (uint8_t)((v >> 16) & 0xFF);
        } else {
            /* 'xxx=' -> 2 bytes output; enforce canonical unused bits: low 2 bits of q2 must be 0 */
            if (q[2] < 0) return -13;
            if ((q[2] & 0x03) != 0) return -14;
            uint32_t v = ((uint32_t)q[0] << 18) |
                         ((uint32_t)q[1] << 12) |
                         ((uint32_t)q[2] <<  6);
            out[o++] = (uint8_t)((v >> 16) & 0xFF);
            out[o++] = (uint8_t)((v >>  8) & 0xFF);
        }
    }

    /* If we ended mid-quad, it's invalid (but should be caught by n%4). */
    if (qi != 0) return -15;

    *out_len = o;
    return 0;
}
