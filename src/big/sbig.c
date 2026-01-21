/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#include "cfx/sbig.h"
#include <stdlib.h>
#include <string.h>

void cfx_sbig_init(cfx_sbig_t* s) {
    cfx_big_init(&s->mag);
    s->sign = 0;
}

void cfx_sbig_free(cfx_sbig_t* s) {
    cfx_big_free(&s->mag);
    s->sign = 0;
}

void cfx_sbig_clear(cfx_sbig_t* s) {
    cfx_big_clear(&s->mag);
    s->sign = 0;
}

void cfx_sbig_normalize(cfx_sbig_t* s) {
    if (cfx_big_is_zero(&s->mag)) {
        s->sign = 0;
    }
}

void cfx_sbig_assign(cfx_sbig_t* dst, const cfx_sbig_t* src) {
    if (dst == src) return;
    cfx_big_assign(&dst->mag, &src->mag);
    dst->sign = src->sign;
}

void cfx_sbig_assign_big(cfx_sbig_t* dst, const cfx_big_t* mag, int8_t sign) {
    cfx_big_assign(&dst->mag, mag);
    if (cfx_big_is_zero(mag)) {
        dst->sign = 0;
    } else {
        dst->sign = (sign < 0) ? -1 : 1;
    }
}

void cfx_sbig_from_i64(cfx_sbig_t* s, int64_t val) {
    if (val == 0) {
        cfx_big_from_u64(&s->mag, 0);
        s->sign = 0;
    } else if (val > 0) {
        cfx_big_from_u64(&s->mag, (uint64_t)val);
        s->sign = 1;
    } else {
        /* INT64_MIN can't be negated directly */
        if (val == INT64_MIN) {
            cfx_big_from_u64(&s->mag, (uint64_t)INT64_MAX + 1);
        } else {
            cfx_big_from_u64(&s->mag, (uint64_t)(-val));
        }
        s->sign = -1;
    }
}

void cfx_sbig_from_u64(cfx_sbig_t* s, uint64_t val) {
    cfx_big_from_u64(&s->mag, val);
    s->sign = (val == 0) ? 0 : 1;
}

void cfx_sbig_neg_eq(cfx_sbig_t* s) {
    s->sign = -s->sign;
}

void cfx_sbig_abs_eq(cfx_sbig_t* s) {
    if (s->sign < 0) s->sign = 1;
}

void cfx_sbig_neg(cfx_sbig_t* out, const cfx_sbig_t* a) {
    cfx_sbig_assign(out, a);
    cfx_sbig_neg_eq(out);
}

void cfx_sbig_abs(cfx_sbig_t* out, const cfx_sbig_t* a) {
    cfx_sbig_assign(out, a);
    cfx_sbig_abs_eq(out);
}

int cfx_sbig_cmp(const cfx_sbig_t* a, const cfx_sbig_t* b) {
    if (a->sign != b->sign) {
        return (a->sign < b->sign) ? -1 : 1;
    }
    if (a->sign == 0) {
        return 0;
    }
    int mag_cmp = cfx_big_cmp(&a->mag, &b->mag);
    /* negative: larger magnitude = smaller value */
    return (a->sign > 0) ? mag_cmp : -mag_cmp;
}

void cfx_sbig_swap(cfx_sbig_t* a, cfx_sbig_t* b) {
    cfx_big_swap(&a->mag, &b->mag);
    int8_t tmp = a->sign;
    a->sign = b->sign;
    b->sign = tmp;
}

void cfx_sbig_add(cfx_sbig_t* out, const cfx_sbig_t* a, const cfx_sbig_t* b) {
    if (a->sign == 0) {
        cfx_sbig_assign(out, b);
        return;
    }
    if (b->sign == 0) {
        cfx_sbig_assign(out, a);
        return;
    }

    if (a->sign == b->sign) {
        cfx_big_add(&out->mag, &a->mag, &b->mag);
        out->sign = a->sign;
        return;
    }

    /* different signs: subtract magnitudes */
    int cmp = cfx_big_cmp(&a->mag, &b->mag);
    if (cmp == 0) {
        cfx_sbig_clear(out);
        return;
    }
    if (cmp > 0) {
        cfx_big_sub(&out->mag, &a->mag, &b->mag);
        out->sign = a->sign;
    } else {
        cfx_big_sub(&out->mag, &b->mag, &a->mag);
        out->sign = b->sign;
    }

    cfx_sbig_normalize(out);
}

void cfx_sbig_sub(cfx_sbig_t* out, const cfx_sbig_t* a, const cfx_sbig_t* b) {
    cfx_sbig_t neg_b;
    cfx_sbig_init(&neg_b);
    cfx_sbig_neg(&neg_b, b);
    cfx_sbig_add(out, a, &neg_b);
    cfx_sbig_free(&neg_b);
}

void cfx_sbig_mul(cfx_sbig_t* out, const cfx_sbig_t* a, const cfx_sbig_t* b) {
    if (a->sign == 0 || b->sign == 0) {
        cfx_sbig_clear(out);
        return;
    }
    cfx_big_mul(&out->mag, &a->mag, &b->mag);
    out->sign = (a->sign == b->sign) ? 1 : -1;
    cfx_sbig_normalize(out);
}

void cfx_sbig_divrem(cfx_sbig_t* q, cfx_sbig_t* rem,
                     const cfx_sbig_t* a, const cfx_sbig_t* b) {
    if (b->sign == 0) {
        if (q) cfx_sbig_clear(q);
        if (rem) cfx_sbig_clear(rem);
        return;
    }
    if (a->sign == 0) {
        if (q) cfx_sbig_clear(q);
        if (rem) cfx_sbig_clear(rem);
        return;
    }

    cfx_big_t q_mag, r_mag;
    cfx_big_init(&q_mag);
    cfx_big_init(&r_mag);
    cfx_big_divrem(&q_mag, &r_mag, &a->mag, &b->mag);

    if (q) {
        cfx_big_assign(&q->mag, &q_mag);
        q->sign = cfx_big_is_zero(&q_mag) ? 0 : ((a->sign == b->sign) ? 1 : -1);
    }

    /* remainder has same sign as dividend (C99) */
    if (rem) {
        cfx_big_assign(&rem->mag, &r_mag);
        if (cfx_big_is_zero(&r_mag)) {
            rem->sign = 0;
        } else {
            rem->sign = a->sign;
        }
    }

    cfx_big_free(&q_mag);
    cfx_big_free(&r_mag);
}

void cfx_sbig_add_eq(cfx_sbig_t* a, const cfx_sbig_t* b) {
    cfx_sbig_t tmp;
    cfx_sbig_init(&tmp);
    cfx_sbig_add(&tmp, a, b);
    cfx_sbig_swap(a, &tmp);
    cfx_sbig_free(&tmp);
}

void cfx_sbig_sub_eq(cfx_sbig_t* a, const cfx_sbig_t* b) {
    cfx_sbig_t tmp;
    cfx_sbig_init(&tmp);
    cfx_sbig_sub(&tmp, a, b);
    cfx_sbig_swap(a, &tmp);
    cfx_sbig_free(&tmp);
}

void cfx_sbig_mul_eq(cfx_sbig_t* a, const cfx_sbig_t* b) {
    cfx_sbig_t tmp;
    cfx_sbig_init(&tmp);
    cfx_sbig_mul(&tmp, a, b);
    cfx_sbig_swap(a, &tmp);
    cfx_sbig_free(&tmp);
}

char* cfx_sbig_dec_alloc(const cfx_sbig_t* s, size_t* len_out) {
    if (s->sign == 0) {
        char* str = (char*)malloc(2);
        if (!str) return NULL;
        str[0] = '0';
        str[1] = '\0';
        if (len_out) *len_out = 1;
        return str;
    }

    size_t mag_len;
    char* mag_str = cfx_big_dec_alloc(&s->mag, &mag_len);
    if (!mag_str) return NULL;

    if (s->sign > 0) {
        if (len_out) *len_out = mag_len;
        return mag_str;
    }

    char* str = (char*)malloc(mag_len + 2);
    if (!str) {
        free(mag_str);
        return NULL;
    }
    str[0] = '-';
    memcpy(str + 1, mag_str, mag_len + 1);
    free(mag_str);

    if (len_out) *len_out = mag_len + 1;
    return str;
}

int cfx_sbig_from_dec(cfx_sbig_t* out, const char* str) {
    if (!str || !*str) {
        cfx_sbig_clear(out);
        return -1;
    }

    int8_t sign = 1;
    if (*str == '-') {
        sign = -1;
        str++;
    } else if (*str == '+') {
        str++;
    }

    if (!*str) {
        cfx_sbig_clear(out);
        return -1;
    }

    int ret = cfx_big_from_dec(&out->mag, str);
    if (ret != 0) {
        cfx_sbig_clear(out);
        return ret;
    }

    if (cfx_big_is_zero(&out->mag)) {
        out->sign = 0;
    } else {
        out->sign = sign;
    }
    return 0;
}

int cfx_sbig_from_hex(cfx_sbig_t* out, const char* str) {
    if (!str || !*str) {
        cfx_sbig_clear(out);
        return -1;
    }

    int8_t sign = 1;
    if (*str == '-') {
        sign = -1;
        str++;
    } else if (*str == '+') {
        str++;
    }

    if (!*str) {
        cfx_sbig_clear(out);
        return -1;
    }

    int ret = cfx_big_from_hex(&out->mag, str);
    if (ret != 0) {
        cfx_sbig_clear(out);
        return ret;
    }

    if (cfx_big_is_zero(&out->mag)) {
        out->sign = 0;
    } else {
        out->sign = sign;
    }
    return 0;
}

void cfx_sbig_xgcd(cfx_sbig_t* g, cfx_sbig_t* x, cfx_sbig_t* y,
                   const cfx_sbig_t* a, const cfx_sbig_t* b) {
    /* gcd(a, b) = gcd(|a|, |b|)
     * If |a|*x' + |b|*y' = g, then a*x + b*y = g where:
     *   x = x' * sign(a)  (or x' if a >= 0)
     *   y = y' * sign(b)  (or y' if b >= 0)
     */
    cfx_sbig_t x_tmp, y_tmp;
    cfx_sbig_init(&x_tmp);
    cfx_sbig_init(&y_tmp);

    cfx_big_t g_tmp;
    cfx_big_init(&g_tmp);

    /* Call the unsigned version with magnitudes */
    cfx_big_xgcd(&g_tmp, &x_tmp, &y_tmp, &a->mag, &b->mag);

    /* Adjust signs of coefficients based on input signs */
    if (a->sign < 0) {
        cfx_sbig_neg_eq(&x_tmp);
    }
    if (b->sign < 0) {
        cfx_sbig_neg_eq(&y_tmp);
    }

    /* Output results */
    if (g) {
        cfx_sbig_assign_big(g, &g_tmp, 1);
    }
    if (x) {
        cfx_sbig_swap(x, &x_tmp);
    }
    if (y) {
        cfx_sbig_swap(y, &y_tmp);
    }

    cfx_big_free(&g_tmp);
    cfx_sbig_free(&x_tmp);
    cfx_sbig_free(&y_tmp);
}
