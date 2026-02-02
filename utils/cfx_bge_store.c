/* cfx_bge_store.c -- KV store manipulation for BGE */

#include "cfx_bge_internal.h"

/* resolve 1-based index to entry name. NULL if out of range. */
const uint8_t *store_name_by_index(const uint8_t *pt, size_t pt_len,
                                   unsigned idx, uint16_t *name_len) {
    const uint8_t *p = pt;
    const uint8_t *end = pt + pt_len;
    unsigned cur = 0;

    while (p + 2 <= end) {
        uint16_t klen = cfx_load16_le(p);
        p += 2;
        if (p + klen > end) break;
        const uint8_t *kptr = p;
        p += klen;
        if (p + 4 > end) break;
        uint32_t vl = cfx_load32_le(p);
        p += 4;
        if (p + vl > end) break;
        p += vl;

        if (++cur == idx) {
            *name_len = klen;
            return kptr;
        }
    }
    return NULL;
}

/* find entry by name, return ptr to value + length. NULL if missing.
   vlen can be null if you don't care about the length */
const uint8_t *store_get(const uint8_t *pt, size_t pt_len,
                         const char *name, size_t *vlen) {
    size_t nlen = strlen(name);
    const uint8_t *p = pt;
    const uint8_t *end = pt + pt_len;

    while (p + 2 <= end) {
        uint16_t klen = cfx_load16_le(p);
        p += 2;
        if (p + klen > end) break;
        const uint8_t *kptr = p;
        p += klen;

        if (p + 4 > end) break;
        uint32_t val_len = cfx_load32_le(p);
        p += 4;
        if (p + val_len > end) break;
        const uint8_t *vptr = p;
        p += val_len;

        if (klen == nlen && memcmp(kptr, name, nlen) == 0) {
            if (vlen) *vlen = val_len;
            return vptr;
        }
    }
    return NULL;
}

/* upsert entry. returns new malloc'd store buf, sets *new_len. */
uint8_t *store_set(const uint8_t *pt, size_t pt_len, const char *name,
                   const uint8_t *val, size_t val_len, size_t *new_len) {
    size_t nlen = strlen(name);

    /* size of new entry: 2 + nlen + 4 + val_len */
    size_t entry_len = 2 + nlen + 4 + val_len;

    /* scan for existing entry */
    const uint8_t *p = pt;
    const uint8_t *end = pt + pt_len;
    const uint8_t *found_start = NULL;
    size_t found_entry_len = 0;

    while (p + 2 <= end) {
        const uint8_t *entry_start = p;
        uint16_t klen = cfx_load16_le(p);
        p += 2;
        if (p + klen > end) break;
        const uint8_t *kptr = p;
        p += klen;

        if (p + 4 > end) break;
        uint32_t vl = cfx_load32_le(p);
        p += 4;
        if (p + vl > end) break;
        p += vl;

        if (klen == nlen && memcmp(kptr, name, nlen) == 0) {
            found_start = entry_start;
            found_entry_len = (size_t)(p - entry_start);
            break;
        }
    }

    size_t out_len;
    uint8_t *out;

    if (found_start) {
        /* replace in-place */
        size_t prefix = (size_t)(found_start - pt);
        size_t suffix_off = prefix + found_entry_len;
        size_t suffix_len = pt_len - suffix_off;
        out_len = prefix + entry_len + suffix_len;
        out = malloc(out_len);
        if (!out) return NULL;

        memcpy(out, pt, prefix);
        uint8_t *w = out + prefix;
        cfx_store16_le(w, (uint16_t)nlen); w += 2;
        memcpy(w, name, nlen); w += nlen;
        cfx_store32_le(w, (uint32_t)val_len); w += 4;
        memcpy(w, val, val_len); w += val_len;
        memcpy(w, pt + suffix_off, suffix_len);
    } else {
        /* append */
        out_len = pt_len + entry_len;
        out = malloc(out_len);
        if (!out) return NULL;

        memcpy(out, pt, pt_len);
        uint8_t *w = out + pt_len;
        cfx_store16_le(w, (uint16_t)nlen); w += 2;
        memcpy(w, name, nlen); w += nlen;
        cfx_store32_le(w, (uint32_t)val_len); w += 4;
        memcpy(w, val, val_len);
    }

    *new_len = out_len;
    return out;
}

/* remove entry by name. returns new malloc'd buf or NULL if not found. */
uint8_t *store_rm(const uint8_t *pt, size_t pt_len, const char *name, size_t *new_len) {
    size_t nlen = strlen(name);
    const uint8_t *p = pt;
    const uint8_t *end = pt + pt_len;
    const uint8_t *found_start = NULL;
    size_t found_entry_len = 0;

    while (p + 2 <= end) {
        const uint8_t *entry_start = p;
        uint16_t klen = cfx_load16_le(p);
        p += 2;
        if (p + klen > end) break;
        const uint8_t *kptr = p;
        p += klen;

        if (p + 4 > end) break;
        uint32_t vl = cfx_load32_le(p);
        p += 4;
        if (p + vl > end) break;
        p += vl;

        if (klen == nlen && memcmp(kptr, name, nlen) == 0) {
            found_start = entry_start;
            found_entry_len = (size_t)(p - entry_start);
            break;
        }
    }

    if (!found_start) {
        fprintf(stderr, "error: '%s' not found\n", name);
        return NULL;
    }

    size_t prefix = (size_t)(found_start - pt);
    size_t suffix_off = prefix + found_entry_len;
    size_t suffix_len = pt_len - suffix_off;
    size_t out_len = prefix + suffix_len;

    uint8_t *out = malloc(out_len + 1);
    if (!out) return NULL;

    memcpy(out, pt, prefix);
    memcpy(out + prefix, pt + suffix_off, suffix_len);
    *new_len = out_len;
    return out;
}

/* for dump - render as "[name]\nvalue\n" pairs. malloc'd, caller frees. */
uint8_t *store_to_text(const uint8_t *pt, size_t pt_len, size_t *text_len) {
    /* pass 1: size */
    size_t total = 0;
    const uint8_t *p = pt;
    const uint8_t *end = pt + pt_len;

    while (p + 2 <= end) {
        uint16_t klen = cfx_load16_le(p);
        p += 2;
        if (p + klen > end) break;
        p += klen;
        if (p + 4 > end) break;
        uint32_t vl = cfx_load32_le(p);
        p += 4;
        if (p + vl > end) break;
        p += vl;

        /* [name]\n + value + \n */
        total += 1 + klen + 2 + vl + 1;
    }

    uint8_t *out = malloc(total + 1);
    if (!out) return NULL;

    /* pass 2: write */
    uint8_t *w = out;
    p = pt;
    while (p + 2 <= end) {
        uint16_t klen = cfx_load16_le(p);
        p += 2;
        if (p + klen > end) break;
        const uint8_t *kptr = p;
        p += klen;
        if (p + 4 > end) break;
        uint32_t vl = cfx_load32_le(p);
        p += 4;
        if (p + vl > end) break;
        const uint8_t *vptr = p;
        p += vl;

        *w++ = '[';
        memcpy(w, kptr, klen); w += klen;
        *w++ = ']';
        *w++ = '\n';
        memcpy(w, vptr, vl); w += vl;
        *w++ = '\n';
    }

    *w = '\0';
    *text_len = (size_t)(w - out);
    return out;
}
