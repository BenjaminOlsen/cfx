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

/* swap two entries by 1-based index. returns new malloc'd buf, or NULL on error. */
uint8_t *store_swap(const uint8_t *pt, size_t pt_len,
                    unsigned idx_a, unsigned idx_b, size_t *new_len) {
    if (idx_a == idx_b || idx_a == 0 || idx_b == 0) {
        fprintf(stderr, "error: swap requires two distinct positive indices\n");
        return NULL;
    }

    /* ensure idx_a < idx_b for simpler reconstruction */
    if (idx_a > idx_b) { unsigned t = idx_a; idx_a = idx_b; idx_b = t; }

    /* scan to locate both entries */
    const uint8_t *p = pt;
    const uint8_t *end = pt + pt_len;
    unsigned cur = 0;

    const uint8_t *a_start = NULL, *a_end = NULL;
    const uint8_t *b_start = NULL, *b_end = NULL;

    while (p + 2 <= end) {
        const uint8_t *entry_start = p;
        uint16_t klen = cfx_load16_le(p);
        p += 2;
        if (p + klen > end) break;
        p += klen;
        if (p + 4 > end) break;
        uint32_t vl = cfx_load32_le(p);
        p += 4;
        if (p + vl > end) break;
        p += vl;

        ++cur;
        if (cur == idx_a) { a_start = entry_start; a_end = p; }
        if (cur == idx_b) { b_start = entry_start; b_end = p; break; }
    }

    if (!a_start || !b_start) {
        fprintf(stderr, "error: index out of range\n");
        return NULL;
    }

    /* rebuild: [prefix][B][middle][A][suffix] */
    size_t a_len = (size_t)(a_end - a_start);
    size_t b_len = (size_t)(b_end - b_start);
    size_t prefix_len = (size_t)(a_start - pt);
    size_t middle_len = (size_t)(b_start - a_end);
    size_t suffix_len = pt_len - (size_t)(b_end - pt);

    uint8_t *out = malloc(pt_len);
    if (!out) return NULL;

    uint8_t *w = out;
    memcpy(w, pt, prefix_len);               w += prefix_len;
    memcpy(w, b_start, b_len);               w += b_len;
    memcpy(w, a_end, middle_len);             w += middle_len;
    memcpy(w, a_start, a_len);               w += a_len;
    memcpy(w, b_end, suffix_len);             w += suffix_len;

    *new_len = pt_len;
    return out;
}

/* --- sort support --- */

typedef struct {
    const uint8_t *start;    /* pointer to entry start */
    size_t         len;      /* total entry length */
    const uint8_t *key;      /* pointer to key bytes */
    size_t         klen;     /* key length */
} store_entry_ref;

/* comparator of keys, case insensitive */
static int entry_cmp(const void *a, const void *b) {
    const store_entry_ref *ea = a;
    const store_entry_ref *eb = b;
    size_t minlen = ea->klen < eb->klen ? ea->klen : eb->klen;
    for (size_t i = 0; i < minlen; i++) {
        unsigned char ca = ea->key[i];
        unsigned char cb = eb->key[i];
        /* to lowercase */
        if (ca >= 'A' && ca <= 'Z') ca += 32;
        if (cb >= 'A' && cb <= 'Z') cb += 32;
        if (ca != cb) return (int)ca - (int)cb;
    }
    return (int)ea->klen - (int)eb->klen;
}

/* sort entries alphabetically returns new malloc'd buf. */
uint8_t *store_sort(const uint8_t *pt, size_t pt_len, size_t *new_len) {
    if (pt_len == 0) {
        *new_len = 0;
        return malloc(1);
    }

    /* count entries */
    unsigned count = 0;
    const uint8_t *p = pt;
    const uint8_t *end = pt + pt_len;
    while (p + 2 <= end) {
        size_t klen = cfx_load16_le(p); p += 2;
        if (p + klen > end) break;
        p += klen;
        if (p + 4 > end) break;
        uint32_t vl = cfx_load32_le(p); p += 4;
        if (p + vl > end) break;
        p += vl;
        count++;
    }

    if (count <= 1) {
        uint8_t *out = malloc(pt_len);
        if (!out) return NULL;
        memcpy(out, pt, pt_len);
        *new_len = pt_len;
        return out;
    }

    /* collect entry pointers */
    store_entry_ref *entries = malloc(count * sizeof(store_entry_ref));
    if (!entries) return NULL;

    p = pt;
    unsigned idx = 0;
    while (p + 2 <= end && idx < count) {
        const uint8_t *entry_start = p;
        size_t klen = cfx_load16_le(p); p += 2;
        if (p + klen > end) break;
        const uint8_t *kptr = p;
        p += klen;
        if (p + 4 > end) break;
        uint32_t vl = cfx_load32_le(p); p += 4;
        if (p + vl > end) break;
        p += vl;

        entries[idx].start = entry_start;
        entries[idx].len   = (size_t)(p - entry_start);
        entries[idx].key   = kptr;
        entries[idx].klen  = klen;
        idx++;
    }

    qsort(entries, idx, sizeof(store_entry_ref), entry_cmp);

    /* rebuild buffer in sorted order */
    uint8_t *out = malloc(pt_len);
    if (!out) { free(entries); return NULL; }

    uint8_t *w = out;
    for (unsigned i = 0; i < idx; i++) {
        memcpy(w, entries[i].start, entries[i].len);
        w += entries[i].len;
    }

    free(entries);
    *new_len = pt_len;
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
