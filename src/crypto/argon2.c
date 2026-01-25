/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

/* argon2.c - Argon2 password hashing (RFC 9106) */

#include "cfx/argon2.h"
#include "cfx/blake2.h"
#include "cfx/base64.h"
#include "cfx/memory.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define ARGON2_BLOCK_SIZE 1024
#define ARGON2_QWORDS_IN_BLOCK 128
#define ARGON2_SYNC_POINTS 4

typedef struct {
    uint64_t v[ARGON2_QWORDS_IN_BLOCK];
} block_t;

typedef struct {
    block_t *memory;
    uint32_t m_cost;
    uint32_t t_cost;
    uint32_t lanes;
    uint32_t lane_length;
    uint32_t segment_length;
    int type;
} argon2_ctx_t;

static inline uint64_t rotr64(uint64_t x, unsigned n) {
    return (x >> n) | (x << (64 - n));
}

/* H' variable length hash */
static void hash_long(uint8_t *out, size_t outlen, const uint8_t *in, size_t inlen) {
    uint8_t outlen_bytes[4];
    cfx_store32_le(outlen_bytes, (uint32_t)outlen);

    if (outlen <= 64) {
        cfx_blake2b_ctx_t ctx;
        cfx_blake2b_init(&ctx, outlen);
        cfx_blake2b_update(&ctx, outlen_bytes, 4);
        cfx_blake2b_update(&ctx, in, inlen);
        cfx_blake2b_final(&ctx, out);
    } else {
        uint8_t v[64];
        cfx_blake2b_ctx_t ctx;

        cfx_blake2b_init(&ctx, 64);
        cfx_blake2b_update(&ctx, outlen_bytes, 4);
        cfx_blake2b_update(&ctx, in, inlen);
        cfx_blake2b_final(&ctx, v);

        memcpy(out, v, 32);
        out += 32;
        outlen -= 32;

        while (outlen > 64) {
            cfx_blake2b_init(&ctx, 64);
            cfx_blake2b_update(&ctx, v, 64);
            cfx_blake2b_final(&ctx, v);
            memcpy(out, v, 32);
            out += 32;
            outlen -= 32;
        }

        size_t final_len = (outlen <= 64) ? outlen : 64;
        cfx_blake2b_init(&ctx, final_len);
        cfx_blake2b_update(&ctx, v, 64);
        cfx_blake2b_final(&ctx, out);
    }
}

/* blake2b mixing function G (modified for argon2) */
#define G(a, b, c, d) do { \
            a = a + b + 2 * (uint64_t)(uint32_t)a * (uint64_t)(uint32_t)b; \
            d = rotr64(d ^ a, 32); \
            c = c + d + 2 * (uint64_t)(uint32_t)c * (uint64_t)(uint32_t)d; \
            b = rotr64(b ^ c, 24); \
            a = a + b + 2 * (uint64_t)(uint32_t)a * (uint64_t)(uint32_t)b; \
            d = rotr64(d ^ a, 16); \
            c = c + d + 2 * (uint64_t)(uint32_t)c * (uint64_t)(uint32_t)d; \
            b = rotr64(b ^ c, 63); \
} while (0)

static void blake2b_round(uint64_t *v) {
    G(v[0], v[4], v[8],  v[12]);
    G(v[1], v[5], v[9],  v[13]);
    G(v[2], v[6], v[10], v[14]);
    G(v[3], v[7], v[11], v[15]);
    G(v[0], v[5], v[10], v[15]);
    G(v[1], v[6], v[11], v[12]);
    G(v[2], v[7], v[8],  v[13]);
    G(v[3], v[4], v[9],  v[14]);
}

static void P(uint64_t *v) {
    blake2b_round(v);
    blake2b_round(v);
}

static void fill_block(const block_t *prev, const block_t *ref, block_t *next, int with_xor) {
    block_t tmp;
    uint64_t r[ARGON2_QWORDS_IN_BLOCK];

    for (int i = 0; i < ARGON2_QWORDS_IN_BLOCK; i++) {
        r[i] = prev->v[i] ^ ref->v[i];
    }
    memcpy(&tmp, r, sizeof(tmp));

    for (int i = 0; i < 8; i++) {
        P(&tmp.v[i * 16]);
    }

    for (int i = 0; i < 8; i++) {
        uint64_t col[16];
        for (int j = 0; j < 16; j++) {
            col[j] = tmp.v[j * 8 + i];
        }
        P(col);
        for (int j = 0; j < 16; j++) {
            tmp.v[j * 8 + i] = col[j];
        }
    }

    if (with_xor) {
        for (int i = 0; i < ARGON2_QWORDS_IN_BLOCK; i++) {
            next->v[i] ^= r[i] ^ tmp.v[i];
        }
    } else {
        for (int i = 0; i < ARGON2_QWORDS_IN_BLOCK; i++) {
            next->v[i] = r[i] ^ tmp.v[i];
        }
    }
}

static uint32_t index_alpha(const argon2_ctx_t *ctx, uint32_t pass, uint32_t slice,
    uint32_t index, uint64_t pseudo_rand, int same_lane) {
    uint32_t ref_area_size;
    uint32_t start_position;

    if (pass == 0) {
        if (slice == 0) {
            ref_area_size = index - 1;
        } else {
            ref_area_size = slice * ctx->segment_length + index - 1;
            if (same_lane) {
                ref_area_size += ctx->segment_length;
            }
        }
        start_position = 0;
    } else {
        ref_area_size = ctx->lane_length - ctx->segment_length + index - 1;
        if (same_lane) {
            ref_area_size += ctx->segment_length;
        }
        start_position = (slice + 1) * ctx->segment_length;
        if (start_position >= ctx->lane_length) {
            start_position = 0;
        }
    }

    uint64_t rel_pos = pseudo_rand & 0xFFFFFFFF;
    rel_pos = (rel_pos * rel_pos) >> 32;
    rel_pos = ref_area_size - 1 - ((ref_area_size * rel_pos) >> 32);

    return (start_position + (uint32_t)rel_pos) % ctx->lane_length;
}

static void fill_segment(argon2_ctx_t *ctx, uint32_t pass, uint32_t lane, uint32_t slice) {
    uint32_t start_index = (pass == 0 && slice == 0) ? 2 : 0;
    uint32_t curr_offset = lane * ctx->lane_length + slice * ctx->segment_length + start_index;
    uint32_t prev_offset = curr_offset - 1;

    if (curr_offset % ctx->lane_length == 0) {
        prev_offset += ctx->lane_length;
    }

    int data_independent = (ctx->type == CFX_ARGON2I) ||
                           (ctx->type == CFX_ARGON2ID && pass == 0 && slice < ARGON2_SYNC_POINTS / 2);

    block_t addr_block, zero_block, input_block;
    memset(&zero_block, 0, sizeof(zero_block));
    memset(&input_block, 0, sizeof(input_block));

    if (data_independent) {
        input_block.v[0] = pass;
        input_block.v[1] = lane;
        input_block.v[2] = slice;
        input_block.v[3] = ctx->lane_length * ctx->lanes;
        input_block.v[4] = ctx->t_cost;
        input_block.v[5] = ctx->type;
    }

    uint32_t addr_index = 0;

    for (uint32_t i = start_index; i < ctx->segment_length; i++, curr_offset++, prev_offset++) {
        if (curr_offset % ctx->lane_length == 1) {
            prev_offset = curr_offset - 1;
        }

        uint64_t pseudo_rand;
        if (data_independent) {
            if (addr_index == 0) {
                input_block.v[6]++;
                fill_block(&zero_block, &input_block, &addr_block, 0);
                fill_block(&zero_block, &addr_block, &addr_block, 0);
            }
            pseudo_rand = addr_block.v[addr_index];
            addr_index = (addr_index + 1) % ARGON2_QWORDS_IN_BLOCK;
        } else {
            pseudo_rand = ctx->memory[prev_offset].v[0];
        }

        uint32_t ref_lane = (uint32_t)(pseudo_rand >> 32) % ctx->lanes;
        if (pass == 0 && slice == 0) {
            ref_lane = lane;
        }

        int same_lane = (ref_lane == lane);
        uint32_t ref_index = index_alpha(ctx, pass, slice, i, pseudo_rand, same_lane);
        uint32_t ref_offset = ref_lane * ctx->lane_length + ref_index;

        fill_block(&ctx->memory[prev_offset], &ctx->memory[ref_offset],
            &ctx->memory[curr_offset], pass > 0);
    }
}

static void initial_hash(uint8_t *h0, const uint8_t *pwd, size_t pwdlen,
    const uint8_t *salt, size_t saltlen,
    uint32_t m_cost, uint32_t t_cost, uint32_t p,
    size_t outlen, int type) {
    cfx_blake2b_ctx_t ctx;
    uint8_t buf[4];

    cfx_blake2b_init(&ctx, 64);

    cfx_store32_le(buf, p);
    cfx_blake2b_update(&ctx, buf, 4);

    cfx_store32_le(buf, (uint32_t)outlen);
    cfx_blake2b_update(&ctx, buf, 4);

    cfx_store32_le(buf, m_cost);
    cfx_blake2b_update(&ctx, buf, 4);

    cfx_store32_le(buf, t_cost);
    cfx_blake2b_update(&ctx, buf, 4);

    cfx_store32_le(buf, CFX_ARGON2_VERSION);
    cfx_blake2b_update(&ctx, buf, 4);

    cfx_store32_le(buf, (uint32_t)type);
    cfx_blake2b_update(&ctx, buf, 4);

    cfx_store32_le(buf, (uint32_t)pwdlen);
    cfx_blake2b_update(&ctx, buf, 4);
    if (pwdlen > 0) {
        cfx_blake2b_update(&ctx, pwd, pwdlen);
    }

    cfx_store32_le(buf, (uint32_t)saltlen);
    cfx_blake2b_update(&ctx, buf, 4);
    cfx_blake2b_update(&ctx, salt, saltlen);

    /* secret and AD not supported, length = 0 */
    cfx_store32_le(buf, 0);
    cfx_blake2b_update(&ctx, buf, 4);
    cfx_blake2b_update(&ctx, buf, 4);

    cfx_blake2b_final(&ctx, h0);
}

int cfx_argon2(uint8_t *out, size_t outlen,
    const uint8_t *pwd, size_t pwdlen,
    const uint8_t *salt, size_t saltlen,
    uint32_t m_cost, uint32_t t_cost, uint32_t p,
    int type) {
    if (out == NULL) return CFX_ARGON2_ERR_OUTPUT;
    if (outlen < CFX_ARGON2_MIN_OUTLEN) return CFX_ARGON2_ERR_OUTPUT;
    if (salt == NULL || saltlen < CFX_ARGON2_MIN_SALT_LEN) return CFX_ARGON2_ERR_SALT;
    if (t_cost < CFX_ARGON2_MIN_TIME) return CFX_ARGON2_ERR_PARAM;
    if (p < CFX_ARGON2_MIN_LANES) return CFX_ARGON2_ERR_PARAM;
    if (type < CFX_ARGON2D || type > CFX_ARGON2ID) return CFX_ARGON2_ERR_PARAM;

    uint32_t memory_blocks = m_cost;
    if (memory_blocks < 8 * p) {
        memory_blocks = 8 * p;
    }
    uint32_t segment_length = memory_blocks / (4 * p);
    memory_blocks = segment_length * 4 * p;

    argon2_ctx_t ctx;
    ctx.memory = (block_t *)calloc(memory_blocks, sizeof(block_t));
    if (ctx.memory == NULL) {
        return CFX_ARGON2_ERR_MEMORY;
    }

    ctx.m_cost = m_cost;
    ctx.t_cost = t_cost;
    ctx.lanes = p;
    ctx.lane_length = memory_blocks / p;
    ctx.segment_length = segment_length;
    ctx.type = type;

    uint8_t h0[64];
    initial_hash(h0, pwd, pwdlen, salt, saltlen, m_cost, t_cost, p, outlen, type);

    uint8_t block_hash_input[72];
    memcpy(block_hash_input, h0, 64);

    for (uint32_t lane = 0; lane < p; lane++) {
        cfx_store32_le(block_hash_input + 64, 0);
        cfx_store32_le(block_hash_input + 68, lane);
        hash_long((uint8_t *)ctx.memory[lane * ctx.lane_length].v, ARGON2_BLOCK_SIZE,
            block_hash_input, 72);

        cfx_store32_le(block_hash_input + 64, 1);
        hash_long((uint8_t *)ctx.memory[lane * ctx.lane_length + 1].v, ARGON2_BLOCK_SIZE,
            block_hash_input, 72);
    }

    for (uint32_t pass = 0; pass < t_cost; pass++) {
        for (uint32_t slice = 0; slice < ARGON2_SYNC_POINTS; slice++) {
            for (uint32_t lane = 0; lane < p; lane++) {
                fill_segment(&ctx, pass, lane, slice);
            }
        }
    }

    block_t final_block;
    memcpy(&final_block, &ctx.memory[ctx.lane_length - 1], sizeof(block_t));
    for (uint32_t lane = 1; lane < p; lane++) {
        uint32_t last_idx = lane * ctx.lane_length + ctx.lane_length - 1;
        for (int i = 0; i < ARGON2_QWORDS_IN_BLOCK; i++) {
            final_block.v[i] ^= ctx.memory[last_idx].v[i];
        }
    }

    hash_long(out, outlen, (uint8_t *)final_block.v, ARGON2_BLOCK_SIZE);

    cfx_memzero_s(ctx.memory, memory_blocks * sizeof(block_t));
    cfx_memzero_s(h0, sizeof(h0));
    cfx_memzero_s(&final_block, sizeof(final_block));
    free(ctx.memory);

    return CFX_ARGON2_OK;
}

int cfx_argon2id(uint8_t *out, size_t outlen,
    const uint8_t *pwd, size_t pwdlen,
    const uint8_t *salt, size_t saltlen,
    uint32_t m_cost, uint32_t t_cost, uint32_t p) {
    return cfx_argon2(out, outlen, pwd, pwdlen, salt, saltlen,
        m_cost, t_cost, p, CFX_ARGON2ID);
}

int cfx_argon2d(uint8_t *out, size_t outlen,
    const uint8_t *pwd, size_t pwdlen,
    const uint8_t *salt, size_t saltlen,
    uint32_t m_cost, uint32_t t_cost, uint32_t p) {
    return cfx_argon2(out, outlen, pwd, pwdlen, salt, saltlen,
        m_cost, t_cost, p, CFX_ARGON2D);
}

int cfx_argon2i(uint8_t *out, size_t outlen,
    const uint8_t *pwd, size_t pwdlen,
    const uint8_t *salt, size_t saltlen,
    uint32_t m_cost, uint32_t t_cost, uint32_t p) {
    return cfx_argon2(out, outlen, pwd, pwdlen, salt, saltlen,
        m_cost, t_cost, p, CFX_ARGON2I);
}

/* PHC string encoding */

static const char * type_name(int type) {
    switch (type) {
    case CFX_ARGON2D:  return "argon2d";
    case CFX_ARGON2I:  return "argon2i";
    case CFX_ARGON2ID: return "argon2id";
    default: return "argon2";
    }
}

int cfx_argon2_encode(char *out, size_t outlen,
    const uint8_t *hash, size_t hashlen,
    const uint8_t *salt, size_t saltlen,
    uint32_t m_cost, uint32_t t_cost, uint32_t p,
    int type) {
    size_t salt_b64_len = cfx_base64_enc_len(saltlen);
    size_t hash_b64_len = cfx_base64_enc_len(hashlen);
    size_t needed = 1 + strlen(type_name(type)) + 1 + 5 + 3 +
                    2 + 10 + 1 + 2 + 10 + 1 + 2 + 10 + 1 +
                    salt_b64_len + 1 + hash_b64_len + 1;

    if (outlen < needed) {
        return -1;
    }

    char salt_b64[256], hash_b64[256];
    size_t slen = sizeof(salt_b64), hlen = sizeof(hash_b64);
    cfx_base64_encode(salt_b64, &slen, salt, saltlen);
    cfx_base64_encode(hash_b64, &hlen, hash, hashlen);

    /* remove padding (PHC uses unpadded base64) */
    while (slen > 0 && salt_b64[slen - 1] == '=') slen--;
    while (hlen > 0 && hash_b64[hlen - 1] == '=') hlen--;
    salt_b64[slen] = '\0';
    hash_b64[hlen] = '\0';

    int len = snprintf(out, outlen, "$%s$v=%d$m=%u,t=%u,p=%u$%s$%s",
        type_name(type), CFX_ARGON2_VERSION,
        m_cost, t_cost, p, salt_b64, hash_b64);

    return len;
}

static int b64_decode_unpadded(uint8_t *out, size_t *outlen, const char *in, size_t inlen) {
    char padded[512];
    if (inlen >= sizeof(padded) - 4) return -1;

    memcpy(padded, in, inlen);
    size_t pad = (4 - (inlen % 4)) % 4;
    for (size_t i = 0; i < pad; i++) {
        padded[inlen + i] = '=';
    }
    padded[inlen + pad] = '\0';

    return cfx_base64_decode(out, outlen, padded, inlen + pad);
}

int cfx_argon2_verify(const char *encoded, const uint8_t *pwd, size_t pwdlen) {
    if (encoded == NULL || encoded[0] != '$') return CFX_ARGON2_ERR_PARAM;

    int type;
    if (strncmp(encoded + 1, "argon2id$", 9) == 0) {
        type = CFX_ARGON2ID;
        encoded += 10;
    } else if (strncmp(encoded + 1, "argon2i$", 8) == 0) {
        type = CFX_ARGON2I;
        encoded += 9;
    } else if (strncmp(encoded + 1, "argon2d$", 8) == 0) {
        type = CFX_ARGON2D;
        encoded += 9;
    } else {
        return CFX_ARGON2_ERR_PARAM;
    }

    int version;
    if (sscanf(encoded, "v=%d$", &version) != 1) {
        return CFX_ARGON2_ERR_PARAM;
    }
    encoded = strchr(encoded, '$');
    if (!encoded) return CFX_ARGON2_ERR_PARAM;
    encoded++;

    uint32_t m_cost, t_cost, p;
    if (sscanf(encoded, "m=%u,t=%u,p=%u$", &m_cost, &t_cost, &p) != 3) {
        return CFX_ARGON2_ERR_PARAM;
    }
    encoded = strchr(encoded, '$');
    if (!encoded) return CFX_ARGON2_ERR_PARAM;
    encoded++;

    const char *salt_end = strchr(encoded, '$');
    if (!salt_end) return CFX_ARGON2_ERR_PARAM;
    size_t salt_b64_len = salt_end - encoded;

    uint8_t salt[256];
    size_t saltlen = sizeof(salt);
    if (b64_decode_unpadded(salt, &saltlen, encoded, salt_b64_len) != 0) {
        return CFX_ARGON2_ERR_SALT;
    }

    encoded = salt_end + 1;
    size_t hash_b64_len = strlen(encoded);

    uint8_t stored_hash[256];
    size_t hashlen = sizeof(stored_hash);
    if (b64_decode_unpadded(stored_hash, &hashlen, encoded, hash_b64_len) != 0) {
        return CFX_ARGON2_ERR_PARAM;
    }

    uint8_t computed_hash[256];
    int rc = cfx_argon2(computed_hash, hashlen, pwd, pwdlen, salt, saltlen,
        m_cost, t_cost, p, type);
    if (rc != CFX_ARGON2_OK) {
        return rc;
    }

    /* constant-time comparison */
    uint8_t diff = 0;
    for (size_t i = 0; i < hashlen; i++) {
        diff |= stored_hash[i] ^ computed_hash[i];
    }

    cfx_memzero_s(computed_hash, sizeof(computed_hash));

    return (diff == 0) ? 0 : 1;
}

const char * cfx_argon2_strerror(int err) {
    switch (err) {
    case CFX_ARGON2_OK:         return "OK";
    case CFX_ARGON2_ERR_MEMORY: return "memory allocation failed";
    case CFX_ARGON2_ERR_PARAM:  return "invalid parameter";
    case CFX_ARGON2_ERR_SALT:   return "salt too short";
    case CFX_ARGON2_ERR_PWD:    return "password error";
    case CFX_ARGON2_ERR_OUTPUT: return "invalid output length";
    default:                    return "unknown error";
    }
}
