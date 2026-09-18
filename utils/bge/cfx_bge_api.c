#include "bge.h"
#include "cfx_bge_internal.h"

#include <limits.h>

int bge_is_armored(const uint8_t *buf, size_t len) {
    return len >= 27 && memcmp(buf, BGE_ARMOR_HEADER, 27) == 0;
}

int bge_armor_decode(const uint8_t *text, size_t text_len,
                     uint8_t **out, size_t *out_len) {
    const char *start = (const char *)text;
    const char *body = memchr(start, '\n', text_len);
    if (!body) return -1;
    body++;

    size_t body_len = text_len - (size_t)(body - start);
    size_t footer_len = strlen(BGE_ARMOR_FOOTER);
    const char *footer = NULL;
    for (size_t i = 0; i + footer_len <= body_len; i++) {
        if (memcmp(body + i, BGE_ARMOR_FOOTER, footer_len) == 0) {
            footer = body + i;
            break;
        }
    }
    if (!footer) return -1;

    size_t encoded_len = (size_t)(footer - body);
    size_t decoded_len = cfx_base64_dec_max_len(encoded_len);
    uint8_t *decoded = malloc(decoded_len ? decoded_len : 1);
    if (!decoded) return -1;
    if (cfx_base64_decode(decoded, &decoded_len, body, encoded_len) != 0) {
        free(decoded);
        return -1;
    }

    *out = decoded;
    *out_len = decoded_len;
    return 0;
}

static int bge_derive_key(const bge_header *header,
                          const uint8_t *passphrase, size_t passphrase_len,
                          uint8_t key[48]) {
    
    uint32_t m = cfx_load32_le(&header->m_cost);
    uint32_t t = cfx_load32_le(&header->t_cost);
    uint32_t p = cfx_load32_le(&header->p_cost);

    if (m < BGE_MIN_M || m > BGE_MAX_M ||
        t < 1 || t > BGE_MAX_T || p < 1 || p > BGE_MAX_P) {
        return -2;
    }

    
    int rc = cfx_argon2id(key, 48, passphrase, passphrase_len, header->salt, sizeof(header->salt), m, t, p) == 0;
    return rc ? 0 : -1;
}

void cfx_bge_free(void *buffer, size_t buffer_len) {
    if (!buffer) return;
    cfx_memzero_s(buffer, buffer_len);
    free(buffer);
}


#define MARK_FAIL_AND_CLEANUP(x) ret = x; goto cleanup;

int cfx_bge_decrypt_stream(FILE* input, FILE* output, const uint8_t *passphrase, size_t passphrase_len) {

    if (!input || !output || !passphrase || passphrase_len == 0) {
        return -1;
    }

    int ret = 0;
    bge_header header;
    uint8_t key[48] = {0};
    uint8_t *plaintext  = (uint8_t*)malloc(CFX_STREAM_CHUNK_SIZE);
    uint8_t *current    = (uint8_t*)malloc(CFX_STREAM_CHUNK_SIZE + CFX_STREAM_TAG_SIZE);
    uint8_t *next       = (uint8_t*)malloc(CFX_STREAM_CHUNK_SIZE + CFX_STREAM_TAG_SIZE);

    if (!plaintext || !current || !next) {
        MARK_FAIL_AND_CLEANUP(-1);
    }

    if (fread(&header, 1, sizeof(header), input) != sizeof(header)) {
        MARK_FAIL_AND_CLEANUP(-1);
    }
    if (header.version != BGE_STREAM_VERSION) {
        MARK_FAIL_AND_CLEANUP(-1);
    }
    if (memcmp(header.magic, BGE_MAGIC, 3) != 0) {
        MARK_FAIL_AND_CLEANUP(-1);
    }
    
    if (bge_derive_key(&header, passphrase, passphrase_len, key) != 0) {
        MARK_FAIL_AND_CLEANUP(-1);
    }


    uint8_t verifier_in[BGE_VERIFIER_LEN] = {0};
    if (fread(verifier_in, 1, sizeof(verifier_in), input) != sizeof(verifier_in)) {
        ret = -1;
        goto cleanup;
    }
    unsigned int diff = 0;
    for (size_t i = 0; i < BGE_VERIFIER_LEN; ++i) {
        diff |= key[32 + i] ^ verifier_in[i];
    }
    if (diff != 0) {
        ret = -1;
        goto cleanup;
    }

    size_t current_len = fread(current, 1, CFX_STREAM_CHUNK_SIZE + CFX_STREAM_TAG_SIZE, input);
    if (ferror(input)) {
        ret = -1;
        goto cleanup;
    }
    uint64_t chunk_cnt = 0;
    while (1) {
        uint8_t *src = current;
        uint8_t *dst = plaintext;
        
        size_t next_len = fread(next, 1, CFX_STREAM_CHUNK_SIZE + CFX_STREAM_TAG_SIZE, input);
        if (ferror(input)) {
            ret = -1;
            goto cleanup;
        }
        int final = (next_len == 0 && feof(input));
        if (current_len < CFX_STREAM_TAG_SIZE) {
            ret = -1;
            goto cleanup;
        }
        size_t len = current_len - CFX_STREAM_TAG_SIZE;
        int rc = cfx_stream_xchacha20_poly1305_decrypt_chunk(
            dst, src, len, src + len,
            chunk_cnt, final, key, header.nonce);
        if (rc != 0) {
            ret = -1;
            goto cleanup;
        }
        if (fwrite(plaintext, 1, len, output) != len) {
            ret = -1;
            goto cleanup;
        }
        uint8_t *tmp = current;
        current = next;
        next = tmp;
        current_len = next_len;
        if (final) break;
        ++chunk_cnt;
    }

    if (fflush(output) == EOF) {
        ret = -1;
        goto cleanup;
    }



cleanup:
    cfx_memzero_s(&header, sizeof(header));
    cfx_memzero_s(key, sizeof(key));
    cfx_bge_free(current, CFX_STREAM_CHUNK_SIZE + CFX_STREAM_TAG_SIZE);
    cfx_bge_free(next, CFX_STREAM_CHUNK_SIZE + CFX_STREAM_TAG_SIZE);
    cfx_bge_free(plaintext, CFX_STREAM_CHUNK_SIZE);
    return ret;
}


        


int cfx_bge_encrypt_stream(FILE *input, FILE *output, const uint8_t *passphrase, size_t passphrase_len) {

    if (!input || !output || !passphrase || passphrase_len == 0) {
        return -1;
    }

    int ret = 0;
    bge_header header;
    uint8_t key[48] = {0};
    uint8_t *current    = (uint8_t*)malloc(CFX_STREAM_CHUNK_SIZE);
    uint8_t *next       = (uint8_t*)malloc(CFX_STREAM_CHUNK_SIZE);
    uint8_t *cipher     = (uint8_t*)malloc(CFX_STREAM_CHUNK_SIZE + CFX_STREAM_TAG_SIZE);
    if (!current || !next || !cipher) {
        ret = -1;
        goto cleanup;
    }


    memcpy(header.magic, BGE_MAGIC, 3);
    header.version = BGE_STREAM_VERSION;
    cfx_store32_le(&header.m_cost, BGE_DEFAULT_M);
    cfx_store32_le(&header.t_cost, BGE_DEFAULT_T);
    cfx_store32_le(&header.p_cost, BGE_DEFAULT_P);
    cfx_rand_bytes_os(header.salt, sizeof(header.salt));
    cfx_rand_bytes_os(header.nonce, sizeof(header.nonce));

    int rc = bge_derive_key(&header, passphrase, passphrase_len, key);
    if (rc != 0) {
        ret = rc;
        goto cleanup;
    }

    if (fwrite(&header, 1, sizeof(header), output) != sizeof(header)) {
        ret = -1;
        goto cleanup;
    }
    if (fwrite(key + 32, 1, BGE_VERIFIER_LEN, output) != BGE_VERIFIER_LEN) {
        ret = -1;
        goto cleanup;
    }


    uint64_t chunk_cnt = 0;
    size_t current_len = fread(current, 1, CFX_STREAM_CHUNK_SIZE, input);
    if (ferror(input)) {
        ret = -1;
        goto cleanup;
    }

    while (1) {
        uint8_t *src = current;
        uint8_t *dst = cipher;
        uint8_t *tag = dst + current_len;
        size_t next_len = fread(next, 1, CFX_STREAM_CHUNK_SIZE, input);
        if (ferror(input)) {
            ret = -1;
            goto cleanup;
        }
        int final = (next_len == 0 && feof(input));
        if (chunk_cnt >= (UINT64_C(1) << 31)) {
            ret = -1;
            goto cleanup;
        }
        rc = cfx_stream_xchacha20_poly1305_encrypt_chunk(
                dst, tag, src, current_len, chunk_cnt, final, key, header.nonce);
        if (rc != 0) {
            ret = -1;
            goto cleanup;
        }

        size_t output_len = current_len + CFX_STREAM_TAG_SIZE;
        if (fwrite(cipher, 1, output_len, output) != output_len) {
            ret = -1;
            goto cleanup;
        }
        uint8_t *tmp = current;
        current = next;
        next = tmp;
        current_len = next_len;
        if (final) break;
        ++chunk_cnt;
    }
    if (fflush(output) == EOF) {
        ret = -1;
        goto cleanup;
    }

cleanup:
    cfx_memzero_s(&header, sizeof(header));
    cfx_memzero_s(key, sizeof(key));
    cfx_bge_free(current, CFX_STREAM_CHUNK_SIZE);
    cfx_bge_free(next, CFX_STREAM_CHUNK_SIZE);
    cfx_bge_free(cipher, CFX_STREAM_CHUNK_SIZE + CFX_STREAM_TAG_SIZE);
    return ret;
}

int cfx_bge_encrypt(const uint8_t *plaintext, size_t plaintext_len,
                    const uint8_t *passphrase, size_t passphrase_len,
                    uint8_t **output, size_t *output_len) {
    if (!output || !output_len ||
        (!plaintext && plaintext_len != 0) ||
        (!passphrase && passphrase_len != 0)) {
        return -1;
    }
    *output = NULL;
    *output_len = 0;

    size_t chunks = plaintext_len / CFX_STREAM_CHUNK_SIZE + 1;
    if (chunks > (SIZE_MAX - BGE_AAD_LEN - plaintext_len) / CFX_STREAM_TAG_SIZE) {
        return -1;
    }

    size_t total = BGE_AAD_LEN + plaintext_len + (chunks * CFX_STREAM_TAG_SIZE);

    uint8_t *result = malloc(total);
    if (!result) return -1;

    bge_header header;
    memcpy(header.magic, BGE_MAGIC, 3);
    header.version = BGE_STREAM_VERSION;
    cfx_store32_le(&header.m_cost, BGE_DEFAULT_M);
    cfx_store32_le(&header.t_cost, BGE_DEFAULT_T);
    cfx_store32_le(&header.p_cost, BGE_DEFAULT_P);
    cfx_rand_bytes_os(header.salt, sizeof(header.salt));
    cfx_rand_bytes_os(header.nonce, sizeof(header.nonce));

    uint8_t key[48] = {0};
    int rc = bge_derive_key(&header, passphrase, passphrase_len, key);
    if (rc != 0) {
        cfx_bge_free(result, total);
        return rc;
    }

    memcpy(result, &header, sizeof(header));
    memcpy(result + BGE_HEADER_LEN, key + 32, BGE_VERIFIER_LEN);

    const uint8_t *src = plaintext;
    uint8_t *dst = result + BGE_AAD_LEN;
    size_t remaining = plaintext_len;
    uint64_t counter = 0;
    for (;;) {
        size_t len = remaining > CFX_STREAM_CHUNK_SIZE ? CFX_STREAM_CHUNK_SIZE : remaining;
        int final = remaining <= CFX_STREAM_CHUNK_SIZE;
        uint8_t *tag = dst + len;

        rc = cfx_stream_xchacha20_poly1305_encrypt_chunk(
            dst, tag, src, len, counter, final, key, header.nonce);
        if (rc != 0) {
            cfx_memzero_s(key, sizeof(key));
            cfx_bge_free(result, total);
            return -1;
        }

        dst += len + CFX_STREAM_TAG_SIZE;
        if (src) {
            src += len;
        }
        remaining -= len;
        counter++;
        if (final) break;
    }

    cfx_memzero_s(key, sizeof(key));
    *output = result;
    *output_len = total;
    return 0;
}

static int bge_decrypt_binary(const uint8_t *input, size_t input_len,
                              const uint8_t *passphrase, size_t passphrase_len,
                              uint8_t **plaintext, size_t *plaintext_len) {
    if (input_len < BGE_MIN_FILE || memcmp(input, BGE_MAGIC, 3) != 0) {
        return -2;
    }

    bge_header header;
    memcpy(&header, input, sizeof(header));
    if (header.version != BGE_FILE_VERSION &&
        header.version != BGE_STREAM_VERSION)
        return -2;

    uint8_t key[48] = {0};
    int rc = bge_derive_key(&header, passphrase, passphrase_len, key);
    if (rc != 0) return rc;

    unsigned diff = 0;
    for (size_t i = 0; i < BGE_VERIFIER_LEN; i++)
        diff |= input[BGE_HEADER_LEN + i] ^ key[32 + i];
    if (diff != 0) {
        cfx_memzero_s(key, sizeof(key));
        return -3;
    }

    const uint8_t *payload = input + BGE_AAD_LEN;
    size_t payload_len = input_len - BGE_AAD_LEN;

    if (header.version == BGE_FILE_VERSION) {
        if (payload_len < BGE_TAG_LEN) {
            cfx_memzero_s(key, sizeof(key));
            return -2;
        }
        size_t len = payload_len - BGE_TAG_LEN;
        uint8_t *out = malloc(len ? len : 1);
        if (!out) {
            cfx_memzero_s(key, sizeof(key));
            return -1;
        }
        rc = cfx_xchacha20_poly1305_decrypt(
            out, payload, len, input, BGE_AAD_LEN, key, header.nonce,
            payload + len);
        cfx_memzero_s(key, sizeof(key));
        if (rc != 0) {
            cfx_bge_free(out, len);
            return -3;
        }
        *plaintext = out;
        *plaintext_len = len;
        return 0;
    }

    uint8_t *out = malloc(payload_len ? payload_len : 1);
    if (!out) {
        cfx_memzero_s(key, sizeof(key));
        return -1;
    }

    const uint8_t *p = payload;
    const uint8_t *end = input + input_len;
    size_t written = 0;
    uint64_t counter = 0;
    while (p < end) {
        size_t remaining = (size_t)(end - p);
        if (remaining < CFX_STREAM_TAG_SIZE) {
            cfx_memzero_s(key, sizeof(key));
            cfx_bge_free(out, payload_len);
            return -2;
        }
        int final = remaining <= CFX_STREAM_CHUNK_SIZE + CFX_STREAM_TAG_SIZE;
        size_t len = final ? remaining - CFX_STREAM_TAG_SIZE
                           : CFX_STREAM_CHUNK_SIZE;
        rc = cfx_stream_xchacha20_poly1305_decrypt_chunk(
            out + written, p, len, p + len, counter, final,
            key, header.nonce);
        if (rc != 0) {
            cfx_memzero_s(key, sizeof(key));
            cfx_bge_free(out, payload_len);
            return -3;
        }
        written += len;
        p += len + CFX_STREAM_TAG_SIZE;
        counter++;
    }

    cfx_memzero_s(key, sizeof(key));
    *plaintext = out;
    *plaintext_len = written;
    return 0;
}

int cfx_bge_decrypt(const uint8_t *input, size_t input_len,
                    const uint8_t *passphrase, size_t passphrase_len,
                    uint8_t **plaintext, size_t *plaintext_len) {
    if (!plaintext || !plaintext_len || !input ||
        (!passphrase && passphrase_len != 0))
        return -1;
    *plaintext = NULL;
    *plaintext_len = 0;

    uint8_t *decoded = NULL;
    size_t decoded_len = 0;
    if (bge_is_armored(input, input_len)) {
        if (bge_armor_decode(input, input_len, &decoded, &decoded_len) != 0)
            return -2;
        input = decoded;
        input_len = decoded_len;
    }

    int rc = bge_decrypt_binary(input, input_len, passphrase, passphrase_len,
                                plaintext, plaintext_len);
    if (decoded) cfx_bge_free(decoded, decoded_len);
    return rc;
}

