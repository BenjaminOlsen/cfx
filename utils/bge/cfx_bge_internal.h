#ifndef CFX_BGE_INTERNAL_H
#define CFX_BGE_INTERNAL_H

#include "cfx/argon2.h"
#include "cfx/aead_chacha20_poly1305.h"
#include "cfx/base64.h"
#include "cfx/macros.h"
#include "cfx/memory.h"
#include "cfx/rand.h"
#include "cfx_utils_common.h"

#include <errno.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BGE_MAGIC          "BGE"
#define BGE_VERSION_STR    "3.1.1"
#define BGE_TAG_LEN        16
#define BGE_DEFAULT_M      0x10000
#define BGE_DEFAULT_T      3
#define BGE_DEFAULT_P      4
#define BGE_MIN_M          64
#define BGE_MAX_M          4194304
#define BGE_MAX_T          100
#define BGE_MAX_P          16
#define BGE_ARMOR_HEADER   "-----BEGIN BGE MESSAGE-----"
#define BGE_ARMOR_FOOTER   "-----END BGE MESSAGE-----"
#define BGE_FILE_VERSION   5
#define BGE_STREAM_VERSION 3
#define BGE_HEADER_LEN     56
#define BGE_VERIFIER_LEN   16
#define BGE_AAD_LEN        (BGE_HEADER_LEN + BGE_VERIFIER_LEN)
#define BGE_MIN_FILE       (BGE_AAD_LEN + BGE_TAG_LEN)

typedef struct {
    uint8_t  magic[3];
    uint8_t  version;
    uint32_t m_cost;
    uint32_t t_cost;
    uint32_t p_cost;
    uint8_t  salt[16];
    uint8_t  nonce[24];
} bge_header;
CFX_STATIC_ASSERT(sizeof(bge_header) == BGE_HEADER_LEN, bge_header_packing);

int bge_is_armored(const uint8_t *buf, size_t len);
int bge_armor_encode(const uint8_t *bin, size_t bin_len,
                     uint8_t **out, size_t *out_len);
int bge_armor_decode(const uint8_t *text, size_t text_len,
                     uint8_t **out, size_t *out_len);
int bge_encrypt_file(int argc, char **argv);
int bge_decrypt_file(int argc, char **argv);

#endif
