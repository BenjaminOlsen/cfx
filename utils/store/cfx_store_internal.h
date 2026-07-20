/* cfx_bge_internal.h -- shared types, constants, and declarations for BGE */

#ifndef CFX_STORE_INTERNAL_H
#define CFX_STORE_INTERNAL_H

#include "cfx/argon2.h"
#include "cfx/aead_chacha20_poly1305.h"
#include "cfx/sha256.h"
#include "cfx/rand.h"
#include "cfx/memory.h"
#include "cfx/base64.h"
#include "cfx/macros.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <errno.h>
#include <time.h>
#include <sys/stat.h>

#ifdef _WIN32
#include <io.h>
#include <fcntl.h>
#include <direct.h>
#include <shlobj.h>
#define mkdir_p(path) _mkdir(path)
#ifndef S_ISDIR
#define S_ISDIR(m) (((m) & _S_IFMT) == _S_IFDIR)
#endif
#else
#include <unistd.h>
#include <fcntl.h>
#include <sys/file.h>
#include <pwd.h>
#include <termios.h>
#define mkdir_p(path) mkdir(path, 0700)
#endif

#include "cfx_cmd.h"
#include "cfx_keyfile.h"
#include "cfx_utils_common.h"

/*  common constants  */

#define BGE_MAGIC         "BGE"
#define BGE_VERSION_STR   "3.1.0"
#define BGE_TAG_LEN       16

#define BGE_DEFAULT_M     0x10000  /* 64 MB (in KB) */
#define BGE_DEFAULT_T     3
#define BGE_DEFAULT_P     4

/* m=8 is RFC-minimum but produces collisions; 64 KB is the practical floor */
#define BGE_MIN_M         64
#define BGE_MAX_M         4194304
#define BGE_MAX_T         100
#define BGE_MAX_P         16

#define BGE_ARMOR_HEADER  "-----BEGIN BGE MESSAGE-----"
#define BGE_ARMOR_FOOTER  "-----END BGE MESSAGE-----"

/*  file encryption constants (used by cfx_bge_file.c)  */

#define BGE_FILE_VERSION  5       /* single-shot file encryption format */
#define BGE_STREAM_VERSION 3      /* streaming file encryption */
#define BGE_HEADER_LEN    56
#define BGE_VERIFIER_LEN  16
#define BGE_AAD_LEN       (BGE_HEADER_LEN + BGE_VERIFIER_LEN)  /* 72 */
#define BGE_MIN_FILE      (BGE_AAD_LEN + BGE_TAG_LEN)          /* 88 */

typedef struct {
    uint8_t  magic[3];
    uint8_t  version;
    uint32_t m_cost;    /* little-endian on disk */
    uint32_t t_cost;    /* little-endian on disk */
    uint32_t p_cost;    /* little-endian on disk */
    uint8_t  salt[16];
    uint8_t  nonce[24];
} bge_header;
CFX_STATIC_ASSERT(sizeof(bge_header) == BGE_HEADER_LEN, bge_header_packing);

/*  store constants  */

#define BGE_STORE_VERSION     4
#define BGE_STORE_HDR_LEN    8
#define BGE_STORE_SLOT_LEN   116
#define BGE_STORE_MAX_SLOTS  8
#define BGE_STORE_DEK_LEN    32
#define BGE_STORE_MIN_FILE   (BGE_STORE_HDR_LEN + BGE_STORE_SLOT_LEN + 24 + BGE_TAG_LEN)

/*  store structs  */

typedef struct {
    uint8_t magic[3];
    uint8_t version;
    uint8_t slot_count;
    uint8_t reserved[3];
} bge_v4_fixed_header;
CFX_STATIC_ASSERT(sizeof(bge_v4_fixed_header) == BGE_STORE_HDR_LEN, bge_v4_hdr_packing);

typedef struct {
    uint32_t m_cost;
    uint32_t t_cost;
    uint32_t p_cost;
    uint8_t  salt[16];
    uint8_t  slot_nonce[24];
    uint8_t  verifier[16];
    uint8_t  wrapped_dek[32];
    uint8_t  wrap_tag[16];
} bge_v4_slot;
CFX_STATIC_ASSERT(sizeof(bge_v4_slot) == BGE_STORE_SLOT_LEN, bge_v4_slot_packing);

typedef struct {
    bge_v4_fixed_header hdr;
    bge_v4_slot         slots[BGE_STORE_MAX_SLOTS];
    uint8_t             data_nonce[24];
    uint8_t             dek[BGE_STORE_DEK_LEN];
    int                 matched_slot;
    uint8_t            *file_buf;
    size_t              file_len;
} bge_v4_store;

/*  unified store — all stores are now v4 multi-slot  */

typedef struct {
    bge_v4_store v4;
} bge_ustore;

/*  global  */

extern const char *g_passphrase_arg;

/*  cfx_bge_crypto.c  */

void bge_v4_store_wipe(bge_v4_store *s);
void bge_ustore_wipe(bge_ustore *us);

/* slot serialization */
void bge_v4_slot_from_buf(bge_v4_slot *slot, const uint8_t buf[BGE_STORE_SLOT_LEN]);
void bge_v4_slot_to_buf(uint8_t buf[BGE_STORE_SLOT_LEN], const bge_v4_slot *slot);

/* v4 crypto */
int bge_v4_wrap_dek(const char *pwd, size_t pwd_len,
                    uint32_t m, uint32_t t, uint32_t p,
                    const uint8_t dek[BGE_STORE_DEK_LEN], bge_v4_slot *slot);
int bge_v4_safe_write(const char *path,
                      const bge_v4_fixed_header *hdr,
                      const bge_v4_slot *slots, int slot_count,
                      const uint8_t data_nonce[24],
                      const uint8_t *ct, size_t ct_len,
                      const uint8_t tag[BGE_TAG_LEN]);
int bge_v4_encrypt_write(const char *path,
                         const uint8_t *pt, size_t pt_len,
                         const bge_v4_fixed_header *hdr,
                         const bge_v4_slot *slots, int slot_count,
                         const uint8_t dek[BGE_STORE_DEK_LEN]);
int bge_v4_init_write(const char *path, const uint8_t *pt, size_t pt_len,
                      const char *pwd, size_t pwd_len,
                      uint32_t m, uint32_t t, uint32_t p);
int bge_v4_init_to_buf(uint8_t **out, size_t *out_len,
                       const uint8_t *pt, size_t pt_len,
                       const char *pwd, size_t pwd_len,
                       uint32_t m, uint32_t t, uint32_t p);
int bge_v4_authenticate_buf(const uint8_t *file_buf, size_t file_len,
                            const char *pwd, size_t pwd_len,
                            bge_v4_store *store);
int bge_v4_decrypt_store(const bge_v4_store *store,
                         uint8_t **pt_out, size_t *pt_len);

/* unified dispatch */
int bge_uauthenticate(const char *path, const char *pwd, size_t pwd_len,
                      bge_ustore *us);
int bge_udecrypt(const bge_ustore *us, uint8_t **pt_out, size_t *pt_len);
int bge_uwrite(const char *path, const uint8_t *pt, size_t pt_len,
               const bge_ustore *us);

/*  cfx_bge_store.c  */

const uint8_t *store_name_by_index(const uint8_t *pt, size_t pt_len,
                                   unsigned idx, uint16_t *name_len);
const uint8_t *store_get(const uint8_t *pt, size_t pt_len,
                         const char *name, size_t *vlen);
uint8_t *store_set(const uint8_t *pt, size_t pt_len,
                   const char *name, const uint8_t *val, size_t val_len,
                   size_t *new_len);
uint8_t *store_rm(const uint8_t *pt, size_t pt_len,
                  const char *name, size_t *new_len);
uint8_t *store_swap(const uint8_t *pt, size_t pt_len,
                    unsigned idx_a, unsigned idx_b, size_t *new_len);
uint8_t *store_sort(const uint8_t *pt, size_t pt_len, size_t *new_len);
uint8_t *store_to_text(const uint8_t *pt, size_t pt_len, size_t *text_len);

/*  cfx_bge_file.c  */

int bge_is_armored(const uint8_t *buf, size_t len);
int bge_armor_encode(const uint8_t *bin, size_t bin_len,
                     uint8_t **out, size_t *out_len);
int bge_armor_decode(const uint8_t *text, size_t text_len,
                     uint8_t **out, size_t *out_len);
int bge_encrypt_file(int argc, char **argv);
int bge_decrypt_file(int argc, char **argv);

/*  cfx_bge_grace.c  */

int  grace_check(const char *store_path, char *pwd_out, size_t pwd_bufsz);
void grace_stamp(const char *store_path, const char *pwd, size_t pwd_len);
void grace_delete(void);

/*  cfx_bge_common.c (I/O, path, and utility helpers)  */

int bge_read_all(FILE *f, uint8_t **out, size_t *out_len);
int bge_read_secret(const char *prompt, char *buf, size_t bufsz);
int bge_read_passphrase(const char *prompt, char *buf, size_t bufsz);
int bge_read_visible(const char *prompt, char *buf, size_t bufsz);
int prompt_passphrase(char *pwd, size_t pwdsz);

int get_cfx_dir(char *buf, size_t bufsz);
int bge_default_path(char *buf, size_t bufsz);
int ensure_cfx_dir(void);

int ct_pwd_match(const char *pw1, int pw1_len, size_t pw1_bufsz,
                 const char *pw2, int pw2_len, size_t pw2_bufsz);
int is_numeric(const char *s);
void store_print_names(const uint8_t *pt, size_t pt_len);
int clip_copy(const uint8_t *data, size_t len);

#endif /* CFX_STORE_INTERNAL_H */
