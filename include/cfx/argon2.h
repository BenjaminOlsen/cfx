/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

/*
 * argon2.h - Argon2 password hashing (RFC 9106)
 *
 * Argon2d:  data-dependent (fastest, side-channel vulnerable)
 * Argon2i:  data-independent (side-channel resistant)
 * Argon2id: hybrid (recommended)
 */

#ifndef CFX_ARGON2_H
#define CFX_ARGON2_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define CFX_ARGON2D  0
#define CFX_ARGON2I  1
#define CFX_ARGON2ID 2

#define CFX_ARGON2_VERSION 0x13

#define CFX_ARGON2_MIN_OUTLEN   4
#define CFX_ARGON2_MIN_SALT_LEN 8
#define CFX_ARGON2_MIN_TIME     1
#define CFX_ARGON2_MIN_MEMORY   8
#define CFX_ARGON2_MIN_LANES    1

#define CFX_ARGON2_OK           0
#define CFX_ARGON2_ERR_MEMORY  -1
#define CFX_ARGON2_ERR_PARAM   -2
#define CFX_ARGON2_ERR_SALT    -3
#define CFX_ARGON2_ERR_PWD     -4
#define CFX_ARGON2_ERR_OUTPUT  -5

/* recommended: m=65536 (64MB), t=3, p=4 */
int cfx_argon2id(uint8_t* out, size_t outlen,
                 const uint8_t* pwd, size_t pwdlen,
                 const uint8_t* salt, size_t saltlen,
                 uint32_t m_cost, uint32_t t_cost, uint32_t p);

int cfx_argon2d(uint8_t* out, size_t outlen,
                const uint8_t* pwd, size_t pwdlen,
                const uint8_t* salt, size_t saltlen,
                uint32_t m_cost, uint32_t t_cost, uint32_t p);

int cfx_argon2i(uint8_t* out, size_t outlen,
                const uint8_t* pwd, size_t pwdlen,
                const uint8_t* salt, size_t saltlen,
                uint32_t m_cost, uint32_t t_cost, uint32_t p);

int cfx_argon2(uint8_t* out, size_t outlen,
               const uint8_t* pwd, size_t pwdlen,
               const uint8_t* salt, size_t saltlen,
               uint32_t m_cost, uint32_t t_cost, uint32_t p,
               int type);

/* PHC string format: $argon2id$v=19$m=65536,t=3,p=4$<salt>$<hash> */
int cfx_argon2_encode(char* out, size_t outlen,
                      const uint8_t* hash, size_t hashlen,
                      const uint8_t* salt, size_t saltlen,
                      uint32_t m_cost, uint32_t t_cost, uint32_t p,
                      int type);

/* returns 0 if password matches, 1 if not, <0 on error */
int cfx_argon2_verify(const char* encoded, const uint8_t* pwd, size_t pwdlen);

const char* cfx_argon2_strerror(int err);

#ifdef __cplusplus
}
#endif

#endif /* CFX_ARGON2_H */
