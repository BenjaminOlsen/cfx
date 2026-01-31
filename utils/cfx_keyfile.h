/* cfx_keyfile.h - encrypted key file format (CFX\x01) shared between utilities */

#ifndef CFX_KEYFILE_H
#define CFX_KEYFILE_H

#include <stdint.h>
#include <stddef.h>

/*
 * Encrypted key file layout (104 bytes):
 *
 *   Offset  Size  Field
 *   0       4     magic: "CFX\x01"
 *   4       4     argon2 m_cost  (LE32)
 *   8       4     argon2 t_cost  (LE32)
 *   12      4     argon2 parallelism (LE32)
 *   16      16    salt
 *   32      24    nonce
 *   56      32    ciphertext (encrypted seed)
 *   88      16    auth tag (Poly1305)
 *
 * The header (bytes 0-55) is used as AAD for the AEAD.
 */
#define CFX_KEY_MAGIC      "CFX\x01"
#define CFX_KEY_HEADER_LEN 56
#define CFX_KEY_FILE_LEN   104

/* default argon2id parameters */
#define CFX_KEY_DEFAULT_M  0x10000  /* 64 MB (in KB) */
#define CFX_KEY_DEFAULT_T  3       /* 3 iterations */
#define CFX_KEY_DEFAULT_P  4       /* 4 lanes */

/* read a passphrase from the terminal with echo disabled.
 * returns length of password, or 0 if empty. */
int cfx_key_read_secret_console(const char *prompt, char *buf, size_t bufsz);

/* decrypt a CFX\x01 encrypted key file buffer (104 bytes).
 * prompts for passphrase on the terminal.
 * returns 0 on success, -1 on error. */
int cfx_key_decrypt(const uint8_t *file_buf, uint8_t *seed_out);

/* encrypt seed and write to path as CFX\x01 file.
 * returns 0 on success, -1 on error. */
int cfx_key_write_encrypted(const char *path, const uint8_t *seed, size_t seed_len, const char *pwd, size_t pwd_len);

#endif /* CFX_KEYFILE_H */
