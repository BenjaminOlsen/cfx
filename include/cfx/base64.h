#ifndef CFX_BASE64_H
#define CFX_BASE64_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

size_t cfx_base64_enc_len(size_t bin_len);
size_t cfx_base64_dec_max_len(size_t b64_len); /* upper bound, ignoring whitespace */

/*
 * Encode canonical Base64 (RFC 4648).
 * out may be NULL to query required length via *out_len.
 * If out is non-NULL, *out_len is in/out: capacity on input, bytes written on success.
 * Returns 0 on success, -1 on insufficient output buffer.
 */
int cfx_base64_encode(char* out, size_t* out_len, const uint8_t* in, size_t in_len);

/*
 * Decode canonical Base64, ignoring ASCII whitespace.
 * out may be NULL to query exact decoded length via *out_len.
 * If out is non-NULL, *out_len is in/out: capacity on input, bytes written on success.
 *
 * Returns:
 *  0  success
 * -1  insufficient output buffer
 * -2  invalid input (bad char, bad padding, bad trailing chars, non-canonical tail bits)
 */
int cfx_base64_decode(uint8_t* out, size_t* out_len, const char* in, size_t in_len);

#ifdef __cplusplus
}
#endif

#endif
