#ifndef CFX_MISC_H
#define CFX_MISC_H

#include <stdint.h>
#include <stddef.h>
#include <string.h>

enum cfx_str_format {
    CFX_STR_FMT_AUTO,
    CFX_STR_FMT_HEX,
    CFX_STR_FMT_ASCII,
    CFX_STR_FMT_BASE64,
    CFX_STR_FMT_BINARY
};

int hexval(int c);

/* parse hex string into exactly outlen bytes. returns 0 on success, -1 on error */
int cfx_parse_hex(const char* s, uint8_t* out, size_t outlen);
void cfx_print_hex(const uint8_t *buf, size_t len);

/* read a line from stdin, trim whitespace. caller must free() result. */
char* cfx_read_line_stdin(void);

/* check if string looks like hex (all hex chars, even length, optionally with 0x prefix) */
int cfx_looks_like_hex(const char* s);

/* parse hex string (with or without 0x prefix) into out buffer. returns bytes written, or -1 on error */
int cfx_parse_hex_auto(const char* s, uint8_t* out, size_t outlen);

/* parse string into out buffer based on mode. returns bytes written, or -1 on error */
int cfx_parse_str(const char* s, uint8_t* out, size_t outlen, enum cfx_str_format mode);

void cfx_printf_output(const uint8_t* data, size_t len, enum cfx_str_format fmt);

/* read text from file (hex/base64/ascii), securely clear temp buffer.
 * returns bytes written to out, or -1 on error. max file size 4KB. */
int cfx_read_file_text(const char* path, uint8_t* out, size_t outlen, enum cfx_str_format mode);

/* read raw binary from file. expects file to be exactly outlen bytes.
 * returns 0 on success, -1 on error. */
int cfx_read_file_bin(const char* path, uint8_t* out, size_t outlen);

/* detect file format based on contents.
 * if out_len is not NULL, writes the decoded data size in bytes. */
enum cfx_str_format cfx_detect_file_format(const char* path, size_t* out_len);

#endif  /* CFX_MISC_H */
