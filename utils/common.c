#include "common.h"
#include "cfx/base64.h"
#include "cfx/memory.h"

#include <stdlib.h>
#include <ctype.h>

int hexval(int c) {
    if ('0' <= c && c <= '9') return c - '0';
    if ('a' <= c && c <= 'f') return 10 + c - 'a';
    if ('A' <= c && c <= 'F') return 10 + c - 'A';
    return -1;
}

/* parse hex string into exactly outlen bytes. returns 0 on success, -1 on error */
int cfx_parse_hex(const char* s, uint8_t* out, size_t outlen) {
    if (s[0]=='0' && (s[1]=='x' || s[1]=='X')) s += 2;

    size_t n = strlen(s);
    if (n != outlen*2) return -1;

    for (size_t i = 0; i < outlen; ++i) {
        int hi = hexval(s[2*i]);
        int lo = hexval(s[2*i+1]);
        if (hi < 0 || lo < 0) return -1;
        out[i] = (uint8_t)((hi << 4) | lo);
    }
    return 0;
}

int cfx_parse_str_exact(const char* s, uint8_t* out, size_t outlen) {
    if (!s) return -1;
    if (strncmp(s, "b64:", 4) == 0) {
        size_t dec_len = outlen;
        if (cfx_base64_decode(out, &dec_len, s + 4, strlen(s + 4)) != 0 ||
            dec_len != outlen)
            return -1;
        return 0;
    }
    return cfx_parse_hex(s, out, outlen);
}

void cfx_print_hex(const uint8_t *buf, size_t len) {
    for (size_t i = 0; i < len; ++i) {
        printf("%02x", (unsigned)buf[i]);
    }
}

char* cfx_read_line_stdin(void) {
    size_t cap = 1024;
    size_t len = 0;
    char* buf = (char*)malloc(cap);
    if (!buf) return NULL;

    int c;
    while ((c = getchar()) != EOF && c != '\n' && c != '\r') {
        if (len + 1 >= cap) {
            cap *= 2;
            char* newbuf = (char*)realloc(buf, cap);
            if (!newbuf) { free(buf); return NULL; }
            buf = newbuf;
        }
        buf[len++] = (char)c;
    }
    buf[len] = '\0';

    char* start = buf;
    while (*start && isspace((unsigned char)*start)) start++;
    char* end = start + strlen(start);
    while (end > start && isspace((unsigned char)end[-1])) end--;
    *end = '\0';

    if (start != buf) {
        memmove(buf, start, strlen(start) + 1);
    }
    return buf;
}

int cfx_looks_like_hex(const char* s) {
    if (!s || !*s) return 0;
    if (s[0] == '0' && (s[1] == 'x' || s[1] == 'X')) s += 2;
    size_t len = strlen(s);
    if (len == 0 || len % 2 != 0) return 0;
    for (size_t i = 0; i < len; ++i) {
        if (hexval(s[i]) < 0) return 0;
    }
    return 1;
}

int cfx_parse_hex_auto(const char* s, uint8_t* out, size_t outlen) {
    if (!s || !*s) return -1;
    if (s[0] == '0' && (s[1] == 'x' || s[1] == 'X')) s += 2;
    size_t slen = strlen(s);
    if (slen % 2 != 0) return -1;
    size_t nbytes = slen / 2;
    if (nbytes > outlen) return -1;
    for (size_t i = 0; i < nbytes; ++i) {
        int hi = hexval(s[2*i]);
        int lo = hexval(s[2*i+1]);
        if (hi < 0 || lo < 0) return -1;
        out[i] = (uint8_t)((hi << 4) | lo);
    }
    return (int)nbytes;
}

int cfx_parse_str(const char* s, uint8_t* out, size_t outlen, enum cfx_str_format mode) {
    if (!s) return -1;
    memset(out, 0, outlen);

    if (mode == CFX_STR_FMT_BASE64) {
        size_t out_len = outlen;
        if (cfx_base64_decode(out, &out_len, s, strlen(s)) != 0) {
            return -1;
        }
        return (int)out_len;
    }

    int use_hex = (mode == CFX_STR_FMT_HEX) || (mode == CFX_STR_FMT_AUTO && cfx_looks_like_hex(s));
    if (mode == CFX_STR_FMT_ASCII || !use_hex) {
        size_t kl = strlen(s);
        if (kl > outlen) kl = outlen;
        memcpy(out, s, kl);
        return (int)kl;
    }

    return cfx_parse_hex_auto(s, out, outlen);
}


void cfx_printf_output(const uint8_t* data, size_t len, enum cfx_str_format fmt) {
    if (fmt == CFX_STR_FMT_BASE64) {
        size_t b64_len = cfx_base64_enc_len(len);
        char* b64 = (char*)malloc(b64_len + 1);
        if (b64) {
            cfx_base64_encode(b64, &b64_len, data, len);
            b64[b64_len] = '\0';
            printf("%s", b64);
            free(b64);
        }
    } else {
        cfx_print_hex(data, len);
    }
}

#define CFX_FILE_READ_MAX (256 * 1024)

int cfx_read_file_text(const char* path, uint8_t* out, size_t outlen, enum cfx_str_format mode) {
    FILE* f = fopen(path, "rb");
    if (!f) return -1;

    char buf[CFX_FILE_READ_MAX + 1];
    size_t n = fread(buf, 1, CFX_FILE_READ_MAX, f);
    int ferr = ferror(f);
    fclose(f);

    if (ferr || n == 0) {
        cfx_memzero_s(buf, sizeof(buf));
        return -1;
    }
    buf[n] = '\0';

    /* trim leading whitespace */
    char* start = buf;
    while (*start && isspace((unsigned char)*start)) start++;

    /* trim trailing whitespace */
    char* end = start + strlen(start);
    while (end > start && isspace((unsigned char)end[-1])) end--;
    *end = '\0';

    int result = cfx_parse_str(start, out, outlen, mode);

    cfx_memzero_s(buf, sizeof(buf));
    return result;
}

int cfx_read_file_bin(const char* path, uint8_t* out, size_t outlen) {
    FILE* f = fopen(path, "rb");
    if (!f) return -1;

    size_t n = fread(out, 1, outlen, f);
    int ferr = ferror(f);
    fclose(f);

    if (ferr || n != outlen) {
        cfx_memzero_s(out, outlen);
        return -1;
    }
    return 0;
}

static int is_base64_char(int c) {
    return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') ||
           (c >= '0' && c <= '9') || c == '+' || c == '/' || c == '=';
}

enum cfx_str_format cfx_detect_file_format(const char* path, size_t* out_len) {
    FILE* f = fopen(path, "rb");
    if (!f) return CFX_STR_FMT_BINARY;

    /* get file size */
    fseek(f, 0, SEEK_END);
    long fsize = ftell(f);
    fseek(f, 0, SEEK_SET);

    if (fsize <= 0 || fsize > CFX_FILE_READ_MAX) {
        fclose(f);
        if (out_len) *out_len = (fsize > 0) ? (size_t)fsize : 0;
        return CFX_STR_FMT_BINARY;
    }

    /* read file */
    char buf[CFX_FILE_READ_MAX + 1];
    size_t n = fread(buf, 1, (size_t)fsize, f);
    fclose(f);
    if (n != (size_t)fsize) {
        cfx_memzero_s(buf, sizeof(buf));
        if (out_len) *out_len = 0;
        return CFX_STR_FMT_BINARY;
    }
    buf[n] = '\0';

    /* check for null bytes or non-printable chars */
    int has_nonprint = 0;
    for (size_t i = 0; i < n; i++) {
        unsigned char c = (unsigned char)buf[i];
        if (c == 0) {
            cfx_memzero_s(buf, sizeof(buf));
            if (out_len) *out_len = n;
            return CFX_STR_FMT_BINARY;
        }
        if (!isprint(c) && !isspace(c)) {
            has_nonprint = 1;
        }
    }

    if (has_nonprint) {
        cfx_memzero_s(buf, sizeof(buf));
        if (out_len) *out_len = n;
        return CFX_STR_FMT_BINARY;
    }

    /* trim whitespace for format detection */
    char* start = buf;
    while (*start && isspace((unsigned char)*start)) start++;
    char* end = start + strlen(start);
    while (end > start && isspace((unsigned char)end[-1])) end--;
    *end = '\0';
    size_t len = (size_t)(end - start);

    /* check if hex */
    if (cfx_looks_like_hex(start)) {
        const char* hex = start;
        if (hex[0] == '0' && (hex[1] == 'x' || hex[1] == 'X')) hex += 2;
        if (out_len) *out_len = strlen(hex) / 2;
        cfx_memzero_s(buf, sizeof(buf));
        return CFX_STR_FMT_HEX;
    }

    /* check if base64 */
    int all_b64 = 1;
    for (size_t i = 0; i < len && all_b64; i++) {
        if (!is_base64_char(start[i])) all_b64 = 0;
    }
    if (all_b64 && len > 0 && len % 4 == 0) {
        /* base64 decoded size: 3 bytes per 4 chars, minus padding */
        size_t decoded = (len / 4) * 3;
        if (len >= 1 && start[len - 1] == '=') decoded--;
        if (len >= 2 && start[len - 2] == '=') decoded--;
        if (out_len) *out_len = decoded;
        cfx_memzero_s(buf, sizeof(buf));
        return CFX_STR_FMT_BASE64;
    }

    if (out_len) *out_len = len;
    cfx_memzero_s(buf, sizeof(buf));
    return CFX_STR_FMT_ASCII;
}


int cfx_read_all_file(FILE* f, uint8_t** out, size_t* out_len) {
    uint8_t* buf;
    size_t cap;
    size_t len;

    if (!f || !out || !out_len) return -1;

    *out = NULL;
    *out_len = 0;

    buf = NULL;
    cap = 0;
    len = 0;

    for (;;) {
        size_t r;
        size_t newcap;
        uint8_t* tmp;

        if (len == cap) {
            newcap = cap ? (cap * 2) : 4096;

            if (newcap > CFX_FILE_READ_MAX) newcap = CFX_FILE_READ_MAX;
            if (newcap <= cap) {  /* cap reached but still more input -> fail */
                if (buf) { cfx_memzero_s(buf, cap); free(buf); }
                return -1;
            }

            tmp = (uint8_t*)realloc(buf, newcap);
            if (!tmp) {
                if (buf) { cfx_memzero_s(buf, cap); free(buf); }
                return -1;
            }
            buf = tmp;
            cap = newcap;
        }

        r = fread(buf + len, 1, cap - len, f);
        len += r;

        if (r == 0) {
            if (ferror(f)) {
                cfx_memzero_s(buf, cap);  /* or use len? */
                free(buf);
                return -1;
            }
            break; /* EOF */
        }
    }

    *out = buf;
    *out_len = len;
    return 0;
}

