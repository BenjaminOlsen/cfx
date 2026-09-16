#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cfx/hmac.h"
#include "cfx/argon2.h"
#include "cfx/memory.h"
#include "cfx_utils_common.h"

#define MAX_OUTLEN 2048

int cfx_passwd_run(int argc, char** argv) {

    char pw[128];
    char info[128];
    int result = 0;
    size_t outlen = 20;
    uint8_t user_salt[32] = {0xC0, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xEE};

    uint8_t master_key[32];
    const uint32_t m_cost = 65536;  /* 64 MB */
    const uint32_t t_cost = 5;
    const uint32_t p_cost = 4;

    static const uint8_t hkdf_salt[] = "cfx_passwd hkdf v1";
    uint8_t prk[32];

    uint8_t *indices = NULL;
    uint8_t *outbuf = NULL;

    for (int i = 1; i < argc; ++i) {
        char *arg = argv[i];
        if ((strcmp("-l", arg) == 0) || (strcmp("--length", arg) == 0)) {
            
            if (i + 1 >= argc) {
                fprintf(stderr, "need a length argument with -l\n");
                return 1;
            }

            char *end;
            unsigned long value = strtoul(argv[i+1], &end, 10);
            if (value == 0 || value > MAX_OUTLEN || *end != '\0' ) {
                fprintf(stderr, "outlen should be between 1 and %d\n", MAX_OUTLEN);
                return 1;
            }
            outlen = (size_t)value;
        } else if (strcmp("-s", arg) == 0 || strcmp("--salt", arg) == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "need a value for salt!\n");
                return 1;
            }
            int chars_read = cfx_parse_str(argv[i+1], user_salt, sizeof(user_salt), CFX_STR_FMT_AUTO);
            if (chars_read < 0) {
                fprintf(stderr, "problem reading a value for salt!\n");
                return 1;
            }
            printf("read salt: ");
            for (int k = 0; k < chars_read; ++k) {
                printf("%02x ", user_salt[k]);
            }
            printf("\n");

        }
    }
    
    int pwlen = cfx_read_secret("enter password:", pw, sizeof(pw));
    if (pwlen < 0) {
        fprintf(stderr, "pw reading failed\n");
        result = 1;
        goto cleanup;        
    }
    int infolen = cfx_read_visible("enter info:", info, sizeof(info));
    
    int rc = cfx_argon2id(master_key, sizeof(master_key), pw, pwlen,
                user_salt, sizeof(user_salt), m_cost, t_cost, p_cost);

    if (rc != 0) {
        fprintf(stderr, "Error: argon2 failed\n");
        result = 1;
        goto cleanup;
    }


    cfx_hmac_sha256_ctx extract_ctx;
    cfx_hmac_sha256_init(&extract_ctx, hkdf_salt, sizeof(hkdf_salt));
    cfx_hmac_sha256_update(&extract_ctx, master_key, sizeof(master_key));
    /* this is to prevent the output to be the same for different outlens: */
    cfx_hmac_sha256_final(&extract_ctx, prk);

    // expand
    size_t written = 0;
    indices = (uint8_t*)malloc(outlen);
    if (!indices) {
        result = 1;
        goto cleanup;
    }
    uint8_t cnt = 0;
    size_t tlen = 0;
    uint8_t t[32];   /* rfc5869 T(n)*/
    cfx_hmac_sha256_ctx expand_ctx;
    const char* alphabet = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789!@#*?:.";
    const size_t upper_bound = strlen(alphabet);
    do {
        cfx_hmac_sha256_init(&expand_ctx, prk, sizeof(prk));
        cfx_hmac_sha256_update(&expand_ctx, t, tlen);
        cfx_hmac_sha256_update(&expand_ctx, info, infolen);
        cfx_hmac_sha256_update(&expand_ctx, (uint8_t*)&outlen, sizeof(outlen));
        if (cnt == UINT8_MAX) {
            /* No more blocks available. */
            printf("no more block available\n");
            goto cleanup;  /* Set an error result first. */
        }
        ++cnt;
        cfx_hmac_sha256_update(&expand_ctx, &cnt, sizeof(cnt));
        cfx_hmac_sha256_final(&expand_ctx, t);
        tlen = 32;

        for (int j = 0; j < sizeof(t) && written < outlen; ++j) {
            if (t[j] < (uint8_t)upper_bound) {
                indices[written] = t[j];
                ++written;
            }
        }
    } while (written < outlen);

    outbuf = (char*)malloc(outlen+1);
    if (!outbuf) {
        result = 1;
        goto cleanup;
    }
    for (int i = 0; i < outlen; ++i) {
        outbuf[i] = alphabet[indices[i]];
    }
    outbuf[outlen] = '\0';
    printf("%s\n", outbuf);

cleanup:
    cfx_memzero_s(pw, sizeof(pw));
    cfx_memzero_s(prk, sizeof(prk));
    cfx_memzero_s(master_key, sizeof(master_key));
    if (outbuf) {
        cfx_memzero_s(outbuf, outlen+1);
        free(outbuf);
    }
    if (indices) {
        cfx_memzero_s(indices, outlen);
        free(indices);
    }
    cfx_memzero_s(t, sizeof(t));

    return result;
}