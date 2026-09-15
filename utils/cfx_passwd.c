#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cfx/argon2.h"
#include "cfx/hmac.h"
#include "cfx_utils_common.h"

#define MAX_OUTLEN 2048

int cfx_passwd_run(int argc, char** argv) {

    char pw[128];
    char info[128];
    int result = 0;
    size_t outlen = 12;

    for (int i = 1; i < argc; ++i) {
        if (strcmp("-l", argv[i]) == 0) {
            
            if (i + 1 >= argc) {
                return -1;
            }
            outlen = (size_t)strtoul(argv[i+1], NULL, 10);

            if ( (outlen == 0) || (outlen > MAX_OUTLEN) ) {
                fprintf(stderr, "outlen invalid\n");
                outlen = 12;
            }

        }
    }
    
    int pwlen = cfx_read_secret("enter password:", pw, sizeof(pw));
    if (pwlen < 0) {
        fprintf(stderr, "pw reading failed\n");
        result = -1;
        goto cleanup;        
    }
    cfx_read_visible("enter info:", info, sizeof(info));
    uint8_t salt[] = {0xC0, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xEE};

    uint8_t master_key[32];
    uint32_t m_cost = 65536;  /* 64 MB */
    uint32_t t_cost = 3;
    uint32_t p_cost = 4;
    int rc = cfx_argon2id(master_key, sizeof(master_key), pw, pwlen,
                salt, sizeof(salt), m_cost, t_cost, p_cost);
    

    if (rc != 0) {
        fprintf(stderr, "Error: argon2 failed\n");
        result = 1;
        goto cleanup;
    }


    cfx_hmac_sha256_ctx extract_ctx;
    cfx_hmac_sha256_init(&extract_ctx, master_key, sizeof(master_key));
    uint8_t hkdf_salt[] = {0x00, 0x01, 0x02, 0x03, 0x04, 0x05};
    cfx_hmac_sha256_update(&extract_ctx, hkdf_salt, sizeof(hkdf_salt));
    /* this is to prevent the output to be the same for different outlens: */
    cfx_hmac_sha256_update(&extract_ctx, (uint8_t*)&outlen, sizeof(outlen));
    uint8_t prk[32];
    cfx_hmac_sha256_final(&extract_ctx, prk);

    printf("prk: ");
    for (int k = 0; k < sizeof(prk); ++k) {
        printf("%02x ", prk[k]);
    }
    printf("\n");
    const char* a = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789!@#)(*?:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";
    const size_t alen = strlen(a);
    size_t upper_bound = alen;

    // expand
    const int hashlen = 32;  // sha256
    size_t written = 0;
    uint8_t *indices = (uint8_t*)malloc(outlen);
    uint8_t cnt = 0;
    size_t tlen = 0;
    uint8_t t[32];
    cfx_hmac_sha256_ctx expand_ctx;
    do {
        cfx_hmac_sha256_init(&expand_ctx, prk, sizeof(prk));
        cfx_hmac_sha256_update(&expand_ctx, t, tlen);
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

    char *outbuf = (char*)malloc(outlen+1);
    for (int i = 0; i < outlen; ++i) {
        outbuf[i] = a[indices[i]];
    }
    outbuf[outlen] = '\0';
    printf("%s\n", outbuf);

    free(outbuf);
    free(indices);
cleanup:
    return result;
}