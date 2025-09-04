// cfx_error.h
#ifndef CFX_ERROR_H
#define CFX_ERROR_H

#include <errno.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum cfx_err {
    CFX_OK = 0,
    /* Map to system errno where available (so strerror() works). */
    CFX_EINVAL    = EINVAL,
    CFX_ENOMEM    = ENOMEM,
#ifdef EOVERFLOW
    CFX_EOVERFLOW = EOVERFLOW,
#else
    CFX_EOVERFLOW = ERANGE,
#endif
    /* ---- Start of library-private range ---- */
    CFX_EBASE_PRIVATE = 100000,
    CFX_EPARSE    = CFX_EBASE_PRIVATE + 1,
} cfx_err_t;

#define CFX_IS_SYSTEM_ERROR(code) ((int)(code) >= 0 && (int)(code) < CFX_EBASE_PRIVATE)

const char* cfx_strerror(cfx_err_t e);

#ifdef __cplusplus
}
#endif
#endif