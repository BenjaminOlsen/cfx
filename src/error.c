#include "cfx/error.h"
#include <string.h>

const char* cfx_strerror(cfx_err_t e) {
    if (CFX_IS_SYSTEM_ERROR(e)) {
        const char* s = strerror((int)e);
        return s ? s : "unknown system error";
    }
    switch (e) {
        case CFX_EPARSE:    return "parse error";
        default:            return "unknown cfx error";
    }
}