/* ---- cfx_dc.c - Desktop calculator for big integers ---- */
#include "cfx/big.h"
#include "cfx/compat.h"
#include "cfx_cmd.h"

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>
#include <limits.h>

/* each token in the string is either a numerical, operator, left paren or right paren */
/* TODO: allow funciton calls */
typedef enum {
    T_NUM,
    T_BIN_OP,  /* binary */
    T_UN_OP,   /* unary */
    T_LP,
    T_RP
} token_type_t;

typedef enum {
    BASE_DEC,
    BASE_HEX,
    BASE_BIN,
    BASE_64
} base_t;

typedef enum {
    MODE_INFIX,
    MODE_RPN
} input_mode_t;

typedef enum {
    ADD,    /* '+' */
    SUB,    /* '-' */
    MUL,    /* '*' */
    DIV,    /* '/' */
    POW,    /* '**' */
    SHR,    /* '>>' */
    SHL,    /* '<<' */
    FAC,    /* '!' */
    MOD,    /* '%'*/
    XOR,    /* '^' */
    AND,    /* '&' */
    OR,     /* '|' */
    ROTR,   /* '>>>' */
    ROTL    /* '<<<' */
} op_t;

typedef struct {
    token_type_t t;
    op_t op;          /* for when token type = T_BIN_OP */
    char *num;        /* malloc'd string for when token type = T_NUM */
} token_t;

typedef struct {
    token_t* v;
    size_t n, cap;
} token_vec_t;

static void tv_init(token_vec_t* t) {
    memset(t, 0, sizeof(*t));
}

static void tv_push(token_vec_t* tv, token_t tk) {
    /* printf("tv push: tk.t: %d, tk.op: %c, tk.num: %s\n", tk.t, tk.op, tk.num); */
    if (tv->n == tv->cap) {
        tv->cap = tv->cap ? tv->cap * 2 : 32;
        tv->v = realloc(tv->v, tv->cap * sizeof(token_t));
    }
    tv->v[tv->n++] = tk;
}

static int precedence(op_t op) {
    switch(op) {
        case FAC: return 4;
        case POW: return 3;
        case MUL: case DIV: case MOD: return 2;
        case ADD: case SUB: return 1;
        default: return -1;
    }
}

static int right_assoc(op_t op) { return op == POW; }

static int is_num_char(char c, base_t b) {
    if (b == BASE_DEC) return (c >= '0' && c <= '9');
    /* if (b == BASE_BIN) return (c == '0' || c == '1'); */
    /* hex */
    return (c >= '0' && c <= '9') ||
           (c >= 'a' && c <= 'f') ||
           (c >= 'A' && c <= 'F');
}

typedef struct {
    const char *str;
    op_t        op;
} op_map_t;

/* important: longest operator first */
static const op_map_t op_map[] = {
    { ">>>", ROTR },
    { "<<<", ROTL },
    { "**",  POW  },
    { ">>",  SHR  },
    { "<<",  SHL  },
    { "^",   XOR  },
    { "&",   AND  },
    { "|",   OR   },
    { "+",   ADD  },
    { "-",   SUB  },
    { "*",   MUL  },
    { "/",   DIV  },
    { "%",   MOD  },
};

/* returns num characters decodes, 0 if no match */
static int to_bin_op(op_t* op, const char* s) {
    for (size_t i = 0; i < sizeof(op_map)/sizeof(op_map[0]); ++i) {
        const char* opstr = op_map[i].str;
        size_t len = strlen(opstr);

        if (strncmp(s, opstr, len) == 0) {
            if (op) *op = op_map[i].op;
            return (int)len;
        }
    }
    return 0;
}

#define TOK_DBG 0
static token_vec_t tokenize(const char* s, base_t inb) {
    token_vec_t tv;
    tv_init(&tv);
    size_t len = strlen(s);
    for (size_t i = 0; s[i]; ) {

        if (isspace((unsigned char)s[i])) {
            ++i;
            continue;
        }
        /* try to parse a number */
        #if TOK_DBG
        printf("TOK i=%zu char='%c'\n", i, s[i]);
        #endif
        size_t ncons = 0;
        cfx_big_t tmp;
        cfx_big_init(&tmp);

        if (cfx_big_scan_num_n(&tmp, (const uint8_t*)s + i, len - i, &ncons) == 0) {
            token_t tk = {.t = T_NUM, .num=cfx_strndup(s+i, ncons)};
            tv_push(&tv, tk);
            i += ncons;
            const char* st = cfx_big_dec_alloc(&tmp, NULL);
            #if TOK_DBG
            printf("NUM %s\n", st);
            #endif
            cfx_big_free(&tmp);
            continue;
        }
        cfx_big_free(&tmp);
        /* binary op*/
        op_t op;
        int oplen = to_bin_op(&op, &s[i]);
        if (oplen > 0) {
            token_t tk = {
                .t   = T_BIN_OP,
                .op  = op,
                .num = NULL
            };
            #if TOK_DBG
            printf("- BINOP %.*s\n", oplen, &s[i]);
            #endif
            tv_push(&tv, tk);
            i += (size_t)oplen;
            continue;
        }
        if (strchr("!", s[i])) { /* unary op */
            token_t tk = {
                .t = T_UN_OP,
                .op = s[i],
                .num = NULL
            };
            #if TOK_DBG
            printf("+ UNOP !\n");
            #endif
            tv_push(&tv, tk);
            ++i;
            continue;
        }
        if (s[i] == '(') {
            token_t tk = {.t = T_LP};
            tv_push(&tv, tk);
            ++i;
            continue;
        }
        if (s[i] == ')') {
            token_t tk = {.t = T_RP};
            tv_push(&tv, tk);
            ++i;
            continue;
        }
        fprintf(stderr,"syntax error near '%c'\n", s[i]); exit(1);
    }
    return tv;
}

static token_vec_t tokenize_rpn(const char* s, base_t inb) {
    (void)inb; /* todo! */
    token_vec_t tv;
    tv_init(&tv);
    const char* p = s;
    while (*p) {
        while (isspace((unsigned char)*p)) ++p;
        if (!*p) break;
        const char* start = p;
        while (*p && !isspace((unsigned char)*p)) ++p;
        size_t len = (size_t)(p - start);
        op_t binop;
        if ((size_t)to_bin_op(&binop, start) == len) {
            token_t tk = {.t = T_BIN_OP, .op = binop, .num = NULL};
            tv_push(&tv, tk);
        } else if (len == 1 && strchr("!", start[0])) {
            token_t tk = {.t = T_UN_OP, .op = start[0], .num = NULL};
            tv_push(&tv, tk);
        } else {
            /* treat as number; validation happens in cfx parse anyway */
            token_t tk = {.t = T_NUM, .op = 0, .num = cfx_strndup(start, len)};
            tv_push(&tv, tk);
        }
    }
    return tv;
}

/* convert infix tokens to RPN using shunting-yard */
static token_vec_t to_rpn(const token_vec_t* in) {
    token_vec_t out = {0};
    token_t* opstack = NULL;
    size_t top = 0, cap = 0;

    #define PUSH_OP(tk) \
        do { \
            if (top == cap) { \
                cap = cap ? (cap*2) : 32; \
                opstack = realloc(opstack, cap * sizeof(token_t)); \
            } \
            opstack[top++] = (tk); \
        } while (0)

    #define POP_OP() (opstack[--top])

    for (size_t i = 0; i < in->n; i++) {
        token_t tk = in->v[i];
        if (tk.t == T_NUM) {
            tv_push(&out, tk);
            continue;
        }
        if (tk.t == T_BIN_OP) {
            while (top > 0 && opstack[top-1].t == T_BIN_OP) {
                char o2 = opstack[top-1].op;
                char o1 = tk.op;
                if ( (right_assoc(o1) ? precedence(o1) < precedence(o2)
                                      : precedence(o1) <= precedence(o2)) ) {
                    tv_push(&out, POP_OP());
                } else break;
            }
            PUSH_OP(tk);
        } else if (tk.t == T_UN_OP) {
            PUSH_OP(tk);
        } else if (tk.t == T_LP) {
            PUSH_OP(tk);
        } else if (tk.t == T_RP) {
            int found = 0;
            while (top > 0) {
                token_t t2 = POP_OP();
                if (t2.t == T_LP) {
                    found = 1;
                    break;
                }
                tv_push(&out, t2);
            }
            if (!found) {
                fprintf(stderr,"mismatched parentheses\n");
                exit(1);
            }
        }
    }
    while (top > 0) {
        token_t t2 = POP_OP();
        if (t2.t == T_LP) {
            fprintf(stderr,"mismatched parentheses\n");
            exit(1);
        }
        tv_push(&out, t2);
    }
    free(opstack);
    return out;

    #undef PUSH_OP
    #undef POP_OP
}

static void apply_un_op(cfx_big_t* out, const cfx_big_t* A, char op) {
    if (op == '!') {
        if (cfx_big_is_zero(A)) {
            cfx_big_assign_one(out);
            return;
        }
        if (A->n > 1) {
            char* s = cfx_big_dec_alloc(A, NULL);
            printf("You are trying to take the factorial of %s, rethink things...\n", s);
            free(s);
            return;
        }
        cfx_limb_t n = A->limb[0];
        cfx_vec_t primes = cfx_sieve_primes(n);
        cfx_fac_t fac;
        cfx_fac_init(&fac);
        cfx_fac_factorial(&fac, n, &primes);
        cfx_big_from_fac_faster(out, &fac);
        cfx_fac_free(&fac);
        cfx_vec_free(&primes);
    }
}

static void apply_bin_op(cfx_big_t* out, const cfx_big_t* A, const cfx_big_t* B, op_t op) {

    switch (op) {
        case ADD: {
            cfx_big_t tmp;
            cfx_big_init(&tmp);
            cfx_big_copy(&tmp, A);
            cfx_big_add_eq(&tmp, B);
            cfx_big_move(out, &tmp);
            return;
        }
        case SUB: {
            cfx_big_t tmp;
            cfx_big_init(&tmp);
            cfx_big_copy(&tmp, A);
            cfx_big_sub_eq(&tmp, B);
            cfx_big_move(out, &tmp);
            return;
        }
        case MUL: {
            cfx_big_t tmp;
            cfx_big_init(&tmp);
            cfx_big_copy(&tmp, A);
            cfx_big_mul_auto(&tmp, B);
            cfx_big_move(out, &tmp);
            return;
        }
        case DIV: {
            cfx_big_t q, r;
            cfx_big_init(&q);
            cfx_big_init(&r);
            cfx_big_divrem(&q, &r, A, B);
            cfx_big_move(out, &q);
            cfx_big_free(&r);
            return;
        }
        case MOD: {
            cfx_big_t q, r;
            cfx_big_init(&q);
            cfx_big_init(&r);
            cfx_big_divrem(&q, &r, A, B);
            cfx_big_move(out, &r);
            cfx_big_free(&q);
            return;
        }
        case POW: {
            /* plain pow (no mod).  */
            /* Here we implement a^b with b as non-negative. */
            cfx_big_t e;
            cfx_big_init(&e);
            cfx_big_copy(&e,B);
            cfx_big_t res, base;
            cfx_big_init(&res);
            cfx_big_init(&base);
            cfx_big_from_limb(&res, 1);
            cfx_big_copy(&base, A);
            cfx_big_t zero, two, q, r;
            cfx_big_init(&zero);
            cfx_big_init(&two);
            cfx_big_init(&q);
            cfx_big_init(&r);
            cfx_big_from_limb(&zero, 0);
            cfx_big_from_limb(&two, 2);
            while (cfx_big_cmp(&e, &zero) > 0) {
                cfx_big_divrem(&q, &r, &e, &two);
                if (!cfx_big_is_zero(&r)) {
                    cfx_big_mul_auto(&res, &base);
                }
                cfx_big_mul_auto(&base,&base);
                cfx_big_copy(&e,&q);
            }
            cfx_big_copy(out,&res);

            cfx_big_free(&e);
            cfx_big_free(&res);
            cfx_big_free(&base);
            cfx_big_free(&zero);
            cfx_big_free(&two);
            cfx_big_free(&q);
            cfx_big_free(&r);
            return;
        }
        case OR: {
            cfx_big_or(out, A, B);
            return;
        }
        case AND: {
            cfx_big_and(out, A, B);
            return;
        }
        case XOR: {
            cfx_big_xor(out, A, B);
            return;
        }
        case SHL: {
            if (B->n == 1) {
                cfx_big_shl_bits(out, A, (int)B->limb[0]);
            } else {
                const char* s = cfx_big_dec_alloc(B, NULL);
                fprintf(stderr, "you're trying to shift left by %s\n, rethink things...\n", s);
            }
            return;
        }
        case SHR: {
            if (B->n == 1) {
                cfx_big_shr_bits(out, A, (int)B->limb[0]);
            } else {
                const char* s = cfx_big_dec_alloc(B, NULL);
                fprintf(stderr, "you're trying to shift right by %s\n, rethink things...\n", s);
            }
            return;
        }
        case ROTR: {
            if (B->n == 1 && B->limb[0] < UINT_MAX) {
                cfx_big_rotr(out, A, (unsigned)B->limb[0]);
            } else {
                const char* s = cfx_big_dec_alloc(B, NULL);
                fprintf(stderr, "you're trying to rotate right by %s\n, rethink things...\n", s);
            }
            return;
        }
        case ROTL: {
            if (B->n == 1&& B->limb[0] < UINT_MAX) {
                cfx_big_rotl(out, A, (unsigned)B->limb[0]);
            } else {
                const char* s = cfx_big_dec_alloc(B, NULL);
                fprintf(stderr, "you're trying to rotate left by %s\n, rethink things...\n", s);
            }
            return;
        }
        default:
            fprintf(stderr,"unknown binary op %d '%c'\n", (int)op, op); exit(1);
    }
}

static void eval_rpn(const token_vec_t* rpn, cfx_big_t* result, base_t inb) {
    cfx_big_t* st = NULL;
    size_t n = 0, cap = 0;

    #define PUSH() \
        do { \
            if (n == cap) { \
                cap = cap ? (cap * 2) : 32; \
                st = realloc(st, cap * sizeof(cfx_big_t)); \
            } \
            cfx_big_init(&st[n]); \
            ++n; \
        } while (0)

    #define POP(dst) \
        do { \
            if (n == 0) { fprintf(stderr, "stack underflow\n"); exit(1); } \
            --n; \
            cfx_big_copy((dst), &st[n]); \
            cfx_big_free(&st[n]); \
        } while (0)

    for (size_t i = 0; i < rpn->n; i++) {
        token_t tk = rpn->v[i];
        if (tk.t == T_NUM) {
            PUSH();
            cfx_big_from_str(&st[n-1], tk.num);
            free(tk.num);
        } else if (tk.t == T_BIN_OP) {
            cfx_big_t b, a, res;
            cfx_big_init(&a);
            cfx_big_init(&b);
            cfx_big_init(&res);
            POP(&b);
            POP(&a);
            apply_bin_op(&res, &a, &b, tk.op);
            cfx_big_free(&a);
            cfx_big_free(&b);
            PUSH();
            cfx_big_copy(&st[n - 1], &res);
            cfx_big_free(&res);
        } else if (tk.t == T_UN_OP) {
            cfx_big_t a, res;
            cfx_big_init(&a);
            cfx_big_init(&res);
            POP(&a);
            apply_un_op(&res, &a, tk.op);
            cfx_big_free(&a);
            PUSH();
            cfx_big_copy(&st[n - 1], &res);
            cfx_big_free(&res);
        }
    }
    if (n != 1) {
        fprintf(stderr,"invalid expression\n");
        exit(1);
    }
    cfx_big_copy(result, &st[0]);
    cfx_big_free(&st[0]); free(st);

    #undef PUSH
    #undef POP
}

static void usage(const char* prog) {
    fprintf(stderr,
        "usage: %s [opts] EXPR...\n"
        "  input base : -id (dec, default) | -ix (hex)\n" /* | -ib (bin)\n" */
        "  output base: -od (dec, default) | -ox (hex) | -ob (bin) | -ob64 (base64) \n"
        "  mode       : -rpn (treat input as Reverse Polish Notation)\n"
        "examples:\n"
        "  %s \"(1237+1238)^3298 %% 237\"\n"
        "  %s -rpn 1237 1238 + 3298 ^ 237 %%\n"
        "  %s -ix -ox FF * F  (hex in, hex out)\n",
        prog, prog, prog, prog);
}



static char* join_args(int argc, char** argv, int start) {
    size_t len = 0;
    for (int i = start; i < argc; i++) len += strlen(argv[i]) + 1;
    char* s = malloc(len + 1);
    if (!s) return NULL;
    s[0] = '\0';
    for (int i = start; i < argc; i++) {
        if (i > start) strcat(s, " ");
        strcat(s, argv[i]);
    }
    return s;
}

static int parse_flags(int argc, char** argv, base_t* inb, base_t* outb,
                       input_mode_t* mode, int* expr_index) {
    *inb = BASE_DEC;
    *outb = BASE_DEC;
    *mode = MODE_INFIX;

    int i = 1;
    for (; i < argc; ++i) {
        const char* a = argv[i];
        if (a[0] != '-') break;
        if (strcmp(a,"-id") == 0) *inb = BASE_DEC;
        else if (strcmp(a,"-ix") == 0) *inb = BASE_HEX;
        /* else if (strcmp(a,"-ib") == 0) *inb = BASE_BIN; */
        else if (strcmp(a,"-od") == 0) *outb = BASE_DEC;
        else if (strcmp(a,"-ox") == 0) *outb = BASE_HEX;
        else if (strcmp(a,"-ob") == 0) *outb = BASE_BIN;
        else if (strcmp(a,"-ob64") == 0) *outb = BASE_64;
        else if (strcmp(a,"-rpn") == 0) *mode = MODE_RPN;
        else if (strcmp(a,"-h") == 0 || strcmp(a,"--help") == 0) { usage(argv[0]); return 0; }
        else { fprintf(stderr,"unknown option: %s\n", a); usage(argv[0]); return 0; }
    }
    if (i >= argc) { usage(argv[0]); return 0; }
    *expr_index = i;
    return 1;
}



int cfx_dc_run(int argc, char** argv) {

    base_t inb, outb;
    input_mode_t mode;
    int idx;

    if (!parse_flags(argc, argv, &inb, &outb, &mode, &idx)) return 2;

    char* expr = join_args(argc, argv, idx);
    if (!expr) { fprintf(stderr,"OOM\n"); return 1; }

    token_vec_t tv = (mode == MODE_RPN) ? tokenize_rpn(expr, inb)
                                        : tokenize(expr, inb);
    token_vec_t rpn = (mode == MODE_RPN) ? tv : to_rpn(&tv);

    cfx_big_t res;
    cfx_big_init(&res);
    eval_rpn(&rpn, &res, inb);

    char* out = NULL;
    char* prefix = NULL;
    size_t digits;
    switch (outb) {
        case BASE_DEC:
            out = cfx_big_dec_alloc(&res, &digits);
            prefix = "";
            break;
        case BASE_HEX:
            out = cfx_big_hex_alloc(&res, &digits);
            prefix = "0x";
            break;
        case BASE_BIN:
            out = cfx_big_bin_alloc(&res, &digits);
            prefix = "0b";
            break;
        case BASE_64:
            out = cfx_big_b64_alloc(&res, &digits);
            prefix = "b64:";
            break;
    }
    printf("%s%s\n", prefix, out);

    if (mode != MODE_RPN) free(tv.v);  /* tv == rpn when MODE_RPN */
    free(rpn.v);
    cfx_big_free(&res);
    free(out);
    free(expr);
    return 0;
}
