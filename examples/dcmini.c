// ---- dcmini.c ----
#include "cfx/cfx.h"

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>


typedef enum {
    T_NUM,
    T_OP,
    T_LP,
    T_RP
} token_type_t;

typedef enum {
    BASE_DEC,
    BASE_HEX,
    BASE_BIN
} base_t;

typedef enum {
    MODE_INFIX,
    MODE_RPN
} input_mode_t;

typedef struct {
    token_type_t t;
    char op;          // for when token type = T_OP
    char *num;        // malloc'd string for when token type = T_NUM
} token_t;

typedef struct {
    token_t* v;
    size_t n, cap;
} token_vec_t;

static void tv_init(token_vec_t* t) {
    memset(t, 0, sizeof(*t));
}

static void tv_push(token_vec_t* tv, token_t tk) { 
    if (tv->n == tv->cap) {
        tv->cap = tv->cap ? tv->cap * 2 : 32;
        tv->v = realloc(tv->v, tv->cap * sizeof(token_t));
    }
    tv->v[tv->n++] = tk;
}

static int precedence(char op) {
    switch(op) {
        case '^': return 3;
        case '*': case '/': case '%': return 2;
        case '+': case '-': return 1;
        default: return -1;
    }
}

static int right_assoc(char op) { return op=='^'; }

static int is_num_char(char c, base_t b){
    if(b == BASE_DEC) return (c>='0' && c<='9');
    // if(b == BASE_BIN) return (c=='0' || c=='1');
    /* hex */
    return (c>='0' && c<='9') ||
           (c>='a' && c<='f') ||
           (c>='A' && c<='F');
}

static token_vec_t tokenize(const char* s, base_t inb) {
    token_vec_t tv;
    tv_init(&tv);
    for (size_t i = 0; s[i]; ) {
        if (isspace((unsigned char)s[i])) {
            ++i;
            continue;
        }
        if (is_num_char(s[i], inb)) {
            size_t j = i;
            while (is_num_char(s[j], inb)) ++j;
            token_t tk = {.t = T_NUM, .op = 0, .num = strndup(s+i, j-i) };
            tv_push(&tv, tk); i = j;
            continue;
        }
        if (strchr("+-*/%^", s[i])) {
            token_t tk = {.t = T_OP, .op = s[i], .num = NULL};
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

static token_vec_t tokenize_rpn(const char* s, base_t inb){
    (void)inb; // todo!
    token_vec_t tv;
    tv_init(&tv);
    const char* p = s;
    while (*p) {
        while (isspace((unsigned char)*p)) ++p;
        if (!*p) break;
        const char* start = p;
        while (*p && !isspace((unsigned char)*p)) ++p;
        size_t len = (size_t)(p - start);
        if (len == 1 && strchr("+-*/%^", start[0])) {
            token_t tk = {.t = T_OP, .op = start[0], .num = NULL};
            tv_push(&tv, tk);
        } else {
            /* treat as number; validation happens in cfx parse anyway */
            token_t tk = {.t = T_NUM, .op = 0, .num = strndup(start, len)};
            tv_push(&tv, tk);
        }
    }
    return tv;
}

// convert infix tokens to RPN using shunting-yard
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
        } while(0)

    #define POP_OP() (opstack[--top])
    
    for (size_t i = 0; i < in->n; i++) {
        token_t tk = in->v[i];
        if (tk.t == T_NUM) {
            tv_push(&out, tk);
            continue;
        }
        if (tk.t == T_OP) {
            while (top>0 && opstack[top-1].t==T_OP) {
                char o2=opstack[top-1].op, o1=tk.op;
                if ( (right_assoc(o1) ? precedence(o1) < precedence(o2)
                                      : precedence(o1) <= precedence(o2)) ) {
                    tv_push(&out, POP_OP());
                } else break;
            }
            PUSH_OP(tk);
        } else if (tk.t==T_LP) {
            PUSH_OP(tk);
        } else if (tk.t==T_RP) {
            int found=0;
            while (top>0) {
                token_t t2=POP_OP();
                if (t2.t==T_LP) { found=1; break; }
                tv_push(&out, t2);
            }
            if (!found) { fprintf(stderr,"mismatched parentheses\n"); exit(1); }
        }
    }
    while (top>0) {
        token_t t2=POP_OP();
        if (t2.t==T_LP) { fprintf(stderr,"mismatched parentheses\n"); exit(1); }
        tv_push(&out, t2);
    }
    free(opstack);
    return out;

    #undef PUSH_OP
    #undef POP_OP
}

static void apply_op(cfx_big_t* out, const cfx_big_t* A, const cfx_big_t* B, char op) {

    if (op=='+') { 
        cfx_big_t tmp;
        cfx_big_init(&tmp);
        cfx_big_copy(&tmp, A);
        cfx_big_add(&tmp, B);
        cfx_big_move(out, &tmp);
        return;
    }
    if (op=='-') { 
        cfx_big_t tmp;
        cfx_big_init(&tmp);
        cfx_big_copy(&tmp, A);
        cfx_big_sub(&tmp, B);
        cfx_big_move(out, &tmp);
        return;
    }
    if (op=='*') {
        cfx_big_t tmp;
        cfx_big_init(&tmp);
        cfx_big_copy(&tmp, A);
        cfx_big_mul_auto(&tmp, B);
        cfx_big_move(out, &tmp);
        return;
    }
    if (op=='/') {
        cfx_big_t q, r;
        cfx_big_init(&q);
        cfx_big_init(&r);
        cfx_big_divrem(&q, &r, A, B);
        cfx_big_move(out, &q);
        cfx_big_free(&r);
        return;
    }
    if (op == '%') {
        cfx_big_t q, r;
        cfx_big_init(&q);
        cfx_big_init(&r);
        cfx_big_divrem(&q, &r, A, B);
        cfx_big_move(out, &r);
        cfx_big_free(&q);
        return;
    }
    if (op == '^') {
        // plain pow (no mod). 
        // Here we implement a^b with b as non-negative.
        // For huge exponents you’ll almost always want (a^b) % m written with % in expr.
        cfx_big_t e;
        cfx_big_init(&e);
        cfx_big_copy(&e,B);
        cfx_big_t res, base;
        cfx_big_init(&res);
        cfx_big_init(&base);
        cfx_big_from_u64(&res, 1);
        cfx_big_copy(&base, A);
        cfx_big_t zero, two, q, r; 
        cfx_big_init(&zero);
        cfx_big_init(&two);
        cfx_big_init(&q);
        cfx_big_init(&r);
        cfx_big_from_u64(&zero, 0);
        cfx_big_from_u64(&two, 2);
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
    fprintf(stderr,"unknown op '%c'\n", op); exit(1);
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
            switch (inb) {
                case BASE_DEC: cfx_big_from_str(&st[n-1], tk.num); break;
                case BASE_HEX: cfx_big_from_hex(&st[n-1], tk.num); break;
                case BASE_BIN: break; /* TODO: cfx_big_from_bin(&st[n-1], tk.num); */ 
                default: break;
            }
            free(tk.num);
        } else if (tk.t == T_OP) {
            cfx_big_t b, a, res;
            cfx_big_init(&a);
            cfx_big_init(&b);
            cfx_big_init(&res);
            POP(&b);
            POP(&a);
            apply_op(&res, &a, &b, tk.op);
            cfx_big_free(&a);
            cfx_big_free(&b);
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
        "  input base : -id (dec, default) | -ix (hex)\n" // | -ib (bin)\n"
        "  output base: -od (dec, default) | -ox (hex) | -ob (bin)\n"
        "  mode       : -rpn (treat input as Reverse Polish Notation)\n"
        "examples:\n"
        "  %s \"(1237+1238)^3298 %% 237\"\n"
        "  %s -rpn 1237 1238 + 3298 ^ 237 %%\n"
        "  %s -ix -ox FF * F  (hex in, hex out)\n",
        prog, prog, prog, prog);
}



static char* join_args(int argc, char** argv, int start) {
    size_t len = 0;
    for(int i=start;i<argc;i++) len += strlen(argv[i]) + 1;
    char* s = malloc(len + 1);
    if(!s) return NULL;
    s[0] = '\0';
    for(int i=start;i<argc;i++){
        if(i>start) strcat(s, " ");
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
        if(a[0] != '-') break;
        if(strcmp(a,"-id")==0) *inb = BASE_DEC;
        else if(strcmp(a,"-ix")==0) *inb = BASE_HEX;
        // else if(strcmp(a,"-ib")==0) *inb = BASE_BIN;
        else if(strcmp(a,"-od")==0) *outb = BASE_DEC;
        else if(strcmp(a,"-ox")==0) *outb = BASE_HEX;
        else if(strcmp(a,"-ob")==0) *outb = BASE_BIN;
        else if(strcmp(a,"-rpn")==0) *mode = MODE_RPN;
        else if(strcmp(a,"-h")==0 || strcmp(a,"--help")==0){ usage(argv[0]); return 0; }
        else { fprintf(stderr,"unknown option: %s\n", a); usage(argv[0]); return 0; }
    }
    if (i >= argc) { usage(argv[0]); return 0; }
    *expr_index = i;
    return 1;
}



int main(int argc, char** argv) {
    
    base_t inb, outb;
    input_mode_t mode;
    int idx;
    
    if(!parse_flags(argc, argv, &inb, &outb, &mode, &idx)) return 2;

    char* expr = join_args(argc, argv, idx);
    if(!expr){ fprintf(stderr,"OOM\n"); return 1; }

    token_vec_t tv = (mode == MODE_RPN) ? tokenize_rpn(expr, inb)
                                        : tokenize(expr, inb);
    token_vec_t rpn = (mode == MODE_RPN) ? tv : to_rpn(&tv);

    cfx_big_t res;
    cfx_big_init(&res);
    eval_rpn(&rpn, &res, inb);

    char* out = NULL;
    size_t digits;
    switch(outb){
        case BASE_DEC: out = cfx_big_to_str(&res, &digits); break;
        case BASE_HEX: out = cfx_big_to_hex(&res, &digits); break;
        case BASE_BIN: out = cfx_big_to_bin(&res, &digits); break;
    }
    puts(out);
    printf("digits: %zu\n", digits);

    if(mode != MODE_RPN) free(tv.v);  // tv == rpn when MODE_RPN
    free(rpn.v);
    cfx_big_free(&res);
    free(out);
    free(expr);
    return 0;
}
