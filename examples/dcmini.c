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

typedef struct {
    token_type_t t;
    char op;          // for T_OP
    char *num;        // malloc'd string for T_NUM
} token_t;

typedef struct {
    token_t* v;
    size_t n, cap;
} token_vec_t;

static void token_vec_init(token_vec_t* t) {
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

static token_vec_t tokenize(const char* s) {
    token_vec_t tv;
    token_vec_init(&tv);

    for (size_t i = 0; s[i]; ) {
        if (isspace((unsigned char)s[i])) {
            ++i;
            continue;
        }
        if (isdigit((unsigned char)s[i])) {
            size_t j = i;
            while (isdigit((unsigned char)s[j])) ++j;
            token_t tk = {.t = T_NUM, .op = 0, .num = strndup(s + i, j - i) };
            tv_push(&tv, tk);
            i = j;
            continue;
        }
        if (strchr("+-*/%^", s[i])) {
            token_t tk = {.t = T_OP, .op=s[i], .num = NULL};
            tv_push(&tv, tk);
            ++i;
            continue;
        }
        if (s[i]=='(') {
            token_t tk = {.t = T_LP};
            tv_push(&tv, tk);
            ++i;
            continue;
        }
        if (s[i]==')') {
            token_t tk = {.t = T_RP};
            tv_push(&tv, tk);
            ++i;
            continue;
        }
        fprintf(stderr,"syntax error near '%c'\n", s[i]); exit(1);
    }
    return tv;
}

#define PUSH_OP(tk) do{if(top==cap){cap=cap?cap*2:32;opstack=realloc(opstack,cap*sizeof(token_t));}opstack[top++]=(tk);}while(0)
#define POP_OP()(opstack[--top])

// convert infix tokens to RPN using shunting-yard
static token_vec_t to_rpn(const token_vec_t* in) {
    token_vec_t out = {0};
    token_t* opstack = NULL;
    size_t top = 0, cap = 0;

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
}

static void apply_op(cfx_big_t* out,
                     const cfx_big_t* A,
                     const cfx_big_t* B,
                     char op)
{
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
        cfx_big_mul(&tmp, B);
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
        // plain pow (no mod). if you want ^ to be powmod only when followed by % m, keep this.
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
                cfx_big_mul(&res, &base);
            }
            cfx_big_mul(&base,&base);
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

static void eval_rpn(const token_vec_t* rpn, cfx_big_t* result) {
    cfx_big_t* st = NULL;
    size_t n = 0, cap = 0;

    #define PUSH()do{if(n==cap){cap=cap?cap*2:32; st=realloc(st,cap*sizeof(cfx_big_t));}cfx_big_init(&st[n]);++n;}while(0)
    #define POP(dst)do{if(n==0){fprintf(stderr,"stack underflow\n");exit(1);}--n;cfx_big_copy((dst),&st[n]);cfx_big_free(&st[n]);}while(0)

    for (size_t i = 0; i < rpn->n; i++) {
        token_t tk = rpn->v[i];
        if (tk.t == T_NUM) {
            PUSH();
            cfx_big_from_str(&st[n - 1], tk.num);
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
}

int main(int argc, char** argv) {
    
    if (argc < 2) {
        fprintf(stderr,"usage: %s \"EXPR\"\n", argv[0]);
        fprintf(stderr,"ops: + - * / %% ^   parentheses ()   decimal integers\n");
        return 2;
    }

    token_vec_t tv = tokenize(argv[1]);
    token_vec_t rpn = to_rpn(&tv);

    cfx_big_t res;
    cfx_big_init(&res);
    eval_rpn(&rpn, &res);

    char *s = cfx_big_to_str(&res, NULL);
    printf("%s\n", s);
    
    free(tv.v);
    free(rpn.v);
    cfx_big_free(&res);
    free(s);
    
    return 0;
}
