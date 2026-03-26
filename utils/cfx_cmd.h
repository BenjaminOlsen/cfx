/* cfx_cmd.h - Run function declarations for cfx utilities */

#ifndef CFX_UTILS_H
#define CFX_UTILS_H

int cfx_aead_run(int argc, char** argv);
int cfx_choose_run(int argc, char** argv);
int cfx_dc_run(int argc, char** argv);
int cfx_div_run(int argc, char** argv);
int cfx_chacha20_run(int argc, char** argv);
int cfx_factor_run(int argc, char** argv);
int cfx_factorial_run(int argc, char** argv);
int cfx_fread_run(int argc, char** argv);
int cfx_hash_run(int argc, char** argv);
int cfx_hmac_run(int argc, char** argv);
int cfx_primegen_run(int argc, char** argv);
int cfx_keygen_run(int argc, char** argv);
int cfx_mac_run(int argc, char** argv);
int cfx_modexp_run(int argc, char** argv);
int cfx_pow_run(int argc, char** argv);
int cfx_isprime_run(int argc, char** argv);
int cfx_primes_run(int argc, char** argv);
int cfx_primes_near_pow2_run(int argc, char** argv);
int cfx_rand_run(int argc, char** argv);
int cfx_ulam_spiral_run(int argc, char** argv);
int cfx_pi_run(int argc, char** argv);
int cfx_x25519_run(int argc, char** argv);
int cfx_ed25519_run(int argc, char** argv);
int cfx_sign_run(int argc, char** argv);
int cfx_verify_run(int argc, char** argv);
int cfx_argon2_run(int argc, char** argv);
int cfx_xgcd_run(int argc, char** argv);
int cfx_bge_run(int argc, char** argv);
int cfx_store_run(int argc, char** argv);
int cfx_randart_run(int argc, char** argv);

typedef struct {
    const char* name;
    const char* description;
    int (*run)(int, char**);
} cfx_cmd_t;

extern const cfx_cmd_t cfx_commands[];
extern const int cfx_commands_count;

#endif
