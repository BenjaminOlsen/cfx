/* cfx.c - Unified dispatcher for cfx utilities */

#include <stdio.h>
#include <string.h>
#include "cfx_cmd.h"

const cfx_cmd_t cfx_commands[] = {
    {"bge",             "Encrypted secrets store (BGE)",                    cfx_bge_run},
    {"choose",          "Compute binomial coefficients C(n,k)",             cfx_choose_run},
    {"dc",              "Desk calculator for big integers",                 cfx_dc_run},
    {"div",             "Divide two big integers",                          cfx_div_run},
    {"chacha20",        "Encrypt text with ChaCha20",                       cfx_chacha20_run},
    {"factor",          "Prime factorization of integers",                  cfx_factor_run},
    {"factorial",       "Compute n! (factorial)",                           cfx_factorial_run},
    {"fread",           "Read big integer from file",                       cfx_fread_run},
    {"hash",            "Compute cryptographic hash (SHA256/SHA3/BLAKE2)",  cfx_hash_run},
    {"primegen",        "Generate random N-bit prime",                      cfx_primegen_run},
    {"keygen",          "Generate random key bytes",                        cfx_keygen_run},
    {"mac",             "Compute Poly1305 MAC",                             cfx_mac_run},
    {"modexp",          "Modular exponentiation (base^exp mod m)",          cfx_modexp_run},
    {"pow",             "Compute n^p",                                      cfx_pow_run},
    {"isprime",      "Miller-Rabin primality test",                      cfx_isprime_run},
    {"primes",          "List primes up to n (sieve)",                      cfx_primes_run},
    {"primes_near_pow2","Find primes near powers of 2",                     cfx_primes_near_pow2_run},
    {"rand",            "Generate random bytes",                            cfx_rand_run},
    {"ulam_spiral",     "Generate Ulam spiral image",                       cfx_ulam_spiral_run},
    {"pi",              "Compute digits of pi",                             cfx_pi_run},
    {"x25519",          "X25519 key exchange (RFC 7748)",                   cfx_x25519_run},
    {"ed25519",         "Ed25519 digital signatures (RFC 8032)",            cfx_ed25519_run},
    {"sign",            "Sign a file with Ed25519",                         cfx_sign_run},
    {"verify",          "Verify an Ed25519 signature",                      cfx_verify_run},
    {"aead",            "ChaCha20-Poly1305 AEAD (RFC 8439)",                cfx_aead_run},
    {"argon2",          "Password hashing with Argon2id",                   cfx_argon2_run},
    {"xgcd",            "Extended GCD (Bezout coefficients)",               cfx_xgcd_run}
};

const int cfx_commands_count = sizeof(cfx_commands) / sizeof(cfx_commands[0]);

static void print_help(const char* prog) {
    printf("cfx - Arbitrary precision math utilities\n\n");
    printf("Usage: %s <command> [args...]\n\n", prog);
    printf("Available commands:\n");
    for (int i = 0; i < cfx_commands_count; i++) {
        printf("  %-16s  %s\n", cfx_commands[i].name, cfx_commands[i].description);
    }
    printf("\nRun '%s <command> --help' for more information on a command.\n", prog);
}

int main(int argc, char** argv) {
    if (argc < 2) {
        print_help(argv[0]);
        return 1;
    }

    const char* cmd = argv[1];

    if (strcmp(cmd, "help") == 0 || strcmp(cmd, "--help") == 0 || strcmp(cmd, "-h") == 0) {
        print_help(argv[0]);
        return 0;
    }

    for (int i = 0; i < cfx_commands_count; i++) {
        if (strcmp(cmd, cfx_commands[i].name) == 0) {
            return cfx_commands[i].run(argc - 1, argv + 1);
        }
    }

    fprintf(stderr, "Unknown command: %s\n", cmd);
    fprintf(stderr, "Run '%s help' for a list of available commands.\n", argv[0]);
    return 1;
}
