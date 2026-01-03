//! Raw FFI bindings to the cfx C library.
//!
//! These are low-level bindings. For safe Rust APIs, use the wrapper types
//! in the parent module.

#![allow(non_camel_case_types)]

use std::os::raw::{c_char, c_int, c_uint, c_void};

// ============================================================================
// Types
// ============================================================================

/// The limb type - matches the C library's cfx_limb_t
#[cfg(target_pointer_width = "64")]
pub type cfx_limb_t = u64;
#[cfg(target_pointer_width = "32")]
pub type cfx_limb_t = u32;

/// Arbitrary precision positive integer
#[repr(C)]
#[derive(Debug)]
pub struct cfx_big_t {
    pub limb: *mut cfx_limb_t,
    pub n: usize,
    pub cap: usize,
}

/// SHA-256 context (opaque)
#[repr(C)]
#[derive(Copy, Clone)]
pub struct cfx_sha256_ctx {
    pub opaque: [u8; 128],
}

/// BLAKE2b context (opaque)
#[repr(C)]
#[derive(Copy, Clone)]
pub struct cfx_blake2b_ctx_t {
    pub opaque: [u8; 256],
}

/// BLAKE2s context (opaque)
#[repr(C)]
#[derive(Copy, Clone)]
pub struct cfx_blake2s_ctx_t {
    pub opaque: [u8; 128],
}

/// ChaCha20 RNG context (opaque, 32-byte aligned)
#[repr(C, align(32))]
#[derive(Copy, Clone)]
pub struct cfx_rng_ctx_t {
    pub opaque: [u8; 736], // max size with AVX2
}

// ============================================================================
// Big Integer Functions
// ============================================================================

extern "C" {
    // --- Lifecycle ---
    pub fn cfx_big_init(b: *mut cfx_big_t);
    pub fn cfx_big_free(b: *mut cfx_big_t);
    pub fn cfx_big_copy(dst: *mut cfx_big_t, src: *const cfx_big_t) -> c_int;

    // --- Conversion from ---
    pub fn cfx_big_from_u64(b: *mut cfx_big_t, v: u64) -> c_int;
    pub fn cfx_big_from_limb(b: *mut cfx_big_t, v: cfx_limb_t) -> c_int;
    pub fn cfx_big_from_dec(b: *mut cfx_big_t, s: *const c_char) -> c_int;
    pub fn cfx_big_from_hex(b: *mut cfx_big_t, s: *const c_char) -> c_int;
    pub fn cfx_big_from_bytes_be(out: *mut cfx_big_t, be: *const u8, len: usize) -> c_int;

    // --- Conversion to ---
    pub fn cfx_big_dec_alloc(b: *const cfx_big_t, len_out: *mut usize) -> *mut c_char;
    pub fn cfx_big_hex_alloc(b: *const cfx_big_t, len_out: *mut usize) -> *mut c_char;

    // --- Predicates ---
    pub fn cfx_big_is_zero(b: *const cfx_big_t) -> c_int;
    pub fn cfx_big_is_one(b: *const cfx_big_t) -> c_int;
    pub fn cfx_big_is_even(b: *const cfx_big_t) -> c_int;
    pub fn cfx_big_is_prime(b: *const cfx_big_t) -> c_int;

    // --- Comparison ---
    pub fn cfx_big_eq(a: *const cfx_big_t, b: *const cfx_big_t) -> c_int;
    pub fn cfx_big_cmp(a: *const cfx_big_t, b: *const cfx_big_t) -> c_int;
    pub fn cfx_big_cmp_sm(a: *const cfx_big_t, n: cfx_limb_t) -> c_int;

    // --- Arithmetic ---
    pub fn cfx_big_add_eq(b: *mut cfx_big_t, a: *const cfx_big_t);
    pub fn cfx_big_add_sm_eq(b: *mut cfx_big_t, n: cfx_limb_t);
    pub fn cfx_big_sub_eq(a: *mut cfx_big_t, b: *const cfx_big_t);
    pub fn cfx_big_sub_sm_eq(b: *mut cfx_big_t, n: cfx_limb_t);
    pub fn cfx_big_mul(out: *mut cfx_big_t, a: *const cfx_big_t, b: *const cfx_big_t);
    pub fn cfx_big_mul_eq(b: *mut cfx_big_t, m: *const cfx_big_t);
    pub fn cfx_big_mul_sm_eq(b: *mut cfx_big_t, m: cfx_limb_t);
    pub fn cfx_big_sq_eq(b: *mut cfx_big_t);
    pub fn cfx_big_divrem_eq(b: *mut cfx_big_t, d: *const cfx_big_t, r: *mut cfx_big_t) -> c_int;
    pub fn cfx_big_div_sm_eq(b: *mut cfx_big_t, d: cfx_limb_t) -> cfx_limb_t;
    pub fn cfx_big_mod(out: *mut cfx_big_t, n: *const cfx_big_t, m: *const cfx_big_t) -> c_int;
    pub fn cfx_big_mod_sm(b: *const cfx_big_t, m: cfx_limb_t) -> cfx_limb_t;

    // --- Exponentiation ---
    pub fn cfx_big_exp(out: *mut cfx_big_t, n: *const cfx_big_t, p: *const cfx_big_t);
    pub fn cfx_big_exp_u64(out: *mut cfx_big_t, n: *const cfx_big_t, p: u64);
    pub fn cfx_big_mod_exp(
        out: *mut cfx_big_t,
        n: *const cfx_big_t,
        p: *const cfx_big_t,
        m: *const cfx_big_t,
    );

    // --- Bit operations ---
    pub fn cfx_big_bitlen(b: *const cfx_big_t) -> usize;
    pub fn cfx_big_shl_bits_eq(b: *mut cfx_big_t, s: c_uint);
    pub fn cfx_big_shr_bits_eq(b: *mut cfx_big_t, s: c_uint);
}

// ============================================================================
// Hash Functions
// ============================================================================

extern "C" {
    // --- SHA-256 ---
    pub fn cfx_sha256_init(ctx: *mut cfx_sha256_ctx);
    pub fn cfx_sha256_update(ctx: *mut cfx_sha256_ctx, data: *const u8, len: usize);
    pub fn cfx_sha256_final(ctx: *mut cfx_sha256_ctx, out: *mut u8);

    // --- BLAKE2b ---
    pub fn cfx_blake2b_init(ctx: *mut cfx_blake2b_ctx_t, outlen: usize) -> c_int;
    pub fn cfx_blake2b_init_key(
        ctx: *mut cfx_blake2b_ctx_t,
        outlen: usize,
        key: *const c_void,
        keylen: usize,
    ) -> c_int;
    pub fn cfx_blake2b_update(ctx: *mut cfx_blake2b_ctx_t, data: *const c_void, len: usize) -> c_int;
    pub fn cfx_blake2b_final(ctx: *mut cfx_blake2b_ctx_t, out: *mut c_void) -> c_int;
    pub fn cfx_blake2b(
        out: *mut c_void,
        outlen: usize,
        data: *const c_void,
        datalen: usize,
        key: *const c_void,
        keylen: usize,
    ) -> c_int;

    // --- BLAKE2s ---
    pub fn cfx_blake2s_init(ctx: *mut cfx_blake2s_ctx_t, outlen: usize) -> c_int;
    pub fn cfx_blake2s_init_key(
        ctx: *mut cfx_blake2s_ctx_t,
        outlen: usize,
        key: *const c_void,
        keylen: usize,
    ) -> c_int;
    pub fn cfx_blake2s_update(ctx: *mut cfx_blake2s_ctx_t, data: *const c_void, len: usize) -> c_int;
    pub fn cfx_blake2s_final(ctx: *mut cfx_blake2s_ctx_t, out: *mut c_void) -> c_int;
    pub fn cfx_blake2s(
        out: *mut c_void,
        outlen: usize,
        data: *const c_void,
        datalen: usize,
        key: *const c_void,
        keylen: usize,
    ) -> c_int;
}

// ============================================================================
// Ed25519 Digital Signatures
// ============================================================================

extern "C" {
    pub fn cfx_ed25519_create_keypair(pk: *mut u8, sk: *mut u8, seed: *const u8);
    pub fn cfx_ed25519_sign(sig: *mut u8, msg: *const u8, msg_len: usize, sk: *const u8);
    pub fn cfx_ed25519_verify(sig: *const u8, msg: *const u8, msg_len: usize, pk: *const u8) -> c_int;
    pub fn cfx_ed25519_get_public_key(pk: *mut u8, sk: *const u8);
}

// ============================================================================
// X25519 Key Exchange
// ============================================================================

extern "C" {
    pub fn cfx_x25519(shared: *mut u8, scalar: *const u8, point: *const u8) -> c_int;
    pub fn cfx_x25519_base(public_key: *mut u8, private_key: *const u8);
    pub fn cfx_x25519_scalarmult(out: *mut u8, scalar: *const u8, point: *const u8);
}

// ============================================================================
// Random Number Generation
// ============================================================================

extern "C" {
    pub fn cfx_chacha20_rng_init(st: *mut cfx_rng_ctx_t, seed: u32);
    pub fn cfx_chacha20_rng(ctx: *mut c_void, out: *mut u8, len: usize) -> c_int;

    pub fn cfx_srand(seed: u32);
    pub fn cfx_rand() -> c_int;
    pub fn cfx_urand() -> u32;
    pub fn cfx_rand_bytes(buf: *mut c_void, len: usize);

    pub fn cfx_srand_os();
    pub fn cfx_rand_os() -> u32;
    pub fn cfx_rand_bytes_os(buf: *mut c_void, len: usize);
}

// ============================================================================
// Memory (for freeing allocated strings)
// ============================================================================

#[cfg(unix)]
pub unsafe fn libc_free(ptr: *mut c_void) {
    extern "C" {
        fn free(ptr: *mut c_void);
    }
    free(ptr);
}

#[cfg(windows)]
pub unsafe fn libc_free(ptr: *mut c_void) {
    extern "C" {
        fn free(ptr: *mut c_void);
    }
    free(ptr);
}
