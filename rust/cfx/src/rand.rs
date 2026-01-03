//! Random number generation.
//!
//! This module provides ChaCha20-based cryptographically secure random
//! number generation, as well as access to the OS RNG.
//!
//! # Example
//! ```
//! use cfx::rand::{ChaChaRng, random_bytes};
//!
//! // Seeded deterministic RNG (for reproducibility)
//! let mut rng = ChaChaRng::from_seed(12345);
//! let mut buf = [0u8; 32];
//! rng.fill_bytes(&mut buf);
//!
//! // OS-seeded random bytes
//! random_bytes(&mut buf);
//! ```

use crate::ffi;
use std::os::raw::c_void;

/// ChaCha20-based random number generator.
///
/// This RNG uses the ChaCha20 stream cipher to generate random bytes.
/// It is deterministic given the same seed.
pub struct ChaChaRng {
    ctx: ffi::cfx_rng_ctx_t,
}

impl ChaChaRng {
    /// Creates a new RNG from a 32-bit seed.
    ///
    /// The same seed always produces the same sequence.
    pub fn from_seed(seed: u32) -> Self {
        let mut ctx = ffi::cfx_rng_ctx_t { opaque: [0u8; 736] };
        unsafe {
            ffi::cfx_chacha20_rng_init(&mut ctx, seed);
        }
        ChaChaRng { ctx }
    }

    /// Fills a buffer with random bytes.
    pub fn fill_bytes(&mut self, dest: &mut [u8]) {
        unsafe {
            ffi::cfx_chacha20_rng(&mut self.ctx as *mut _ as *mut c_void, dest.as_mut_ptr(), dest.len());
        }
    }

    /// Returns a random u32.
    pub fn next_u32(&mut self) -> u32 {
        let mut buf = [0u8; 4];
        self.fill_bytes(&mut buf);
        u32::from_le_bytes(buf)
    }

    /// Returns a random u64.
    pub fn next_u64(&mut self) -> u64 {
        let mut buf = [0u8; 8];
        self.fill_bytes(&mut buf);
        u64::from_le_bytes(buf)
    }
}

// ============================================================================
// Global RNG functions (use internal global state)
// ============================================================================

/// Seeds the global RNG with a 32-bit value.
///
/// Call this before using `rand()` or `rand_bytes()`.
pub fn seed(s: u32) {
    unsafe { ffi::cfx_srand(s) };
}

/// Returns a random integer in range [0, 0x7FFFFFFF].
pub fn rand() -> i32 {
    unsafe { ffi::cfx_rand() }
}

/// Returns a random u32 in full range [0, 0xFFFFFFFF].
pub fn urand() -> u32 {
    unsafe { ffi::cfx_urand() }
}

/// Fills a buffer with random bytes using the global RNG.
///
/// Uses the seed set by `seed()`.
pub fn rand_bytes(dest: &mut [u8]) {
    unsafe { ffi::cfx_rand_bytes(dest.as_mut_ptr() as *mut c_void, dest.len()) };
}

// ============================================================================
// OS RNG functions
// ============================================================================

/// Seeds the global RNG from the OS entropy source.
///
/// Uses /dev/urandom on Unix or CryptGenRandom on Windows.
pub fn seed_from_os() {
    unsafe { ffi::cfx_srand_os() };
}

/// Returns a random u32 from the OS entropy source.
pub fn rand_os() -> u32 {
    unsafe { ffi::cfx_rand_os() }
}

/// Fills a buffer with random bytes from the OS entropy source.
///
/// This is the recommended function for cryptographic purposes.
pub fn random_bytes(dest: &mut [u8]) {
    unsafe { ffi::cfx_rand_bytes_os(dest.as_mut_ptr() as *mut c_void, dest.len()) };
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_chacha_deterministic() {
        let mut rng1 = ChaChaRng::from_seed(42);
        let mut rng2 = ChaChaRng::from_seed(42);

        let mut buf1 = [0u8; 64];
        let mut buf2 = [0u8; 64];

        rng1.fill_bytes(&mut buf1);
        rng2.fill_bytes(&mut buf2);

        assert_eq!(buf1, buf2);
    }

    #[test]
    fn test_chacha_different_seeds() {
        let mut rng1 = ChaChaRng::from_seed(1);
        let mut rng2 = ChaChaRng::from_seed(2);

        let mut buf1 = [0u8; 32];
        let mut buf2 = [0u8; 32];

        rng1.fill_bytes(&mut buf1);
        rng2.fill_bytes(&mut buf2);

        assert_ne!(buf1, buf2);
    }

    #[test]
    fn test_next_u32_u64() {
        let mut rng = ChaChaRng::from_seed(123);
        let a = rng.next_u32();
        let b = rng.next_u32();
        assert_ne!(a, b); // Extremely unlikely to be equal

        let c = rng.next_u64();
        let d = rng.next_u64();
        assert_ne!(c, d);
    }

    #[test]
    fn test_os_random_not_all_zeros() {
        let mut buf = [0u8; 32];
        random_bytes(&mut buf);
        // Extremely unlikely to be all zeros
        assert!(buf.iter().any(|&b| b != 0));
    }
}
