//! Cryptographic hash functions.
//!
//! This module provides SHA-256, BLAKE2b, and BLAKE2s hash functions.
//!
//! # Example
//! ```
//! use cfx::hash::{Sha256, Blake2b};
//!
//! // SHA-256 (returns 32 bytes)
//! let hash = Sha256::hash(b"hello world");
//! assert_eq!(hash.len(), 32);
//!
//! // BLAKE2b (one-shot, returns 64 bytes)
//! let hash = Blake2b::hash(b"hello world");
//! assert_eq!(hash.len(), 64);
//! ```

use crate::error::{Error, Result};
use crate::ffi;
use std::ptr;

// SHA-256

/// SHA-256 hasher with streaming support.
///
/// # Example
/// ```
/// use cfx::hash::Sha256;
///
/// let mut hasher = Sha256::new();
/// hasher.update(b"hello ");
/// hasher.update(b"world");
/// let hash = hasher.finalize();
/// ```
pub struct Sha256 {
    ctx: ffi::cfx_sha256_ctx,
}

impl Sha256 {
    /// Creates a new SHA-256 hasher.
    pub fn new() -> Self {
        let mut ctx = ffi::cfx_sha256_ctx { opaque: [0u8; 128] };
        unsafe { ffi::cfx_sha256_init(&mut ctx) };
        Sha256 { ctx }
    }

    /// Feeds data into the hasher.
    pub fn update(&mut self, data: &[u8]) {
        unsafe {
            ffi::cfx_sha256_update(&mut self.ctx, data.as_ptr(), data.len());
        }
    }

    /// Finalizes the hash and returns the 32-byte digest.
    ///
    /// Consumes the hasher. To hash more data, create a new instance.
    pub fn finalize(mut self) -> [u8; 32] {
        let mut out = [0u8; 32];
        unsafe { ffi::cfx_sha256_final(&mut self.ctx, out.as_mut_ptr()) };
        out
    }

    /// One-shot hash function.
    pub fn hash(data: &[u8]) -> [u8; 32] {
        let mut hasher = Sha256::new();
        hasher.update(data);
        hasher.finalize()
    }
}

impl Default for Sha256 {
    fn default() -> Self {
        Self::new()
    }
}

// BLAKE2b

/// BLAKE2b hasher (64-bit optimized, up to 64-byte output).
///
/// BLAKE2b is optimized for 64-bit platforms and supports:
/// - Variable output length (1-64 bytes, default 64)
/// - Optional key for MAC mode (up to 64 bytes)
///
/// # Example
/// ```
/// use cfx::hash::Blake2b;
///
/// // Default 64-byte output
/// let hash = Blake2b::hash(b"hello world");
///
/// // Custom output length
/// let mut hasher = Blake2b::new(32).unwrap();
/// hasher.update(b"hello world");
/// let hash = hasher.finalize();
/// assert_eq!(hash.len(), 32);
/// ```
pub struct Blake2b {
    ctx: ffi::cfx_blake2b_ctx_t,
    outlen: usize,
}

impl Blake2b {
    /// Creates a new BLAKE2b hasher with specified output length (1-64 bytes).
    pub fn new(outlen: usize) -> Result<Self> {
        if outlen == 0 || outlen > 64 {
            return Err(Error::InvalidOutputLength);
        }
        let mut ctx = ffi::cfx_blake2b_ctx_t { opaque: [0u8; 256] };
        let rc = unsafe { ffi::cfx_blake2b_init(&mut ctx, outlen) };
        if rc != 0 {
            return Err(Error::HashError);
        }
        Ok(Blake2b { ctx, outlen })
    }

    /// Creates a new keyed BLAKE2b hasher (MAC mode).
    ///
    /// Key must be 1-64 bytes.
    pub fn with_key(outlen: usize, key: &[u8]) -> Result<Self> {
        if outlen == 0 || outlen > 64 {
            return Err(Error::InvalidOutputLength);
        }
        if key.is_empty() || key.len() > 64 {
            return Err(Error::InvalidKeyLength);
        }
        let mut ctx = ffi::cfx_blake2b_ctx_t { opaque: [0u8; 256] };
        let rc = unsafe {
            ffi::cfx_blake2b_init_key(&mut ctx, outlen, key.as_ptr() as *const _, key.len())
        };
        if rc != 0 {
            return Err(Error::HashError);
        }
        Ok(Blake2b { ctx, outlen })
    }

    /// Feeds data into the hasher.
    pub fn update(&mut self, data: &[u8]) {
        unsafe {
            ffi::cfx_blake2b_update(&mut self.ctx, data.as_ptr() as *const _, data.len());
        }
    }

    /// Finalizes and returns the hash.
    pub fn finalize(mut self) -> Vec<u8> {
        let mut out = vec![0u8; self.outlen];
        unsafe { ffi::cfx_blake2b_final(&mut self.ctx, out.as_mut_ptr() as *mut _) };
        out
    }

    /// One-shot hash with default 64-byte output.
    pub fn hash(data: &[u8]) -> [u8; 64] {
        let mut out = [0u8; 64];
        unsafe {
            ffi::cfx_blake2b(
                out.as_mut_ptr() as *mut _,
                64,
                data.as_ptr() as *const _,
                data.len(),
                ptr::null(),
                0,
            );
        }
        out
    }

    /// One-shot keyed hash (MAC).
    pub fn mac(key: &[u8], data: &[u8], outlen: usize) -> Result<Vec<u8>> {
        if outlen == 0 || outlen > 64 {
            return Err(Error::InvalidOutputLength);
        }
        if key.is_empty() || key.len() > 64 {
            return Err(Error::InvalidKeyLength);
        }
        let mut out = vec![0u8; outlen];
        let rc = unsafe {
            ffi::cfx_blake2b(
                out.as_mut_ptr() as *mut _,
                outlen,
                data.as_ptr() as *const _,
                data.len(),
                key.as_ptr() as *const _,
                key.len(),
            )
        };
        if rc != 0 {
            return Err(Error::HashError);
        }
        Ok(out)
    }
}

// BLAKE2s

/// BLAKE2s hasher (32-bit optimized, up to 32-byte output).
///
/// BLAKE2s is optimized for 32-bit platforms and embedded systems.
/// Supports variable output length (1-32 bytes) and optional keying.
pub struct Blake2s {
    ctx: ffi::cfx_blake2s_ctx_t,
    outlen: usize,
}

impl Blake2s {
    /// Creates a new BLAKE2s hasher with specified output length (1-32 bytes).
    pub fn new(outlen: usize) -> Result<Self> {
        if outlen == 0 || outlen > 32 {
            return Err(Error::InvalidOutputLength);
        }
        let mut ctx = ffi::cfx_blake2s_ctx_t { opaque: [0u8; 128] };
        let rc = unsafe { ffi::cfx_blake2s_init(&mut ctx, outlen) };
        if rc != 0 {
            return Err(Error::HashError);
        }
        Ok(Blake2s { ctx, outlen })
    }

    /// Creates a new keyed BLAKE2s hasher (MAC mode).
    pub fn with_key(outlen: usize, key: &[u8]) -> Result<Self> {
        if outlen == 0 || outlen > 32 {
            return Err(Error::InvalidOutputLength);
        }
        if key.is_empty() || key.len() > 32 {
            return Err(Error::InvalidKeyLength);
        }
        let mut ctx = ffi::cfx_blake2s_ctx_t { opaque: [0u8; 128] };
        let rc = unsafe {
            ffi::cfx_blake2s_init_key(&mut ctx, outlen, key.as_ptr() as *const _, key.len())
        };
        if rc != 0 {
            return Err(Error::HashError);
        }
        Ok(Blake2s { ctx, outlen })
    }

    /// Feeds data into the hasher.
    pub fn update(&mut self, data: &[u8]) {
        unsafe {
            ffi::cfx_blake2s_update(&mut self.ctx, data.as_ptr() as *const _, data.len());
        }
    }

    /// Finalizes and returns the hash.
    pub fn finalize(mut self) -> Vec<u8> {
        let mut out = vec![0u8; self.outlen];
        unsafe { ffi::cfx_blake2s_final(&mut self.ctx, out.as_mut_ptr() as *mut _) };
        out
    }

    /// One-shot hash with default 32-byte output.
    pub fn hash(data: &[u8]) -> [u8; 32] {
        let mut out = [0u8; 32];
        unsafe {
            ffi::cfx_blake2s(
                out.as_mut_ptr() as *mut _,
                32,
                data.as_ptr() as *const _,
                data.len(),
                ptr::null(),
                0,
            );
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sha256_empty() {
        let hash = Sha256::hash(b"");
        // SHA-256 of empty string
        let expected = [
            0xe3, 0xb0, 0xc4, 0x42, 0x98, 0xfc, 0x1c, 0x14,
            0x9a, 0xfb, 0xf4, 0xc8, 0x99, 0x6f, 0xb9, 0x24,
            0x27, 0xae, 0x41, 0xe4, 0x64, 0x9b, 0x93, 0x4c,
            0xa4, 0x95, 0x99, 0x1b, 0x78, 0x52, 0xb8, 0x55,
        ];
        assert_eq!(hash, expected);
    }

    #[test]
    fn test_sha256_streaming() {
        let mut hasher = Sha256::new();
        hasher.update(b"hello ");
        hasher.update(b"world");
        let hash1 = hasher.finalize();

        let hash2 = Sha256::hash(b"hello world");
        assert_eq!(hash1, hash2);
    }

    #[test]
    fn test_blake2b_lengths() {
        // Valid lengths
        assert!(Blake2b::new(1).is_ok());
        assert!(Blake2b::new(32).is_ok());
        assert!(Blake2b::new(64).is_ok());

        // Invalid lengths
        assert!(Blake2b::new(0).is_err());
        assert!(Blake2b::new(65).is_err());
    }

    #[test]
    fn test_blake2s_lengths() {
        assert!(Blake2s::new(1).is_ok());
        assert!(Blake2s::new(32).is_ok());
        assert!(Blake2s::new(0).is_err());
        assert!(Blake2s::new(33).is_err());
    }

    #[test]
    fn test_blake2b_keyed() {
        let key = b"supersecretkey!!";
        let mac = Blake2b::mac(key, b"hello world", 32).unwrap();
        assert_eq!(mac.len(), 32);

        // Different key should give different MAC
        let mac2 = Blake2b::mac(b"differentkey1234", b"hello world", 32).unwrap();
        assert_ne!(mac, mac2);
    }
}
