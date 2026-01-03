//! Ed25519 digital signatures (RFC 8032).
//!
//! Ed25519 is a deterministic signature scheme providing 128-bit security.
//! All operations are constant-time with respect to secret data.
//!
//! # Key Sizes
//! - Seed: 32 bytes (random)
//! - Secret key: 64 bytes (seed + public key bundled)
//! - Public key: 32 bytes
//! - Signature: 64 bytes
//!
//! # Example
//! ```
//! use cfx::ed25519::{Ed25519Keypair, Ed25519Signature};
//!
//! // Generate keypair from a 32-byte seed
//! let seed = [0u8; 32]; // Use proper random bytes in production!
//! let keypair = Ed25519Keypair::from_seed(&seed);
//!
//! // Sign a message
//! let message = b"Hello, Ed25519!";
//! let signature = keypair.sign(message);
//!
//! // Verify the signature
//! assert!(signature.verify(message, keypair.public_key()));
//! ```

use crate::ffi;

/// Ed25519 keypair (secret key + public key).
///
/// The secret key is 64 bytes: the 32-byte seed concatenated with the
/// 32-byte public key. This bundling prevents a class of attacks where
/// signing with mismatched keys leaks the private scalar.
#[derive(Clone)]
pub struct Ed25519Keypair {
    secret_key: [u8; 64],
    public_key: [u8; 32],
}

impl Ed25519Keypair {
    /// Creates a keypair from a 32-byte seed.
    ///
    /// The seed should be generated using a cryptographically secure RNG.
    pub fn from_seed(seed: &[u8; 32]) -> Self {
        let mut pk = [0u8; 32];
        let mut sk = [0u8; 64];
        unsafe {
            ffi::cfx_ed25519_create_keypair(pk.as_mut_ptr(), sk.as_mut_ptr(), seed.as_ptr());
        }
        Ed25519Keypair {
            secret_key: sk,
            public_key: pk,
        }
    }

    /// Returns the 32-byte public key.
    #[inline]
    pub fn public_key(&self) -> &[u8; 32] {
        &self.public_key
    }

    /// Returns the 64-byte secret key (seed || public_key).
    ///
    /// # Security
    /// Handle with care - this contains the private key material.
    #[inline]
    pub fn secret_key(&self) -> &[u8; 64] {
        &self.secret_key
    }

    /// Signs a message, returning a 64-byte signature.
    pub fn sign(&self, message: &[u8]) -> Ed25519Signature {
        let mut sig = [0u8; 64];
        unsafe {
            ffi::cfx_ed25519_sign(
                sig.as_mut_ptr(),
                message.as_ptr(),
                message.len(),
                self.secret_key.as_ptr(),
            );
        }
        Ed25519Signature(sig)
    }

    /// Reconstructs keypair from a 64-byte secret key.
    ///
    /// The secret key format is seed (32 bytes) || public_key (32 bytes).
    pub fn from_secret_key(sk: &[u8; 64]) -> Self {
        let mut pk = [0u8; 32];
        unsafe {
            ffi::cfx_ed25519_get_public_key(pk.as_mut_ptr(), sk.as_ptr());
        }
        Ed25519Keypair {
            secret_key: *sk,
            public_key: pk,
        }
    }
}

/// Ed25519 signature (64 bytes).
#[derive(Clone, Copy, PartialEq, Eq)]
pub struct Ed25519Signature(pub [u8; 64]);

impl Ed25519Signature {
    /// Creates a signature from raw bytes.
    pub fn from_bytes(bytes: [u8; 64]) -> Self {
        Ed25519Signature(bytes)
    }

    /// Returns the signature as a byte array.
    #[inline]
    pub fn as_bytes(&self) -> &[u8; 64] {
        &self.0
    }

    /// Verifies this signature against a message and public key.
    ///
    /// Returns `true` if the signature is valid.
    pub fn verify(&self, message: &[u8], public_key: &[u8; 32]) -> bool {
        let rc = unsafe {
            ffi::cfx_ed25519_verify(
                self.0.as_ptr(),
                message.as_ptr(),
                message.len(),
                public_key.as_ptr(),
            )
        };
        rc == 0 // 0 = valid
    }
}

impl std::fmt::Debug for Ed25519Signature {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Ed25519Signature({:02x?}...)", &self.0[..8])
    }
}

/// Standalone signature verification.
///
/// Convenience function that doesn't require constructing an Ed25519Signature.
pub fn verify(signature: &[u8; 64], message: &[u8], public_key: &[u8; 32]) -> bool {
    Ed25519Signature(*signature).verify(message, public_key)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sign_verify() {
        let seed = [42u8; 32];
        let keypair = Ed25519Keypair::from_seed(&seed);
        let message = b"test message";

        let sig = keypair.sign(message);
        assert!(sig.verify(message, keypair.public_key()));
    }

    #[test]
    fn test_wrong_message_fails() {
        let seed = [42u8; 32];
        let keypair = Ed25519Keypair::from_seed(&seed);

        let sig = keypair.sign(b"correct message");
        assert!(!sig.verify(b"wrong message", keypair.public_key()));
    }

    #[test]
    fn test_wrong_pubkey_fails() {
        let keypair1 = Ed25519Keypair::from_seed(&[1u8; 32]);
        let keypair2 = Ed25519Keypair::from_seed(&[2u8; 32]);
        let message = b"test";

        let sig = keypair1.sign(message);
        assert!(!sig.verify(message, keypair2.public_key()));
    }

    #[test]
    fn test_from_secret_key() {
        let seed = [99u8; 32];
        let keypair1 = Ed25519Keypair::from_seed(&seed);
        let keypair2 = Ed25519Keypair::from_secret_key(keypair1.secret_key());

        assert_eq!(keypair1.public_key(), keypair2.public_key());

        // Both should produce same signature
        let sig1 = keypair1.sign(b"test");
        let sig2 = keypair2.sign(b"test");
        assert_eq!(sig1.as_bytes(), sig2.as_bytes());
    }

    #[test]
    fn test_deterministic() {
        let seed = [7u8; 32];
        let keypair = Ed25519Keypair::from_seed(&seed);
        let message = b"deterministic test";

        // Same seed + message should always produce same signature
        let sig1 = keypair.sign(message);
        let sig2 = keypair.sign(message);
        assert_eq!(sig1.as_bytes(), sig2.as_bytes());
    }
}
