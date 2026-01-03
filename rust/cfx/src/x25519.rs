//! X25519 Diffie-Hellman key exchange (RFC 7748).
//!
//! X25519 provides 128-bit security for key agreement. All operations
//! are constant-time.
//!
//! # Example
//! ```
//! use cfx::x25519::X25519PrivateKey;
//!
//! // Alice generates her keypair
//! let alice_private = X25519PrivateKey::from_bytes([1u8; 32]); // Use random bytes!
//! let alice_public = alice_private.public_key();
//!
//! // Bob generates his keypair
//! let bob_private = X25519PrivateKey::from_bytes([2u8; 32]);
//! let bob_public = bob_private.public_key();
//!
//! // Both compute the same shared secret
//! let alice_shared = alice_private.diffie_hellman(&bob_public).unwrap();
//! let bob_shared = bob_private.diffie_hellman(&alice_public).unwrap();
//! assert_eq!(alice_shared.as_bytes(), bob_shared.as_bytes());
//! ```

use crate::error::{Error, Result};
use crate::ffi;

/// X25519 private key (32 bytes).
///
/// The private key is clamped internally during operations:
/// - Bits 0, 1, 2 are cleared
/// - Bit 254 is set
/// - Bit 255 is cleared
#[derive(Clone)]
pub struct X25519PrivateKey([u8; 32]);

impl X25519PrivateKey {
    /// Creates a private key from 32 random bytes.
    ///
    /// The bytes should be generated using a CSPRNG.
    /// Clamping is applied automatically during DH operations.
    pub fn from_bytes(bytes: [u8; 32]) -> Self {
        X25519PrivateKey(bytes)
    }

    /// Returns the raw private key bytes.
    #[inline]
    pub fn as_bytes(&self) -> &[u8; 32] {
        &self.0
    }

    /// Computes the corresponding public key.
    ///
    /// This performs scalar multiplication: public = private * basepoint
    pub fn public_key(&self) -> X25519PublicKey {
        let mut pk = [0u8; 32];
        unsafe {
            ffi::cfx_x25519_base(pk.as_mut_ptr(), self.0.as_ptr());
        }
        X25519PublicKey(pk)
    }

    /// Performs Diffie-Hellman key exchange.
    ///
    /// Computes: shared = self * their_public
    ///
    /// Returns an error if the result is all zeros (indicates a low-order point,
    /// which would be a security issue).
    pub fn diffie_hellman(&self, their_public: &X25519PublicKey) -> Result<X25519SharedSecret> {
        let mut shared = [0u8; 32];
        let rc = unsafe {
            ffi::cfx_x25519(shared.as_mut_ptr(), self.0.as_ptr(), their_public.0.as_ptr())
        };
        if rc != 0 {
            return Err(Error::AllZeroResult);
        }
        Ok(X25519SharedSecret(shared))
    }
}

/// X25519 public key (32 bytes).
///
/// This is a point on Curve25519 in Montgomery form (x-coordinate only).
#[derive(Clone, Copy, PartialEq, Eq)]
pub struct X25519PublicKey(pub [u8; 32]);

impl X25519PublicKey {
    /// Creates a public key from raw bytes.
    pub fn from_bytes(bytes: [u8; 32]) -> Self {
        X25519PublicKey(bytes)
    }

    /// Returns the public key as bytes.
    #[inline]
    pub fn as_bytes(&self) -> &[u8; 32] {
        &self.0
    }
}

impl std::fmt::Debug for X25519PublicKey {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "X25519PublicKey({:02x?}...)", &self.0[..8])
    }
}

/// X25519 shared secret (32 bytes).
///
/// The result of a successful Diffie-Hellman exchange. This should typically
/// be passed through a KDF (like HKDF) before use as a symmetric key.
#[derive(Clone)]
pub struct X25519SharedSecret([u8; 32]);

impl X25519SharedSecret {
    /// Returns the shared secret as bytes.
    #[inline]
    pub fn as_bytes(&self) -> &[u8; 32] {
        &self.0
    }
}

impl Drop for X25519SharedSecret {
    fn drop(&mut self) {
        // Zero out the secret on drop
        self.0.iter_mut().for_each(|b| *b = 0);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_key_agreement() {
        // Alice
        let alice_priv = X25519PrivateKey::from_bytes([
            0x77, 0x07, 0x6d, 0x0a, 0x73, 0x18, 0xa5, 0x7d,
            0x3c, 0x16, 0xc1, 0x72, 0x51, 0xb2, 0x66, 0x45,
            0xdf, 0x4c, 0x2f, 0x87, 0xeb, 0xc0, 0x99, 0x2a,
            0xb1, 0x77, 0xfb, 0xa5, 0x1d, 0xb9, 0x2c, 0x2a,
        ]);
        let alice_pub = alice_priv.public_key();

        // Bob
        let bob_priv = X25519PrivateKey::from_bytes([
            0x5d, 0xab, 0x08, 0x7e, 0x62, 0x4a, 0x8a, 0x4b,
            0x79, 0xe1, 0x7f, 0x8b, 0x83, 0x80, 0x0e, 0xe6,
            0x6f, 0x3b, 0xb1, 0x29, 0x26, 0x18, 0xb6, 0xfd,
            0x1c, 0x2f, 0x8b, 0x27, 0xff, 0x88, 0xe0, 0xeb,
        ]);
        let bob_pub = bob_priv.public_key();

        // Both compute same shared secret
        let alice_shared = alice_priv.diffie_hellman(&bob_pub).unwrap();
        let bob_shared = bob_priv.diffie_hellman(&alice_pub).unwrap();

        assert_eq!(alice_shared.as_bytes(), bob_shared.as_bytes());
    }

    #[test]
    fn test_public_key_derivation() {
        // From RFC 7748 test vector
        let private = X25519PrivateKey::from_bytes([
            0x77, 0x07, 0x6d, 0x0a, 0x73, 0x18, 0xa5, 0x7d,
            0x3c, 0x16, 0xc1, 0x72, 0x51, 0xb2, 0x66, 0x45,
            0xdf, 0x4c, 0x2f, 0x87, 0xeb, 0xc0, 0x99, 0x2a,
            0xb1, 0x77, 0xfb, 0xa5, 0x1d, 0xb9, 0x2c, 0x2a,
        ]);
        let public = private.public_key();

        let expected = [
            0x85, 0x20, 0xf0, 0x09, 0x89, 0x30, 0xa7, 0x54,
            0x74, 0x8b, 0x7d, 0xdc, 0xb4, 0x3e, 0xf7, 0x5a,
            0x0d, 0xbf, 0x3a, 0x0d, 0x26, 0x38, 0x1a, 0xf4,
            0xeb, 0xa4, 0xa9, 0x8e, 0xaa, 0x9b, 0x4e, 0x6a,
        ];
        assert_eq!(public.as_bytes(), &expected);
    }
}
