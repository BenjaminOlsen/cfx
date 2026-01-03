//! Error types for cfx operations.

use std::fmt;

/// Errors that can occur in cfx operations.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Error {
    /// The provided key length is invalid.
    InvalidKeyLength,
    /// The requested output length is invalid (e.g., BLAKE2 outlen > 64).
    InvalidOutputLength,
    /// Signature verification failed.
    VerificationFailed,
    /// X25519 produced an all-zero result (indicates a low-order point).
    AllZeroResult,
    /// A hash operation failed.
    HashError,
    /// Memory allocation failed.
    AllocationFailed,
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::InvalidKeyLength => write!(f, "invalid key length"),
            Error::InvalidOutputLength => write!(f, "invalid output length"),
            Error::VerificationFailed => write!(f, "signature verification failed"),
            Error::AllZeroResult => write!(f, "operation produced all-zero result"),
            Error::HashError => write!(f, "hash operation failed"),
            Error::AllocationFailed => write!(f, "memory allocation failed"),
        }
    }
}

impl std::error::Error for Error {}

/// Result type for cfx operations.
pub type Result<T> = std::result::Result<T, Error>;
